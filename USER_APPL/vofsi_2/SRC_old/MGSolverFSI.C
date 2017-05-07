// ===============================================================
// --------------   NAVIER-STOKES system [FSI_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef FSI_EQUATIONS
// ==============================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverFSI.h"       // Navier-Stokes class header file
#ifdef TBK_EQUATIONS
#include "MGSolverTBK.h" 
#endif

// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif
// ==============================================================
#ifdef HAVE_MED
#include "vof_config.h"
#endif
// #define NON_LINEAR_ITER (0)

// void  my1_add(DenseMatrixM &);
//  void bc_Dirichlet(DenseMatrixM KeM, DenseVectorM &FeM,int _FSI_idx,int  _nFSIdim,int dir_maxnormal,int dir,
//    int bc_in[], int bc_vol[],int sur_toply[], int el_ndof[],int nvars[],double ub[],
//     double u_old[], double normal[],double Ipenalty
//   );
//  void MGSolFSI::bc_Dirichlet(/*DenseMatrixM KeM, DenseVectorM &FeM,int _FSI_idx,int  _nFSIdim,int dir_maxnormal,int dir,*/
//    int bc_in[]/*, int bc_vol[],int sur_toply[], int el_ndof[]*/,int nvars[]  /*,double ub[],*/
//     double u_old[], double normal[],double Ipenalty
//   );

// ==================================================================
/// This routine constructs the FSI class:
MGSolFSI::MGSolFSI(
  MGEquationsSystem& mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
//
/*===================================================================================================*/
/*                                    A) Reading parameters                                          */
/*===================================================================================================*/
//
  _offset(_mgmesh._NoNodes[_NoLevels-1]),             // mesh nodes (top level)
  _dt(_mgutils.get_par("dt")),                        // parameter  dt
  _uref(mg_equations_map_in.get_par("Uref")),         // parameter  u reference
  _lref(mg_equations_map_in.get_par("Lref")),         // parameter  l reference
  _rhof(mg_equations_map_in.get_par("rho0")),         // parameter density
  _T_nc(mg_equations_map_in.get_par("T_nc")),
  _muf(mg_equations_map_in.get_par("mu0")),
   _rhos(mg_equations_map_in.get_par("rhos")),    // parameter density
  _ni(mg_equations_map_in.get_par("nis")),
  _Emod(mg_equations_map_in.get_par("Es")),
  _hs(mg_equations_map_in.get_par("hs"))
  {          // parameter viscosity
//
  /*===================================================================================================*/
  /*                                    B) Setting class variables                                     */
  /*===================================================================================================*/


  _nFSIdim=DIMENSION;
  _dir=0;
 _pres_order=(_nvars[0]>0)? 0:1;
  // class equation ---------------------------------------------------------------------
  _FSI_idx=-1; //  navier-stokes-equations
  _T_idx=-1;   //  energy-equations
  _K_idx=-1;   //  kappa-equations
  _DS_idx=-1; 
  // class variable names ------------------------------------------------------------
#if FSI_EQUATIONS==2       //   segregated ( P in NSP_EQUATIONS)
  if(!varname_in.compare("u")) { _dir=0; }   // u-equation
  if(!varname_in.compare("v")) { _dir=1; }   // v-equation
  if(!varname_in.compare("w")) { _dir=2; }   // w-equation
  _var_names[0]=varname_in;   _refvalue[0]=_uref;
#else
  _var_names[_nFSIdim-1]="w"; _refvalue[_nFSIdim-1]=_uref; // velocity 3D
  _var_names[0]="u";   _refvalue[0]=_uref; // velocity 2D
  _var_names[1]="v";   _refvalue[1]=_uref; // velocity 2D
#if FSI_EQUATIONS==1       //  coupled  (+P)
  _var_names[_nFSIdim]="p";
  _refvalue[_nFSIdim]=_rhof*_uref*_uref;  // pressure
#endif
#endif
//
  /*===================================================================================================*/
  /*                                    C) Setting solver type                                         */
  /*===================================================================================================*/
//
  for(int l=0; l<_NoLevels; l++) { _solver[l]->set_solver_type(SOLVER_FSI); }
//
  /*===================================================================================================*/
  /*                                    D) Setting non_dimensional parameters                           */
  /*===================================================================================================*/
//
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
  _IFr=9.81*_lref/(_uref*_uref);          // Froud number
  _dirg[0] = mg_equations_map_in.get_par("dirgx");      // x-gravity
  _dirg[1] = mg_equations_map_in.get_par("dirgy");      // y-gravity
  _dirg[2] = mg_equations_map_in.get_par("dirgz");      // z-gravity
//   _beta_mat = mg_equations_map_in.get_par("beta_mat");
//   _grav_mod = mg_equations_map_in.get_par("grav");      // gravity modulus for buoyancy term

FSI_parameter.set_param();
  _y_bcout=sqrt(2*0.09*_IRe/80.); 
   _lambda=(_ni*_Emod)/(_rhos*(1+_ni)*(1-2*_ni));//.001;
  _mus=_Emod/(_rhos*2*(1+_ni)*_lref*_uref);

//   std::cout << "lambda  " << _lambda << " mus " <<  _mus << "\n";
  
  return;
} /******************************************************************************************************/
//
//

// #include "MGsolverFSI_vol.h"

// ====================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================
void  MGSolFSI::GenMatRhs(const double time, const int
                         Level,const  int mode) {
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

  int iaxisim=0;// axysimmetry
#ifdef  AXISYM
  iaxisim=1; // axysimmetry
#endif
  const double euler_impl=1;  // 1 for full implicit, 0.5 for Crank-Nicholson, 0 for full explicit
//
  /*===================================================================================================*/
  /*                                    A) Set up                                                      */
  /*===================================================================================================*/
//
  /*-------------------------- Geometry -------------------------------------------------------------*/
//
//   const int   _nFSIdim = DIMENSION;                                           //dimension
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
//   double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  double     normal[DIMENSION];                                          // normal to the boundary
  double     mu_m;
//
  /*-------------------------- Gauss integration -----------------------------------------------------*/
//
  int  el_ngauss = _fe[2]->_NoGauss1[ _nFSIdim-1];                //elem gauss points
  const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points
  double det2,det1,JxW_g2,JxW_g1;           // Jac, Jac*w Jacobean
  double dphijdx_g2[DIMENSION],dphijdx_g1[DIMENSION];
  double dphiidx_g2[DIMENSION],dphiidx_g1[DIMENSION];
//   double dphiidx_g[3][DIMENSION]; // global derivatives at g point

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] -----------------------------------
  int el_ndof[3];  el_ndof[0]=NDOF_K;  int elb_ndof[3];  elb_ndof[0]=0;
  if(_nvars[0]>0) elb_ndof[0]=1; // number of el dofs
  int el_mat_nrows =0;                                            // number of mat rows (dofs)
  for(int ideg=1; ideg<3; ideg++) {                               //     ...
    el_ndof[ideg]=((_nvars[ideg]>0)?    _fe[ideg]->_NoShape[ _nFSIdim-1]:0);                    //   computing
    elb_ndof[ideg]=((_nvars[ideg]>0)?_fe[ideg]->_NoShape[ _nFSIdim-2]:0);                  //     ...
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  }
#if FSI_EQUATIONS%2==0
  el_ndof[1]= _fe[1]->_NoShape[ _nFSIdim-1];
#endif


  el_mat_nrows +=  el_ndof[0]*_nvars[0];
  int el_mat_ncols = el_mat_nrows;                                //square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);                  // element dof vector

  // fields -> Navier-Stokes ----------------------------------------------------------------------

//   int _FSI_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[FSI_F]];      // FSI equation [FSI_F]]
  double u_nlg[DIMENSION];
  double u_old[DIMENSION*NDOF_FEM];
  double u_oold[DIMENSION*NDOF_FEM];
  double u_nl[DIMENSION*NDOF_FEM];       /* velocity vector for non linear terms -> it contains all the velocity components */
  double p_proj[DIMENSION*NDOF_FEM];
  double dp_proj[DIMENSION*NDOF_FEM];
  double ds_old[DIMENSION*NDOF_FEM];
  double val_tbg[30];
//   double h_eff;               // turbulence, h_eff
//   double Pe_h;
//   double f_upwind,Phi_supg;            // local Peclet, upwind
  double src_value[DIMENSION];
  double bc_value[DIMENSION];
  double x_m[DIMENSION];


//   int T_idx = _data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];       // Equazione temperatura
// coupling  basic system fields -------------------------------------------------------------------------------
  int idx_t= _data_eq[2].tab_eqs[T_F]; // Temperature ------------------------------------
  _T_idx=(idx_t>=0)?_data_eq[2].indx_ub[idx_t]:-1;
  int idx_ns= _data_eq[2].tab_eqs[FS_F];// Navier-Stokes ---------------------------------
  _FSI_idx=(idx_ns>=0)? _data_eq[2].indx_ub[idx_ns]:-1;
  double  val_ve[DIMENSION];   double  val_vg[DIMENSION];
   int idx_ds= _data_eq[2].tab_eqs[SDSX_F];// Displacement ---------------------------------
  _DS_idx=(idx_ds>=0)? _data_eq[2].indx_ub[idx_ds]:-1;
  for(int idim=0; idim<DIMENSION; idim++) {val_vg[idim] =0.; val_ve[idim] =0.; }
  int idx_k= _data_eq[2].tab_eqs[K_F];// Turbulence --------------------------------------
  _K_idx=(idx_k>=0)? _data_eq[2].indx_ub[idx_k]:-1;
  double kappa_mg[2];
  for(int ieq=0; ieq<2; ieq++) { kappa_mg[ieq] =0.; }

  /* ============== Projection and Penalty ============== */
  int  el_ndofp=  (el_ndof[0]>0)? _fe[0]->_NoShape[ _nFSIdim-1]   :_fe[1]->_NoShape[ _nFSIdim-1];
#if FSI_EQUATIONS==1
  double penalty_f2=0.e+0;
  double impl_pen_f=0.e+0;//implict term of penalty method fluid
#else
  double proj_f=1.e+0;
  double penalty_f2=0.e0;
  double impl_pen_f=1.e+0;//implict term of penalty method fluid
//   int  el_ndofp=_fe[1]->_NoShape[ _nFSIdim-1];
#endif
  double val_matp[27*DIMENSION];
//   int imesh=atoi(_mgutils.get_file("MESHNUMBER").c_str());

  MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);

#ifdef  HAVE_MED
    EquationSystemsExtendedM *ext_es=dynamic_cast<EquationSystemsExtendedM *>(&_mgeqnmap);
    int interface_id2=3;    double src_value2[1];
    InterfaceFunctionM * source2= ext_es->get_interface_fun(interface_id2);
    int * map_mg2 ;
    int * map_med2;
    int n_map2;
    map_mg2  = source2->get_map_mg();
    map_med2 = source2->get_map_med();
    n_map2 = source2->get_n();

// 	double cc1=0;
// 	    for (int iz=0; iz<n_map; iz++) {
//
// 	 source->eval( iz, 1, src_value);
//                    cc1 += src_value[0];
// 		   std::cout <<  iz << "  "<< src_value[0]<< "  ";
// 	    }
// 	std::cout <<  "cc1"<< cc1 << "\n ";

#endif

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
  A[Level]->zero();  if(mode ==1) { b[Level]->zero(); }                 // global matrix+rhs
  DenseMatrixM KeM;  DenseVectorM FeM;                     // local  matrix+rhs

  KeM.resize(el_mat_nrows,el_mat_ncols);  FeM.resize(el_mat_nrows); // resize  local  matrix+rhs

  int ndof_lev=0;
  for(int pr=0; pr <_mgmesh._iproc; pr++) {
    int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    ndof_lev +=delta;
  }
 int flag_interface=0;
  /*===================================================================================================*/
  /*                                    B) Element  Loop over the volume (n_elem)                      */
  /*===================================================================================================*/
//
  const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    // set to zero matrix and rhs and center
    KeM.zero();        FeM.zero();

    // geometry and element  fields ------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,_xx_qnds);
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);
    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
    // field data  ------------------------------------------------------
    get_el_field_data(iel, Level,
                      el_conn,
                      offset,el_ndof,ndof_lev,
                      u_old,u_oold, u_nl,p_proj,dp_proj);


    //  Average vel and init flag bc , ---------------------------------------
    for(int idim=0; idim< _nFSIdim; idim++) {// based on scalar product between v and each element side
      u_nlg[idim] =0.; for(int d=0; d< NDOF_FEM; d++) {
        const int  dnode=idim*NDOF_FEM+d;    // index local points
        u_nlg[idim] += u_nl[dnode]/NDOF_FEM; // non linear solution average vel
        _bc_el[dnode]=1;                     // Neumann flag to all points
      }
    }
     flag_interface=0;
     for(int d=0; d< NDOF_FEM; d++) {
       _bc_el[_nFSIdim*NDOF_FEM+d]=1;
       if(ext_mesh->_bc_id[el_conn[d]] >999) flag_interface=1;
     }
       
    //  h_eff computation for SUPG, ---------------------------------------
   double  h_eff=1.e-20; double vdothmax_old=1.e-20; for(int kdim=0; kdim< _nFSIdim; kdim++) { x_m[kdim]=0.; }
    for(int d=0; d< NDOF_P; d++) {
      double hcurr=0, vdothmax=0.;
      for(int idim=0; idim<_nFSIdim; idim++) {
        x_m[idim] +=_xx_qnds[idim*NDOF_FEM+d]/NDOF_P;
        const double dist_idim=(_xx_qnds[idim*NDOF_FEM+(d+1)%NDOF_P]-_xx_qnds[idim*NDOF_FEM+d]);
        hcurr += dist_idim*dist_idim;  vdothmax += abs(u_nlg[idim]*dist_idim); // ????????????????????
      }
      if(vdothmax > vdothmax_old) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;}
      else if(vdothmax > vdothmax_old-1.e-10 && h_eff>sqrt(hcurr)) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;}
    }
    // end h_eff computation--------------------------------------------------
    // // // //     VOF
            double xm[3];
        xm[2]=0.;
        for (int idim=0; idim<DIMENSION; idim++) xm[idim]=0.;
        for (int idim=0; idim<DIMENSION; idim++) {
            for (int d=0; d< NDOF_FEM; d++) xm[idim] += _xx_qnds[idim*NDOF_FEM+d]/NDOF_FEM;
        }
        double  hx[3];
        hx[2]=0;
        int ixyz[3];
        ixyz[2]=0;
        int nx[3]= {NX,NY,NZ};
        for (int idim=0; idim<_nFSIdim; idim++) {
//            nx[idim]=16;
            hx[idim] = 1./nx[idim];
            ixyz[idim]=(int)((xm[idim]-0.25*hx[idim])*nx[idim]);
        }

        double cc=0; double cc2=0;
        int count=0;
        for (int iz=0; iz<_nFSIdim-1; iz++) {
            for (int ix=0; ix<2; ix++) {
                for (int iy=0; iy<2; iy++) {
                    int indx = ixyz[0]+ix+( nx[0]+1)*(ixyz[1]+iy)+( nx[0]+1)*( nx[1]+1)*(ixyz[2]+iz);
//                     source->eval(indx, 1, src_value);
//                     cc += src_value[0];
		    source2->eval(indx, 1, src_value2);
                    cc2 += src_value2[0];
// 		     if(fabs(src_value[0])>0.000001){std::cout << " cccccccccccccccccccccccccc"<<src_value[0];}
                    count++;
                }
            }
        }
// // //         END VOF
    // element fields ----------------------------------
          _vof_phase= cc2/4;
    
    double phase=1.; double rho =phase;// _y_bcout=0.001;
    if(_K_idx>=0) {
#ifdef DIST_FIX
      _y_bcout=DIST_FIX;
//       _y_bcout=(x_m[0]>0.15)? 0.3-x_m[0] : x_m[0] ; //distance from wall [0,0.3]
 #else
     _y_bcout=_mgmesh._dist[ iel+nel_b];
#endif
    }
    kappa_mg[0] =0; kappa_mg[1] =0;  mu_m =0.;

    // MED functions ----------------------------------------------------------------
    for(int i=0 ; i< _nFSIdim; i++) { src_value[i]=0.; }
// #ifdef  HAVE_MED
//     if(imesh==1) {
//       if(mat[iel+nel_b]==interface_id) {
//         int i_mg =0;  int i_med=-1;
//         for(int i_mg =0; i_mg<n_map; i_mg++) if(map_mg[i_mg] ==  iel+nel_b) {i_med=i_mg; break;}
//         if(i_med >-1) {
//           int elem_med= map_med[i_med];
//           source->eval(elem_med,  _nFSIdim, src_value);
//         }
//       }
//     }
// #endif
    // ------------------ Calculation of the cell velocity divergence 'sum' ----------------------------------
//     double sum=0.;
// #if FSI_EQUATIONS%2==0
//     for(int qp=0; qp<  el_ngauss; qp++) {
//       // shape functions at gaussian points --------------------------------------------------------------
//       // quadratic continuous
//       det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);     // Jacobian
//       JxW_g2 =det2*_fe[2]->_weight1[ _nFSIdim-1][qp];       // weight
//       _fe[2]->get_phi_gl_g(_nFSIdim,qp,_phi_g[2]);
//       _fe[2]->get_dphi_gl_g(_nFSIdim,qp,_InvJac2,_dphi_g[2]);  // global coord deriv
// 
//       det1      = _fe[1]->Jac(qp,_xx_qnds,_InvJac1);     // Jacobian
//       _fe[1]->get_phi_gl_g(_nFSIdim,qp,_phi_g[1]);
//       _fe[1]->get_dphi_gl_g(_nFSIdim,qp,_InvJac1,_dphi_g[1]);  // global coord deriv
//       // Assemblying Matrix quad ---------------------------------
// 
//       for(int jvar=0; jvar<DIMENSION; jvar++) {
//         for(int j=0; j<el_ndof[2]; j++) {
//           double vv=JxW_g2/**_bc_vol[j+jvar*NDOF_FEM]*/*_dphi_g[2][j+jvar*el_ndof[2]];
//           sum += vv*u_oold[j+jvar*NDOF_FEM];
//         }
//       }
// 
//     }
// //          std::cout<<"iel "<< iel<<" sum "<< sum <<" sum2 "<< sum2 <<std::endl;
// #endif

// ----------------------------------------------------------------------------------
//                                       Boundary
// ----------------------------------------------------------------------------------
//     double pressure[NDOF_FEMB];


    for(int  iside=0; iside< el_sides; iside++) {
      if(el_neigh[iside] == -1) {
        // setup boundary element -> connectivity+coordinates
        for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
          sur_toply[lbnode]=lnode;          // lbnode -> lnode
          elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
          for(int idim=0; idim< _nFSIdim; idim++) { // coordinates
            _xxb_qnds[idim*NDOF_FEMB+lbnode]=_xx_qnds[idim*NDOF_FEM+lnode];
          }
//           pressure[lbnode]=_data_eq[1].ub[lnode];
        }
        // u turbulent
//         double u_tau[elb_ndof[2]];    double vel_bound[1];  double y1_plus[elb_ndof[2]];
//         for(int i=0; i<elb_ndof[2]; i++) {
//           vel_bound[0] = fabs(_data_eq[2].ub[(_FSI_idx+AX_DIR)*NDOF_FEM+sur_toply[i]]);
//           u_tau[i]=(_K_idx == -1)?  0.: eval_var2(vel_bound);
//           y1_plus[i]   = _y_bcout*u_tau[i]/_IRe;
//         }

        // det hyperplane
        det2= _fe[2]->JacSur(elb_ngauss-1,_xxb_qnds, _InvJac2);// jacobian
        double Ipenalty=det2/_dt;                            // Dirichlet bc flag
        // normal
        _fe[2]->normal_g(_xxb_qnds,x_m,normal);
        int dir_maxnormal = (fabs(normal[0])>fabs(normal[1]))?0:1 ;
        dir_maxnormal= (fabs(normal[dir_maxnormal])>fabs(normal[DIMENSION-1]))? dir_maxnormal:DIMENSION-1;
//         if(x_m[0]< -0.139 && x_m[1]< -0.19 ){
//          double s=0.; 
//         }
        //   normal boundary condition (Diriclet normal=0; Neumann=1) ----------------
        // enum bound_cond{simm=0,velocity_in0=1,velocity_tg0=2,wall=3,simm1=4,velocity_in=5,velocity_tg=6,velocity=7,
        //                 pressure_outlet0=10,interior=11,pressure_inlet=14)
        set_bc_matrix(KeM,FeM,dir_maxnormal,sur_toply,
                      el_ndof,elb_ndof, elb_ngauss,
                      u_old, normal,Ipenalty);
//            cout << "KeM" << endl; cout << KeM << endl;cout << "FeM" << endl; cout << FeM << endl;// check
      } // iside
    }  //
 
  
    // ------------------------------ Volume --------------------------------------------

  if(ext_mesh->_mat_id[ iel+nel_b]  !=4){  // liquid region
    matrixrhsvol_liq(KeM,FeM,
                 el_ndof,mode,FSI_parameter.CRANK_NICK,h_eff,
                 u_old,u_oold, u_nl,
                 p_proj,dp_proj,flag_interface);
    
  }
    else{  // solid region ==> 4
     
        matrixrhsvol_sol(KeM,FeM,
                 el_ndof,mode,FSI_parameter.CRANK_NICK,h_eff,
                 u_old,u_oold, u_nl,
                 p_proj,dp_proj,flag_interface);
    }

    // ---------------------- end volume (element) --------------------------------------

    // ----------------------------------------------------------------------------------
    //   E) Add them to the global matrix
    // ----------------------------------------------------------------------------------
//           cout << "KeM" << endl; cout << KeM << endl;cout << "FeM" << endl; cout << FeM << endl; // check 
    A[Level]->add_matrix(KeM,el_dof_indices);
    if(mode == 1) { b[Level]->add_vector(FeM,el_dof_indices); }

  } //  =============== End of element loop =============================================

// ===========================================================================================================

  el_dof_indices.clear();
  A[Level]->close();
  if(mode == 1) { b[Level]->close(); }
//   A[Level]->print();
//   b[Level]->print();
// ----------------------------------------------------------
#ifdef PRINT_INFO
  std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

  return;
} /******************************************************************************************************/
//

// void my1_add(DenseMatrixM & KeM) {
//
// //     std::vector<int> el_dof_indices(1);
//   for(int j=0; j<9; j++)  {KeM(j,j) +=1.e-12;}
// //     el_dof_indices[0]=0;
// //     Ke
// //     A[0]->add_matrix(KeM,el_dof_indices);
// //     el_dof_indices.clear();
//   return;
// }

//
// void MGSolFSI::my_add() {
//   DenseMatrixM KeM;
//   std::vector<int> el_dof_indices(1);
//   el_dof_indices[0]=0;
//   KeM.resize(1,1);
//   A[0]->add_matrix(KeM,el_dof_indices);
//   el_dof_indices.clear();
// }

// =========================================================================================
/// This function controls the assembly and the solution of the FSI_equation system:
void MGSolFSI::MGTimeStep(
  const double time,  // time
  const int iter  // Number of max inter
) {
// =========================================================================================

  /* ========================================================================================= */
  /*              A) Set up the time step                                                      */
  /* ========================================================================================= */
//
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
  x[_NoLevels-1]->localize(*x_nonl[_NoLevels-1]);


  /* ========================================================================================= */
  /*              B) Assemblying of the Matrix-Rhs                                             */
  /* ========================================================================================= */

#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                              // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); }  // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
//
  /* ========================================================================================= */
  /*              C) Solution of the linear MGsystem (MGSolFSI::MGSolve)                       */
  /* ========================================================================================= */
//
  MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
//   x_oold[_NoLevels-1]->close();
//   x_oold[_NoLevels-1]->zero();
//   double norm_check=  x_oold[_NoLevels-1]->l2_norm();
//   if(norm_check>2.e+5) { x_oold[_NoLevels-1]->scale(1.e-3); }
//   double norm_diff=1.;
//
  /* ========================================================================================= */
  /*              D) Non Linear Iterations                                                     */
  /* ========================================================================================= */
  
  int iter_nl=FSI_parameter.NL_ITER;if(time < FSI_parameter.NL_TIME0) iter_nl=FSI_parameter.NL_ITER0;
  double penalty =1.e+1; double factor =1.e+0;
  for(int iter=0; iter<iter_nl; iter++) {
    x[_NoLevels-1]->localize(*x_nonl[_NoLevels-1]);
//     double norm_in= x_nonl[_NoLevels-1]->l2_norm();
//     std::cout <<" norm xnl*********** " << norm_in  << std::endl;
// #if FSI_EQUATIONS%2==0
//     x_oold[_NoLevels-1]->add(penalty,*x_nonl[_NoLevels-1]); // penalty
// //     penalty*=factor;
// #endif
//     std::cout <<" norm x_oold--------- " <<  x_oold[_NoLevels-1]->l2_norm()  << std::endl;
//
    /* ----------- Non linear solution ------------------------------------------------------- */
//
    std::cout  << std::endl << " === NON LINEAR ITERATOR: "<< iter+1 << "- " << _eqname.c_str() << " solution "  << std::endl;
//     A[_NoLevels-1]->close();
    GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
    for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); } // matrix
    MGSolve(1.e-6,40);                                                            // solve
//
    /* ----------- Check error --------------------------------------------------------------- */
//
//     x[_NoLevels-1]->localize(*disp[_NoLevels-1]);
//     disp[_NoLevels-1]->add(-1.e+0,*x_nonl[_NoLevels-1]);
//     norm_diff=  disp[_NoLevels-1]->l2_norm();
//     std::cout << " Check:  "<< " norm_in "<< norm_in <<
//               "  norm_diff/-norm_in " <<norm_diff/norm_in <<std::endl
//               <<" === END NON LINEAR ITERATOR"<<std::endl;
//     if(norm_diff/norm_in<1.e-6) { break; }

  }
  x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);

  return;
} /******************************************************************************************************/
//
//



// ================================================================
// ====================== END velocity ============================

// ================================================================



// ================================================================



#endif  //ENDIF FSI_EQUATIONS
// #endif  // FSI_equation is personal

