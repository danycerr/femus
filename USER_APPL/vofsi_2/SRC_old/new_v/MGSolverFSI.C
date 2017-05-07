// ===============================================================
// --------------   NAVIER-STOKES system [FS_F] ------------------
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



// ==================================================================
/// This routine constructs the FSI class:
MGSolFSI::MGSolFSI(
  MGEquationsSystem &mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
//
//===================================================================================================//
//                                    A) Reading parameters                                          //
//===================================================================================================//
//
//   _offset(_mgmesh._NoNodes[_NoLevels-1]),             // mesh nodes (top level)
//   _dt(_mgutils.get_par("dt")),                        // parameter  dt
//   _uref(mg_equations_map_in.get_par("Uref")),         // parameter  u reference
//   _lref(mg_equations_map_in.get_par("Lref")),         // parameter  l reference
//   _rhof(mg_equations_map_in.get_par("rho0")),         // parameter density
//   _muf(mg_equations_map_in.get_par("mu0")) {          // parameter viscosity
    
  _offset(_mgmesh._NoNodes[_NoLevels-1]),             // mesh nodes (top level)
  _dt(_mgutils.get_par("dt")),                        // parameter  dt
  _uref(mg_equations_map_in.get_par("Uref")),         // parameter  u reference
  _lref(mg_equations_map_in.get_par("Lref")),         // parameter  l reference
  _rhof(mg_equations_map_in.get_par("rho0")),         // parameter density
  _muf(mg_equations_map_in.get_par("mu0")),
   _rhos(mg_equations_map_in.get_par("rhos")),    // parameter density
  _ni(mg_equations_map_in.get_par("nis")),
  _Emod(mg_equations_map_in.get_par("Es")),
  _hs(mg_equations_map_in.get_par("hs"))    {
    
  //===================================================================================================//
  //                                    B) Setting class variables                                     //
  //===================================================================================================//

  _nNSdim=DIMENSION;
  // class equation ---------------------------------------------------------------------
  for(int k_index=0; k_index<30; k_index++) { _FF_idx[k_index]=-1; }
  _dir=0;
  _pres_order=(_nvars[0]>0)? 0:1;
  // class variable names ------------------------------------------------------------
#if FSI_EQUATIONS==2       //   segregated ( P in NSP_EQUATIONS)
  if(!varname_in.compare("u")) { _dir=0; }   // u-equation
  if(!varname_in.compare("v")) { _dir=1; }   // v-equation
  if(!varname_in.compare("w")) { _dir=2; }   // w-equation
  _var_names[0]=varname_in;   _refvalue[0]=_uref;
#else
  _var_names[_nNSdim-1]="w"; _refvalue[_nNSdim-1]=_uref; // velocity 3D
  _var_names[0]="u";   _refvalue[0]=_uref; // velocity 2D
  _var_names[1]="v";   _refvalue[1]=_uref; // velocity 2D
#if FSI_EQUATIONS==1       //  coupled  (+P)
  _var_names[_nNSdim]="p";
  _refvalue[_nNSdim]=_rhof*_uref*_uref;  // pressure
#endif
#endif

  //===================================================================================================//
  //                                    C) Setting solver type                                         //
  //===================================================================================================//

  for(int l=0; l<_NoLevels; l++) { _solver[l]->set_solver_type(SOLVER_FSI); }

  //===================================================================================================//
  //                                    D) Setting non_dimensional parameters                          //
  //===================================================================================================//

  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
  _IFr=9.81*_lref/(_uref*_uref);          // Froud number
  _dirg[0] = mg_equations_map_in.get_par("dirgx");      // x-gravity
  _dirg[1] = mg_equations_map_in.get_par("dirgy");      // y-gravity
  _dirg[2] = mg_equations_map_in.get_par("dirgz");      // z-gravity

  _y_bcout=sqrt(2*0.09*_IRe/80.);
  _lambda=(_ni*_Emod)/(_rhos*(1+_ni)*(1-2*_ni));//.001;
  _mus=_Emod/(_rhos*2*(1+_ni)*_lref*_uref);
  _control;

  return;
} //****************************************************************************************************//
//
//

// #include "MGsolverNS_vol.h"

// ====================================================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================================================
void  MGSolFSI::GenMatRhs(const double/* time*/, const int
                         Level,const  int mode) {
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// =====================================================================================================

  // ===================================================================================================
  //                                    A) Set up
  // ===================================================================================================

  // NS parameters
  const int unsteady_flag=1;//_FSI_parameter.UNSTEADY;
  const int iaxisim=_FSI_parameter.AXISYM;
  const double les=_FSI_parameter.LES;
  // -------------------------- Geometry -------------------------------------------------------------
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     normal[DIMENSION];                                          // normal to the boundary

  // -------------------------- Gauss integration -----------------------------------------------------

  const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] ---------------------------------------
  int el_ndof[3];  el_ndof[0]=NDOF_K;  int elb_ndof[3];
  if(el_ndof[0]>0) elb_ndof[0]=1; // number of el dofs
  int el_mat_nrows =0;                                            // number of mat rows (dofs)
  for(int ideg=1; ideg<3; ideg++) {                               //     ...
    el_ndof[ideg]=((_nvars[ideg]>0)?  _fe[ideg]->_NoShape[ _nNSdim-1]:0);                    //   computing
    elb_ndof[ideg]=((_nvars[ideg]>0)?_fe[ideg]->_NoShape[ _nNSdim-2]:0);                  //     ...
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  }
#if FSI_EQUATIONS%2==0
  el_ndof[1]= _fe[1]->_NoShape[ _nNSdim-1];
#endif
  el_mat_nrows +=  el_ndof[0]*_nvars[0];
  int el_mat_ncols = el_mat_nrows;                                //square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);                  // element dof vector

  // fields -> Navier-Stokes ----------------------------------------------------------------------
  double u_nlg[DIMENSION];
  double u_old[DIMENSION*NDOF_FEM];
  double u_oold[DIMENSION*NDOF_FEM];
  double u_nl[DIMENSION*NDOF_FEM];       // velocity vector for non linear terms -> it contains all the velocity components //
  double p_proj[DIMENSION*NDOF_FEM];
  double dp_proj[DIMENSION*NDOF_FEM];
//   double x_m[DIMENSION];
  int flag_group[NDOF_FEM];
 MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
// coupling  basic system fields -------------------------------------------------------------------------------



  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
  A[Level]->zero();  if(mode ==1) { b[Level]->zero(); }             // global matrix+rhs
  DenseMatrixM KeM;  DenseVectorM FeM;                              // local  matrix+rhs
  KeM.resize(el_mat_nrows,el_mat_ncols);  FeM.resize(el_mat_nrows); // resize  local  matrix+rhs

  int ndof_lev=0;
  for(int pr=0; pr <_mgmesh._iproc; pr++) {
    int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    ndof_lev +=delta;
  }

  // ===================================================================================================
  //                                    B) Element  Loop over the volume (n_elem)
  // ===================================================================================================

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
    for(int idim=0; idim< _nNSdim; idim++) {// based on scalar product between v and each element side
      u_nlg[idim] =0.;_x_m[idim]=0.; for(int d=0; d< NDOF_FEM; d++) {
        _x_m[idim] +=_xx_qnds[idim*NDOF_FEM+d]/NDOF_FEM;
        const int  dnode=idim*NDOF_FEM+d;    // index local points
          flag_group[d]=fabs(ext_mesh->_bc_id[el_conn[d]]);
        u_nlg[idim] += u_nl[dnode]/NDOF_FEM; // non linear solution average vel
        _bc_el[dnode]=1;                     // Neumann flag to all points
      }
    }
    for(int d=0; d< NDOF_FEM; d++) {_bc_el[_nNSdim*NDOF_FEM+d]=1;}
    // set up boundary conditions (wall(5)-> all Dirichlet)
  
//  for(int d=0; d< NDOF_P; d++) {_bc_el[_nNSdim*NDOF_FEM+d]=(_bc_vol[d]/10 ==1)?0:1; }

    // element fields ----------------------------------
    //     double phase=1.; //double rho =phase;// _y_bcout=0.001;
    if(_FF_idx[K_F]>=0) { _y_bcout=_FSI_parameter.DIST_FIX;}//  _y_bcout=_mgmesh._dist[ iel+nel_b];

    // -------------------------------- Boundary ----------------------------------------------
    for(int  iside=0; iside< el_sides; iside++) {
      if(el_neigh[iside] == -1) {
        // setup boundary element -> connectivity+coordinates
        for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
          sur_toply[lbnode]=lnode;          // lbnode -> lnode
          elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
            _bc_el[lnode] *=-1;  
          for(int idim=0; idim< _nNSdim; idim++) { // coordinates
            _xxb_qnds[idim*NDOF_FEMB+lbnode]=_xx_qnds[idim*NDOF_FEM+lnode];
          }
        }
        // normal
        _fe[2]->normal_g(_xxb_qnds,_x_m,normal);
        int dir_maxnormal = (fabs(normal[0])>fabs(normal[1]))?0:1 ;
        dir_maxnormal= (fabs(normal[dir_maxnormal])>fabs(normal[DIMENSION-1]))? dir_maxnormal:DIMENSION-1;
        // boundary condition
//         if(ext_mesh->_mat_id[ iel+nel_b]  !=4){
//          set_bc_matrix(KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof, elb_ngauss,normal,iaxisim);
         
          if(ext_mesh->_mat_id[ iel+nel_b]  !=4){  // liquid region
    set_liq_bc_matrix(KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof, elb_ngauss,normal,iaxisim);
   
   }
   else{
      set_sol_bc_matrix(KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof, elb_ngauss,normal,iaxisim);
     
   }
         
         
//  if( x_m[1]>0.15 && x_m[0]<0.2){ 
//   std::cout << x_m[0] << " " << x_m[1] << " " <<  x_m[2] << " "  << "bound Kem \n "<< KeM << "  \n" <<  iel << "\n bc rhs \n" <<  FeM;
//          }
          
        } // iside
    }  // -----------------------------  End Boundary -------------------------------------

    // ------------------------------ Volume --------------------------------------------
    
   if(ext_mesh->_mat_id[ iel+nel_b]  !=4){  // liquid region
  
    matrixrhs_liq_vol(KeM,FeM,el_ndof,
                 u_old,u_oold, u_nl, p_proj,dp_proj,
                 unsteady_flag,iaxisim,les,mode,flag_group);
   }
   else{
      matrixrhs_sol_vol(KeM,FeM,el_ndof,
                 u_old,u_oold, u_nl, p_proj,dp_proj,
                 unsteady_flag,iaxisim,les,mode,flag_group);
   }
     
//      if( x_m[1]>0.15 && x_m[0]<0.2) { 
//     std::cout << x_m[0] << " " << x_m[1] << " " <<  x_m[2] << " "   << "vol KeM \n" << KeM << " ielem \n" <<  iel << " rhs \n" <<  FeM;
      
//      }
    // ---------------------- end volume (element) --------------------------------------
double max_diag=0.; 
   int count=0;
    if(ext_mesh->_mat_id[ iel+nel_b]  !=4){ 
     for(int id=0; id< DIMENSION*NDOF_FEM; id++) {
       if(_bc_el[id] !=0 && fabs(flag_group[id])<999) if(max_diag<fabs(KeM(id,id))) { max_diag=fabs(KeM(id,id)); }
     }
     for(int id=0; id< DIMENSION*NDOF_FEM; id++) if(_bc_el[id]==0) {
         const double pivot=fabs(KeM(id,id));
          if(fabs(pivot)<1.e-10) {
//           if( x_m[0]>0.1971 && x_m[0]<0.2019 && x_m[1]<-0.23){
//            double sss=1.; 
//           }
                std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr vol KeM  ielem " <<  iel << " row " << id<< " \n";
           
         }
         FeM(id) *=  max_diag/pivot;
         for(int jd=0; jd< el_mat_nrows; jd++) { KeM(id,jd) *= max_diag/pivot; }
       }
    }
    else{
        for(int id=0; id< DIMENSION*NDOF_FEM; id++) {
        if(max_diag<fabs(KeM(id,id))) { max_diag=fabs(KeM(id,id)); }
     }
     for(int id=0; id<  DIMENSION*NDOF_FEM; id++)  {
         const double pivot=fabs(KeM(id,id));
          if(fabs(pivot)<1.e-10) {
//           if( x_m[0]>0.1971 && x_m[0]<0.2019 && x_m[1]<-0.23){
//            double sss=1.; 
//           }
                std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr vol KeM  ielem " <<  iel << " row " << id<< " \n";
           
         }
         FeM(id) *=  max_diag/pivot;
         for(int jd=0; jd< el_mat_nrows; jd++) { KeM(id,jd) *= max_diag/pivot; }
       }
      
      
    }
//    if( x_m[0]>0.48)  
//  if( x_m[0]>0.1971 && x_m[0]<0.2019&& x_m[1]<-0.23 ) 
//       std::cout << " renorm KeM \n" << KeM << " ielem \n" <<  iel << " rhs \n" <<  FeM;
    // ----------------------------------------------------------------------------------
    //   E) Add them to the global matrix
    // ----------------------------------------------------------------------------------
    //   check -----------------------------------------
    A[Level]->add_matrix(KeM,el_dof_indices);
    if(mode == 1) { b[Level]->add_vector(FeM,el_dof_indices); }

  } //  =============== End of element loop =============================================

// ===========================================================================================================

  el_dof_indices.clear();
  A[Level]->close();  if(mode == 1) { b[Level]->close(); }
//    A[Level]->print();  b[Level]->print();
// ----------------------------------------------------------
#ifdef PRINT_INFO
  std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

  return;
} //****************************************************************************************************//
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
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolFSI::MGTimeStep(
  const double time,  // time
  const int iter  // Number of max inter
) {
// =========================================================================================

  // ========================================================================================= //
  //              A) Set up the time step                                                      //
  // ========================================================================================= //
//
  
    for(int k=0; k<30; k++) {// coupling  basic system fields
    const int idx= _data_eq[2].tab_eqs[k];
    _FF_idx[k]=(idx>=0)?_data_eq[2].indx_ub[idx]:-1;
  }
  for(int kdim=0;kdim<_nNSdim;kdim++) {
     const int num=_data_eq[2].tab_eqs[SDSX_F+kdim];
    _mgmesh.MoveMesh(kdim,(*_data_eq[2].mg_eqs[num]->x_old[_NoLevels-1]));
  }
  
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
  x[_NoLevels-1]->localize(*x_nonl[_NoLevels-1]);

// _mgmesh.MoveMesh(_NoLevels-1,_dir,*x_old[_NoLevels-1]);
  // ========================================================================================= //
  //              B) Assemblying of the Matrix-Rhs                                             //
  // ========================================================================================= //

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
  // ========================================================================================= //
  //              C) Solution of the linear MGsystem (MGSolFSI::MGSolve)                       //
  // ========================================================================================= //
//
  if(_mgutils.get_name() != 1){
  MGSolve(1.e-6,40);
  }
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
  // ========================================================================================= //
  //              D) Non Linear Iterations                                                     //
  // ========================================================================================= //

  int iter_nl=_FSI_parameter.NL_ITER; if(time < _FSI_parameter.NL_TIME0) { iter_nl=_FSI_parameter.NL_ITER0; }
//   double penalty =1.e+1; double factor =1.e+0;
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
    // ----------- Non linear solution ------------------------------------------------------- //
//
    std::cout  << std::endl << " === NON LINEAR ITERATOR: "<< iter+1 << "- " << _eqname.c_str() << " solution "  << std::endl;
//     A[_NoLevels-1]->close();
    GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
    for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); } // matrix
    MGSolve(1.e-6,40);                                                            // solve
//
    // ----------- Check error --------------------------------------------------------------- //
//
//     x[_NoLevels-1]->localize(*disp[_NoLevels-1]);
//     disp[_NoLevels-1]->add(-1.e+0,*x_nonl[_NoLevels-1]);
//     norm_diff=  disp[_NoLevels-1]->l2_norm();
//     std::cout << " Check:  "<< " norm_in "<< norm_in <<
//               "  norm_diff/-norm_in " <<norm_diff/norm_in <<std::endl
//               <<" === END NON LINEAR ITERATOR"<<std::endl;
//     if(norm_diff/norm_in<1.e-6) { break; }

  }

 
  
//   x_oold[_NoLevels-1]->zero();
//   x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);
   x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
// 
// 
//   const int flag_moving_mesh = _mgutils.get_par("moving_mesh");
//   const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
//   const int n_elem=_mgmesh._NoElements[0][_NoLevels-1];
//   const int offsetp=_dir*n_nodes;
// //   if(flag_moving_mesh) {
// // #ifdef FINE_MOVEMENT
//     /// E) mesh update
// //     const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
// //     const int offsetp=_dir*n_nodes;
//     for(int inode=0; inode<n_nodes; inode++) {
//       
//        for(int kdim=0; kdim<DIMENSION; kdim++) {
//       double disp=(*x_old[_NoLevels-1])(inode);//-(*x_oold[_NoLevels-1])(inode);
//       if(_mgmesh._xyzo[inode+0*kdim*n_nodes] >0.19999  /*&& kdim==0*/){
//           _mgmesh._xyz[inode+kdim*n_nodes]    =_mgmesh._xyzo[inode+kdim*n_nodes] +5.*disp;
//           _mgmesh._dxdydz[inode+kdim*n_nodes]  = _mgmesh._xyzo[inode+kdim*n_nodes]+5.*disp;
//       }
//      }
       
//     }
// #endif
// #ifdef COARSE_MOVEMENT
//     MoveMesh(_NoLevels-1);
// // disp[_NoLevels-1]->zero();
//     for(int inode=0; inode<n_nodes; inode++) {
//       _mgmesh._xyz[inode+offsetp] += (*disp[_NoLevels-1])(inode);
//       _mgmesh._dxdydz[inode+offsetp]= (*disp[_NoLevels-1])(inode);
//     }
//     MoveMesh(_NoLevels-1);
// #endif
//   }
  // ==============================================================
  
  

  return;
} //****************************************************************************************************//
//
//



// ================================================================
// ====================== END velocity ============================

// ================================================================



// ================================================================



#endif  //ENDIF FSI_EQUATIONS
// #endif  // NS_equation is personal

