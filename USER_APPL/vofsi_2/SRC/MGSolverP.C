// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
#include "MGSolverP.h"       // Navier-Stokes class header file

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

#define NON_LINEAR_ITER (0)

// ================================================================


#if (NS_EQUATIONS%2==0)
// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// ==================================================================
/// This function constructs the 3d-2D MGSolNSP class I
MGSolNSP::MGSolNSP(
  MGEquationsSystem& mg_equations_map_in,
  int             nvars_in[],    ///
  std::string     eqname_in,    ///< name equation (in)
  std::string     varname_in    ///< name variable (in)
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
//
/*===================================================================================================*/
/*                                    A) Reading parameters                                          */
/*===================================================================================================*/
//
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
  _dt(_mgutils.get_par("dt")),            // parameter  dt
  _uref(mg_equations_map_in.get_par("Uref")),         // parameter  u reference
  _lref(mg_equations_map_in.get_par("Lref")),         // parameter  l reference
  _rhof(mg_equations_map_in.get_par("rho0")),         // parameter density
  _muf(mg_equations_map_in.get_par("mu0")) {          // parameter viscosity
//
  _nPdim=DIMENSION;
  _NS_idx=-1;
  /*===================================================================================================*/
  /*                                    B) Setting class variables                                     */
  /*===================================================================================================*/
  _var_names[0]="p";  _refvalue[0]=_rhof*_uref*_uref;    // class variable names
//
//
  /*===================================================================================================*/
  /*                                    C) Setting solver type                                         */
  /*===================================================================================================*/
  for(int  l=0; l<_NoLevels; l++) { _solver[l]->set_solver_type(SOLVERNS); }
//
//
  /*===================================================================================================*/
  /*                                    D) Setting no _nNSdimensional parameters                           */
  /*===================================================================================================*/
  return;
}
//
//  ====================================================
//  This function assembles the matrix and the rhs:
//  ====================================================
//
void  MGSolNSP::GenMatRhs(
  const double /* time*/, // time  <-
  const int  Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

  /*===================================================================================================*/
  /*                                    A) Set up                                                      */
  /*===================================================================================================*/
  /*----------------------------------- Geometry  -----------------------------------------------------*/
  const int    _nNSdim = DIMENSION;                      				// dimension
  const int   offset = _mgmesh._NoNodes[_NoLevels-1];				// mesh nodes
  double      xx_qnds[DIMENSION*NDOF_FEM];           				// element node coords
  //   int         el_conn[NDOF_FEM];                    			// element connectivity
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];              	        // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];             	        // element connectivity
  int        el_neigh[NDOF_FEM];                                 	        // element connectivity
  double uvw_dxg[DIMENSION];
  /*----------------------------------- Gauss integration  --------------------------------------------*/

  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           		// Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];  double dphiidx_g[3][DIMENSION]; 		// global derivatives at g point
  const int  el_ngauss = _fe[2]->_NoGauss1[ _nNSdim-1];                		//elem gauss points
  double u_nl[DIMENSION*NDOF_FEM];  double u_old[DIMENSION*NDOF_FEM];
  // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
  int el_ndof[3];  el_ndof[0]=1;                  // number of el dofs
  int el_mat_nrows =0;                            // number of mat rows (dofs)
  for(int ideg=1; ideg<3; ideg++) {
    el_ndof[ideg]=_fe[ideg]->_NoShape[ _nNSdim-1];
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  };
  int el_mat_ncols = el_mat_nrows;                     //square matrix
  std::vector<int > el_dof_indices(el_mat_ncols);      // element dof vector
//
  /*---------------------------- Fields -> Navier-Stokes  [NS_F]  -------------------------------------*/
//
  int idx_ns= _data_eq[2].tab_eqs[NS_F];// Navier-Stokes ---------------------------------
  _NS_idx=(idx_ns>=0)? _data_eq[2].indx_ub[idx_ns]:-1;

//
  /* --------------------------- Element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --*/
//
  A[Level]->zero();                       if(mode ==1) { b[Level]->zero(); }   // global matrix+rhs
  DenseMatrixM KeM;                       DenseVectorM FeM;                // local  matrix+rhs
  KeM.resize(el_mat_nrows,el_mat_ncols);  FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs
//
  int ndof_lev=0;
  for(int pr=0; pr <_mgmesh._iproc; pr++) {
    int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    ndof_lev +=delta;
  }
//
  const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    // ===================================================================================================
    //                          B) Element  Loop over the volume (n_elem)
    // ===================================================================================================

    // ------------------------ Set to zero matrix and rhs and center -----------------------------------
    KeM.zero();         FeM.zero();
    // ------------------------ Geometry and element fields ---------------------------------------------
    // ------------------------ Element Connectivity (el_conn) and coordinates (xx_qnds) ----------------
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);
    get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

    // element nodes coordinates
    for(int  d=0; d< NDOF_FEM; d++)
      for(int  idim=0; idim<DIMENSION; idim++) {
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
      }

#if NS_EQUATIONS==2 //
    for(int idim=0; idim<_nPdim; idim++) {
      //       _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F+idim]]->get_el_nonl_sol(0,1,el_ndof[2],el_conn, offset,idim,u_nl);
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F+idim]]->get_el_sol(0,1,el_ndof[2],el_conn, offset,idim,u_old);
    }
#endif
#if NS_EQUATIONS==0
    //     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_nonl_sol(0,DIMENSION,el_ndof[2],el_conn, offset,0,u_nl);
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol(0,DIMENSION,el_ndof[2],el_conn, offset,0,u_old);
#endif
    //
    // ===================================================================================================
    //                          C) Gaussian integration loop (n_gauss)
    // ===================================================================================================

    for(int  qp=0; qp< el_ngauss; qp++) {
      // ------------------------- Shape functions at gaussian points --------------------------------------
      for(int  ideg=1; ideg<3; ideg++) {  				// linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);    	// Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ _nNSdim-1][qp];      	// weight
        _fe[ideg]->get_phi_gl_g(_nNSdim,qp,_phi_g[ideg]);                	// shape funct
        _fe[ideg]->get_dphi_gl_g(_nNSdim,qp,InvJac[ideg],_dphi_g[ideg]); 	// global coord deriv
      }
#if NS_EQUATIONS==2
      for(int kdim=0; kdim<_nPdim; kdim++) {
        interp_el_gdx(_data_eq[2].ub,_NS_idx+kdim,1,_dphi_g[2],el_ndof[2],uvw_dxg);          /* Funzione l'interpolazione delle derivate delle funzioni di forma sui nodi di gauss */
        for(int jdim=0; jdim<DIMENSION; jdim++) { _ub_dxg[jdim+kdim*DIMENSION]=uvw_dxg[jdim]; } /* _ub_dxg e' vettore che memorizza le derivate delle componenti di velocita' */
      }
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof[2],_ub_g[2]);
#else
      // interp_el_sol(_data_eq[2].ub,0,_nPdim,_phi_g[2],el_ndof[2],_ub_g[2]);
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof[2],_ub_g[2]);
      //external fields (quad only)
      interp_el_gdx(_data_eq[2].ub,_NS_idx,_nPdim,_dphi_g[2],el_ndof[2],_ub_dxg);
#endif
#ifdef AXISYM
      JxW_g2  *=_ub_g[2][0];  JxW_g1  *=_ub_g[2][0];//  printf("  %f \n",_ub_g[2][0]);
#endif

      // ===================================================================================================
      //                          D) Local (element) assemblying pressure equation
      // ===================================================================================================
      for(int  i=0; i<el_ndof[1]; i++)     {
        const double phii_g=_phi_g[1][i];
        double Div_g=0.;
        for(int  idim=0; idim<  _nPdim; idim++)  {
          dphiidx_g[1][idim]=_dphi_g[1][i+idim*el_ndof[1]];
          //        Div_g +=_ub_dxg[idim+idim*DIMENSION];      //divergence of velocity field
        }
        // =================================================
        //  enum bound_cond_p {outflowp0=0,outflowp=4,vel_fix=10,interiorp=11}
        // =================================================
        int bc_DN=_bc_vol[i]/10;
        double dtxJxW_g=JxW_g[2]*bc_DN;
        // ------------------------- Rhs Assemblying  -------------------------------------------------------
        //   if(mode == 1)  FeM(i) += -1.*dtxJxW_g*Div_g*phii_g/_dt;
        for(int  j=0; j<el_ndof[2]; j++) {
          const double phij_g=_phi_g[1][j];
          for(int  idim=0; idim<_nPdim; idim++) {
            dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
            FeM(i) += -1.*dtxJxW_g*(u_old[NDOF_FEM*idim+j])* dphijdx_g[2][idim]*phii_g/_dt;
#ifdef AXISYM
            FeM(i) += -1.*(1-(idim))*dtxJxW_g*(u_old[NDOF_FEM*idim+j])*phii_g*phij_g/_ub_g[2][0]/_dt ;
#endif
          }
        }
        // ------------------------- Matrix Assemblying  ----------------------------
        for(int  j=0; j<el_ndof[1]; j++) {
          const double phij_g=_phi_g[1][j];
          double Lap=0.;
          for(int  idim=0; idim<_nPdim; idim++) {
#ifdef  AXISYM
            Lap +=(1-(idim))*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
            dphijdx_g[1][idim]=_dphi_g[1][j+idim*el_ndof[1]];
            Lap += dphijdx_g[1][idim]*dphiidx_g[1][idim]; // Laplacian
          }
          // Pressure matrix assembling ---------------------------------------------
          KeM(i,i) += JxW_g[2]*(1-bc_DN);
          KeM(i,j) += dtxJxW_g*Lap* bc_DN;
        } // ---------------------------------------------
      }
    } // end of the quadrature point qp-loop +++++++++++++++++++++++++

    // ===================================================================================================
    //                          E) Global assemblying pressure equation
    // ===================================================================================================

    A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
    if(mode == 1) { b[Level]->add_vector(FeM,el_dof_indices); } // global rhs

  } // end of element loop
  // clean and close
  el_dof_indices.clear();
  A[Level]->close();    if(mode == 1) { b[Level]->close(); }

#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(P)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

  return;
} /******************************************************************************************************/
//
//

// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolNSP::MGTimeStep(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================
//
  /* ========================================================================================= */
  /*              A) Set up the time step                                                      */
  /* ========================================================================================= */
//
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /* ========================================================================================= */
  /*              B) Assemblying of the Matrix-Rhs                                             */
  /* ========================================================================================= */

#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); } // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
//
  /* ========================================================================================= */
  /*               C) Solution of the linear MGsystem (MGSolP::MGSolve)                        */
  /* ========================================================================================= */
//
  MGSolve(1.e-6,50);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
//
  /* ========================================================================================= */
  /*               D) Update of the old solution at the top Level		                     */
  /* ========================================================================================= */
//
  x[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);  // dp
  x_old[_NoLevels-1]->add(.99,*x_oold[_NoLevels-1]);// p

  return;
} /******************************************************************************************************/
//
//

/// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
/// =========================================================================================
void MGSolNSP::MGTimeStep_nl_setup(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================

  return;
} /******************************************************************************************************/
//
//

/// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
/// =========================================================================================

void MGSolNSP::MGTimeStep_nl_sol_up(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {

  /* ========================================================================================= */
  /*              A) Set up the time step                                                      */
  /* ========================================================================================= */
//
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

//
  /* ========================================================================================= */
  /*              B) Assemblying of the Matrix-Rhs                                             */
  /* ========================================================================================= */
//
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                               // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); } // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
//
  /* ========================================================================================= */
  /*               C) Solution of the linear MGsystem (MGSolP::MGSolve)                        */
  /* ========================================================================================= */
//
  MGSolve(1.e-6,50);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
//
  /* ========================================================================================= */
  /*               D) Update of the old solution at the top Level		                     */
  /* ========================================================================================= */
//
  x[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);  // dp
  x_old[_NoLevels-1]->add(1.,*x_oold[_NoLevels-1]);// p
  return;
}// =======================================================================================

/// ======================================================
/// This function controls the time step operations:
/// ======================================================
int MGSolNSP::MGTimeStep_nl_iter(const double time, int) {

  return 0;
} /******************************************************************************************************/
//
//


#endif  //ENDIF NS_EQUATIONS%2==0
// // *************************************************

#endif  //ENDIF NS_EQUATIONS
// #endif  // NS_equation is personal
// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 

