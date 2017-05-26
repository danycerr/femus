

#include "Equations_conf.h"

// ============================================
#ifdef COLOR_EQUATIONS // 3D-2D Energy equation
// ============================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverCOL.h"       // Navier-Stokes class header file
#include "UserCOL.h"

// config file --------------------------------------------------------------------------
#include "MGGeomEl.h"        // Geometrical element
#include "MGFE_conf.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ======================================================================================

/**     \addtogroup user_function User functions
 * @{
 */

// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
/// \ingroup  user_function ic_read
// =================================================
void MGSolCOL::ic_read(int bc_gam,int bc_mat, double xp[],int iel,double u_value[]) {

// =================================================
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(u,v,p)
  
//    boundary conditions box plane_ch_rotate.med ==================================
  u_value[0] = 0.;
//   if(xp[1]>0.05 && xp[1]<0.15 &&xp[0]>0.15 && xp[0]<0.20) { u_value[0] = 200.; }  
// =================================================

#else
// =================================================
  // xp[]=(xp,ypzp) u_value[]=(u,v,w,p)
  u_value[0] = 0.;
#endif
}



// ========================================================
/// This function  defines the boundary conditions for the energy equation:
/// \ingroup  user_function bc_read
void MGSolCOL::bc_read(int bc_gam,int bc_mat, double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================

// enum bound_cond{
//      Twall0=0,Twall=1,
//      insulation=10,heat_flux=11,robinT=12,ext_conv=13
// };

  double ILref = 1./_lref;
#if DIMENSION==2
// // xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
  
  
//    boundary conditions box plane_ch_rotate.med ==================================

  bc_Neum[0]=Twallc;
//   if(bc_gam==11) {bc_Neum[0]=insulation;}    // top
//   if(bc_gam==12) {bc_Neum[0]=heat_flux;/*simmetry;*/}    // left
//   if(bc_gam==13) {bc_Neum[0]=heat_flux;}    // right
//   if(bc_gam==10) {bc_Neum[0]=Twall;}           // bottom
//   if(xp[0]<0.00001) {bc_Neum[0]=Twall/*Twall*/;}        // bottom
//    if(xp[1]< -0.2499) {bc_Neum[0]=Twall/*Twall*/;}
// ======================= end  plane_ch.med ======================================
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
  bc_flag[0]=insulation;
#endif

  return;
} // end boundary conditions ==========================
/** @} */

#endif
