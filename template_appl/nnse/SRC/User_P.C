// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"
#if (NS_EQUATIONS%2==0)
// ======================================================================================
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
#include "MGSolverP.h"

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

// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================
// ==============================================================
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================



// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for 
/// the pressure equation in the pressure projection method
// =================================================
void MGSolP::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]) {
// xp[]=(xp,yp) u_value[]=(u,v,p)


// ============================================================================
//  Initial pressure conditions for mesh in MESH/test_3d.med (split unplit)
// ============================================================================
  u_value[0] =1.;
// ============================================================================
// ============================================================================
}

// ============================================================================
/// This function  defines the boundary conditions  for 
/// the pressure equation in the pressure projection method
void MGSolP::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],int bc_Neum[],int bc_flag[]) {
  // ==========================================================================
// outflowp0 = 0= Dirichlet homogeneus        (p=0)
// outflowp  = 4= Dirichlet nonhomogeneus    (p=p_0)
//
// vel_fix   =10= Neuman homogeneus or simmetry   (dp.n=0)
// interiorp =11= Neuman nonhomogeneus   (dp.n=dp0)
// ============================================================================


#if DIMENSION==2

 bc_Neum[0]=outflowp0;
// ============================================================================
//  Boundary conditions for mesh in MESH/test_ch.med (split unplit)
// ============================================================================
  if(bc_gam==12) {bc_Neum[0]=vel_fix;/*simmetry;*/}    // left
  if(bc_gam==13) {bc_Neum[0]=vel_fix;}                 // right
  if(bc_gam==10) {bc_Neum[0]=vel_fix;}                 // bottom
  if(bc_gam==11) {bc_Neum[0]=outflowp0;}                // top
// ============================================================================
// ============================================================================




#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann


// ============================================================================
//  Boundary conditions for mesh in MESH/test_3d.med (split unplit)
// ============================================================================
  bc_Neum[0]=vel_fix;
  if(bc_gam==12) {bc_Neum[0]=vel_fix;}
  if(bc_gam==13) {bc_Neum[0]=outflowp0;}    // z top
  if(bc_gam==10) {bc_Neum[0]=vel_fix;}
  if(bc_gam==11) {bc_Neum[0]=vel_fix;}
  if(bc_gam==14) {bc_Neum[0]=vel_fix;}
  if(bc_gam==15) {bc_Neum[0]=vel_fix;}
  if(xp[1]<0.001) {bc_Neum[0]=vel_fix;}
  if(xp[0]<0.001) {bc_Neum[0]=vel_fix;}
  if(xp[0]<0.25001 && xp[0]>0.24999 && xp[2]>0.4330 && xp[2]< 0.4331) {bc_Neum[0]=vel_fix;}
// ============================================================================
// ============================================================================

#endif
  return;
} // end boundary conditions ==========================



#endif // ENDIF NS_EQUATIONS==0

