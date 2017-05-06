// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS
// #if NS_EQUATIONS==1
// ======================================================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file


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

/// \ingroup  user_function ic_read
// ======================================================================================
/// This function generates the initial conditions for the NS system:
void MGSolNS::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]
) {// ===================================================================================
  

   double ILref = 1./_lref;   
#if DIMENSION==2
  
  #if NS_EQUATIONS==2
 
    //    boundary conditions box file plane_ch.med ==================================
    if (_dir==0)  u_value[0] =0.;    if (_dir==1)  u_value[0] = 2.; 
    // ================================================================================

  #else
    //    boundary conditions box file plane_ch.med ==================================
    u_value[0] = 0.;  u_value[1] = 2.;
   // ================================================================================ 
    
  #endif
  
  #if NS_EQUATIONS==1 //coupled
    double pref = _refvalue[DIMENSION];
    u_value[DIMENSION] =10.;// (LYE*ILref - xp[1])/pref;
  #endif

#endif

#if DIMENSION==3


  #if NS_EQUATIONS==1
  u_value[0]=0.;  u_value[1]=0.;  u_value[2]=1.;  u_value[3]=0.;
  #endif//NSEQUATION==1

  #if NS_EQUATIONS==2
  if (_dir==0) u_value[0]=0.;//x direction
  if (_dir==1) u_value[0]=0.;
  if (_dir==2) u_value[0]=0.17;
  #endif//NSEQUATION==2

#endif
  return;
}

/// \ingroup  user_function bc_read
// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolNS::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], 	// normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0  ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 ->  tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;

#if DIMENSION==2  // ----- 2D boundary conditions ----------

  bc_flag[0]=0;
  //    boundary conditions box file plane_ch.med ==================================
  bc_Neum[0]=wall;
  if(bc_gam==11) {bc_Neum[0]=pressure_inlet;}                          // top
  if(bc_gam==12) {bc_Neum[0]=wall/*simmx*/;}                           // left side
  if(bc_gam==13) {bc_Neum[0]=wall;}                                    // right side
  if(bc_gam==10) {bc_Neum[0]=velocity_in;}                             // bottom
//   if(xp[0]<0.00001 && xp[1]<0.00001) {bc_Neum[0]=wall/*velocity_in*/;} // pt bottom
  
// ======================= end  plane_ch.med ======================================
  
  
  

#endif // //---------------------------------------------

#if DIMENSION==3 // -------- 3D boundary conditions ---------



#endif  //  end DIMENSION ==3 ---------------------------

  return;
}



// double MGSolNS::musker (double dist, double utau, double nu){
//     double yplus = dist*utau/nu;
//     double vel = 5.424*atan((2.*yplus-8.15)/16.7) 
// 	       + 4.1693*log(yplus + 10.6) 
// 	       - 0.8686*log(yplus*yplus - 8.15*yplus+86)
// 	       - 3.52;
//     return vel; 
//    }
// #endif  //
#endif  //ENDIF NS_EQUATIONS
