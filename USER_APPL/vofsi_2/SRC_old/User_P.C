// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ======================================================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
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
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#if (NS_EQUATIONS%2==0)


// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolNSP::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]) {
// xp[]=(xp,yp) u_value[]=(u,v,p)
   u_value[0] =(0.9-xp[0])*(0.9-xp[1])*1000.;
//    u_value[0] =1.e+0;
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolNSP::bc_read(
    int bc_gam,
    int bc_mat,
  double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================
  //  enum bound_cond_p { outflowp0=0,outflowp=4, vel_fix=10,interiorp=11}
 // =================================================
  
//   const double Lref = _mgphys.get_par("Lref");
  double ILref = 1.;
//
#if DIMENSION==2
      
     
  
//   if ( bc_gam == 11 ) {
//     bc_Neum[0] =vel_fix; 
//     
//   } 
//   if ( bc_gam == 13 ) { 
//       bc_Neum[0] =outflowp0;
//     }
//      if ( bc_gam == 14 ) {
//      bc_Neum[0] = vel_fix;
//      
//   } 
//   if ( bc_gam == 12 ) {
//       bc_Neum[0] =vel_fix;  
//       
//     } 
  
   
  
  
//   
//     if(xp[0]< 0.00001 && xp[1] < 0.00001) { bc_Neum[0] =vel_fix;} // left bottom corner
//         if(xp[0]>1.99 && xp[1] < 0.01) { bc_Neum[0] =velocity;} // right bottom corner
//          if(xp[0]>1.9999 && xp[1] > 0.9999) { bc_Neum[0] =wall;} // right top corner
//     if(xp[0]< 1.0001 && xp[1] > 0.999) { bc_Neum[0] =wall;} // left top corner
//   if(xp[0]< 0.00001 && xp[1] < 0.00001) { bc_Neum[0] =wall;} // left bottom corner
     // =======================================================================
//  

 bc_Neum[0] =vel_fix;
 if (bc_gam == 14) {bc_Neum[0]=outflowp0;}    
  if (xp[0]>0.9-BDRY_TOLL) {bc_Neum[0]=outflowp0;}  // top

#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
//   bc_Neum[0]=1;    bc_flag[0]=0;
//   if (xp[2]> LZE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }
//   if (xp[2]< LZB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // left
//   if (xp[0]> LXE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // right
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }  // bottom
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }  // top
//
//
//     /* Spigolo 1 */
//     if (xp[1] > LYE*ILref - BDRY_TOLL && xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }
//     /* Spigolo 2 */
//     if (xp[1] < LYB*ILref + BDRY_TOLL && xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }
//     /* Spigolo 3 */
//     if (xp[1] < LYB*ILref + BDRY_TOLL && xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }
//     /* Spigolo 4 */
//     if (xp[1] > LYE*ILref - BDRY_TOLL && xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }
//
//     /* Spigolo 5 */
//     if (xp[0] > LXE*ILref - BDRY_TOLL && xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }
//     /* Spigolo 6 */
//     if (xp[0] > LXE*ILref - BDRY_TOLL && xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }
//     /* Spigolo 7 */
//     if (xp[0] > LXE*ILref - BDRY_TOLL && xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }
//     /* Spigolo 8 */
//     if (xp[0] > LXE*ILref - BDRY_TOLL && xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }

//     bc_Neum[0]=0; bc_flag[0]=0;
//    if( bc_gam==15 || bc_gam==16){
//      bc_Neum[0]=1; bc_flag[0]=0;
//     }    // outlet
//    if( bc_gam==50){
//      bc_Neum[0]=0; bc_flag[0]=0;
//     }
//   if( bc_gam==100 || bc_gam==101){
//      bc_Neum[0]=1; bc_flag[0]=0;
//     }

#endif
  return;
} // end boundary conditions ==========================



#endif // ENDIF NS_EQUATIONS==0

// double MGSolNS::musker (double dist, double utau, double nu){
//     double yplus = dist*utau/nu;
//     double vel = 5.424*atan((2.*yplus-8.15)/16.7) 
// 	       + 4.1693*log(yplus + 10.6) 
// 	       - 0.8686*log(yplus*yplus - 8.15*yplus+86)
// 	       - 3.52;
//     return vel; 
//    }
#endif  //ENDIF NS_EQUATIONS
