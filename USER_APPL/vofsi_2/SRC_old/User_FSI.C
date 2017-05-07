// ======================================================================================
// ----------------------- Fuild Structure Interaction  [FSI_F] -------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef FSI_EQUATIONS
// #if FSI_EQUATIONS==1
// ======================================================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in FSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in FSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"         // Fuild Structure Interaction class conf file
#include "MGSolverFSI.h"           // Fuild Structure Interaction class header file


// config file --------------------------------------------------------------------------
#include "MGGeomEl.h"              // Geometrical element
#include "MGFE_conf.h"             // FEM approximation
#include "Printinfo_conf.h"        // Print options

// local Femus class include ------------------------------------------------------------
#include "MGMesh.h"                // Mesh class
#include "MGSystem.h"              // System class
#include "MGEquationsSystem.h"     // Equation map class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>                // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"         // algebra dense matrices
#include "sparse_matrixM.h"        // algebra sparse matrices
#include "dense_vectorM.h"         // algebra dense vectors
#include "numeric_vectorM.h"       // algebra numerical vectors
#include "linear_solverM.h"        // algebra solvers
// ======================================================================================



// ======================================================================================
/// This function generates the initial conditions for the FSI system:
void MGSolFSI::ic_read (
     int bc_gam,
     int bc_mat,
     double xp[],
     int iel,
     double u_value[]
)  // ===================================================================================
{
     double ILref = 1./_lref;
// ======================================================================================
#if DIMENSION==2
     u_value[0] = 0.;  u_value[1] = 0.;
//      if ( bc_gam==11 ) u_value[1] = 2.;
//      if ( xp[0]<0.0001 && xp[1]<-0.24999) u_value[1] = 2.;

#if FSI_EQUATIONS==1 //coupled
     double pref = _refvalue[DIMENSION];
     u_value[DIMENSION] = 0.*1000. /pref;
#endif
#endif
// ======================================================================================
#if DIMENSION==3
     u_value[0] = 0.;
     u_value[1] = 0.;
     u_value[2] = 0.;
#if FSI_EQUATIONS==1 //coupled
     double pref = _refvalue[DIMENSION];
     u_value[DIMENSION] =( 20+1.* ( 0.5-xp[2] ) ) *0.6;
#endif
#endif
// ======================================================================================
     return;
}
// ======================================================================================
// This function  defines the boundary conditions for the NS system:
// ======================================================================================
void MGSolFSI::bc_read (
     int bc_gam,
     int bc_mat,
     double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
     int bc_Neum[], 	// normal
     int bc_flag[]         // boundary condition flag
)
{
//  ============================================================================================
// number refer to tg components
// simm=0,velocity_in0=1,velocity_tg0=2,wall=3,
// simm1=4,velocity_in=5,velocity_tg=6,velocity=7,
// outflow=10,pressure_outlet=12,interior=11,
// wall_turb=13,outflow_p=14,pressure_inlet=16
//   
//  =====================================================================================
   double ILref = 1./_lref;

#if DIMENSION==2  // --------------------- 2D boundary conditions -----------------------
     bc_Neum[0] =wall;
     bc_flag[0]=0;
// mesh: channel_fsi1t.med                                                               
//      if ( bc_gam == 11 )      bc_Neum[0]=velocity;          //inlet                      
//      if ( bc_gam == 15 )      bc_Neum[0]=simm;              //left wall                  
//      if ( bc_gam == 13 )      bc_Neum[0]=outflow;           //outlet                     
//      if ( bc_gam == 1000 )    bc_Neum[0]=wall;              //FSI interface (top bottom) 
//      if ( bc_gam == 23 )      bc_Neum[0]=wall;              //top solid face             
//      if ( bc_gam == 21 )      bc_Neum[0]=wall;              //bot solid face             
//      if ( bc_gam == 20 )      bc_Neum[0]=outflow;           //right lateral wall         
//                                                                                          
//      if ( xp[0] < 0.0001 &&
//                xp[1] < -0.2499 )   bc_Neum[0]=velocity;       // left bottom point
  if( bc_gam==21 ){bc_Neum[0]=outflow;   } //bott solid
  if( bc_gam==22 ){bc_Neum[0]=wall;   } // lateral solid     
  if( bc_gam==12 ){bc_Neum[0]=outflow;   } //laterl fluid
  if( bc_gam==10 ){bc_Neum[0]=outflow;  }// top 
  if (bc_gam==1000){bc_Neum[0]=outflow;  }
               
               
               
#endif // -------------------------------- 2D boundary conditions -----------------------

#if DIMENSION==3 // ---------------------- 3D boundary conditions -----------------------

// FSI
     bc_Neum[0]=wall;
//   if (xp[2]< 0.0001) {bc_Neum[0]=velocity_in;} //inlet
//   if (xp[2]> .499999) {bc_Neum[0]=outflow;} //inlet
     if ( bc_gam == 13 ) {
          bc_Neum[0]=pressure_inlet;   //inlet
     }
     if ( xp[2]<0.0001 && xp[0]>0.1999 && xp[0]<2.00001 ) {
          bc_Neum[0]=wall;   //inlet
     }
     if ( bc_gam == 23 ) {
          bc_Neum[0]=wall;   //inlet
     }

     if ( bc_gam == 20 ) {
          bc_Neum[0]=outflow;   //wall
     }
//  if (bc_gam == 15) {bc_Neum[0]=simm;}        //wall
//
     if ( bc_gam == 11 ) {
          bc_Neum[0]=outflow_p;   //outlet
     }
     if ( xp[2]>0.4999999 && xp[0]>0.1999 && xp[0]<2.00001 ) {
          bc_Neum[0]=wall;   //inlet
     }
     if ( bc_gam == 21 ) {
          bc_Neum[0]=wall;   //wall
     }
//
//   if (bc_gam == 12) {bc_Neum[0]=wall;}     //outlet
//  if (bc_gam == 22) {bc_Neum[0]=wall;}        //wall
//    if (bc_gam == 14) {bc_Neum[0]=wall;}     //outlet
//  if (bc_gam == 24) {bc_Neum[0]=wall;}        //wall
     if ( xp[0]<0.001 ) bc_Neum[0]=simm;

     if ( xp[2]<0.0001 && xp[0]<.0001 )  bc_Neum[0]=pressure_inlet;
     if ( xp[1]>0.499999 ) bc_Neum[0]=wall;
     if ( xp[1]<0.0000001 ) bc_Neum[0]=wall;
#endif  //  // -------------------------------- 3D boundary conditions -----------------------
return;
}
// ============================================================================================
#endif  //ENDIF FSI_EQUATIONS
// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 
