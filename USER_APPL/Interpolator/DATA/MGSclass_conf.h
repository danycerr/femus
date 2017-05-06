#include <Equations_conf.h>

//
/*===================================== Navier Stokes Equation ============================================*/
//
#ifdef NS_EQUATIONS
//

#ifndef __mgsnsc_h__
#define __mgsnsc_h__

// ==============================
// femlcore
#define D_EQ  0.0129
#define R_CORE  0.5408

//  ---------------------
// MULTIGRID PARAMETERS

// #define SOLVERNS VANKANSM
#define SOLVERNS GMRESM
// #define SOLVERNS BICGSTABM
#define SOLVERNSP CGM
// #define SOLVERNSP GMRESM


//  ---------------------------------
// 3D NAVIER-STOKES ADVECTION term
//
// A) Stokes flow    ADVPIC_NS 0. ADVNEW_NS=0
// B)  Navier-Stokes ADVPIC_NS 1. ADVNEW_NS={0,1}
#define ADVPIC_NS 0.
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_NS=1, ADVNEW_NS=0
// B) Newton iteration  ADVPIC_NS=1, ADVNEW_NS=1
#define ADVNEW_NS 0.


// -----------  stabilization ------------------
// --------------------------------------------
//  Navier-Stokes  stab=STAB_NS*0.5*(div u,u.v)
#define STAB_NS 0.
//  Navier-Stokes  compressibility=KOMP_NS*dp/dt
#define KOMP_NS 1.e-20
//  Navier-Stokes  compressibility  upwind=Re+UP_WIND_NS*v^2
#define UP_WIND_NS (1.)
//  Navier-Stokes  penalty =LAMBDA grad div u (solid)
#define LAMBDA (1.)
// Navier-Stokes supg
#define SUPG_NS (0.)
// Navier-Stokes  c -> antisymmetric (0.5)
#define ADV_ASYM 0.

// Navier-Stokes non linear iteration
// #define N_NNL_ITER 0

/* For Inclined Geometry */
#define INCLINED
#define OR_ANGLE 30



// Crank-Nicolson first order 0. 2nd order 0.5
#define CN_TIME 0.

//  boundary integral
// #define P_0 (0.)



#define SQCMU (0.3)
#define YPLUS (1.)
#define ALPHA0 (1.)

#endif


#endif /* End Navier Stokes Equation ----------------------------------------------------------------------*/
//
//
/*===================================== Energy Equation ===================================================*/
//
#ifdef T_EQUATIONS

#define V_MID (0.12085553719)  //(1.75867107438) //1.0234 // 1.06

#ifndef  __mgstconf_h__
#define  __mgstconf_h__

// turbulence  ----------------------
#define ADVE 1.

// Turbulence Prandl number
#define PRT (0.85)

// #define SOLVERT VANKATM -----------------
#define SOLVERT GMRESM



// temperature lows -------------------------
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)
#define  UP_WIND_T (1.)
// ------------------------
// #define  LINEAR 1 quad (2)
#define LQ_T (2)

// boundary flux
// #define Q0 (0.)

#endif

#endif /* End Energy Equation -----------------------------------------------------------------------------*/
//
//
/*===================================== DA Equation =======================================================*/
//
#ifdef DA_EQUATIONS
// ===============================->0.4268 ==============
#ifndef __mgsnsc_h__
#define __mgsnsc_h__



//  -----------------------------
// NAVIER-STOKES ADVECTION term
// --------------------------------
// A) Stokes flow   ADV 0. B)  Navier-Stokes ADV 1.
#define ADV 1.
//  Navier-Stokes  A) Picard iteration  ADV1 (0.) B) Newton ADV 1.
#define ADV1 0.
//  Navier-Stokes  stab  +0.5*(div u,u.v)
#define STAB 1.
//


#endif

#endif /* End DA Equation ---------------------------------------------------------------------------------*/
//
//
/*===================================== FSI Equation ======================================================*/
//
#ifdef FSI_EQUATIONS
// =============================================

#ifndef __mgsnscfsi_h__
#define __mgsnscfsi_h__


//
// MULTIGRID PARAMETERS
// -------------------------------
// Navier-Stokes solver type
#define SOLVER_FSI GMRESM  // options -> GMRESM BICGSTABM

// Pressure solver type (projection method only)
#define SOLVER_FSIP CGM  // options -> GMRESM CGM BICGSTABM


// ------------------------
//   SOLID MODEL
// ----------------------------
// geometric non-linearity
 #define NL_GEOM  (1)

// define penalty  (only with FSIP_EQUATIONS==1)
#define  PENALTY_FSI (100.)



 // ---------------------------
//   SOLID-FLUID REGIONS
// -----------------------------
//  #define SOLID 0
//  #define STIFF 10



//  --------------------------
// 3D NAVIER-STOKES ADVECTION term
// ---------------------------------
// A) Stokes flow    ADVPIC_SM 0. ADVNEW_SM=0
// B)  Navier-Stokes ADVPIC_SM 1. ADVNEW_SM={0,1}
#define ADVPIC_FSI 0.
//  Navier-Stokes  nonlinear iterations
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_SM=1, ADVNEW_SM=0
// B) Newton iteration  ADVPIC_SM=1, ADVNEW_SM=1
#define ADVNEW_FSI 0.


// -----------  stabilization ------------------
// --------------------------------------------
//  Navier-Stokes  stab=STAB_SM*0.5*(div u,u.v)
// #define STAB_FSI 0.
//  Navier-Stokes  compressibility=KOMP_SM*dp/dt
#define KOMP_FSI 1.e-20
//  Navier-Stokes  compressibility  upwind=Re+UP_WIND_SM*v^2
#define UP_WIND_FSI (1.)
//  Navier-Stokes  penalty =LAMBDA grad div u (solid)
// #define LAMBDA (2100)
// Navier-Stokes supg
// #define SUPG_FSI (0.)
// Navier-Stokes  c -> antisymmetric (0.5)
// #define ADV_ASYM 0.


// #define LQ (2)


// Crank-Nicolson first order 0. 2nd order 0.5
// #define CN_TIME 0.

//  boundary integral
// #define P_0 (0.)

// turbulent
// #define MU_LOW (1.e-12)
// #define MU_TOP (1.e+12)


// #define SQCMU (0.3)
// #define YPLUS (1.)
// #define ALPHA0 (1.)



// P solver ------------------------------------

// #define SOLVERT VANKATM ---------------------
// #define SOLVERFSIP GMRESM


// temperature lows -----------------------------
//  #define  CONST 1
// // #define densityT(x) (1.)
// // #define cipiT(x) (1.)
// // #define kappa(x)  (1.)
//
// // quadratic QL=2 linear QL=1 -------------------
// #define LQ (2)
// #define LQ_P (1)
// // boundary pressure --------------------------
// #define P_BD (1.)


#endif
#endif /* End FSI Equation --------------------------------------------------------------------------------*/
//
//
/*===================================== DS Equation =======================================================*/
//
#ifdef DS_EQUATIONS
// =============================================


#endif /* End DS Equation ---------------------------------------------------------------------------------*/
//
//
/*===================================== Dynamic Turbulence ================================================*/
//
#ifdef TBK_EQUATIONS
// =============================================
#ifndef  __mgstbkconf_h__
#define  __mgstbkconf_h__

// turbulence  -----------------------------
#define ADVE 1.

// Turbulence Prandl number
#define PRT (0.85)

// #define SOLVERT VANKATM-----------------------------
#define SOLVER_TBK GMRESM

// #define DIST_FIX (1.e-6)  // pay attention, only on the boundary!
// #define LES (100.)


#define MU_LOW  (1.e-6)                 /* Limit for minimum turbulent viscosity */
#define MU_TOP  (1.e+5)                 /* Limit for maximum turbulent viscosity */

#define MAX_TAU (1.e+10)                 /* Limit for maximum dynamic characteristic time */
#define MIN_TAU (1.e-10)                 /* Limit for minimum dynamic characteristic time */

#define UP_WIND_K (1.)

#define BETASTAR (0.09)
#define CMU (0.09)

//boundary
#define Y_WALL (0.00005)
#define SQCMU (0.3)
// #define KAPPA_VK (0.4)
//
/*-------------------------------------- The Dynamic K-E turbulence models ---------------------------------*/
//
#if ((TBK_EQUATIONS/2)==1)
//
/*-------------------------------------- For Nagano K-E ----------------------------------------------------*/
//
// #define NAGANO (1)
#ifdef NAGANO
 #define SIGMA_K (0.714285714286)
 #define SIGMA_E (0.714285714286)
 #define C2O (1.9)
 #define C1O (1.5)
#endif /* End Nagano K-E -----------------------------------------------------------------------------------*/
//
/*-------------------------------------- For LAUNDER -------------------------------------------------------*/
//
// #define LAUNDER (1)
#ifdef LAUNDER
 #define SIGMA_K (1.)
 #define SIGMA_E (0.769230769231)
 #define C2O (1.92)
 #define C1O (1.44)
#endif /* End LAUNDER --------------------------------------------------------------------------------------*/

#endif /* End Dynamic K-E turbulence models ----------------------------------------------------------------*/
//
/*-------------------------------------- The Dynamic K-W turbulence models ---------------------------------*/
//
#if ((TBK_EQUATIONS/2)==2)
//
/*-------------------------------------- For SST -----------------------------------------------------------*/
//
#ifdef SST
#define BETASTAR (0.09)
#define SIGMA_K (.5)
#define SIGMA_W (.5)
#endif /* End SST ------------------------------------------------------------------------------------------*/
//
/*-------------------------------------- For Nagano K-w ----------------------------------------------------*/
//
#define NAGANO (1)
#ifdef NAGANO
 #define CMU (0.09)
 #define SIGMA_W (0.714285714286)
 #define SIGMA_K (0.714285714286)
 #define C2O (1.9)
 #define C1O (1.5)
#endif /* End Nagano K-W -----------------------------------------------------------------------------------*/

// #define WILCOX (1)
#ifdef WILCOX
 #define CMU (0.09)
 #define BETASTAR (0.09)
 #define SIGMA_K (.5)
 #define SIGMA_W (.5)
 #define LOWRE
#endif /* End Wilcox K-W -----------------------------------------------------------------------------------*/


#endif
#endif
#endif /* End TBK_EQUATIONS -------------------------------------------------------------------------------*/
//
//
/*===================================== Thermal Turbulence ================================================*/
//
#ifdef TTBK_EQUATIONS

#ifndef  __mgsttbkconf_h__
#define  __mgsttbkconf_h__

// turbulence  ==========================
#define ADVE 1.

// #define SOLVERT VANKATM =======================
#define SOLVERTBK GMRESM

// temperature lows ================================
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)

// #define  LINEAR 1 quad (2)
//  #define LQ_TB  (2)

// #define LES (100.)

#define UP_WIND_TTK (1.)

//boundary
// #define SQCMU (0.3)
// #define KAPPA_VK (0.4)
//
/*-------------------------------------- For Nagano K-E ----------------------------------------------------*/
//
#if ((TTBK_EQUATIONS/2)==1)

#define MAX_TAUT (1.e+8)                 /* Limit for maximum thermal characteristic time */
#define MIN_TAUT (1.e-8)                 /* Limit for minimum thermal characteristic time */

#define SIGMA_EH (0.714285714286)
#define SIGMA_KH (0.714285714286)

#define CD1 (1.0)
#define CP1 (0.925)
#define CD2 (1.9)
#define CP2 (0.9)

#define CMU (0.09)

#endif /* End Nagano K-E -----------------------------------------------------------------------------------*/
//
/*-------------------------------------- The Thermal K-W turbulence models ---------------------------------*/
//
#if ((TTBK_EQUATIONS/2)==2)
//
#define MAX_TAUT (1.e+14)                /* Limit for maximum thermal characteristic time */
#define MIN_TAUT (1.e-14)                /* Limit for minimum thermal characteristic time */
#define CMU (0.09)
#define BETA (0.09)
//
/*-------------------------------------- For SST -----------------------------------------------------------*/
//
#ifdef SST
#define SIGMA_WH (0.714285714286)
#define SIGMA_KH (0.714285714286)
#define CD1 (1.4)
#define CP1 (1.1)
#define CD2 (0.8)
#define CP2 (0.6)
#endif /* End SST ------------------------------------------------------------------------------------------*/
//
/*-------------------------------------- For Nagano K-w ----------------------------------------------------*/
//
#define NAGANO (1)
#ifdef NAGANO
#define SIGMA_WH (0.714285714286)
#define SIGMA_KH (0.714285714286)
#define CD1 (0.1)
#define CP1 (0.025)
#define CD2 (1.9)
#define CP2 (0.9)
#endif /* End Nagano K-W -----------------------------------------------------------------------------------*/

#endif /* End K-W models -----------------------------------------------------------------------------------*/

#endif
#endif /* End TTBK_EQUATIONS ------------------------------------------------------------------------------*/
