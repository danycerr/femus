

#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================
#include <sstream>
#include "MGGeomEl.h"
// configuration files -----------
#include "MGSclass_conf.h"
#include "MGFE_conf.h"
// #include "MGSTBKconf.h"
// #include "MGSTTBKconf.h"   // turbulent energy
#include "Printinfo_conf.h"
#include "MGEquationsSystem.h"

// class local configuration -------
#include "MGSolverT.h"

// local include -------------------
#include "MGMesh.h"
#include "MGSystem.h"
#include "MGUtils.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"

#include "parallelM.h"

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolT::bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
)   // ===================================
{
// ============================================================================
  if( mat_gam==4) {  bc_flag[0]=0;    bc_Neum[0]=3 ;  }
  else            { bc_flag[0]=0;     bc_Neum[0]=1 ;  }
  if (bc_gam==1000)  {    bc_flag[0]=0;  bc_Neum[0]=5 ;  }

  return;
}


// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolT::ic_read(int bc_gam,int bc_mat, double xp[],double u_value[])
{
// =================================================
  u_value[0] = 0.;
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(u,v,p)
//   u_value[0] =573.15;

  u_value[0] =  20.;//*(1.-xp[1]);
  if( (xp[0]-0.5)*(xp[0]-0.5)+(xp[1]-0.7)*(xp[1]-0.7)<0.2*0.2 )  u_value[0] =  800.;//*(1.-xp[1]);
//  if( bc_gam==10 ) u_value[0] =  800.;//*(1.-xp[1]);
//  if( bc_gam==21 ||bc_gam==22 ) u_value[0] =  20.;//*(1.-xp[1]);
//    if (xp[1]<0.001) u_value[0] = 573.;

#else
// =================================================
  // xp[]=(xp,ypzp) u_value[]=(u,v,w,p)
  u_value[0] = 20.;
  double xc=0.5; double yc=0.5; double zc=1; 
//   double zc=0.7; 
  double r=0.3;
 if( bc_gam==10 || xp[2]>0.999) u_value[0] =  800.;
  if((xp[0]-xc)*(xp[0]-xc)+(xp[1]-yc)*(xp[1]-yc)+(xp[2]-zc)*(xp[2]-zc)<(r)*(r)*(r))u_value[0] =  800.;
#endif
//   u_value[0] =10.;
}


// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolT::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int bc_Neum[],
  int bc_flag[]
)
{
  // =================================================
  const double Lref = 1;//*this->get_par("Lref");
  double ILref = 1./Lref;
  bc_flag[0]=0;     bc_Neum[0]=1 ;
  if( bc_gam==22 ) { //lateral solid
    bc_flag[0]=0;     bc_Neum[0]=3 ;
  }
  if( bc_gam==21 ) { //bottom
    bc_flag[0]=0;     bc_Neum[0]=3 ;
  }

  if( bc_gam==10  || xp[2]>0.999) { //top
    bc_flag[0]=1;     bc_Neum[0]=0;
  }
  if (bc_gam==1000) {
    bc_flag[0]=0;     bc_Neum[0]=5;
  }

  return;
} // end boundary conditions ==========================


#endif

