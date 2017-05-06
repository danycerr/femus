

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




// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolT::ic_read(int bc_gam,int bc_mat, double xp[],double u_value[]) {
// =================================================
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(u,v,p)
//   u_value[0] =573.15;
  
  u_value[0] = 0.;//xp[0]*(2.-xp[0])*xp[1]*(1.-xp[1]);
  
//   -1.*(1.5+0.5*xp[0]-xp[1])*(xp[0]-1.)/**(1.-xp[0])*(1.5+0.5*xp[0]-xp[1])*/;
  
//   (0.5+2.*xp[0]*xp[0]);//*(xp[1]*xp[1]*xp[1]+2);
  
//   if(xp[0]>0.5) u_value[0] =  1.;//
//    if (xp[1]<0.001) u_value[0] = 573.;
     
#else
// =================================================
  // xp[]=(xp,ypzp) u_value[]=(u,v,w,p)
 u_value[0] = 0.;
#endif
}


// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolT::bc_read(
  int bc_gam,
  int bc_mat, 
  double xp[],
  int bc_Neum[],
  int bc_flag[]
) {
  // =================================================
  const double Lref = 1;//*this->get_par("Lref");
  double ILref = 1./Lref;
  int imesh=atoi(_mgutils.get_file("MESHNUMBER").c_str());
#if DIMENSION==2
// // xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
//   //    boundary conditions box
//   bc_Neum[0]=1;       bc_flag[0]=0;
//    // mesh1 simple_hex27.med
//   if(imesh==1){
  if( bc_gam==11){ 
    bc_Neum[0]=0; bc_flag[0]=0;
  }   //  lateral
    if( bc_gam==12){
    bc_Neum[0]=1; bc_flag[0]=0;
  }  // bottom
  if( bc_gam==13){ 
    bc_Neum[0]=0; bc_flag[0]=1;
  }   //  lateral
    if( bc_gam==21){ 
    bc_Neum[0]=0; bc_flag[0]=1;
  }

 
    bc_Neum[0]=0;       bc_flag[0]=0;
    if (xp[0]>LXE-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}
    if (xp[0]<LXB+BDRY_TOLL) {bc_Neum[0]=0; bc_flag[0]=0;}      
    if (xp[1]>LYE-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}
    if (xp[1]<LYB+BDRY_TOLL) {bc_Neum[0]=0; bc_flag[0]=1;} 
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box

    bc_Neum[0]=0;       bc_flag[0]=0;
    if (xp[0]>LXE-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}
    if (xp[0]<LXB+BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}      
    if (xp[1]>LYE-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}
    if (xp[1]<LYB+BDRY_TOLL) {bc_Neum[0]=0; bc_flag[0]=1;}
    if (xp[2]>LZE-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}
    if (xp[2]<LZB+BDRY_TOLL) {bc_Neum[0]=0; bc_flag[0]=0;}
#endif

  return;
} // end boundary conditions ==========================


#endif

