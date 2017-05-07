
// C++ include files that we need
#include "Equations_conf.h"

#ifdef TWO_PHASE

// laspack includes
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <map>
#include "itersolv.h"

// vof includes and config
#include "vof_config.h"
#include "MGSolverCC.h"


// #include "MGFE_conf.h"
#include "Domain_conf.h"
#include "Solverlib_conf.h"
// MED includes
#ifdef HAVE_MED

#include "InterfaceFunctionDD.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif


// ---------------------------------
// --------- Constructors -----------
// ---------------------------------
/// Constructor
MGSolCC::MGSolCC(const unsigned int NoLevels_in,int name) :
  _dim(DIMENSION),_Nlev_cc(NoLevels_in) {
  _name=name;
  V_Constr(&_uvw, (char *) "uvw",1,Normal, _LPTrue);
  V_Constr(&_dsvw, (char *) "dsvw",1,Normal, _LPTrue);
  /* allocation of dynamic variables */
  for(int idim=0;idim<_dim;idim++) _nxyz[idim] = (size_t *) malloc(NoLevels_in * sizeof(size_t));
 
//   _nx = (size_t *) malloc(NoLevels_in * sizeof(size_t));
//   _ny = (size_t *) malloc(NoLevels_in * sizeof(size_t));
  _invnode_dof= (int *) malloc(sizeof(int));
  c1 = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
  c1_old = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
  _mx1 = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
  _my1 = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
  _mx1_old = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
  _my1_old = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
#if DIMENSION==3   
//    _nz = (size_t *)malloc(NoLevels_in * sizeof(size_t));
  _mz1_old = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
   _mz1 = (Matrix *) malloc(NoLevels_in * sizeof(Matrix));
#endif   
  return;
}

// ------------------------------------------------
/// Initial function (<- inital level)
void MGSolCC::init(
   const int ne_xyz[] ///< xyz line elements
//   const int nex_in,const int ney_in,const int nez_in
) {
  // initialization
  int npt_xyz[DIMENSION];int nvi=1;
  for(int idim=0;idim<_dim;idim++) {
    npt_xyz[idim]=ne_xyz[idim]+1;
     nvi *=npt_xyz[idim];
  }
//   const unsigned int nx=ne_xyz[0]+1;  const unsigned int ny=ne_xyz[1]+1;
//   const unsigned int nz=ne_xyz[2]+1;   // number of nodes
  fprintf(stderr,"-\n Solution: ");

  // Init level 0 of cc cc_old
//   int nvi=nx*ny*nz;
  V_Constr(&cc, (char *) "cc",nvi,Normal, _LPTrue);
  V_Constr(&cc_old, (char *) "cc_old",nvi,Normal, _LPTrue);
  V_Constr(&_mx, (char *) "mx",nvi,Normal, _LPTrue);
  V_Constr(&_my, (char *) "my",nvi,Normal, _LPTrue);
 #if DIMENSION==3 
  V_Constr(&_mz, (char *) "mz",nvi,Normal, _LPTrue);
#endif  
  init_level(0,npt_xyz);//nx,ny,nz);
  init_dofCA(0/*,mgmesh1,Levelmesh*/);

  // Init level matrix c1[Level] and c1[Level]_old
  int fcl=1;
  for(int level=1; level<_Nlev_cc; level++) {
    fcl *=2;
      for(int idim=0;idim<_dim;idim++) npt_xyz[idim]=ne_xyz[idim]*fcl+1;
//     init_level(level,(nx-1)*fcl+1,(ny-1) *fcl+1,(nz-1) *fcl+1);
    init_level(level,npt_xyz);
  }
  // Aux storage
//   int dimt=(_nx[_Nlev_cc-1]+1)*9;
  int dimt=( _nxyz[0][_Nlev_cc-1]+1)*9;
  _tempc1=new double[dimt];
  _tempc2=new double[dimt];
  _tempmx1=new double[dimt];
  _tempmy1=new double[dimt];
 #if DIMENSION==3 
  _tempmz1=new double[dimt];
#endif  
  GenSol(CCLEV-1);

  fprintf(stderr," - \n");
  return;
}
// ------------------------------------------------
/// Initial level
void MGSolCC::init_level(
  const int Level,    ///< cc level
  const int npt_xyz[] ///< xyz line points
) {
  // initialization
  for(int idim=0;idim<_dim;idim++) _nxyz[idim][Level]=npt_xyz[idim]; // storage class
   int nyz=1; for(int idim=1;idim<_dim;idim++) nyz *= npt_xyz[idim];
  
  char Name[30];
  sprintf(Name, "c1[%d]", Level % 1000);
  M_Constr(&c1[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "c1_old[%d]", Level % 1000);
  M_Constr(&c1_old[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "mx1[%d]", Level % 1000);
  M_Constr(&_mx1[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "my1[%d]", Level % 1000);
  M_Constr(&_my1[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "mx1_old[%d]", Level % 1000);
  M_Constr(&_mx1_old[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "my1_old[%d]", Level % 1000);
  M_Constr(&_my1_old[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
#if DIMENSION==3
  sprintf(Name, "mz1[%d]", Level % 1000);
  M_Constr(&_mz1[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
  sprintf(Name, "mz1_old[%d]", Level % 1000);
  M_Constr(&_mz1_old[Level],Name,nyz,npt_xyz[0],Rowws,Normal,_LPTrue);
#endif
  return;
}

// --------------------------------------------------------------------------------------
void MGSolCC::init_dofCA(
  const unsigned int Level
//   const MGMeshC & mesh_uw,
//   const int unsigned meshLev
) {
  //
  int n_nodes1=1; double   hxyz[DIMENSION];
   for(int idim=0;idim<_dim;idim++) {
     hxyz[idim]=1./(_nxyz[idim][Level]-1); // space grid
     n_nodes1 *=_nxyz[idim][Level];        // number of nodes
   }

  _invnode_dof =(int *) malloc((n_nodes1+1)*sizeof(int));
  // Translating into a Cartesian Grid
//   int icount=0;
//     for (unsigned int ix=0;ix< nx-1;ix++){
//         for (unsigned int iy=0;iy< ny-1;iy++){
//            for (unsigned int iz=0;iz< nz-1;iz++){
//          int index=ix+iy*ny+iz*nx*ny;
//           _invnode_dof[index] =icount+1;
//          icount++;
//         }
//         }
//     }
  for(unsigned int iu=0; iu< n_nodes1; iu++) {
//     double xp=mesh_uw._xyz[iu]; double yp=mesh_uw._xyz[iu+offset];
//     double zp=mesh_uw._xyz[iu+offset2];
//     int ix=(int)(xp/hx+0.5); int iy=(int)(yp/hy+0.5);
//     int iz=(int)(zp/hz+0.5);
//     int ind=ix+(iy+iz*ny)*nx;
//     _invnode_dof[ind] =ind+1;
     _invnode_dof[iu] =iu+1;
  }
  return;
}
// -----------

// --------------------------------------
// ----------- Destructors --------------
// --------------------------------------
/// Destructor (<-clear)
MGSolCC::~MGSolCC() {
  /* calling of destructors for matrices, vectors of unknowns,
     vectors of right hand side, restriction and prolongation operators */
  clear();
  /* release of dynamic variables */
  for(int idim=0;idim<_dim;idim++) {if(_nxyz[idim] != NULL) free(_nxyz[idim]); }
  
  if(c1 != NULL)  free(c1);  if(c1_old != NULL)  free(c1_old);
  if(_mx1 != NULL)  free(_mx1);  if(_my1 != NULL)  free(_my1);
  if(_mx1_old != NULL)  free(_mx1_old);  if(_my1_old != NULL)  free(_my1_old);
#if DIMENSION==3  
  if(_mz1_old != NULL)  free(_mz1_old);delete [] _tempmz1;
#endif
  delete [] _invnode_dof;
  delete [] _tempc1; delete [] _tempc2;
  delete [] _tempmx1; delete [] _tempmy1; 
  return;
}
// --------------------------------------------------------------------------------------
/// Clear
void MGSolCC::clear() {
  V_Destr(&cc);  V_Destr(&cc_old); V_Destr(&_mx); V_Destr(&_my); 
#if DIMENSION==3   
  V_Destr(&_mz);
#endif   
  for(unsigned int Level =0; Level< _Nlev_cc; Level++) {
    M_Destr(&c1[Level]);  M_Destr(&c1_old[Level]);
    M_Destr(&_my1[Level]);  M_Destr(&_mx1[Level]); 
    M_Destr(&_my1_old[Level]);  M_Destr(&_mx1_old[Level]); 
 #if DIMENSION==3   
   M_Destr(&_mz1[Level]); M_Destr(&_mz1_old[Level]);
#endif   
  }
  return;
}
// ======================================================================================
//  IO functions
// ======================================================================================

//---------------------------------------------------------------------------------------
/// write solution
void MGSolCC::WriteSol(const int Level,const std::string& name, const int var) {
  QVector *sol;
  if(var == 0) sol=&cc;
  else if(var == 1) sol=&_mx;
  else if(var == 2) sol=&_my;
#if DIMENSION==3  
  else if(var == 3) sol=&_mz;
#endif  
  std::ofstream outfile(name.c_str());
  // header  std::ofstream in(name.c_str());
  outfile << " Solution Level " << Level << std::endl;
  outfile << var << "  " << V_GetDim(sol) <<  _nxyz[0][0]   << std::endl;
  // data
  for(unsigned int i=1; i<=V_GetDim(sol); i++)   outfile  << V_GetCmp(sol,i) << "  ";
  outfile.close();
  std::cout << "   write solution cc  " << std::endl;
  return;
}



//--------------------------------------------------------------------------------------
/// read solution c1[Level]
void MGSolCC::ReadSol(const int Level,const std::string& name, const int var) {

  std::ifstream infile(name.c_str());
  QVector *sol;
  if(var == 0) sol=&cc; else if(var == 1) sol=&_mx;
  else if(var == 2) sol=&_my;
#if DIMENSION==3  
  else if(var == 3) sol=&_mz;
#endif  
  // header  std::ofstream in(name.c_str());
  const int  bufLen = 256; char  buf[bufLen+1];
  int dummy,dummy1; double ddummy;
  while(strncmp(buf,"Level",5) != 0) infile >> buf; infile >> dummy;
  if(Level != dummy) std::cout << "error input sol" << name << std::endl;
  // reshape
  infile >> buf  >> dummy >> dummy1;  _nxyz[0][0]=dummy1;
  // data
  for(unsigned int i=1; i<= (unsigned  int) dummy; i++) {
    infile >> ddummy;  V__SetCmp(sol,i,ddummy);
  }
  infile.close();
  std::cout << " Sol(cc) " << std::endl;
}










#endif // ======= TWO_PHASE =============================
