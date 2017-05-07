#include "vof_config.h"
// C++ include files that we need
#include "Equations_conf.h"
#include "MGFE_conf.h"
#include "Domain_conf.h"
#include "Solverlib_conf.h"
#ifdef TWO_PHASE
#if DIMENSION==3

#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <map>
#include "itersolv.h"
//#include "itersolv2.h
#include "MGSolverCC.h"
// #include "MGSolver.h"
//  #include "MGSolverT.h"
#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionDD.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif


#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#undef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN_VAL 1.e-12
#define PI 3.14159265358979
//#define PRINT_STEP 10
//#define ALPHA 1
//#define ELVIRA 1
#define YOUNG 1



// inline double PSI(double xi,double yj)  {
//   // single vortex
//   return (1.*sin(PI*xi) *sin(PI*xi) *sin(PI*yj) *sin(PI*yj) /PI);
//   // 4 vortices
//   // return (.25*sin(4.*pi*(xi+0.5))*cos(4*pi*(yj+0.5))/pi);
// }


// --------------------------------
//  IO functions
// --------------------------------
// /// write solution
// void MGSolCC::WriteSol(const int Level,const std::string& name, const int var) {
//   QVector *sol;
//   if(var == 0) sol=&cc;
//   else if(var == 1) sol=&_mx;
//   else if(var == 2) sol=&_my;
//   else if(var == 3) sol=&_mz;
//   std::ofstream outfile(name.c_str());
//   // header  std::ofstream in(name.c_str());
//   outfile << " Solution Level " << Level << std::endl;
//   outfile << var << "  " << V_GetDim(sol) <<  _nxyz[0][0]   << std::endl;
//   // data
//   for(unsigned int i=1; i<=V_GetDim(sol); i++)   outfile  << V_GetCmp(sol,i) << "  ";
//   outfile.close();
//   std::cout << "   write solution cc  " << std::endl;
//   return;
// }
// -------------------------------------------
/// read solution c1[Level]
// void MGSolCC::ReadSol(const int Level,const std::string& name, const int var) {
// 
//   std::ifstream infile(name.c_str());
//   QVector *sol;
//   if(var == 0) sol=&cc; else if(var == 1) sol=&_mx;
//   else if(var == 2) sol=&_my; else if(var == 3) sol=&_mz;
//   // header  std::ofstream in(name.c_str());
//   const int  bufLen = 256; char  buf[bufLen+1];
//   int dummy,dummy1; double ddummy;
//   while(strncmp(buf,"Level",5) != 0) infile >> buf; infile >> dummy;
//   if(Level != dummy) std::cout << "error input sol" << name << std::endl;
//   // reshape
//   infile >> buf  >> dummy >> dummy1;  _nxyz[0][0]=dummy1;
//   // data
//   for(unsigned int i=1; i<= (unsigned  int) dummy; i++) {
//     infile >> ddummy;  V__SetCmp(sol,i,ddummy);
//   }
//   infile.close();
//   std::cout << " Sol(cc) " << std::endl;
// }

// --------------------------
// Solution
// -----------------------------------
/// cc color function  update
void MGSolCC::CCSol_update() {
  const unsigned int nxl=_nxyz[0][0]; const unsigned int ny=_nxyz[1][0]-1;
  const unsigned int nz=_nxyz[2][0]-1;
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int jy=0; jy<ny; jy++) {
      // expanding c1_old -> tempc1
      ExpRow(&c1[0],_tempc1,0,jy+kz*ny,0);
      // back to cc
      int indyz=(jy+kz*(ny+1))*nxl;
      for(unsigned int ix=1; ix<nxl; ix++) cc.Cmp[ix+indyz]=_tempc1[ix];
    }// global grid jy
  return;
}
// ------------------------------------------
/// Generating solution at level Level
void MGSolCC::GenSol(const int Level) {

  unsigned int nx,ny,nz,dummy;
  std::ifstream infile("./RESU/c0.in");
  // reading from function init_circle(Level,0.5,.75,0.,0.15);
  if(!infile) {
    std::cout <<"Input file "<<infile<<" not opened." << std::endl;
    double xc[3];// drop center
    xc[0]=0.5; xc[1]=0.5; xc[2]=1.;
  
    init_drop(Level,xc,0.20);  
  } else {
    // reading from file "../RESU/c0.in"
    const int  bufLen = 256; char  buf[bufLen+1]; buf[0]='0';
    while(strncmp(buf,"Level",5) != 0)    infile >> buf;
    infile >> dummy;
    if(Level+1 != dummy) std::cout << " _NoLevels is " << _Nlev_cc << std::endl;
    // Cartesian dimension (nx,ny,nz)
    infile >> nx >> ny >> nz;
    if(nx != _nxyz[0][Level]-1 || ny != _nxyz[1][Level]-1|| nz != _nxyz[2][Level]-1) {
      std::cout << " Input file error  " << nx << " != " << _nxyz[0][Level]-1 << std::endl; exit(3);
    }
    // reading
    for(unsigned int k=0; k<nz; k++) for(unsigned int j=0; j<ny; j++) {
        for(unsigned int i=1; i<=nx; i++)  infile >> _tempc2[i];
        CtrRow(_tempc2,&c1[Level],Level,j+k*ny,0);
      }
  }
  OldSol_update(Level,c1[Level],c1_old[Level]);
  GenNorm(Level,&MGSolCC::rec_Young);
  // Restriction -----------------------------
  for(int level=Level-1; level>=0; level--) {
    RestSol(level); OldSol_update(level,c1[level],c1_old[level]);
    //   RestNormal(level);
  } // ---------------------------------
  // Projection -----------------------------
//   for (int level=Level+1;level<_Nlev_cc;level++) {
  //   ProjSol(level); OldSol_update(level,c1[level],c1_old[level]);
//   } // ------------------------------------
  // Update
  CCSol_update();
//  Normal_update();
  return;
}
// --------------  Old Solution cc ----------------------
/// Operator
/// \f$ Matrix^L(nx\times ny) \rightarrow Matrix^L(nx\times ny) :\f$
/// \f$ \quad c1[L] = c1_{old}[L]  \f$
void MGSolCC::GenOldSol(const int Level) {
  for(unsigned int ind=1; ind<_nxyz[1][Level]*_nxyz[2][Level]; ind++) {
    int len=M_GetLen(&c1[Level],ind);   M_SetLen(&c1_old[Level],ind,len);
    for(unsigned int k=0; k<len; k++)      M_SetEntry(&c1_old[Level],ind,k,k+1,0.);
  }
  OldSol_update(Level,c1[Level],c1_old[Level]);
  // std::cout << "   gen old solution  " << std::endl;
}
// -------------------------------------------------------
/// Update old solution
void MGSolCC::OldSol_update(const int Level,Matrix & sol,Matrix & sol_old) {
  for(unsigned int iy=1; iy<= (_nxyz[1][Level]-1)*(_nxyz[2][Level]-1); iy++) {
    int len=M__GetLen(&sol,iy); M_SetLen(&sol_old,iy,len);
    for(unsigned int kx=0; kx<len; kx++) M__SetEntry(&sol_old,iy,kx,M__GetPos(&sol,iy,kx),M__GetVal(&sol,iy,kx));
  }
  // std::cout << "  update old sol " << std::endl;
}


// ----------------------------------------------
void MGSolCC::init_drop(unsigned int Level,double xc[], double r) {
  // initilize a circle with (x0,y0) = center; r = radius
  double x0=xc[0]; double y0=xc[1]; double z0=xc[2];
  const unsigned int nx=_nxyz[0][Level]-1;   double hx=1./nx;
  const unsigned int ny=_nxyz[1][Level]-1;   double hy=1./ny;
  const unsigned int nz=_nxyz[2][Level]-1;   double hz=1./nz;
  // number of points in cell
    for(unsigned int k=0; k<nz; k++) 
      for(unsigned int j=0; j<ny; j++) {
        for(unsigned int i=1; i<=nx; i++) { 
	    unsigned int ind = i+(j+k*(ny+1))*(nx+1)+1;
	    const double xl = (i+.5)*hx;  const double yl = (j+0.5)*hy;
	    const double zl = (k+0.5)*hz;
	      _tempc2[i]=0.;
	  
// // // 	    drop  find_me
// 	    if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < (r+hx)*(r+hx)) _tempc2[i]=.5;
// 	    if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < r*r) _tempc2[i]=1.;
	    
// 	      	  if( (xl-0.5)*(xl-0.5)+ (yl-0.5)*(yl-0.5)<0.04 &&
// 	      zl>1.-1.*hz 
// 	  ){
// 	   _tempc2[i]=1.0;
// 	  }
	   
// // // //     //custom
//             if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < (r+hx)*(r+hx)) _tempc2[i]=.5;
//             if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < r*r) _tempc2[i]=1.;
//             if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < (r+hx)*(r+hx)) _tempc2[i]=.5;
//             if((xl-x0)*(xl-x0)+(yl-y0)*(yl-y0)+(zl-z0)*(zl-z0) < r*r) _tempc2[i]=1.;
	      
	    //custom
    //         if((xl-x0)*(xl-x0)+(yl-y0*0.5)*(yl-y0*0.5)+(zl-z0)*(zl-z0) < (r+hx)*(r+hx)) _tempc2[i]=.5;
    //         if((xl-x0)*(xl-x0)+(yl-y0*0.5)*(yl-y0*0.5)+(zl-z0)*(zl-z0) < r*r) _tempc2[i]=1.;
    //         if((xl-x0)*(xl-x0)+(yl-y0*0.25)*(yl-y0*0.25)+(zl-z0)*(zl-z0) < (r+hx)*(r+hx)) _tempc2[i]=.5;
    //         if((xl-x0)*(xl-x0)+(yl-y0*0.25)*(yl-y0*0.25)+(zl-z0)*(zl-z0) < r*r) _tempc2[i]=1.;
	      
// 	      double x1=hx,y1=hy,z1=hz;
// 	      double x2=0.35,y2=1-hy,z2=0.26;
	      
// 	      double xc=0.5,zc=0.7,r=0.025;
// 	      if(
// 		yl<3*hx+0.00001 && ((xl-xc)*(xl-xc)+(zl-zc)*(zl-zc)<(r+0.1)*(r+0.1))
// 	      ){
// 		_tempc2[i]=0.5;
// 	      }
// 	      if(
// 		yl<2*hx+0.00001 && ((xl-xc)*(xl-xc)+(zl-zc)*(zl-zc)<r*r)
// 	      ){
// 		_tempc2[i]=1.;
// 	      }
	      
	      
// 	    if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)){
// 	      _tempc2[i]=1.;
// 	      if((yl<y1+hy && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)) _tempc2[i]=.5;
// 	      if((yl<y2 && yl>y2-hy) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)) _tempc2[i]=.5;
// 	      if((yl<y2 && yl>y1) && (xl>x1 && xl<x1+hx) && (zl>z1 && zl<z2))_tempc2[i]=.5;
// 	      if((yl<y2 && yl>y1) && (xl>x2-hx && xl<x2) && (zl>z1 && zl<z2)) _tempc2[i]=.5;
// 	      if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z1+hz))_tempc2[i]=.5;
// 	      if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z2-hz && zl<z2)) _tempc2[i]=.5;
// 	    }
	    
	    
	      ic_read(_tempc2,i,j,k,hx,hy,hz,0);
	    
        }
        CtrRow(_tempc2,&c1[Level],Level,j+k*ny,0);
      }
//   std::cout<< "Initializing c function on a circle of center (" << x0 << ","
//            << y0 << ") and radius " << r << std::endl;
//   std::cout<< nx << " x " << nx << "  "<< hx
//            << " h " << hx << std::endl;
}

// ---------------------------------
// Normal functions
// -----------------------------------
// ------------------------------------
/// Normal Generation from color function
// void MGSolCC::GenNorm(const int Level) {
void MGSolCC::GenNorm(const unsigned int Level,
                      void (MGSolCC::*rec1)(const unsigned int Level,int ixyz[], double mxyz[])) {

  int nxl=_nxyz[0][Level]; int ny=_nxyz[1][Level]-1; int nz=_nxyz[2][Level]-1;
  int ixyz[3];
  // Set up
  for(unsigned int ix=0; ix<9*(nxl+1); ix++) {
    _tempc1[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
  }             //line 2

  for(int kk=0; kk<2; kk++)
    for(int jj=0; jj<2; jj++) {
      int indk=jj+kk*ny; int indkl=(jj+kk*3)*(nxl+1);
      ExpRow(&c1_old[Level],_tempc1,Level,indk,indkl);
      ExpRow(&_mx1_old[Level],_tempmx1,Level,indk,indkl);
      ExpRow(&_my1_old[Level],_tempmy1,Level,indk,indkl);
      ExpRow(&_mz1_old[Level],_tempmz1,Level,indk,indkl);
    }

// loop --------------------------------------
  for(unsigned int kz=0; kz<nz; kz++) {
    // boundary condition




    for(unsigned int jy=0; jy<ny; jy++) {
      int ind=jy+kz*ny;
      int len=M__GetLen(&c1_old[Level],ind+1);
      M_SetLen(&_mx1[Level],ind+1,len);
      M_SetLen(&_my1[Level],ind+1,len);
      M_SetLen(&_mz1[Level],ind+1,len);

      // loop over columns
      for(unsigned int kx=0; kx<len; kx++) {
        double mxyz[4];
        int ix=M__GetPos(&c1_old[Level],ind+1,kx)-1;
        ixyz[0]=ix; ixyz[1]=jy; ixyz[2]=kz;

        double val=M__GetVal(&c1_old[Level],ind+1,kx);
        // computing normal ( 0<c<1)
        if(val <1.)(this->*rec1)(Level,ixyz,mxyz);

        M__SetEntry(&_mx1[Level],ind+1,kx,ix+1,mxyz[0]);
        M__SetEntry(&_my1[Level],ind+1,kx,ix+1,mxyz[1]);
        M__SetEntry(&_mz1[Level],ind+1,kx,ix+1,mxyz[2]);
      }

      // Row update
      for(int kkup=-1; kkup<2; kkup++) {
        int indjt=jy+2+(kz+kkup)*ny; int indjl=((jy-1+3)%3+((kkup+kz+3)%3)*3)*(nxl+1);
        ExpRow(&c1_old[Level],_tempc1,Level,indjt,indjl);
        ExpRow(&_mx1_old[Level],_tempmx1,Level,indjt,indjl);
        ExpRow(&_my1_old[Level],_tempmy1,Level,indjt,indjl);
        ExpRow(&_mz1_old[Level],_tempmz1,Level,indjt,indjl);
      }


    }
    // Initialize new k line
    for(unsigned int ix=0; ix<9*(nxl+1); ix++) {
      _tempc1[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
    }             //line 2
    for(int k0=kz; k0<kz+3; k0++)
      for(int j0=0; j0<2; j0++) {
        int indkt=k0*ny+j0; int indkl=((j0%3)+(k0%3)*3)*(nxl+1);
        ExpRow(&c1_old[Level],_tempc1,Level,indkt,indkl);
        ExpRow(&_mx1_old[Level],_tempmx1,Level,indkt,indkl);
        ExpRow(&_my1_old[Level],_tempmy1,Level,indkt,indkl);
        ExpRow(&_mz1_old[Level],_tempmz1,Level,indkt,indkl);
      }
//           printf(" %d %d %d \n ",jy+2,0,Level);
//      for (int kxc=0; kxc<nxl; kxc++) printf(" %e ",_tempc1[kxc+((jy+2)%3)*nxl]);
//        printf(" \n \n");
//
//          for (int kxc=0; kxc<M__GetLen(&c1_old[Level],jy+3); kxc++) printf(" %e ",M__GetVal(&c1_old[Level],jy+3,kxc));
//        printf(" \n \n");
//



    //  }
  }
  return;
}
// ------------------------------------------------------------------
/// unit normal update
void MGSolCC::Normal_update() {
  const unsigned int nxl=_nxyz[0][0]; const unsigned int ny=_nxyz[1][0]-1;
  const unsigned int nz=_nxyz[2][0]-1;
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int jy=0; jy<ny; jy++)  {
      int ind=jy+kz*ny;
      int indyz=(jy+kz*(ny+1))*(nxl);
      // expanding mx1 -> tempmx1
      for(unsigned int ix=1; ix<nxl; ix++) {
        _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
      }
      for(unsigned int kx=0; kx<M_GetLen(&c1[0],ind+1); kx++)  {
        int pos=M__GetPos(&c1[0],ind+1,kx);
        _tempmx1[pos]=M__GetVal(&_mx1[0],ind+1,kx);
        _tempmy1[pos]=M__GetVal(&_my1[0],ind+1,kx);
        _tempmz1[pos]=M__GetVal(&_mz1[0],ind+1,kx);
      }
      // back to cc
      for(unsigned int ix=1; ix<nxl; ix++) {
        _mx.Cmp[ix+indyz+1]=_tempmx1[ix];
        _my.Cmp[ix+indyz+1]=_tempmy1[ix];
        _mz.Cmp[ix+indyz+1]=_tempmz1[ix];
      }
    }// global grid jy
  return;
}


// ----------------------------------------------------

// ------------------------------------
// Multigrid  Operators
// ------------------------------------

// -------------------------------------------------------
/// 3D solver function
void MGSolCC::MGSolve(/*MGMeshC &mgmesh,*/
                      unsigned int nc_step,unsigned int itime,double dt) {



  for(unsigned int is=0; is<nc_step; is++) {
    // advection
    Adv(_Nlev_cc-1,dt/(nc_step), itime);
    OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
    // normal
#ifdef ALPHA
    GenNorm(_Nlev_cc-1,&MGSolCC::rec_melv1);
#else
    GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
#endif
    OldSol_update(_Nlev_cc-1,_mx1[_Nlev_cc-1],_mx1_old[_Nlev_cc-1]);
    OldSol_update(_Nlev_cc-1,_my1[_Nlev_cc-1],_my1_old[_Nlev_cc-1]);
    OldSol_update(_Nlev_cc-1,_mz1[_Nlev_cc-1],_mz1_old[_Nlev_cc-1]);
    // Restriction
    for(int level=_Nlev_cc-2; level>=0; level--) {
      RestSol(level); RestNormal(level);
    }
    CCSol_update();
    Normal_update();
  }
#ifndef NS_EQUATIONS
//  print(mgmesh,nc_step,_Nlev_cc-1);
#endif
  return;
}
// --------------------------------------------------------
///  Operator for color function from Coarse Level  C1[Levc] to
///  Fine Level  c1[Levf]   \f$ C_1[L_f] = {\bf P}  C_1[L_c] \f$
void MGSolCC::ProjSol(const int Levf) {
  return;
}
// -----------------------------------
/// restriction operator for color function c1[Level fine]-> c1[Level coarse]
void MGSolCC::RestSol(const int Levc) {
  int Levf=Levc+1;
  const unsigned int nxlc=_nxyz[0][Levc]; const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxlf=_nxyz[0][Levf]; const unsigned int nyf=_nxyz[1][Levf]-1;
  const unsigned int nzc=_nxyz[2][Levc]-1; const unsigned int nzf=_nxyz[2][Levf]-1;

  //OldSol_update(Levf);
  for(unsigned int kzc=0; kzc<nzc; kzc++)  {
    int kzf=kzc*2;
    for(unsigned int jyc=0; jyc<nyc; jyc++)  {
      int jyf=jyc*2;
      //for (unsigned int ixf=0; ixf<=(nxf+1)*2; ixf++){_tempc1[ixf]=0.;_tempc2[ixf]=0.;}

//      Expanding matrix c1_old -> new line _tempc1
      ExpRow(&c1[Levf],_tempc1,Levf,jyf+kzf*nyf,0); // line 1
      ExpRow(&c1[Levf],_tempc1,Levf,jyf+1+kzf*nyf,nxlf+1); // line 2

      ExpRow(&c1[Levf],_tempc1,Levf,jyf+(kzf+1)*nyf,2*(nxlf+1)); // line 1
      ExpRow(&c1[Levf],_tempc1,Levf,jyf+1+(kzf+1)*nyf,3*(nxlf+1)); // line 2
      // Area restriction
      for(unsigned int ixc=1; ixc<nxlc; ixc++) {
        int ixf=ixc*2-1;
        _tempc2[ixc]=0.125* (
                       _tempc1[ixf]+_tempc1[ixf+1]+_tempc1[ixf+nxlf+1]+_tempc1[ixf+nxlf+2]+
                       _tempc1[ixf+2*(nxlf+1)]+_tempc1[ixf+1+2*(nxlf+1)]+_tempc1[ixf+nxlf+1+2*(nxlf+1)]+_tempc1[ixf+nxlf+2+2*(nxlf+1)]);
      }
      // Contracting matrix tempc2 -> c1
      CtrRow(_tempc2,&c1[Levc],Levc,jyc+kzc*nyc,0);

    }// global grid jy
  }// global grid kz

  //CtrRow(_tempc2,&c1[Levc],Levc,8+8*nyc,0);
  return;
}
// -----------------------------------
/// for normal _mxy1[Level c]->_mxy1[Level f]
void MGSolCC::ProjNormal(const int Levf) {

  return;
}
// -----------------------------------
/// restriction operator for color function mxy1Level fine]-> mxy1[Level coarse]
void MGSolCC::RestNormal(const int Levc) {

  int Levf=Levc+1;
  const unsigned int nxlc=_nxyz[0][Levc]; const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxlf=_nxyz[0][Levf]+1; const unsigned int nyf=_nxyz[1][Levf]-1;
  const unsigned int nzc=_nxyz[2][Levc]-1; const unsigned int nzf=_nxyz[2][Levf]-1;
  //OldSol_update(Levf);

//OldSol_update(Levf);
  for(unsigned int kzc=0; kzc<nzc; kzc++)  {
    int kzf=kzc*2;
    for(unsigned int jyc=0; jyc<nyc; jyc++)  {
      int jyf=jyc*2;

//      Expanding matrix  mx1 -> new line _tempmx1
      ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf+kzf*nyf,0); // line 1
      ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf+1+kzf*nyf,nxlf); // line 2
      ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf+(kzf+1)*nyf,2*nxlf); // line 1
      ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf+1+(kzf+1)*nyf,3*nxlf);  // line 2
      // Expanding matrix my1 -> new line _tempmy1
      ExpRow(&_my1[Levf],_tempmy1,Levf,jyf+kzf*nyf,0); // line 1
      ExpRow(&_my1[Levf],_tempmy1,Levf,jyf+1+kzf*nyf,nxlf); // line
      ExpRow(&_my1[Levf],_tempmy1,Levf,jyf+(kzf+1)*nyf,2*nxlf); // line 1
      ExpRow(&_my1[Levf],_tempmy1,Levf,jyf+1+(kzf+1)*nyf,3*nxlf);  // line
      // Expanding matrix mz1 -> new line _tempmz1
      ExpRow(&_mz1[Levf],_tempmz1,Levf,jyf+kzf*nyf,0); // line 1
      ExpRow(&_mz1[Levf],_tempmz1,Levf,jyf+1+kzf*nyf,nxlf); // line
      ExpRow(&_mz1[Levf],_tempmz1,Levf,jyf+(kzf+1)*nyf,2*nxlf); // line 1
      ExpRow(&_mz1[Levf],_tempmz1,Levf,jyf+1+(kzf+1)*nyf,3*nxlf);  // line 2

      // Area restriction
      int lenc=M__GetLen(&c1[Levc],jyc+kzc*nyc+1);
      M_SetLen(&_mx1[Levc],jyc+kzc*nyc+1,lenc);
      M_SetLen(&_my1[Levc],jyc+kzc*nyc+1,lenc);
      M_SetLen(&_mz1[Levc],jyc+kzc*nyc+1,lenc);
      for(unsigned int kxc=0; kxc<lenc; kxc++) {
        int ixc= M__GetPos(&c1[Levc],jyc+kzc*nyc+1,kxc);
        int ixf=ixc*2-1;
        double mx=_tempmx1[ixf]+_tempmx1[ixf+1]+_tempmx1[ixf+nxlf]+_tempmx1[ixf+nxlf+1]+ _tempmx1[ixf+2*nxlf]+_tempmx1[ixf+1+2*nxlf]+_tempmx1[ixf+nxlf+2*nxlf]+_tempmx1[ixf+nxlf+1+2*nxlf];
        double my=_tempmy1[ixf]+_tempmy1[ixf+1]+_tempmy1[ixf+nxlf]+_tempmy1[ixf+nxlf+1]+ _tempmy1[ixf+2*nxlf]+_tempmy1[ixf+1+2*nxlf]+_tempmy1[ixf+nxlf+2*nxlf]+_tempmy1[ixf+nxlf+1+2*nxlf];
        double mz=_tempmz1[ixf]+_tempmz1[ixf+1]+_tempmz1[ixf+nxlf]+_tempmz1[ixf+nxlf+1]+ _tempmz1[ixf+2*nxlf]+_tempmz1[ixf+1+2*nxlf]+_tempmz1[ixf+nxlf+2*nxlf]+_tempmz1[ixf+nxlf+1+2*nxlf];
        // normalization
        double mm = fabs(mx)+ fabs(my)+ fabs(mz);
        if(mm > 0) {mx = mx/mm; my = my/mm; mz = mz/mm;  }
        M__SetEntry(&_mx1[Levc],jyc+kzc*nyc+1,kxc,ixc,mx);
        M__SetEntry(&_my1[Levc],jyc+kzc*nyc+1,kxc,ixc,my);
        M__SetEntry(&_mz1[Levc],jyc+kzc*nyc+1,kxc,ixc,mz);
      }

    }// global grid jy
  }// global grid kz



  return;
}
// -------------------------------------------------------------
/// The 3d advection function
void MGSolCC::Adv(unsigned int Level, double dt, unsigned int itime) {


  const unsigned int nx=_nxyz[0][Level]-1; const unsigned int ny=_nxyz[1][Level]-1; const unsigned int nz=_nxyz[2][Level]-1;

  // Volume computation
  double sum=0.;
  for(unsigned int i=0; i<ny*nz; i++) {
    int len=M__GetLen(&c1[Level],i+1);
    for(unsigned int kx=0; kx<len; kx++) {
      sum +=M__GetVal(&c1[Level],i+1,kx);
    }
  }
  std::cout << "VOF SOLVING...\n";
  std::cout << nx << " x " << ny << " x " << nz << " cc= " << sum << std::endl;
  
  
//   // Split volume advection
  lagrangeX(Level,dt/3.,itime);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   //fprintf(stderr,"\n Split: x ");
  lagrangeY(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
  //fprintf(stderr,"- y ");
  lagrangeZ(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"- z direction ");


  lagrangeY(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"\n Split: y ");
  lagrangeZ(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"- z ");
  lagrangeX(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"- x direction ");


  lagrangeZ(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"\n Split: z ");
  lagrangeX(Level,dt/3.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_Young);
//   fprintf(stderr,"- x ");
  lagrangeY(Level,dt/3.);
//   fprintf(stderr,"- y direction ");

  return;
}

// ---------------------
// RECONSTRUCTION
// ---------------------
// -------------------------------------
/// reconstruction by using Young stensil,
///.brutal finite difference across the cell with weight
/// {[1,2,1][2,4,2][1,2,1]} over each plane
void MGSolCC::rec_Young(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nxyz[0][Level]+1; const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nz=_nxyz[2][Level]-1; double min_err = 16.;
  const int ix=ixyz[0]; const int jy=ixyz[1]; const int kz=ixyz[2];


  Matrix & mcc=c1[Level];
  // index def
  int ind9[3][3][3];
  for(int kk=0; kk<3; kk++)
    for(int jj=0; jj<3; jj++)
      for(int ii=0; ii<3; ii++)
        ind9[ii][jj][kk]=((ix+1)+ii-1+nxl)%nxl +((jy+jj-1)%3+((kz+kk-1)%3)*3)*nxl;

  // unit normal
  double m1,m2,mxt,myt,mzt,alphat;
  if(kz !=0) {
    // z direction
    m1=(_tempc1[ind9[0][0][0]]+_tempc1[ind9[0][2][0]]+_tempc1[ind9[2][0][0]]+_tempc1[ind9[2][2][0]])+
       (_tempc1[ind9[0][1][0]]+_tempc1[ind9[2][1][0]]+_tempc1[ind9[1][0][0]]+_tempc1[ind9[1][2][0]])*2.+
       (_tempc1[ind9[1][1][0]])*4.;
    m2=(_tempc1[ind9[0][0][2]]+_tempc1[ind9[0][2][2]]+_tempc1[ind9[2][0][2]]+_tempc1[ind9[2][2][2]])+
       (_tempc1[ind9[0][1][2]]+_tempc1[ind9[2][1][2]]+_tempc1[ind9[1][0][2]]+_tempc1[ind9[1][2][2]])*2.+
       (_tempc1[ind9[1][1][2]])*4.;
    mzt=m1-m2;
    // y-direction
    m1=(_tempc1[ind9[0][0][0]]+_tempc1[ind9[0][0][2]]+_tempc1[ind9[2][0][0]]+_tempc1[ind9[2][0][2]])+
       (_tempc1[ind9[0][0][1]]+_tempc1[ind9[2][0][1]]+_tempc1[ind9[1][0][0]]+_tempc1[ind9[1][0][2]])*2.+
       (_tempc1[ind9[1][0][1]])*4.;
    m2=(_tempc1[ind9[0][2][0]]+_tempc1[ind9[0][2][2]]+_tempc1[ind9[2][2][0]]+_tempc1[ind9[2][2][2]])+
       (_tempc1[ind9[0][2][1]]+_tempc1[ind9[2][2][1]]+_tempc1[ind9[1][2][0]]+_tempc1[ind9[1][2][2]])*2.+
       (_tempc1[ind9[1][2][1]])*4.;
    myt=m1-m2;
    // x-direction
    m1=(_tempc1[ind9[0][0][0]]+_tempc1[ind9[0][0][2]]+_tempc1[ind9[0][2][0]]+_tempc1[ind9[0][2][2]])+
       (_tempc1[ind9[0][0][1]]+_tempc1[ind9[0][2][1]]+_tempc1[ind9[0][1][2]]+_tempc1[ind9[0][1][0]])*2.+
       (_tempc1[ind9[0][1][1]])*4.;
    m2=(_tempc1[ind9[2][0][0]]+_tempc1[ind9[2][0][2]]+_tempc1[ind9[2][2][0]]+_tempc1[ind9[2][2][2]])+
       (_tempc1[ind9[2][0][1]]+_tempc1[ind9[2][2][1]]+_tempc1[ind9[2][1][2]]+_tempc1[ind9[2][1][0]])*2.+
       (_tempc1[ind9[2][1][1]])*4.;
    mxt=m1-m2;
  } else {
    // boundary
    mzt=0.;
    mxt=(_tempc1[ind9[0][2][1]]+2.*_tempc1[ind9[0][1][1]]+_tempc1[ind9[0][0][1]])   -(_tempc1[ind9[2][2][1]]+2.*_tempc1[ind9[2][1][1]]+_tempc1[ind9[2][0][1]]);
    myt=(_tempc1[ind9[2][0][1]]+2.*_tempc1[ind9[1][0][1]]+_tempc1[ind9[0][0][1]])   -(_tempc1[ind9[2][2][1]]+2.*_tempc1[ind9[1][2][1]]+_tempc1[ind9[0][2][1]]);
  }
  // normalization
  double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
  if(mm >0.) {   mxt /= mm; myt /= mm; mzt /= mm; }
  mxyz[0]=mxt; mxyz[1]=myt;; mxyz[2]=mzt;
  return;
}


// ---------------------------------------------------
// 3D alpha reconstruction
///  alpha reconstruction (elvira modified for low resolutions)
///  by using contracing matrices c_old and mx1,my1,mz1
// --------------------------------------------------
void MGSolCC::rec_melv1(unsigned int Level,int ixyz[],double mxyz[]) {

//                int ind = i+(j+k*ny)*nx; const double  c_0=cc[Level].Cmp[ind+1];
//
//         /* only if  0.< c <1.*/
//         if (c_0 > 0. && c_0 < 1. ) {
//           double mmx; double mmy; double mmz;
//
//           mmx=_mx1[Level].Cmp[ind+1];  mmy=_my1[Level].Cmp[ind+1];
//           mmz=_mz1[Level].Cmp[ind+1];
//
//           for (int i1=-1;i1<2;i1++) {
//             for (int j1=-1;j1<2;j1++) {
//               for (int k1=-1;k1<2;k1++) {
//                 int indx=ind+i1+(j1+k1*ny)*nx;
//                 cl[i1+1][j1+1][k1+1]=p.Cmp[indx+1];
//               }
//             }
//           }
//
//           if (fabs(mmx)>fabs(mmy)) {
//             if (fabs(mmx)>fabs(mmz)) {
//               int i12 = (mmx > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=i12 || j1 !=1 || k1 !=1) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mx1[Level].Cmp[indx+1]*mmx <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(i1-i12) < 2) {
//                           int indx1=ind+i1-i12+(j1+k1*ny)*nx;
//                           if (_mx1[Level].Cmp[indx1+1]*mmx <0) cl[i1-i12+1][j1+1][k1+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//             else {
//               int k12 = (mmz > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1  || j1 !=1 || k1 !=k12) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mz1[Level].Cmp[indx+1]*mmz <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(k1-k12) < 2) {
//                           int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
//                           if (_mz1[Level].Cmp[indx1+1]*mmz <0) cl[i1+1][j1+1][k1-k12+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//           }
//           else {
//             if (fabs(mmy)>fabs(mmz)) {
//               int j12 = (mmy > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1 || j1 !=j12 || k1 !=1) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_my1[Level].Cmp[indx+1]*mmy <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(j1-j12) < 2) {
//                           int indx1=ind+i1+(j1-j12+k1*ny)*nx;
//                           if (_my1[Level].Cmp[indx1+1]*mmy <0) cl[i1+1][j1-j12+1][k1+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//             else {
//               int k12 = (mmz > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1 || j1 !=1 || k1 !=k12) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mz1[Level].Cmp[indx+1]*mmz <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(k1-k12) < 2) {
//                           int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
//                           if (_mz1[Level].Cmp[indx1+1]*mmz <0) cl[i1+1][j1+1][k1-k12+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//           }

  const int ix=ixyz[0];
  const int jy=ixyz[1];
  const int kz=ixyz[2];


  double mmx,mmy,mmz,alp; double min_err = 16.;
  unsigned int nxl=_nxyz[0][Level]+1; unsigned int ny=_nxyz[1][Level]-1;
  unsigned int nz=_nxyz[2][Level]-1;
  int ind[3][3][3]; double ccl[3][3][3];

  // initial crw and ccl
  for(unsigned int k0=0; k0<3; k0++)
    for(unsigned int j0=0; j0<3; j0++)
      for(unsigned int i0=0; i0<3; i0++) {
        ind[i0][j0][k0]=(ix+1)+i0-1+ ((jy+j0-1+3) %3+((kz+k0-1+3) %3)*3)*(nxl);
        ccl[i0][j0][k0]=_tempc1[ind[i0][j0][k0]];
      }



  double ccc=_tempc1[ind[1][1][1]];
  double amx=_tempmx1[ind[1][1][1]];
  double amy=_tempmy1[ind[1][1][1]];
  double amz=_tempmz1[ind[1][1][1]];

  if(fabs(amx)>fabs(amy)) {
    // amx > (amy,amz)
    if(fabs(amx)>fabs(amz)) {
      int i12 = (amx > 0.)? 1:-1;
      for(int i1=-1; i1<2; i1++)
        for(int j1=-1; j1<2; j1++)
          for(int k1=-1; k1<2; k1++) {
            if(i1 !=i12 || j1 !=1 || k1 !=1) {
              //int indx=ind+i1+(j1+k1*ny)*nx;
              double mxl=_tempmx1[ind[i1+1][j1+1][k1+1]];
              if(mxl*amx < 0.) {
                ccl[i1+1][j1+1][k1+1]=1.;
                if(fabs(i1-i12) < 2) {
                  //int indx1=ind+i1-i12+(j1+k1*ny)*nx;
                  double mxl2=_tempmx1[ind[i1-i12+1][j1+1][k1+1]];
                  if(mxl2*amx < 0.) ccl[i1-i12+1][j1+1][k1+1]=1.;
                }
              }
            }
          }
    } // end amx > amy &&    amx > amz)
    else {
      int k12 = (amz > 0.)? 1:-1;
      for(int i1=-1; i1<2; i1++)
        for(int j1=-1; j1<2; j1++)
          for(int k1=-1; k1<2; k1++) {
            if(i1 !=1  || j1 !=1 || k1 !=k12) {
              //int indx=ind+i1+(j1+k1*ny)*nx;
              double mzl=_tempmz1[ind[i1+1][j1+1][k1+1]];
              if(mzl*amz < 0.) {
                ccl[i1+1][j1+1][k1+1]=1.;
                if(fabs(k1-k12) < 2) {
                  //int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
                  double mzl2=_tempmx1[ind[i1+1][j1+1][k1-k12+1]];
                  if(mzl2*amz <0) ccl[i1+1][j1+1][k1-k12+1]=1.;
                }
              }
            }
          }
    }
  } else {
    if(fabs(amy)>fabs(amz)) {
      int j12 = (amy > 0.)? 1:-1;
      for(int i1=-1; i1<2; i1++)
        for(int j1=-1; j1<2; j1++)
          for(int k1=-1; k1<2; k1++) {
            if(i1 !=1 || j1 !=j12 || k1 !=1) {
              //  int indx=ind+i1+(j1+k1*ny)*nx;
              double myl=_tempmy1[ind[i1+1][j1+1][k1+1]];
              if(myl*mmy <0) {
                ccl[i1+1][j1+1][k1+1]=1.;
                if(fabs(j1-j12) < 2) {
                  // int indx1=ind+i1+(j1-j12+k1*ny)*nx;
                  double myl2=_tempmy1[ind[i1+1][j1-j12+1][k1+1]];
                  if(myl2*mmy <0) ccl[i1+1][j1-j12+1][k1+1]=1.;
                }
              }
            }
          }
    } else {
      int k12 = (amz > 0.)? 1:-1;
      for(int i1=-1; i1<2; i1++)
        for(int j1=-1; j1<2; j1++)
          for(int k1=-1; k1<2; k1++) {
            if(i1 !=1 || j1 !=1 || k1 !=k12) {
              // int indx=ind+i1+(j1+k1*ny)*nx;
              double mzl=_tempmz1[ind[i1+1][j1+1][k1+1]];
              if(mzl*amz <0) {
                ccl[i1+1][j1+1][k1+1]=1.;
                if(fabs(k1-k12) < 2) {
                  // int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
                  double mzl2=_tempmz1[ind[i1+1][j1+1][k1-k12+1]];
                  if(mzl2*amz <0) ccl[i1+1][j1+1][k1-k12+1]=1.;
                }
              }
            }
          }
    }
  }
  double m1=ccl[0][0][0]+ccl[0][2][0]+ ccl[2][0][0]+ccl[2][2][0]+
            2*(ccl[0][1][0]+ccl[2][1][0]+ ccl[1][0][0]+ccl[1][2][0])+4*(ccl[1][1][0]);
  double m2=ccl[0][0][2]+ccl[0][2][2]+ ccl[2][0][2]+ccl[2][2][2]+
            2*(ccl[0][1][2]+ccl[2][1][2]+ ccl[1][0][2]+ccl[1][2][2])+4*(ccl[1][1][2]);
  mmz=m1-m2;
  m1=ccl[0][0][0]+ccl[2][0][0]+ ccl[0][0][2]+ccl[2][0][2]+
     2*(ccl[1][0][0]+ccl[2][0][1]+ ccl[1][0][2]+ccl[0][0][1])+4*(ccl[1][0][1]);
  m2=ccl[0][2][0]+ccl[2][2][0]+ ccl[0][2][2]+ccl[2][2][2]+
     2*(ccl[1][2][0]+ccl[2][2][1]+ ccl[1][2][2]+ccl[0][2][1])+4*(ccl[1][2][1]);
  mmy=m1-m2;
  m1=ccl[0][0][0]+ccl[0][2][0]+ ccl[0][0][2]+ccl[0][2][2]+
     2*(ccl[0][1][0]+ccl[0][2][1]+ ccl[0][1][2]+ccl[0][0][1])+4*(ccl[0][1][1]);
  m2=ccl[2][0][0]+ccl[2][2][0]+ ccl[2][0][2]+ccl[2][2][2]+
     2*(ccl[2][1][0]+ccl[2][2][1]+ ccl[2][1][2]+ccl[2][0][1])+4*(ccl[2][1][1]);
  mmx=m1-m2;

  // rec_Young(Level,ind,&mmx,&mmy,&mmz);
// double mmm=sqrt(mmx*mmx+mmy*mmy+mmz*mmz+1.e-10);
  double mmm=(fabs(mmx)+fabs(mmy)+fabs(mmz)+1.e-16);
  //	  printf("\n rec  mmx %e mm %e mmz %e \n",mmx,mmy,mmz);
  mxyz[0] = mmx/mmm;   mxyz[1] = mmy/mmm; mxyz[2] = mmz/mmm;

  return;
}

// // -----------------------------------------------------------------
// /** geometrical unsplit */
// void MGSolCC::rec_melv(unsigned int Level) {
//
//
//   const unsigned int nx=_nxyz[0][Level];const unsigned int ny=_nxyz[1][Level];
//   const unsigned int nz=_nxyz[2][Level];const unsigned int nvi=nx*ny*nz;
//   const double hx=1./nx;const   double hy=1./ny;const double hz=1./nz;
//   const int *dofn=_invnode_dof[Level];   const unsigned int NM=nx;
//   QVector &p=cc[Level];double cl[3][3][3];
//
//   // normal
//   for (unsigned int i=1 ; i<nx-1 ; i++) {
//     for (unsigned int j=1 ; j<ny-1 ; j++) {
//       for (unsigned int k=1 ; k<nz-1 ; k++) {
//         int ind = i+(j+k*ny)*nx; const double  c_0=cc[Level].Cmp[ind+1];
//
//         /* only if  0.< c <1.*/
//         if (c_0 > 0. && c_0 < 1. ) {
//           double mmx; double mmy; double mmz;
//
//           mmx=_mx1[Level].Cmp[ind+1];  mmy=_my1[Level].Cmp[ind+1];
//           mmz=_mz1[Level].Cmp[ind+1];
//
//           for (int i1=-1;i1<2;i1++) {
//             for (int j1=-1;j1<2;j1++) {
//               for (int k1=-1;k1<2;k1++) {
//                 int indx=ind+i1+(j1+k1*ny)*nx;
//                 cl[i1+1][j1+1][k1+1]=p.Cmp[indx+1];
//               }
//             }
//           }
//
//           if (fabs(mmx)>fabs(mmy)) {
//             if (fabs(mmx)>fabs(mmz)) {
//               int i12 = (mmx > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=i12 || j1 !=1 || k1 !=1) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mx1[Level].Cmp[indx+1]*mmx <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(i1-i12) < 2) {
//                           int indx1=ind+i1-i12+(j1+k1*ny)*nx;
//                           if (_mx1[Level].Cmp[indx1+1]*mmx <0) cl[i1-i12+1][j1+1][k1+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//             else {
//               int k12 = (mmz > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1  || j1 !=1 || k1 !=k12) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mz1[Level].Cmp[indx+1]*mmz <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(k1-k12) < 2) {
//                           int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
//                           if (_mz1[Level].Cmp[indx1+1]*mmz <0) cl[i1+1][j1+1][k1-k12+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//           }
//           else {
//             if (fabs(mmy)>fabs(mmz)) {
//               int j12 = (mmy > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1 || j1 !=j12 || k1 !=1) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_my1[Level].Cmp[indx+1]*mmy <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(j1-j12) < 2) {
//                           int indx1=ind+i1+(j1-j12+k1*ny)*nx;
//                           if (_my1[Level].Cmp[indx1+1]*mmy <0) cl[i1+1][j1-j12+1][k1+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//             else {
//               int k12 = (mmz > 0.)? 1:-1;
//               for (int i1=-1;i1<2;i1++)
//                 for (int j1=-1;j1<2;j1++)
//                   for (int k1=-1;k1<2;k1++) {
//                     if (i1 !=1 || j1 !=1 || k1 !=k12) {
//                       int indx=ind+i1+(j1+k1*ny)*nx;
//                       if (_mz1[Level].Cmp[indx+1]*mmz <0) {
//                         cl[i1+1][j1+1][k1+1]=1.;
//                         if (fabs(k1-k12) < 2) {
//                           int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
//                           if (_mz1[Level].Cmp[indx1+1]*mmz <0) cl[i1+1][j1+1][k1-k12+1]=1.;
//                         }
//                       }
//                     }
//                   }
//             }
//           }
//
//           double m1=cl[0][0][0]+cl[0][2][0]+ cl[2][0][0]+cl[2][2][0]+
//                     2*(cl[0][1][0]+cl[2][1][0]+ cl[1][0][0]+cl[1][2][0])+4*(cl[1][1][0]);
//           double m2=cl[0][0][2]+cl[0][2][2]+ cl[2][0][2]+cl[2][2][2]+
//                     2*(cl[0][1][2]+cl[2][1][2]+ cl[1][0][2]+cl[1][2][2])+4*(cl[1][1][2]);
//           mmz=m1-m2;
//           m1=cl[0][0][0]+cl[2][0][0]+ cl[0][0][2]+cl[2][0][2]+
//              2*(cl[1][0][0]+cl[2][0][1]+ cl[1][0][2]+cl[0][0][1])+4*(cl[1][0][1]);
//           m2=cl[0][2][0]+cl[2][2][0]+ cl[0][2][2]+cl[2][2][2]+
//              2*(cl[1][2][0]+cl[2][2][1]+ cl[1][2][2]+cl[0][2][1])+4*(cl[1][2][1]);
//           mmy=m1-m2;
//           m1=cl[0][0][0]+cl[0][2][0]+ cl[0][0][2]+cl[0][2][2]+
//              2*(cl[0][1][0]+cl[0][2][1]+ cl[0][1][2]+cl[0][0][1])+4*(cl[0][1][1]);
//           m2=cl[2][0][0]+cl[2][2][0]+ cl[2][0][2]+cl[2][2][2]+
//              2*(cl[2][1][0]+cl[2][2][1]+ cl[2][1][2]+cl[2][0][1])+4*(cl[2][1][1]);
//           mmx=m1-m2;
//
//           // rec_Young(Level,ind,&mmx,&mmy,&mmz);
//           double mmm=sqrt(mmx*mmx+mmy*mmy+mmz*mmz+1.e-10);
//           //	  printf("\n rec  mmx %e mm %e mmz %e \n",mmx,mmy,mmz);
//           _mx1[Level].Cmp[ind+1] = mmx/mmm;   _my1[Level].Cmp[ind+1] = mmy/mmm;
//           _mz1[Level].Cmp[ind+1] = mmz/mmm;
//         }
//         else {
//           _mx1[Level].Cmp[ind+1] = 0.;   _my1[Level].Cmp[ind+1] = 0.;
//           _mz1[Level].Cmp[ind+1] = 0.;
//         }
//       }
//     }
//   }
//   return;
// }
//
//


// ---------------------------------------------------------------
// Split advection along the x-direction:
/// Standard advection in one direction by using the normal from
/// (mx1 my1 mz1) previously computed by Young function
void  MGSolCC::lagrangeX(const unsigned int Level,const double dt, unsigned int itime) {



  const unsigned int nxl=_nxyz[0][Level]+1;  const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nz=_nxyz[2][Level]-1;
  const double hx=1./(nxl-2); const double hy=1./ny; const double hz=1./nz;
  unsigned int nxc=_nxyz[0][0]-1;  int fd=(nxl-2)/nxc;
// double *c1;c1=new double[n_nodes];
  double vof1,vof2,vof3,mmx,mmy,mmz,tnot;
  double mxt,myt,mzt,alphat;
  double a1,a2; int ixyz[3]; double utemp[6];double dstemp[6];
  double st=dt/(hx);
  double time=itime*dt;
  // First line
  for(unsigned int ix=0; ix< nxl*9; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;
    _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
  }
  for(unsigned int k0=0; k0<2; k0++) for(unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }

  /*c1=0.*/
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for(int k=0; k<nz ; k++) {
    for(int j=0; j<ny ; j++) {

      for(int i=0; i<nxl-2 ; i++) {
        ixyz[0]=i; ixyz[1]=j; ixyz[2]=k;
        //  int indv = i+(j+k*(ny+1))*nxl;
	get_disp(ixyz,fd,dstemp);
        get_vel(ixyz,fd,utemp);
        a1=utemp[0]*st; a2=utemp[1]*st;

//         double ix3=i*hx;double iy3=j*hy;
//         a1 =-(PSI(ix3,(iy3+hy))-PSI(ix3,iy3)) *dt/(hy*hx);
//         a2 =-(PSI(ix3+hx,iy3+hy)-PSI(ix3+hx,iy3)) *dt/(hy*hx);
        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        vof1 = 0.; vof2 = 0.; vof3 = 0.;
	
	
	
	//////////////////////////////////////////////////

	//internal boundary
	double x1=0.1,y1=hy,z1=hz;
	double x2=0.35,y2=1-hy,z2=0.26;
	double xl=i*hx, yl= j*hy, zl=k*hz;
// 	if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)){
// 	  _tempc1[indl+1]=0.9999999;
// 	}
// 	if( xl>0.45 && xl<0.55 && 
// 	    yl>0.45 && yl<0.55 && 
// 	    zl>1-2*hz-0.000001 && zl<1-hz){
// 	  _tempc1[indl+1]=1.;
// 	}
	
// // // 	find_me2 here we can put sources
	
//     double xc=0.5;double yc=0.5;double zc=1.;double r=0.3;
//     if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)<(r+2*hx)*(r+2*hy) && zl>1-2*hz){
//       _tempc1[indl+1]=0.8;
//      if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)<(r)*(r)) _tempc1[indl+1]=0.99;
// 
//     }
	
	
// 	if(
// 	      zl>0.59 && zl<0.71 
// 	    && yl>0.44 && yl<0.56
// 	    && xl<1*hx+0.00001
// 	  
// 	){
// 	  _tempc1[indl+1]=1.;
// 	}
	
// 	double yc=0.5,zc=0.6,r=0.05;
// 	if(
// 	  xl<2*hx+0.00001 && ((yl-yc)*(yl-yc)+(zl-zc)*(zl-zc)<(r+0.1)*(r+0.1))
// 	){
// 	  _tempc1[indl+1]=0.5;
// 	}
// 	if(
// 	  xl<1*hx+0.00001 && ((yl-yc)*(yl-yc)+(zl-zc)*(zl-zc)<r*r)
// 	){
// 	  _tempc1[indl+1]=1.;
// 	}
	/////////////Neuman outlet
// 	double xc=0.5,zc=0.7,r=0.025;
// 	if(
// 	  yl<3*hx+0.00001 && ((xl-xc)*(xl-xc)+(zl-zc)*(zl-zc)<(r+0.1)*(r+0.1))
// 	){
// 	   _tempc1[indl+1]=0.9;
// 	}
// 	if(
// 	  yl<2*hx+0.00001 && ((xl-xc)*(xl-xc)+(zl-zc)*(zl-zc)<r*r)
// 	){
// 	  _tempc1[indl+1]=1.;
// 	}
	
	////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////
	//Boundary check LUCA //check if cc=1 in boundary cells
// 	if(
// 	  i*hx<0.+1.*hx || i*hx>1-1.*hx ||
// 	  j*hy<0.+1.*hy || j*hy>1-1.*hy ||
// 	  k*hz<0.+1.*hz || k*hz>1-1.*hz 
// 
// 	){
// 	  if(_tempc1[indl+1] > 0.9999999) { //=1
// 	    _tempc1[indl+1]=0.9999999;
// 	  }
// 	}
	//Boundary check END
	//////////////////////////////////////////////////

	////////////////////////////////////////////////////
	//Outlet condition LUCA
// 	if(_tempc1[indl+1] > 0.0) {
	  //2nd cell near boundary
	  //cells on boundary
// 	  if( (xl-0.5)*(xl-0.5)+ (yl-0.5)*(yl-0.5)>0.05 &&
// 	      k*hz>1.-1.*hz -0.0001
// // 	      && !(
// // 	      i*hx>0.44 && i*hx<0.56 && 
// // 	      j*hy>0.44 && j*hy>0.44 
// // 	    )	//top no rubinetto
// 	  ){
// // 	   _tempc1[indl+1]=0.0000001;
// 	  }
//         }
// 	Outlet condition END
	/////////////////////////////////////////////////
	
        const double  c_0=_tempc1[indl+1];
	
        // c=1.
        if(c_0 > 1.-1e-16) {
          vof1 =c_0* MAX(-a1,0.);
          vof2 =c_0* (1.- MAX(a1,0.)-MAX(-a2,0.));
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if(c_0 > 0.) {
          // unit normal
          mxt=_tempmx1[indl+1]; myt=_tempmy1[indl+1]; mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          if(mm >0.) {   mxt /= mm; myt /= mm; mzt /= mm; }
          // alpha
          alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          mmx=mxt; mmy=myt;; mmz=mzt; tnot=alphat;
          // splitting correction factor
          mmx=mmx/(1.-a1+a2); tnot=tnot+mmx*a1;
          if(a1<0.) vof1=get_vol3D(mmx,mmy,mmz,tnot,a1,-a1);
          if(a2>0.) vof3=get_vol3D(mmx,mmy,mmz,tnot,1.,a2);
          vof2=get_vol3D(mmx,mmy,mmz,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // temporary storage
        _tempc2[indl+1]+=vof2; _tempc2[indl]+=vof1; _tempc2[indl+2]+=vof3;
	    double xc=0.5;double yc=0.5;double zc=1.;double r=0.3*0.5;
    if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)<(r+1*hx)*(r+1*hy) && zl>1-2*hz){
      _tempc2[indl+1]=((time>1.)?0:1)*0.8;
     if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)<(r)*(r)) _tempc2[indl+1]=((time>1.)?0:1)*0.99;

    }
      }
      // Contracting _tempc2 -> c1
      CtrRow(_tempc2,&c1[Level],Level,j+k*ny,((j%3)+(k%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for(int k2=-1; k2<2; k2++) {
        int indj=j+2+(k+k2)*ny; int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
      }

    }
    // new block for new k line
    // storage fo the last line in j
    // CtrRow(_tempc2,c1[Level],Level,ny-1+(k)*ny,((ny-1)%3+((k)%3)*3)*nxl);
//     // zero the temp storage
    for(unsigned int ix=0; ix< nxl *3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
    }
//     // set up the temp stprage for new k line
    for(unsigned int k0=k; k0<k+3; k0++) for(unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}

// ---------------------------------------------------------------
// Split advection along the y-direction
/// Standard advection in one direction by using the normal from
/// (mx1 my1 mz1) previously computed by Young function.
/// The matrices mx1,my1,mz1 are in contracted form
void  MGSolCC::lagrangeY(const unsigned int Level,const double dt, unsigned int itime) {

  const unsigned int nxl=_nxyz[0][Level]+1;  const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nz=_nxyz[2][Level]-1;
  const double hx=1./(nxl-2); const double hy=1./ny; const double hz=1./nz;
  unsigned int nxc=_nxyz[0][0]-1; int fd=(nxl-2)/nxc;
  double vof1,vof2,vof3,mmx,mmy,mmz,tnot;
  double mxt,myt,mzt,alphat;
  double a1,a2; int ixyz[3]; double utemp[6];double dstemp[6];
  double st=dt/(hx);
  // First line
  for(unsigned int ix=0; ix< nxl *3*3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;
    _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
  }
  for(unsigned int k0=0; k0<2; k0++) for(unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }

  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for(int k=0 ; k<nz ; k++) {
    for(int j=0 ; j<ny ; j++) {
      for(int i=0 ; i<nxl-2 ; i++) {

        //int indv = i+(j+k*(ny+1))*(nxl-1);
        ixyz[0]=i; ixyz[1]=j; ixyz[2]=k;
	get_disp(ixyz,fd,dstemp);
        get_vel(ixyz,fd,utemp);
        a1=utemp[2]*st; a2=utemp[3]*st;

        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        vof1 = 0.; vof2 = 0.; vof3 = 0.;
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if(c_0 >= 1.) {
          vof1 =c_0* MAX(-a1,0.);
          vof2 =c_0* (1. - MAX(a1,0.) - MAX(-a2,0.));
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if(c_0 > 0.) {
          mxt=_tempmx1[indl+1]; myt=_tempmy1[indl+1]; mzt=_tempmz1[indl+1];
          // mxt=0.;myt=0.;mzt=1.;
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          if(mm >0.) { mxt /= mm; myt /= mm; mzt /= mm;  }
          alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          mmx=mxt; mmy=myt;; mmz=mzt; tnot=alphat;

          mmy=mmy/(1.-a1+a2); tnot=tnot+mmy*a1;
          if(a1<0.) vof1=get_vol3D(mmy,mmx,mmz,tnot,a1,-a1);
          if(a2>0.)  vof3=get_vol3D(mmy,mmx,mmz,tnot,1.,a2);
          vof2=get_vol3D(mmy,mmx,mmz,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // new values of c1
        _tempc2[indl+1]+=vof2;
        _tempc2[i+((j-1+3)%3+(k%3)*3)*nxl+1]+=vof1;  _tempc2[i+((j+1)%3+(k%3)*3)*nxl+1] +=vof3;
        //if(fabs(vof1)+fabs(vof3) > 1.e-15) printf("\n vof vof! \n");
      }

      // contracting  tempc2 -> c2
      if(j>0) CtrRow(_tempc2,&c1[Level],Level,j-1+k*ny,((j-1)%3+(k%3)*3)*nxl);

      //Expanding matrix c1_old ->line  _tempc1
      for(int k2=-1; k2<2; k2++) {
        int indjt=j+2+(k+k2)*ny;   int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indjt,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indjt,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indjt,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indjt,indjl);
      }

    }
    // new block for new k line
    // storage fo the last line in j
    CtrRow(_tempc2,&c1[Level],Level,ny-1+k*ny,((ny-1)%3+(k%3)*3)*nxl);
    // zero the temp storage
    for(unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for(unsigned int k0=k; k0<k+3; k0++) for(unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}
// ---------------------------------------------------------------
// Split advection along the z-direction
/// Standard advection in one direction by using the normal from
/// (mx1 my1 mz1) previously computed by Young function.
/// The matrices mx1,my1,mz1 are in contracted form
// ---------------------------------------------------------------
void  MGSolCC::lagrangeZ(const unsigned int Level,
                         const double dt, unsigned int itime) {
  // Set up
  const unsigned int nxl=_nxyz[0][Level]+1;  const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nz=_nxyz[2][Level]-1;
  const double hx=1./(nxl-2); const double hy=1./ny; const double hz=1./nz;
  unsigned int nxc=_nxyz[0][0]-1;  int fd=(nxl-2)/nxc;
  double vof1,vof2,vof3,mmx,mmy,mmz,tnot;
  double mxt,myt,mzt,alphat;
  double a1,a2; int ixyz[3]; double utemp[6];double dstemp[6];
  double st=dt/(hx);
  // First line
  for(unsigned int ix=0; ix< nxl*3*3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;
    _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
  }
  for(unsigned int k0=0; k0<2; k0++) for(unsigned int j0=0; j0<2; j0++) {
      int ind=k0*ny+j0; int indl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,ind,indl);
      ExpRow(&_mx1[Level],_tempmx1,Level,ind,indl);
      ExpRow(&_my1[Level],_tempmy1,Level,ind,indl);
      ExpRow(&_mz1[Level],_tempmz1,Level,ind,indl);
    }

  // Domain cartesian loop
  for(int j=0 ; j<ny ; j++) {
    for(int k=0 ; k<nz ; k++) {
      for(int i=0 ; i<nxl-2 ; i++) {
        ixyz[0]=i; ixyz[1]=j; ixyz[2]=k;

        //a1 =-0.;
        //a2 =-0.;
	get_disp(ixyz,fd,dstemp);
        get_vel(ixyz,fd,utemp);
        a1=utemp[4]*st; a2=utemp[5]*st;

        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        vof1 = 0.; vof2 = 0.; vof3 = 0.;
// // // 	aaaaaaa
	double xl=i*hx, yl= j*hy, zl=k*hz;
// 		if((xl-0.5)*(xl-0.5)+ (yl-0.5)*(yl-0.5)<(0.2+hx)*(0.2+hx) && 
// 	    zl> 1-2*hz-0.000001 
// 	  )
// 	{
// 	  _tempc2[indl+1]=1.;
// 	}
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if(c_0 >= 1.) {
          vof1 =c_0*MAX(-a1,0.);
          vof2 =c_0*(1.-MAX(a1,0.)-MAX(-a2,0.));
          vof3 =c_0*MAX(a2,0.);
        }
        //
        else if(c_0 > 0.) {
          // unit normal
          mxt=_tempmx1[indl+1]; myt=_tempmy1[indl+1]; mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          if(mm >0.) {
            mxt /= mm; myt /= mm; mzt /= mm;
          }
          alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          mmx=mxt; mmy=myt;; mmz=mzt; tnot=alphat;

          mmz=mmz/(1.-a1+a2); tnot=tnot+mmz*a1;
          if(a1<0.)   vof1=get_vol3D(mmz,mmx,mmy,tnot,a1,-a1);
          if(a2>0.)   vof3=get_vol3D(mmz,mmx,mmy,tnot,1.,a2);
          vof2=get_vol3D(mmz,mmx,mmy,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // new values of c1 (in temp2)
        _tempc2[indl+1]+=vof2;
        _tempc2[i+((j%3)+((k-1+3)%3)*3)*nxl+1]+=vof1; _tempc2[i+(j%3+((k+1)%3)*3)*nxl+1]+=vof3;
//     if((xl-0.5)*(xl-0.5)+ (yl-0.5)*(yl-0.5)<(0.2+4*hx)*(0.2+4*hx) && zl> 1-2*hz-0.000001 )
// 	{
// 	 _tempc2[indl+1]=0.2;
// 	}
//       if((xl-0.5)*(xl-0.5)+ (yl-0.5)*(yl-0.5)<(0.2)*(0.2) && zl> 1-2*hz-0.000001 )
// 	{
// 	 _tempc2[indl+1]=0.7;
// 	}
      }

      // contraction temp2 -> c1
      if(k>0) CtrRow(_tempc2,&c1[Level],Level,j+(k-1)*ny,((j)%3+((k-1)%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for(int j2=-1; j2<2; j2++) {
        int indk=j+j2+(k+2)*ny; int indkl=((j+j2+3)%3+((k+2)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indkl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indkl);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indkl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indk,indkl);
      }

    }
    // new block for new j line
    // storage fo the last line in k
    CtrRow(_tempc2,&c1[Level],Level,j+(nz-1)*ny,(j%3+((nz-1)%3)*3)*nxl);
    // zero the temp storage
    for(unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.; _tempmy1[ix]=0.; _tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for(unsigned int k0=0; k0<2; k0++) for(unsigned int j0=j; j0<j+3; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0%3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}

// -------------------------------------------------
double  MGSolCC::get_vol3D(double m1,double m2,double m3,double alpha,double r0,
                           double dr0) {
  double al0,b1,b2,b3,b12,bm,tmp,pr,vol;

  /* (1) move origin to r0 along r; (2) reflect parallelepiped;
     (3) limit alpha (0<=al0<=0.5); (4) order coefficients b1<b2<b3;
     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)   */
  /* (1) */
  alpha = alpha - m1*r0;
  /* (2) */
  alpha = alpha + MAX(0.0,-m1*dr0) + MAX(0.0,-m2) + MAX(0.0,-m3);
  tmp = fabs(m1)*dr0 + fabs(m2) + fabs(m3);
  m1 = fabs(m1)/tmp;
  m2 = fabs(m2)/tmp;
  m3 = fabs(m3)/tmp;
  alpha = MAX(0.0,MIN(1.0,alpha/tmp));
  /* (3) */
  al0 = MIN(alpha,1.0-alpha);
  /* (4) */
  b1 = MIN(m1*dr0,m2);
  b3 = MAX(m1*dr0,m2);
  b2 = m3;
  if(b2 < b1) {
    tmp = b1;
    b1  = b2;
    b2  = tmp;
  } else if(b2 > b3) {
    tmp = b3;
    b3  = b2;
    b2  = tmp;
  }
  b12 = b1 + b2;
  bm = MIN(b12,b3);
  pr = MAX(6.*b1*b2*b3,1.0e-50);
  /* (5) */
  if(al0 < b1)
    tmp = al0*al0*al0/pr;
  else if(al0 < b2)
    tmp = 0.5*al0*(al0-b1)/(b2*b3) + b1*b1*b1/pr;
  else if(al0 < bm)
    tmp = (al0*al0*(3.0*b12-al0) + b1*b1*(b1-3.0*al0) +
           b2*b2*(b2-3.0*al0))/pr;
  else if(bm == b12)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3.0-2.0*al0) + b1*b1*(b1-3.0*al0) +
           b2*b2*(b2-3.0*al0) + b3*b3*(b3-3.0*al0))/pr;

  if(alpha <= 0.5) return  tmp*dr0;
  return (1.0-tmp)*dr0;
}
// ---------------------------------------------------------
// void  MGSolCC::print_cc(MGMeshC &mgmesh,std::ofstream& out,const unsigned int offset,
//                         const unsigned int Level_vof,const unsigned int mode) {
//   const unsigned int level=mgmesh._NoLevels-1-Level_vof; int nx=_nxyz[0][Level_vof];
//   const unsigned int goffset=mgmesh._NoNodes[mgmesh._NoLevels-1];
//   const unsigned int goffset2=goffset*2;
//   double hx=1./(nx-1);
//   if(mode == 1) {
//     int conn[8][8];
//     conn[0][0]=0; conn[0][1]=8; conn[0][2]=20; conn[0][3]=11; conn[0][4]=12; conn[0][5]=21; conn[0][6]=26; conn[0][7]=24;
//     conn[1][0]=8; conn[1][1]=1; conn[1][2]=9; conn[1][3]=20; conn[1][4]=21; conn[1][5]=13; conn[1][6]=22; conn[1][7]=26;
//     conn[2][0]=11; conn[2][1]=20; conn[2][2]=10; conn[2][3]=3; conn[2][4]=24; conn[2][5]=26; conn[2][6]=23; conn[2][7]=15;
//     conn[3][0]=20; conn[3][1]=9; conn[3][2]=2; conn[3][3]=10; conn[3][4]=26; conn[3][5]=22; conn[3][6]=14; conn[3][7]=23;
// 
//     conn[4][0]=12; conn[4][1]=21; conn[4][2]=26; conn[4][3]=24; conn[4][4]=4; conn[4][5]=16; conn[4][6]=25; conn[4][7]=19;
//     conn[5][0]=21; conn[5][1]=13; conn[5][2]=22; conn[5][3]=26; conn[5][4]=16; conn[5][5]=5; conn[5][6]=17; conn[5][7]=25;
//     conn[6][0]=24; conn[6][1]=26; conn[6][2]=23; conn[6][3]=15; conn[6][4]=19; conn[6][5]=25; conn[6][6]=18; conn[6][7]=7;
//     conn[7][0]=26; conn[7][1]=22; conn[7][2]=14; conn[7][3]=23; conn[7][4]=25; conn[7][5]=17; conn[7][6]=6; conn[7][7]=18;
// 
//     for(unsigned int el=0; el<offset; el++)
//       for(unsigned int se=0; se<8; se++)      {
//         int kn=mgmesh._elem_map[level][conn[se][0]+el*27];
//         double xp=mgmesh._xyz[kn]; double yp=mgmesh._xyz[kn+goffset];
//         double zp=mgmesh._xyz[kn+goffset2];
//         int  ix=(int)(xp/hx+0.5); int  iy=(int)(yp/hx+0.5);
//         int  iz=(int)(zp/hx+0.5);
//         int ind=ix+(iy+iz*nx)*nx;
// 
//         out <<  cc.Cmp[ind+1] << " ";
//       }
//     out <<" \n";
// 
//   }
// }
// // ------------------------------------------
// void MGSolCC::print(MGMeshC &mgmesh,
//                     const unsigned int flag_print,const unsigned int Level) {
// 
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   // label
//   sprintf(Named,"%d",flag_print);
//   // velocity and pressure
//   sprintf(outf2d,"./RESU/cc.%s.vtu",Named);  std::ofstream out(outf2d);
// 
//   unsigned int n_nodes=mgmesh._NoNodes[mgmesh._NoLevels -1]; unsigned int n_elements=mgmesh._NoElements[mgmesh._NoLevels -1];
//   out << "<VTKFile type=\"UnstructuredGrid\"  byte_type=\"LittleEndian\">\n";
//   out << "<UnstructuredGrid>\n";
//   out << "<Piece NumberOfPoints=\""<<n_nodes << "\" NumberOfCells=\"" << n_elements*NDOF_P << "\">\n";
//   // write the nodes --------------------------------------
//   out << "<Points> \n";
//   out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//   mgmesh.print_nodes(out,mgmesh._NoLevels -1,0);
//   out << "</DataArray>\n";
//   out<<  "</Points>\n";
//   // write the connectivity
//   out << "<Cells>\n";
//   out << "<DataArray type=\"Int32\" Name=\"connectivity\"  format=\"ascii\">\n";
//   mgmesh.print_conn(out,mgmesh._NoLevels -1,NDOF_FEM);
//   out << "</DataArray>\n";
//   out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
//   int offsets=0;
//   for(unsigned int el=0; el<n_elements; el++) {
//     for(unsigned int se=0; se<NDOF_P; se++)      {
//       offsets +=NDOF_P;      out << offsets << " ";
//     }
//   }
//   out << "\n";
//   out << "</DataArray>\n";
//   out << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
//   for(unsigned int el=0; el<n_elements; el++) {
//     for(unsigned int se=0; se<NDOF_P; se++)   out << 9+3*(DIMENSION-2) << " ";
//   }
//   out << "\n";
//   out << "</DataArray>\n";
//   out << "</Cells>\n";
//   // color function
//   out << "<CellData Scalars=\"CC\"> \n";
//   out << "<DataArray type=\"Float64\" Name=\"CC\" NumberOfComponents=\"1\" format=\"ascii\"> \n";
//   print_cc(mgmesh,out,n_elements,0,1);
//   out <<  "\n</DataArray>\n";
//   out<< "</CellData>\n";
//   out << "</Piece>\n";
//   out << "</UnstructuredGrid>\n";
//   out << "</VTKFile>\n";
//   out.close();
// 
//   // printing interface reconstruction
//   print_fine(flag_print,Level);
//   return;
// }


// ------------------------------------------
void MGSolCC::print_fine(const unsigned int flag_print,const unsigned int Level) {

  char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
  // label 00*,0**,***
  sprintf(Named,"%d",flag_print);
  // velocity and pressure
  sprintf(outf2d,"./RESU/cf.%s.vtu",Named);  std::ofstream out(outf2d);

  const unsigned int nx =_nxyz[0][Level]-1; const unsigned int ny =_nxyz[1][Level]-1;
  const unsigned int nz =_nxyz[2][Level]-1;
  double hx=1./nx; double hy=1./ny; double hz=1./nz;
  int fc=1; for(int lk=0; lk<Level; lk++) fc *=2;
  int n_cell=0;
  for(int ik=1; ik<=ny*nz; ik++) n_cell += M__GetLen(&c1[Level],ik);

  out << "<VTKFile type=\"UnstructuredGrid\"  byte_type=\"LittleEndian\">\n";
  out << "<UnstructuredGrid>\n";
  out << "<Piece NumberOfPoints=\""<< 8*n_cell << "\" NumberOfCells=\"" << n_cell << "\">\n";
  // write the nodes --------------------------------------
  out << "<Points> \n";
  out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ny+1);
      for(unsigned int kx=0; kx<len; kx++) {
        int ix=M__GetPos(&c1[Level],iy+kz*ny+1,kx)-1;
        double  ix0= hx*ix; double  iy0= hy*iy;
        double  iz0= hz*kz;
        out << ix0 << " " << iy0 << " " << iz0 << " ";
        out << ix0+hx << " " << iy0 << " " << iz0 << " ";
        out << ix0+hx << " " << iy0+hy << " " << iz0 << " ";
        out << ix0 << " " << iy0+hy << " " << iz0 << " ";
        out << ix0 << " " << iy0 << " " << iz0+hz << " ";
        out << ix0+hx << " " << iy0 << " " << iz0+hz << " ";
        out << ix0+hx << " " << iy0+hy << " " << iz0+hz << " ";
        out << ix0 << " " << iy0+hy << " " << iz0+hz << std::endl;
      }
    }
  out << "</DataArray>\n";
  out<<  "</Points>\n";
  // write the connectivity
  out << "<Cells>\n";
  out << "<DataArray type=\"Int32\" Name=\"connectivity\"  format=\"ascii\">\n";
  for(int  ik=0; ik<8*n_cell; ik++)  out << ik << " "; out << "\n";
  out << "</DataArray>\n";
  out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(int  ik=1; ik<=n_cell; ik++)  out << ik*8 << " ";
  out << "\n";
  out << "</DataArray>\n";
  out << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(int ik=1; ik<=n_cell; ik++)   out << 9+3*(DIMENSION-2) << " ";
  out << "\n";
  out << "</DataArray>\n";
  out << "</Cells>\n";

  // color function
  out << "<CellData Scalars=\"CC\" Vectors=\"Norm\"> \n";
  out << "<DataArray type=\"Float64\" Name=\"CC\" NumberOfComponents=\"1\" format=\"ascii\"> \n";
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ny+1);
      for(unsigned int kx=0; kx<len; kx++) {
        double val=M__GetVal(&c1[Level],iy+kz*ny+1,kx); if(val>1.) val=1.;
        out << val << " " ;
      }
    }
  out <<  "\n</DataArray>\n";
  // normal
  out << "<DataArray type=\"Float64\" Name=\"Norm\" NumberOfComponents=\"3\" format=\"ascii\"> \n";
  double mx,my,mz;
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ny+1);
      for(unsigned int kx=0; kx<len; kx++) {
        double val=M__GetVal(&c1[Level],iy+kz*ny+1,kx);
        mx=0.; my=0.; mz=0.;
        if(val<1.) {
          mx=M__GetVal(&_mx1[Level],iy+kz*ny+1,kx);
          my=M__GetVal(&_my1[Level],iy+kz*ny+1,kx);
          mz=M__GetVal(&_mz1[Level],iy+kz*ny+1,kx);
        }
        out << mx << " " << my << " "<< mz << " ";
      }
    }
  out <<  "\n</DataArray>\n";
  out<< "</CellData>\n";
  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  // print for restarting ------------------------
  out << "\n<!-- RestartData \n";
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ny+1);
      out << len;
      for(unsigned int kx=0; kx<len; kx++) {
        int pos=M__GetPos(&c1[Level],iy+kz*ny+1,kx);
        out << " " << pos;
        double val=M__GetVal(&c1[Level],iy+kz*ny+1,kx);
        out << " " << val;
      }
      out << "\n";
    }
  out << "-->\n";
  // ------------------------------------------

  out.close();

  return;
}
// ------------------
// ------------------------------------------
void MGSolCC::read_fine(const unsigned int flag_print,const unsigned int Level) {

  const unsigned int ny =_nxyz[1][Level]-1; const unsigned int nz =_nxyz[2][Level]-1;
  char *buf; buf=new char[50];
  char *Named; Named=new char[30]; char *inf2d; inf2d=new char[30];
  sprintf(Named,"%d",flag_print);
  // velocity and pressure
  sprintf(inf2d,"./RESU/cf.%s.vtu",Named);
  std::ifstream in(inf2d);
  if(!in) {printf("input file cf.%d.vtu not found\n",flag_print); exit(3);}
  while(strncmp(buf,"RestartData",11) != 0) in >> buf;
  int len,pos; double value;
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      in >> len; M_SetLen(&c1[Level],iy+kz*ny+1,len);
      for(unsigned int i=0; i<len; i++) {
        in >> pos;
        in >> value;
        M_SetEntry(&c1[Level],iy+kz*ny+1,i,pos,value);
      }
    }
  in.close();
  // Restriction --------------------------------
  OldSol_update(Level,c1[Level],c1_old[Level]);
  GenNorm(Level,&MGSolCC::rec_Young);
  for(int level=Level-1; level>=0; level--) {
    RestSol(level); OldSol_update(level,c1[level],c1_old[level]);
    RestNormal(level);
  }// ----------------------------------------
  // Projection --------------------------
  for(int level=Level+1; level<_Nlev_cc; level++) {
    ProjSol(level); OldSol_update(level,c1[level],c1_old[level]);
  } // -----------------------------------
  // Update ---------------------------
  CCSol_update();
  Normal_update();
  delete []buf; delete []Named; delete []inf2d;
  return;
}
// ------------------
// ------------------------------------------------------
// Utility for contracting and expanding compressed Matrix
// ------------------------------------------------------
// ----------------------------------
/// Row expansion  c1_old -> new line _tempc1
void MGSolCC::ExpRow(Matrix * mtrx,double tempc1[],unsigned int Level,unsigned int jy,unsigned int rowl) {
  int nxl=_nxyz[0][Level]; int ny=_nxyz[1][Level]-1; int nz=_nxyz[2][Level]-1;
  // Expanding matrix c1 -> new line _tempc1
  for(unsigned int ix=1; ix<nxl; ix++) tempc1[ix+rowl]=0.;
  if(jy%ny < ny && jy/ny <nz && jy%ny >=0 && jy/ny >=0) {  // no ghost cell
    for(unsigned int kx=0; kx<M__GetLen(mtrx,jy+1); kx++)  {
      int ix=M__GetPos(mtrx,jy+1,kx);
      double val=M__GetVal(mtrx,jy+1,kx);
      if(val <1.)  tempc1[ix+rowl]=val;
      else for(int ki=0; ki< (int) val; ki++) tempc1[ix+rowl+ki]=1.;
    }
  }
  return;
}
// -------------------------------------------
/// Row contraction  line _tempc2 -> c1 matrix
void MGSolCC::CtrRow(double tempc2[],Matrix * mtrx,unsigned int Level,unsigned int jy,unsigned int rowl) {
  int nxl=_nxyz[0][Level];
  int icount=0; int icount1=0;  int flag1=0;
  for(unsigned int ixf=1; ixf<nxl; ixf++) {
    int ind=ixf+rowl;
    if(tempc2[ind] > 1.-1.e-12) {
      if(flag1 == 0) icount1++;
      flag1 = 1; tempc2[ind]=1.;
    } else {
      flag1 = 0; if(tempc2[ind]> 1.e-12) icount++;
      else tempc2[ind]=0.;
    }
  }
  // compressed storage matrix c1
  M_SetLen(mtrx,jy+1,icount+icount1);
  int irow1=0; int count1=0;
  for(unsigned int ixf=1; ixf<nxl; ixf++) {
    int ind=ixf+rowl;
    if(tempc2[ind]> 0.) {
      if(tempc2[ind]< 1.) { // _tempc2[ind] < 1.
        M__SetEntry(mtrx,jy+1,irow1,ixf,tempc2[ind]);
        irow1++;
      } // <1
      else {// _tempc2[ind]== 1.
        count1++;
        if(tempc2[ind+1]<1.) {
          M__SetEntry(mtrx,jy+1,irow1,M__GetPos(mtrx,jy+1,irow1-1) +1,count1);
          irow1++;
          count1=0;
        }
      } // else
      tempc2[ind]=0.;
    } // .>0
  }
  return;
}

// --------------------------------
// Cartesian cell evaluation
// ----------------------------------------
/// Return cell velocity field
void MGSolCC::get_vel(int ixyz[],int fd,double u[]) {
  // geometry level
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1; int nzc=_nxyz[2][0]-1;
  int n_nodes= (nxc+1)*(nyc+1)*(nzc+1); int n_nodes2= 2*n_nodes;
  // velocity field  ul(left),ur(right),vt(top),vb(bot)
  int ixc=ixyz[0]/fd; int ixf=ixyz[0]%fd;
  int jyc=ixyz[1]/fd; int jyf=ixyz[1]%fd;
  int kzc=ixyz[2]/fd; int kzf=ixyz[2]%fd;
  int indv = ixc+(jyc+kzc*(nyc+1))*(nxc+1);
  const unsigned int indy = indv+nxc+1;
  const unsigned int indx = indv+1;
  const unsigned int indz = indv+(nxc+1)*(nyc+1);
  // u-component
  double u0 =_uvw.Cmp[_invnode_dof[indv]]; // v-right
  double u1 =_uvw.Cmp[_invnode_dof[indv+1]]; // v-right
  double u2 =_uvw.Cmp[_invnode_dof[indv+1+(nxc+1)]]; // v-right
  double u3 =_uvw.Cmp[_invnode_dof[indv+(nxc+1)]]; // v-right
  double u4 =_uvw.Cmp[_invnode_dof[indz]]; // v-right
  double u5 =_uvw.Cmp[_invnode_dof[indz+1]]; // v-right
  double u6 =_uvw.Cmp[_invnode_dof[indz+1+(nxc+1)]]; // v-right
  double u7 =_uvw.Cmp[_invnode_dof[indz+(nxc+1)]]; // v-right

  // v-component
  double v0 =_uvw.Cmp[_invnode_dof[indv]+n_nodes]; // v-right
  double v1 =_uvw.Cmp[_invnode_dof[indv+1]+n_nodes]; // v-right
  double v2 =_uvw.Cmp[_invnode_dof[indv+1+(nxc+1)]+n_nodes]; // v-right
  double v3 =_uvw.Cmp[_invnode_dof[indv+(nxc+1)]+n_nodes]; // v-right
  double v4 =_uvw.Cmp[_invnode_dof[indz]+n_nodes]; // v-right
  double v5 =_uvw.Cmp[_invnode_dof[indz+1]+n_nodes]; // v-right
  double v6 =_uvw.Cmp[_invnode_dof[indz+1+(nxc+1)]+n_nodes]; // v-right
  double v7 =_uvw.Cmp[_invnode_dof[indz+(nxc+1)]+n_nodes]; // v-right

  // w-component
  double w0 =_uvw.Cmp[_invnode_dof[indv]+n_nodes2]; // v-right
  double w1 =_uvw.Cmp[_invnode_dof[indv+1]+n_nodes2]; // v-right
  double w2 =_uvw.Cmp[_invnode_dof[indv+1+(nxc+1)]+n_nodes2]; // v-right
  double w3 =_uvw.Cmp[_invnode_dof[indv+(nxc+1)]+n_nodes2]; // v-right
  double w4 =_uvw.Cmp[_invnode_dof[indz]+n_nodes2]; // v-right
  double w5 =_uvw.Cmp[_invnode_dof[indz+1]+n_nodes2]; // v-right
  double w6 =_uvw.Cmp[_invnode_dof[indz+1+(nxc+1)]+n_nodes2]; // v-right
  double w7 =_uvw.Cmp[_invnode_dof[indz+(nxc+1)]+n_nodes2]; // v-right

  u[0] =(u0*(fd-ixf)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u1*(ixf)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u2*(ixf)*(jyf+0.5)*(fd-kzf-0.5)+
         u3*(fd-ixf)*(jyf+0.5)*(fd-kzf-0.5)+
         u4*(fd-ixf)*(fd-jyf-0.5)*(kzf+0.5)+
         u5*(ixf)*(fd-jyf-0.5)*(kzf+0.5)+
         u6*(ixf)*(jyf+0.5)*(kzf+0.5)+
         u7*(fd-ixf)*(jyf+0.5)*(kzf+0.5)
        )/(fd*fd*fd);
  u[1] =(u0*(fd-ixf-1)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u1*(ixf+1)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u2*(ixf+1)*(jyf+0.5)*(fd-kzf-0.5)+
         u3*(fd-ixf-1)*(jyf+0.5)*(fd-kzf-0.5)+
         u4*(fd-ixf-1)*(fd-jyf-0.5)*(kzf+0.5)+
         u5*(ixf+1)*(fd-jyf-0.5)*(kzf+0.5)+
         u6*(ixf+1)*(jyf+0.5)*(kzf+0.5)+
         u7*(fd-ixf-1)*(jyf+0.5)*(kzf+0.5)
        )/(fd*fd*fd);
  u[2] =(v0*(fd-ixf-0.5)*(fd-jyf)*(fd-kzf-0.5)+
         v1*(ixf+0.5)*(fd-jyf)*(fd-kzf-0.5)+
         v2*(ixf+0.5)*(jyf)*(fd-kzf-0.5)+
         v3*(fd-ixf-0.5)*(jyf)*(fd-kzf-0.5)+
         v4*(fd-ixf-0.5)*(fd-jyf)*(kzf+0.5)+
         v5*(ixf+0.5)*(fd-jyf)*(kzf+0.5)+
         v6*(ixf+0.5)*(jyf)*(kzf+0.5)+
         v7*(fd-ixf-0.5)*(jyf)*(kzf+0.5)
        )/(fd*fd*fd);
  u[3] =(v0*(fd-ixf-0.5)*(fd-jyf-1.)*(fd-kzf-0.5)+
         v1*(ixf+0.5)*(fd-jyf-1.)*(fd-kzf-0.5)+
         v2*(ixf+0.5)*(jyf+1.)*(fd-kzf-0.5)+
         v3*(fd-ixf-0.5)*(jyf+1.)*(fd-kzf-0.5)+
         v4*(fd-ixf-0.5)*(fd-jyf-1.)*(kzf+0.5)+
         v5*(ixf+0.5)*(fd-jyf-1.)*(kzf+0.5)+
         v6*(ixf+0.5)*(jyf+1.)*(kzf+0.5)+
         v7*(fd-ixf-0.5)*(jyf+1.)*(kzf+0.5)
        )/(fd*fd*fd);

  u[4] =(w0*(fd-ixf-0.5)*(fd-jyf-0.5)*(fd-kzf)+
         w1*(ixf+0.5)*(fd-jyf-0.5)*(fd-kzf)+
         w2*(ixf+0.5)*(jyf+0.5)*(fd-kzf)+
         w3*(fd-ixf-0.5)*(jyf+0.5)*(fd-kzf)+
         w4*(fd-ixf-0.5)*(fd-jyf-0.5)*(kzf)+
         w5*(ixf+0.5)*(fd-jyf-0.5)*(kzf)+
         w6*(ixf+0.5)*(jyf+0.5)*(kzf)+
         w7*(fd-ixf-0.5)*(jyf+0.5)*(kzf)
        )/(fd*fd*fd);
  u[5] =(w0*(fd-ixf-0.5)*(fd-jyf-0.5)*(fd-kzf-1.)+
         w1*(ixf+0.5)*(fd-jyf-0.5)*(fd-kzf-1.)+
         w2*(ixf+0.5)*(jyf+0.5)*(fd-kzf-1.)+
         w3*(fd-ixf-0.5)*(jyf+0.5)*(fd-kzf-1.)+
         w4*(fd-ixf-0.5)*(fd-jyf-0.5)*(kzf+1.)+
         w5*(ixf+0.5)*(fd-jyf-0.5)*(kzf+1.)+
         w6*(ixf+0.5)*(jyf+0.5)*(kzf+1.)+
         w7*(fd-ixf-0.5)*(jyf+0.5)*(kzf+1.)
        )/(fd*fd*fd);

  // Analytic vortex test
  //   double ix3=ixyz[0]*hx;double iy3=ixyz[1]*hy;
  //   double ul =-(PSI(ix3,(iy3+hy))-PSI(ix3,iy3))/hy;
  //         double ur =-(PSI(ix3+hx,iy3+hy)-PSI(ix3+hx,iy3))/hy;
  //         double vb = (PSI(ix3+hx,iy3)-PSI(ix3,iy3))/hx;
  //         double vt = (PSI(ix3+hx,iy3+hy)-PSI(ix3,iy3+hy))/hx;
  return;
}
// --------------------------------
// Cartesian cell evaluation
// ----------------------------------------
/// Return cell velocity field
void MGSolCC::get_disp(int ixyz[],int fd,double u[]) {
  // geometry level
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1; int nzc=_nxyz[2][0]-1;
  int n_nodes= (nxc+1)*(nyc+1)*(nzc+1); int n_nodes2= 2*n_nodes;
  // velocity field  ul(left),ur(right),vt(top),vb(bot)
  int ixc=ixyz[0]/fd; int ixf=ixyz[0]%fd;
  int jyc=ixyz[1]/fd; int jyf=ixyz[1]%fd;
  int kzc=ixyz[2]/fd; int kzf=ixyz[2]%fd;
  int indv = ixc+(jyc+kzc*(nyc+1))*(nxc+1);
  const unsigned int indy = indv+nxc+1;
  const unsigned int indx = indv+1;
  const unsigned int indz = indv+(nxc+1)*(nyc+1);
  // u-component
  double u0 =_dsvw.Cmp[_invnode_dof[indv]]; // v-right
  double u1 =_dsvw.Cmp[_invnode_dof[indv+1]]; // v-right
  double u2 =_dsvw.Cmp[_invnode_dof[indv+1+(nxc+1)]]; // v-right
  double u3 =_dsvw.Cmp[_invnode_dof[indv+(nxc+1)]]; // v-right
  double u4 =_dsvw.Cmp[_invnode_dof[indz]]; // v-right
  double u5 =_dsvw.Cmp[_invnode_dof[indz+1]]; // v-right
  double u6 =_dsvw.Cmp[_invnode_dof[indz+1+(nxc+1)]]; // v-right
  double u7 =_dsvw.Cmp[_invnode_dof[indz+(nxc+1)]]; // v-right

  // v-component
  double v0 =_dsvw.Cmp[_invnode_dof[indv]+n_nodes]; // v-right
  double v1 =_dsvw.Cmp[_invnode_dof[indv+1]+n_nodes]; // v-right
  double v2 =_dsvw.Cmp[_invnode_dof[indv+1+(nxc+1)]+n_nodes]; // v-right
  double v3 =_dsvw.Cmp[_invnode_dof[indv+(nxc+1)]+n_nodes]; // v-right
  double v4 =_dsvw.Cmp[_invnode_dof[indz]+n_nodes]; // v-right
  double v5 =_dsvw.Cmp[_invnode_dof[indz+1]+n_nodes]; // v-right
  double v6 =_dsvw.Cmp[_invnode_dof[indz+1+(nxc+1)]+n_nodes]; // v-right
  double v7 =_dsvw.Cmp[_invnode_dof[indz+(nxc+1)]+n_nodes]; // v-right

  // w-component
  double w0 =_dsvw.Cmp[_invnode_dof[indv]+n_nodes2]; // v-right
  double w1 =_dsvw.Cmp[_invnode_dof[indv+1]+n_nodes2]; // v-right
  double w2 =_dsvw.Cmp[_invnode_dof[indv+1+(nxc+1)]+n_nodes2]; // v-right
  double w3 =_dsvw.Cmp[_invnode_dof[indv+(nxc+1)]+n_nodes2]; // v-right
  double w4 =_dsvw.Cmp[_invnode_dof[indz]+n_nodes2]; // v-right
  double w5 =_dsvw.Cmp[_invnode_dof[indz+1]+n_nodes2]; // v-right
  double w6 =_dsvw.Cmp[_invnode_dof[indz+1+(nxc+1)]+n_nodes2]; // v-right
  double w7 =_dsvw.Cmp[_invnode_dof[indz+(nxc+1)]+n_nodes2]; // v-right

  u[0] =(u0*(fd-ixf)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u1*(ixf)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u2*(ixf)*(jyf+0.5)*(fd-kzf-0.5)+
         u3*(fd-ixf)*(jyf+0.5)*(fd-kzf-0.5)+
         u4*(fd-ixf)*(fd-jyf-0.5)*(kzf+0.5)+
         u5*(ixf)*(fd-jyf-0.5)*(kzf+0.5)+
         u6*(ixf)*(jyf+0.5)*(kzf+0.5)+
         u7*(fd-ixf)*(jyf+0.5)*(kzf+0.5)
        )/(fd*fd*fd);
  u[1] =(u0*(fd-ixf-1)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u1*(ixf+1)*(fd-jyf-0.5)*(fd-kzf-0.5)+
         u2*(ixf+1)*(jyf+0.5)*(fd-kzf-0.5)+
         u3*(fd-ixf-1)*(jyf+0.5)*(fd-kzf-0.5)+
         u4*(fd-ixf-1)*(fd-jyf-0.5)*(kzf+0.5)+
         u5*(ixf+1)*(fd-jyf-0.5)*(kzf+0.5)+
         u6*(ixf+1)*(jyf+0.5)*(kzf+0.5)+
         u7*(fd-ixf-1)*(jyf+0.5)*(kzf+0.5)
        )/(fd*fd*fd);
  u[2] =(v0*(fd-ixf-0.5)*(fd-jyf)*(fd-kzf-0.5)+
         v1*(ixf+0.5)*(fd-jyf)*(fd-kzf-0.5)+
         v2*(ixf+0.5)*(jyf)*(fd-kzf-0.5)+
         v3*(fd-ixf-0.5)*(jyf)*(fd-kzf-0.5)+
         v4*(fd-ixf-0.5)*(fd-jyf)*(kzf+0.5)+
         v5*(ixf+0.5)*(fd-jyf)*(kzf+0.5)+
         v6*(ixf+0.5)*(jyf)*(kzf+0.5)+
         v7*(fd-ixf-0.5)*(jyf)*(kzf+0.5)
        )/(fd*fd*fd);
  u[3] =(v0*(fd-ixf-0.5)*(fd-jyf-1.)*(fd-kzf-0.5)+
         v1*(ixf+0.5)*(fd-jyf-1.)*(fd-kzf-0.5)+
         v2*(ixf+0.5)*(jyf+1.)*(fd-kzf-0.5)+
         v3*(fd-ixf-0.5)*(jyf+1.)*(fd-kzf-0.5)+
         v4*(fd-ixf-0.5)*(fd-jyf-1.)*(kzf+0.5)+
         v5*(ixf+0.5)*(fd-jyf-1.)*(kzf+0.5)+
         v6*(ixf+0.5)*(jyf+1.)*(kzf+0.5)+
         v7*(fd-ixf-0.5)*(jyf+1.)*(kzf+0.5)
        )/(fd*fd*fd);

  u[4] =(w0*(fd-ixf-0.5)*(fd-jyf-0.5)*(fd-kzf)+
         w1*(ixf+0.5)*(fd-jyf-0.5)*(fd-kzf)+
         w2*(ixf+0.5)*(jyf+0.5)*(fd-kzf)+
         w3*(fd-ixf-0.5)*(jyf+0.5)*(fd-kzf)+
         w4*(fd-ixf-0.5)*(fd-jyf-0.5)*(kzf)+
         w5*(ixf+0.5)*(fd-jyf-0.5)*(kzf)+
         w6*(ixf+0.5)*(jyf+0.5)*(kzf)+
         w7*(fd-ixf-0.5)*(jyf+0.5)*(kzf)
        )/(fd*fd*fd);
  u[5] =(w0*(fd-ixf-0.5)*(fd-jyf-0.5)*(fd-kzf-1.)+
         w1*(ixf+0.5)*(fd-jyf-0.5)*(fd-kzf-1.)+
         w2*(ixf+0.5)*(jyf+0.5)*(fd-kzf-1.)+
         w3*(fd-ixf-0.5)*(jyf+0.5)*(fd-kzf-1.)+
         w4*(fd-ixf-0.5)*(fd-jyf-0.5)*(kzf+1.)+
         w5*(ixf+0.5)*(fd-jyf-0.5)*(kzf+1.)+
         w6*(ixf+0.5)*(jyf+0.5)*(kzf+1.)+
         w7*(fd-ixf-0.5)*(jyf+0.5)*(kzf+1.)
        )/(fd*fd*fd);

  // Analytic vortex test
  //   double ix3=ixyz[0]*hx;double iy3=ixyz[1]*hy;
  //   double ul =-(PSI(ix3,(iy3+hy))-PSI(ix3,iy3))/hy;
  //         double ur =-(PSI(ix3+hx,iy3+hy)-PSI(ix3+hx,iy3))/hy;
  //         double vb = (PSI(ix3+hx,iy3)-PSI(ix3,iy3))/hx;
  //         double vt = (PSI(ix3+hx,iy3+hy)-PSI(ix3,iy3+hy))/hx;
  return;
}
// -------------------------------------------------------------------
// Volume
// -------------------------------------------------------------------
// -----------------------------------------------------------------
double MGSolCC::get_area(double mx,double my,double alpha,double x0,double y0) {
  double area,tmp;
  /* move origin to (x0,y0) */
  alpha = alpha - my*y0 - mx*x0;
  /* rotate figure in the proper way */
  alpha = alpha + MAX(0.,-mx) + MAX(0.,-my);
  mx = fabs(mx);  my = fabs(my);
  /* limit alpha so that 0. <= area <= 1. */
  alpha = MAX(0.,MIN(alpha,mx+my));
  if(MIN(mx,my) < MIN_VAL)   area = alpha;
  else    area = 0.5/(mx*my)*(alpha*alpha - MAX(0.0,alpha-mx)*
                                (alpha-mx) - MAX(0.0,alpha-my)*(alpha-my));
  return area;
}


// -------------------------------------------------------------------
// ALPHA reconstruction
// -------------------------------------------------------------------
double MGSolCC::get_alpha(double m1, double m2,  double m3, double v) {

  double v1,v2,v3,alpha,m,m12;
  double a3,a2,a1,a0;
  double teta,p,q,t;
  int k=1,l=1;

  double temp;
  if(m2> m3) {temp= m3; m3= m2; m2=temp;}
  if(m1> m2) {temp= m2; m2= m1; m1=temp;}
  if(m2> m3) {temp= m3; m3= m2; m2=temp;}

  if(m2==0) alpha=v;
  else if(m1==0) {
    m=m2;     v1=m/(2*(1-m));
    if(0.<v && v<v1) alpha=sqrt(2*m*(1-m)*v);
    else if(v1<=v && v<=1-v1) alpha=v*(1-m)+m/2;
    else if(1-v1<=v && v<=1.) alpha=1-sqrt(2*m*(1-m)*(1-v));
  } else {
    m12=m1+m2;
    if(v>0.5) { v=1.-v; k=2; }
    if(m12<m3) l=2;

    v1=m1*m1/(6*m2*m3);    v2=v1+(m2-m1)/(2*m3);
    v3=(m3*m3*(3*m12-m3)+m1*m1*(m1-3*m3)+m2*m2*(m2-3*m3))/(6*m1*m2*m3);
    if(l==2) v3=m12/(2*m3);

    if(0<=v && v<v1) {
      t=6*m1*m2*m3*v; alpha=pow(t,(1./3));
    } else if(v1<=v && v<v2) alpha=0.5*(m1+sqrt(m1*m1+8*m2*m3*(v-v1)));
    else if(v2<=v && v<v3) {
      a3=-1.;       a2=3*m12/a3;
      a1=-3*(m1*m1+m2*m2)/a3;  a0=(m1*m1*m1+m2*m2*m2-6*m1*m2*m3*v)/a3;

      p=a2*a2/9.-a1/3.; q=(a1*a2-3.*a0)/6.-a2*a2*a2/27.;
      teta=acos(q/sqrt(p*p*p))/3.;
      alpha=sqrt(p)*(sqrt(3.)*sin(teta)-cos(teta))-a2/3.;
    } else if(v>=v3) {
      if(l==1) {
        a3=-2.;  a0=(m1*m1*m1+m2*m2*m2+m3*m3*m3-6*m1*m2*m3*v)/a3;
        a2=3/a3; a1=-3*(m1*m1+m2*m2+m3*m3)/a3;
        p=a2*a2/9.-a1/3.; q=(a1*a2-3*a0)/6.-a2*a2*a2/27.;
        teta=acos(q/sqrt(p*p*p))/3.;
        alpha=sqrt(p)*(sqrt(3.)*sin(teta)-cos(teta))-a2/3.;
      } else if(l==2) alpha=m3*v+m12/2;
    }
    if(k==2) alpha=1-alpha;
  }
  return alpha;
}




/// return cell phase
double MGSolCC::get_2phase(int blck,double xp[]) {
// Set up
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1; int nzc=_nxyz[2][0]-1;
  int  ix=(int)(xp[0]*nxc+0.5); int  iy=(int)(xp[1]*nyc+0.5);
  int  iz=(int)(xp[2]*nzc+0.5);
  const int ind=ix+(iy+iz*(nyc+1))*(nxc+1);
  double phase=0.;
  // average over neighbouring cells
  for(int kk=-blck; kk<blck; kk++)
    for(int jj=-blck; jj<blck; jj++)
      for(int ii=-blck; ii<blck; ii++)
        phase +=cc.Cmp[ind+(jj+kk*(nyc+1))*(nxc+1)+ii+1];
  return phase/(8.*blck*blck*blck);
}




// ----------------------------------------------------
double shape2(int i,double xi) {
  switch(i) {
  case 0:  return .5*xi*(xi-1.);
  case 2:  return .5*xi*(xi+1.);
  case 1:  return (1.-xi*xi);
  }
  return 0.;
}
// ----------------------------------------------------
/// Return surface force tension
void MGSolCC::get_2surten(double xp[],double ff[],int ord[]) {

  for(int kf=0; kf<100; kf++)  ff[kf]=0.;
  int nxf=_nxyz[0][_Nlev_cc-1]-1;
  int nyf=_nxyz[1][_Nlev_cc-1]-1; int nzf=_nxyz[2][_Nlev_cc-1]-1;
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1; int nzc=_nxyz[2][0]-1;
  double hxc=1./nxc;
  double hxf=1./nxf; double hyf=1./nyf; double hzf=1./nzf;
  int fd=nxf/nxc;
  double pt1[3],pt2[3],vl[3],vt[3];
  double m1,m2,ma;
  double ww[3]= {0.555555555556,0.888888888889,0.555555555556};
  double xgw[3]= {-0.774596669241,0.,0.774596669241};
  double nns[3][3];
  int  ixc=(int)(xp[0]*nxc+0.5); int  jyc=(int)(xp[1]*nyc+0.5);
  int  kzc=(int)(xp[2]*nzc+0.5);
  const int indc=ixc+(jyc+kzc*(nyc+1))*(nxc+1);




  for(int kk=-fd; kk<fd; kk++)
    for(int jj=-fd; jj<fd; jj++) {
      int indf=(jyc*fd+jj)+(kzc*fd+kk)*nyf;

      ExpRow(&c1[_Nlev_cc-1],_tempc1,_Nlev_cc-1,indf,0);
      ExpRow(&_mx1[_Nlev_cc-1],_tempmx1,_Nlev_cc-1,indf,0);
      ExpRow(&_my1[_Nlev_cc-1],_tempmy1,_Nlev_cc-1,indf,0);
      ExpRow(&_mz1[_Nlev_cc-1],_tempmz1,_Nlev_cc-1,indf,0);

      for(int ii=-fd; ii<fd; ii++) {
        int ixf=(ixc*fd+ii)+1;   double ccc=_tempc1[ixf];
        if(ccc >0. && ccc<1.) {
          double mx=_tempmx1[ixf]; double my=_tempmy1[ixf];
          double mz=_tempmz1[ixf];

          double  alpha=get_alpha(fabs(mx),fabs(my),fabs(mz),ccc);
          alpha += MIN(mx,0.)+MIN(my,0.)+MIN(mz,0.);

          // double alpha=get_alpha(mx,my,mz,ccc);

          for(int iface=0; iface<6; iface++) {
            int dir=iface%3; int pos=iface/3;
            double length=0.;
            // direction  x
            if(dir == 0) {
              ma=1./(fabs(my)+fabs(mz)+1.e-6);
              m1=my*ma; m2=mz*ma;
              ma=(alpha-mx*pos)*ma;
              get_pts(m1,m2,ma,pt1,pt2);
              length=hxf*sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+
                              (pt2[1]-pt1[1])*(pt2[1]-pt1[1]));
              // interface reconstruction points
              pt1[2] = (pt1[1]+kk)/(double)fd;
              pt2[2] = (pt2[1]+kk)/(double)fd;
              pt1[1] = (pt1[0]+jj)/(double)fd;
              pt2[1] = (pt2[0]+jj)/(double)fd;
              pt1[0] = (pos+ii)/(double)fd;
              pt2[0] = (pos+ii)/(double)fd;
            }
            // direction  y
            if(dir == 1) {
              ma=1./(fabs(mz)+fabs(mx)+1.e-6);
              m1=mz*ma; m2=mx*ma;
              ma=(alpha-my*pos)*ma;
              get_pts(m1,m2,ma,pt1,pt2);
              length=hyf*sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+
                              (pt2[1]-pt1[1])*(pt2[1]-pt1[1]));
              // interface reconstruction points
              pt1[2] = (pt1[0]+kk)/(double)fd;
              pt2[2] = (pt2[0]+kk)/(double)fd;
              pt1[0] = (pt1[1]+ii)/(double)fd;
              pt2[0] = (pt2[1]+ii)/(double)fd;
              pt1[1] = (pos+jj)/(double)fd;
              pt2[1] = (pos+jj)/(double)fd;
            }
            // direction  z
            if(dir == 2) {
              ma=1./(fabs(mx)+fabs(my)+1.e-6);
              m1=mx*ma; m2=my*ma;
              ma=(alpha-mz*pos)*ma;
              get_pts(m1,m2,ma,pt1,pt2);
              length=hzf*sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+
                              (pt2[1]-pt1[1])*(pt2[1]-pt1[1]));
              // interface reconstruction points
              pt1[2] = (pos+kk)/(double)fd;
              pt2[2] = (pos+kk)/(double)fd;
              pt1[1] = (pt1[1]+jj)/(double)fd;
              pt2[1] = (pt2[1]+jj)/(double)fd;
              pt1[0] = (pt1[0]+ii)/(double)fd;
              pt2[0] = (pt2[0]+ii)/(double)fd;
            } // dir =

            if(length >0) {


              vl[0]=pt2[0]-pt1[0]; vl[1]=pt2[1]-pt1[1]; vl[2]=pt2[2]-pt1[2];
              vt[0]=my*vl[2]-vl[1]*mz;
              vt[1]=-(mx*vl[2]-vl[0]*mz);
              vt[2]=mx*vl[1]-vl[0]*my;
              double vtm=sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2])*(1.-2*pos);
              vt[0] /=vtm; vt[1] /=vtm; vt[2] /=vtm;
              if(pos == 0 && vt[dir] <0) printf(" error dir ");
              if(pos == 1 && vt[dir] >0) printf(" error dir ");

              // printf(" vl %d %d %d  %d %d %e %e %e %e \n",ixf,jyc*fd+jj,kzc*fd+kk,dir,pos,vt[0],vt[1],vt[2],length);
              vtm=sqrt(mx*mx+my*my+mz*mz);
              double mx3 =mx/vtm;  double my3 =my/vtm;  double mz3 =mz/vtm;
              nns[0][0]=1.-mx3*mx3; nns[0][1]=-mx3*my3; nns[0][2]=-mx3*mz3;
              nns[1][0]=-mx3*my3; nns[1][1]=1.-my3*my3; nns[1][2]=-my3*mz3;
              nns[2][0]=-mx3*mz3; nns[2][1]=-mz3*my3; nns[2][2]=1-mz3*mz3;
              for(int ik=0; ik<3; ik++)
                vl[ik]= nns[ik][0]*vt[0]+nns[ik][1]*vt[1]+nns[ik][2]*vt[2];
              for(unsigned int ks=0; ks<3; ks++) {
                for(unsigned int js=0; js<3; js++) {
                  for(unsigned int is=0; is<3; is++) {
                    int inds=is+(js+ks*3)*3;
                    double integral=0.;
                    // gauss integration
                    for(int ig=0; ig<3; ig++) {
                      double xg= 0.5*(pt2[0]*(1+xgw[ig])+pt1[0]*(1-xgw[ig]));
                      double yg= 0.5*(pt2[1]*(1+xgw[ig])+pt1[1]*(1-xgw[ig]));
                      double zg= 0.5*(pt2[2]*(1+xgw[ig])+pt1[2]*(1-xgw[ig]));
                      integral +=ww[ig]*shape2(is,xg)*shape2(js,yg)*shape2(ks,zg);
                    } // igauss
                    //integral=fabs(integral);
                    ff[inds]    += 0.5*integral*length*vl[0];
                    ff[inds+27] += 0.5*integral*length*vl[1];
                    ff[inds+54] += 0.5*integral*length*vl[2];



                  }
                }
              }
            } // length

          } // iface


          // }
        } // if ccc
      } //  ii ->  fd
    }

// cartesian ordering
  for(int ksc=-1; ksc<2; ksc++)
    for(int jsc=-1; jsc<2; jsc++)
      for(int isc=-1; isc<2; isc++)
        ord[isc+1+(jsc+1+(ksc+1)*3)*3]=
          _invnode_dof[indc+isc+(jsc+ksc*(nyc+1))*(nxc+1)];

  return;

}

// -------------------------------------------------------------------
// Return interface points
void MGSolCC::get_pts(double mx, double my, double alpha,double *pt1, double *pt2) {
  int invx,invy; double tmp;

  /* rotate figure in the proper way */
  alpha = alpha + MAX(0.,-mx) + MAX(0.,-my);
  invx = invy = 0;
  if(mx < 0.) {invx = 1; mx = fabs(mx); }
  mx = mx + MIN_VAL;
  if(my < 0.) {invy = 1; my = fabs(my); }
  my = my + MIN_VAL;

  /* find two intersections of: mx*x + my*y = alpha with unit square */
  pt2[0] = 0.; pt2[1] = alpha/my;
  if(pt2[1] > 1.) {
    pt2[1] = 1.0; pt2[0] = (alpha-my) /mx;
  }
  pt1[0] = 1.; pt1[1] = (alpha-mx) /my;
  if(pt1[1] < 0.) {
    pt1[1] = 0.;  pt1[0] = alpha/mx;
  }
  /* if necessary rotate around x=0.5 */
  if(invx) {
    tmp = pt2[1]; pt2[1] = pt1[1]; pt1[1] = tmp;
    tmp = pt2[0]; pt2[0] = 1. - pt1[0]; pt1[0] = 1. - tmp;
  }
  /* if necessary rotate around y=0.5 */
  if(invy) {
    tmp = pt2[1]; pt2[1] = 1. - pt1[1]; pt1[1] = 1. - tmp;
    tmp = pt2[0]; pt2[0] = pt1[0]; pt1[0] = tmp;
  }
  if(pt1[0] <0. || pt1[0] >1.) {pt1[0]=0.; pt1[1]=0.; pt2[0]=0.; pt2[1]=0.;}
  if(pt1[1] <0. || pt1[1] >1.) {pt1[0]=0.; pt1[1]=0.; pt2[0]=0.; pt2[1]=0.;}
  if(pt2[0] <0. || pt2[0] >1.) {pt1[0]=0.; pt1[1]=0.; pt2[0]=0.; pt2[1]=0.;}
  if(pt2[1] <0. || pt2[1] >1.) {pt1[0]=0.; pt1[1]=0.; pt2[0]=0.; pt2[1]=0.;}
  return;
}
// void MGSolCC::setFieldSource(
//   int interface_name,
//   int n_cmp,
//   const ParaMEDMEM::MEDCouplingFieldDouble * srcField,
//    InterfaceFunctionDD * fct 
//   
//                             ) {}

// // -------------------------------------------------------
// --------------------------------------------------------

// void MGSolCC::InitVel(
//   const MGMeshC & mgmesh,
// //   const MGSol & mgs,
//   const double dt) {
//  V_Destr(&uvw);
// // #ifdef NS_EQUATIONS
// //   uvw=mgs.uw[mgs._NoLevels-1];
// // #else
//   V_Constr(&uvw,(char *)"uvw",3*_nxyz[0][0]*_nxyz[1][0]*_nxyz[2][0],Normal,_LPTrue);
//   GenVel(mgmesh,uvw,dt);
// // #endif
//   return;
// }
// // ============================================================================
// void MGSolCC::GenVel(
// //   const MGMeshC & mgmesh,
// // QVector &sol,
//   const int Level, ///< Level of the velocity field
//   const double dt  ///< time step
// ) { // ==========================================================================
// 
// // Get a constant reference to the MG solution.
// //   int Level=mgmesh._NoLevels -1;
//   int nxyz[3]; nxyz[0]=_nxyz[0][Level]; nxyz[1]=_nxyz[1][Level]; nxyz[2]=1;
//   double hxyz[3]; hxyz[0]=1./(nxyz[0]-1); hxyz[1]=1./(nxyz[1]-1); hxyz[2]=0;
// #if DIMENSION==3
//   nxyz[2]=_nxyz[2][Level]; hxyz[2]=1./(nxyz[2]-1);
// #endif
//   int  offset=1;
//   for(int idim=0; idim<DIMENSION; idim++) offset *=nxyz[idim];
// //   const  unsigned int  offset2= offset*2;
// 
//   V_Destr(&_uvw);
//   V_Constr(&_uvw,(char *)"uvw",DIMENSION*offset,Normal,_LPTrue);
// //   const unsigned int nx=_nxyz[0][0];const unsigned int ny=_nxyz[1][0];const unsigned int nz=_nxyz[2][0];
// //   const unsigned int nxf=_nxyz[0][_Nlev_cc-1];
// //   const double hx=1./(nx-1);const double hy=1./(ny-1);const double hz=1./(nz-1);
// //   const double hxf=1./(nxf-1);
// 
//   // Get a constant reference to the mesh objects.
// //   const int *map_nodes=mgmesh._elem_map[Level];
// //   const double *xyz_glob=mgmesh._xyz;
// //   const  unsigned int  offset=mgmesh._NoNodes[Level];
// 
// //   const unsigned int  n_elem=mgmesh._NoElements[Level];
// 
//   // Dof
//   //int *idx_dof=_node_dof[Level];
//   double u_value[3]; int dof_u[3];//double Real v_value =0.;
//   for(int idim=0; idim<DIMENSION; idim++) {
//     u_value[idim]=0.;  dof_u[idim]=0;
//   }
// //   const unsigned int n_u_dofs = NDOF_FEM; const unsigned int n_p_dofs = NDOF_P;
//   for(unsigned int ix=0 ; ix <nxyz[0]; ix++) {
//     for(unsigned int iy=0 ; iy <nxyz[1]; iy++) {
//       for(unsigned int iz=0 ; iz <nxyz[2]; iz++) {
//         double xi =ix*hxyz[0];
//         double yi =iy*hxyz[1];
//         double zi =iz*hxyz[2];
// //       const Real yj = xyz_glob[k+offset];
// //       const Real zk = xyz_glob[k+offset2];
//         // dofs
// //         unsigned int dof_u=ix+iy*nxyz[0]+iz*nxyz[0]*nxyz[1]+1;
// //         unsigned int dof_v=dof_u+offset;
// //         unsigned int dof_w=dof_u+offset2;
// //rotation --------------------------
// //         u_value[0]=/*/*3.14159265358979*/*/*(yj-0.5);
// //         u_value[1]=-3.14159265358979*(xi-0.5);
// //         u_value[2]=0.;
// // --------------------------------------------
// // // vortex 2D ---------------------------
// //       u_value[1] =  ((double) nx)*(PSI(xi+hx,yi)-PSI(xi,yi));
// //       u_value[0] = - ((double) nx)*(PSI(xi,yi+hx)-PSI(xi,yi));
// //       u_value[2] = 0;
// // ---------------------------------
// // vortex 3D ---------------------------
//         const double Pigreco=3.14159265358979;
//         u_value[0] =10* sin(Pigreco*xi)*sin(Pigreco*xi)*(sin(Pigreco*(yi-.5))-sin(Pigreco*(zi-.5)));
//         u_value[1] =10*sin(Pigreco*yi)*sin(Pigreco*yi)*(sin(Pigreco*(zi-.5))-sin(Pigreco*(xi-.5)));
//         u_value[2] =10*sin(Pigreco*zi)*sin(Pigreco*zi)*(sin(Pigreco*(xi-.5))-sin(Pigreco*(yi-.5)));
// // -------------------------------------------------
// // translation -------------------------
// //     u_value[0] =10.;u_value[1]=0.;u_value[2]=0.;
// 	
//         //       // set velocity field
//        for(int jdim=0; jdim<DIMENSION; jdim++) {
// // //         dof_u[0]=i+j*nxyz[0]+1+idim*offset;  //dof_u[1]=dof_u[0]+offset;
//          dof_u[jdim]=ix+iy*nxyz[0]+iz*nxyz[0]*nxyz[1]+1+jdim*offset;
// 	    V_SetCmp(&_uvw,dof_u[jdim],u_value[jdim]);
//         }
// //         V_SetCmp(&_uvw,dof_u,u_value);
// //         V_SetCmp(&_uvw,dof_v,v_value);
// //         V_SetCmp(&_uvw,dof_w,w_value);
//       }
//     }
//   }
// 
// 
// 
// 
//   // loop element --------------------------------------
// //   for (unsigned int iel=0 ; iel <n_elem ; ++iel) {
// //     unsigned int elem_gidx=iel*n_u_dofs;
// //     // the local nodes
// //     for (unsigned int i=0; i<n_u_dofs; i++) {
// //       // coordinates
// //       int k=map_nodes[elem_gidx+i];
// //       const Real xi = xyz_glob[k];
// //       const Real yj = xyz_glob[k+offset];
// //       const Real zk = xyz_glob[k+offset2];
// //
// //
// //       // dofs
// //       unsigned int dof_u=k+1; unsigned int dof_v=k+offset+1;
// //       unsigned int dof_w=k+offset2+1;
// //       // Set the initial velocity (note:: nondimensional field)
// //
// // // translation -------------------------
// //     const Real u_value =10.;const Real v_value=-0.;const Real w_value=0.;
// // // ------------------
// // //rotation --------------------------
// //     //const Real u_value=3.14159265358979*(yj-0.5);
// //     //const Real v_value=-3.14159265358979*(xi-0.5);
// //     //const Real w_value=0.;
// // // --------------------------------------------
// // // vortex 2D ---------------------------
// // //       const Real v_value =  ((double) nx)*(PSI(xi+hx,yj)-PSI(xi,yj));
// // //       const Real u_value = - ((double) nx)*(PSI(xi,yj+hx)-PSI(xi,yj));
// // //       const Real w_value = 0;
// // // ---------------------------------
// // // vortex 3D ---------------------------
// // //         const double Pigreco=3.14159265358979;
// // //         const Real u_value =10* sin(Pigreco*xi)*sin(Pigreco*xi)*(sin(Pigreco*(yj-.5))-sin(Pigreco*(zk-.5)));
// // //         const Real v_value =10*sin(Pigreco*yj)*sin(Pigreco*yj)*(sin(Pigreco*(zk-.5))-sin(Pigreco*(xi-.5)));
// // //         const Real w_value =10*sin(Pigreco*zk)*sin(Pigreco*zk)*(sin(Pigreco*(xi-.5))-sin(Pigreco*(yj-.5)));
// // // -------------------------------------------------
// //       // set velocity field
// //       V_SetCmp(&sol,dof_u,u_value);
// //       V_SetCmp(&sol,dof_v,v_value);
// //       V_SetCmp(&sol,dof_w,w_value);
// //
// //     } // --------------------------------------------------------
// //   } // end of element loop
//   return;
// }
// // -------------------------------------------------------





// // ======================================================================================
// //       fine cc from  cc1[Level]
// // ======================================================================================
//
//
//
// /// This function prints the time xml file (to run the single time step file ccf.#.xmf)
// void  MGSolCC::print_time_fine_xmf(
//   const int t_init,       ///<  intial time
//   const int n_time_step,  ///<  number of time steps
//   const int print_step    ///<  print every print_step
// ) {// ===================================================================================
//
//   // file name (time_fine.xmf) ----------------------------------------------------------
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   sprintf(outf2d,"RESU/time_fine.xmf");  std::ofstream out(outf2d);
//
//   //  file xmf text  template -----------------------------------------------------------
//   out << "<?xml version=\"1.0\" ?> \n";
//   out << "<!DOCTYPE Xdmf SYSTEM "  <<  "\"" <<  "../DATA" << "\"" << "[]>\n";
//   out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
//   out << "<Domain> \n";
//   out << "<Grid Name=\" ccf \"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
//   // time loop for grid sequence +++++++++++++++++++++++++++++++++++++++++++++++
//   for(int it=t_init; it<=t_init+n_time_step; it++)
//     if(it%print_step ==0)   {
//       out << "<xi:include href=\""<< "ccf."<< it <<  ".xmf" << "\""
//           << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\" >\n";
//       out << "<xi:fallback />\n";
//       out << " </xi:include>\n";
//     } // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   out << "</Grid> \n"; // Grid Collection end
//   out << "</Domain> \n";
//   out << "</Xdmf> \n";
//   // close and clear --------------------------------------------------------------------
//   out.close();
//   delete [] Named; delete []outf2d;
//   return;
// }
//
// // ======================================================================================
// /// This function prints the fine color file xmf for fine cc (c1[Level])
// /// only for cells 0 < C < 1.
// /// It is called by the print fine function (print_fine_hdf5())
// void MGSolCC::print_fine_xmf(
//   int flag,                    ///<    print flag (time)
//   int n_cell_fine,             ///< number of fine cells with 0< C<1
//   const unsigned int Level     ///< fine Level
// ) {// ===================================================================================
//
//   // file name (ccf(flag).xmf) ----------------------------------------------------------
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   sprintf(Named,"%d",flag); sprintf(outf2d,"./RESU/ccf.%s.xmf",Named);
//   std::ofstream out(outf2d);
//   char *Name_type; Name_type=new char[30];sprintf(Name_type,"%s","Quadrilateral");
//   if(_dim==3)sprintf(Name_type,"%s","Hexahedron");
//   // cell topology (to define cell coordinates) -----------------------------------------
//   int npt_elem=4;  if(_dim==3) npt_elem=8;
//
//   // xmf file template ------------------------------------------------------------------
//   out << "<?xml version=\"1.0\" ?> \n";
//   out << "<Xdmf>  \n";
//   out << "<Domain> \n";
//   out << "<Grid Name=\"Mesh\"> \n";
//   out << "<Time Value =\""<< flag <<"\" />  \n";
//   out << "<Topology Type=\""<<Name_type<<"\"  Dimensions=\""<<n_cell_fine<<"\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""
//       << n_cell_fine<< " "<< npt_elem<<"\" Format=\"HDF\"> \n";
//   out << "ccf."<< flag <<".h5:conn \n";
//   out << "</DataStructure>  \n";
//   out << "</Topology> \n";
//   out << "<Geometry Type=\"X_Y_Z\">  \n";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_cell_fine*npt_elem<< "  1\" Format=\"HDF\">  \n";
//   out << "ccf."<< flag <<".h5:X1 \n";
//   out << "</DataStructure> ";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_cell_fine*npt_elem<< " 1\" Format=\"HDF\"> \n" ;
//   out << "ccf."<< flag <<".h5:X2 \n";
//   out << "</DataStructure>  \n";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_cell_fine*npt_elem<< " 1\" Format=\"HDF\">  \n";
//   out << "ccf."<< flag <<".h5:X3 \n";
//   out << "</DataStructure>  \n";
//   out << " </Geometry> \n";
//   out << " <Attribute Name=\"C\" AttributeType=\"Scalar\" Center=\"Cell\" \n>";
//   out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_cell_fine<< "  1\" Format=\"HDF\"> \n ";
//   out << "ccf."<< flag <<".h5:CCf \n";
//   out << "</DataItem> \n";
//   out << "</Attribute> \n";
//   out << "</Grid> \n";
//   out << "</Domain> \n";
//   out << "</Xdmf> \n";
//   out.close();
//   return;
// }
//
//
// // ======================================================================================
// /// This function prints the fine color function and the file xmf. Only for cells 0< C<1.
// void MGSolCC::print_fine_hdf5(
//   int flag,                   ///<    print flag (time)
//   const unsigned int Level    ///<    fine Level
// ) { // ==================================================================================
//
//   // storage in hf5 file (Xdmf) ---------------------------------------------------------
//   std::ostringstream namefile;  namefile  << "./RESU/ccf."<< flag <<".h5";
//   std::cout << namefile.str() << std::endl;
//   hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
//
//   // element ----------------------------------------------------------------------------
//   int npt_element=4;  int ne_z =1;                            // points for element
//   const int ne_x =_nxyz[0][Level]-1; const int ne_y =_nxyz[1][Level]-1; // number of elements
//   if(_dim==3) {
//     npt_element=8;
//     #if DIMENSION==3
//     ne_z =_nxyz[2][Level]-1;
// #endif
//   }          // 3D
//   double hx=1./ne_x; double hy=1./ne_y; double hz=1./ne_z;    // cell dim (hx,hy)
//   // number of fine cell with 0< C <1
//   int n_cell=0; for (int ik=1;ik<=ne_y*ne_z;ik++) n_cell += M__GetLen(&c1[Level],ik);
// //    for(int ik=1; ik<=ne_y; ik++) n_cell += M__GetLen(&c1[Level],ik);
//
//   print_fine_xmf(flag,n_cell,Level); // print the corresponding xmf file
//
//   // coordinates ------------------------------------------------------------------------
//   double * xxx=new double[n_cell*npt_element];
//   double * yyy=new double[n_cell*npt_element];
//   double * zzz=new double[n_cell*npt_element];
//
//   if(_dim==2) { // coord in 2D ----------------------------------------------------------
//     int icount=0; // getting the coordinate
//     for(unsigned int iy=0; iy<ne_y; iy++) {
//       int len=M_GetLen(&c1[Level],iy+1);
//       for(unsigned int kx=0; kx<len; kx++) {
//         int ix=M__GetPos(&c1[Level],iy+1,kx)-1;
//         double  ix0= hx*ix; double  iy0= hy*iy;
//         xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
//         zzz[icount*npt_element+0]= 0.;
//         xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
//         zzz[icount*npt_element+1]= 0.;
//         xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
//         zzz[icount*npt_element+2]= 0.;
//         xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
//         zzz[icount*npt_element+3]= 0.;
//         icount++;
//       }
//     }
//   }
//
//   //3d
//     if(_dim==3) { // coord in 2D ----------------------------------------------------------
//     int icount=0; // getting the coordinate
//     for (unsigned int kz=0; kz<ne_z; kz++)
//     for (unsigned int iy=0; iy<ne_y; iy++) {
//       int len=M_GetLen(&c1[Level],iy+kz*ne_y+1);
//       for (unsigned int kx=0; kx<len; kx++) {
//         int ix=M__GetPos(&c1[Level],iy+kz*ne_y+1,kx)-1;
//         double  ix0= hx*ix;double  iy0= hy*iy;
//         double  iz0= hz*kz;
//
//         xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
//         zzz[icount*npt_element+0]= iz0;
//         xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
//         zzz[icount*npt_element+1]= iz0;
//         xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
//         zzz[icount*npt_element+2]= iz0;
//         xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
//         zzz[icount*npt_element+3]= iz0;
//
//         xxx[icount*npt_element+4]= ix0;    yyy[icount*npt_element+4]= iy0;
//         zzz[icount*npt_element+4]= iz0+hz;
//         xxx[icount*npt_element+5]= ix0+hx; yyy[icount*npt_element+5]= iy0;
//         zzz[icount*npt_element+5]= iz0+hz;
//         xxx[icount*npt_element+6]= ix0+hx; yyy[icount*npt_element+6]= iy0+hy;
//         zzz[icount*npt_element+6]= iz0+hz;
//         xxx[icount*npt_element+7]= ix0;    yyy[icount*npt_element+7]= iy0+hy;
//         zzz[icount*npt_element+7]= iz0+hz;
//
//         icount++;
//
// //         out << ix0 << " " << iy0 << " " << iz0 << " ";
// //         out << ix0+hx << " " << iy0 << " " << iz0 << " ";
// //         out << ix0+hx << " " << iy0+hy << " " << iz0 << " ";
// //         out << ix0 << " " << iy0+hy << " " << iz0 << " ";
// //         out << ix0 << " " << iy0 << " " << iz0+hz << " ";
// //         out << ix0+hx << " " << iy0 << " " << iz0+hz << " ";
// //         out << ix0+hx << " " << iy0+hy << " " << iz0+hz << " ";
// //         out << ix0 << " " << iy0+hy << " " << iz0+hz << std::endl;
//       }
//     }
//     }
//
//
//
//
//
// //   for(unsigned int iy=0; iy<ne_y; iy++) {
// //       int len=M_GetLen(&c1[Level],iy+1);
// //       for(unsigned int kx=0; kx<len; kx++) {
// //         int ix=M__GetPos(&c1[Level],iy+1,kx)-1;
// //         double  ix0= hx*ix; double  iy0= hy*iy;
// //         xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
// //         zzz[icount*npt_element+0]= 0.;
// //         xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
// //         zzz[icount*npt_element+1]= 0.;
// //         xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
// //         zzz[icount*npt_element+2]= 0.;
// //         xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
// //         zzz[icount*npt_element+3]= 0.;
// //         icount++;
// //       }
// //     }
// //   }
//   // writing _NoNodes -------------------------------------------------------------------
//   hsize_t dimsf[2]; dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;
//   // coord x -> x1 --------------------
//   std::ostringstream Name; Name << "X1";          // name hdf5 dir
//   hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace
//   hid_t dtset=H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                         ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                        );                          // dataset -> xxx
//   H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,xxx);//write
//   H5Sclose(dtsp); H5Dclose(dtset);        // close dataspace dataset
//   delete[]xxx;                            // delete coord  hdf5 print vector
//   // coord x -> x2 --------------------
//   std::ostringstream Name2; Name2 << "X2";         // name hdf5 dir
//   hid_t dtsp2 = H5Screate_simple(2, dimsf, NULL);  // dataspace
//   hid_t  dtset2=H5Dcreate(file,Name2.str().c_str(),H5T_NATIVE_DOUBLE,dtsp2,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                           ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                          );                       // dataset -> yyy
//   H5Dwrite(dtset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,yyy);//write
//   H5Sclose(dtsp2); H5Dclose(dtset2);           // close dataspace dataset
//   delete[]yyy;                                 // delete coord yyy hdf5 print vector
//   // coord x -> x3 -------------------
//   std::ostringstream Name3; Name3 << "X3";                 // name hdf5 dir
//   hid_t dtsp3 = H5Screate_simple(2, dimsf, NULL);          // dataspace
//   hid_t  dtset3=H5Dcreate(file,Name3.str().c_str(),H5T_NATIVE_DOUBLE,dtsp3,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                           ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                          );                                // dataset -> zzz
//   H5Dwrite(dtset3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,zzz);//write
//   H5Sclose(dtsp3); H5Dclose(dtset3);                       // close dataspace dataset
//   delete[]zzz;                                   // delete coord zzz hdf5 print vector
//
//   // connectivity -----------------------------------------------------------------------
//   int *cc_conn=new int[n_cell*npt_element];        // connectivity hdf5 print vector
//   for(int iy=0; iy<n_cell*npt_element; iy++) cc_conn[iy]=iy;
//   dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;     // dimension hdf5 print vector
//   std::ostringstream Nameconn; Nameconn << "conn"; // name hdf5 dir
//   hid_t dtspconn = H5Screate_simple(2, dimsf, NULL);// dataspace
//   hid_t  dtsetconn=H5Dcreate(file,Nameconn.str().c_str(),H5T_NATIVE_INT,dtspconn,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                              ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                             );                 // dataset -> conn
//   H5Dwrite(dtsetconn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_conn);// write
//   H5Sclose(dtspconn); H5Dclose(dtsetconn);    // close dataspace dataset
//   delete[]cc_conn;                            // delete connectivity hdf5 print vector
//
//   // color function --------------------------------------------------------------------
//   double *cc_tmp=new double[n_cell];
//   int icount=0; // getting the color function from c1[Level] ---------
//
//     // color function
//   for (unsigned int kz=0; kz<ne_z; kz++)
//     for (unsigned int iy=0; iy<ne_y; iy++) {
//       int len=M_GetLen(&c1[Level],iy+kz*ne_y+1);
//       for (unsigned int kx=0; kx<len; kx++) {
//         double val=M__GetVal(&c1[Level],iy+kz*ne_y+1,kx);// if (val>1.) val=1.;
//         cc_tmp[icount]= val;icount++;
//       }
//     }
//
//
//
// //   for(unsigned int iy=0; iy<ne_y; iy++) {
// //     int len=M_GetLen(&c1[Level],iy+1);
// //     for(unsigned int kx=0; kx<len; kx++) {
// //       double val=M__GetVal(&c1[Level],iy+1,kx); if(val>1.) val=1.;
// //       cc_tmp[icount]= val; icount++;
// //     }
// //   }
//   // hdf5 storage ----------------------------------------------------
//   dimsf[0] =n_cell;  dimsf[1] = 1;                 // dimension hdf5 print vector
//   std::ostringstream Nameccf; Nameccf << "CCf";     // name hdf5 dir
//   hid_t dtspccf = H5Screate_simple(2, dimsf, NULL); // dataspace
//   hid_t  dtsetccf=H5Dcreate(file,Nameccf.str().c_str(),H5T_NATIVE_DOUBLE,dtspccf,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                             ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                            );                       // dataset -> CCf
//   H5Dwrite(dtsetccf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_tmp);// write
//   H5Sclose(dtspccf); H5Dclose(dtsetccf);            // close dataspace dataset
//   delete[]cc_tmp;                                   // delete connectivity hdf5 print vector
//
//   // close file -------------------------------------------------------------------------
//   H5Fclose(file);
//
//
//   return;
// }
//
//
//
//
// // ======================================================================================
// //       Point color cc from  cc (vector nx_pt*ny_pt)
// // ======================================================================================
//
// // ======================================================================================
// /// This function prints the time xml file (to run the single time step file ccf.#.xmf)
// void  MGSolCC::print_time_cc_xmf(
//   const int t_init,       ///<  intial time
//   const int n_time_step,  ///<  number of time steps
//   const int print_step    ///<  print every print_step
// ) {// ===================================================================================
//
//   // file name (time_fine.xmf) ----------------------------------------------------------
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   sprintf(outf2d,"RESU/time_cc.xmf");  std::ofstream out(outf2d);
//
//   //  file xmf text  template -----------------------------------------------------------
//   out << "<?xml version=\"1.0\" ?> \n";
//   out << "<!DOCTYPE Xdmf SYSTEM "  <<  "\"" <<  "../DATA" << "\"" << "[]>\n";
//   out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
//   out << "<Domain> \n";
//   out << "<Grid Name=\" cc \"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
//   // time loop for grid sequence +++++++++++++++++++++++++++++++++++++++++++++++
//   for(int it=t_init; it<=t_init+n_time_step; it++)
//     if(it%print_step ==0)   {
//       out << "<xi:include href=\""<< "cc."<< it <<  ".xmf" << "\""
//           << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\" >\n";
//       out << "<xi:fallback />\n";
//       out << " </xi:include>\n";
//     } // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   out << "</Grid> \n"; // Grid Collection end
//   out << "</Domain> \n";
//   out << "</Xdmf> \n";
//   // close and clear --------------------------------------------------------------------
//   out.close();
//   delete [] Named; delete []outf2d;
//   return;
// }
//
// // ======================================================================================
// /// This function prints the cc color file xmf for  cc vector
// /// It is called by the print fine function (print_fine_hdf5())
// void MGSolCC::print_cc_xmf(
//   int flag,                    ///<    print flag (time)
//   const unsigned int Level     ///< fine Level
// ) {// ===================================================================================
//
//   // file name (ccf(flag).xmf) ----------------------------------------------------------
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   sprintf(Named,"%d",flag); sprintf(outf2d,"./RESU/cc.%s.xmf",Named);
//   std::ofstream out(outf2d);
//
//    // element ----------------------------------------------------------------------------
//   char *Named_top;Named_top=new char[30];sprintf(Named_top,"Quadrilateral"); // name
//   int npt_element=4;              // points for element
//   const int ne_x=_nxyz[0][Level]-1; const int ne_y=_nxyz[1][Level]-1;int ne_z =1;// n elements
//   int n_cell=ne_x*ne_y;  int n_pts=(ne_x+1)*(ne_y+1); // 2D n elements and points
//   double hx=1./ne_x; double hy=1./ne_y; ;  double hz=0.;    // cell dim (hx,hy)
//   if(_dim==3) { // 3D
//     npt_element=8;sprintf(Named_top,"Hexahedron"); // 3D element name
// #ifndef DIMS2
//     ne_z =_nxyz[2][Level]-1;
// #endif
//      n_cell=ne_x*ne_y*ne_z;n_pts=(ne_x+1)*(ne_y+1)*(ne_z+1); // 3D n element
//      hz=1./ne_z;
//   }
//
//   // xmf file template ------------------------------------------------------------------
//   out << "<?xml version=\"1.0\" ?> \n";
//   out << "<Xdmf>  \n";
//   out << "<Domain> \n";
//   out << "<Grid Name=\"Mesh\"> \n";
//   out << "<Time Value =\""<< flag <<"\" />  \n";
//   out << "<Topology Type=\""<<Named_top<<"\"  Dimensions=\""<<n_cell<<"\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""
//       << n_cell<< " "<< npt_element<<"\" Format=\"HDF\"> \n";
//   out << "mesh_0.h5:MSHCONN \n";//<< flag <<".h5:conn \n";
//   out << "</DataStructure>  \n";
//   out << "</Topology> \n";
//   out << "<Geometry Type=\"X_Y_Z\">  \n";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_pts<< "  1\" Format=\"HDF\">  \n";
//   out << "mesh_0.h5:/NODES/COORD/X1 \n";//cc."<< flag <<".h5:X1 \n";
//   out << "</DataStructure> ";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_pts<< " 1\" Format=\"HDF\"> \n" ;
//   out << "mesh_0.h5:/NODES/COORD/X2 \n";//"cc."<< flag <<".h5:X2 \n";
//   out << "</DataStructure>  \n";
//   out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_pts<< " 1\" Format=\"HDF\">  \n";
//   out << "mesh_0.h5:/NODES/COORD/X3 \n";//"cc."<< flag <<".h5:X3 \n";
//   out << "</DataStructure>  \n";
//   out << " </Geometry> \n";
//   out << " <Attribute Name=\"C\" AttributeType=\"Scalar\" Center=\"Node\" \n>";
//   out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_pts<< "  1\" Format=\"HDF\"> \n ";
//   out << "cc."<< flag <<".h5:CC \n";
//   out << "</DataItem> \n";
//   out << "</Attribute> \n";
//   out << "</Grid> \n";
//   out << "</Domain> \n";
//   out << "</Xdmf> \n";
//   // end xmf file template --------------------------------------------------------------
//   out.close();
//   return;
// }
//
//
//
// // ======================================================================================
// /// This function prints the fine color function and the file xmf. Only for cells 0< C<1.
// void MGSolCC::print_cc_hdf5(
//   int flag,                   ///<    print flag (time)
//   const unsigned int Level    ///<    fine Level
// ) { // ==================================================================================
//
//   // storage in hf5 file (Xdmf) ---------------------------------------------------------
//   std::ostringstream namefile;  namefile  << "./RESU/cc."<< flag <<".h5";
//   std::cout << namefile.str() << std::endl;
//   hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
//
//   // element ----------------------------------------------------------------------------
//   int npt_element=4;  int ne_z =1;                            // points for element
//   const int ne_x =_nxyz[0][Level]-1; const int ne_y =_nxyz[1][Level]-1; // number of elements
//   int n_cell=ne_x*ne_y; int n_pts=(ne_x+1)*(ne_y+1);
//   if(_dim==3) { npt_element=8;
// #ifndef DIMS2
//     ne_z =_nxyz[2][Level]-1;
// #endif
//     n_cell=ne_x*ne_y*ne_z;n_pts=(ne_x+1)*(ne_y+1)*(ne_z+1); }          // 3D
//   double hx=1./ne_x; double hy=1./ne_y; double hz=1./ne_z;    // cell dim (hx,hy)
//   // number of fine cell with 0< C <1
//
//
//   print_cc_xmf(flag,Level); // print the corresponding xmf file
//
// //   // coordinates ------------------------------------------------------------------------
// //   double * xxx=new double[n_cell*npt_element];
// //   double * yyy=new double[n_cell*npt_element];
// //   double * zzz=new double[n_cell*npt_element];
// //
// //   if(_dim==2) { // coord in 2D ----------------------------------------------------------
// //     int icount=0; // getting the coordinate
// //     for(unsigned int iy=0; iy<ne_y; iy++) {
// //       int len=M_GetLen(&c1[Level],iy+1);
// //       for(unsigned int kx=0; kx<len; kx++) {
// //         int ix=M__GetPos(&c1[Level],iy+1,kx)-1;
// //         double  ix0= hx*ix; double  iy0= hy*iy;
// //         xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
// //         zzz[icount*npt_element+0]= 0.;
// //         xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
// //         zzz[icount*npt_element+1]= 0.;
// //         xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
// //         zzz[icount*npt_element+2]= 0.;
// //         xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
// //         zzz[icount*npt_element+3]= 0.;
// //         icount++;
// //       }
// //     }
// //   }
// //   // writing _NoNodes -------------------------------------------------------------------
// //   hsize_t dimsf[2]; dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;
// //   // coord x -> x1 --------------------
// //   std::ostringstream Name; Name << "X1";          // name hdf5 dir
// //   hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace
// //   hid_t dtset=H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
// // #if HDF5_VERSIONM == 1812
// //                         ,H5P_DEFAULT,H5P_DEFAULT
// // #endif
// //                        );                          // dataset -> xxx
// //   H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,xxx);//write
// //   H5Sclose(dtsp); H5Dclose(dtset);        // close dataspace dataset
// //   delete[]xxx;                            // delete coord  hdf5 print vector
// //   // coord x -> x2 --------------------
// //   std::ostringstream Name2; Name2 << "X2";         // name hdf5 dir
// //   hid_t dtsp2 = H5Screate_simple(2, dimsf, NULL);  // dataspace
// //   hid_t  dtset2=H5Dcreate(file,Name2.str().c_str(),H5T_NATIVE_DOUBLE,dtsp2,H5P_DEFAULT
// // #if HDF5_VERSIONM == 1812
// //                           ,H5P_DEFAULT,H5P_DEFAULT
// // #endif
// //                          );                       // dataset -> yyy
// //   H5Dwrite(dtset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,yyy);//write
// //   H5Sclose(dtsp2); H5Dclose(dtset2);           // close dataspace dataset
// //   delete[]yyy;                                 // delete coord yyy hdf5 print vector
// //   // coord x -> x3 -------------------
// //   std::ostringstream Name3; Name3 << "X3";                 // name hdf5 dir
// //   hid_t dtsp3 = H5Screate_simple(2, dimsf, NULL);          // dataspace
// //   hid_t  dtset3=H5Dcreate(file,Name3.str().c_str(),H5T_NATIVE_DOUBLE,dtsp3,H5P_DEFAULT
// // #if HDF5_VERSIONM == 1812
// //                           ,H5P_DEFAULT,H5P_DEFAULT
// // #endif
// //                          );                                // dataset -> zzz
// //   H5Dwrite(dtset3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,zzz);//write
// //   H5Sclose(dtsp3); H5Dclose(dtset3);                       // close dataspace dataset
// //   delete[]zzz;                                   // delete coord zzz hdf5 print vector
// //
// //   // connectivity -----------------------------------------------------------------------
// //   int *cc_conn=new int[n_cell*npt_element];        // connectivity hdf5 print vector
// //   for(int iy=0; iy<n_cell*npt_element; iy++) cc_conn[iy]=iy;
// //   dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;     // dimension hdf5 print vector
// //   std::ostringstream Nameconn; Nameconn << "conn"; // name hdf5 dir
// //   hid_t dtspconn = H5Screate_simple(2, dimsf, NULL);// dataspace
// //   hid_t  dtsetconn=H5Dcreate(file,Nameconn.str().c_str(),H5T_NATIVE_INT,dtspconn,H5P_DEFAULT
// // #if HDF5_VERSIONM == 1812
// //                              ,H5P_DEFAULT,H5P_DEFAULT
// // #endif
// //                             );                 // dataset -> conn
// //   H5Dwrite(dtsetconn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_conn);// write
// //   H5Sclose(dtspconn); H5Dclose(dtsetconn);    // close dataspace dataset
// //   delete[]cc_conn;                            // delete connectivity hdf5 print vector
//
//   // color function --------------------------------------------------------------------
//
// //     unsigned  int nx=_nxyz[0][Level_vof];
// //   const unsigned int level=mgmesh._NoLevels-1-Level_vof;
// //   const unsigned int goffset=mgmesh._NoNodes[mgmesh._NoLevels-1];
// //   double hx=1./ (nx-1);
// //   // std::cout << nx << " ";
// //   if(mode == 1) {  // for quad 9 --------------------------------------------------------
// //     int conn[4][4];
// //     conn[0][0]=0; conn[0][1]=4; conn[0][2]=8; conn[0][3]=7;// 1-subdivision quad9
// //     conn[1][0]=4; conn[1][1]=1; conn[1][2]=5; conn[1][3]=8;// 2-subdivision quad9
// //     conn[2][0]=8; conn[2][1]=5; conn[2][2]=2; conn[2][3]=6;// 3-subdivision quad9
// //     conn[3][0]=7; conn[3][1]=8; conn[3][2]=6; conn[3][3]=3;// 4-subdivision quad9
// //     for(unsigned int el=0; el<offset; el++)
// //       for(unsigned int se=0; se<4; se++)      {
// //         int kn=mgmesh._elem_map[level][conn[se][0]+el*9];// 0-node  in subdivision
// //         double xp=mgmesh._xyz[kn]; double yp=mgmesh._xyz[kn+goffset]; // 0-node coord
// //         int  ix= (int)(xp/hx+0.5); int  iy= (int)(yp/hx+0.5); // central coordinates
// //         int ind=ix+iy*nx;// cartesian index
// //         out <<   cc.Cmp[ind+1] << " ";
// //       }
// //   } else { // for quad 4 ------------------------------------------------------------------
// // //     int nx=_nxyz[0][Level]; int ny=_nxyz[1][Level];
// // //     double hx=
// //     for(unsigned int el=0; el<offset; el++)  {
// // //       int kn=mgmesh._elem_map[level][0+el*4]; // 0-node ->(-1,-1) in quad4
// // //       double xp=mgmesh._xyz[kn]; double yp=mgmesh._xyz[kn+goffset]; // 0-node coord
// // //       int  ix= (int)(xp/hx+0.5); int  iy= (int)(yp/hx+0.5); // central coordinates
// // //       int ind=ix+iy*nx; // cartesian index
// // //       out <<   cc.Cmp[ind+1] << " ";
// //       out <<   cc.Cmp[el+1] << " ";
// //     }
// //
// //   }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//   double *cc_tmp=new double[n_pts];
//   int icount=0; // getting the color function from c1[Level] ---------
//       for(unsigned int ipt=0; ipt<n_pts; ipt++)  {
// //       int kn=mgmesh._elem_map[level][0+el*4]; // 0-node ->(-1,-1) in quad4
// //       double xp=mgmesh._xyz[kn]; double yp=mgmesh._xyz[kn+goffset]; // 0-node coord
// //       int  ix= (int)(xp/hx+0.5); int  iy= (int)(yp/hx+0.5); // central coordinates
// //       int ind=ix+iy*nx; // cartesian index
// //       out <<   cc.Cmp[ind+1] << " ";
//       cc_tmp[ipt]=  cc.Cmp[ipt+1];
//     }
//
// /*
//   for(unsigned int iy=0; iy<ne_y; iy++) {
//     int len=M_GetLen(&c1[Level],iy+1);
//     for(unsigned int kx=0; kx<len; kx++) {
//       double val=M__GetVal(&c1[Level],iy+1,kx); if(val>1.) val=1.;
//       cc_tmp[icount]= val; icount++;
//     }
//   }*/
//   // hdf5 storage ----------------------------------------------------
//   hsize_t dimsf[2]; dimsf[0] =n_pts;  dimsf[1] = 1;                 // dimension hdf5 print vector
//   std::ostringstream Nameccf; Nameccf << "CC";     // name hdf5 dir
//   hid_t dtspccf = H5Screate_simple(2, dimsf, NULL); // dataspace
//   hid_t  dtsetccf=H5Dcreate(file,Nameccf.str().c_str(),H5T_NATIVE_DOUBLE,dtspccf,H5P_DEFAULT
// #if HDF5_VERSIONM == 1812
//                             ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                            );                       // dataset -> CCf
//   H5Dwrite(dtsetccf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_tmp);// write
//   H5Sclose(dtspccf); H5Dclose(dtsetccf);            // close dataspace dataset
//   delete[]cc_tmp;                                   // delete connectivity hdf5 print vector
//
//   // close file -------------------------------------------------------------------------
//   H5Fclose(file);
//
//
//   return;
// }




#endif
#endif
