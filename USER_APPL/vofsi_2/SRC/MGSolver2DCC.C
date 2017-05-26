#include "vof_config.h"
// C++ include files that we need
#include "Equations_conf.h"
#include <sstream>

#ifdef TWO_PHASE
#if DIMENSION==2
#include "DrawShape2D.C"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <map>
#include "itersolv.h"
//#include "itersolv2.h"
#include "MGSolverCC.h"
//#include "MGMesh.h"

#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#undef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN_VAL 1.e-12
#define PI 3.14159265358979
//#define ALPHA 1
 #define ELVIRA 1
// #define SPLIT  1

// #define YOUNG 1


#define CRLIMIT 20.

// inline double PSI(double xi,double yj)  {
//   // single vortex
//   return (1.*sin(PI*xi) *sin(PI*xi) *sin(PI*yj) *sin(PI*yj) /PI);
//   // 4 vortices
//   // return (.25*sin(4.*pi*(xi+0.5))*cos(4*pi*(yj+0.5))/pi);
// 
// }


// --------------------------------
//  IO functions
// --------------------------------

/// write solution
// void MGSolCC::WriteSol(const int Level,const std::string& name, const int var) {
//   QVector *sol;
//   if(var == 0) sol=&cc;
//   else if(var == 1) sol=&_mx;
//   else if(var == 2) sol=&_my;
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



// ---------------------------------------------------
/// read solution c1[Level]
// void MGSolCC::ReadSol(const int Level,const std::string& name, const int var) {
// 
//   std::ifstream infile(name.c_str());
//   QVector *sol;
//   if(var == 0) sol=&cc; else if(var == 1) sol=&_mx;
//   else if(var == 2) sol=&_my;
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
// ------------------------------------
// Multigrid  Operators
// ------------------------------------
// -------------------------------------------------------
/// solution of discret problem with multigrid solver
void MGSolCC::MGSolve(/*MGMeshC &mgmesh,*/
                      unsigned int nc_step,unsigned int itime,double dt) {


  for(unsigned int is=0; is<nc_step; is++) {
    // advection
    Adv(_Nlev_cc-1,dt/ (nc_step));
    OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
    // normal
#ifdef ALPHA
    GenNorm(_Nlev_cc-1,&MGSolCC::rec_melv1);
#else
    GenNorm(_Nlev_cc-1,&MGSolCC::rec_elv1);
#endif
    OldSol_update(_Nlev_cc-1,_mx1[_Nlev_cc-1],_mx1_old[_Nlev_cc-1]);
    OldSol_update(_Nlev_cc-1,_my1[_Nlev_cc-1],_my1_old[_Nlev_cc-1]);
    // Restriction
    for(int level=_Nlev_cc-2; level>=0; level--) {
      RestSol(level); OldSol_update(level,c1[level],c1_old[level]);
      RestNormal(level);
    }
    CCSol_update();
    Normal_update();
  }
  return;
}
// --------------------------------------------------------
/// Projection operator for color function c1[Level c]->c1[Level f]
void MGSolCC::ProjSol(const int Levf) {
  int Levc=Levf-1;
  const unsigned int nxc=_nxyz[0][Levc]-1;  const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxf=_nxyz[0][Levf]-1;  const unsigned int nyf=_nxyz[1][Levf]-1;

  // ghost cells (bandh)
  int memh=3/2+1; int ixy[2]; double mxy[3];
  for(unsigned int ixc=0; ixc<= (nxc+1) * (memh-1); ixc++) _tempc1[ixc]=0.;
  // first lines
  for(unsigned int kc=0; kc<3-memh+1; kc++) ExpRow(&c1_old[Levf],_tempc1,Levc,kc,kc*(nxc+1));

  // domain cells
  for(int jyc=0; jyc<nyc; jyc++) {  // line jy
    for(int ixc=0; ixc<nxc; ixc++) {    // column ix
      ixy[0]=ixc; ixy[1]=jyc;
      int indc=(jyc%3)*(nxc+1)+ixc;
      double valc=_tempc1[indc];
      // add refined values
      if(valc <1. && valc >0.) {  // Area is split into 4 pieces
        // normal
        rec_elv1(Levc,ixy,mxy);
        double mmx=mxy[0]; double mmy=mxy[1]; double alp=mxy[2];
        // Area is split into 4 pieces
        for(int jj1=0; jj1<2; jj1++) for(int ii1=0; ii1<2; ii1++)
            _tempc2[2*ixc+ii1+jj1*(nxf+1)] =4.*get_vol(mmx,mmy,0.,alp,ii1*0.5,jj1*0.5,.5,.5);
      } else { // Area is 1
        for(int jj1=0; jj1<2; jj1++) for(int ii1=0; ii1<2; ii1++)
            _tempc2[2*ixc+ii1+jj1*(nxf+1)]=valc;
      }
    } // global grid ix
    // line _tempc2  contracted into row matrix c1
    CtrRow(_tempc2,&c1[Levf],Levf,2*jyc,0);  //line 1
    CtrRow(_tempc2,&c1[Levf],Levf,2*jyc+1,nxf+1);  //line 2
    // line expand row matrix c1_old -> tempc1
    ExpRow(&c1_old[Levc],_tempc1,Levc,jyc+memh, ((jyc+memh) %3) * (nxc+1));
  }// global grid jy
  return;
}
// -----------------------------------
/// restriction operator for color function c1[Level fine]-> c1[Level coarse]
void MGSolCC::RestSol(const int Levc) {
  int Levf=Levc+1;
  const unsigned int nxlc=_nxyz[0][Levc]-1; const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxfl=_nxyz[0][Levf]+1; const unsigned int nyf=_nxyz[1][Levf]-1;
  //OldSol_update(Levf);
  for(unsigned int jyc=0; jyc<nyc; jyc++)  {
    int jyf=jyc*2;
    //for (unsigned int ixf=0; ixf<=(nxf+1)*2; ixf++){_tempc1[ixf]=0.;_tempc2[ixf]=0.;}

    // Expanding matrix c1_old -> new line _tempc1
    ExpRow(&c1_old[Levf],_tempc1,Levf,jyf,0); // line 1
    ExpRow(&c1_old[Levf],_tempc1,Levf,jyf+1,nxfl); // line 2
    // Area restriction
    for(unsigned int ixc=1; ixc<nxlc; ixc++) {
      int ixf=ixc*2-1;
      _tempc2[ixc]=0.25* (
                     _tempc1[ixf]+_tempc1[ixf+1]+_tempc1[ixf+nxfl]+_tempc1[ixf+nxfl+1]);
    }
    // Contracting matrix tempc2 -> c1
    CtrRow(_tempc2,&c1[Levc],Levc,jyc,0);

  }// global grid jy

  return;
}
/// Projection operator for normal _mxy1[Level c]->_mxy1[Level f]
void MGSolCC::ProjNormal(const int Levf) {
  int Levc=Levf-1;
  const unsigned int nxc=_nxyz[0][Levc]-1;  const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxf=_nxyz[0][Levf]-1;  const unsigned int nyf=_nxyz[1][Levf]-1;

// domain cells
  for(int jyf=0; jyf<nyf; jyf++) {    // line jy
    int jyc=jyf/2;
    int len=M__GetLen(&c1[Levf],jyf+1);
    M_SetLen(&_mx1[Levf],jyf+1,len); M_SetLen(&_my1[Levf],jyf+1,len);
    for(int kxf=0; kxf<len; kxf++) {    // column ix
      double valc=M__GetVal(&c1[Levf],jyf+1,kxf);
      // add refined values
      if(valc <1. && valc >0.) {    // Area is split into 4 pieces


        int ixf= M__GetPos(&c1[Levf],jyf+1,kxf)-1;
        int ixc=ixf/2;
        // unit norm in the coarse level
        double mx=M__GetVal(&_mx1[Levc],jyc+1,ixc);
        double my=M__GetVal(&_my1[Levc],jyc+1,ixc);
        // unit norm in the fine level
        M__SetEntry(&_mx1[Levf],jyf+1,kxf,ixf+1,mx);
        M__SetEntry(&_my1[Levf],jyf+1,kxf,ixf+1,my);

      }
    } // global grid ix
  }// global grid jy
  return;
}
// -----------------------------------
/// restriction operator for color function mxy1Level fine]-> mxy1[Level coarse]
void MGSolCC::RestNormal(const int Levc) {

  int Levf=Levc+1;
  const unsigned int nxc=_nxyz[0][Levc]-1; const unsigned int nyc=_nxyz[1][Levc]-1;
  const unsigned int nxf=_nxyz[0][Levf]-1; const unsigned int nyf=_nxyz[1][Levf]-1;
  //OldSol_update(Levf);
  for(unsigned int jyc=0; jyc<nyc; jyc++)  {
    int jyf=jyc*2;

    // Expanding matrix mx1 -> new line _tempmx1
    ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf,0); // line 1
    ExpRow(&_mx1[Levf],_tempmx1,Levf,jyf+1,nxf+1); // line 2

    // Expanding matrix my1 -> new line _tempmy1
    ExpRow(&_my1[Levf],_tempmy1,Levf,jyf,0); // line 1
    ExpRow(&_my1[Levf],_tempmy1,Levf,jyf+1,nxf+1); // line 2

    // Area restriction
    int lenc=M__GetLen(&c1_old[Levc],jyc+1);
    M_SetLen(&_mx1[Levc],jyc+1,lenc);
    M_SetLen(&_my1[Levc],jyc+1,lenc);
    for(unsigned int kxc=0; kxc<lenc; kxc++) {
      int ixc= M__GetPos(&c1_old[Levc],jyc+1,kxc)-1;

      int ixf=ixc*2;
      double mx=(_tempmx1[ixf]+_tempmx1[ixf+1]+_tempmx1[ixf+nxf+1]+_tempmx1[ixf+nxf+2]);
      double my=(_tempmy1[ixf]+_tempmy1[ixf+1]+_tempmy1[ixf+nxf+1]+_tempmy1[ixf+nxf+2]);
      double mm = fabs(mx)+ fabs(my);
      if(mm > 0) {
        mx = mx/mm; my = my/mm;
      }
      M__SetEntry(&_mx1[Levc],jyc+1,kxc,ixc+1,mx);
      M__SetEntry(&_my1[Levc],jyc+1,kxc,ixc+1,my);
    }

  }// global grid jy
  return;
}

// --------------------------
// Solution
// -----------------------------------

// ======================================================================================
/// This function update the cc color function
void MGSolCC::CCSol_update(
) {// ==================================================================================
  const unsigned int nxl=_nxyz[0][0]; const unsigned int ny=_nxyz[1][0]-1;
  for(unsigned int jy=0; jy<ny; jy++)  {
    // expanding c1_old -> tempc1
    ExpRow(&c1[0],_tempc1,0,jy,0);
    // back to cc
    int indy=jy*nxl;
    for(unsigned int ix=1; ix<nxl; ix++) cc.Cmp[ix+indy]=_tempc1[ix];
  }// global grid jy
  return;
}
// ----------------------------------------
/// Generating solution at all levels < Level
void MGSolCC::GenSol(
  const int Level
) {

  unsigned int nx,ny,nz,dummy;
  std::ifstream infile("../RESU/c0.in");
//   // reading -------------------------------------------
  if(!infile) {
    std::cout <<"Input file "<<infile<<" not opened." << std::endl;
    double xc[2]; // drop center coordinates
    xc[0]=0.25; xc[1]=0.25;
    init_drop(Level,xc,0.15,0.1495);

  } else {
    std::cout << " reading from file ../RESU/c0.in" << std::endl;
    // reading from file "../RESU/c0.in"
    const int  bufLen = 256; char  buf[bufLen+1];
    while(strncmp(buf,"Level",5) != 0)    infile >> buf;
    infile >> dummy;
    if(_Nlev_cc != dummy) std::cout << " _NoLevels is " << _Nlev_cc << std::endl;
    infile >> nx >> ny;
    for(unsigned int i=0; i<nx; i++)  for(unsigned int j=0; j<ny; j++)
        infile >> cc.Cmp[i+j* (nx+1) +1];
  } // ----------------------------------------------
  // Restriction --------------------------------
  OldSol_update(Level,c1[Level],c1_old[Level]);
  GenNorm(Level,&MGSolCC::rec_elv1);
 unsigned int rlen0=M_GetRowDim(&(c1[Level]));
//   unsigned int rlen1=M_GetRowDim(&(c1_old[Level]));
  for(int level=Level-1; level>=0; level--) {
    RestSol(level); OldSol_update(level,c1[level],c1_old[level]);
    RestNormal(level);
//     unsigned int rlen2=M_GetRowDim(&(c1[level]));
//     unsigned int rlen3=M_GetRowDim(&(c1_old[level]));
//     double a=0.;
  }// ----------------------------------------
  // Projection --------------------------
  for(int level=Level+1; level<_Nlev_cc; level++) {
    ProjSol(level); OldSol_update(level,c1[level],c1_old[level]);
  } // -----------------------------------
  // Update ---------------------------
  CCSol_update();   // update the cc color vector (transfer cc vector)
  Normal_update();  // update the normal vector (transfer normal vector _mx _my)
//   unsigned int rlen=M_GetRowDim(&(c1[Level]));
  return;
}

// -------------------------------------------------------

// ======================================================================================
/// This function initialize circle 2D in c1 Matrix at level Level
void MGSolCC::init_drop(
  unsigned int Level,
  double xc[],
  double r,
  double min_y
) {
  // initilize a circle with (x0,y0) = center; r = radius

  const unsigned int nx=_nxyz[0][Level]-1;
  const unsigned int ny=_nxyz[1][Level]-1;
  double hx=1./nx; double hy=1./ny;
  double x0=xc[0]; double y0=xc[1];
  // number of points in cell
  const unsigned int nadd = 20;
  for(unsigned int ix=1; ix<nx+1; ix++)  _tempc2[ix]=0.;

  for(unsigned int j=0; j<ny; j++) {
    for(unsigned int i=0; i<nx; i++) {

      unsigned int ind = i+j* (nx+1);
      double ic = i*hx; double jc = j*hy;
      if(j*hy>min_y)
      ic_read(_tempc2,i,j,0,hx,hy,0,nadd); //this read the inizial conditions for CC

    } // i
    // numbering nonzero elements
    CtrRow(_tempc2,&c1[Level],Level,j,0);

  } // j
//    unsigned int rlen=M_GetRowDim(&(c1[Level]));
//     unsigned int clen=M_GetClmDim(&(c1[Level]));
//   std::cout<< "Initializing c function on a circle of center (" << x0 << ","
//            << y0 << ") and radius " << r << std::endl;
//            " on matrix " << rlen << " x " << clen << std::endl;

  return;
}
// --------------  Old Solution cc ----------------------
/// Operator
/// \f$ Matrix^L(nx\times ny) \rightarrow Matrix^L(nx\times ny) :\f$
/// \f$ \quad c1[L] = c1_{old}[L]  \f$
void MGSolCC::GenOldSol(const int Level) {
  for(unsigned int ind=1; ind<_nxyz[1][Level]; ind++) {
    int len=M_GetLen(&c1[Level],ind);   M_SetLen(&c1_old[Level],ind,len);
    for(unsigned int k=0; k<len; k++)      M_SetEntry(&c1_old[Level],ind,k,k+1,0.);
  }
  OldSol_update(Level,c1[Level],c1_old[Level]);
  // std::cout << "   gen old solution  " << std::endl;
}
// ======================================================================================
/// This function updates the old solution c_old withthe new one c
void MGSolCC::OldSol_update(const int Level_cc,Matrix & c,Matrix & c_old) {
  for(unsigned int iy=1; iy<_nxyz[1][Level_cc]; iy++) {
    int len=M__GetLen(&c,iy); M_SetLen(&c_old,iy,len);
    for(unsigned int kx=0; kx<len; kx++)
      M__SetEntry(&c_old,iy,kx,M__GetPos(&c,iy,kx),M__GetVal(&c,iy,kx));
  }
  std::cout << "  update old sol Level  " << Level_cc << std::endl;
}
// ---------------------------------
// Normal functions
// -----------------------------------
// ------------------------------------
/// Normal Generation from color function
void MGSolCC::GenNorm(const unsigned int Level,
                      void (MGSolCC::*rec1)(const unsigned int,int[], double [])) {

  int nxl=_nxyz[0][Level]+1; int ny=_nxyz[1][Level]-1;
  int ixy[3]; double mxy[3];

  ExpRow(&c1[Level],_tempc1,Level,0,0);// line 0
  ExpRow(&c1[Level],_tempc1,Level,1,nxl);// line 1
  ExpRow(&_mx1_old[Level],_tempmx1,Level,0,0);// line 0
  ExpRow(&_mx1_old[Level],_tempmx1,Level,1,nxl);// line 1
  ExpRow(&_my1_old[Level],_tempmy1,Level,0,0);// line 0
  ExpRow(&_my1_old[Level],_tempmy1,Level,1,nxl);// line 1
  for(unsigned int ix=2*nxl; ix<3*nxl; ix++) {
    _tempc1[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
  }             //line 2


  for(unsigned int jy=0; jy<ny; jy++) {
    int len=M__GetLen(&c1[Level],jy+1);
    M_SetLen(&_mx1[Level],jy+1,len);
    M_SetLen(&_my1[Level],jy+1,len);

    // loop over columns
    for(unsigned int kx=0; kx<len; kx++) {

      int ix=M__GetPos(&c1[Level],jy+1,kx);
      ixy[0]=ix; ixy[1]=jy;
      double val=M__GetVal(&c1[Level],jy+1,kx);
      // computing normal ( 0<c<1)
      if(val <1.)(this->*rec1)(Level,ixy,mxy);
      double mx=mxy[0]; double my=mxy[1]; double alpha=mxy[2];
      M__SetEntry(&_mx1[Level],jy+1,kx,ix,mx);
      M__SetEntry(&_my1[Level],jy+1,kx,ix,my);
    }
    // Row update
    int indl=((jy+2)%3)*nxl;
    ExpRow(&c1[Level],_tempc1,Level,jy+2,indl);
    ExpRow(&_mx1_old[Level],_tempmx1,Level,jy+2,indl);
    ExpRow(&_my1_old[Level],_tempmy1,Level,jy+2,indl);
  }
  return;
}
// ======================================================================================
/// This function updates the normal vector _mx _my (transfer vector)
void MGSolCC::Normal_update(
) {// ===================================================================================
  const unsigned int nxl0=_nxyz[0][0]; const unsigned int ny=_nxyz[1][0]-1;
  for(unsigned int jy=0; jy<ny; jy++)  {
    int indy=jy*nxl0;
    // expanding c1_old -> tempc1
    for(unsigned int ix=1; ix<nxl0; ix++) {
      _tempmx1[ix]=0.; _tempmy1[ix]=0.;
    }
    for(unsigned int lx=0; lx<M_GetLen(&c1[0],jy+1); lx++)  {
      int pos=M__GetPos(&c1[0],jy+1,lx);
      _tempmx1[pos]=M__GetVal(&_mx1[0],jy+1,lx);
      _tempmy1[pos]=M__GetVal(&_my1[0],jy+1,lx);
    }
    // back to cc
    for(unsigned int ix=1; ix<nxl0; ix++) {
      _mx.Cmp[ix+indy]=_tempmx1[ix];
      _my.Cmp[ix+indy]=_tempmy1[ix];
    }
  }// global grid jy
  return;
}

// -----------------------------------------
// RECONSTRUCTION
// ------------------------------------------


// --------------------------------------------------------
/// Metodo Centrale Modificato
void MGSolCC::rec_Cent(const unsigned int Level,int indxy[],
                       double mxy[]) {
  const unsigned int nx=_nxyz[0][Level]-1; const unsigned int ny=_nxyz[0][Level]-1;
  const QVector cc0=cc;
  int ind=indxy[0];
  double hx=1./ (nx); double hy=1./ (ny);
  const int indt=ind+nx+1; const int indb=ind-nx-1;
  double m1,m2;
  // mx
  m2= 0.5* ((cc0.Cmp[indt]+cc0.Cmp[ind] +cc0.Cmp[indb])-
            (cc0.Cmp[indt+2]+cc0.Cmp[ind+2]+cc0.Cmp[indb+2]));
  double signm1 = (m2 >0) ? -1:1;
  // my
  m1= 0.5* ((cc0.Cmp[indb]+cc0.Cmp[indb+1]+cc0.Cmp[indb+2])-
            + (cc0.Cmp[indt]+cc0.Cmp[indt+1]+cc0.Cmp[indt+2]));
  double signm2 = (m1 >0) ? -1:1;

  // minimization
  if(fabs(m1) < fabs(m2)) {
    mxy[1] =m1/ (fabs(m1) +fabs(signm1)); mxy[0] =-signm1/ (fabs(m1) +fabs(signm1));
  } else {
    mxy[1] =-signm2/ (fabs(m2) +fabs(signm2)); mxy[0] =m2/ (fabs(m2) +fabs(signm2));
  }
  return;
}
// ----------------------------------------------------------
/** Youngs local reconstruction */
// ---------------------------------------------------------------
void MGSolCC::rec_Young(const unsigned int Level ,int ixy[],
                        double mxy[]) {
  const unsigned int nx=_nxyz[0][Level]-1; const unsigned int ny=_nxyz[0][Level]-1;
  double hx=1./ (nx); double hy=1./ (ny);
  int ix=ixy[0]; int iy=ixy[1];

  int ind=ix+iy*(nx+1);
  const int indt=ind+nx+1; const int indb=ind-nx-1;
  double mx =- (cc.Cmp[indt+2]+2.* cc.Cmp[ind+2] +cc.Cmp[indb+2] -
                cc.Cmp[indt]-2.* cc.Cmp[ind] -cc.Cmp[indb]) /hx;
  double my =- (cc.Cmp[indt+2]+2.*cc.Cmp[indt+1]+cc.Cmp[indt] -
                cc.Cmp[indb]-2.*cc.Cmp[indb+1]-cc.Cmp[indb+2])  /hy;

  // rec_Cent(Level,ind,mx,my,alpha);
  int invx = 1; int invy = 1;
  if(mx < 0.) {invx = -1; mx = - mx; }
  if(my < 0.) {invy = -1; my = - my;}
  double mm = mx+my;  if(mm != 0.) { mx = mx/mm;  my = my/mm;}
  /* get alpha for the equation of the interface */
  mm = MIN(mx,my);  double alpha = get_alpha(mm,1.-mm,0.,cc.Cmp[ind+1]);
  /* now back to the original line */
  mxy[0] = invx*mx; mxy[1]= invy*my; mxy[1]=alpha +MIN(0.,mx)+MIN(0.,my);
}

// --------------------------------------------------
// 2D Elvira reconstruction
/// elvira reconstruction by using contracing matrices
/// c_old and mx1,my1,mz1
// --------------------------------------------------
void MGSolCC::rec_elv1(unsigned int Level,int ixyz[],double mxy[]) {
  // set up
  int ix=ixyz[0]; int jy=ixyz[1];
  const unsigned int nxl=_nxyz[0][Level]+1;
  double m[2][3],mmx,mmy,alp; double min_err = 16.;

  // index def
  const int ind= ix+ ((jy) %3)*nxl;
  const int indt=ix+((jy+1) %3)*nxl;
  const int indb=ix+((jy+2) %3)*nxl;

  //  horizontal height function	(c_t,c_1,c_b)
  double c_t = _tempc1[indt-1]  + _tempc1[indt] + _tempc1[indt+1];
  double c_1 = _tempc1[ind-1]   + _tempc1[ind]  + _tempc1[ind+1];
  double c_b = _tempc1[indb-1]  + _tempc1[indb] + _tempc1[indb+1];
//  vertical  height function	(c_r,c_2,c_l)
  double c_r = _tempc1[indb+1]+ _tempc1[ind+1]  + _tempc1[indt+1];
  double c_2 = _tempc1[indb]+ _tempc1[ind]  + _tempc1[indt];
  double c_l = _tempc1[indb-1]  + _tempc1[ind-1]    + _tempc1[indt-1];

  double ccc=_tempc1[ind];

  /* angular coefficients from backward centered and forward finite
     differences */
  m[0][0] = c_l - c_2; m[0][1] = 0.5* (c_l-c_r);  m[0][2] = c_2 - c_r;
  m[1][0] = c_b - c_1; m[1][1] = 0.5* (c_b-c_t);  m[1][2] = c_1 - c_t;
  int iyx = 1; if(c_t > c_b)    iyx = -1;
  int ixy = 1; if(c_r > c_l)    ixy = -1;

  /* k=0: y=y(x), k=1: x = x(y) */
  for(int k=0; k<=1; k++) {
    int invy = (1-k) *iyx + k*ixy;
    for(int kk=0; kk<=2; kk++) {
      int invx =1; mmy = 1.; mmx = fabs(m[k][kk]);
      if(m[k][kk] < 0.) invx = -1;
      /* mx and my are now positive, set  mx + my = 1 */
      double mm = mmx+ mmy; double mmx1 = mmx/mm; double mmy1 = mmy/mm;
      /* get alpha for the equation of the interface */
      mm = MIN(mmx1,mmy1); alp = get_alpha(mm,1.-mm,0.,ccc);
      mmx = (k) *invy*mmy1 + (1-k) *invx*mmx1;
      mmy = (1-k) *invy*mmy1 + (k) *invx*mmx1;
      alp += MIN(0.,mmx) + MIN(0.,mmy);

      /* get area diff. with the surrounding cells */
      double sum_area_diff = 0.;
      for(int jj=-1; jj<=1; jj++) {
        for(int ii=-1; ii<=1; ii++) {
          int indl=(ix+ii)+((jy+jj+3)%3)*nxl;
          double area = get_area(mmx,mmy,alp, (double) ii,jj);
          sum_area_diff +=(area-_tempc1[indl])*(area-_tempc1[indl]);
        }
      }
      /* now update the minimum value and final line */
      if(sum_area_diff < min_err) {
        min_err = sum_area_diff; mxy[2] = alp; mxy[0] = mmx; mxy[1] = mmy;
      }
    }
  }
  return;
}
// --------------------------------------------------
// 2D alpha reconstruction
///  alpha reconstruction (elvira modified for low resolutions)
///  by using contracing matrices c_old and mx1,my1,mz1
// --------------------------------------------------
void MGSolCC::rec_melv1(unsigned int Level,int ixyz[],
                        double mxy[]) {

  int ix=ixyz[0]; int jy=ixyz[1];
  double m[2][3],mmx,mmy,alp; double min_err = 16.;
  unsigned int nxl=_nxyz[0][Level]+1; unsigned int ny=_nxyz[1][Level]-1;
  int ind[3][3]; double crw[3][3];  double ccl[3][3];

//   int tab1[9];int invtab1[9]; tab1[0]=0; tab1[1]=1; tab1[2]=2;
//   tab1[3]=7; tab1[4]=8; tab1[5]=3;tab1[6]=6; tab1[7]=5; tab1[8]=4;
//   invtab1[0]=0; invtab1[1]=1; invtab1[2]=2;
//   invtab1[7]=3; invtab1[8]=4; invtab1[3]=5; invtab1[6]=6; invtab1[5]=7; invtab1[4]=8;


  // initial crw and ccl
  for(unsigned int j0=0; j0<3; j0++)
    for(unsigned int i0=0; i0<3; i0++) {
      ind[i0][j0]=ix+i0-1+((jy+j0+2)%3) *nxl;
      crw[i0][j0]= ccl[i0][j0]=_tempc1[ind[i0][j0]];
    }

  double ccc=_tempc1[ind[1][1]];
  double amx=_tempmx1[ind[1][1]];
  double amy=_tempmy1[ind[1][1]];

// correction to color function *********************************
  int kl[3]; kl[1]=1;
  double valfk[3]= {1.,1.,0.};
  if(3.*fabs(amx) < fabs(amy)) {
    kl[0]= (amy >=0) ? 0:2;  kl[2]=2-kl[0];
    // corrected ccl matrix (h column function)
    for(int j1=0; j1<3; j1++) {
      for(int i1=0; i1<3; i1++) {
        double myl=_tempmy1[ind[i1][kl[j1]]]; double mxl=_tempmx1[ind[i1][kl[j1]]];
        if(myl*amy+ mxl*amx<0.) {
          ccl[i1][kl[j1]]=valfk[j1];
          if(j1 == 1 && amy*_tempmy1[ind[i1][kl[0]]]+_tempmx1[ind[i1][kl[0]]]*amx < MIN_VAL)
            ccl[i1][kl[0]]=1.;
        }
      }
    }
  } else if(fabs(amx) >fabs(amy) *3.) {
    kl[0]= (amx >=0) ? 0:2; kl[1]=1; kl[2]=2-kl[0];
    //corrected crw matrix (h row function)
    for(int j1=0; j1<3; j1++) {
      for(int i1=0; i1<3; i1++) {
        double myl=_tempmy1[ind[kl[i1]][j1]];
        double mxl=_tempmx1[ind[kl[i1]][j1]];
        if(mxl*amx+myl*amy < 0.) {
          crw[kl[i1]][j1]=valfk[i1];
          if(i1 == 1 && amy*_tempmy1[ind[kl[0]][j1]]+_tempmx1[ind[kl[0]][j1]]*amx < MIN_VAL)
            crw[kl[0]][j1]=1.;
        }
      }
    }
  }






//   else {
//     int j1= (amx >=0) ? 0:2; int i1= (amy >=0) ? 0:2;
//     int k1=j1+3*i1; int kp=invtab1[(tab1[k1]+1) %8];
//     int  km=invtab1[(tab1[k1]+7) %8];
//
//     // bottom  line
//     int ind0=ind-1+j1+ (i1-1) * (nx+1);
//     double mx0=_mx.Cmp[ind0+1];double my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy < 0.) {
//       crw[i1][j1]=1.; ccl[i1][j1]=1.;
//     }
//
//     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[kp/3][kp%3]=1.; ccl[kp/3][kp%3]=1.;
//     }
//
//     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[km/3][km%3]=1.; ccl[km/3][km%3]=1.;
//     }
//
//     // center line
//     kp=invtab1[(tab1[k1]+2) %8];km=invtab1[(tab1[k1]+6) %8];
//     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[kp/3][kp%3]=1.; ccl[kp/3][kp%3]=1.;
//     }
//
//     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[km/3][km%3]=1.; ccl[km/3][km%3]=1.;
//     }
//
//     // top line
//     int ind2=ind-1+2-j1+ (2-i1-1) * (nx+1);
//     double mx2=_mx.Cmp[ind2+1];double my2=_my.Cmp[ind2+1];
//     if (mx2*amx +my2*amy < 0.) {
//       crw[2-i1][2-j1]=0.;ccl[2-i1][2-j1]=0.;
//     }
//
//     kp=invtab1[(tab1[k1]+5) %8];km=invtab1[(tab1[k1]+3) %8];
//     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[kp/3][kp%3]=0.; ccl[kp/3][kp%3]=0.;
//     }
//
//     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
//     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
//     if (mx0*amx+my0*amy <0.) {
//       crw[km/3][km%3]=0.;ccl[km/3][km%3]=0.;
//     }
//   }
  // *****************************

  // height function
  double c_r = ccl[2][0]+ccl[2][1]+ccl[2][2];
  double c_2 = ccl[1][0]+ccl[1][1]+ccl[1][2];
  double c_l = ccl[0][0]+ccl[0][1]+ccl[0][2];
  double c_t = crw[0][2]+crw[1][2]+crw[2][2];
  double c_1 = crw[0][1]+crw[1][1]+crw[2][1];
  double c_b = crw[0][0]+crw[1][0]+crw[2][0];

  /* angular coefficients from backward centered and forward finite  differences */
  m[0][0] = c_l-c_2; m[0][1] = 0.5* (c_l-c_r); m[0][2] = c_2 - c_r;
  m[1][0] = c_b - c_1; m[1][1] = 0.5* (c_b-c_t);  m[1][2] = c_1 - c_t;
//   int iyx = 1; if (c_t > c_b)    iyx = -1;
//   int ixy = 1; if (c_r > c_l)    ixy = -1;
  int iyx = 1; if(ccl[0][2]+ccl[1][2]+ccl[2][2] > ccl[0][0]+ccl[1][0]+ccl[2][0])    iyx = -1;
  int ixy=  1; if(crw[2][0]+crw[2][1]+crw[2][2] > crw[0][0]+crw[0][1]+crw[0][2])    ixy = -1;

// for (int k=0; k<=1; k++) {
  int k=1; if(fabs(amx) <fabs(amy)) k=0;   {

    // unit normal
    int invy = (1-k) *iyx + k*ixy;
    for(int kk=0; kk<=2; kk++) {
      int invx = 1;   mmy = 1.;     mmx = fabs(m[k][kk]);
      if(m[k][kk] < 0.) invx = -1;
      /* mx and my are now positive, set  mx + my = 1 */
      double mm = mmx+ mmy; double mmx1 = mmx/mm; double mmy1 = mmy/mm;
      /* get alpha for the equation of the interface */
      mm = MIN(mmx1,mmy1); alp = get_alpha(mm,1.-mm,0.,ccc);
      /* now back to the original line */
      mmx = (k) *invy*mmy1 + (1-k) *invx*mmx1; mmy = (1-k) *invy*mmy1 + (k) *invx*mmx1;
      alp += MIN(0.,mmx) + MIN(0.,mmy);

      /* get area diff. with the surrounding cells */
      double sum_area_diff = 0.;     double sum1=0.;
      if(k==0) {
        for(int ii=-1; ii<=1; ii++) for(int jj=-1; jj<=1; jj++) {
            sum1 += ccl[ii+1][jj+1]*ccl[ii+1][jj+1];
            double area = get_area(mmx,mmy,alp, (double) ii, (double) jj);
            sum_area_diff += (area-ccl[ii+1][jj+1]) * (area - ccl[ii+1][jj+1]);
          }
      } else { // k !=0
        for(int ii=-1; ii<=1; ii++) for(int jj=-1; jj<=1; jj++) {
            sum1 +=crw[ii+1][jj+1]*crw[ii+1][jj+1];
            double area = get_area(mmx,mmy,alp, (double) ii, (double) jj);
            sum_area_diff += (area-crw[ii+1][jj+1]) * (area-crw[ii+1][jj+1]);
          }
      }
      /* now update the minimum value and final line */
      sum_area_diff =  sum_area_diff;
      if(sum_area_diff < min_err) {
        min_err = sum_area_diff;
        mxy[2] = alp; mxy[0] = mmx; mxy[1] = mmy;
      }
    }
  }

  return;
}

// --------------------------------
// Cartesian cell evaluation
// ----------------------------------------
/// Return cell velocity field
void MGSolCC::get_vel(int ixy[],int fd,double u[]) {
  // geometry level
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[0][0]-1;
  int n_nodes= (nxc+1)*(nyc+1);
  // velocity field  ul(left),ur(right),vt(top),vb(bot)
  int ixc=ixy[0]/fd; int ixf=ixy[0]%fd;
  int jyc=ixy[1]/fd; int jyf=ixy[1]%fd;
  int indv = ixc+jyc*(nxc+1);
  const unsigned int indy = indv+nxc+1;
  const unsigned int indx = indv+1;
  // u-component
  double u0 =_uvw.Cmp[_invnode_dof[indv]]; // v-right
  double u1 =_uvw.Cmp[_invnode_dof[indx]]; // v-right
  double u2 =_uvw.Cmp[_invnode_dof[indx+(nxc+1)]]; // v-right
  double u3 =_uvw.Cmp[_invnode_dof[indv+(nxc+1)]]; // v-right

  // v-component
  double v0 =_uvw.Cmp[_invnode_dof[indv]+n_nodes]; // v-b left
  double v3 =_uvw.Cmp[_invnode_dof[indy]+n_nodes]; // v- t left
  double v2 =_uvw.Cmp[_invnode_dof[indy+1]+n_nodes]; // v- t right
  double v1 =_uvw.Cmp[_invnode_dof[indv+1]+n_nodes]; // v-b right

  u[0] = (u0*(fd-ixf)*(fd-jyf-0.5)+u1*(ixf)*(fd-jyf-0.5)+
          u2*(ixf)*(jyf+0.5)+u3*(fd-ixf)*(jyf+0.5))/(fd*fd);
  u[1] = (u0*(fd-ixf-1)*(fd-jyf-0.5)+u1*(ixf+1)*(fd-jyf-0.5)+
          u2*(ixf+1)*(jyf+0.5)+u3*(fd-ixf-1)*(jyf+0.5))/(fd*fd);
  u[2] = (v0*(fd-ixf-0.5)*(fd-jyf)+v1*(ixf+0.5)*(fd-jyf)+
          v2*(ixf+0.5)*(jyf)+v3*(fd-ixf-0.5)*(jyf))/(fd*fd);
  u[3] = (v0*(fd-ixf-0.5)*(fd-jyf-1)+v1*(ixf+0.5)*(fd-jyf-1)+
          v2*(ixf+0.5)*(jyf+1)+v3*(fd-ixf-0.5)*(jyf+1))/(fd*fd);
	  
	  
  // Analytic vortex test
  //   double ix3=ix*hx;double iy3=jy*hy;
  //   double ul =-(PSI(ix3,(iy3+hy))-PSI(ix3,iy3)) *dt/(hy*hx);
  //         double ur =-(PSI(ix3+hx,iy3+hy)-PSI(ix3+hx,iy3)) *dt/(hy*hx);
  //         double vb = (PSI(ix3+hx,iy3)-PSI(ix3,iy3)) *dt/(hy*hx);
  //         double vt = (PSI(ix3+hx,iy3+hy)-PSI(ix3,iy3+hy)) *dt/(hy*hx);
	   
// 	  std::cout<<"not Reducing velocity  "<< ixy[1] <<std::endl;
	  
// if (ixy[1]<30){ 
// //   std::cout<<"Reducing velocity  "<< ixy[1] <<std::endl;
//   u[0] = 0; u[1] = 0; u[2] = 0; u[3] = 0;
//   
// }
	  
  return;
}
/// Return cell displacement field
void MGSolCC::get_disp(int ixy[],int fd,double u[]) {
  // geometry level
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[0][0]-1;
  int n_nodes= (nxc+1)*(nyc+1);
  // velocity field  ul(left),ur(right),vt(top),vb(bot)
  int ixc=ixy[0]/fd; int ixf=ixy[0]%fd;
  int jyc=ixy[1]/fd; int jyf=ixy[1]%fd;
  int indv = ixc+jyc*(nxc+1);
  const unsigned int indy = indv+nxc+1;
  const unsigned int indx = indv+1;
  // u-component
  double u0 =_dsvw.Cmp[_invnode_dof[indv]]; // v-right
  double u1 =_dsvw.Cmp[_invnode_dof[indx]]; // v-right
  double u2 =_dsvw.Cmp[_invnode_dof[indx+(nxc+1)]]; // v-right
  double u3 =_dsvw.Cmp[_invnode_dof[indv+(nxc+1)]]; // v-right

  // v-component
  double v0 =_dsvw.Cmp[_invnode_dof[indv]+n_nodes]; // v-b left
  double v3 =_dsvw.Cmp[_invnode_dof[indy]+n_nodes]; // v- t left
  double v2 =_dsvw.Cmp[_invnode_dof[indy+1]+n_nodes]; // v- t right
  double v1 =_dsvw.Cmp[_invnode_dof[indv+1]+n_nodes]; // v-b right

  u[0] = (u0*(fd-ixf)*(fd-jyf-0.5)+u1*(ixf)*(fd-jyf-0.5)+
          u2*(ixf)*(jyf+0.5)+u3*(fd-ixf)*(jyf+0.5))/(fd*fd);
  u[1] = (u0*(fd-ixf-1)*(fd-jyf-0.5)+u1*(ixf+1)*(fd-jyf-0.5)+
          u2*(ixf+1)*(jyf+0.5)+u3*(fd-ixf-1)*(jyf+0.5))/(fd*fd);
  u[2] = (v0*(fd-ixf-0.5)*(fd-jyf)+v1*(ixf+0.5)*(fd-jyf)+
          v2*(ixf+0.5)*(jyf)+v3*(fd-ixf-0.5)*(jyf))/(fd*fd);
  u[3] = (v0*(fd-ixf-0.5)*(fd-jyf-1)+v1*(ixf+0.5)*(fd-jyf-1)+
          v2*(ixf+0.5)*(jyf+1)+v3*(fd-ixf-0.5)*(jyf+1))/(fd*fd);
	  
	  
  // Analytic vortex test
  //   double ix3=ix*hx;double iy3=jy*hy;
  //   double ul =-(PSI(ix3,(iy3+hy))-PSI(ix3,iy3)) *dt/(hy*hx);
  //         double ur =-(PSI(ix3+hx,iy3+hy)-PSI(ix3+hx,iy3)) *dt/(hy*hx);
  //         double vb = (PSI(ix3+hx,iy3)-PSI(ix3,iy3)) *dt/(hy*hx);
  //         double vt = (PSI(ix3+hx,iy3+hy)-PSI(ix3,iy3+hy)) *dt/(hy*hx);
	   
// 	  std::cout<<"not Reducing velocity  "<< ixy[1] <<std::endl;
	  
// if (ixy[1]<30){ 
// //   std::cout<<"Reducing velocity  "<< ixy[1] <<std::endl;
//   u[0] = 0; u[1] = 0; u[2] = 0; u[3] = 0;
//   
// }
	  
  return;
}
// ------------------------------------------------
///  Return line intersect
double MGSolCC::get_alpha(double m1, double m2, double m3, double cc) {
  double alpha;    double V1 = .5*m1/m2;  double V2 = 1. - V1; double mm=m1*m2;
  if(cc < V1)    alpha = sqrt(2.*mm*cc);
  else if(cc <= V2)    alpha = m2*cc + 0.5*m1;
  else    alpha =  1. - sqrt(2.* (1.-cc) *mm);
  return alpha;
}
//     ----------------------------
/// Return line intersect
double MGSolCC::get_alpha2(double mx, double my, double mz, double cc0) {

  /* mx and my are now positive, set  mx + my = 1 */
// double mm = fabs(mx)+ fabs(my); double mmx1 = fabs(mx)/mm; double mmy1 = fabs(my)/mm;
  /* get alpha for the equation of the interface */
  double mm = MIN(fabs(mx),fabs(my));
  double V1 = .5*mm/(1-mm); double V2 =1.-V1;
  double alpha;
  if(cc0 < V1)    alpha = sqrt(2.*mm*(1.-mm)*cc0);
  else if(cc0 <= V2)    alpha = (1.-mm)*cc0 + 0.5*mm;
  else    alpha =  1. - sqrt(2.* (1.-cc0) *mm*(1.-mm));
  /* now back to the original line */
  alpha += MIN(0.,mx) + MIN(0.,my);
  return alpha;
}
//---------------------------------------------------------
// Volume
//  -------------------------
// -------------------------------------
/// Return Area over the [x0,1]x[y0,1] from (mx,my,alpha)
double MGSolCC::get_area(double mx,double my,double alpha,double x0,double y0) {
  alpha = alpha - my*y0 - mx*x0;  // move origin to (x0,y0)
  alpha = alpha + MAX(0.,-mx) + MAX(0.,-my);// rotate the figure
  mx = fabs(mx);  my = fabs(my);
  /* limit alpha so that 0. <= area <= 1. */
  alpha = MAX(0.,MIN(alpha,mx+my));
  if(MIN(mx,my) < MIN_VAL)  return alpha;
  return 0.5*(alpha*alpha-MAX(0.0,alpha-mx)*(alpha-mx)-
              MAX(0.0,alpha-my)*(alpha-my))/(mx*my);
}

// --------------------------------------------
// Volume computation
/// Area computation  over the [x0,dx]x[y0,dy] from (mx,my,alpha)
///
/// (1) move origin to (x0,y0)   ; (2) reflect rectangle;
/// (3) limit alpha (0<=al0<=0.5); (4) order coefficients (b1<=b2);
/// (5) get volume V
double MGSolCC::get_vol(double mx,double my,double mz,double alpha,double x0,double y0,double dx, double dy) {

  double al,b1,b2,V,tmp,al0;
  al = alpha-mx*x0-my*y0;//move origin to (x0,y0)
  // reflect rectangle
  b1 = mx*dx; b2=my*dy; al=al+MAX(0.,-b1) +MAX(0.,-b2);
  // limit alpha (0<=al0<=0.5)
  b1 = fabs(b1); b2=fabs(b2);  tmp = b1 + b2;
  al = MAX(0.,MIN(1.,al/tmp));  al0 = MIN(al,1.-al);
  // order coefficients (b1<=b2)
  b1 = MIN(b1,b2) /tmp; b2=1.-b1;
  // get volume V
  if(al0 < b1) tmp=0.5*al0*al0/ (b1*b2);
  else tmp = 0.5* (2.*al0-b1) /b2;
  if(al<0.5) return tmp*dx*dy; return (1.-tmp) *dx*dy;
}

double shape2(int i,double xi) {
  switch(i) {
  case 0:
    return .5*xi*(xi - 1.);
  case 2:
    return .5*xi*(xi + 1);
  case 1:
    return (1. - xi*xi);
  }
}
// ----------------------------------------------------
/// Return surface force tension
void MGSolCC::get_2surten(double xp[],double ff[],int ord[]) {

  for(int kf=0; kf<18; kf++)  ff[kf]=0.;
  int nxf=_nxyz[0][_Nlev_cc-1]-1; int nyf=_nxyz[1][_Nlev_cc-1]-1;
  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1;
  int fd=nxf/nxc;
  double pt1[3],pt2[3];
  int  ixc=(int)(xp[0]*nxc+0.5); int  jyc=(int)(xp[1]*nyc+0.5);
  const int indc=ixc+jyc*(nxc+1);

  for(int jj=-fd; jj<fd; jj++) {
    int jyf=(jyc*fd+jj);

    ExpRow(&c1[_Nlev_cc-1],_tempc1,_Nlev_cc-1,jyf,0);
    ExpRow(&_mx1[_Nlev_cc-1],_tempmx1,_Nlev_cc-1,jyf,0);
    ExpRow(&_my1[_Nlev_cc-1],_tempmy1,_Nlev_cc-1,jyf,0);

    for(int ii=-fd; ii<fd; ii++) {
      int ixf=(ixc*fd+ii);   double ccc=_tempc1[ixf+1];
      if(ccc >1.e-15  && ccc<1.-1.e-15) {
        double mx=_tempmx1[ixf+1]; double my=_tempmy1[ixf+1];
        double alpha=get_alpha2(mx,my,0.,ccc);
        get_pts(mx,my,alpha,pt1,pt2);
        double mom=sqrt(mx*mx+my*my);
        mx /=mom; my /=mom;
        // interface reconstruction points
        pt1[0] = (pt1[0]+ii)/(double)fd; pt1[1] = (pt1[1]+jj)/(double)fd;
        pt2[0] = (pt2[0]+ii)/(double)fd; pt2[1] = (pt2[1]+jj)/(double)fd;

        for(unsigned int js=0; js<3; js++) {
          for(unsigned int is=0; is<3; is++) {
            int inds=is+js*3;
            ff[inds] +=-my*(shape2(is,pt1[0])*shape2(js,pt1[1])-
                            shape2(is,pt2[0])*shape2(js,pt2[1]));
            ff[inds+9] +=mx*(shape2(is,pt1[0])*shape2(js,pt1[1])-
                             shape2(is,pt2[0])*shape2(js,pt2[1]));
          }
        }

      }
    }
  }
  // cartesian ordering
  for(int jsc=-1; jsc<2; jsc++) {
    for(int isc=-1; isc<2; isc++) {
      int inds=isc+1+(jsc+1)*3;
      ord[inds]=_invnode_dof[indc+isc+jsc*(nxc+1)];
    }
  }

  return;

}
// ----------------------------------------------------
/// Return phase
double MGSolCC::get_2phase(int blck,double xp[]) {

  int nxc=_nxyz[0][0]-1; int nyc=_nxyz[1][0]-1;
  int  ix=(int)(xp[0]*nxc+0.5); int  iy=(int)(xp[1]*nyc+0.5);
  const int ind=ix+iy*(nxc+1);
  double phase=0.;
  for(int ii=-blck; ii<blck; ii++)
    for(int jj=-blck; jj<blck; jj++) {
      int idl=ind+jj*(nxc+1)+ii+1;
      phase +=cc.Cmp[idl];
    }
//     if(phase/(4.*blck*blck) != phase/(4.*blck*blck))
//     {
//       int a=12;
//     }
  return phase/(4.*blck*blck);

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
  if(pt2[1] > 1.) { pt2[1] = 1.0; pt2[0] = (alpha-my) /mx;}
  pt1[0] = 1.; pt1[1] = (alpha-mx) /my;
  if(pt1[1] < 0.) { pt1[1] = 0.;  pt1[0] = alpha/mx;}

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

  return;
}

// --------------------------------------
// Advection
// --------------------------------------

// ======================================================================================
// Advection driver function
/// Driver for 2D split and unsplit algorithm.
/// The unsplit conserve mass exactly if div=0 exactly
void MGSolCC::Adv(
  unsigned int Level, ///< cc level
  double dt ,          ///< time step
  unsigned int itime
) {
  int ny=_nxyz[1][Level]-1;
  // Initial volume
  double sumc=0.;
  for(unsigned int iy=0; iy<ny; iy++)  {
    int len=M_GetLen(&c1[Level],iy+1);
    for(unsigned int kx=0; kx<len; kx++)  {
      sumc +=M__GetVal(&c1[Level],iy+1,kx);
    }
  }
  std::cout << _nxyz[0][Level]-1 << " x " << ny << " cc= " << sumc << std::endl;
#ifndef SPLIT
  // (babble)Unsplit volume advection
  AdvVelLXY(Level,dt/2.); // xy
  fprintf(stderr,"\n Bub UnSplit: xy ");
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_elv1);
  AdvVelLYX1(Level,dt/2.);
  fprintf(stderr,"\n Bub UnSplit: yx \n");
#else
  // Split volume advection
  lagrangeX(Level,dt/4.,itime);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_elv1);
  lagrangeY(Level,dt/4.);
  fprintf(stderr,"\n Split xy \n");
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_elv1);
  lagrangeX(Level,dt/4.);
  OldSol_update(_Nlev_cc-1,c1[_Nlev_cc-1],c1_old[_Nlev_cc-1]);
  GenNorm(_Nlev_cc-1,&MGSolCC::rec_elv1);
  lagrangeY(Level,dt/4.);
  fprintf(stderr,"\n  Split yx \n");
#endif



  return;
}

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// xy direction
void MGSolCC::AdvVelLXY(unsigned int Level, double dt) {
  // geometry level
  unsigned int nxl=_nxyz[0][Level]+1; unsigned int ny=_nxyz[0][Level]-1;
  unsigned int nxc=_nxyz[0][0]-1; unsigned int nyc=_nxyz[0][0]-1;
  double hx=1./(nxl-2);  double hy=1./ny;
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof;
  double uu[3]; double utemp[4]; int ixy[3];double disptemp[4];

  // First line
//   for(unsigned int ix=0; ix< nxl *3; ix++) {
//     _tempc1[ix]=0.; _tempc2[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
//   }
  for(unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }

  // Cartesian domain loop
  for(int jy=0; jy<ny; jy++) {

    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if(M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) >-1.e-10) {
      // column ix

      for(int ix=0; ix<nxl-2; ix++)  {
        ixy[0]=ix; ixy[1]=jy;
	get_disp(ixy,fd,disptemp);
//        std::cout<<"  "<<disptemp[0]<<"  "<<disptemp[1]<<"  "<<disptemp[2]<<"  "<<disptemp[3]<<std::endl;
        get_vel(ixy,fd,utemp);
        uu[0] =utemp[0]*stx; uu[2] =-utemp[1]*stx; // ul and -ur
        double vb=utemp[2]*sty; double vt =utemp[3]*sty;
// 	vb=0;vt=0;
        // xy  Advection and volume computation *********************
        // geometric parameters
        double du = -uu[2]-uu[0]; double aa = 1.-du; double ai=1./aa;
        double dyt = MAX(0.,vt); double dyb = -MIN(0.,vb);
        double y0 = MAX(0.,vb);  double  y1 = MIN(1.,1.+vt);
        double dyc = y1 - y0;
        double x0 = MAX(0.,-uu[0]);   double x1 = MIN(1.,1.+uu[2]);
        uu[1] = x1 - x0;
        double x6 = MAX(0.0,uu[0]*ai); //double x7 = MIN(1.0,1.0 + (-uu[2]*ai));

        // loop left(0) center(1) right(2) cell
        for(int ii=0; ii<3; ii++) {
          // color function
          int indl= ix+ii-1+(jy%3)*nxl; double cc1=_tempc1[indl+1];
          if(cc1 > 0. && uu[ii]>0.) {
            // set up volume computations
            double vcold = uu[ii]; double dx = vcold*ai;
            double vct = dx*dyt; double vcb = dx*dyb;
            double vcc = dx*dyc;
	    
	    
//               //LUCA OUTLET BOUNDARY
//               if(
// 		bc_outlet_read(ix,jy,indl,hx,hy,cc1,0) //2nd cell near boundary
// 	      ){
// 		  if(cc1 > 0.) {
// 			cc1=0.01;
// 			_tempc1[indl+1]=0.01;
// 			if(
// 			  bc_outlet_read(ix,jy,indl,hx,hy,cc1,1) //last cell near boundary
// 			){
// 			  cc1=0.;
// 			  _tempc1[indl+1]=0.;
// 			}
// 		      
// 			//normali manuali //need tests
// // 			if(ix*hx <0.0+1*hx ||
// // 			  ix*hx >1.0-1*hx		    
// // 			){
// // 			  _tempmx1[indl+1]=0.;
// // 			  _tempmy1[indl+1]=1.;
// // 			}
// // 			if( jy*hy <0.0+1*hy ||
// // 			  jy*hy >1.0-1*hy
// // 			){
// // 			  _tempmx1[indl+1]=1.;
// // 			  _tempmy1[indl+1]=0.;
// // 			}
// 		  }
// 		  
// 		}
// 	      //LUCA END
// 	      //LUCA INTERNAL BOUNDARY CONDITIONS
// // 	      internal_read(_tempc1,ix, jy,indl, hx, hy, cc1);
// 	    0.5,0.35,0.1
	    
// 	      if((ix*hx-0.5)*(ix*hx-0.5)+(jy*hy-0.85)*(jy*hy-0.85)<0.1*0.1){
//   // 		
// 		      cc1=1.;
// 		      _tempc1[indl+1]=1.;
// 		
// 	      }
	      //LUCA END
            if(cc1 < 1.) {
              // unit normal and alpha
              double mx=_tempmx1[indl+1];  double my=_tempmy1[indl+1];
              double alpha=get_alpha2(mx,my,0.,cc1);
              // computation of the volume fraction
              double w1=0.5*(1.-uu[0])*(2-ii)*(1-ii)+x0*(2-ii)*(ii);
              vcold = get_vol(mx,my,0.,alpha,w1,0.,uu[ii],1.0);
              vct = vcb = vcc = 0.;
              if(vcold > 0.) {
                double w2=0.5*(1.-dx)*(ii-1)*(ii)+x6*(2-ii)*(ii);
                my *= ai;  alpha += mx*(uu[0]+ii-1.) + my*vb;  mx *= aa;
                if(vt > 0.) vct += get_vol(mx,my,0.,alpha,w2,1.,dx,dyt);
                if(vb < 0.) vcb += get_vol(mx,my,0.,alpha,w2,vb,dx,dyb);
                vcc +=  get_vol(mx,my,0.,alpha,w2,y0,dx,dyc);
              }
            }
            
	      
	      
            // storage of the volume fraction
            _tempc2[ix+((jy+1)%3)*nxl+1] += vct; _tempc2[ix+((jy)%3)*nxl+1] += vcc;
            _tempc2[ix+((jy+2)%3)*nxl+1] += vcb;
          }
        }

      } // global grid ix
    } // if len>0
    // Contracting tempc2 -> c1
    if(jy>1) CtrRow(_tempc2,&c1[Level],Level,jy-1,indup);
    // Expanding matrix c1_old jy+2 ->line indup  _tempc1
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // Last line contraction
  CtrRow(_tempc2,&c1[Level],Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// xy direction
void MGSolCC::AdvVelLXY1(unsigned int Level, double dt) {
  // geometry level
  unsigned int nxl=_nxyz[0][Level]+1; unsigned int ny=_nxyz[0][Level]-1;
  unsigned int nxc=_nxyz[0][0]-1; unsigned int nyc=_nxyz[0][0]-1;
  double hx=1./(nxl-2);  double hy=1./ny;
  int n_nodes= (nxc+1) * (nyc+1);
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof; double utemp[4]; int ixy[3];double disptemp[4];

//   // First line
//   for(unsigned int ix=0; ix< nxl *3; ix++) {
//     _tempc1[ix]=0.; _tempc2[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
//   }
  for(unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }
  // Cartesian domain loop
  for(int jy=0; jy<ny; jy++) {

    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if(M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) >-1.e-10) {
      // column ix
      for(int ix=0; ix<nxl-2; ix++)  {
        // Color function stensil  left top center bottom right
        int indc=ix+(jy%3)*nxl;    double ccc=_tempc1[indc+1];
        int indt=ix+((jy+1)%3)*nxl; int indb=ix+indup; // not use
        int indl=ix-1+(jy%3)*nxl;  double ccl=_tempc1[indl+1];
        int indr=ix+1+(jy%3)*nxl;  double ccr=_tempc1[indr+1];

        // velocity field  ul(left),ur(right),vt(top),vb(bot)
        ixy[0]=ix; ixy[1]=jy;
     get_disp(ixy,fd,disptemp);
       std::cout<<"  "<<disptemp[0]<<"  "<<disptemp[1]<<"  "<<disptemp[2]<<"  "<<disptemp[3]<<std::endl;
        get_vel(ixy,fd,utemp);
        double ul =utemp[0]*stx; double ur =utemp[1]*stx; // ul and ur
        double vb=utemp[2]*sty; double vt =utemp[3]*sty;

        // xy  Advection and volume computation *********************
        double vcold, dx,vct,vcb,vcc, x0,x1,dxold ;
        double du = ur - ul;  double    aa = 1. - du;   double   ai = 1./aa;
        double dyt = MAX(0.,vt);  double    dyb = - MIN(0.,vb);
        double  y0 = MAX(0.,vb);   double    y1 = MIN(1.,1.+vt);     double  dyc = y1 - y0;

        // left side ----------------------------
        if(ccl > 0. && ul > 0.) {
          vcold = ul;	 dx = ul*ai;
          vct = dx*dyt; vcb = dx*dyb; vcc = dx*dyc;
          if(ccl < 1.) {
            // unit normal
            double mx=_tempmx1[indl+1]; double my=_tempmy1[indl+1];
            double alpha=get_alpha2(mx,my,0.,ccl);
            // Volume
            vcold = get_vol(mx,my,0.,alpha,1.-ul,0.,ul,1.);
            vct = vcb = vcc = 0.;
            if(vcold > 0.) {
              my *= ai;  alpha += mx*(ul-1.) + my*vb;  mx *= aa;
              if(vt > 0.)  vct = get_vol(mx,my,0.,alpha,0.,1.,dx,dyt);
              if(vb < 0.)  vcb = get_vol(mx,my,0.,alpha,0.,vb,dx,dyb);
              vcc =  get_vol(mx,my,0.,alpha,0.,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct; _tempc2[indc+1] += vcc; _tempc2[indb+1] += vcb;
        }
        // right side ---------------------------------
        if(ccr > 0. && ur < 0.) {
          vcold = -ur;	dx = -ur*ai;
          vct = dx*dyt;	vcb = dx*dyb;	vcc = dx*dyc;
          if(ccr < 1.) {
            double mx=_tempmx1[indr+1]; double my=_tempmy1[indr+1];
            double alpha=get_alpha2(mx,my,0.,ccr);
            // Volume
            vcold = get_vol(mx,my,0.,alpha,0.,0.,-ur,1.);
            vct = vcb = vcc = 0.;
            if(vcold > 0.) {
              my *= ai; alpha += mx*(ul+1.) + my*vb;   mx *= aa;
              if(vt > 0.) vct = get_vol(mx,my,0.,alpha,1.-dx,1.,dx,dyt);
              if(vb < 0.) vcb = get_vol(mx,my,0.,alpha,1.-dx,vb,dx,dyb);
              vcc =  get_vol(mx,my,0.,alpha,1.-dx,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct; _tempc2[indc+1] += vcc; _tempc2[indb+1] += vcb;
        }
        // center ----------------------------------------
        if(ccc > 0.) {  // c > 0.
          x0 = MAX(0.,-ul); x1 = MIN(1.,1.-ur);
          dxold = x1 - x0; vcold = dxold;	dx = dxold*ai;
          vct = dyt*dx;	vcb = dyb*dx;	vcc = dyc*dx;
          if(ccc < 1.) {  // c < 1.
            double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
            double alpha=get_alpha2(mx,my,0.,ccc);
            // Volume
            vcold = get_vol(mx,my,0.,alpha,x0,0.0,dxold,1.0);
            vct = vcb = vcc = 0.0;
            if(vcold > 0.0) {
              my *= ai; alpha += mx*ul + my*vb; mx *= aa;
              x0 = MAX(0.0,ul*ai); x1 = MIN(1.0,1.0 + (ur*ai));
              if(vt > 0.) vct = get_vol(mx,my,0.,alpha,x0,1.0,dx,dyt);
              if(vb < 0.0) vcb = get_vol(mx,my,0.,alpha,x0,vb,dx,dyb);
              vcc =  get_vol(mx,my,0.,alpha,x0,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct; _tempc2[indc+1] += vcc; _tempc2[indb+1] += vcb;
        }
      } // global grid ix
    } // if
    // Contraction and expansion  in  indup line
    if(jy>1) CtrRow(_tempc2,&c1[Level],Level,jy-1,indup);
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // final line
  CtrRow(_tempc2,&c1[Level],Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// yx direction
void MGSolCC::AdvVelLYX1(unsigned int Level, double dt) {
  // geometry level
  unsigned int nxl=_nxyz[0][Level]+1; unsigned int ny=_nxyz[0][Level]-1;
  unsigned int nxc=_nxyz[0][0]-1; unsigned int nyc=_nxyz[0][0]-1;
  double hx=1./(nxl-2);  double hy=1./ny;
  int n_nodes= (nxc+1) * (nyc+1);
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof; double utemp[4]; int ixy[3];

  // First line
  for(unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
  }
  for(unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }

  // Cartesian domain
  for(int jy=0; jy<ny; jy++) {
    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if(M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) >-1.e-10) {
      // column ix
      for(int ix=0; ix<nxl-2; ix++)  {
        // Color function stensil  left top center bottom right
        int indc=ix+(jy%3)*nxl;    double ccc=_tempc1[indc+1];
        int indt=ix+((jy+1)%3)*nxl; double cct=_tempc1[indt+1];
        int indb=ix+indup;         double ccb=_tempc1[indb+1];
        int indl=ix-1+(jy%3)*nxl;  int indr=ix+1+(jy%3)*nxl;

        // velocity field  ul(left),ur(right),vt(top),vb(bot)
        ixy[0]=ix; ixy[1]=jy;
        get_vel(ixy,fd,utemp);
        double ul =utemp[0]*stx; double ur =utemp[1]*stx; // ul and -ur
        double vb=utemp[2]*sty; double vt =utemp[3]*sty;

        // Volume computations
        double du,dv,aa,ai,y0,y1,x0,x1,dx,dy;
        double dyc,dxc,dyt,dyb,dxr,dxl,vcold,vct,vcb,vcr,vcl,vcc;
        double dxold,dyold;

        // Advection yx
        dv = vt - vb;      aa = 1. - dv;     ai = 1./aa;
        dxr = MAX(0.0,ur);     dxl = - MIN(0.0,ul);
        x0 = MAX(0.0,ul);      x1 = MIN(1.0,1.0+ur);    dxc = x1 - x0;
        // bottom side -------------------------------
        if(ccb > 0. && vb > 0.) {
          vcold = vb;	dy = vb*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;

          if(ccb < 1.) {
            // unit normal
            double mx=_tempmx1[indb+1]; double my=_tempmy1[indb+1];
            double alpha=get_alpha2(mx,my,0.,ccb);
            // Volume
            vcold = get_vol(mx,my,0.,alpha,0.,1.-vb,1.,vb);
            vcr = vcl = vcc = 0.;
            if(vcold > 0.) {
              mx *= ai;  alpha += my*(vb-1.) + mx*ul;  my *= aa;
              if(ur > 0.) vcr = get_vol(mx,my,0.,alpha,1.,0.,dxr,dy);
              if(ul < 0.) vcl = get_vol(mx,my,0.,alpha,ul,0.,dxl,dy);
              vcc =  get_vol(mx,my,0.,alpha,x0,0.,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc;
          _tempc2[indr+1] += vcr;
        }
        // top -----------------------------
        if(cct > 0. && vt < 0.) {
          vcold = -vt;	dy = -vt*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;
          if(cct < 1.) {
            // unit normal
            double mx=_tempmx1[indt+1]; double my=_tempmy1[indt+1];
            double alpha=get_alpha2(mx,my,0.,cct);
            // Volume
            vcold = get_vol(mx,my,0.,alpha,0.,0.,1.,-vt);
            vcr = vcl = vcc = 0.;
            if(vcold > 0.) {
              mx *= ai;  alpha += my*(vb+1.0) + mx*ul;   my *= aa;
              if(ur > 0.) vcr = get_vol(mx,my,0.,alpha,1.,1.-dy,dxr,dy);
              if(ul < 0.) vcl = get_vol(mx,my,0.,alpha,ul,1.-dy,dxl,dy);
              vcc =  get_vol(mx,my,0.,alpha,x0,1.-dy,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc; _tempc2[indr+1] += vcr;
        }
        // center  --------------------------------------
        if(ccc > 0.) {
          y0 = MAX(0.,-vb);	y1 = MIN(1.,1.-vt);
          dyold = y1 - y0;	vcold = dyold;	dy = dyold*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;
          if(ccc < 1.) {
            // unit normal
            double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
            double alpha=get_alpha2(mx,my,0.,ccc);
            // volume
            vcold = get_vol(mx,my,0.,alpha,0.,y0,1.,dyold);
            vcr = vcl = vcc = 0.;
            if(vcold > 0.) {
              mx *= ai;	    alpha += my*vb + mx*ul;    my *= aa;
              y0 = MAX(0.,vb*ai);	    y1 = MIN(1.,1. + (vt*ai));
              if(ur > 0.)     vcr = get_vol(mx,my,0.,alpha,1.,y0,dxr,dy);
              if(ul < 0.)    vcl = get_vol(mx,my,0.,alpha,ul,y0,dxl,dy);
              vcc =  get_vol(mx,my,0.,alpha,x0,y0,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc; _tempc2[indr+1] += vcr;
        }
      }
    }
    // storage in indup line
    if(jy>1) CtrRow(_tempc2,&c1[Level],Level,jy-1,indup);
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // last line
  CtrRow(_tempc2,&c1[Level],Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ---------------------------------------------------------------
// Split advection along the x-direction:
/// Standard advection in one direction by using the normal from
/// (mx1 my1) previously computed by Young function
void  MGSolCC::lagrangeX(const unsigned int Level,
                         const double dt, 
			 unsigned int itime) {

  // Set up
  const unsigned int nxl=_nxyz[0][Level]+1;  const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nxc=_nxyz[0][0]-1;  const unsigned int nyc=_nxyz[1][0]-1;
  const unsigned int n_nodes=(nxl-1)*(ny+1);
  const double hx=1./(nxl-2); const double hy=1./ny;
  const int fd=(nxl-2)/nxc;
  double stx=dt/hx; double utemp[4]; int ixy[2];

  // First line
  for(unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
  }
  for(unsigned int j0=0; j0<2; j0++) {
    int indl=j0*nxl;
    ExpRow(&c1_old[Level],_tempc1,Level,j0,indl);
    ExpRow(&_mx1[Level],_tempmx1,Level,j0,indl);
    ExpRow(&_my1[Level],_tempmy1,Level,j0,indl);
  }

  // Domain loop (i,j)=[0,nz]x[0,ny]
  for(int j=0; j<ny ; j++) {

    // Update temporary line -> j-1 && solved line -> j
    int indup=((j+2)%3)*nxl; int indout=(j%3)*nxl;
    // index i loop
    for(int i=0; i<nxl-1 ; i++) {
      // color function
      int indc=i+(j%3)*nxl; const double  ccc=_tempc1[indc+1];

      // velocity field  ul(left),ur(right),vt(top),vb(bot)
      ixy[0]=i; ixy[1]=j;
      get_vel(ixy,fd,utemp);
      double a1 =utemp[0]*stx; double a2 =utemp[1]*stx; // ul and -ur

      // Volume computation ******************************
      double vof1 = 0.; double  vof2 = 0.; double vof3 = 0.; // c=0.
      if(ccc > 1.-1e-16) {   // c=1.
        vof1 =ccc* MAX(-a1,0.);    vof3 =ccc* MAX(a2,0.);
        vof2 =ccc* (1.- MAX(a1,0.)-MAX(-a2,0.));
      } else if(ccc > 0.) {  // 0.< c <1.
        // unit normal
        double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
        // intersection alpha
        double alpha=get_alpha2(mx,my,0.,ccc);
        // Advected Volume
        mx=mx/(1.-a1+a2); alpha +=mx*a1;
        if(a1<0.) vof1=get_vol(mx,my,0.,alpha,a1,0.,-a1,1.);
        if(a2>0.) vof3=get_vol(mx,my,0.,alpha,1.,0.,a2,1.);
        vof2=get_vol(mx,my,0.,alpha,MAX(0.,a1),0.,1-MAX(0.,-a2)-MAX(0.,a1),1.);
      }
      // temporary storage
      _tempc2[indc+1]+=vof2; _tempc2[indc]+=vof1; _tempc2[indc+2]+=vof3;
    }
    // Contracting _tempc2 -> c1
    CtrRow(_tempc2,&c1[Level],Level,j,indout);
    //Expanding matrix c1_old j+2->line indup _tempc1
    ExpRow(&c1_old[Level],_tempc1,Level,j+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,j+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,j+2,indup);
  }
  return;
}

// ---------------------------------------------------------------
// Split advection along the y-direction
/// Standard advection in one direction by using the normal from
/// (mx1,my1) previously computed by Young function.
/// The matrices mx1,my1 are in contracted form
void  MGSolCC::lagrangeY(const unsigned int Level,
                         const double dt, 
			 unsigned int itime) {

  const unsigned int nxl=_nxyz[0][Level]+1;  const unsigned int ny=_nxyz[1][Level]-1;
  const unsigned int nxc=_nxyz[0][0]-1;  const unsigned int nyc=_nxyz[1][0]-1;
  const unsigned int n_nodes=(nxc+1)*(nyc+1);
  const double hx=1./(nxl-2); const double hy=1./ny;
  const int fd=(nxl-2)/nxc;
  double sty=dt/(hy); double utemp[4]; int ixy[2];

  // First line
  for(unsigned int ix=0; ix< nxl*3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.; _tempmx1[ix]=0.; _tempmy1[ix]=0.;
  }
  for(unsigned int j0=0; j0<2; j0++) {
    ExpRow(&c1_old[Level],_tempc1,Level,j0,j0*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,j0,j0*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,j0,j0*nxl);
  }

  // Cartesian domain loop  (i,j)=[0,nx]x[0,ny]
  for(int j=0 ; j<ny ; j++) {

    // Update temporary line -> j-1 && solved line -> j
    int indup=((j+2)%3)*nxl;  int indout=((j-1)%3)*nxl;
    // index i loop
    for(int i=0; i<nxl-2 ; i++) {
      // color function
      const int indc=i+(j%3)*nxl; const double  ccc=_tempc1[indc+1];
      const int indb=i+((j+2)%3)*nxl; const int indt=i+((j+1)%3)*nxl;
      ixy[0]=i; ixy[1]=j;
      // velocity field  ul(left),ur(right),vt(top),vb(bot) -> utemp
      get_vel(ixy,fd,utemp);
      double a1=utemp[2]*sty; double a2 =utemp[3]*sty;  //vb vt

      // advected volume ---------------------------------
      double vof1=0.; double vof2=0.; double vof3=0.; // c=0.
      if(ccc >= 1.) {       // c=1.
        vof1 =ccc* MAX(-a1,0.); vof3 =ccc* MAX(a2,0.);
        vof2 =ccc* (1. - MAX(a1,0.) - MAX(-a2,0.));
      } else if(ccc > 0.) { // 0.< c <1.
        // unit normal
        double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
        // intersect alpha
        double alpha=get_alpha2(mx,my,0.,ccc) ;
        // advected volume
        my=my/(1.-a1+a2); alpha +=my*a1;
        if(a1<0.)  vof1=get_vol(my,mx,0.,alpha,a1,0.,-a1,1.);
        if(a2>0.)  vof3=get_vol(my,mx,0.,alpha,1.,0.,a2,1.);
        vof2=get_vol(my,mx,0.,alpha,MAX(0.,a1),0.,1-MAX(0.,-a2)-MAX(0.,a1),1.);
      }
      // new values of c1
      _tempc2[indb+1]+=vof1; _tempc2[indc+1]+=vof2; _tempc2[indt+1] +=vof3;
    }

    // Contracting  tempc2[((j-1)%3)*nxl] -> c2[j-1] --------------
    if(j>0) CtrRow(_tempc2,&c1[Level],Level,j-1,indout);
    // Expanding matrix c1_old[j+2] ->line  _tempc1[((j+2)%3)*nxl]
    ExpRow(&c1_old[Level],_tempc1,Level,j+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,j+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,j+2,indup);
  }
  // storage of the last line
  CtrRow(_tempc2,&c1[Level],Level,ny-1,((ny-1)%3)*nxl);
  return;
}



// ------------------------------------------------------
// Utility for contracting and expanding compressed Matrix
// ------------------------------------------------------
// ----------------------------------
/// Row expansion  c1_old -> new line _tempc1
void MGSolCC::ExpRow(Matrix * mtrx,double tempc1[],unsigned int Level,unsigned int jy,unsigned int rowl) {
  int nxl=_nxyz[0][Level]; int ny=_nxyz[1][Level]-1;
  // Expanding matrix c1 -> new line _tempc1
  for(unsigned int ix=1; ix<nxl; ix++) tempc1[ix+rowl]=0.;
  if(jy < ny && jy>= 0) {  // no ghost cell
    for(unsigned int kx=0; kx<M_GetLen(mtrx,jy+1); kx++)  {
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
    if(tempc2[ind] > 1.-MIN_VAL) {
      if(flag1 == 0) icount1++;
      flag1 = 1; tempc2[ind]=1.;
    } else {
      flag1 = 0; if(tempc2[ind]> MIN_VAL) icount++;
      else tempc2[ind]=0.;
    }
  }
  // compressed storage matrix c1
  M_SetLen(mtrx,jy+1,icount+icount1);
  int irow1=0; int count1=0;
  for(unsigned int ixf=1; ixf<nxl; ixf++) {
    int ind=ixf+rowl;
    if(tempc2[ind]> 0.) {
      if(tempc2[ind]< 1.) {
        M__SetEntry(mtrx,jy+1,irow1,ixf,tempc2[ind]);
        irow1++;
      } // <1
      else {// _tempc2[ind]== 1.
        count1++;
        if(tempc2[ind+1]<1.) {
          M__SetEntry(mtrx,jy+1,irow1,M__GetPos(mtrx,jy+1,irow1-1)+1,count1);
          irow1++;
          count1=0;
        }
      } // else
      tempc2[ind]=0.;
    } // .>0
  }
}
// --------------------------------------------------------

// void MGSolCC::InitVel(const MGMeshC & mgmesh/*,const MGSol & mgs*/,const double dt) {
//   V_Destr(&_uvw);
// // #ifdef NS_EQUATIONS
// //   _uvw=mgs.uw[mgs._NoLevels-1];
// // #else
//   int Level=mgmesh._NoLevels -1;
//   V_Constr(&_uvw,(char *)"_uvw",2*_nxyz[0][Level]*_nxyz[1][Level],Normal,_LPTrue);
//   GenVel(mgmesh,_uvw,dt);
// // #endif
//   return;
// }

// // ======================================================================================
// /// This function generates the velocity fields (Cartesian mesh)
// void MGSolCC::GenVel(
//    const int Level,
// //   const MGMeshC & mgmesh,   ///<   mesh
// //   QVector &sol,             ///<   velocity field
//   const double dt           ///<   time step
// ) {
// 
//   // Box elements and space steps
//   int nxyz[3];nxyz[0]=_nxyz[0][Level];nxyz[1]=_nxyz[1][Level];nxyz[2]=0;
//   double hxyz[3];hxyz[0]=1./(nxyz[0]-1);hxyz[1]=1./(nxyz[1]-1);hxyz[2]=0;
// #if DIMENSION==3
//    nxyz[2]=_nxyz[2][Level];hxyz[2]=1./(nxyz[2]-1)];
// #endif
//     int  offset=1; 
//     for(int idim=0;idim<DIMENSION;idim++) offset *=nxyz[idim];
//     // velocity vector _uvw
//     V_Destr(&_uvw);
//     V_Constr(&_uvw,(char *)"uvw",2*offset,Normal,_LPTrue);
// 
//   // Dof
//   const unsigned int n_u_dofs = NDOF_P;//NDOF_FEM; const unsigned int n_p_dofs = NDOF_P;
// 
// 
// //   double xi[NDOF_P];  double yj[NDOF_P]; int dofu[NDOF_P];
//   double u_value[3]; int dof_u[3];//double Real v_value =0.;
//   for(int idim=0;idim<DIMENSION;idim++) {
//     u_value[idim]=0.;
//     dof_u[idim]=0;
//   }
//   for(int j=0; j<nxyz[1]; j++) {
//     for(int i=0; i<nxyz[0]; i++) {
// 
//       // vortex
//     
//       u_value[1] = 1*(PSI(i*hxyz[0]+hxyz[0],j*hxyz[1])-PSI(i*hxyz[0],j*hxyz[1]))/hxyz[0];
//       u_value[0] = -1*(PSI(i*hxyz[0],j*hxyz[1]+hxyz[1])-PSI(i*hxyz[0],j*hxyz[1]))/hxyz[1];
//       //rotation --------------------------
// //          int    dof_u=i+j*nx+1;       int  dof_v=dof_u+offset;
// //   const Real u_value=-3.14159265358979*(j*hy-0.5); const Real v_value=3.14159265358979*(i*hx-0.5);
// //   const Real norm_u=sqrt(u_value*u_value+v_value*v_value);
//       // translation -------------------------
// //         u_value =0.5; v_value=0.;
//       for(int idim=0;idim<DIMENSION;idim++) {
//         dof_u[idim]=i+j*nxyz[0]+1+idim*offset;  //dof_u[1]=dof_u[0]+offset;
//       V_SetCmp(&_uvw,dof_u[idim],u_value[idim]);
// //       V_SetCmp(&_uvw,dof_u[1],u_value[0]);
//       }
//     }
//   }
// //   }
// 
// //   // loop element --------------------------------------
// //   for (unsigned int iel=0 ; iel <n_elem ; ++iel) {
// //     unsigned int elem_gidx=iel*n_u_dofs;
// //     // the local nodes
// //     for (unsigned int i=0; i<n_u_dofs; i++) {
// //       // coordinates
// //       int k=map_nodes[elem_gidx+i];
// //       const Real xi = xyz_glob[k];  const Real yj = xyz_glob[k+offset];
// //       // dofs
// //       unsigned int dof_u=k+1; unsigned int dof_v=k+offset+1;
// //       // Set the initial velocity (note:: nondimensional field)
// //
// // // translation -------------------------
// //       //   const Real u_value =0.;const Real v_value=1.;
// // // ------------------
// // //rotation --------------------------
// // // const Real u_value=3.14159265358979*(yj-0.5); const Real v_value=-3.14159265358979*(xi-0.5);
// // // --------------------------------------------
// // // vortex ---------------------------
// //       const Real v_value =  ((double) nx)*(PSI(xi+hx,yj)-PSI(xi,yj));
// //       const Real u_value = - ((double) nx)*(PSI(xi,yj+hx)-PSI(xi,yj));
// // // ---------------------------------
// //
// // // -------------------------------------------------
// //       // set velocity field
// //       V_SetCmp(&sol,dof_u,u_value);   V_SetCmp(&sol,dof_v,v_value);
// //
// //     } // --------------------------------------------------------
// //   } // end of element loop
//   return;
// }


void MGSolCC::SourceVel_ext(
//   const MGMeshC & mgmesh,
  double _uvw_field[],
  const double dt) {

// mesh structure --------------------------------------------------------------------
  double  hxyz[3]; hxyz[2]=0;  // space step
  int nxyz[3]= {NX,NY,NZ};     // number of divisions
  int n_cell=1; int nNodes=1;
  for(int idim=0; idim<DIMENSION; idim++) {
    hxyz[idim] = 1./nxyz[idim]; n_cell *=nxyz[idim];
    nNodes *=(nxyz[idim]+1);
  }
  
//   int Level=mgmesh._NoLevels -1;
  V_Destr(&_uvw); V_Constr(&_uvw,(char *)"_uvw",DIMENSION*nNodes,Normal,_LPTrue);
//   int offset=_nxyz[0][Level]*_nxyz[1][Level];
  for(int idim=0;idim<DIMENSION;idim++){
    for(int inode=0;inode<nNodes;inode++){
      V_SetCmp(&_uvw,inode,_uvw_field[inode+idim*nNodes]); 
    }
  }
//     V_SetCmp(&sol,dof_u,u_value); V_SetCmp(&sol,dof_v,v_value);
//   GenVel(mgmesh,_uvw,dt);
// // #endif
//   return;
// }
// 
// // ======================================================================================
// /// This function generates the velocity fields (Cartesian mesh)
// void MGSolCC::GenVel(
//   const MGMeshC & mgmesh,   ///<   mesh
//   QVector &sol,             ///<   velocity field
//   const double dt           ///<   time step
// ) {

// Get a constant reference to the MG solution
//   int Level=mgmesh._NoLevels -1;
// //   const unsigned int nx=_nxyz[0][0]; const unsigned int ny=_nxyz[1][0]; // number of nodes
// //   const unsigned int nxf=_nxyz[0][_Nlev_cc-1];                     //
// //   const double hx=1./(nx-1); const double hy=1./(ny-1);       //  grid dx dy
// //   const double hxf=1./(nxf-1);
// 
//   const unsigned int nx=_nxyz[0][Level]; const unsigned int ny=_nxyz[1][Level]; // number of nodes
// //   const unsigned int nxf=_nxyz[0][_Nlev_cc-1];                     //
//   const double hx=1./(nx-1); const double hy=1./(ny-1);       //  grid dx dy
// //   const double hxf=1./(nxf-1);
//   // Get a constant reference to the mesh objects.
//   const int *map_nodes=mgmesh._elem_map[Level];
//   const double *xyz_glob=mgmesh._xyz;
//   const  unsigned int  offset=mgmesh._NoNodes[Level];
//   const unsigned int  n_elem=mgmesh._NoElements[Level];
// 
//   // Dof
//   //int *idx_dof=_node_dof[Level];
//   const unsigned int n_u_dofs = NDOF_P;//NDOF_FEM; const unsigned int n_p_dofs = NDOF_P;
// 
// 
//   double xi[NDOF_P];  double yj[NDOF_P]; int dofu[NDOF_P];
//   Real u_value =0.; Real v_value =0.;
//   for(int j=0; j<ny; j++) {
//     for(int i=0; i<nx; i++) {
// //       int ielem= i+j*(nx-1);
//       // connectivity cartesian cell    quad4
// //       xi[0]=i*hx;      yj[0]= j*hy;     dofu[0]=i+j*nx;
// //       xi[1]=(i+1)*hx;  yj[1] = j*hy;    dofu[1]=(i+1)+j*nx;
// //       xi[2]=(i+1)*hx;  yj[2] = (j+1)*hy;dofu[2]=(i+1)+(j+1)*nx;
// //       xi[3]=(i)*hx;    yj[3] = (j+1)*hy;dofu[3]=i+j*nx;
// 
// 
// //       for(unsigned int i=0; i<n_u_dofs; i++) {
// 
// //         const Real    dof_u=dofu[i]+1;      const Real  dof_v=dof_u+offset;
// 
//       // vortex
//       int    dof_u=i+j*nx+1;       int  dof_v=dof_u+offset;
//       const Real v_value = 1*((double)(nx-1))*(PSI(i*hx+hx,j*hy)-PSI(i*hx,j*hy));
//       const Real u_value = -1*((double)(ny-1))*(PSI(i*hx,j*hy+hy)-PSI(i*hx,j*hy));
//       //rotation --------------------------
// //          int    dof_u=i+j*nx+1;       int  dof_v=dof_u+offset;
// //   const Real u_value=-3.14159265358979*(j*hy-0.5); const Real v_value=3.14159265358979*(i*hx-0.5);
// //   const Real norm_u=sqrt(u_value*u_value+v_value*v_value);
//       // translation -------------------------
// //         u_value =0.5; v_value=0.;
//       V_SetCmp(&sol,dof_u,u_value); V_SetCmp(&sol,dof_v,v_value);
// 
//     }
//   }
// //   }
// 
// //   // loop element --------------------------------------
// //   for (unsigned int iel=0 ; iel <n_elem ; ++iel) {
// //     unsigned int elem_gidx=iel*n_u_dofs;
// //     // the local nodes
// //     for (unsigned int i=0; i<n_u_dofs; i++) {
// //       // coordinates
// //       int k=map_nodes[elem_gidx+i];
// //       const Real xi = xyz_glob[k];  const Real yj = xyz_glob[k+offset];
// //       // dofs
// //       unsigned int dof_u=k+1; unsigned int dof_v=k+offset+1;
// //       // Set the initial velocity (note:: nondimensional field)
// //
// // // translation -------------------------
// //       //   const Real u_value =0.;const Real v_value=1.;
// // // ------------------
// // //rotation --------------------------
// // // const Real u_value=3.14159265358979*(yj-0.5); const Real v_value=-3.14159265358979*(xi-0.5);
// // // --------------------------------------------
// // // vortex ---------------------------
// //       const Real v_value =  ((double) nx)*(PSI(xi+hx,yj)-PSI(xi,yj));
// //       const Real u_value = - ((double) nx)*(PSI(xi,yj+hx)-PSI(xi,yj));
// // // ---------------------------------
// //
// // // -------------------------------------------------
// //       // set velocity field
// //       V_SetCmp(&sol,dof_u,u_value);   V_SetCmp(&sol,dof_v,v_value);
// //
// //     } // --------------------------------------------------------
// //   } // end of element loop
  return;
}



#endif
#endif
