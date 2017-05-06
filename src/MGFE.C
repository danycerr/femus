// std library -----------------
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// class headers ----------
#include "MGFE.h"           // class header
#include "MGFE_conf.h"       // configuration

// local ------------------
#include "MGUtils.h"        // data and path to file
#include "Printinfo_conf.h" // to print info


// ==============================
/// This function constructs the class MGFE:
MGFE::MGFE(
//   const MGUtils& mgutils_in, // file names and path input
  const int order_in,        // Fem order input
  const int fem_type         // Fem type: 27(Quad,Hex) or 10(Tri,Tet)
) : //_mgutils(mgutils_in),    // file names and path
    _dim(DIMENSION),         // space dimension
    _order(order_in)         // Fem order
{// =============================

  int n_shape[3];     
  int deg[3];          
  int n_gauss[3];

if (order_in == 2)
switch(fem_type) {
  case 27:
   n_shape[0]=3; n_shape[1]=9; n_shape[2]=27;   // hex quadratic  shapes
   deg[0]=2;     deg[1]=4;     deg[2]=8;        // degree of the shape
   n_gauss[0]=3; n_gauss[1]=9; n_gauss[2]=27;   // 3x3x3 gaussian points
  break;
  case 10:
   n_shape[0]=3; n_shape[1]=6; n_shape[2]=10;   // tet quadratic  shapes
   deg[0]=2;     deg[1]=2;     deg[2]=3;        // degree of the shape
   n_gauss[0]=3; n_gauss[1]=7; n_gauss[2]=14;   // 3x3x3 gaussian points
  break;
  default :
    std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
  break;
}
else if (order_in == 1)
switch(fem_type) {
  case 27:
   n_shape[0]=2; n_shape[1]=4; n_shape[2]=8;   // hex linear  shapes
   deg[0]=1;     deg[1]=2;     deg[2]=3;        // degree of the shape
   n_gauss[0]=3; n_gauss[1]=9; n_gauss[2]=27;   // 3x3x3 gaussian points
  break;
  case 10:
   n_shape[0]=2; n_shape[1]=3; n_shape[2]=4;   // tet linear  shapes
   deg[0]=1;     deg[1]=2;     deg[2]=3;        // degree of the shape
   n_gauss[0]=1; n_gauss[1]=7; n_gauss[2]=14;   // 3x3x3 gaussian points
  break;
  default :
    std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
  break;
}
else if (order_in == 0)
switch(fem_type) {
  case 27:
   n_shape[0]=((NDOF_K>1)? 2  : 1); n_shape[1]=((NDOF_K>1)?3 : 1); n_shape[2]=((NDOF_K>1)? 4  : 1);   // hex piecewise shapes
   deg[0]=((NDOF_K>0)? 1 : 0);     deg[1]=((NDOF_K>0)? 1 : 0);     deg[2]=((NDOF_K>0)? 1 : 0);        // degree of the shape
   n_gauss[0]=3; n_gauss[1]=9; n_gauss[2]=27;   // 3x3x3 gaussian points
  break;
  case 10:
   n_shape[0]=1; n_shape[1]=1; n_shape[2]=1;   // tet piecewise shapes
   deg[0]=0;     deg[1]=0;     deg[2]=0;        // degree of the shape
   n_gauss[0]=3; n_gauss[1]=7; n_gauss[2]=14;   // 3x3x3 gaussian points
  break;
  default :
    std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
  break;
}
else std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";

// Fem -------------------------
  for (int idim=0;idim<3;idim++) {
    _NoGauss1[idim]=n_gauss[idim];// # of Gaussian points
    _NoShape[idim]=n_shape[idim]; //  # of shape functions 
    _deg[idim]=deg[idim];
  }

// weights, shapes and shape derivatives -----------------------------
  for (int idim=0;idim<3;idim++) {
    _weight1[idim]= new double [_NoGauss1[idim]];
    _phi_map1[idim]=new double[_NoShape[idim]*_NoGauss1[idim]];
    _dphidxez_map1[idim]=new double[(idim+1)*_NoShape[idim]*_NoGauss1[idim]];
    _dphidxez_map1_nodes[idim]=new double[(idim+1)*_NoShape[idim]*_NoGauss1[idim]];
    _dphidxx_map1[idim]=new double[(idim+1)*(idim+1)*_NoShape[idim]*_NoGauss1[idim]];
  }

#ifdef PRINT_INFO
  std::cout << "\n MGFE::MGFE: gaussian points: "
            << _NoGauss1[DIMENSION-1] << " \n";
  std::cout << " MGFE::MGFE: polynomial order " << _order
            << " - # shape func=(";
  for(int idim=0;idim<3;idim++) std::cout <<	_NoShape[idim] << ",";    
      std::cout    << ") \n";
#endif
  return;
}

// ==============================
/// This function destroys the MGFE class
MGFE::~MGFE(
) {// ==============================
  clear();
}

// ================================
/// This function clears the substructures
void MGFE::clear(
) {// ========================
  for (int idim=0;idim<_dim;idim++) {
    delete []_weight1[idim];       // weight
    delete []_phi_map1[idim];      // shapes at g.p.
    delete []_dphidxez_map1[idim]; // derivatives at g.p.
    delete []_dphidxez_map1_nodes[idim]; // derivatives on nodes
    delete []_dphidxx_map1[idim]; // derivatives on nodes
  }
  return;
}

// // ==========================================
// /// This function initializes the MGFE class from file:
// void MGFE::init(
// ) {// ==========================================
// // #ifdef PRINT_INFO
// //   std::cout << "\n Reading FE:: =============== " << std::endl;
// // #endif
// // 
// //   std::string    f_shape = _mgutils.get_file("F_SHAPE");
// //   std::ostringstream file;
// // 
// //   switch (_order) {
// //     /// A) _order =2 family Hex27(3D) Quad9(2D) Lin 3(1D)
// // //   if (_order == 2) {
// //   case 2:
// //     // a) for  3D functions ("fem/shape3D_2727.in")
// //     // for  2D functions ("fem/shape2D_0909.in")
// //     // for  1D functions ("fem/shape1D_0303.in")
// //     for (int idim=1;idim<= DIMENSION;idim++) {
// //       file.str("");
// //       file << _mgutils._fem_dir  << f_shape << idim<< "D_"
// //       << _NoGauss1[idim-1] << "." << _NoShape[idim-1] << ".in";
// //       std::ifstream in(file.str().c_str()); read_c(in,idim);
// //     }
// //     break;
// // //   }
// //     /// B) mode =8 family Hex8(3D) Quad4(2D) Lin 2(1D)
// // //   else if (_order == 1) {
// //   case 1:
// //     // a) for  3D functions  ("fem/shape3D_2708.in")
// //     // b) for  2D functions ("fem/shape2D_0904.in")
// //     // c) for  1D functions ("fem/shape1D_0302.in")
// //     for (int idim=1;idim<= DIMENSION;idim++) {
// //       file.str("");
// //       file << _mgutils._fem_dir
// //       << f_shape << idim<< "D_"  << _NoGauss1[idim-1] << "." << _NoShape[idim-1]<< ".in";
// //       std::ifstream in(file.str().c_str()); read_c(in,idim);
// //     }
// //     break;
// // 
// //     /// c) const
// // //   else if (_order == 1) {
// //   case 0:
// //     // a) for  3D functions  ("fem/shape3D_2708.in")
// //     // b) for  2D functions ("fem/shape2D_0904.in")
// //     // c) for  1D functions ("fem/shape1D_0302.in")
// //     for (int idim=1;idim<= DIMENSION;idim++) {
// //       file.str("");
// //       file << _mgutils._fem_dir
// //       << f_shape << idim<< "D_"  << _NoGauss1[idim-1] << "." <<_order<< "." << _NoShape[idim-1]<< ".in";
// //       std::cout << file.str() << " file \n";
// //       std::ifstream in(file.str().c_str());read_c(in,idim);
// //    
// //     }
// //     break;
// // 
// //   default:
// //     std::cerr << "MGFE::init: FE order " << _order << " not supported" << std::endl; abort();
// //   break;
// //  }
//   return;
// }


/// /// This function generates the Lagrangian quad shape functions
void MGFE::init_qua(
) {// ================================

#if ELTYPE==27

//               ********************************************
//                               QUAD 9  
// 		 ********************************************
//                              
// 			       3 ______6_____ 2
// 				|            |
// 			        |            |
//			       7|      8     |5
//				|            |
// 			        |____________|
//			       0       4      1

//gaussian coordinates
 const double a=-sqrt(3./5.); const double b=0.; const double c=-a;
 const double a_n=-1.; const double b_n=0.; const double c_n=1.;
 
 const double x1D[3]={a,b,c};
 const double x1D_n[3]={a_n,b_n,c_n};
 
 const double x2D[9]={a,a,a,b,b,b,c,c,c};
 const double y2D[9]={a,b,c,a,b,c,a,b,c};
 const double x2D_n[9]={a_n,c_n,c_n,a_n,b_n,c_n,b_n,a_n,b_n};
 const double y2D_n[9]={a_n,a_n,c_n,c_n,a_n,b_n,c_n,b_n,b_n};
 
 const double x3D[27]={a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
 const double y3D[27]={a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
 const double z3D[27]={a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};
 const double x3D_n[27]={a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n};
 const double y3D_n[27]={a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n};
 const double z3D_n[27]={a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n};
 
 
//gaussian weights
// 1D --------------------------------
 _weight1[0][0]= 5./9.; _weight1[0][1]= 8./9.; _weight1[0][2]= 5./9.;
 
// 2D --------------------------------
 const double m=_weight1[0][0]*_weight1[0][1]; const double l=_weight1[0][0]*_weight1[0][0];
 const double h=_weight1[0][1]*_weight1[0][1];
 const double weight[9]={l,m,l,m,h,m,l,m,l};
 for (int i=0;i<_NoGauss1[1];i++) _weight1[1][i]= weight[i];
 
// 3D --------------------------------
 const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0]; const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
 const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1]; const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
 const double weight_3[27]={w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
 for (int i=0;i<_NoGauss1[2];i++) _weight1[2][i]= weight_3[i];

 
// ****QUADRATIC SHAPES AND DERIVATIVES****

  //1D------------------------------------ 
  for (int i=0;i<_NoGauss1[0];i++) {
    const double xx=x1D[i]; const double xx_n=x1D_n[i];
    //shape functions
   _phi_map1[0][i]               =.5*xx*(xx - 1.);
   _phi_map1[0][i+_NoGauss1[0]]  =.5*xx*(xx + 1);
   _phi_map1[0][i+2*_NoGauss1[0]]=(1. - xx*xx);
   //d/dx gaussian points
   _dphidxez_map1[0][i]                = xx-.5;
   _dphidxez_map1[0][i+_NoGauss1[0]]   =  xx+.5;
   _dphidxez_map1[0][i+2*_NoGauss1[0]] = -2.*xx;
   //     nodal points
   _dphidxez_map1_nodes[0][i]                 = xx_n-.5;
   _dphidxez_map1_nodes[0][i+_NoGauss1[0]]    =  xx_n+.5;
   _dphidxez_map1_nodes[0][i+2*_NoGauss1[0]]  = -2.*xx_n;
   //d2/dx2 gaussian points
   _dphidxx_map1[0][i]                = 1.;
   _dphidxx_map1[0][i+_NoGauss1[0]]   = 1.;
   _dphidxx_map1[0][i+2*_NoGauss1[0]] = -2.;
  }
    
  // 2D --------------------------------
  for (int i=0;i<_NoGauss1[1];i++) { 
    const double xx=x2D[i]; const double yy=y2D[i];
    const double xx_n=x2D_n[i]; const double yy_n=y2D_n[i];
    //shape functions
    _phi_map1[1][i+0*_NoGauss1[1]]= .5*xx*(xx - 1.)*.5*yy*(yy - 1.);
    _phi_map1[1][i+1*_NoGauss1[1]]= .5*xx*(xx + 1.)*.5*yy*(yy - 1.);
    _phi_map1[1][i+2*_NoGauss1[1]]= .5*xx*(xx + 1)*.5*yy*(yy + 1);
    _phi_map1[1][i+3*_NoGauss1[1]]= .5*xx*(xx - 1.)*.5*yy*(yy + 1);
    _phi_map1[1][i+4*_NoGauss1[1]]= (1. - xx*xx)*.5*yy*(yy - 1.);
    _phi_map1[1][i+5*_NoGauss1[1]]= .5*xx*(xx + 1)*(1. - yy*yy);
    _phi_map1[1][i+6*_NoGauss1[1]]= (1. - xx*xx)*.5*yy*(yy + 1);
    _phi_map1[1][i+7*_NoGauss1[1]]= .5*xx*(xx - 1.)* (1. - yy*yy);
    _phi_map1[1][i+8*_NoGauss1[1]]= (1. - xx*xx)*(1. - yy*yy);   
    
    //d/dx gaussian points
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]=  (xx - 0.5)*.5*yy*(yy - 1.);
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]=  (xx + 0.5)*.5*yy*(yy - 1.);
    _dphidxez_map1[1][i+(2)*_NoGauss1[1]]=  (xx + 0.5)*.5*yy*(yy + 1);
    _dphidxez_map1[1][i+(3)*_NoGauss1[1]]=  (xx - 0.5)*.5*yy*(yy + 1);
    _dphidxez_map1[1][i+(4)*_NoGauss1[1]]=  (-2.*xx)  *.5*yy*(yy - 1.);
    _dphidxez_map1[1][i+(5)*_NoGauss1[1]]=  (xx + 0.5)*(1. - yy*yy);
    _dphidxez_map1[1][i+(6)*_NoGauss1[1]]=  (-2.*xx)* .5*yy*(yy + 1);
    _dphidxez_map1[1][i+(7)*_NoGauss1[1]]=  (xx - 0.5)*(1. - yy*yy);
    _dphidxez_map1[1][i+(8)*_NoGauss1[1]]=  (-2.*xx)*(1. - yy*yy); 
    //     nodal points
    _dphidxez_map1_nodes[1][i+(0)*_NoGauss1[1]]=  (xx_n - 0.5)*.5*yy_n*(yy_n - 1.);
    _dphidxez_map1_nodes[1][i+(1)*_NoGauss1[1]]=  (xx_n + 0.5)*.5*yy_n*(yy_n - 1.);
    _dphidxez_map1_nodes[1][i+(2)*_NoGauss1[1]]=  (xx_n + 0.5)*.5*yy_n*(yy_n +  1);
    _dphidxez_map1_nodes[1][i+(3)*_NoGauss1[1]]=  (xx_n - 0.5) *  .5*yy_n*(yy_n+ 1);
    _dphidxez_map1_nodes[1][i+(4)*_NoGauss1[1]]=  (-2.*xx_n) *  .5*yy_n*(yy_n - 1.);
    _dphidxez_map1_nodes[1][i+(5)*_NoGauss1[1]]=  (xx_n + 0.5)*(1. - yy_n*yy_n);
    _dphidxez_map1_nodes[1][i+(6)*_NoGauss1[1]]=  (-2.*xx_n)*.5*yy_n*(yy_n + 1);
    _dphidxez_map1_nodes[1][i+(7)*_NoGauss1[1]]=  (xx_n - 0.5)* (1. - yy_n*yy_n);
    _dphidxez_map1_nodes[1][i+(8)*_NoGauss1[1]]=  (-2.*xx_n)*(1. - yy_n*yy_n);
    
    //d/dy gaussian points
    _dphidxez_map1[1][i+(9)*_NoGauss1[1]]=  .5*xx*(xx - 1.)*(yy - 0.5);
    _dphidxez_map1[1][i+(10)*_NoGauss1[1]]= .5*xx*(xx + 1)*(yy - 0.5);
    _dphidxez_map1[1][i+(11)*_NoGauss1[1]]= .5*xx*(xx + 1)*(yy + 0.5);
    _dphidxez_map1[1][i+(12)*_NoGauss1[1]]= .5*xx*(xx - 1.) * (yy + 0.5);
    _dphidxez_map1[1][i+(13)*_NoGauss1[1]]= (1. - xx*xx)*(yy - 0.5);
    _dphidxez_map1[1][i+(14)*_NoGauss1[1]]= .5*xx*(xx + 1)*(-2.*yy);
    _dphidxez_map1[1][i+(15)*_NoGauss1[1]]= (1. - xx*xx)*(yy + 0.5);
    _dphidxez_map1[1][i+(16)*_NoGauss1[1]]= .5*xx*(xx - 1.)* (-2.*yy);
    _dphidxez_map1[1][i+(17)*_NoGauss1[1]]= (1. - xx*xx)*(-2.*yy);   
    //     nodal points
    _dphidxez_map1_nodes[1][i+(9)*_NoGauss1[1]]=  .5*xx_n*(xx_n - 1.)*(yy_n - 0.5);
    _dphidxez_map1_nodes[1][i+(10)*_NoGauss1[1]]= .5*xx_n*(xx_n + 1.)*(yy_n - 0.5);
    _dphidxez_map1_nodes[1][i+(11)*_NoGauss1[1]]= .5*xx_n*(xx_n + 1.)*(yy_n + 0.5);
    _dphidxez_map1_nodes[1][i+(12)*_NoGauss1[1]]= .5*xx_n*(xx_n - 1.)*(yy_n + 0.5);
    _dphidxez_map1_nodes[1][i+(13)*_NoGauss1[1]]= (1.-xx_n*xx_n)*(yy_n - 0.5);
    _dphidxez_map1_nodes[1][i+(14)*_NoGauss1[1]]= .5*xx_n*(xx_n + 1.)*(-2.*yy_n);
    _dphidxez_map1_nodes[1][i+(15)*_NoGauss1[1]]= (1.-xx_n*xx_n)*(yy_n + 0.5);
    _dphidxez_map1_nodes[1][i+(16)*_NoGauss1[1]]= .5*xx_n*(xx_n - 1.)*(-2.*yy_n);
    _dphidxez_map1_nodes[1][i+(17)*_NoGauss1[1]]= (1.-xx_n*xx_n)*(-2.*yy_n); 
    
    
    //d2/dx2 gaussian points
    _dphidxx_map1[1][i+(0)*_NoGauss1[1]]=  .5*yy*(yy - 1.);
    _dphidxx_map1[1][i+(1)*_NoGauss1[1]]=  .5*yy*(yy - 1.);
    _dphidxx_map1[1][i+(2)*_NoGauss1[1]]=  .5*yy*(yy + 1);
    _dphidxx_map1[1][i+(3)*_NoGauss1[1]]=  .5*yy*(yy + 1);
    _dphidxx_map1[1][i+(4)*_NoGauss1[1]]=  -1.*yy*(yy - 1.);
    _dphidxx_map1[1][i+(5)*_NoGauss1[1]]=  (1. - yy*yy);
    _dphidxx_map1[1][i+(6)*_NoGauss1[1]]=  -1.*yy*(yy + 1);
    _dphidxx_map1[1][i+(7)*_NoGauss1[1]]=  (1. - yy*yy);
    _dphidxx_map1[1][i+(8)*_NoGauss1[1]]=  -2.*(1. - yy*yy);
    
     //d2/dxdy gaussian points
    _dphidxx_map1[1][i+(9)*_NoGauss1[1]]=   (xx - 0.5)*(yy - 0.5);
    _dphidxx_map1[1][i+(10)*_NoGauss1[1]]=  (xx + 0.5)*(yy - 0.5);
    _dphidxx_map1[1][i+(11)*_NoGauss1[1]]=  (xx + 0.5)*(yy + 0.5);
    _dphidxx_map1[1][i+(12)*_NoGauss1[1]]=  (xx - 0.5)*(yy + 0.5);
    _dphidxx_map1[1][i+(13)*_NoGauss1[1]]=  -2.*xx*(yy - 0.5);
    _dphidxx_map1[1][i+(14)*_NoGauss1[1]]=  (xx + 0.5)*(-2.*yy);
    _dphidxx_map1[1][i+(15)*_NoGauss1[1]]=  -2.*xx*(yy + 0.5);
    _dphidxx_map1[1][i+(16)*_NoGauss1[1]]=  (xx - 0.5)*(-2.*yy);
    _dphidxx_map1[1][i+(17)*_NoGauss1[1]]=  4.*xx*yy; 
     //d2/dxdy gaussian points
    _dphidxx_map1[1][i+(18)*_NoGauss1[1]]=  (xx - 0.5)*(yy - 0.5);
    _dphidxx_map1[1][i+(19)*_NoGauss1[1]]=  (xx + 0.5)*(yy - 0.5);
    _dphidxx_map1[1][i+(20)*_NoGauss1[1]]=  (xx + 0.5)*(yy + 0.5);
    _dphidxx_map1[1][i+(21)*_NoGauss1[1]]=  (xx - 0.5)*(yy + 0.5);
    _dphidxx_map1[1][i+(22)*_NoGauss1[1]]=  -2.*xx*(yy - 0.5);
    _dphidxx_map1[1][i+(23)*_NoGauss1[1]]=  (xx + 0.5)*(-2.*yy);
    _dphidxx_map1[1][i+(24)*_NoGauss1[1]]=  -2.*xx*(yy + 0.5);
    _dphidxx_map1[1][i+(25)*_NoGauss1[1]]=  (xx - 0.5)*(-2.*yy);
    _dphidxx_map1[1][i+(26)*_NoGauss1[1]]=  4.*xx*yy; 

    //d2/dy2 gaussian points
    _dphidxx_map1[1][i+(27)*_NoGauss1[1]]= xx*(0.5*xx - 0.5);
    _dphidxx_map1[1][i+(28)*_NoGauss1[1]]= xx*(0.5*xx + 0.5);
    _dphidxx_map1[1][i+(29)*_NoGauss1[1]]= xx*(0.5*xx + 0.5);
    _dphidxx_map1[1][i+(30)*_NoGauss1[1]]= xx*(0.5*xx - 0.5);
    _dphidxx_map1[1][i+(31)*_NoGauss1[1]]= (1. - xx*xx);
    _dphidxx_map1[1][i+(32)*_NoGauss1[1]]= -1.*xx*(xx + 1);
    _dphidxx_map1[1][i+(33)*_NoGauss1[1]]= (1. - xx*xx);
    _dphidxx_map1[1][i+(34)*_NoGauss1[1]]= -1.*xx*(xx - 1.);
    _dphidxx_map1[1][i+(35)*_NoGauss1[1]]= -2.*(1. - xx*xx);
  }

  // 3D -----------------------------------------------
  if (_dim == 3) {
    for (int i=0;i<_NoGauss1[2];i++) {
      const double xx=x3D[i]; const double yy=y3D[i]; const double zz=z3D[i];
      const double xx_n=x3D_n[i]; const double yy_n=y3D_n[i]; const double zz_n=z3D_n[i];
      //shape functions
      _phi_map1[2][i+0*_NoGauss1[2]]= .5*xx*(xx - 1.)*.5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _phi_map1[2][i+1*_NoGauss1[2]]=  .5*xx*(xx + 1)*.5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _phi_map1[2][i+2*_NoGauss1[2]]= .5*xx*(xx + 1) *.5*yy*(yy + 1)  * .5*zz*(zz - 1.);
      _phi_map1[2][i+3*_NoGauss1[2]]= .5*xx*(xx - 1.)*.5*yy*(yy + 1)  * .5*zz*(zz - 1.);
      _phi_map1[2][i+4*_NoGauss1[2]]= .5*xx*(xx - 1.)*.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _phi_map1[2][i+5*_NoGauss1[2]]= .5*xx*(xx + 1) *.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _phi_map1[2][i+6*_NoGauss1[2]]=  .5*xx*(xx + 1)* .5*yy*(yy + 1) * .5*zz*(zz + 1) ;
      _phi_map1[2][i+7*_NoGauss1[2]]= .5*xx*(xx - 1.)* .5*yy*(yy + 1) * .5*zz*(zz + 1) ;
      _phi_map1[2][i+8*_NoGauss1[2]]= (1. - xx*xx)   * .5*yy*(yy - 1.)* .5*zz*(zz - 1.);
      _phi_map1[2][i+9*_NoGauss1[2]]= .5*xx*(xx + 1) * (1. - yy*yy)   * .5*zz*(zz - 1.);
      _phi_map1[2][i+10*_NoGauss1[2]]= (1. - xx*xx)  * .5*yy*(yy + 1) * .5*zz*(zz - 1.);
      _phi_map1[2][i+11*_NoGauss1[2]]=.5*xx*(xx - 1.)*   (1. - yy*yy) * .5*zz*(zz - 1.);
      _phi_map1[2][i+12*_NoGauss1[2]]=.5*xx*(xx - 1.)* .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _phi_map1[2][i+13*_NoGauss1[2]]= .5*xx*(xx + 1)* .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _phi_map1[2][i+14*_NoGauss1[2]]= .5*xx*(xx + 1)* .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _phi_map1[2][i+15*_NoGauss1[2]]=.5*xx*(xx - 1.)* .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _phi_map1[2][i+16*_NoGauss1[2]]=  (1. - xx*xx) *.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _phi_map1[2][i+17*_NoGauss1[2]]= .5*xx*(xx + 1)* (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _phi_map1[2][i+18*_NoGauss1[2]]=  (1. - xx*xx) *.5*yy*(yy + 1)  * .5*zz*(zz + 1) ;
      _phi_map1[2][i+19*_NoGauss1[2]]=.5*xx*(xx - 1.)* (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _phi_map1[2][i+20*_NoGauss1[2]]= (1. - xx*xx)  *  (1. - yy*yy)  * .5*zz*(zz - 1.);
      _phi_map1[2][i+21*_NoGauss1[2]]=  (1. - xx*xx) *.5*yy*(yy - 1.) *  (1. - zz*zz)  ;
      _phi_map1[2][i+22*_NoGauss1[2]]= .5*xx*(xx + 1)*  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _phi_map1[2][i+23*_NoGauss1[2]]=  (1. - xx*xx) * .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _phi_map1[2][i+24*_NoGauss1[2]]=.5*xx*(xx - 1.)*  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _phi_map1[2][i+25*_NoGauss1[2]]=  (1. - xx*xx) *  (1. - yy*yy) * .5*zz*(zz + 1)  ;
      _phi_map1[2][i+26*_NoGauss1[2]]=  (1. - xx*xx) *  (1. - yy*yy)  * (1. - zz*zz)   ;
      
       // df/dx gaussian points
      _dphidxez_map1[2][i+(0)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(1)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(2)*_NoGauss1[2]]=  (xx + 0.5) *.5*yy*(yy + 1)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(3)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy + 1)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(4)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(5)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(6)*_NoGauss1[2]]=  (xx + 0.5)* .5*yy*(yy + 1) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(7)*_NoGauss1[2]]=  (xx - 0.5)* .5*yy*(yy + 1) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(8)*_NoGauss1[2]]=  (- 2.*xx)   * .5*yy*(yy - 1.)* .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(9)*_NoGauss1[2]]=  (xx + 0.5) * (1. - yy*yy)   * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(10)*_NoGauss1[2]]=  (- 2.*xx)   * .5*yy*(yy + 1) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(11)*_NoGauss1[2]]= (xx - 0.5)*   (1. - yy*yy) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(12)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(13)*_NoGauss1[2]]=  (xx + 0.5)* .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(14)*_NoGauss1[2]]=  (xx + 0.5)* .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(15)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(16)*_NoGauss1[2]]=   (- 2.*xx)  *.5*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(17)*_NoGauss1[2]]=  (xx + 0.5)* (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(18)*_NoGauss1[2]]=   (- 2.*xx)  *.5*yy*(yy + 1)  * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(19)*_NoGauss1[2]]= (xx - 0.5)* (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(20)*_NoGauss1[2]]=  (- 2.*xx)   *  (1. - yy*yy)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(21)*_NoGauss1[2]]=  (- 2.*xx)  *.5*yy*(yy - 1.) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(22)*_NoGauss1[2]]=  (xx + 0.5)*  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(23)*_NoGauss1[2]]=   (- 2.*xx)  * .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(24)*_NoGauss1[2]]= (xx - 0.5)*  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(25)*_NoGauss1[2]]=   (- 2.*xx)  *  (1. - yy*yy) * .5*zz*(zz + 1)  ;
      _dphidxez_map1[2][i+(26)*_NoGauss1[2]]=  (- 2.*xx)  *  (1. - yy*yy)  * (1. - zz*zz)   ;
      //    nodal points
      _dphidxez_map1_nodes[2][i+(0)*_NoGauss1[2]]=  (xx_n - 0.5)*.5*yy_n*(yy_n - 1.) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(1)*_NoGauss1[2]]=  (xx_n + 0.5)*.5*yy_n*(yy_n - 1.) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(2)*_NoGauss1[2]]=  (xx_n + 0.5) *.5*yy_n*(yy_n + 1)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(3)*_NoGauss1[2]]=  (xx_n - 0.5)*.5*yy_n*(yy_n + 1)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(4)*_NoGauss1[2]]=  (xx_n - 0.5)*.5*yy_n*(yy_n - 1.) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(5)*_NoGauss1[2]]=  (xx_n + 0.5)*.5*yy_n*(yy_n - 1.) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(6)*_NoGauss1[2]]=  (xx_n + 0.5)* .5*yy_n*(yy_n + 1) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(7)*_NoGauss1[2]]=  (xx_n - 0.5)* .5*yy_n*(yy_n + 1) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(8)*_NoGauss1[2]]=  (- 2.*xx_n)   * .5*yy_n*(yy_n - 1.)* .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(9)*_NoGauss1[2]]=  (xx_n + 0.5) * (1. - yy_n*yy_n)   * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(10)*_NoGauss1[2]]=  (- 2.*xx_n)   * .5*yy_n*(yy_n + 1) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(11)*_NoGauss1[2]]= (xx_n - 0.5)*   (1. - yy_n*yy_n) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(12)*_NoGauss1[2]]= (xx_n - 0.5)* .5*yy_n*(yy_n - 1.)*  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(13)*_NoGauss1[2]]=  (xx_n + 0.5)* .5*yy_n*(yy_n - 1.)*  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(14)*_NoGauss1[2]]=  (xx_n + 0.5)* .5*yy_n*(yy_n + 1) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(15)*_NoGauss1[2]]= (xx_n - 0.5)* .5*yy_n*(yy_n + 1) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(16)*_NoGauss1[2]]=   (- 2.*xx_n)  *.5*yy_n*(yy_n - 1.) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(17)*_NoGauss1[2]]=  (xx_n + 0.5)* (1. - yy_n*yy_n)   * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(18)*_NoGauss1[2]]=   (- 2.*xx_n)  *.5*yy_n*(yy_n + 1)  * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(19)*_NoGauss1[2]]= (xx_n - 0.5)* (1. - yy_n*yy_n)   * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(20)*_NoGauss1[2]]=  (- 2.*xx_n)   *  (1. - yy_n*yy_n)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(21)*_NoGauss1[2]]=  (- 2.*xx_n)  *.5*yy_n*(yy_n - 1.) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(22)*_NoGauss1[2]]=  (xx_n + 0.5)*  (1. - yy_n*yy_n)  *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(23)*_NoGauss1[2]]=   (- 2.*xx_n)  * .5*yy_n*(yy_n + 1) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(24)*_NoGauss1[2]]= (xx_n - 0.5)*  (1. - yy_n*yy_n)  *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(25)*_NoGauss1[2]]=   (- 2.*xx_n)  *  (1. - yy_n*yy_n) * .5*zz_n*(zz_n + 1)  ;
      _dphidxez_map1_nodes[2][i+(26)*_NoGauss1[2]]=  (- 2.*xx_n)  *  (1. - yy_n*yy_n)  * (1. - zz_n*zz_n)   ;
      
      
      // df/dy    gaussian points                             
      _dphidxez_map1[2][i+(27)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(28)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(29)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy + 0.5)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(30)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(31)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(32)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy - 0.5) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(33)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(34)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(35)*_NoGauss1[2]]=   (1. - xx*xx)   * (yy - 0.5)* .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(36)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (-2.*yy)   * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(37)*_NoGauss1[2]]=    (1. - xx*xx)  * (yy + 0.5) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(38)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*   (-2.*yy) * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(39)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5)*  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(40)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5)*  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(41)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(42)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(43)*_NoGauss1[2]]=     (1. - xx*xx) *(yy - 0.5) * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(44)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (-2.*yy)   * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(45)*_NoGauss1[2]]=     (1. - xx*xx) *(yy + 0.5)  * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(46)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (-2.*yy)   * .5*zz*(zz + 1) ;
      _dphidxez_map1[2][i+(47)*_NoGauss1[2]]=    (1. - xx*xx)  *  (-2.*yy)  * .5*zz*(zz - 1.);
      _dphidxez_map1[2][i+(48)*_NoGauss1[2]]=     (1. - xx*xx) * (yy - 0.5) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(49)*_NoGauss1[2]]=    .5*xx*(xx + 1)*  (-2.*yy)  *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(50)*_NoGauss1[2]]=     (1. - xx*xx) * (yy + 0.5) *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(51)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*  (-2.*yy)  *  (1. - zz*zz)  ;
      _dphidxez_map1[2][i+(52)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy) * .5*zz*(zz + 1)  ;
      _dphidxez_map1[2][i+(53)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy)  * (1. - zz*zz)   ;
      //    nodal points
      _dphidxez_map1_nodes[2][i+(27)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n - 0.5) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(28)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (yy_n - 0.5) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(29)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) * (yy_n + 0.5)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(30)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n + 0.5)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(31)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n - 0.5) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(32)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) * (yy_n - 0.5) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(33)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (yy_n + 0.5) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(34)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n + 0.5) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(35)*_NoGauss1[2]]=   (1. - xx_n*xx_n)   * (yy_n - 0.5)* .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(36)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) * (-2.*yy_n)   * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(37)*_NoGauss1[2]]=    (1. - xx_n*xx_n)  * (yy_n + 0.5) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(38)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*   (-2.*yy_n) * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(39)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n - 0.5)*  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(40)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (yy_n - 0.5)*  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(41)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (yy_n + 0.5) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(42)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (yy_n + 0.5) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(43)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *(yy_n - 0.5) * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(44)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (-2.*yy_n)   * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(45)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *(yy_n + 0.5)  * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(46)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (-2.*yy_n)   * .5*zz_n*(zz_n + 1) ;
      _dphidxez_map1_nodes[2][i+(47)*_NoGauss1[2]]=    (1. - xx_n*xx_n)  *  (-2.*yy_n)  * .5*zz_n*(zz_n - 1.);
      _dphidxez_map1_nodes[2][i+(48)*_NoGauss1[2]]=     (1. - xx_n*xx_n) * (yy_n - 0.5) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(49)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)*  (-2.*yy_n)  *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(50)*_NoGauss1[2]]=     (1. - xx_n*xx_n) * (yy_n + 0.5) *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(51)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*  (-2.*yy_n)  *  (1. - zz_n*zz_n)  ;
      _dphidxez_map1_nodes[2][i+(52)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *  (-2.*yy_n) * .5*zz_n*(zz_n + 1)  ;
      _dphidxez_map1_nodes[2][i+(53)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *  (-2.*yy_n)  * (1. - zz_n*zz_n)   ;
      
      // df/dz gaussian points               
      _dphidxez_map1[2][i+(54)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy - 1.) * (zz - 0.5);
      _dphidxez_map1[2][i+(55)*_NoGauss1[2]]=    .5*xx*(xx + 1)*.5*yy*(yy - 1.) * (zz - 0.5);
      _dphidxez_map1[2][i+(56)*_NoGauss1[2]]=   .5*xx*(xx + 1) *.5*yy*(yy + 1)  * (zz - 0.5);
      _dphidxez_map1[2][i+(57)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy + 1)  * (zz - 0.5);
      _dphidxez_map1[2][i+(58)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy - 1.) * (zz + 0.5) ;
      _dphidxez_map1[2][i+(59)*_NoGauss1[2]]=   .5*xx*(xx + 1) *.5*yy*(yy - 1.) * (zz + 0.5) ;
      _dphidxez_map1[2][i+(60)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy + 1) * (zz + 0.5) ;
      _dphidxez_map1[2][i+(61)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy + 1) * (zz + 0.5) ;
      _dphidxez_map1[2][i+(62)*_NoGauss1[2]]=   (1. - xx*xx)   * .5*yy*(yy - 1.)* (zz - 0.5);
      _dphidxez_map1[2][i+(63)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (1. - yy*yy)   * (zz - 0.5);
      _dphidxez_map1[2][i+(64)*_NoGauss1[2]]=    (1. - xx*xx)  * .5*yy*(yy + 1) * (zz - 0.5);
      _dphidxez_map1[2][i+(65)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*   (1. - yy*yy) * (zz - 0.5);
      _dphidxez_map1[2][i+(66)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy - 1.)* (-2.* zz)  ;
      _dphidxez_map1[2][i+(67)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy - 1.)* (-2.* zz)  ;
      _dphidxez_map1[2][i+(68)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy + 1) * (-2.* zz)  ;
      _dphidxez_map1[2][i+(69)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy + 1) * (-2.* zz)  ;
      _dphidxez_map1[2][i+(70)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy - 1.) * (zz + 0.5) ;
      _dphidxez_map1[2][i+(71)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (1. - yy*yy)   * (zz + 0.5) ;
      _dphidxez_map1[2][i+(72)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy + 1)  * (zz + 0.5) ;
      _dphidxez_map1[2][i+(73)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (1. - yy*yy)   * (zz + 0.5) ;
      _dphidxez_map1[2][i+(74)*_NoGauss1[2]]=    (1. - xx*xx)  *  (1. - yy*yy)  * (zz - 0.5);
      _dphidxez_map1[2][i+(75)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy - 1.) * (-2.* zz)   ;
      _dphidxez_map1[2][i+(76)*_NoGauss1[2]]=    .5*xx*(xx + 1)*  (1. - yy*yy)  * (-2.* zz)   ;
      _dphidxez_map1[2][i+(77)*_NoGauss1[2]]=     (1. - xx*xx) * .5*yy*(yy + 1) * (-2.* zz)   ;
      _dphidxez_map1[2][i+(78)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*  (1. - yy*yy)  * (-2.* zz)   ;
      _dphidxez_map1[2][i+(79)*_NoGauss1[2]]=     (1. - xx*xx) *  (1. - yy*yy) * (zz + 0.5)  ;
      _dphidxez_map1[2][i+(80)*_NoGauss1[2]]=     (1. - xx*xx) *  (1. - yy*yy)  * (-2.* zz)   ;
      //    nodal points
      _dphidxez_map1_nodes[2][i+(54)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*.5*yy_n*(yy_n - 1.) * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(55)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)*.5*yy_n*(yy_n - 1.) * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(56)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) *.5*yy_n*(yy_n + 1)  * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(57)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*.5*yy_n*(yy_n + 1)  * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(58)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*.5*yy_n*(yy_n - 1.) * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(59)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) *.5*yy_n*(yy_n - 1.) * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(60)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* .5*yy_n*(yy_n + 1) * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(61)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* .5*yy_n*(yy_n + 1) * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(62)*_NoGauss1[2]]=   (1. - xx_n*xx_n)   * .5*yy_n*(yy_n - 1.)* (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(63)*_NoGauss1[2]]=   .5*xx_n*(xx_n + 1) * (1. - yy_n*yy_n)   * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(64)*_NoGauss1[2]]=    (1. - xx_n*xx_n)  * .5*yy_n*(yy_n + 1) * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(65)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*   (1. - yy_n*yy_n) * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(66)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* .5*yy_n*(yy_n - 1.)* (-2.* zz_n)  ;
      _dphidxez_map1_nodes[2][i+(67)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* .5*yy_n*(yy_n - 1.)* (-2.* zz_n)  ;
      _dphidxez_map1_nodes[2][i+(68)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* .5*yy_n*(yy_n + 1) * (-2.* zz_n)  ;
      _dphidxez_map1_nodes[2][i+(69)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* .5*yy_n*(yy_n + 1) * (-2.* zz_n)  ;
      _dphidxez_map1_nodes[2][i+(70)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *.5*yy_n*(yy_n - 1.) * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(71)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)* (1. - yy_n*yy_n)   * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(72)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *.5*yy_n*(yy_n + 1)  * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(73)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)* (1. - yy_n*yy_n)   * (zz_n + 0.5) ;
      _dphidxez_map1_nodes[2][i+(74)*_NoGauss1[2]]=    (1. - xx_n*xx_n)  *  (1. - yy_n*yy_n)  * (zz_n - 0.5);
      _dphidxez_map1_nodes[2][i+(75)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *.5*yy_n*(yy_n - 1.) * (-2.* zz_n)   ;
      _dphidxez_map1_nodes[2][i+(76)*_NoGauss1[2]]=    .5*xx_n*(xx_n + 1)*  (1. - yy_n*yy_n)  * (-2.* zz_n)   ;
      _dphidxez_map1_nodes[2][i+(77)*_NoGauss1[2]]=     (1. - xx_n*xx_n) * .5*yy_n*(yy_n + 1) * (-2.* zz_n)   ;
      _dphidxez_map1_nodes[2][i+(78)*_NoGauss1[2]]=   .5*xx_n*(xx_n - 1.)*  (1. - yy_n*yy_n)  * (-2.* zz_n)   ;
      _dphidxez_map1_nodes[2][i+(79)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *  (1. - yy_n*yy_n) * (zz_n + 0.5)  ;
      _dphidxez_map1_nodes[2][i+(80)*_NoGauss1[2]]=     (1. - xx_n*xx_n) *  (1. - yy_n*yy_n)  * (-2.* zz_n)   ;
      
      
      //--------------- Second derivatives ---------------------
       // d2f/dx2 gaussian points
      _dphidxx_map1[2][i+(0)*_NoGauss1[2]]=  .5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(1)*_NoGauss1[2]]=  .5*yy*(yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(2)*_NoGauss1[2]]=  .5*yy*(yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(3)*_NoGauss1[2]]=  .5*yy*(yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(4)*_NoGauss1[2]]=  .5*yy*(yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(5)*_NoGauss1[2]]=  .5*yy*(yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(6)*_NoGauss1[2]]=  .5*yy*(yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(7)*_NoGauss1[2]]=  .5*yy*(yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(8)*_NoGauss1[2]]=  -1.*yy*(yy - 1.)* .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(9)*_NoGauss1[2]]=  (1. - yy*yy)   * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(10)*_NoGauss1[2]]= -1.*yy*(yy + 1) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(11)*_NoGauss1[2]]=  (1. - yy*yy) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(12)*_NoGauss1[2]]= .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(13)*_NoGauss1[2]]= .5*yy*(yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(14)*_NoGauss1[2]]= .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(15)*_NoGauss1[2]]= .5*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(16)*_NoGauss1[2]]= -1.*yy*(yy - 1.) * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(17)*_NoGauss1[2]]=  (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(18)*_NoGauss1[2]]=  -1.*yy*(yy + 1)  * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(19)*_NoGauss1[2]]= (1. - yy*yy)   * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(20)*_NoGauss1[2]]=  -1. *  (1. - yy*yy)  * zz*(zz - 1.);
      _dphidxx_map1[2][i+(21)*_NoGauss1[2]]=  -1.*yy*(yy - 1.) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(22)*_NoGauss1[2]]=  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(23)*_NoGauss1[2]]=   -1.*yy*(yy + 1) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(24)*_NoGauss1[2]]=  (1. - yy*yy)  *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(25)*_NoGauss1[2]]=   -1.* (1. - yy*yy) * zz*(zz + 1)  ;
      _dphidxx_map1[2][i+(26)*_NoGauss1[2]]=  - 2. *  (1. - yy*yy)  * (1. - zz*zz)   ;
      
      // d2f/dxy gaussian points
      _dphidxx_map1[2][i+(27)*_NoGauss1[2]]=  (xx - 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(28)*_NoGauss1[2]]=  (xx + 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(29)*_NoGauss1[2]]=  (xx + 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(30)*_NoGauss1[2]]=  (xx - 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(31)*_NoGauss1[2]]=  (xx - 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(32)*_NoGauss1[2]]=  (xx + 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(33)*_NoGauss1[2]]=  (xx + 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(34)*_NoGauss1[2]]=  (xx - 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(35)*_NoGauss1[2]]=  (- 2.*xx) *.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(36)*_NoGauss1[2]]=  (xx + 0.5)*    (-2.*yy)    * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(37)*_NoGauss1[2]]=  (- 2.*xx) *.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(38)*_NoGauss1[2]]=  (xx - 0.5)*  (-2.*yy)      * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(39)*_NoGauss1[2]]=  (xx - 0.5)* .5*(2.*yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(40)*_NoGauss1[2]]=  (xx + 0.5)* .5*(2.*yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(41)*_NoGauss1[2]]=  (xx + 0.5)* .5*(2.*yy + 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(42)*_NoGauss1[2]]=  (xx - 0.5)* .5*(2.*yy + 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(43)*_NoGauss1[2]]=  (- 2.*xx) * .5*(2.*yy - 1.)* .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(44)*_NoGauss1[2]]=  (xx + 0.5)* (-2.*yy)       * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(45)*_NoGauss1[2]]=  (- 2.*xx) * .5*(2.*yy + 1) * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(46)*_NoGauss1[2]]=  (xx - 0.5)* (-2.*yy)       * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(47)*_NoGauss1[2]]=  (- 2.*xx) * (-2.*yy)       * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(48)*_NoGauss1[2]]=  (- 2.*xx) *.5*(2.*yy - 1.) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(49)*_NoGauss1[2]]=  (xx + 0.5)*  (-2.*yy)      *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(50)*_NoGauss1[2]]=  (- 2.*xx) * .5*(2.*yy + 1) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(51)*_NoGauss1[2]]=  (xx - 0.5)*  (-2.*yy)      *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(52)*_NoGauss1[2]]=  (- 2.*xx) *  (-2.*yy)      * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(53)*_NoGauss1[2]]=  (- 2.*xx) *  (-2.*yy)      * (1. - zz*zz)   ;
      
      // d2f/dxz gaussian points
      _dphidxx_map1[2][i+(54)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(55)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(56)*_NoGauss1[2]]=  (xx + 0.5) *.5*yy*(yy + 1) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(57)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy + 1)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(58)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(59)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(60)*_NoGauss1[2]]=  (xx + 0.5)* .5*yy*(yy + 1) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(61)*_NoGauss1[2]]=  (xx - 0.5)* .5*yy*(yy + 1) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(62)*_NoGauss1[2]]=  (- 2.*xx) * .5*yy*(yy - 1.)* .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(63)*_NoGauss1[2]]=  (xx + 0.5)* (1. - yy*yy)   * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(64)*_NoGauss1[2]]= (- 2.*xx) * .5*yy*(yy + 1) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(65)*_NoGauss1[2]]= (xx - 0.5)*   (1. - yy*yy) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(66)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy - 1.)*  (-2.*zz)  ;
      _dphidxx_map1[2][i+(67)*_NoGauss1[2]]= (xx + 0.5)* .5*yy*(yy - 1.)*  (-2.*zz)  ;
      _dphidxx_map1[2][i+(68)*_NoGauss1[2]]= (xx + 0.5)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(69)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(70)*_NoGauss1[2]]= (- 2.*xx) *.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(71)*_NoGauss1[2]]= (xx + 0.5)* (1. - yy*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(72)*_NoGauss1[2]]= (- 2.*xx) *.5*yy*(yy + 1)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(73)*_NoGauss1[2]]= (xx - 0.5)* (1. - yy*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(74)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(75)*_NoGauss1[2]]=  (- 2.*xx)*.5*yy*(yy - 1.) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(76)*_NoGauss1[2]]= (xx + 0.5)*  (1. - yy*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(77)*_NoGauss1[2]]=  (- 2.*xx)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(78)*_NoGauss1[2]]= (xx - 0.5)*  (1. - yy*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(79)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  * .5*(2.*zz + 1)  ;
      _dphidxx_map1[2][i+(80)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  *  (-2.*zz)   ;
      
      // d2f/dyx gaussian points
      _dphidxx_map1[2][i+(81)*_NoGauss1[2]]=   (xx - 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(82)*_NoGauss1[2]]=   (xx + 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(83)*_NoGauss1[2]]=   (xx + 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(84)*_NoGauss1[2]]=   (xx - 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(85)*_NoGauss1[2]]=   (xx - 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(86)*_NoGauss1[2]]=   (xx + 0.5)*.5*(2.*yy - 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(87)*_NoGauss1[2]]=   (xx + 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(88)*_NoGauss1[2]]=   (xx - 0.5)*.5*(2.*yy + 1.) * .5*zz*(zz + 1.);
      _dphidxx_map1[2][i+(89)*_NoGauss1[2]]=   (- 2.*xx) *.5*(2.*yy - 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(90)*_NoGauss1[2]]=   (xx + 0.5)*    (-2.*yy)    * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(91)*_NoGauss1[2]]=   (- 2.*xx) *.5*(2.*yy + 1.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(92)*_NoGauss1[2]]=   (xx - 0.5)*  (-2.*yy)      * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(93)*_NoGauss1[2]]=   (xx - 0.5)* .5*(2.*yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(94)*_NoGauss1[2]]=   (xx + 0.5)* .5*(2.*yy - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(95)*_NoGauss1[2]]=   (xx + 0.5)* .5*(2.*yy + 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(96)*_NoGauss1[2]]=   (xx - 0.5)* .5*(2.*yy + 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(97)*_NoGauss1[2]]=   (- 2.*xx) * .5*(2.*yy - 1.)* .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(98)*_NoGauss1[2]]=   (xx + 0.5)* (-2.*yy)       * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(99)*_NoGauss1[2]]=   (- 2.*xx) * .5*(2.*yy + 1) * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(100)*_NoGauss1[2]]=  (xx - 0.5)* (-2.*yy)       * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(101)*_NoGauss1[2]]=  (- 2.*xx) * (-2.*yy)       * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(102)*_NoGauss1[2]]=  (- 2.*xx) *.5*(2.*yy - 1.) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(103)*_NoGauss1[2]]=  (xx + 0.5)*  (-2.*yy)      *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(104)*_NoGauss1[2]]=  (- 2.*xx) * .5*(2.*yy + 1) *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(105)*_NoGauss1[2]]=  (xx - 0.5)*  (-2.*yy)      *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(106)*_NoGauss1[2]]=  (- 2.*xx) *  (-2.*yy)      * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(107)*_NoGauss1[2]]=  (- 2.*xx) *  (-2.*yy)      * (1. - zz*zz)   ;
      
      // d2f/dy2    gaussian points                             
      _dphidxx_map1[2][i+(108)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(109)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)*  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(110)*_NoGauss1[2]]=0.*   .5*xx*(xx + 1) *  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(111)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(112)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(113)*_NoGauss1[2]]=0.*   .5*xx*(xx + 1) *  .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(114)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)*  .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(115)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(116)*_NoGauss1[2]]=0.*   (1. - xx*xx)   *  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(117)*_NoGauss1[2]]=0.*   .5*xx*(xx + 1) * (-2.)   * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(118)*_NoGauss1[2]]=0.*    (1. - xx*xx)  *  .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(119)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*   (-2.) * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(120)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(121)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(122)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(123)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(124)*_NoGauss1[2]]=0.*     (1. - xx*xx) * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(125)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)* (-2.)   * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(126)*_NoGauss1[2]]=0.*     (1. - xx*xx) *   .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(127)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)* (-2.)   * .5*zz*(zz + 1) ;
      _dphidxx_map1[2][i+(128)*_NoGauss1[2]]=0.*    (1. - xx*xx)  *  (-2.)  * .5*zz*(zz - 1.);
      _dphidxx_map1[2][i+(129)*_NoGauss1[2]]=0.*     (1. - xx*xx) *   (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(130)*_NoGauss1[2]]=0.*    .5*xx*(xx + 1)*  (-2.)  *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(131)*_NoGauss1[2]]=0.*     (1. - xx*xx) *   (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(132)*_NoGauss1[2]]=0.*   .5*xx*(xx - 1.)*  (-2.)  *  (1. - zz*zz)  ;
      _dphidxx_map1[2][i+(133)*_NoGauss1[2]]=0.*     (1. - xx*xx) *  (-2.) * .5*zz*(zz + 1)  ;
      _dphidxx_map1[2][i+(134)*_NoGauss1[2]]=0.*     (1. - xx*xx) *  (-2.)  * (1. - zz*zz)   ;
      
      // d2f/dyz    gaussian points                             
      _dphidxx_map1[2][i+(135)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(136)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(137)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(138)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(139)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(140)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy - 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(141)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(142)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(143)*_NoGauss1[2]]=   (1. - xx*xx)   * (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(144)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (-2.*yy)   * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(145)*_NoGauss1[2]]=    (1. - xx*xx)  * (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(146)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*   (-2.*yy) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(147)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(148)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(149)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(150)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(151)*_NoGauss1[2]]=     (1. - xx*xx) *(yy - 0.5)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(152)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (-2.*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(153)*_NoGauss1[2]]=     (1. - xx*xx) *(yy + 0.5)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(154)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (-2.*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(155)*_NoGauss1[2]]=    (1. - xx*xx)  *  (-2.*yy)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(156)*_NoGauss1[2]]=     (1. - xx*xx) * (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(157)*_NoGauss1[2]]=    .5*xx*(xx + 1)*  (-2.*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(158)*_NoGauss1[2]]=     (1. - xx*xx) * (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(159)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*  (-2.*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(160)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy)  * .5*(2.*zz + 1)  ;
      _dphidxx_map1[2][i+(161)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy)  *  (-2.*zz)   ;
      
      // d2f/dzx gaussian points
      _dphidxx_map1[2][i+(162)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(163)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(164)*_NoGauss1[2]]=  (xx + 0.5) *.5*yy*(yy + 1) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(165)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy + 1)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(166)*_NoGauss1[2]]=  (xx - 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(167)*_NoGauss1[2]]=  (xx + 0.5)*.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(168)*_NoGauss1[2]]=  (xx + 0.5)* .5*yy*(yy + 1) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(169)*_NoGauss1[2]]=  (xx - 0.5)* .5*yy*(yy + 1) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(170)*_NoGauss1[2]]=  (- 2.*xx) * .5*yy*(yy - 1.)* .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(171)*_NoGauss1[2]]=  (xx + 0.5)* (1. - yy*yy)   * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(172)*_NoGauss1[2]]= (- 2.*xx) * .5*yy*(yy + 1) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(173)*_NoGauss1[2]]= (xx - 0.5)*   (1. - yy*yy) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(174)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy - 1.)*  (-2.*zz)  ;
      _dphidxx_map1[2][i+(175)*_NoGauss1[2]]= (xx + 0.5)* .5*yy*(yy - 1.)*  (-2.*zz)  ;
      _dphidxx_map1[2][i+(176)*_NoGauss1[2]]= (xx + 0.5)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(177)*_NoGauss1[2]]= (xx - 0.5)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(178)*_NoGauss1[2]]= (- 2.*xx) *.5*yy*(yy - 1.) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(179)*_NoGauss1[2]]= (xx + 0.5)* (1. - yy*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(180)*_NoGauss1[2]]= (- 2.*xx) *.5*yy*(yy + 1)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(181)*_NoGauss1[2]]= (xx - 0.5)* (1. - yy*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(182)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(183)*_NoGauss1[2]]=  (- 2.*xx)*.5*yy*(yy - 1.) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(184)*_NoGauss1[2]]= (xx + 0.5)*  (1. - yy*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(185)*_NoGauss1[2]]=  (- 2.*xx)* .5*yy*(yy + 1) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(186)*_NoGauss1[2]]= (xx - 0.5)*  (1. - yy*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(187)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  * .5*(2.*zz + 1)  ;
      _dphidxx_map1[2][i+(188)*_NoGauss1[2]]=  (- 2.*xx)*  (1. - yy*yy)  *  (-2.*zz)   ;
      
      // d2f/dzy    gaussian points                             
      _dphidxx_map1[2][i+(189)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(190)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(191)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(192)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(193)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(194)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (yy - 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(195)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(196)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(197)*_NoGauss1[2]]=   (1. - xx*xx)   * (yy - 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(198)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (-2.*yy)   * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(199)*_NoGauss1[2]]=    (1. - xx*xx)  * (yy + 0.5) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(200)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*   (-2.*yy) * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(201)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(202)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(203)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(204)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(205)*_NoGauss1[2]]=     (1. - xx*xx) *(yy - 0.5)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(206)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (-2.*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(207)*_NoGauss1[2]]=     (1. - xx*xx) *(yy + 0.5)  * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(208)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (-2.*yy)   * .5*(2.*zz + 1) ;
      _dphidxx_map1[2][i+(209)*_NoGauss1[2]]=    (1. - xx*xx)  *  (-2.*yy)  * .5*(2.*zz - 1.);
      _dphidxx_map1[2][i+(210)*_NoGauss1[2]]=     (1. - xx*xx) * (yy - 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(211)*_NoGauss1[2]]=    .5*xx*(xx + 1)*  (-2.*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(212)*_NoGauss1[2]]=     (1. - xx*xx) * (yy + 0.5) *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(213)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*  (-2.*yy)  *  (-2.*zz)  ;
      _dphidxx_map1[2][i+(214)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy)  * .5*(2.*zz + 1)  ;
      _dphidxx_map1[2][i+(215)*_NoGauss1[2]]=     (1. - xx*xx) *  (-2.*yy)  *  (-2.*zz)   ;
     
      // d2f/dz2 gaussian points               
      _dphidxx_map1[2][i+(216)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy - 1.) ;
      _dphidxx_map1[2][i+(217)*_NoGauss1[2]]=    .5*xx*(xx + 1)*.5*yy*(yy - 1.) ;
      _dphidxx_map1[2][i+(218)*_NoGauss1[2]]=   .5*xx*(xx + 1) *.5*yy*(yy + 1)  ;
      _dphidxx_map1[2][i+(219)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy + 1)  ;
      _dphidxx_map1[2][i+(220)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*.5*yy*(yy - 1.) ;
      _dphidxx_map1[2][i+(221)*_NoGauss1[2]]=   .5*xx*(xx + 1) *.5*yy*(yy - 1.) ;
      _dphidxx_map1[2][i+(222)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy + 1) ;
      _dphidxx_map1[2][i+(223)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy + 1) ;
      _dphidxx_map1[2][i+(224)*_NoGauss1[2]]=   (1. - xx*xx)   * .5*yy*(yy - 1.);
      _dphidxx_map1[2][i+(225)*_NoGauss1[2]]=   .5*xx*(xx + 1) * (1. - yy*yy)   ;
      _dphidxx_map1[2][i+(226)*_NoGauss1[2]]=    (1. - xx*xx)  * .5*yy*(yy + 1) ;
      _dphidxx_map1[2][i+(227)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*   (1. - yy*yy) ;
      _dphidxx_map1[2][i+(228)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy - 1.)* (-2.)  ;
      _dphidxx_map1[2][i+(229)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy - 1.)* (-2.)  ;
      _dphidxx_map1[2][i+(230)*_NoGauss1[2]]=    .5*xx*(xx + 1)* .5*yy*(yy + 1) * (-2.)  ;
      _dphidxx_map1[2][i+(231)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* .5*yy*(yy + 1) * (-2.)  ;
      _dphidxx_map1[2][i+(232)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy - 1.) ;
      _dphidxx_map1[2][i+(233)*_NoGauss1[2]]=    .5*xx*(xx + 1)* (1. - yy*yy)   ;
      _dphidxx_map1[2][i+(234)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy + 1)  ;
      _dphidxx_map1[2][i+(235)*_NoGauss1[2]]=   .5*xx*(xx - 1.)* (1. - yy*yy)   ;
      _dphidxx_map1[2][i+(236)*_NoGauss1[2]]=    (1. - xx*xx)  *  (1. - yy*yy)  ;
      _dphidxx_map1[2][i+(237)*_NoGauss1[2]]=     (1. - xx*xx) *.5*yy*(yy - 1.) * (-2.)   ;
      _dphidxx_map1[2][i+(238)*_NoGauss1[2]]=    .5*xx*(xx + 1)*  (1. - yy*yy)  * (-2.)   ;
      _dphidxx_map1[2][i+(239)*_NoGauss1[2]]=     (1. - xx*xx) * .5*yy*(yy + 1) * (-2.)   ;
      _dphidxx_map1[2][i+(240)*_NoGauss1[2]]=   .5*xx*(xx - 1.)*  (1. - yy*yy)  * (-2.)   ;
      _dphidxx_map1[2][i+(241)*_NoGauss1[2]]=     (1. - xx*xx) *  (1. - yy*yy)  ;
      _dphidxx_map1[2][i+(242)*_NoGauss1[2]]=     (1. - xx*xx) *  (1. - yy*yy)  * (-2.)  ;

      
    }//end for loop
  }//end if 3D
#endif                                        
                   
#if ELTYPE==10     
                   
/*               ********************************************
//                               TRI 6  
// 		 ********************************************
//                              2
// 				|\
// 				| \
// 				|  \
// 			   r=l1 |   \ 1-s-r=l0
//				|    \
//				|     \
// 			        |______\
//			       0  s=l2 1
*/
//gaussian coordinates
// 1D --------------------------------
 const double x[3]={-sqrt(3./5.),0.,sqrt(3./5.)};
 
// 2D --------------------------------
 const double alphabeta[4]={1.-2.*( 2./7. + sqrt(15.)/21.),  1.-2.*( 2./7. - sqrt(15.)/21.), 2./7. + sqrt(15.)/21., 2./7. - sqrt(15.)/21.};
 const double r[7]={1./3.,alphabeta[2],alphabeta[0],alphabeta[2],alphabeta[3],alphabeta[1],alphabeta[3]};
 const double s[7]={1./3.,alphabeta[2],alphabeta[2],alphabeta[0],alphabeta[3],alphabeta[3],alphabeta[1]};
 
//3D --------------------------------
// const double p[5]={0.25, 0.5,   1./6., 1./6., 1./6.};   
// const double q[5]={0.25, 1./6., 0.5,   1./6., 1./6.};  
// const double t[5]={0.25, 1./6., 1./6., 0.5,   1./6.};
 
 const double p[14]={0.31088591926330060980, 0.31088591926330060980,   1.-3.*0.31088591926330060980, 0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402,   1.-3.*0.092735250310891226402, 0.092735250310891226402,0.5-0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492};   
 const double q[14]={0.31088591926330060980, 1.-3.*0.31088591926330060980, 0.31088591926330060980,   0.31088591926330060980, 0.092735250310891226402, 1.-3.*0.092735250310891226402, 0.092735250310891226402,   0.092735250310891226402,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492};  
 const double t[14]={0.31088591926330060980, 0.31088591926330060980, 0.31088591926330060980,  1.-3.*0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402,  1.-3.*0.092735250310891226402,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.5-0.045503704125649649492};

//gaussian weights
// 1D --------------------------------
 _weight1[0][0]= 5./9.; _weight1[0][1]= 8./9.; _weight1[0][2]= 5./9.;
 
// 2D --------------------------------
  const double weight[7]={9./80., 31./480. + sqrt(15.)/2400., 31./480. + sqrt(15.)/2400., 
    31./480. + sqrt(15.)/2400.,31./480. - sqrt(15.)/2400., 31./480. - sqrt(15.)/2400., 
    31./480. - sqrt(15.)/2400.}; 
              
 for (int i=0;i<_NoGauss1[1];i++) _weight1[1][i]= weight[i];
 
// 3D --------------------------------
//  const double weight_3[5]={-16./120., 9./120., 9./120., 9./120., 9./120.};
  const double weight_3[14]={0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730};
 for (int i=0;i<_NoGauss1[2];i++) _weight1[2][i]= weight_3[i];


// ****QUADRATIC SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates                  derivatives
//                                           |_______r_______|_______s_______
// phi0= 2.*(1-r-s)*((1.-r-s)-0.5)           | -3.+4.s+4.r   |  -3.+4.*s+4.*r
// phi1= 2.*r*(r-0.5)                        |    4.*r-1     |      0.
// phi2= 2.*s*(s-0.5)                        |    0.         |   4.*s-1
// phi3= 4.*r*(1.-r-s)                       |  4.-4.*s-8.*r |   -4.*r
// phi4= 4.*r*s                              |   4.*s        |   4.*r
// phi5= 4.*s*(1.-r-s)                       |  -4.*s        |   4-4.*r-8.*s                                  

  // 1D --------------------------------
  for (int i=0;i<_NoGauss1[0];i++) {
    const double xi=x[i];
    _phi_map1[0][i]               =.5*xi*(xi - 1.);
    _phi_map1[0][i+_NoGauss1[0]]  =.5*xi*(xi + 1);
    _phi_map1[0][i+2*_NoGauss1[0]]=(1. - xi*xi);
    _dphidxez_map1[0][i]                = xi-.5;
    _dphidxez_map1[0][i+_NoGauss1[0]]   =  xi+.5;
    _dphidxez_map1[0][i+2*_NoGauss1[0]] = -2.*xi;
    //d2/dx2 gaussian points
   _dphidxx_map1[0][i]                = 1.;
   _dphidxx_map1[0][i+_NoGauss1[0]]   = 1.;
   _dphidxx_map1[0][i+2*_NoGauss1[0]] = -2.;
  }
    
  // 2D --------------------------------
  for (int i=0;i<_NoGauss1[1];i++) {
    const double l1=r[i]; const double l2=s[i];
    const double l0 = 1. - l1 - l2;  
    //shape functions
    _phi_map1[1][i+(0)*_NoGauss1[1]]= 2.*l0*(l0-0.5);
    _phi_map1[1][i+(1)*_NoGauss1[1]]= 2.*l1*(l1-0.5); 
    _phi_map1[1][i+(2)*_NoGauss1[1]]= 2.*l2*(l2-0.5);
    _phi_map1[1][i+(3)*_NoGauss1[1]]= 4.*l1*l0;
    _phi_map1[1][i+(4)*_NoGauss1[1]]= 4.*l1*l2;
    _phi_map1[1][i+(5)*_NoGauss1[1]]= 4.*l2*l0;
    
    //d/dr
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]=  -3.+4.*l2+4.*l1; 
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]=  4.*l1-1.; 
    _dphidxez_map1[1][i+(2)*_NoGauss1[1]]=  0.; 
    _dphidxez_map1[1][i+(3)*_NoGauss1[1]]=  4.-4.*l2-8.*l1;
    _dphidxez_map1[1][i+(4)*_NoGauss1[1]]=  4.*l2; 
    _dphidxez_map1[1][i+(5)*_NoGauss1[1]]=  -4.*l2;
    //d/ds
    _dphidxez_map1[1][i+(6)*_NoGauss1[1]]=  -3.+4.*l2+4.*l1; 
    _dphidxez_map1[1][i+(7)*_NoGauss1[1]]=  0.; 
    _dphidxez_map1[1][i+(8)*_NoGauss1[1]]=  4.*l2-1.; 
    _dphidxez_map1[1][i+(9)*_NoGauss1[1]]=  -4.*l1;
    _dphidxez_map1[1][i+(10)*_NoGauss1[1]]= 4.*l1; 
    _dphidxez_map1[1][i+(11)*_NoGauss1[1]]= 4.-4.*l1-8.*l2; 
    
    //d2/dr2
    _dphidxx_map1[1][i+(0)*_NoGauss1[1]]=  4.; 
    _dphidxx_map1[1][i+(1)*_NoGauss1[1]]=  4.; 
    _dphidxx_map1[1][i+(2)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(3)*_NoGauss1[1]]=  -8.;
    _dphidxx_map1[1][i+(4)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(5)*_NoGauss1[1]]=  0.;
    //d2/dsdr
    _dphidxx_map1[1][i+(6)*_NoGauss1[1]]=  4.; 
    _dphidxx_map1[1][i+(7)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(8)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(9)*_NoGauss1[1]]=  -4.;
    _dphidxx_map1[1][i+(10)*_NoGauss1[1]]= 4.; 
    _dphidxx_map1[1][i+(11)*_NoGauss1[1]]= -4.; 
    //d2/drds
    _dphidxx_map1[1][i+(12)*_NoGauss1[1]]=  4.;  
    _dphidxx_map1[1][i+(13)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(14)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(15)*_NoGauss1[1]]=  -4.;
    _dphidxx_map1[1][i+(16)*_NoGauss1[1]]=  4.; 
    _dphidxx_map1[1][i+(17)*_NoGauss1[1]]=  -4.; 
    //d2/ds2
    _dphidxx_map1[1][i+(18)*_NoGauss1[1]]=   4.; 
    _dphidxx_map1[1][i+(19)*_NoGauss1[1]]=   0.; 
    _dphidxx_map1[1][i+(20)*_NoGauss1[1]]=   4.; 
    _dphidxx_map1[1][i+(21)*_NoGauss1[1]]=   0.;
    _dphidxx_map1[1][i+(22)*_NoGauss1[1]]=  0.; 
    _dphidxx_map1[1][i+(23)*_NoGauss1[1]]= -8.; 
  }

  // 3D -----------------------------------------------
  if (_dim == 3) {
    
  for (int i=0;i<_NoGauss1[2];i++) {
    const double l1=p[i]; const double l2=q[i]; const double l3=t[i];
    const double l0 = 1. - l1 - l2 - l3;
    //shape functions
    _phi_map1[2][i+(0)*_NoGauss1[2]]=  l0*(2.*l0 - 1.);
    _phi_map1[2][i+(1)*_NoGauss1[2]]=  l1*(2.*l1 - 1.);
    _phi_map1[2][i+(2)*_NoGauss1[2]]=  l2*(2.*l2 - 1.);
    _phi_map1[2][i+(3)*_NoGauss1[2]]=  l3*(2.*l3 - 1.);
    _phi_map1[2][i+(4)*_NoGauss1[2]]=  4.*l0*l1;
    _phi_map1[2][i+(5)*_NoGauss1[2]]=  4.*l1*l2;
    _phi_map1[2][i+(6)*_NoGauss1[2]]=  4.*l2*l0;
    _phi_map1[2][i+(7)*_NoGauss1[2]]=  4.*l0*l3;
    _phi_map1[2][i+(8)*_NoGauss1[2]]=  4.*l1*l3;
    _phi_map1[2][i+(9)*_NoGauss1[2]]=  4.*l2*l3;
    
    //df/dp
    _dphidxez_map1[2][i+(0)*_NoGauss1[2]]=  -(4.*l0 - 1.);
    _dphidxez_map1[2][i+(1)*_NoGauss1[2]]=  (4.*l1 - 1.);
    _dphidxez_map1[2][i+(2)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(3)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(4)*_NoGauss1[2]]=  4.*(l0 - l1);
    _dphidxez_map1[2][i+(5)*_NoGauss1[2]]=  4.*l2;
    _dphidxez_map1[2][i+(6)*_NoGauss1[2]]=  -4.*l2;
    _dphidxez_map1[2][i+(7)*_NoGauss1[2]]=  -4.*l3;
    _dphidxez_map1[2][i+(8)*_NoGauss1[2]]=  4.*l3;
    _dphidxez_map1[2][i+(9)*_NoGauss1[2]]=  0.;
    //df/dq
    _dphidxez_map1[2][i+(10)*_NoGauss1[2]]=  1.-4.*l0;
    _dphidxez_map1[2][i+(11)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(12)*_NoGauss1[2]]=  4.*l2-1.;
    _dphidxez_map1[2][i+(13)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(14)*_NoGauss1[2]]=  -4.*l1;
    _dphidxez_map1[2][i+(15)*_NoGauss1[2]]=  4.*l1;
    _dphidxez_map1[2][i+(16)*_NoGauss1[2]]=  4.*(l0-l2);
    _dphidxez_map1[2][i+(17)*_NoGauss1[2]]=  -4.*l3;
    _dphidxez_map1[2][i+(18)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(19)*_NoGauss1[2]]=  4.*l3;
    //df/dt
    _dphidxez_map1[2][i+(20)*_NoGauss1[2]]=  1.-4.*l0;
    _dphidxez_map1[2][i+(21)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(22)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(23)*_NoGauss1[2]]=  4.*l3-1.;
    _dphidxez_map1[2][i+(24)*_NoGauss1[2]]=  -4.*l1;
    _dphidxez_map1[2][i+(25)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(26)*_NoGauss1[2]]=  -4.*l2;
    _dphidxez_map1[2][i+(27)*_NoGauss1[2]]=  4.*(l0-l3);
    _dphidxez_map1[2][i+(28)*_NoGauss1[2]]=  4.*l1;
    _dphidxez_map1[2][i+(29)*_NoGauss1[2]]=  4.*l2;
    
    //d2f/dp2
    _dphidxx_map1[2][i+(0)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(1)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(2)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(3)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(4)*_NoGauss1[2]]=  -8.;
    _dphidxx_map1[2][i+(5)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(6)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(7)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(8)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(9)*_NoGauss1[2]]=   0.;
     //d2f/dpq
    _dphidxx_map1[2][i+(10)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(11)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(12)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(13)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(14)*_NoGauss1[2]]= -4.;
    _dphidxx_map1[2][i+(15)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(16)*_NoGauss1[2]]= -4.;
    _dphidxx_map1[2][i+(17)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(18)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(19)*_NoGauss1[2]]=  0.;
    //d2f/dpt
    _dphidxx_map1[2][i+(20)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(21)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(22)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(23)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(24)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(25)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(26)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(27)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(28)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(29)*_NoGauss1[2]]=   0.;
    //d2f/dqp
    _dphidxx_map1[2][i+(30)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(31)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(32)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(33)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(34)*_NoGauss1[2]]= -4.;
    _dphidxx_map1[2][i+(35)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(36)*_NoGauss1[2]]= -4.;
    _dphidxx_map1[2][i+(37)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(38)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(39)*_NoGauss1[2]]=  0.;
    //d2f/dq2
    _dphidxx_map1[2][i+(40)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(41)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(42)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(43)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(44)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(45)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(46)*_NoGauss1[2]]= -8.;
    _dphidxx_map1[2][i+(47)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(48)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(49)*_NoGauss1[2]]=  0.;
    //d2f/dqt
    _dphidxx_map1[2][i+(50)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(51)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(52)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(53)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(54)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(55)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(56)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(57)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(58)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(59)*_NoGauss1[2]]=   4.;
    //d2f/dtp
    _dphidxx_map1[2][i+(60)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(61)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(62)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(63)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(64)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(65)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(66)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(67)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(68)*_NoGauss1[2]]=   4.;
    _dphidxx_map1[2][i+(69)*_NoGauss1[2]]=   0.;
    //d2f/dtq
    _dphidxx_map1[2][i+(70)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(71)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(72)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(73)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(74)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(75)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(76)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(77)*_NoGauss1[2]]=  -4.;
    _dphidxx_map1[2][i+(78)*_NoGauss1[2]]=   0.;
    _dphidxx_map1[2][i+(79)*_NoGauss1[2]]=   4.;
    //d2f/dt2
    _dphidxx_map1[2][i+(80)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(81)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(82)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(83)*_NoGauss1[2]]=  4.;
    _dphidxx_map1[2][i+(84)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(85)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(86)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(87)*_NoGauss1[2]]= -8.;
    _dphidxx_map1[2][i+(88)*_NoGauss1[2]]=  0.;
    _dphidxx_map1[2][i+(89)*_NoGauss1[2]]=  0.;

  }//for loop 3D
  
 }//if 3D  
 
#endif 

  return;
}

/// /// This function generates the Lagrangian linear shape functions
void MGFE::init_lin(
) {// ================================

#if ELTYPE==27

//               ********************************************
//                               QUAD 4  
// 		 ********************************************
//                              
// 			       3 ___________2
// 				|           |
// 			        |           |
//			        |           |
//				|           |
// 			        |___________|
//			       0             1

//gaussian coordinates
 const double a=-sqrt(3./5.); const double b=0.; const double c=-a;
 
 const double x1D[3]={a,b,c};
 
 const double x2D[9]={a,a,a,b,b,b,c,c,c};
 const double y2D[9]={a,b,c,a,b,c,a,b,c};
 
 const double x3D[27]={a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
 const double y3D[27]={a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
 const double z3D[27]={a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};
 
 
//gaussian weights
// 1D --------------------------------
 _weight1[0][0]= 5./9.; _weight1[0][1]= 8./9.; _weight1[0][2]= 5./9.;
 
// 2D --------------------------------
 const double m=_weight1[0][0]*_weight1[0][1]; const double l=_weight1[0][0]*_weight1[0][0];
 const double h=_weight1[0][1]*_weight1[0][1];
 const double weight[9]={l,m,l,m,h,m,l,m,l};
 for (int i=0;i<_NoGauss1[1];i++) _weight1[1][i]= weight[i];
 
// 3D --------------------------------
 const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0]; const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
 const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1]; const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
 const double weight_3[27]={w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
 for (int i=0;i<_NoGauss1[2];i++) _weight1[2][i]= weight_3[i];


// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D                                            derivatives
//                                           |_______x_______|_______y_______
// phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// phi1= x*(1.-y)                            |     1.-y      |      -x
// phi2= x*y                                 |       y       |       x
// phi3= y*(1.-x)                            |      -y       |     1.-x
  
//1D --------------------------------
  for (int i=0;i<_NoGauss1[0];i++) {
   _phi_map1[0][i]=0.5*(1.-x1D[i]);
   _phi_map1[0][i+_NoGauss1[0]]=0.5*(1.+ x1D[i]);
   _dphidxez_map1[0][i]=-0.5;
   _dphidxez_map1[0][i+_NoGauss1[0]]=0.5; 
  }
    
  // 2D --------------------------------
  for (int i=0;i<_NoGauss1[1];i++) {
    const double xx=x2D[i]; const double yy=y2D[i];
    //shape functions
    _phi_map1[1][i+0*_NoGauss1[1]]= 0.5*(1.- xx)*0.5*(1.- yy);
    _phi_map1[1][i+1*_NoGauss1[1]]= 0.5*(1.+ xx)*0.5*(1.- yy);
    _phi_map1[1][i+2*_NoGauss1[1]]= 0.5*(1.+ xx)*0.5*(1.+ yy);
    _phi_map1[1][i+3*_NoGauss1[1]]= 0.5*(1.- xx)*0.5*(1.+ yy);
    
    //d/dx
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]= -0.25*(1.- yy); 
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]=  0.25*(1.- yy);
    _dphidxez_map1[1][i+(2)*_NoGauss1[1]]=  0.25*(1.+ yy);
    _dphidxez_map1[1][i+(3)*_NoGauss1[1]]= -0.25*(1.+ yy);
    //d/dy
    _dphidxez_map1[1][i+(4)*_NoGauss1[1]]= -0.25*(1.- xx); 
    _dphidxez_map1[1][i+(5)*_NoGauss1[1]]= -0.25*(1.+ xx);
    _dphidxez_map1[1][i+(6)*_NoGauss1[1]]=  0.25*(1.+ xx);
    _dphidxez_map1[1][i+(7)*_NoGauss1[1]]=  0.25*(1.- xx);
  }

  // 3D -----------------------------------------------
  if (_dim == 3) {
      
    for (int i=0;i<_NoGauss1[2];i++) { 
     const double xx=x3D[i]; const double yy=y3D[i]; const double zz=z3D[i];
     
     //shape functions
    _phi_map1[2][i+(0)*_NoGauss1[2]]= 0.5*(1.- xx)*0.5*(1.- yy)*0.5*(1.-zz);
    _phi_map1[2][i+(1)*_NoGauss1[2]]= 0.5*(1.+ xx)*0.5*(1.- yy)*0.5*(1.-zz);
    _phi_map1[2][i+(2)*_NoGauss1[2]]= 0.5*(1.+ xx)*0.5*(1.+ yy)*0.5*(1.-zz);
    _phi_map1[2][i+(3)*_NoGauss1[2]]= 0.5*(1.- xx)*0.5*(1.+ yy)*0.5*(1.-zz);
    _phi_map1[2][i+(4)*_NoGauss1[2]]= 0.5*(1.- xx)*0.5*(1.- yy)*0.5*(1.+zz);  
    _phi_map1[2][i+(5)*_NoGauss1[2]]= 0.5*(1.+ xx)*0.5*(1.- yy)*0.5*(1.+zz); 
    _phi_map1[2][i+(6)*_NoGauss1[2]]= 0.5*(1.+ xx)*0.5*(1.+ yy)*0.5*(1.+zz);
    _phi_map1[2][i+(7)*_NoGauss1[2]]= 0.5*(1.- xx)*0.5*(1.+ yy)*0.5*(1.+zz);
    
    //d/dx
    _dphidxez_map1[2][i+(0)*_NoGauss1[2]]= -0.25*(1.- yy)*0.5*(1.-zz); 
    _dphidxez_map1[2][i+(1)*_NoGauss1[2]]=  0.25*(1.- yy)*0.5*(1.-zz); 
    _dphidxez_map1[2][i+(2)*_NoGauss1[2]]=  0.25*(1.+ yy)*0.5*(1.-zz); 
    _dphidxez_map1[2][i+(3)*_NoGauss1[2]]= -0.25*(1.+ yy)*0.5*(1.-zz);
    _dphidxez_map1[2][i+(4)*_NoGauss1[2]]= -0.25*(1.- yy)*0.5*(1.+zz); 
    _dphidxez_map1[2][i+(5)*_NoGauss1[2]]=  0.25*(1.- yy)*0.5*(1.+zz); 
    _dphidxez_map1[2][i+(6)*_NoGauss1[2]]=  0.25*(1.+ yy)*0.5*(1.+zz); 
    _dphidxez_map1[2][i+(7)*_NoGauss1[2]]= -0.25*(1.+ yy)*0.5*(1.+zz);
    //d/dy
    _dphidxez_map1[2][i+(8)*_NoGauss1[2]]= -0.25*(1.- xx)*0.5*(1.-zz);
    _dphidxez_map1[2][i+(9)*_NoGauss1[2]]= -0.25*(1.+ xx)*0.5*(1.-zz);
    _dphidxez_map1[2][i+(10)*_NoGauss1[2]]= 0.25*(1.+ xx)*0.5*(1.-zz);
    _dphidxez_map1[2][i+(11)*_NoGauss1[2]]= 0.25*(1.- xx)*0.5*(1.-zz);
    _dphidxez_map1[2][i+(12)*_NoGauss1[2]]=-0.25*(1.- xx)*0.5*(1.+zz);
    _dphidxez_map1[2][i+(13)*_NoGauss1[2]]=-0.25*(1.+ xx)*0.5*(1.+zz);
    _dphidxez_map1[2][i+(14)*_NoGauss1[2]]= 0.25*(1.+ xx)*0.5*(1.+zz);
    _dphidxez_map1[2][i+(15)*_NoGauss1[2]]= 0.25*(1.- xx)*0.5*(1.+zz);
    //d/dz
    _dphidxez_map1[2][i+(16)*_NoGauss1[2]]= -0.25*(1.- xx)*0.5*(1.- yy);
    _dphidxez_map1[2][i+(17)*_NoGauss1[2]]= -0.25*(1.+ xx)*0.5*(1.- yy);
    _dphidxez_map1[2][i+(18)*_NoGauss1[2]]= -0.25*(1.+ xx)*0.5*(1.+ yy);
    _dphidxez_map1[2][i+(19)*_NoGauss1[2]]= -0.25*(1.- xx)*0.5*(1.+ yy);
    _dphidxez_map1[2][i+(20)*_NoGauss1[2]]=  0.25*(1.- xx)*0.5*(1.- yy);
    _dphidxez_map1[2][i+(21)*_NoGauss1[2]]=  0.25*(1.+ xx)*0.5*(1.- yy);
    _dphidxez_map1[2][i+(22)*_NoGauss1[2]]=  0.25*(1.+ xx)*0.5*(1.+ yy);
    _dphidxez_map1[2][i+(23)*_NoGauss1[2]]=  0.25*(1.- xx)*0.5*(1.+ yy);
  }
  
 }//end if 3D
#endif

#if ELTYPE==10
/* //               ********************************************
//                               TRI 3  
// 		 ********************************************
//                              2
// 				|\
// 				| \
// 				|  \
// 			   r=l1 |   \ 1-s-r=l0
//				|    \
//				|     \
// 			        |______\
//			       0  s=l2 1
*/
//gaussian coordinates
// 1D --------------------------------
 const double x[3]={-sqrt(3./5.),0.,sqrt(3./5.)};
 
// 2D --------------------------------
 const double r[4]={1./3., 1./5., 1./5., 3./5.};
 const double s[4]={1./3., 3./5., 1./5., 1./5.};
 
//3D --------------------------------
 const double p[5]={0.25, 0.5,   1./6., 1./6., 1./6.};   
 const double q[5]={0.25, 1./6., 0.5,   1./6., 1./6.};  
 const double t[5]={0.25, 1./6., 1./6., 0.5,   1./6.};

//gaussian weights
// 1D --------------------------------
 _weight1[0][0]= 5./9.; _weight1[0][1]= 8./9.; _weight1[0][2]= 5./9.;
 
// 2D --------------------------------
 const double weight[4]={-27./48., 25./48., 25./48., 25./48.};
 for (int i=0;i<_NoGauss1[1];i++) _weight1[1][i]= weight[i];
 
// 3D --------------------------------
 const double weight_3[5]={-4./5., 9./20., 9./20., 9./20., 9./20.};
 for (int i=0;i<_NoGauss1[2];i++) _weight1[2][i]= weight_3[i];

 
// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates              derivatives
//                                      |_______r_______|_______s_______
// phi0= 1-r-s                          |      -1.      |      -1.
// phi1= r                              |       1.      |       0.
// phi2= s                              |       0.      |       1.                                 

  // 1D --------------------------------
  for (int i=0;i<_NoGauss1[0];i++) {
    const double xi=x[i];
    _phi_map1[0][i]=1.-xi;
    _phi_map1[0][i+_NoGauss1[0]]= 1.+ xi;
    _dphidxez_map1[0][i]=-1.;
    _dphidxez_map1[0][i+_NoGauss1[0]]=1.;  
  }
    
  // 2D --------------------------------
  for (int i=0;i<_NoGauss1[1];i++) {
    const double l0=1.-r[i]-s[i];
    _phi_map1[1][i+(0)*_NoGauss1[1]]= l0;
    _phi_map1[1][i+(1)*_NoGauss1[1]]= r[i]; 
    _phi_map1[1][i+(2)*_NoGauss1[1]]= s[i];
    
    //d/dr
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]=  -1.; 
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]=  1.; 
    _dphidxez_map1[1][i+(2)*_NoGauss1[1]]=  0.; 
    //d/ds
    _dphidxez_map1[1][i+(3)*_NoGauss1[1]]=  -1.;
    _dphidxez_map1[1][i+(4)*_NoGauss1[1]]=  0.; 
    _dphidxez_map1[1][i+(5)*_NoGauss1[1]]=  1.;
  }

  // 3D -----------------------------------------------
  if (_dim == 3) {
    for (int i=0;i<_NoGauss1[2];i++) {
    const double l0=1.-p[i]-q[i]-t[i]; 
    _phi_map1[2][i+(0)*_NoGauss1[2]]= l0;
    _phi_map1[2][i+(1)*_NoGauss1[2]]= p[i]; 
    _phi_map1[2][i+(2)*_NoGauss1[2]]= q[i];
    _phi_map1[2][i+(3)*_NoGauss1[2]]= t[i];
    
    //d/dp
    _dphidxez_map1[2][i+(0)*_NoGauss1[2]]= -1.; 
    _dphidxez_map1[2][i+(1)*_NoGauss1[2]]=  1.; 
    _dphidxez_map1[2][i+(2)*_NoGauss1[2]]=  0.; 
    _dphidxez_map1[2][i+(3)*_NoGauss1[2]]=  0.;
    //d/dq
    _dphidxez_map1[2][i+(4)*_NoGauss1[2]]= -1.; 
    _dphidxez_map1[2][i+(5)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(6)*_NoGauss1[2]]=  1.; 
    _dphidxez_map1[2][i+(7)*_NoGauss1[2]]=  0.; 
    //d/dt
    _dphidxez_map1[2][i+(8)*_NoGauss1[2]]= -1.; 
    _dphidxez_map1[2][i+(9)*_NoGauss1[2]]=  0.;
    _dphidxez_map1[2][i+(10)*_NoGauss1[2]]= 0.; 
    _dphidxez_map1[2][i+(11)*_NoGauss1[2]]= 1.;
  }
                                   
 }//if 3D
#endif
  return;
}


/// /// This function generates the Lagrangian piecewise shape functions
void MGFE::init_pie(
) {// ================================

  
 //               ********************************************
//                               S_0 
//               ********************************************
//                              
//                               ___________
//                              |           |
//                              |     0     |
//                              |     x     |
//                              |           |
//                              |___________|
//                                          

//gaussian coordinates
 const double a=-sqrt(3./5.); const double b=0.; const double c=-a;
 
 const double x1D[3]={a,b,c};
 
 const double x2D[9]={a,a,a,b,b,b,c,c,c};
 const double y2D[9]={a,b,c,a,b,c,a,b,c};
 
 const double x3D[27]={a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
 const double y3D[27]={a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
 const double z3D[27]={a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};
 
 
//gaussian weights
// 1D --------------------------------
 _weight1[0][0]= 5./9.; _weight1[0][1]= 8./9.; _weight1[0][2]= 5./9.;
 
// 2D --------------------------------
 const double m=_weight1[0][0]*_weight1[0][1]; const double l=_weight1[0][0]*_weight1[0][0];
 const double h=_weight1[0][1]*_weight1[0][1];
 const double weight[9]={l,m,l,m,h,m,l,m,l};
 for (int i=0;i<_NoGauss1[1];i++) _weight1[1][i]= weight[i];
 
// 3D --------------------------------
 const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0]; const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
 const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1]; const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
 const double weight_3[27]={w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
 for (int i=0;i<_NoGauss1[2];i++) _weight1[2][i]= weight_3[i];


// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D                                            derivatives
//                                           |_______x_______|_______y_______
// phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// phi1= x*(1.-y)                            |     1.-y      |      -x
// phi2= x*y                                 |       y       |       x
// phi3= y*(1.-x)                            |      -y       |     1.-x
  
//1D --------------------------------
  for (int i=0;i<_NoGauss1[0];i++) {
   _phi_map1[0][i]=1.;
   _dphidxez_map1[0][i]=0.;
   
   if(NDOF_K>1) {
     _phi_map1[0][i+_NoGauss1[0]]= x1D[i]; 
     _dphidxez_map1[0][i+_NoGauss1[0]]=1.; 
   }
  }
    
  // 2D --------------------------------
  for (int i=0;i<_NoGauss1[1];i++) {
    const double xx=x2D[i]; const double yy=y2D[i];
    //shape functions
    _phi_map1[1][i+0*_NoGauss1[1]]=1.;
    // derivatives
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]= 0.;  //d/dx
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]= 0.;  //d/dy
    
     if(NDOF_K>1) {
    _phi_map1[1][i+0*_NoGauss1[1]]=1.;
    _phi_map1[1][i+1*_NoGauss1[1]]=xx;
    _phi_map1[1][i+2*_NoGauss1[1]]=yy;
    
    //d/dx
    _dphidxez_map1[1][i+(0)*_NoGauss1[1]]= 0.; 
    _dphidxez_map1[1][i+(1)*_NoGauss1[1]]= 1.;
    _dphidxez_map1[1][i+(2)*_NoGauss1[1]]= 0.;
    //d/dy
    _dphidxez_map1[1][i+(3)*_NoGauss1[1]]= 0.; 
    _dphidxez_map1[1][i+(4)*_NoGauss1[1]]= 0.;
    _dphidxez_map1[1][i+(5)*_NoGauss1[1]]= 1.; 
       
     }
    
  }

  // 3D -----------------------------------------------
  if (_dim == 3) {
      
    for (int i=0;i<_NoGauss1[2];i++) { 
     const double xx=x3D[i]; const double yy=y3D[i]; const double zz=z3D[i];
     
     //shape functions
    _phi_map1[2][i+(0)*_NoGauss1[2]]= 1.;
    // derivatives
    _dphidxez_map1[2][i+(0)*_NoGauss1[2]]= 0.;  //d/dx 
    _dphidxez_map1[2][i+(1)*_NoGauss1[2]]= 0.;  //d/dy
    _dphidxez_map1[2][i+(2)*_NoGauss1[2]]= 0.; //d/dz
    
    
     if(NDOF_K>1) {
         //shape functions
    _phi_map1[2][i+(0)*_NoGauss1[2]]= 1.;
    _phi_map1[2][i+(1)*_NoGauss1[2]]= xx;
    _phi_map1[2][i+(2)*_NoGauss1[2]]= yy;
    _phi_map1[2][i+(3)*_NoGauss1[2]]= zz;
    
    //d/dx
    _dphidxez_map1[2][i+(0)*_NoGauss1[2]]= 0.; 
    _dphidxez_map1[2][i+(1)*_NoGauss1[2]]= 1.; 
    _dphidxez_map1[2][i+(2)*_NoGauss1[2]]= 0.; 
    _dphidxez_map1[2][i+(3)*_NoGauss1[2]]= 0.;
    //d/dy
    _dphidxez_map1[2][i+(4)*_NoGauss1[2]]= 0.;
    _dphidxez_map1[2][i+(5)*_NoGauss1[2]]= 0.;
    _dphidxez_map1[2][i+(6)*_NoGauss1[2]]= 1.;
    _dphidxez_map1[2][i+(7)*_NoGauss1[2]]= 0.;
    //d/dz
    _dphidxez_map1[2][i+(8)*_NoGauss1[2]]= 0.;
    _dphidxez_map1[2][i+(9)*_NoGauss1[2]]= 0.;
    _dphidxez_map1[2][i+(10)*_NoGauss1[2]]=0.;
    _dphidxez_map1[2][i+(11)*_NoGauss1[2]]= 1.;
       
     }
    
    
  }
  
 }//end if 3D 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
// //gaussian weights
//  for(int i=0;i<_NoGauss1[0];i++) _weight1[0][i] = 2.;
//  for(int i=0;i<_NoGauss1[1];i++) _weight1[1][i] = 4.;
//  for(int i=0;i<_NoGauss1[2];i++) _weight1[2][i] = 8.;
// 
// //shape functions and derivatives
//  for(int j=0;j<_dim;j++) for(int i=0;i<_NoGauss1[j];i++)           _phi_map1[j][i]= 1.;     //set all to zero
//  for(int j=0;j<_dim;j++) for(int i=0;i<_NoGauss1[j]*(j+1);i++) _dphidxez_map1[j][i]=0.;

  return;
}




// ================================
/// This function writes FEM
void MGFE::write(
  const std::string& name,// file <-
  const int kdim         // dimension <-
) {// ================================
  std::ofstream in(name.c_str());
  this->write_c(in,kdim);
}

// ====================================================
/// This function writes shape and derivative values at the gaussian points
void MGFE::write_c(
  std::ostream& out,    // file <-
  const int kdim       // dimension  <-
) { // ================================================

  if (!out) {std::cout<<" Gauss Outfile "<<out<<" not opened."<<std::endl;exit(3);}

  // heading
  out<< "gpoints" <<  _NoGauss1[kdim-1]  << std::endl  ;
  // max_nodes (decreasing level order)
  for (int k=0;k<_NoGauss1[kdim-1];k++) {
    out<< " weight " << std::setprecision(20)<<  _weight1[kdim-1][k] << std::endl;
    for (int s=0;s<_NoShape[kdim-1];s++) {
      out<<std::setprecision(20)<<_phi_map1[kdim-1][s*_NoGauss1[kdim-1]+k]<<" ";
      for (int idim=0;idim<3;idim++) {
        out<<std::setprecision(20)<<
        _dphidxez_map1[kdim-1][(idim*_NoShape[kdim-1]+s)*_NoGauss1[kdim-1]+k] << "  ";
      }
      out<< "  \n";
    }
    out<< "  \n";
  }
#ifdef PRINT_INFO
  std::cout << " MGFE::write_c:  dim " << kdim<< " fem with " << _NoGauss1[kdim-1] << " gaussian points and "
            << _NoShape[kdim-1] << " shape functions \n" << std::endl;
#endif
  return;
}


// ===========================================================
/// This function read the gaussian points
void MGFE::read_c(
  std::istream& infile, // file <-
  const int kdim       // dimension <-
) {// =========================================================


  // file
  if (!infile) {std::cout << " Read_c3:Gauss Input file not opened \n";std::exit(3);}

  const int  bufLen = 30; char  buf[bufLen+1];sprintf(buf,"0");
  int ng;double dummy;
  int ngss=_NoGauss1[kdim-1]; int nsh=_NoShape[kdim-1];

  // Reading file
  while (strncmp(buf,"gpoints",7) != 0)    infile >> buf;
  // # of gaussian points
  infile >> ng;
  if (_NoGauss1[kdim-1] != ng) std::cout <<
                                           " MGFE::read_c: error _NoGauss is " << _NoGauss1[kdim-1] << std::endl;

  // Reading 3d shape functions at guassian points 
					  
  for (int k=0;k<ngss;k++) {             // gauss points
    infile >> buf; infile >> dummy; _weight1[kdim-1][k]=dummy;
//         std::cout<<" "<< kdim-1 << " " << k << " " << dummy<<"\n";
    for (int s=0;s< nsh;s++) {           // shape functions
      infile >> dummy; _phi_map1[kdim-1][s*ngss+k]=dummy;

      for (int idim=0;idim<kdim;idim++) { // derivatives
        infile >> dummy;  _dphidxez_map1[kdim-1][(idim*nsh+s)*ngss+k]=dummy;
      }
    }
  }
  std::cout<<"\n\n**************************/*phi*/*******************************\n\n";
  for (int k=0;k<ngss*nsh;k++) {
//   if (k%27==0) std::cout<< "phi" << k/27  << "\n\n";
    std::cout/*<< k<< " "*/ << _phi_map1[2][k] << " \n ";
//     std::cout <<  _weight1[2][k] << "\n";
  }
  
//   std::cout<<"\n\n**************************/*DERIVATA*/*******************************\n\n";
//   for (int i=0;i<ngss*nsh*3;i++) {	
//        if (i%27==0) std::cout<< "dphi" << i/27  << "\n\n";
//     std::cout<<" "<< i<< " "<< " " <<_dphidxez_map1[2][i]<<"\n";
//   }

#ifdef PRINT_INFO
  std::cout << " MGFE::read_c:  dim " << kdim<< " fem with " << ngss << " gaussian points and "
            << nsh << " shape functions " << std::endl;
#endif

  return;
}

// ========================================
/// This function computes the  derivatives at
/// the gaussian point ng:
double MGFE::Jac(//  jacobean ->
  const int ng, // gaussian point <-
  double x[],     // coordinates  <-
  double InvJac[] // Inverted Jacobean ->
) {// =====================================
  double Jac=0.;
#if  DIMENSION==1
  Jac = Jac1D(ng,x,InvJac);
#endif
#if  DIMENSION==2
  Jac = Jac2D(ng,x,InvJac);
#endif
#if  DIMENSION==3
  Jac = Jac3D(ng,x,InvJac);
#endif
  return Jac;
}

// ========================================
/// This function computes the 3D derivatives at
/// the gaussian point ng:
double MGFE::Jac3D(// 3D  jacobean ->
  const int ng,   // gaussian point <-
  double xyz[],     // coordinates <-
  double InvJac[] // Jacobean ->
) { // =====================================

  double x_xi=0., x_eta =0.,x_zeta =0.;
  double y_xi =0.,y_eta=0., y_zeta =0.;
  double z_xi=0., z_eta =0.,z_zeta =0.;
  int nshape=_NoShape[2];
  int offset=nshape*_NoGauss1[2];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[2]+ng;
    double dphidx=_dphidxez_map1[2][sng];
    double dphideta=_dphidxez_map1[2][sng+offset];
    double dphidzeta=_dphidxez_map1[2][sng+2*offset];

    x_xi  += xyz[s]*dphidx;
    x_eta += xyz[s]*dphideta;
    x_zeta +=xyz[s]*dphidzeta;
    y_xi  += xyz[s+NDOF_FEM]*dphidx;
    y_eta += xyz[s+NDOF_FEM]*dphideta;
    y_zeta +=xyz[s+NDOF_FEM]*dphidzeta;
    z_xi  += xyz[s+2*NDOF_FEM]*dphidx;
    z_eta += xyz[s+2*NDOF_FEM]*dphideta;
    z_zeta +=xyz[s+2*NDOF_FEM]*dphidzeta;

  }
  double det=x_xi*(y_eta*z_zeta-y_zeta*z_eta)-
             x_eta*(y_xi*z_zeta-y_zeta*z_xi)+
             x_zeta*(y_xi*z_eta-y_eta*z_xi);
  double invdet=1./det;
  InvJac[0]=(y_eta*z_zeta-y_zeta*z_eta)*invdet;
  InvJac[3]=-(x_eta*z_zeta-x_zeta*z_eta)*invdet;//
  InvJac[6]=(y_zeta*x_eta-x_zeta*y_eta)*invdet;
  InvJac[1]=-(y_xi*z_zeta-y_zeta*z_xi)*invdet;//
  InvJac[4]=(x_xi*z_zeta-x_zeta*z_xi)*invdet;
  InvJac[7]=-(y_zeta*x_xi-y_xi*x_zeta)*invdet;//
  InvJac[2]=(z_eta*y_xi-y_eta*z_xi)*invdet;
  InvJac[5]=-(x_xi*z_eta-x_eta*z_xi)*invdet;//
  InvJac[8]=(x_xi*y_eta-y_xi*x_eta)*invdet;


  return(det);
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac1D(// 2D  jacobean ->
  const int ng,   // Gaussian point <-
  double x[],       // coordinates  <-
  double InvJac[] // Jacobean
)  {// ===================================================

  double x_xi=0.;
  int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss1[0];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[0]+ng;
    x_xi  +=x[s]*_dphidxez_map1[0][sng];
//     x_eta +=x[s]*_dphidxez_map1[1][sng+offset];
//     y_xi  +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng];
//     y_eta +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng+offset];
  }
  double det=x_xi;double idet=1./det;
//   InvJac[0]=y_eta*idet;    // dxi dx
//   InvJac[1]=-y_xi*idet;    //deta dx
//   InvJac[2]=-x_eta*idet;   // dxi dy
  InvJac[0]=x_xi*idet;     //deta dy

  return(det);
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac2D(// 2D  jacobean ->
  const int ng,   // Gaussian point <-
  double x[],       // coordinates  <-
  double InvJac[] // Inverted Jacobean
)  {// ===================================================

  double x_xi=0., x_eta =0.,y_xi =0., y_eta=0.;
  int nshape=_NoShape[1];
  int offset=_NoShape[1]*_NoGauss1[1];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[1]+ng;
    x_xi  +=x[s]*_dphidxez_map1[1][sng];
    x_eta +=x[s]*_dphidxez_map1[1][sng+offset];
    y_xi  +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng];
    y_eta +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng+offset];
  }
  double det=x_xi*y_eta-y_xi*x_eta;double idet=1./det;
  InvJac[0]=y_eta*idet;    // dxi dx
  InvJac[1]=-y_xi*idet;    //deta dx
  InvJac[2]=-x_eta*idet;   // dxi dy
  InvJac[3]=x_xi*idet;     //deta dy

  return(det);
}

// ========================================
/// This function computes the  derivatives at
/// the nodal point:
double MGFE::Jac_nodes(//  jacobean ->
  const int ng, // nodal point <-
  double x[],     // coordinates  <-
  double InvJac[] // Jacobean ->
) {// =====================================
  double Jac=0.;
#if  DIMENSION==1
  Jac = Jac1D_nodes(ng,x,InvJac);
#endif
#if  DIMENSION==2
  Jac = Jac2D_nodes(ng,x,InvJac);
#endif
#if  DIMENSION==3
  Jac = Jac3D_nodes(ng,x,InvJac);
#endif
  return Jac;
}

// ========================================
/// This function computes the 3D derivatives at
/// the nodal point ng:
double MGFE::Jac3D_nodes(// 3D  jacobean ->
  const int ng,   // nodal point <-
  double xyz[],     // coordinates <-
  double InvJac[] // Jacobean ->
) { // =====================================

  double x_xi=0., x_eta =0.,x_zeta =0.;
  double y_xi =0.,y_eta=0., y_zeta =0.;
  double z_xi=0., z_eta =0.,z_zeta =0.;
  int nshape=_NoShape[2];
  int offset=nshape*_NoGauss1[2];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[2]+ng;
    double dphidx=_dphidxez_map1_nodes[2][sng];
    double dphideta=_dphidxez_map1_nodes[2][sng+offset];
    double dphidzeta=_dphidxez_map1_nodes[2][sng+2*offset];

    x_xi  += xyz[s]*dphidx;
    x_eta += xyz[s]*dphideta;
    x_zeta +=xyz[s]*dphidzeta;
    y_xi  += xyz[s+NDOF_FEM]*dphidx;
    y_eta += xyz[s+NDOF_FEM]*dphideta;
    y_zeta +=xyz[s+NDOF_FEM]*dphidzeta;
    z_xi  += xyz[s+2*NDOF_FEM]*dphidx;
    z_eta += xyz[s+2*NDOF_FEM]*dphideta;
    z_zeta +=xyz[s+2*NDOF_FEM]*dphidzeta;

  }
  double det=x_xi*(y_eta*z_zeta-y_zeta*z_eta)-
             x_eta*(y_xi*z_zeta-y_zeta*z_xi)+
             x_zeta*(y_xi*z_eta-y_eta*z_xi);
  double invdet=1./det;
  InvJac[0]=(y_eta*z_zeta-y_zeta*z_eta)*invdet;
  InvJac[3]=-(x_eta*z_zeta-x_zeta*z_eta)*invdet;//
  InvJac[6]=(y_zeta*x_eta-x_zeta*y_eta)*invdet;
  InvJac[1]=-(y_xi*z_zeta-y_zeta*z_xi)*invdet;//
  InvJac[4]=(x_xi*z_zeta-x_zeta*z_xi)*invdet;
  InvJac[7]=-(y_zeta*x_xi-y_xi*x_zeta)*invdet;//
  InvJac[2]=(z_eta*y_xi-y_eta*z_xi)*invdet;
  InvJac[5]=-(x_xi*z_eta-x_eta*z_xi)*invdet;//
  InvJac[8]=(x_xi*y_eta-y_xi*x_eta)*invdet;


  return(det);
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac1D_nodes(// 2D  jacobean ->
  const int ng,   // Gaussian point <-
  double x[],       // coordinates  <-
  double InvJac[] // Jacobean
)  {// ===================================================

  double x_xi=0.;
  int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss1[0];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[0]+ng;
    x_xi  +=x[s]*_dphidxez_map1_nodes[0][sng];
  }
  double det=x_xi;double idet=1./det;
  InvJac[0]=x_xi*idet;     //deta dy
  return(det);
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac2D_nodes(// 2D  jacobean ->
  const int ng,   // Gaussian point <-
  double x[],       // coordinates  <-
  double InvJac[] // Jacobean
)  {// ===================================================

  double x_xi=0., x_eta =0.,y_xi =0., y_eta=0.;
  int nshape=_NoShape[1];
  int offset=_NoShape[1]*_NoGauss1[1];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[1]+ng;
    x_xi  +=x[s]*_dphidxez_map1_nodes[1][sng];
    x_eta +=x[s]*_dphidxez_map1_nodes[1][sng+offset];
    y_xi  +=x[s+NDOF_FEM]*_dphidxez_map1_nodes[1][sng];
    y_eta +=x[s+NDOF_FEM]*_dphidxez_map1_nodes[1][sng+offset];
  }
  double det=x_xi*y_eta-y_xi*x_eta;double idet=1./det;
  InvJac[0]=y_eta*idet;    // dxi dx
  InvJac[1]=-y_xi*idet;    //deta dx
  InvJac[2]=-x_eta*idet;   // dxi dy
  InvJac[3]=x_xi*idet;     //deta dy

  return(det);
}


// ==========================================
double MGFE::JacSur(// surface jacobean ->
  const int ng,    // gaussian point
  double x[],       // coordinates
  double InvJac[]  // Jacobean ->
) const { // ================================

  double JacSur=0.;
#if  DIMENSION==1
  JacSur = JacSur1D(ng,x,InvJac);
#endif
#if  DIMENSION==2
  JacSur = JacSur2D(ng,x,InvJac);
#endif
#if  DIMENSION==3
  JacSur = JacSur3D(ng,x,InvJac);
#endif
  return JacSur;
}
// =====================================
/// This function compute the line Jacobian at the gaussian point ng
double MGFE::JacSur1D( // 1D surface jacobean ->
  const int ng,       // gaussian point  <-
  double x[],       // coordinates <-
  double InvJac[]  // Jacobean ->
) const { // ==========================

  // Values to compute at gaussian points
  double dxdxi=0.;//
//   double dydxi=0.;// d(x,y)d(xi,eta)
  const int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss[0];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[0]+ng;
    double dphidxi=_dphidxez_map1[0][sng];
    dxdxi += x[s]*dphidxi;
//     dydxi += x[s+NDOF_FEMB]*dphidxi;
  }
  //surface weighted jacobean
  double det=sqrt(dxdxi*dxdxi);
  return(det);
}
// // =====================================
// /// This function compute the line Jacobian at the gaussian point ng
// double MGFE::JacSur2D( // 2D surface jacobean ->
//   const int ng,       // gaussian point  <-
//   double x[],       // coordinates <-
//   double InvJac[]  // Jacobean ->
// ) const { // ==========================
// 
//   // Values to compute at gaussian points
//   double dxdxi=0.;//
//   double dydxi=0.;// d(x,y)d(xi,eta)
//   const int nshape=_NoShape[0];
// //   int offset=_NoShape[0]*_NoGauss[0];
//   for (int s=0;s<nshape;s++) {
//     int sng=s*_NoGauss1[0]+ng;
//     double dphidxi=_dphidxez_map1[0][sng];
//     dxdxi += x[s]*dphidxi;
//     dydxi += x[s+NDOF_FEMB]*dphidxi;
//   }
//   //surface weighted jacobean
//   double det=sqrt(dxdxi*dxdxi+dydxi*dydxi);
//   return(det);
// }

// =====================================
/// This function compute the line Jacobian at the gaussian point ng
double MGFE::JacSur2D( // 2D surface jacobean ->
  const int ng,       // gaussian point  <-
  double x[],       // coordinates <-
  double InvJac[]  // Jacobean ->
) const { // ==========================

  // Values to compute at gaussian points
  double dxdxi=0.;//
  double dydxi=0.;// d(x,y)d(xi,eta)
 const int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss[0];
  for (int s=0;s<nshape;s++) {
    int sng=s*_NoGauss1[0]+ng;
    double dphidxi=_dphidxez_map1[0][sng];
    dxdxi += x[s]*dphidxi; 
    dydxi += x[s+NDOF_FEMB]*dphidxi; 
  }
  //It functions only on straight boundary edge!!!!!!!!!!!!!!!!!!!! 
  InvJac[0] = 1./(dxdxi+dydxi); //porcata!!!!!!!!!!!
  //surface weighted jacobean
  double det=sqrt(dxdxi*dxdxi+dydxi*dydxi);
  return(det);
}

// ===============================================
/// This function computes the surface Jacobian at the gaussian point ng:
// double MGFE::JacSur3D(
//   const int ng,// gaussian point <-
//   double x[],    // coordinates <-
//   double InvJac[] // Jacobean
// ) const {// ======================================
//
//   // Values to compute at gaussian points
//   double dxdxi=0.;double dxdeta=0.;//
//   double dydxi=0.;double dydeta=0.;// d(x,y,z)d(xi,eta)
//   double dzdxi=0.;double dzdeta=0.;//
//
//   int nshape=_NoShape[1];
//   int offset=_NoShape[1]*_NoGauss1[1];
//
//   for (int s=0;s<(int)_NoShape[1];s++) {
//
//     int sng=s*_NoGauss1[1]+ng;
//     double dphidxi=_dphidxez_map1[1][sng];
//     double dphideta=_dphidxez_map1[1][sng+offset];
//
//   //  double dphidzeta=_dphidxez3D_map[sng+2*offset];
//   //  dxdxi += x[s]*dphidxi; dxdeta += x[s]*dphideta;
//   //  dydxi += y[s+nshape]*dphidxi; dydeta += y[s+nshape]*dphideta;
//   //  dzdxi += z[s+2*nshape]*dphidxi; dzdeta += z[s+2*nshape]*dphideta;
//
//     dxdxi += x[s]*dphidxi;             dxdeta += x[s]*dphideta;
//     dydxi += x[s+NDOF_FEMB]*dphidxi;   dydeta += x[s+NDOF_FEMB]*dphideta;
//     dzdxi += x[s+2*NDOF_FEMB]*dphidxi; dzdeta += x[s+2*NDOF_FEMB]*dphideta;
//   }
//   //surface weighted jacobean
//
//
//   double det=sqrt((dxdxi*dydeta-dxdeta*dydxi)*(dxdxi*dydeta-dxdeta*dydxi)+
//                   (dydxi*dzdeta-dydeta*dzdxi)*(dydxi*dzdeta-dydeta*dzdxi)+
//                   (dzdxi*dxdeta-dzdeta*dxdxi)*(dzdxi*dxdeta-dzdeta*dxdxi)
//                  );
//   return(det);
// }

// ===============================================
/// This function computes the surface Jacobian at the gaussian point ng:
double MGFE::JacSur3D(
  const int ng,// gaussian point <-
  double x[],       // coordinates <-
  double InvJac[]  // Jacobean ->
) const {// ======================================

  // Values to compute at gaussian points
  double dxdxi=0.;double dxdeta=0.;//
  double dydxi=0.;double dydeta=0.;// d(x,y,z)d(xi,eta)
  double dzdxi=0.;double dzdeta=0.;//

//   int nshape=_NoShape[1];
  int offset=_NoShape[1]*_NoGauss1[1];

  for (int s=0;s<(int)_NoShape[1];s++) {

    int sng=s*_NoGauss1[1]+ng;
    double dphidxi=_dphidxez_map1[1][sng];
    double dphideta=_dphidxez_map1[1][sng+offset];

    //  double dphidzeta=_dphidxez3D_map[sng+2*offset];
    //  dxdxi += x[s]*dphidxi; dxdeta += x[s]*dphideta;
    //  dydxi += y[s+nshape]*dphidxi; dydeta += y[s+nshape]*dphideta;
    //  dzdxi += z[s+2*nshape]*dphidxi; dzdeta += z[s+2*nshape]*dphideta;

    dxdxi += x[s]*dphidxi;             dxdeta += x[s]*dphideta;
    dydxi += x[s+NDOF_FEMB]*dphidxi;   dydeta += x[s+NDOF_FEMB]*dphideta;
    dzdxi += x[s+2*NDOF_FEMB]*dphidxi; dzdeta += x[s+2*NDOF_FEMB]*dphideta;
  }
  //surface weighted jacobean
  double det=sqrt((dxdxi*dydeta-dxdeta*dydxi)*(dxdxi*dydeta-dxdeta*dydxi)+
                  (dydxi*dzdeta-dydeta*dzdxi)*(dydxi*dzdeta-dydeta*dzdxi)+
                  (dzdxi*dxdeta-dzdeta*dxdxi)*(dzdxi*dxdeta-dzdeta*dxdxi)
                 );
  double idet = 1./det;
  InvJac[0]=dydeta*idet;    // dxi dx
  InvJac[1]=-dydxi*idet;    //deta dx
  InvJac[2]=-dxdeta*idet;   // dxi dy
  InvJac[3]=dxdxi*idet;     //deta dy
  return(det);
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
void MGFE::get_phi_gl_g(
  const int kdim, // dimension <-
  const int qp,   // gaussian point <-
  double phi[]     // shape functions ->
) {// =================================================
  for (int ish=0; ish<_NoShape[kdim-1]; ish++) {
    phi[ish] = _phi_map1[kdim-1][ish*_NoGauss1[kdim-1]+qp];
  }
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
void MGFE::get_phi_gl_g(
  const int kdim,                       // dimension        <-
  const int qp,                         // gaussian point   <-
  std::vector<double> & phi             // shape functions  ->
) {// =====================================
  for (int ish=0;ish<_NoShape[kdim-1];ish++) {
    phi[ish] = _phi_map1[kdim-1][ish*_NoGauss1[kdim-1]+qp];
  }
  return;
}

// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ] 
void   MGFE::get_dphi_gl_g(
  const int kdim,// dimension <-
  const int qp,  // gaussian point <
  const double InvJac[], // Jacobean
  double dphi[]   // global derivatives ->
) {// =========================================

  const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
  const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
  const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset

  double gradphi_g[DIMENSION]; // temp grad phi

  for (int eln=0; eln<el_nnodes; eln++)    {
    int lqp=eln* el_ngauss+qp;
    for (int idim=0;idim<_dim; idim++) gradphi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
    for (int idim=0;idim<_dim; idim++) {
      double sum = 0.;
      for (int jdim=0; jdim<_dim; jdim++) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
      dphi[eln+idim* el_nnodes] = sum;
    }
  }
  return;
}
// =============================================================
/// Shape functions 2nd derivatives at the gauss point qp
/// ddphi[id*_dim*_dim+i*_dim+j]= Tij
///    Txx Txy Txz
///    Tyx Tyy Tyz
///    Tzx Tzy Tzz
///   
///     ddphi[Txx(0),Txy(0),Txz(0),Tyx(0),Tyy(0),Tyz(0),Tzx(0),Tzy(0),Tzz(0),
///                 ....................................
///                 Txx(id),Txy(id),Txz,Tyx(id),Tyy,Tyz,Tzx(id),Tzy(id),Tzz(id)]
///   tensor  order  at the gauss point qp, node inode
void   MGFE::get_ddphi_gl_g(
  const int kdim,// dimension <-
  const int qp,  // gaussian point <
  const double InvJac[], // Jacobean
  double ddphi[]   // global derivatives ->
) {// =========================================

  const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
  const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
  const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset
  double ddphi_loc[NDOF_FEM*DIMENSION*DIMENSION];

  for (int eln=0; eln<el_nnodes; eln++) {
    const int lqp=eln*el_ngauss+qp;
    for (int idim=0;idim<_dim*_dim; idim++) 
      ddphi_loc[eln*_dim*_dim+idim] = _dphidxx_map1[kdim-1][lqp+idim*goffset]; 
  
  
  double Hess_tmp[DIMENSION*DIMENSION],  JacTf[DIMENSION*DIMENSION];
  for (int init=0;init<_dim; init++) for (int jnit=0;jnit<_dim; jnit++) {
    Hess_tmp[init*_dim+jnit]=0.; ddphi[eln*_dim*_dim+init*_dim+jnit]=0.;
    JacTf[init*_dim+jnit] = InvJac[jnit*_dim+init];
   }
   
   for(int ii=0;ii<_dim; ii++) for(int jj=0;jj<_dim; jj++)  for(int ss=0;ss<_dim; ss++)
    Hess_tmp[ii*_dim+jj] += ddphi_loc[eln*_dim*_dim+ii*_dim+ss]*JacTf[ss*_dim+jj];
	  
    for(int ii=0;ii<_dim; ii++) for(int jj=0;jj<_dim; jj++) for(int ss=0;ss<_dim; ss++)
      ddphi[eln*_dim*_dim+ii*_dim+jj] += InvJac[ii*_dim+ss]*Hess_tmp[ss*_dim+jj];
  }
  //Pay attention, different order with respect to first order derivatives!
//   delete[]ddphi_loc;
  return;
}
// =============================================================
/// Shape functions derivatives at the nodal points
void   MGFE::get_dphi_node(
  const int kdim,// dimension <-
  const int node,  // nodal point <-
  const double InvJac[], // Jacobean <-
  double dphi[]   // global derivatives at the nodal point ->
) {// =========================================

  const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
  const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
  const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset

  double gradphi_g[DIMENSION]; // temp grad phi

  for (int eln=0; eln<el_nnodes; eln++)    {
    int lqp=eln* el_ngauss+node;
    for (int idim=0;idim<_dim; idim++) gradphi_g[idim] = _dphidxez_map1_nodes[kdim-1][lqp+idim*goffset];
    for (int idim=0;idim<_dim; idim++) {
      double sum = 0.;
      for (int jdim=0; jdim<_dim; jdim++) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
      dphi[eln+idim* el_nnodes] = sum;
    }
  }
  return;
}
// =============================================================
/// Shape functions derivatives at the gauss point qp
void   MGFE::get_dphi_gl_g(
  const int kdim,// dimension <-
  const int qp,  // gaussian point <
  const double InvJac[], // Jacobean
  double dphi[],   // global derivatives ->
  int sdim
) {// =========================================

  const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
  const int el_ngauss =  _NoGauss1[kdim-1];            // # of guass points
  const int goffset   =  el_nnodes*el_ngauss;  //gauss offset

  double gradphi_g[DIMENSION]; // temp grad phi
  for (int idim=0;idim<3; idim++) gradphi_g[idim] =0.;
  for (int eln=0; eln<el_nnodes; eln++)    {
    int lqp=eln* el_ngauss+qp;
    for (int idim=0;idim<sdim; idim++) gradphi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
    for (int idim=0;idim<sdim; idim++) {
      double sum = 0.;
      for (int jdim=0; jdim<sdim; jdim++) sum += InvJac[jdim+idim*sdim]*gradphi_g[jdim];
      dphi[eln+idim* el_nnodes] = sum;
    }
  }
  return;
}

// =============================================================
/// Shape functions derivatives at the gauss point qp
void   MGFE::get_dphi_gl_g(
  const int kdim,           // dimension <-
  const int qp,             // gaussian point <-
  const double InvJac[], // Jacobean
  std::vector<double> & dphi // global derivatives ->
) {// =========================================

  const int el_nnodes =    _NoShape[kdim-1];           // # of shape functions
  const int el_ngauss   =  _NoGauss1[kdim-1];          // # of guass points
  const int goffset   =    el_nnodes*_NoGauss1[kdim-1];   //gauss offset

  double dphidxi_g[DIMENSION];// temp grad phi

  for (int eln=0; eln<el_nnodes; eln++)    {
    int lqp=eln* el_ngauss+qp;
    for (int idim=0; idim<_dim; idim++) dphidxi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
    for (int idim=0; idim<_dim; idim++) {
      double sum = 0.;
      for (int jdim=0; jdim<_dim; jdim++) sum += InvJac[jdim+idim*_dim]*dphidxi_g[jdim];
      dphi[eln+idim* el_nnodes] = sum;
    }
  }
  return;
}

// ======================================
///This function computes the normal at the gauss point
void MGFE::normal_g(
  const double* xx,   // all surface coordinates <-
  double* normal_g    // normal ->
) const {// ======================================

 // coordinates
  double xx3D[3*NDOF_FEMB];
  for (int i=0; i<_dim*NDOF_FEMB; i++) xx3D[i]=xx[i];

  double tg01[DIMENSION];  // tangent  line
 #if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
  for (int i=0; i<2; i++)  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
  normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
#else // 3D ---------------------------------------------------
  double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
  for (int i=0; i<3; i++)  {
    tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
    tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
    tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
  }
//   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
  normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
  normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];
  
#endif
 // normalization -----------------------
  double mm=0.;
  for (int idim=0; idim< _dim; idim++) mm +=normal_g[idim]*normal_g[idim];
  mm=sqrt(mm);  for (int idim=0; idim< _dim; idim++) normal_g[idim] /=mm;

  return;

}


// ======================================
///This function computes the normal at the gauss point
///   x_c[] is an interior point
void MGFE::normal_g(
  const double* xx,   // all surface coordinates <-
  const double x_c[],   // point inside the element <-
  double* normal_g    // normal ->
) const {// ======================================

 // coordinates
  double xx3D[3*NDOF_FEMB];
  for (int i=0; i<_dim*NDOF_FEMB; i++) xx3D[i]=xx[i];

  double tg01[DIMENSION];  // tangent  line
 #if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
  for (int i=0; i<2; i++)  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
  normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
  // if the sign is not correct then reverse it

  if( normal_g[0]*(x_c[0]-xx3D[0] )+normal_g[1]*(x_c[1]-xx3D[0+NDOF_FEMB]) >0.   ){
    normal_g[0] *=-1.;     normal_g[1] *=-1.;  
//    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  
#else // 3D ---------------------------------------------------
  double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
  for (int i=0; i<3; i++)  {
    tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
    tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
    tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
  }
//   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
  normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
  normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];
  
  // if the sign is not correct then reverse it
  if( normal_g[0]*(x_c[0]-xx3D[0])+
      normal_g[1]*(x_c[1]-xx3D[0+NDOF_FEMB])+ 
      normal_g[2]*(x_c[2]-xx3D[0+2*NDOF_FEMB]) >0   ){
    normal_g[0] *=-1;     normal_g[1] *=-1;  normal_g[2] *=-1;
  }
  
#endif
 // normalization -----------------------
  double mm=0.;
  for (int idim=0; idim< _dim; idim++) mm +=normal_g[idim]*normal_g[idim];
  mm=sqrt(mm);  for (int idim=0; idim< _dim; idim++) normal_g[idim] /=mm;

  return;

}


// ======================================
///This function computes the normal at the gauss point
///   x_c[] is an interior point
void MGFE::normal_g(
  const double* xx,   // all surface coordinates <-
  const double x_c[],   // point inside the element <-
  double* normal_g,    // normal ->
  int & sign_normal
) const {// ======================================

 // coordinates
  double xx3D[3*NDOF_FEMB];
  for (int i=0; i<_dim*NDOF_FEMB; i++) xx3D[i]=xx[i];

  double tg01[DIMENSION];  // tangent  line
 #if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
  for (int i=0; i<2; i++)  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
  normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
  // if the sign is not correct then reverse it
  sign_normal=1;
  if( normal_g[0]*(x_c[0]-xx3D[0] )+normal_g[1]*(x_c[1]-xx3D[0+NDOF_FEMB]) >0.   ){
    normal_g[0] *=-1.;     normal_g[1] *=-1.;  sign_normal=-1;
//    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  
#else // 3D ---------------------------------------------------
  double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
  for (int i=0; i<3; i++)  {
    tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
    tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
    tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
  }
//   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
  normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
  normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];
  
  // if the sign is not correct then reverse it
  if( normal_g[0]*(x_c[0]-xx3D[0])+
      normal_g[1]*(x_c[1]-xx3D[0+NDOF_FEMB])+ 
      normal_g[2]*(x_c[2]-xx3D[0+2*NDOF_FEMB]) >0   ){
    normal_g[0] *=-1;     normal_g[1] *=-1;  normal_g[2] *=-1;
  }
  
#endif
 // normalization -----------------------
  double mm=0.;
  for (int idim=0; idim< _dim; idim++) mm +=normal_g[idim]*normal_g[idim];
  mm=sqrt(mm);  for (int idim=0; idim< _dim; idim++) normal_g[idim] /=mm;

  return;

}




