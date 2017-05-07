#include "MGSolverCC.h"

#include "DrawShape2D.C"

/// Color function and interface initialization
void MGSolCC::ic_read(
  double tempc2[],
  unsigned int i, unsigned int j, unsigned int k,
  double hx, double hy,double hz,
  unsigned int nadd
){
  //2D
#if DIMENSION==2
      //bullet impact
      /*if(_name==0){
// 	DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.1,0.5,0.025);
	DrawShape2D::rectangle(_tempc2,i,j,hx,hy,0.06,0.49,0.1,0.51); //bullet
      }
      if(_name==1){
//        DrawShape2D::rectangle(_tempc2,i,j,hx,hy,0.25,0.1,0.26,0.48);
//        DrawShape2D::rectangle(_tempc2,i,j,hx,hy,0.25,0.49,0.26,0.9);
	  DrawShape2D::rectangle(_tempc2,i,j,hx,hy,0.25,0.1,0.26,0.9);
      }*/
	  
	  //cilinder explosion
// 	  if(_name==0){
// 		DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.5,0.5,0.05-hx);
// // 		
// 
// 	  }
// 	  if(_name==1){
// 	    DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.5,0.5,0.06);
// 	    DrawShape2D::circle_hole_raw(_tempc2,i,j,nadd,hx,hy,0.5,0.5,0.05);
// 	  }

	  if(_name==1){
// // 	     DrawShape2D::circle(tempc2,i,j,nadd,hx,hy,xc,yc,r);
	    DrawShape2D::circle(tempc2,i,j,nadd,hx,hy,0.515,0.15,0.14);
// 	    DrawShape2D::circle_hole_raw(_tempc2,i,j,nadd,hx,hy,0.5,0.5,0.05);
	    
// 	    DrawShape2D::rectangle(tempc2,i,j,hx,hy,0.0+hx,0.0+hy,0.41/*-0.0625*/,0.3);
// 	    DrawShape2D::sinwave(_tempc2,i,j,hx,hy,10.,0.,0.7,0.,0.3);
	  }
      
       
//        DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.25,0.3,0.03);
//        DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.75,0.5,0.04);
       
      //box empty
//      DrawShape2D::rectangle(_tempc2,i,j,hx,hy,0.45,0.1,0.55,0.2);
//      DrawShape2D::rectangle_hole(_tempc2,i,j,hx,hy,0.47,0.11,0.53,0.19);
      
      //impact test
      //DrawShape2D::circle(_tempc2,i,j,nadd,hx,hy,0.2,0.5,0.1);
      //DrawShape2D::fracture(_tempc2,i,j,hx,hy,0.1,0.4,0.3,0.6);
  //2D END
#endif

//3D
#if DIMENSION==3
    double xl = (i+.5)*hx;
    double yl = (j+0.5)*hy;
    double zl = (k+0.5)*hz;
	    

    tempc2[i]=0.0;
//     //dam break	///////////////////////////////////////////////////////////////////////// 
//     double x1=6*hx,y1=10*hy,z1=4*4*hz;
//     double x2=0.6,y2=1-10*hy,z2=0.95;
//     if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)){
//       tempc2[i]=1.;
//       if((yl<y1+hy && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)) 	tempc2[i]=.5;
//       if((yl<y2 && yl>y2-hy) && (xl>x1 && xl<x2) && (zl>z1 && zl<z2)) 	tempc2[i]=.5;
//       if((yl<y2 && yl>y1) && (xl>x1 && xl<x1+hx) && (zl>z1 && zl<z2))	tempc2[i]=.5;
//       if((yl<y2 && yl>y1) && (xl>x2-hx && xl<x2) && (zl>z1 && zl<z2))	tempc2[i]=.5;
//       if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z1 && zl<z1+hz))	tempc2[i]=.5;
//       if((yl<y2 && yl>y1) && (xl>x1 && xl<x2) && (zl>z2-hz && zl<z2)) 	tempc2[i]=.5;
//     }
// //     ///////////////////////////////////////////////////////////////////////////////////////  
//     
        //sphere	/////////////////////////////////////////////////////////////////////////
    double xc=0.5;double yc=0.5;double zc=1.;double r=0.3;
//     if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)+(zl-zc)*(zl-zc)<(r+2*hx)*(r+2*hy)*(r+2*hz)){
//       tempc2[i]=0.9;
//      if((xl-xc)*(xl-xc)+(yl-yc)*(yl-yc)+(zl-zc)*(zl-zc)<(r)*(r)*(r)) tempc2[i]=1.;
// 
//     }
    ///////////////////////////////////////////////////////////////////////////////////////  
      
#endif
//3d END
}

/// internal conditions
//NOTE: a interface must be initialized in a phase inlet area
void MGSolCC::internal_read(
    double tempc1[],int ix,  int jy, int kz,
    int indl,
    double hx, double hy,double hz,
    double cc1){
  #if DIMENSION==2
   if((ix*hx-0.2)*(ix*hx-0.2)+(jy*hy-0.6)*(jy*hy-0.5)<0.05*0.05){
	  cc1=1.;
	  tempc1[indl+1]=1.;
    }
  #endif
  
  #if DIMENSION==3

  #endif
}

///set the outlet condition for the VOF phase by setting the phase=0 near the boundary
bool MGSolCC::bc_outlet_read(
  int ix,  int jy,int indl,
  double hx, double hy,
  double cc1,
  int boundary/*, //0 far from boundary, 1 near boundary  
  double tempc1[],
  double tempmx1[],
  double tempmy1[]*/
) {
  bool outlet=false;
  
//2D outlet conditions
#if DIMENSION==2
//       if(boundary==0 &&(	//2nd cell near boundary CC->0.01 (used for interface recostruction)
// // 		ix*hx <0.0+2*hx ||	//Left
// 		ix*hx >1.0-2*hx ||	//right
// // 		jy*hy <0.0+2*hy ||	//bottom
// 		jy*hy >1.0-8*hy 	//top //serve min *8 non so il motivo
// 	      )){
// 		outlet=true;
// 		}
//       if(boundary==1 &&( 	//cell near boundary CC->0.
// // 		ix*hx <0.0+1.*hx ||	//left
// 		ix*hx >1.0-1.*hx ||	//right
// // 		jy*hy <0.0+1.*hy ||	//bottom
// 		jy*hy >1.0-4.*hy 	//top //vedi sopra*/
// 	      )){
// 		outlet=true;
// 		}
//2D END outlet conditions
#endif
 
  return outlet;
}