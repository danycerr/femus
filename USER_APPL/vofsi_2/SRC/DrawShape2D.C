#include <iostream>
#include <math.h>


#define MIN_VAL 1.e-12
using namespace std;

class DrawShape2D {
  private:
    static bool isNear(double val,double to,double h){
	if(val>to-h*0.5 && val<to+h*0.5){
	  return true;
	}else{
	  return false;
	}
      }

    
  public:
    
    
    static void rectangle (	///draw a rectangle
			double* temp,int i, int j,
			double hx,double hy,
			double xa, double ya,///bottom left corner
			double xb, double yb///top right corner
			
		       ) {
      double x1 = xa; double x2 = xb;
      double y1 = ya; double y2 = yb;
      double ic = i*hx; double jc = j*hy;
      
      if( jc>y1 && jc<y2 && ic>x1  && ic<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
      if(  jc>y1 && jc<y2 && isNear(ic,x1,hx))  temp[i+1]=0.5;
      if(  jc>y1 && jc<y2 && isNear(ic,x2,hx))  temp[i+1]=0.5;
      if(  ic>x1 && ic<x2 && isNear(jc,y1,hy))  temp[i+1]=0.5;
      if(  ic>x1 && ic<x2 && isNear(jc,y2,hy))  temp[i+1]=0.5;
      if(  isNear(ic,x1,hx) && isNear(jc,y1,hy))  temp[i+1]=0.25;
      if(  isNear(ic,x1,hx) && isNear(jc,y2,hy))  temp[i+1]=0.25;
      if(  isNear(ic,x2,hx) && isNear(jc,y1,hy))  temp[i+1]=0.25;
      if(  isNear(ic,x2,hx) && isNear(jc,y2,hy))  temp[i+1]=0.25;
      
      if( j>y1 && j<y2 && i*hx>x1  && i*hx<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
      if(  j>y1 && j<y2 && isNear(i*hx,x1,hx))  temp[i+1]=0.5;
      if(  j>y1 && j<y2 && i==x2)  temp[i+1]=0.5;
      if(  i*hx>x1 && i<x2 && j==y1)  temp[i+1]=0.5;
      if(  i*hx>x1 && i<x2 && j==y2)  temp[i+1]=0.5;
      if(  i*hx==x1 && j==y1)  temp[i+1]=0.25;
      if(  i*hx==x1 && j==x2)  temp[i+1]=0.25;
      if(  i*hx==x2 && j==y1)  temp[i+1]=0.25;
      if(  i*hx==x2 && j==x2)  temp[i+1]=0.25;
    };
    
        static void rectangle_filled (	///draw a rectangle
			double* temp,int i, int j,
			double hx,double hy,
			double xa, double ya,///bottom left corner
			double xb, double yb///top right corner
			
		       ) {
      double x1 = xa; double x2 = xb;
      double y1 = ya; double y2 = yb;
      double ic = i*hx; double jc = j*hy;
      
      if( jc>y1 && jc<y2 && ic>x1  && ic<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
//       if(  jc>y1 && jc<y2 && isNear(ic,x1,hx))  temp[i+1]=0.5;
//       if(  jc>y1 && jc<y2 && isNear(ic,x2,hx))  temp[i+1]=0.5;
//       if(  ic>x1 && ic<x2 && isNear(jc,y1,hy))  temp[i+1]=0.5;
//       if(  ic>x1 && ic<x2 && isNear(jc,y2,hy))  temp[i+1]=0.5;
//       if(  isNear(ic,x1,hx) && isNear(jc,y1,hy))  temp[i+1]=0.25;
//       if(  isNear(ic,x1,hx) && isNear(jc,y2,hy))  temp[i+1]=0.25;
//       if(  isNear(ic,x2,hx) && isNear(jc,y1,hy))  temp[i+1]=0.25;
//       if(  isNear(ic,x2,hx) && isNear(jc,y2,hy))  temp[i+1]=0.25;
      
//       if( j>y1 && j<y2 && i*hx>x1  && i*hx<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
//       if(  j>y1 && j<y2 && isNear(i*hx,x1,hx))  temp[i+1]=0.5;
//       if(  j>y1 && j<y2 && i==x2)  temp[i+1]=0.5;
//       if(  i*hx>x1 && i<x2 && j==y1)  temp[i+1]=0.5;
//       if(  i*hx>x1 && i<x2 && j==y2)  temp[i+1]=0.5;
//       if(  i*hx==x1 && j==y1)  temp[i+1]=0.25;
//       if(  i*hx==x1 && j==x2)  temp[i+1]=0.25;
//       if(  i*hx==x2 && j==y1)  temp[i+1]=0.25;
//       if(  i*hx==x2 && j==x2)  temp[i+1]=0.25;
    };
    
    static void pool (	///draw a pool at bottom
			double* temp,int i, int j,
			double hx,double hy,
			double xa, double ya,///bottom left corner
			double xb, double yb///top right corner
			
		       ) {
      double x1 = xa; double x2 = xb;
      double y1 = ya; double y2 = yb;
      double ic = i*hx; double jc = j*hy;
      
      if( jc>y1 && jc<y2 && ic>x1  && ic<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
      if(  ic>x1 && ic<x2 && isNear(jc,y2,hy))  temp[i+1]=0.5;
      
//       if( j>y1 && j<y2 && i*hx>x1  && i*hx<x2  )  temp[i+1]=1.;//cc[Level].Cmp[ind]=1 ;
//       if(  j>y1 && j<y2 && isNear(i*hx,x1,hx))  temp[i+1]=0.5;
//       if(  j>y1 && j<y2 && i==x2)  temp[i+1]=0.5;
//       if(  i*hx>x1 && i<x2 && j==y1)  temp[i+1]=0.5;
//       if(  i*hx>x1 && i<x2 && j==y2)  temp[i+1]=0.5;
//       if(  i*hx==x1 && j==y1)  temp[i+1]=0.25;
//       if(  i*hx==x1 && j==x2)  temp[i+1]=0.25;
//       if(  i*hx==x2 && j==y1)  temp[i+1]=0.25;
//       if(  i*hx==x2 && j==x2)  temp[i+1]=0.25;
    };
    
    static void sinwave (	///draw a sin
			double* temp,int i, int j,
			double hx,double hy,
			double w, double phi,
			double scale,  
			double d1, double d2// dominio
			
		       ) {
      double ic = i*hx; double jc = j*hy;
      
      if(ic>0.0+hx/4 && jc>0.0+hy/4 && ic<1.-hx/4 && jc<1.-hy/4){
	      if(ic>d1 && ic <d2){
		if(jc < sin(ic*w+phi)*scale){
		  temp[i+1]=1.;
		  if(jc<0.+hy/2 || ic<0.0+hx/2){
		    temp[i+1]=0.5;
		  }
		}
  // 	      
		if(jc > (sin(ic*w+phi)*scale-hy) && jc < sin(ic*w+phi)*scale ){
		  temp[i+1]=0.5;
		}
		if(jc > (sin(ic*w+phi-0.1)*scale-hy) && jc < sin(ic*w)*scale ){
		  temp[i+1]=0.5;
		}
		if(jc > (sin(ic*w+phi+0.1)*scale-hy) && jc < sin(ic*w)*scale ){
		  temp[i+1]=0.5;
		}	      
	      }
	      
    }
    };
    
    static void circle (	///draw a circle
			double* temp,int i, int j,unsigned int nadd,
			double hx,double hy,
			double x0, double y0,
			double r
		       ) {
      double pxl = i*hx; double pxr = (i+1) *hx;
      double pyb = j*hy; double pyt = (j+1) *hy;
      unsigned int vin = 0; // polygon vertex index
      double ptx[nadd+6]; double pty[nadd+6];
      
      if((pxl-x0) * (pxl-x0) + (pyb-y0) * (pyb-y0) < r*r) vin++;  // bottom-left
      if((pxr-x0) * (pxr-x0) + (pyb-y0) * (pyb-y0) < r*r) vin++;   // bottom-right
      if((pxr-x0) * (pxr-x0) + (pyt-y0) * (pyt-y0) < r*r) vin++;   // top-right
      if((pxl-x0) * (pxl-x0) + (pyt-y0) * (pyt-y0) < r*r) vin++;   // top-left

      // case vin =4 ********************************
      if(vin == 4)  temp[i+1]=1.; //cc.Cmp[ind]=1.;
      // case vin >0 ********************************
      else if(vin > 0) {  //there's some fluid
        vin = 0;
        // checking intersection with cell
        int int1 = -1; int int2 = -1;
        // bottom-left -----------------------------
        if((pxl-x0) * (pxl-x0) + (pyb-y0) * (pyb-y0) < r*r) {
          ptx[vin] = pxl; pty[vin] = pyb; vin++;
        }
        // bottom edge -----------------------------
        double delta = r*r - (pyb - y0) * (pyb - y0);
        if(delta >= - MIN_VAL) {
          delta = (delta < 0.) ? 0. : delta;
          if(x0+sqrt(delta) > pxl && x0+sqrt(delta) <= pxr) {
            ptx[vin] = x0+sqrt(delta); pty[vin] = pyb;
            int1 = vin; vin++;
          }
          if(x0-sqrt(delta) > pxl && x0-sqrt(delta) <= pxr) {
            ptx[vin] = x0-sqrt(delta); pty[vin] = pyb;
            int1 = vin; vin++;
          }
        }
        // bottom-right --------------------------
        if((pxr-x0) * (pxr-x0) + (pyb-y0) * (pyb-y0) < r*r) {
          ptx[vin] = pxr; pty[vin] = pyb; vin++;
        }
        // right edge ---------------------------
        delta = r*r - (pxr - x0) * (pxr - x0);
        if(delta >= - MIN_VAL) {
          delta = (delta < 0.) ? 0. : delta;
          if(y0+sqrt(delta) > pyb && y0+sqrt(delta) <= pyt) {
            ptx[vin] = pxr; pty[vin] = y0+sqrt(delta);
            if(int1==-1) int1=vin; else int2=vin;
            vin++;
          }
          if(y0-sqrt(delta) > pyb && y0-sqrt(delta) <= pyt) {
            ptx[vin] = pxr; pty[vin] = y0-sqrt(delta);
            if(int1==-1) int1=vin; else int2=vin;
            vin++;
          }
        }
        // top-right --------------------------
        if((pxr-x0) * (pxr-x0) + (pyt-y0) * (pyt-y0) < r*r) {
          ptx[vin] = pxr; pty[vin] = pyt; vin++;
        }
        // top edge ----------------------------
        delta = r*r - (pyt - y0) * (pyt - y0);
        if(delta >= - MIN_VAL) {
          delta = (delta < 0.) ? 0. : delta;
          if(x0+sqrt(delta) >= pxl && x0+sqrt(delta) < pxr) {
            ptx[vin] = x0+sqrt(delta); pty[vin] = pyt;
            if(int1==-1) int1=vin; else int2=vin;
            vin++;
          }
          if(x0-sqrt(delta) >= pxl && x0-sqrt(delta) < pxr) {
            ptx[vin] = x0-sqrt(delta); pty[vin] = pyt;
            if(int1==-1) int1=vin; else int2=vin;
            vin++;
          }
        }
        // top-left
        if((pxl-x0) * (pxl-x0) + (pyt-y0) * (pyt-y0) < r*r) {
          ptx[vin] = pxl; pty[vin] = pyt; vin++;
        }
        // left edge
        delta = r*r - (pxl - x0) * (pxl - x0);
        if(delta >= - MIN_VAL) {
          delta = (delta < 0.) ? 0. : delta;
          if(y0+sqrt(delta) >= pyb && y0+sqrt(delta) < pyt) {
            ptx[vin] = pxl; pty[vin] = y0+sqrt(delta);
            int2=vin; vin++;
          }
          if(y0-sqrt(delta) >= pyb && y0-sqrt(delta) < pyt) {
            ptx[vin] = pxl; pty[vin] = y0-sqrt(delta);
            int2=vin; vin++;
          }
        }
        // case vin > 2 ********************************
        if(vin > 2) {
          // intersection angles
          double piadd = (ptx[int1] < x0) ? acos(-1.) : 0.;
          double theta1 = atan((pty[int1]-y0) / (ptx[int1]-x0+MIN_VAL)) + piadd;
          piadd = (ptx[int2] < x0) ? acos(-1.) : 0.;
          double theta2 = atan((pty[int2]-y0) / (ptx[int2]-x0+MIN_VAL)) + piadd;
          if(fabs(theta2-theta1) > acos(-1.)) theta2 -= 2*acos(-1.);

          // nadd points inside the cell
          if(int2-int1 != 1) {
            int1 = vin-1;
            double tmp=theta1; theta1=theta2; theta2=tmp;
          }
          for(int k=vin-1; k>int1; k--) {
            ptx[k+nadd] = ptx[k]; pty[k+nadd] = pty[k];
          }
          for(unsigned int k=0; k<nadd; k++) {
            double x = (k+1.) / (nadd+1.);
            ptx[int1+1+k] = x0+r*cos((1-x) *theta1+x*theta2);
            pty[int1+1+k] = y0+r*sin((1-x) *theta1+x*theta2);
          }
          vin += nadd;

          // close polygon
          ptx[vin] = ptx[0]; pty[vin] = pty[0];

          // fluid in polygon
          double area=0.;
          for(unsigned int n=0; n<vin; n++) {
            area += ptx[n]*pty[(n+1) %vin] - pty[n]*ptx[(n+1) %vin];
          }
          area /= 2.;
          // circle
          temp[i+1]=fabs(area) / (hx*hy);
          //cc.Cmp[ind]=fabs(area)/(hx*hy);
  // square
    
          // std::cout << "\n area " << V_GetCmp(&cc[Level],ind) << " "<< ind <<" "<< std::endl;
//           if(temp[i+1] > 1.+MIN_VAL) {
//             std::cout << "Color function initialization failed in cell (" << i << ","
//                       << j << ")" << std::endl; exit(4);
//           }
        }
       
         
      }
    };
};

