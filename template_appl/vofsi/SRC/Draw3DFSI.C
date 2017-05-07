#include <iostream>
#include <math.h>


#define  TOLL  1.e-12//tolerance for setting the BCs

using namespace std;

class Draw3DFSI {
private:
    //static double toll=1.e-5; //BDRY_TOLL;

    static bool isNear(double val,double to,double h) {
        if(val>to-h*0.5 && val<to+h*0.5) {
            return true;
        } else {
            return false;
        }
    }


public:

    /////////////////////////////////////////////////////////////
    // User_FSI - internal
    ////////////////////////////////////////////
    static void init (int bc_flag[],int bc_Neum[]) {
        bc_flag[0]=0;
        bc_Neum[0]=1 ;
        bc_flag[1]=0;
        bc_Neum[1]=1 ;
        bc_flag[2]=0;
        bc_Neum[2]=1 ;
    };

    static void internal_rectangle (double xp[],int bc_flag[], int bc_Neum[],double x1, double x2, double y1,double y2,double z1, double z2) {
        double toll=TOLL;
        if(	xp[0] > x1 - toll &&
                xp[0] < x2 + toll  &&
                xp[1] > y1 - toll &&
                xp[1] < y2 + toll &&
                xp[2] > z1 - toll &&
                xp[2] < z2 + toll  )
        {
            bc_flag[0]=0;
            bc_Neum[0]=3 ;
            bc_flag[1]=0;
            bc_Neum[1]=3 ;
            bc_flag[2]=0;
            bc_Neum[2]=3 ;
        }
    };

//     static void rectangle_interface (double xp[],int bc_flag[], int bc_Neum[],double x1, double y1, double z1, double x2, double y2, double z2) {
//       double toll=TOLL;
//       if(	xp[0] > x1 - toll &&
// 		xp[0] < x2 + toll  &&
// 		xp[1] > y1 - toll &&
// 		xp[1] < y2 + toll &&
// 		xp[2] > z1 - toll &&
// 		xp[2] < z2 + toll  )
// 	{
// 	  bc_flag[0]=0;     bc_Neum[0]=3 ;
// 	  bc_flag[1]=0;     bc_Neum[1]=3 ;
// 	}
//     };

    static void plane_interface(double xp[],int bc_flag[], int bc_Neum[],double x1, double x2, double y1,double y2,double z1, double z2) {
        double toll=TOLL;
	      if(xp[0]>x1 && xp[0]< x2 &&
	xp[1]>y1 && xp[1]< y2 &&
	xp[2]>z1 && xp[2]< z2 
      ){bc_flag[0]=0;bc_Neum[0]=5 ;
	bc_flag[1]=0;bc_Neum[1]=5 ;
	bc_flag[2]=0;bc_Neum[2]=5 ;
      } 
	
//         if((x1-x2)<toll) {
//             if(	xp[2] > z1 - toll &&
//                     xp[2] < z2 + toll  &&
//                     xp[1] > y1 - toll &&
//                     xp[1] < y2 + toll )
//             {
//                 bc_flag[0]=0;
//                 bc_Neum[0]=5 ;
//                 bc_flag[1]=0;
//                 bc_Neum[1]=5 ;
//                 bc_flag[2]=0;
//                 bc_Neum[2]=5 ;
//             }
//         } else if((y1-y2)<toll) {
//             if(	xp[2] > z1 - toll &&
//                     xp[2] < z2 + toll  &&
//                     xp[0] > x1 - toll &&
//                     xp[0] < x2 + toll )
//             {
//                 bc_flag[0]=0;
//                 bc_Neum[0]=5 ;
//                 bc_flag[1]=0;
//                 bc_Neum[1]=5 ;
//                 bc_flag[2]=0;
//                 bc_Neum[2]=5 ;
//             }
//         } else if((z1-z2)<toll) {
//             if(	xp[1] > y1 - toll &&
//                     xp[1] < y2 + toll  &&
//                     xp[0] > x1 - toll &&
//                     xp[0] < x2 + toll )
//             {
//                 bc_flag[0]=0;
//                 bc_Neum[0]=5 ;
//                 bc_flag[1]=0;
//                 bc_Neum[1]=5 ;
//                 bc_flag[2]=0;
//                 bc_Neum[2]=5 ;
//             }
//         }


    };



    //////////////////////////////////////////////////
    // User_FSI - bc_read
    ///////////////////////////////////////////////////

    static bool bc_plane(double xp[],int bc_flag[], int bc_Neum[],double x1, double x2, double y1,double y2,double z1, double z2) {
        double toll=TOLL;
        if((x1-x2)<toll) {
            if(	xp[2] > z1 - toll &&
                    xp[2] < z2 + toll  &&
                    xp[1] > y1 - toll &&
                    xp[1] < y2 + toll )
            {
                return true;
            } else {
                return false;
            }
        } else if((y1-y2)<toll) {
            if(	xp[2] > z1 - toll &&
                    xp[2] < z2 + toll  &&
                    xp[0] > x1 - toll &&
                    xp[0] < x2 + toll )
            {
                return true;
            } else {
                return false;
            }
        } else if((z1-z2)<toll) {
            if(	xp[1] > y1 - toll &&
                    xp[1] < y2 + toll  &&
                    xp[0] > x1 - toll &&
                    xp[0] < x2 + toll )
            {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }


    };
    /*
    static bool bc_HLine(double xp[],double y1, double x1, double x2) {
      double toll=TOLL;
      if(	xp[0] > x1 - toll &&
    	xp[0] < x2 + toll  &&
    	xp[1] > y1 - toll &&
    	xp[1] < y1 + toll )
    {
       return true;
    }else{
      return false;
    }
    };

    static bool bc_VLine(double xp[],double x1, double y1, double y2) {
      double toll=TOLL;
      if(	xp[0] > y1 - toll &&
    	xp[0] < y2 + toll  &&
    	xp[1] > x1 - toll &&
    	xp[1] < x1 + toll )
    {
       return true;
    }else{
      return false;
    }
    };*/

    //////////////////////////////////////////////////
    // User_DS - bc_intern_read
    ///////////////////////////////////////////////////

    static void UserDS_init (int bc_flag[],int bc_Neum[])
    {
        bc_flag[0]=0;
        bc_Neum[0]=1 ;
    };

    static void UserDS_internal_rectangle (double xp[],int bc_flag[], int bc_Neum[],
                                           double x1, double x2,
                                           double y1, double y2,
                                           double z1, double z2) {
        double toll=TOLL;
        if(	xp[0] > x1 - toll &&
                xp[0] < x2 + toll  &&
                xp[1] > y1 - toll &&
                xp[1] < y2 + toll &&
                xp[2] > z1 - toll &&
                xp[2] < z2 + toll  )
        {
            bc_flag[0]=0;
            bc_Neum[0]=3;
        }
    };


    static void UserDS_plane_interface(double xp[],int bc_flag[], int bc_Neum[],double x1, double x2, double y1,double y2,double z1, double z2) {
      
      double toll1=1.e-6; 
     
      if(xp[0]>x1 && xp[0]< x2 &&
	xp[1]>y1 && xp[1]< y2 &&
	xp[2]>z1 && xp[2]< z2 
      ){bc_flag[0]=0;bc_Neum[0]=5 ;}

    };





};
