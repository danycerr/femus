#include "Equations_conf.h"  // for  TWO_PHASE
#include "Domain_conf.h"     // for  DIMENSION
#include <iomanip>
#ifdef TWO_PHASE

// C++ stdlib include files -------------------------------
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

// library include ----------------------------------------
#include "hdf5.h"   //hdf5 format
// config file
#include "vof_config.h"

// class includes -----------------------------------------
#include "MGSolverCC.h"
#include "InterfaceFunctionDD.h"

// ======================================================================================
void MGSolCC::setFieldSource(
  double dt,
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * srcField,

 InterfaceFunctionDD * fct 
  
) {

  // mesh structure --------------------------------------------------------------------
  double  hxyz[3]; hxyz[2]=0;  // space step
  int nxyz[3]= {NX,NY,NZ};     // number of divisions
  int n_cell=1; int nNodes=1;
  for(int idim=0; idim<DIMENSION; idim++) {
    hxyz[idim] = 1./nxyz[idim]; n_cell *=nxyz[idim];
    nNodes *=(nxyz[idim]+1);
  }
  
  // velocity field _uvw
  V_Destr(&_uvw);
  V_Constr(&_uvw,(char *)"uvw",n_cmp*nNodes,Normal,_LPTrue);

  
  
  
  
  // mesh structure
//   int nxyz[3]; nxyz[0]=_nxyz[0][Level]; nxyz[1]=_nxyz[1][Level]; nxyz[2]=1;
//   double hxyz[3]; hxyz[0]=1./(nxyz[0]-1); hxyz[1]=1./(nxyz[1]-1); hxyz[2]=0;
// #if DIMENSION==3
//   nxyz[2]=_nxyz[2][Level]; hxyz[2]=1./(nxyz[2]-1);
// #endif
  int  offset=1;
  for(int idim=0; idim<_dim; idim++) offset *=(nxyz[idim]+1);

  // velocity field _uvw
  V_Destr(&_uvw);
  V_Constr(&_uvw,(char *)"uvw",DIMENSION*offset,Normal,_LPTrue);

  // Dof
  double u_value[3]; int dof_u[3];//double Real v_value =0.;
  for(int idim=0; idim<DIMENSION; idim++) {
    u_value[idim]=0.;  dof_u[idim]=0;
  }

  
  
  
  
  
//    int nNodes =offset;
   int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_vof();
  int * map_med = fct->get_map_med();

  ParaMEDMEM::TypeOfField type = srcField->getTypeOfField();
 const ParaMEDMEM::DataArrayDouble * array = srcField->getArray();

  // filling array
  for(int i_mg=0; i_mg < nNodes; i_mg++) {
 int node_mg   = map_mg[i_mg];  // mg  node
  int node_med  = map_med[i_mg];  // med node
//    std::cout<< i_mg << "vof ("<<  node_mg <<","<<  node_med  <<  " ) ";
//     if(node_med>1 && node_mg< nNodes)
     for(int j= 0; j< n_cmp; j++) {
     double valueu=array->getIJ(node_med,j);

//      std::cout<<"     "<<node_mg;
           V_SetCmp(&_uvw,node_mg+j*nNodes+1,valueu/*(valueu*dt)/hxyz[j]*/);
// 	   std::cout<<   valueu <<" ";
     }
//      std::cout<<  std::endl;
//      else{
//         for(int j= 0; j< n_cmp; j++) {
//        double valueu=0.;
//      V_SetCmp(&_uvw,node_mg+j*nNodes+1,valueu/*(valueu*dt)/hxyz[j]*/);
//      } 
//      }
  }
  array->decrRef();     // delete array

  return;
}
// ======================================================================================
void MGSolCC::setFieldSource_disp(
  double dt,
  int interface_name,
  int n_cmp,
  const std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField,
//   const ParaMEDMEM::MEDCouplingFieldDouble * srcField,
 InterfaceFunctionDD * fct 
  
) {

  // mesh structure --------------------------------------------------------------------
  double  hxyz[3]; hxyz[2]=0;  // space step
  int nxyz[3]= {NX,NY,NZ};     // number of divisions
  int n_cell=1; int nNodes=1;
  for(int idim=0; idim<DIMENSION; idim++) {
    hxyz[idim] = 1./nxyz[idim]; n_cell *=nxyz[idim];
    nNodes *=(nxyz[idim]+1);
  }
  
  // velocity field _uvw
  V_Destr(&_dsvw);
  V_Constr(&_dsvw,(char *)"dsvw",n_cmp*nNodes,Normal,_LPTrue);

  
  
  
  
  // mesh structure
//   int nxyz[3]; nxyz[0]=_nxyz[0][Level]; nxyz[1]=_nxyz[1][Level]; nxyz[2]=1;
//   double hxyz[3]; hxyz[0]=1./(nxyz[0]-1); hxyz[1]=1./(nxyz[1]-1); hxyz[2]=0;
// #if DIMENSION==3
//   nxyz[2]=_nxyz[2][Level]; hxyz[2]=1./(nxyz[2]-1);
// #endif
  int  offset=1;
  for(int idim=0; idim<_dim; idim++) offset *=(nxyz[idim]+1);

  // velocity field _uvw
  V_Destr(&_dsvw);
  V_Constr(&_dsvw,(char *)"dsvw",DIMENSION*offset,Normal,_LPTrue);

  // Dof
  double u_value[3]; int dof_u[3];//double Real v_value =0.;
  for(int idim=0; idim<DIMENSION; idim++) {
    u_value[idim]=0.;  dof_u[idim]=0;
  }

  
  
  
  
  
//    int nNodes =offset;
   int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_vof();
  int * map_med = fct->get_map_med();
// 
  for (int ifield=0;ifield<DIMENSION;ifield++){
  
  ParaMEDMEM::TypeOfField type = srcField[ifield]->getTypeOfField();
 const ParaMEDMEM::DataArrayDouble * array = srcField[ifield]->getArray();

  // filling array
  for(int i_mg=0; i_mg < nNodes; i_mg++) {
 int node_mg   = map_mg[i_mg];  // mg  node
  int node_med  = map_med[i_mg];  // med node
//    std::cout<< i_mg << "vof ("<<  node_mg <<","<<  node_med  <<  " ) ";
//     if(node_med>1 && node_mg< nNodes)
     for(int j= 0; j< n_cmp; j++) {
     double valueu=array->getIJ(node_med,j);

//      std::cout<<"     "<<node_mg;
           V_SetCmp(&_dsvw,node_mg+j*nNodes+ (ifield+j)*nNodes +1,valueu/*(valueu*dt)/hxyz[j]*/);
// 	   std::cout<<   valueu <<" ";
     }
//      std::cout<<  std::endl;
//      else{
//         for(int j= 0; j< n_cmp; j++) {
//        double valueu=0.;
//      V_SetCmp(&_uvw,node_mg+j*nNodes+1,valueu/*(valueu*dt)/hxyz[j]*/);
//      } 
//      }
  }
  array->decrRef();     // delete array
  }//end ifield

  return;
}
// ======================================================================================
#define PI 3.14159265358979
inline double PSI(double xi,double yj)  {
  // single vortex
//   return (1.*sin(PI*xi) *sin(PI*xi) *sin(PI*yj) *sin(PI*yj) /PI);
  // 4 vortices
   return (.25*sin(4.*PI*(xi+0.5))*cos(4*PI*(yj+0.5))/PI);

}

// ======================================================================================
void MGSolCC::GenVel(
  const int Level, ///< Level of the velocity field
  const double dt  ///< time step
) { // ==========================================================================

// mesh structure
  int nxyz[3]; nxyz[0]=_nxyz[0][Level]; nxyz[1]=_nxyz[1][Level]; nxyz[2]=1;
  double hxyz[3]; hxyz[0]=1./(nxyz[0]-1); hxyz[1]=1./(nxyz[1]-1); hxyz[2]=0;
#if DIMENSION==3
  nxyz[2]=_nxyz[2][Level]; hxyz[2]=1./(nxyz[2]-1);
#endif
  int  offset=1;
  for(int idim=0; idim<_dim; idim++) offset *=nxyz[idim];

  // velocity field _uvw
  V_Destr(&_uvw);
  V_Constr(&_uvw,(char *)"uvw",DIMENSION*offset,Normal,_LPTrue);

  // Dof
  double u_value[3]; int dof_u[3];//double Real v_value =0.;
  for(int idim=0; idim<DIMENSION; idim++) {
    u_value[idim]=0.;  dof_u[idim]=0;
  }
//   const unsigned int n_u_dofs = NDOF_FEM; const unsigned int n_p_dofs = NDOF_P;
  for(unsigned int ix=0 ; ix <nxyz[0]; ix++) {
    for(unsigned int iy=0 ; iy <nxyz[1]; iy++) {
      for(unsigned int iz=0 ; iz <nxyz[2]; iz++) {
        double xi =ix*hxyz[0];
        double yi =iy*hxyz[1];
        double zi =iz*hxyz[2];
        
        // *************************************************************************
//rotation --------------------------
//       u_value[0]=/*/*3.14159265358979*/*/*(yj-0.5);
//       u_value[1]=-3.14159265358979*(xi-0.5);
//       u_value[2]=0.;
// vortex 2D ---------------------------
//       u_value[1] =  1.*(PSI(xi+hxyz[0],yi)-PSI(xi,yi))/hxyz[0];
//       u_value[0] = -1.*(PSI(xi,yi+hxyz[1])-PSI(xi,yi))/hxyz[1];
//       u_value[2] = 0;
// // vortex 3D ---------------------------
//         const double Pigreco=3.14159265358979;
//         u_value[0] =10* sin(Pigreco*xi)*sin(Pigreco*xi)*(sin(Pigreco*(yi-.5))-sin(Pigreco*(zi-.5)));
//         u_value[1] =10* sin(Pigreco*yi)*sin(Pigreco*yi)*(sin(Pigreco*(zi-.5))-sin(Pigreco*(xi-.5)));
//         u_value[2] =10* sin(Pigreco*zi)*sin(Pigreco*zi)*(sin(Pigreco*(xi-.5))-sin(Pigreco*(yi-.5)));
// translation -------------------------
     u_value[0] =0.1;u_value[1]=0.;u_value[2]=-0.;
        // *****************************************************************************
        
        // set velocity field
        for(int jdim=0; jdim<_dim; jdim++) {
          dof_u[jdim]=ix+iy*nxyz[0]+iz*nxyz[0]*nxyz[1]+1+jdim*offset;
          V_SetCmp(&_uvw,dof_u[jdim],u_value[jdim]);
        }
      }
    }
  }
  return;
}
// -----
// // ------------------------------
// //  Print function
// // ------------------------------





// ======================================================================================
//       fine cc from  cc1[Level]
// ======================================================================================

/// This function prints the time xml file (to run the single time step file ccf.#.xmf)
void  MGSolCC::print_time_fine_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,    ///<  print every print_step
    const int ndigits       ///< number of digit (namefile)
) {// ===================================================================================

  // file name (time_fine.xmf) ----------------------------------------------------------
  char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
  std::stringstream file_name_ccf;
  file_name_ccf << "RESU/time_fine_"<<  _name <<".xmf";
  sprintf(outf2d,file_name_ccf.str().c_str());  std::ofstream out(outf2d);

  //  file xmf text  template -----------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM "  <<  "\"" <<  "../DATA" << "\"" << "[]>\n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  out << "<Grid Name=\" ccf \"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
  // time loop for grid sequence +++++++++++++++++++++++++++++++++++++++++++++++
  for(int it=t_init; it<t_init+n_time_step; it ++)
    if(it%print_step ==0)   {
      out << "<xi:include href=\""<< "ccf_"<<_name<<"."
       << std::setw(ndigits) << std::setfill('0') 
      << it   <<  ".xmf" << "\""
          << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\" >\n";
      out << "<xi:fallback />\n";
      out << " </xi:include>\n";
    } // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  out << "</Grid> \n"; // Grid Collection end
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  // close and clear --------------------------------------------------------------------
  out.close();
  delete [] Named; delete []outf2d;
  return;
}

// ======================================================================================
/// This function prints the fine color file xmf for fine cc (c1[Level])
/// only for cells 0 < C < 1.
/// It is called by the print fine function (print_fine_hdf5())
void MGSolCC::print_fine_xmf(
  double time,
  std::string flag,                    ///<    print flag (time)
  int n_cell_fine,             ///< number of fine cells with 0< C<1
  const unsigned int Level     ///< fine Level
) {// ===================================================================================
  std::ostringstream namefile; namefile<<"./RESU/ccf_"<<_name<< "."<<flag<<".xmf";
  // file name (ccf(flag).xmf) ----------------------------------------------------------
  std::ofstream out(namefile.str().c_str());
  
  // cell topology (to define cell coordinates) -----------------------------------------
  char *Name_type; Name_type=new char[30]; sprintf(Name_type,"%s","Quadrilateral");
  int npt_elem=4; 
  if(_dim==3) {sprintf(Name_type,"%s","Hexahedron"); npt_elem=8;}

  // xmf file template ------------------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<Xdmf>  \n";
  out << "<Domain> \n";
  out << "<Grid Name=\"Mesh\"> \n";
  out << "<Time Value =\""<< time <<"\" />  \n";
  out << "<Topology Type=\""<<Name_type<<"\"  Dimensions=\""<<n_cell_fine<<"\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\""
      << n_cell_fine<< " "<< npt_elem<<"\" Format=\"HDF\"> \n";
  out << "ccf_"<<_name<<"."<< flag <<".h5:conn \n";
  out << "</DataStructure>  \n";
  out << "</Topology> \n";
  out << "<Geometry Type=\"X_Y_Z\">  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_cell_fine*npt_elem<< "  1\" Format=\"HDF\">  \n";
  out << "ccf_"<<_name<<"."<< flag <<".h5:X1 \n";
  out << "</DataStructure> ";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_cell_fine*npt_elem<< " 1\" Format=\"HDF\"> \n" ;
  out << "ccf_"<<_name<<"."<< flag <<".h5:X2 \n";
  out << "</DataStructure>  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_cell_fine*npt_elem<< " 1\" Format=\"HDF\">  \n";
  out << "ccf_"<<_name<<"."<< flag <<".h5:X3 \n";
  out << "</DataStructure>  \n";
  out << " </Geometry> \n";
  out << " <Attribute Name=\"C\" AttributeType=\"Scalar\" Center=\"Cell\" \n>";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_cell_fine<< "  1\" Format=\"HDF\"> \n ";
  out << "ccf_"<<_name<<"."<< flag <<".h5:CCf \n";
  out << "</DataItem> \n";
  out << "</Attribute> \n";
  out << "</Grid> \n";
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();
  return;
}


// ======================================================================================
/// This function prints the fine color function and the file xmf. Only for cells 0< C<1.
void MGSolCC::print_fine_hdf5(
  double time,
  std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level
) { // ==================================================================================

  // storage in hf5 file (Xdmf) ---------------------------------------------------------
//   std::string namefile="./RESU/ccf_"<<_name<<"."+ flag+".h5";
  
  std::ostringstream namefile;  namefile  << "./RESU/ccf_"<<_name<<"." <<flag <<".h5";
//   std::cout << namefile.str() << std::endl;
  hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  // element ----------------------------------------------------------------------------
  int npt_element=4;  int ne_z =1;                            // points for element
  const int ne_x =_nxyz[0][Level]-1; const int ne_y =_nxyz[1][Level]-1; // number of elements
  if(_dim==3) {
    npt_element=8;
#if DIMENSION==3
    ne_z =_nxyz[2][Level]-1;
#endif
  }          // 3D
  double hx=1./ne_x; double hy=1./ne_y; double hz=1./ne_z;    // cell dim (hx,hy)
  // number of fine cell with 0< C <1
  int n_cell=0; for(int ik=1; ik<=ne_y*ne_z; ik++) n_cell += M__GetLen(&c1[Level],ik);
//    for(int ik=1; ik<=ne_y; ik++) n_cell += M__GetLen(&c1[Level],ik);

  print_fine_xmf(time,flag,n_cell,Level); // print the corresponding xmf file

  // coordinates ------------------------------------------------------------------------
  double * xxx=new double[n_cell*npt_element];
  double * yyy=new double[n_cell*npt_element];
  double * zzz=new double[n_cell*npt_element];

  if(_dim==2) { // coord in 2D ----------------------------------------------------------
    int icount=0; // getting the coordinate
    for(unsigned int iy=0; iy<ne_y; iy++) {
      int len=M_GetLen(&c1[Level],iy+1);
      for(unsigned int kx=0; kx<len; kx++) {
        int ix=M__GetPos(&c1[Level],iy+1,kx)-1;
        double  ix0= hx*ix; double  iy0= hy*iy;
        xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
        zzz[icount*npt_element+0]= 0.;
        xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
        zzz[icount*npt_element+1]= 0.;
        xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
        zzz[icount*npt_element+2]= 0.;
        xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
        zzz[icount*npt_element+3]= 0.;
        icount++;
      }
    }
  }

  //3d
  if(_dim==3) { // coord in 2D ----------------------------------------------------------
    int icount=0; // getting the coordinate
    for(unsigned int kz=0; kz<ne_z; kz++)
      for(unsigned int iy=0; iy<ne_y; iy++) {
        int len=M_GetLen(&c1[Level],iy+kz*ne_y+1);
        for(unsigned int kx=0; kx<len; kx++) {
          int ix=M__GetPos(&c1[Level],iy+kz*ne_y+1,kx)-1;
          double  ix0= hx*ix; double  iy0= hy*iy;
          double  iz0= hz*kz;

          xxx[icount*npt_element+0]= ix0;    yyy[icount*npt_element+0]= iy0;
          zzz[icount*npt_element+0]= iz0;
          xxx[icount*npt_element+1]= ix0+hx; yyy[icount*npt_element+1]= iy0;
          zzz[icount*npt_element+1]= iz0;
          xxx[icount*npt_element+2]= ix0+hx; yyy[icount*npt_element+2]= iy0+hy;
          zzz[icount*npt_element+2]= iz0;
          xxx[icount*npt_element+3]= ix0;    yyy[icount*npt_element+3]= iy0+hy;
          zzz[icount*npt_element+3]= iz0;

          xxx[icount*npt_element+4]= ix0;    yyy[icount*npt_element+4]= iy0;
          zzz[icount*npt_element+4]= iz0+hz;
          xxx[icount*npt_element+5]= ix0+hx; yyy[icount*npt_element+5]= iy0;
          zzz[icount*npt_element+5]= iz0+hz;
          xxx[icount*npt_element+6]= ix0+hx; yyy[icount*npt_element+6]= iy0+hy;
          zzz[icount*npt_element+6]= iz0+hz;
          xxx[icount*npt_element+7]= ix0;    yyy[icount*npt_element+7]= iy0+hy;
          zzz[icount*npt_element+7]= iz0+hz;

          icount++;

//         out << ix0 << " " << iy0 << " " << iz0 << " ";
//         out << ix0+hx << " " << iy0 << " " << iz0 << " ";
//         out << ix0+hx << " " << iy0+hy << " " << iz0 << " ";
//         out << ix0 << " " << iy0+hy << " " << iz0 << " ";
//         out << ix0 << " " << iy0 << " " << iz0+hz << " ";
//         out << ix0+hx << " " << iy0 << " " << iz0+hz << " ";
//         out << ix0+hx << " " << iy0+hy << " " << iz0+hz << " ";
//         out << ix0 << " " << iy0+hy << " " << iz0+hz << std::endl;
        }
      }
  }





//   for(unsigned int iy=0; iy<ne_y; iy++) {
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
  // writing _NoNodes -------------------------------------------------------------------
  hsize_t dimsf[2]; dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;
  // coord x -> x1 --------------------
  std::ostringstream Name; Name << "X1";          // name hdf5 dir
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace
  hid_t dtset=H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT,H5P_DEFAULT
#endif
                       );                          // dataset -> xxx
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,xxx);//write
  H5Sclose(dtsp); H5Dclose(dtset);        // close dataspace dataset
  delete[]xxx;                            // delete coord  hdf5 print vector
  // coord x -> x2 --------------------
  std::ostringstream Name2; Name2 << "X2";         // name hdf5 dir
  hid_t dtsp2 = H5Screate_simple(2, dimsf, NULL);  // dataspace
  hid_t  dtset2=H5Dcreate(file,Name2.str().c_str(),H5T_NATIVE_DOUBLE,dtsp2,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          ,H5P_DEFAULT,H5P_DEFAULT
#endif
                         );                       // dataset -> yyy
  H5Dwrite(dtset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,yyy);//write
  H5Sclose(dtsp2); H5Dclose(dtset2);           // close dataspace dataset
  delete[]yyy;                                 // delete coord yyy hdf5 print vector
  // coord x -> x3 -------------------
  std::ostringstream Name3; Name3 << "X3";                 // name hdf5 dir
  hid_t dtsp3 = H5Screate_simple(2, dimsf, NULL);          // dataspace
  hid_t  dtset3=H5Dcreate(file,Name3.str().c_str(),H5T_NATIVE_DOUBLE,dtsp3,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          ,H5P_DEFAULT,H5P_DEFAULT
#endif
                         );                                // dataset -> zzz
  H5Dwrite(dtset3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,zzz);//write
  H5Sclose(dtsp3); H5Dclose(dtset3);                       // close dataspace dataset
  delete[]zzz;                                   // delete coord zzz hdf5 print vector

  // connectivity -----------------------------------------------------------------------
  int *cc_conn=new int[n_cell*npt_element];        // connectivity hdf5 print vector
  for(int iy=0; iy<n_cell*npt_element; iy++) cc_conn[iy]=iy;
  dimsf[0] =n_cell*npt_element;  dimsf[1] = 1;     // dimension hdf5 print vector
  std::ostringstream Nameconn; Nameconn << "conn"; // name hdf5 dir
  hid_t dtspconn = H5Screate_simple(2, dimsf, NULL);// dataspace
  hid_t  dtsetconn=H5Dcreate(file,Nameconn.str().c_str(),H5T_NATIVE_INT,dtspconn,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                             ,H5P_DEFAULT,H5P_DEFAULT
#endif
                            );                 // dataset -> conn
  H5Dwrite(dtsetconn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_conn);// write
  H5Sclose(dtspconn); H5Dclose(dtsetconn);    // close dataspace dataset
  delete[]cc_conn;                            // delete connectivity hdf5 print vector

  // color function --------------------------------------------------------------------
  double *cc_tmp=new double[n_cell];
  int icount=0; // getting the color function from c1[Level] ---------

  // color function
  for(unsigned int kz=0; kz<ne_z; kz++)
    for(unsigned int iy=0; iy<ne_y; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ne_y+1);
      for(unsigned int kx=0; kx<len; kx++) {
        double val=M__GetVal(&c1[Level],iy+kz*ne_y+1,kx);// if (val>1.) val=1.;
        cc_tmp[icount]= val; icount++;
      }
    }



//   for(unsigned int iy=0; iy<ne_y; iy++) {
//     int len=M_GetLen(&c1[Level],iy+1);
//     for(unsigned int kx=0; kx<len; kx++) {
//       double val=M__GetVal(&c1[Level],iy+1,kx); if(val>1.) val=1.;
//       cc_tmp[icount]= val; icount++;
//     }
//   }
  // hdf5 storage ----------------------------------------------------
  dimsf[0] =n_cell;  dimsf[1] = 1;                 // dimension hdf5 print vector
  std::ostringstream Nameccf; Nameccf << "CCf";     // name hdf5 dir
  hid_t dtspccf = H5Screate_simple(2, dimsf, NULL); // dataspace
  hid_t  dtsetccf=H5Dcreate(file,Nameccf.str().c_str(),H5T_NATIVE_DOUBLE,dtspccf,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            ,H5P_DEFAULT,H5P_DEFAULT
#endif
                           );                       // dataset -> CCf
  H5Dwrite(dtsetccf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_tmp);// write
  H5Sclose(dtspccf); H5Dclose(dtsetccf);            // close dataspace dataset
  delete[]cc_tmp;                                   // delete connectivity hdf5 print vector


#ifdef PRINTNORMAL
  // normal --------------------------------------------------------------------
  double *mx_tmp=new double[n_cell];
  double *my_tmp=new double[n_cell];
  double *mz_tmp=new double[n_cell];
  int icount=0; // getting the normal from _mx(y)(z)[Level] --------
  double mx,my,mz;
  for(unsigned int kz=0; kz<nz; kz++)
    for(unsigned int iy=0; iy<ny; iy++) {
      int len=M_GetLen(&c1[Level],iy+kz*ny+1);
      for(unsigned int kx=0; kx<len; kx++) {
        double val=M__GetVal(&c1[Level],iy+kz*ny+1,kx);
        mx=0.; my=0.; mz=0.;
        if(val<1.) {
          mx_tmp[icount]=M__GetVal(&_mx1[Level],iy+kz*ny+1,kx);
          mx_tmp[icount]=M__GetVal(&_my1[Level],iy+kz*ny+1,kx);
          mz=M__GetVal(&_mz1[Level],iy+kz*ny+1,kx);
        }
        out << mx << " " << my << " "<< mz << " ";
      }
    }

#endif

  // close file -------------------------------------------------------------------------
  H5Fclose(file);


  return;
}




// ======================================================================================
//       Point color cc from  cc (vector nx_pt*ny_pt)
// ======================================================================================

// ======================================================================================
/// This function prints the time xml file (to run the single time step file ccf.#.xmf)
void  MGSolCC::print_time_cc_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,   ///<  print every print_step
  const int ndigits       ///< number of digit (namefile)
) {// ===================================================================================

  // file name (time_fine.xmf) ----------------------------------------------------------
  char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
    std::ostringstream file_name; file_name << "RESU/time_cc_"<<_name<<".xmf";
  sprintf(outf2d,file_name.str().c_str());  std::ofstream out(outf2d);

  //  file xmf text  template -----------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM "  <<  "\"" <<  "../DATA" << "\"" << "[]>\n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  out << "<Grid Name=\" cc \"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
  // time loop for grid sequence +++++++++++++++++++++++++++++++++++++++++++++++
  for(int it=t_init; it<t_init+n_time_step; it++)
    if(it%print_step ==0)   {
              out << "<xi:include href=\""
            << "cc_" <<_name<< "."
            << std::setw(ndigits) << std::setfill('0') << it <<  ".xmf" << "\""
            << " xpointer=\"xpointer(//Xdmf/Domain/Grid[1])\" >\n";
        out << "<xi:fallback />\n";
        out << " </xi:include>\n";
      
      
    } // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  out << "</Grid> \n"; // Grid Collection end
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  // close and clear --------------------------------------------------------------------
  out.close();
  delete [] Named; delete []outf2d;
  return;
}

// ======================================================================================
/// This function prints the cc color file xmf for  cc vector
/// It is called by the print fine function (print_fine_hdf5())
void MGSolCC::print_cc_xmf(
  double time,
  std::string &flag,                    ///<    print flag (time)
  const unsigned int Level     ///< fine Level
) {// ===================================================================================

  // file name (ccf(flag).xmf) ----------------------------------------------------------
    
    std::ostringstream namefile; namefile<<"./RESU/cc_"<<_name<<"."<<flag<<".xmf"; std::ofstream out(namefile.str().c_str());

  // mesh initialization  ---------------------------------------------------------------
  char *Named_top; Named_top=new char[30]; sprintf(Named_top,"Quadrilateral"); // name
  int n_cell=1;  int n_pts=1; int npt_element=4; /* double hxyz[DIMENSION];*/
  if(_dim==3) {npt_element=8;  sprintf(Named_top,"Hexahedron");/*hxyz[2]=0.;*/}
  for(int idim=0;idim<_dim;idim++) {
//     _nxyz[idim][Level]=npt_xyz[idim]; // storage class
    n_pts *= _nxyz[idim][Level];    n_cell *= (_nxyz[idim][Level]-1);
  }
//    int nyz=1; for(int idim=1;idim<_dim;idim++) nyz *= npt_xyz[idim];
//   
//   const int ne_x=_nxyz[0][Level]-1; const int ne_y=_nxyz[1][Level]-1; int ne_z =1; // n elements
//   int n_cell=ne_x*ne_y;  int n_pts=(ne_x+1)*(ne_y+1); // 2D n elements and points
//   double hx=1./ne_x; double hy=1./ne_y; ;  double hz=0.;    // cell dim (hx,hy)
//   if(_dim==3) { // 3D
//     npt_element=8; sprintf(Named_top,"Hexahedron"); // 3D element name
// #if DIMENSION==3
//     ne_z =_nxyz[2][Level]-1;
// #endif
//     n_cell=ne_x*ne_y*ne_z; n_pts=(ne_x+1)*(ne_y+1)*(ne_z+1); // 3D n element
//     hz=1./ne_z;
//   }

  // xmf file template ------------------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<Xdmf>  \n";
  out << "<Domain> \n";
  out << "<Grid Name=\"Mesh\"> \n";
  out << "<Time Value =\""<< time <<"\" />  \n";
  out << "<Topology Type=\""<<Named_top<<"\"  Dimensions=\""<<n_cell<<"\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\""
      << n_cell<< " "<< npt_element<<"\" Format=\"HDF\"> \n";
  out << "mesh_0.h5:MSHCONN \n";//<< flag <<".h5:conn \n";
  out << "</DataStructure>  \n";
  out << "</Topology> \n";
  out << "<Geometry Type=\"X_Y_Z\">  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_pts<< "  1\" Format=\"HDF\">  \n";
  out << "mesh_0.h5:/NODES/COORD/X1 \n";//cc."<< flag <<".h5:X1 \n";
  out << "</DataStructure> ";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_pts<< " 1\" Format=\"HDF\"> \n" ;
  out << "mesh_0.h5:/NODES/COORD/X2 \n";//"cc."<< flag <<".h5:X2 \n";
  out << "</DataStructure>  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_pts<< " 1\" Format=\"HDF\">  \n";
  out << "mesh_0.h5:/NODES/COORD/X3 \n";//"cc."<< flag <<".h5:X3 \n";
  out << "</DataStructure>  \n";
  out << " </Geometry> \n";
  out << " <Attribute Name=\"C\" AttributeType=\"Scalar\" Center=\"Cell\" \n>";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      <<n_cell<< "  1\" Format=\"HDF\"> \n ";
  out << "cc_"<<_name<<"."<< flag <<".h5:CC \n";
  out << "</DataItem> \n";
  out << "</Attribute> \n";
  out << "</Grid> \n";
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  // end xmf file template --------------------------------------------------------------
  out.close();
  return;
}

// ======================================================================================
/// This function prints the fine color function and the file xmf. Only for cells 0< C<1.
void MGSolCC::print_cc_hdf5(
  double time,
  std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level
) { // ==================================================================================

  // storage in hf5 file (Xdmf) ---------------------------------------------------------
  
  std::ostringstream namefile; namefile<<"./RESU/cc_"<<_name<<"."<< flag<<".h5"; std::cout << namefile << std::endl;
  hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  
 // mesh initialization  ---------------------------------------------------------------
  int n_cell=1;  int ne_xyz[3]; ne_xyz[2]=1;
  for(int idim=0;idim<_dim;idim++) {  
    n_cell *= (_nxyz[idim][Level]-1);ne_xyz[idim]= _nxyz[idim][Level]-1; 
  }

  // print the corresponding xmf file
  print_cc_xmf(time,flag,Level); 

  double *cc_tmp=new double[n_cell];
  int icount=0; // getting the color function from c1[Level] ---------
    for(unsigned int iz=0; iz<ne_xyz[2]; iz++)  {
   for(unsigned int iy=0; iy<ne_xyz[1]; iy++)  {
      for(unsigned int ix=0; ix<ne_xyz[0]; ix++)  {
//       int kn=mgmesh._elem_map[level][0+el*4]; // 0-node ->(-1,-1) in quad4
//       double xp=mgmesh._xyz[kn]; double yp=mgmesh._xyz[kn+goffset]; // 0-node coord
//       int  ix= (int)(xp/hx+0.5); int  iy= (int)(yp/hx+0.5); // central coordinates
       int ind=ix+iy*(ne_xyz[0]+1)+iz*(ne_xyz[0]+1)*(ne_xyz[1]+1); // cartesian index
       int ind1=ix+iy*(ne_xyz[0])+iz*(ne_xyz[0])*(ne_xyz[1]); // cartesian index
//       out <<   cc.Cmp[ind+1] << " ";
    cc_tmp[ind1]=  cc.Cmp[ind+1];
  }
    }
    }

  // hdf5 storage ----------------------------------------------------
  hsize_t dimsf[2]; dimsf[0] =n_cell;  dimsf[1] = 1;                 // dimension hdf5 print vector
  std::ostringstream Nameccf; Nameccf << "CC";     // name hdf5 dir
  hid_t dtspccf = H5Screate_simple(2, dimsf, NULL); // dataspace
  hid_t  dtsetccf=H5Dcreate(file,Nameccf.str().c_str(),H5T_NATIVE_DOUBLE,dtspccf,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            ,H5P_DEFAULT,H5P_DEFAULT
#endif
                           );                       // dataset -> CCf
  H5Dwrite(dtsetccf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_tmp);// write
  H5Sclose(dtspccf); H5Dclose(dtsetccf);            // close dataspace dataset
  delete[]cc_tmp;                                   // delete connectivity hdf5 print vector

  // close file -------------------------------------------------------------------------
  H5Fclose(file);


  return;
}

// ======================================================================================
// ======================================================================================
//   interface print only 2D
// ======================================================================================
#if DIMENSION==2
// ======================================================================================

// ======================================================================================
/// This function reads the fine color function
void MGSolCC::read_fine(
  const unsigned int flag_print, ///<  file number to print 
  const unsigned int Level       ///<  level
) { // ==================================================================================

  const unsigned int ny =_nxyz[1][Level]-1;
  char *buf; buf=new char[50];
  char *Named; Named=new char[30]; char *inf2d; inf2d=new char[30];
  sprintf(Named,"%d",flag_print);
  // velocity and pressure
  sprintf(inf2d,"./output/cf.%s.vtu",Named);
  std::ifstream in(inf2d);
  if(!in) {printf("input file cf.%d.vtu not found\n",flag_print); exit(3);}
  while(strncmp(buf,"RestartData",11) != 0) in >> buf;
  int len,pos; double value;
  for(unsigned int iy=0; iy<ny; iy++) {
    in >> len; M_SetLen(&c1[Level],iy+1,len);
    for(unsigned int i=0; i<len; i++) {
      in >> pos;
      in >> value;
      M_SetEntry(&c1[Level],iy+1,i,pos,value);
    }
  }
  in.close();
  // Restriction --------------------------------
  OldSol_update(Level,c1[Level],c1_old[Level]);
  GenNorm(Level,&MGSolCC::rec_elv1);
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


// ======================================================================================
//      read cc mesh level from  cc (vector) hdf5 format
// ======================================================================================

// ------------------------------------------
void MGSolCC::read_fine_hdf5(const unsigned int flag_print,const unsigned int Level) {

  const int ny =_nxyz[1][Level]; const int nx =_nxyz[1][Level];
  const double hx =1./(nx-1); const double hy =1./(ny-1);
  int n_pt_elem=4; int *len=new int[ny];
//   char *buf; buf=new char[50];
//   char *Named; Named=new char[30]; char *inf2d; inf2d=new char[30];
//   sprintf(Named,"%d",flag_print);
//   // velocity and pressure
//   sprintf(inf2d,"./output/cf.%s.vtu",Named);
//   std::ifstream in(inf2d);
//   if(!in) {printf("input file cf.%d.vtu not found\n",flag_print); exit(3);}


  // Open an existing file. ---------------
  std::ostringstream meshname;  meshname << "RESU/ccf_"<<_name<<"." << flag_print << ".h5";
  std::cout << " Reading cc from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  // Reading  cc ---------------------------------------

  // Getting dataset
  std::ostringstream Name; Name << "CCf";
  hid_t dtset = H5Dopen(file_id,Name.str().c_str()
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hsize_t dims[2];
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) std::cerr << "MGSolCC::read_fine_hdf5";
  int n_cell=dims[0];  double *cc=new double[n_cell];
  int *ix=new int[n_cell]; int *iy=new int[n_cell];
  status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&cc[0]); // reading
//   H5Fclose(filespace);
//   H5Dclose(dtset);



  // Reading  X ---------------------------------------

  // Getting dataset
  std::ostringstream Name1; Name1 << "X1";
  dtset = H5Dopen(file_id,Name1.str().c_str()
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT
#endif
                 );
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */

  hid_t status2  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status2 ==0) std::cerr << "MGSolCC::read_fine_hdf5";
  int n_coord=n_cell*n_pt_elem;  double *xxx=new double[n_coord];
  status2=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&xxx[0]); // reading
//   H5Fclose(filespace);
//   H5Dclose(dtset);

  // Reading  y ---------------------------------------

  // Getting dataset
  std::ostringstream Name2; Name2 << "X2";
  dtset = H5Dopen(file_id,Name2.str().c_str()
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT
#endif
                 );
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */

  hid_t status3  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status3 ==0) std::cerr << "MGSolCC::read_fine_hdf5";
  double *yyy=new double[n_coord];

  status3=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&yyy[0]); // reading
//   H5Fclose(filespace);
//   H5Dclose(dtset);

  H5Fclose(file_id);

  unsigned int rlen2=M_GetRowDim(&(c1[Level]));
  unsigned int clen2=M_GetClmDim(&(c1[Level]));



  // find ix iy for each cell  and setting len vector
  int iy0=yyy[0]/hy;     int ix_count=0; int iy_count=0; int tot_len=0; // initial values
  for(int ii=0; ii<iy0; ii++) {len[ii]=0; iy_count++;}
  for(int ii=0; ii<n_cell; ii++) {
    iy[ii]=(0.25*(yyy[ii*n_pt_elem]+yyy[ii*n_pt_elem+1]+yyy[ii*n_pt_elem+2]+yyy[ii*n_pt_elem+3]))/hy;
    ix_count++;
    if(iy0!=iy[ii]) {len[iy_count]=ix_count; iy0=iy[ii]; tot_len +=len[iy_count]; ix_count=0; iy_count++;}
    ix[ii]=(0.25*(xxx[ii*n_pt_elem]+xxx[ii*n_pt_elem+1]+xxx[ii*n_pt_elem+2]+xxx[ii*n_pt_elem+3]))/hx;
  }
  for(int ii=iy_count+1; ii<ny-1; ii++) {len[ii]=0;}

  if(iy_count >ny) std::cout << " iy_count >ny ";
  ix_count=0; int icount=0;
  for(unsigned int iy=0; iy<ny-1; iy++) {
    int ilen=len[iy];
    M_SetLen(&c1[Level],iy+1,ilen);

    for(unsigned int i=0; i<ilen; i++) {
      int pos=ix[ix_count+i]+1;
      double value=cc[icount];
      M_SetEntry(&c1[Level],iy+1,i,pos,value);
      icount++;
    }
    ix_count +=ilen;
    unsigned int rleni=M_GetRowDim(&(c1[Level]));
    unsigned int cleni=M_GetClmDim(&(c1[Level]));
  }
  unsigned int rlen=M_GetRowDim(&(c1[Level]));
  unsigned int clen=M_GetClmDim(&(c1[Level]));
  delete []len;
  delete []ix; delete []iy;
  delete []xxx; delete []yyy; delete []cc;


  // Restriction --------------------------------
  OldSol_update(Level,c1[Level],c1_old[Level]);
  GenNorm(Level,&MGSolCC::rec_elv1);
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
//    delete []Named; delete []inf2d;
  return;
}

// ======================================================================================
//       interface cc from  cc vector
// ======================================================================================



/// This function prints the time xml file (to run the single time step file ccf.#.xmf)
void  MGSolCC::print_time_interface_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,   ///<  print every print_step
  const int ndigits       ///< number of digit (namefile)
) {// ===================================================================================

  // file name (time_fine.xmf) ----------------------------------------------------------
  char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
  std::ostringstream file_name; file_name << "RESU/time_interface_"<<_name<<".xmf";
  sprintf(outf2d,file_name.str().c_str());  std::ofstream out(outf2d);

  //  file xmf text  template -----------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM "  <<  "\"" <<  "../DATA" << "\"" << "[]>\n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  out << "<Grid Name=\" c_interface \"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
  // time loop for grid sequence +++++++++++++++++++++++++++++++++++++++++++++++
  for(int it=t_init; it<t_init+n_time_step; it++)
    if(it%print_step ==0)   {
      out << "<xi:include href=\""<< "c_interface_"<<_name<<"."
       << std::setw(ndigits) << std::setfill('0') 
      << it <<  ".xmf" << "\""
          << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\" >\n";
      out << "<xi:fallback />\n";
      out << " </xi:include>\n";
    } // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  out << "</Grid> \n"; // Grid Collection end
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  // close and clear --------------------------------------------------------------------
  out.close();
  delete [] Named; delete []outf2d;
  return;
}

// ======================================================================================
/// This function prints the fine color file xmf for fine cc (c1[Level])
/// only for cells 0 < C < 1.
/// It is called by the print fine function (print_fine_hdf5())
void MGSolCC::print_interface_xmf(
  double time,
  std::string flag,                    ///<    print flag (time)
  int n_lines,             ///< number of lines with cells with 0< C<1
  const unsigned int Level     ///< fine Level
) {// ===================================================================================

  // file name (ccf(flag).xmf) ----------------------------------------------------------
//   char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
//   sprintf(Named,"%d",flag); sprintf(outf2d,"./RESU/c_interface.%s.xmf",Named);
//   std::ofstream out(outf2d);
  std::ostringstream namefile; namefile << "./RESU/c_interface_"<<_name<<"."<<flag<<".xmf";
//   std::string namefile="./RESU/c_interface_"<<_name<<"."+flag+".xmf";
  std::ofstream out(namefile.str().c_str());
  
  // cell topology (to define cell coordinates) -----------------------------------------
  int npt_elem=2;  if(_dim==3) npt_elem=4;

  // xmf file template ------------------------------------------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<Xdmf>  \n";
  out << "<Domain> \n";
  out << "<Grid Name=\"Mesh\"> \n";
  out << "<Time Value =\""<< time <<"\" />  \n";
  out << "<Topology Type=\"Polyline\"  Dimensions=\""<<n_lines<<"\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\""
      << n_lines<< " "<< npt_elem<<"\" Format=\"HDF\"> \n";
  out << "c_interface_"<<_name<<"."<< flag <<".h5:conn \n";
  out << "</DataStructure>  \n";
  out << "</Topology> \n";
  out << "<Geometry Type=\"X_Y_Z\">  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_lines*npt_elem<< "  1\" Format=\"HDF\">  \n";
  out << "c_interface_"<<_name<<"."<< flag <<".h5:X1 \n";
  out << "</DataStructure> ";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_lines*npt_elem<< " 1\" Format=\"HDF\"> \n" ;
  out << "c_interface_"<<_name<<"."<< flag <<".h5:X2 \n";
  out << "</DataStructure>  \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_lines*npt_elem<< " 1\" Format=\"HDF\">  \n";
  out << "c_interface_"<<_name<<"."<< flag <<".h5:X3 \n";
  out << "</DataStructure>  \n";
  out << " </Geometry> \n";
//   out << " <Attribute Name=\"C\" AttributeType=\"Scalar\" Center=\"Cell\" \n>";
//   out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
//       << n_lines<< "  1\" Format=\"HDF\"> \n ";
//   out << "c_interface."<< flag <<".h5:CCf \n";
//   out << "</DataItem> \n";
//   out << "</Attribute> \n";
  out << "</Grid> \n";
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();
  return;
}


// ======================================================================================
/// This function prints the fine color function and the file xmf. Only for cells 0< C<1.
void MGSolCC::print_interface_hdf5(
  double time,
  std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level
) { // ==================================================================================

  // storage in hf5 file (Xdmf) ---------------------------------------------------------
   std::ostringstream namefile;  namefile  << "./RESU/c_interface_"<<_name<<"."<< flag <<".h5";
//     std::string namefile="./RESU/c_interface_"<<_name<<"."+ flag+".h5";
  std::cout << namefile << std::endl;
  hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  // element ----------------------------------------------------------------------------
  int npt_element=2;  int ne_z =1;                            // points for element
  const int ne_x =_nxyz[0][Level]-1; const int ne_y =_nxyz[1][Level]-1; // number of elements
  if(_dim==3) { npt_element=4; /*ne_z =_nxyz[2][Level]-1;*/ }          // 3D
  double hx=1./ne_x; double hy=1./ne_y; double hz=1./ne_z;    // cell dim (hx,hy)

  // evaluating no. of interface segments and points
  int n_lines=0;
  double pt1[3],pt2[3];
  for(unsigned int jy=1; jy<ne_y-1; jy++) {
    int len=M__GetLen(&c1[Level],jy+1);
    for(unsigned int kx=0; kx<len; kx++) {
      double c_ind = M__GetVal(&c1[Level],jy+1,kx);
      if(c_ind > 1.e-10 && c_ind < 1.-1.e-10)  n_lines++;
    }
  }

  print_interface_xmf(time,flag,n_lines,Level); // print the corresponding xmf file

  // coordinates ------------------------------------------------------------------------
  double * xxx=new double[n_lines*npt_element];
  double * yyy=new double[n_lines*npt_element];
  double * zzz=new double[n_lines*npt_element];


  if(_dim==2) { // coord in 2D ----------------------------------------------------------
    int icount=0;
    for(unsigned int jy=1; jy<ne_y; jy++) {
      int len=M__GetLen(&c1[Level],jy+1);
      for(unsigned int kx=0; kx<len; kx++) {
        int ix=M__GetPos(&c1[Level],jy+1,kx)-1;
        unsigned int ind = ix+jy*(ne_x+1);
        double c_ind = M__GetVal(&c1[Level],jy+1,kx);
        if(c_ind > 1.e-10 && c_ind < 1.-1.e-10) {
//      icount++;
          // rec_local_elv(Level,ind,&mx,&my,&alpha);
          double mx=M__GetVal(&_mx1[Level],jy+1,kx);
          double my=M__GetVal(&_my1[Level],jy+1,kx);
//     for (unsigned int ix=1; ix<nx-1; ix++) {
//       unsigned int ind = ix+jy* (nx+1);
//       double c_ind = cc.Cmp[ind+1];
//       if (c_ind > 1.e-10 && c_ind < 1.-1.e-10) {
//       //  rec_local_elv(0,ind,&mx,&my,&alpha);
//        mx=_mx.Cmp[ind+1];my=_my.Cmp[ind+1];
          double alpha=get_alpha2(mx,my,0.,c_ind);
          get_pts(mx,my,alpha,pt1,pt2);
//      if(pt1[0]<0.)pt1[0]=0.;if(pt1[1]<0.) pt1[1]=0.;
//      if(pt2[0]<0.)pt2[0]=0.;if(pt2[1]<0.)pt2[1]=0.;
//      if(pt1[0]>1.)pt1[0]=1.;if(pt1[1]>1.)pt1[1]=1.;
//      if(pt2[0]>1.)pt2[0]=1.;if(pt2[1]>1.)pt2[1]=1.;
          // interface reconstruction points
          xxx[icount*2+0] = (ix+pt1[0]) *hx; yyy[icount*2+0] = (jy+pt1[1]) *hy; zzz[icount*2+0] =0.;
          xxx[icount*2+1] = (ix+pt2[0]) *hx; yyy[icount*2+1] = (jy+pt2[1]) *hy; zzz[icount*2+1] =0.;
//         out << pt1[0] << " " << pt1[1] << " " << 0 << " " << pt2[0] << " " << pt2[1] << " " << 0 <<"\n";

          icount++;
        }
      }
    }

  }
// connectivity




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
  // writing _NoNodes -------------------------------------------------------------------
  hsize_t dimsf[2]; dimsf[0] =n_lines*npt_element;  dimsf[1] = 1;
  // coord x -> x1 --------------------
  std::ostringstream Name; Name << "X1";          // name hdf5 dir
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace
  hid_t dtset=H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT,H5P_DEFAULT
#endif
                       );                          // dataset -> xxx
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,xxx);//write
  H5Sclose(dtsp); H5Dclose(dtset);        // close dataspace dataset
  delete[]xxx;                            // delete coord  hdf5 print vector
  // coord x -> x2 --------------------
  std::ostringstream Name2; Name2 << "X2";         // name hdf5 dir
  hid_t dtsp2 = H5Screate_simple(2, dimsf, NULL);  // dataspace
  hid_t  dtset2=H5Dcreate(file,Name2.str().c_str(),H5T_NATIVE_DOUBLE,dtsp2,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          ,H5P_DEFAULT,H5P_DEFAULT
#endif
                         );                       // dataset -> yyy
  H5Dwrite(dtset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,yyy);//write
  H5Sclose(dtsp2); H5Dclose(dtset2);           // close dataspace dataset
  delete[]yyy;                                 // delete coord yyy hdf5 print vector
  // coord x -> x3 -------------------
  std::ostringstream Name3; Name3 << "X3";                 // name hdf5 dir
  hid_t dtsp3 = H5Screate_simple(2, dimsf, NULL);          // dataspace
  hid_t  dtset3=H5Dcreate(file,Name3.str().c_str(),H5T_NATIVE_DOUBLE,dtsp3,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          ,H5P_DEFAULT,H5P_DEFAULT
#endif
                         );                                // dataset -> zzz
  H5Dwrite(dtset3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,zzz);//write
  H5Sclose(dtsp3); H5Dclose(dtset3);                       // close dataspace dataset
  delete[]zzz;                                   // delete coord zzz hdf5 print vector

  // connectivity -----------------------------------------------------------------------
  int *cc_conn=new int[n_lines*npt_element];        // connectivity hdf5 print vector
  for(int iy=0; iy<n_lines*npt_element; iy++) cc_conn[iy]=iy;
  dimsf[0] =n_lines*npt_element;  dimsf[1] = 1;     // dimension hdf5 print vector
  std::ostringstream Nameconn; Nameconn << "conn"; // name hdf5 dir
  hid_t dtspconn = H5Screate_simple(2, dimsf, NULL);// dataspace
  hid_t  dtsetconn=H5Dcreate(file,Nameconn.str().c_str(),H5T_NATIVE_INT,dtspconn,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                             ,H5P_DEFAULT,H5P_DEFAULT
#endif
                            );                 // dataset -> conn
  H5Dwrite(dtsetconn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_conn);// write
  H5Sclose(dtspconn); H5Dclose(dtsetconn);    // close dataspace dataset
  delete[]cc_conn;                            // delete connectivity hdf5 print vector

//   // color function --------------------------------------------------------------------
//   double *cc_tmp=new double[n_cell];
//   int icount=0; // getting the color function from c1[Level] ---------
//   for(unsigned int iy=0; iy<ne_y; iy++) {
//     int len=M_GetLen(&c1[Level],iy+1);
//     for(unsigned int kx=0; kx<len; kx++) {
//       double val=M__GetVal(&c1[Level],iy+1,kx); if(val>1.) val=1.;
//       cc_tmp[icount]= val; icount++;
//     }
//   }
//   // hdf5 storage ----------------------------------------------------
//    dimsf[0] =n_cell;  dimsf[1] = 1;                 // dimension hdf5 print vector
//   std::ostringstream Nameccf; Nameccf << "CCf";     // name hdf5 dir
//   hid_t dtspccf = H5Screate_simple(2, dimsf, NULL); // dataspace
//   hid_t  dtsetccf=H5Dcreate(file,Nameccf.str().c_str(),H5T_NATIVE_DOUBLE,dtspccf,H5P_DEFAULT
// #if HDF5_VERSIONM != 1808
//                             ,H5P_DEFAULT,H5P_DEFAULT
// #endif
//                            );                       // dataset -> CCf
//   H5Dwrite(dtsetccf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,cc_tmp);// write
//   H5Sclose(dtspccf); H5Dclose(dtsetccf);            // close dataspace dataset
//   delete[]cc_tmp;                                   // delete connectivity hdf5 print vector

  // close file -------------------------------------------------------------------------
  H5Fclose(file);


  return;
}
#endif   // =====================  end 2D ==============================================


// #endif
#endif          // ----------------TWO_PHASE ---------------------- 
