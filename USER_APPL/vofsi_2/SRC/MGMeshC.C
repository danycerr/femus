#include "vof_config.h"
#include "MGMeshC.h"
#include "MGFE_conf.h"
#include "hdf5.h"
#include <sstream>
# include <mpi.h>  // This is needed in the constructor 
// -------------------------------------------
MGMeshC::MGMeshC(const unsigned int dim_in,const unsigned int NoLevels_in):
  _dim(dim_in), _NoLevels(NoLevels_in) {
  _xyz= (double *)malloc(3 * sizeof(double));
  _NoNodes= (int *)malloc(_NoLevels * sizeof(int));
  _NoElements= (int *)malloc(_NoLevels * sizeof(int));
  _offset  = (int **)malloc(_NoLevels * sizeof(int *));
  _elem_map= (int **)malloc(_NoLevels * sizeof(int *));
  _NXNYNZ[0]=1; _NXNYNZ[1]=1; _NXNYNZ[2]=1;
  _aLen[0]=0.; _aLen[1]=0.; _aLen[2]=0.;
  _bLen[0]=1.; _bLen[1]=1.; _bLen[2]=1.;
  int i;
  MPI_Comm_rank(MPI_COMM_WORLD, &i);
  std::cout <<  " ===== Mesh: ===== \n communicator: " << MPI_COMM_WORLD << std::endl;
  _iproc= static_cast< int>(i);
  return;
}





// ---------------------------------
MGMeshC::~MGMeshC() {
  clear();
  if(_elem_map != NULL) free(_elem_map);
  if(_offset != NULL) free(_offset);
  if(_NoElements != NULL)  free(_NoElements);
  if(_NoNodes != NULL)  free(_NoNodes);
  if(_xyz != NULL) free(_xyz);
}
// --------------------------------
void MGMeshC::clear() {
  /* release of dynamic variables */
  for(unsigned int i = 0; i < _NoLevels; i++) {
    if(_elem_map[i] != NULL) free(_elem_map[i]);
    if(_offset[i] != NULL) free(_offset[i]);
  }
}
// ------------------------------------------------
void MGMeshC::init(const std::string& name) {
  std::ifstream in(name.c_str()); this->read_c(in);
}
// ------------------------------------------------
void MGMeshC::write(const std::string& name,int level) {
  std::ofstream in(name.c_str()); this->write_c(in,level);
}
// --------------------------------------------------
void MGMeshC::print_nodes(std::ofstream& out,const unsigned int Level,const unsigned int mode) {
  const unsigned int offset=_NoNodes[_NoLevels-1];
  const unsigned int offset2=2*offset;
  const unsigned int n_nodes=_NoNodes[Level];
  if(mode == 0) {
    if(_dim == 2)      for(unsigned int v=0; v<n_nodes; v++)
        out << _xyz[v] << " " << _xyz[v+offset] << " " << 0  << " ";
    else      for(unsigned int v=0; v<n_nodes; v++)
        out << _xyz[v] << " " << _xyz[v+offset] << " " << _xyz[v+offset2]  << " ";
  } else if(mode == 1) {
    for(unsigned int v=0; v<n_nodes; v++) out << _xyz[v] << " ";
    out << '\n';
    for(unsigned int v=0; v<n_nodes; v++) out <<  _xyz[v+offset] << " ";
    out << '\n';
    if(_dim == 2)  for(unsigned int v=0; v<n_nodes; v++)  out << 0  << " ";
    else  for(unsigned int v=0; v<n_nodes; v++)  out <<  _xyz[v+offset2]  << " ";
  }
  out << "\n";
  return;
}
#ifdef DIM2
// -----------------------------------------------------------------------
void MGMeshC::print_conn(std::ofstream& out,const unsigned int Level,const unsigned int mode) {
  if(mode == 9) {	// Quad9 elements
    int conn[4][4];
    conn[0][0]=0; conn[0][1]=4; conn[0][2]=8; conn[0][3]=7;
    conn[1][0]=4; conn[1][1]=1; conn[1][2]=5; conn[1][3]=8;
    conn[2][0]=8; conn[2][1]=5; conn[2][2]=2; conn[2][3]=6;
    conn[3][0]=7; conn[3][1]=8; conn[3][2]=6; conn[3][3]=3;
    for(unsigned int el=0; el<(unsigned int)_NoElements[Level]; el++) {
      for(unsigned int se=0; se<4; se++)      {
#ifdef OUTGMV
        out << "quad 4 \n";
#endif
        for(unsigned int i=0; i<4; i++) {
          out <<  _elem_map[Level][conn[se][i]+el*9]
#ifdef OUTGMV
              +1
#endif
              << " ";
        }
        out << "\n";
      }
    }
  }
  if(mode == 4) {       // Quad4 elements
    for(unsigned int el=0; el<(unsigned int)_NoElements[Level]; el++) {
      for(unsigned int i=0; i<4; i++) {
        out <<  _elem_map[Level][i+el*4]
            << " ";
      }
      out << "\n";
    }
  }
  return;
}

#else
// ---------------------------------------------------------
void MGMeshC::print_conn(std::ofstream& out,const unsigned int Level,const unsigned int mode) {
  if(mode == 27) {
    int conn[8][8];
    conn[0][0]=0; conn[0][1]=8; conn[0][2]=20; conn[0][3]=11; conn[0][4]=12; conn[0][5]=21; conn[0][6]=26; conn[0][7]=24;
    conn[1][0]=8; conn[1][1]=1; conn[1][2]=9; conn[1][3]=20; conn[1][4]=21; conn[1][5]=13; conn[1][6]=22; conn[1][7]=26;
    conn[2][0]=11; conn[2][1]=20; conn[2][2]=10; conn[2][3]=3; conn[2][4]=24; conn[2][5]=26; conn[2][6]=23; conn[2][7]=15;
    conn[3][0]=20; conn[3][1]=9; conn[3][2]=2; conn[3][3]=10; conn[3][4]=26; conn[3][5]=22; conn[3][6]=14; conn[3][7]=23;

    conn[4][0]=12; conn[4][1]=21; conn[4][2]=26; conn[4][3]=24; conn[4][4]=4; conn[4][5]=16; conn[4][6]=25; conn[4][7]=19;
    conn[5][0]=21; conn[5][1]=13; conn[5][2]=22; conn[5][3]=26; conn[5][4]=16; conn[5][5]=5; conn[5][6]=17; conn[5][7]=25;
    conn[6][0]=24; conn[6][1]=26; conn[6][2]=23; conn[6][3]=15; conn[6][4]=19; conn[6][5]=25; conn[6][6]=18; conn[6][7]=7;
    conn[7][0]=26; conn[7][1]=22; conn[7][2]=14; conn[7][3]=23; conn[7][4]=25; conn[7][5]=17; conn[7][6]=6; conn[7][7]=18;
    for(unsigned int el=0; el<(unsigned int)_NoElements[Level]; el++) {
      for(unsigned int se=0; se<8; se++)      {
#ifdef OUTGMV
        out << " hex 8 \n";
#endif
        for(unsigned int i=0; i<8; i++) {
          out <<  _elem_map[Level][conn[se][i]+el*27]
#ifdef OUTGMV
              +1
#endif
              << " ";
        }
        out << "\n";
      }
    }
    out << "\n";
  }
}
#endif
// --------------------------------
void MGMeshC::print(const unsigned int flag_print,const unsigned int Level) {

  // label 00*,0**,***
  char *Named; Named=new char[30]; char *outf2d; outf2d=new char[30];
  if(flag_print < 10) sprintf(Named,"00%d",flag_print);
  else if(flag_print < 100)  sprintf(Named,"0%d",flag_print);
  else sprintf(Named,"%d",flag_print);
  // file
  sprintf(outf2d,"./output/mesh.%s.vtu",Named);  std::ofstream out(outf2d);
  // set up
  unsigned int n_nodes=_NoNodes[Level]; unsigned int n_elements=_NoElements[Level];
  // file head
  out << "<VTKFile type=\"UnstructuredGrid\"  byte_type=\"LittleEndian\">\n";
  out << "<UnstructuredGrid>\n";
  out << "<Piece NumberOfPoints=\""<<n_nodes << "\" NumberOfCells=\"" << n_elements*NDOF_P << "\">\n";
  // write the nodes --------------------------------------
  out << "<Points> \n";
  out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  print_nodes(out,Level,0);
  out << "</DataArray>\n";
  out<<  "</Points>\n";
  // write the connectivity
  out << "<Cells>\n";
  out << "<DataArray type=\"Int32\" Name=\"connectivity\"  format=\"ascii\">\n";
  print_conn(out,Level,NDOF_FEM);
  out << "</DataArray>\n";
  out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int offsets=0;
  for(unsigned int el=0; el<n_elements; el++) {
    for(unsigned int se=0; se<NDOF_P; se++) {offsets +=NDOF_P; out << offsets << " ";}
  }
  out << "\n";
  out << "</DataArray>\n";
  out << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(unsigned int el=0; el<n_elements; el++) {
    for(unsigned int se=0; se<NDOF_P; se++)   out <<  9+3*(DIMENSION-2) << " ";
  }
  out << "\n";
  out << "</DataArray>\n";
  out << "</Cells>\n";
  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";
  out << "</VTKFile>\n";
  out.close();
  return;
}
// ---------------------------------------------------------
// private
// ---------------------------------------------
void MGMeshC::write_c(std::ostream& out,int level) {
}

// --------------------------------------------
void MGMeshC::read_c(std::ifstream& infile) {
  // file -------------------------------------------------
  if(!infile) {
    std::cout << "Input file "<< infile << " not opened."
              << std::endl;    exit(3);
  }
  const int  bufLen = 256; char  buf[bufLen+1];
  int val; unsigned int level;
  // levels -----------------------------------------------
  while(strncmp(buf,"levels",6) != 0)    infile >> buf;
  infile >> level; fprintf(stderr," Reading %d levels with nodes", level);
  if(_NoLevels != level) std::cout << " while _NoLevels is " << _NoLevels << std::endl;
  // max_nodes (decreasing level order)
  for(unsigned int k=0; k<level; k++) {
    infile >> val; _NoNodes[level-k-1]=val; fprintf(stderr," %d(%d) ",val, level-k-1);
  }
  const unsigned int n_nodes=(unsigned int)_NoNodes[_NoLevels-1];
  // coordinates ---------------------------------
  // memory for coordinates
  if(_xyz != NULL) free(_xyz);
  _xyz= (double *)malloc(n_nodes*_dim * sizeof(double));
  while(strncmp(buf,"nodes",5) != 0) {
    infile >> buf;
  }
  infile >> val;// std::cout << " reading  " << val << " nodes " << std::endl;
  if(_NoNodes[_NoLevels-1] !=val) std::cout << val << "  should be"<< _NoNodes[_NoLevels-1] << std::endl;
  // reading coordinates
  for(unsigned int i=0; i<n_nodes; i++) {infile >> _xyz[i];}
  for(unsigned int i=0; i<n_nodes; i++) {infile >> _xyz[i+n_nodes];}
  if(_dim == 3)
    for(unsigned int i=0; i<n_nodes; i++) {infile >> _xyz[i+2*n_nodes];}
  // cells ---------------------------------
  int nel; unsigned  int mem;
  for(unsigned int k=0; k<_NoLevels; k++) {
    while(strncmp(buf,"level",5) != 0) infile >> buf;
    infile >> level  >> buf  >> nel  >> mem;
    _NoElements[level]=nel;
    _offset[level]=(int *)malloc((_NoElements[level]+1) * sizeof(int));
    _elem_map[level]=(int *)malloc(mem *_NoElements[level]* sizeof(int));
    // connectivity
    for(unsigned int i=0; i<mem *_NoElements[level]; i++) {
      infile >> val; _elem_map[level][i]=val;
    }
    //offset (element nodes till to  offset)
    _offset[level][0]=0;
    for(unsigned int i=1; i<=(unsigned int)_NoElements[level]; i++) {
      _offset[level][i] =mem*i;
    }
  }
  infile.close();
  // std::cout<< " mesh read  "  << std::endl;
}
// ======================================================================================
/// This function generates a Cartesian mesh.
void MGMeshC::gen_c(
//   int  NXYZ0[],   // number of points (each dir) (NX,NY,NZ) in vof system
//   double bxbybz[] // dimension  (bx,by,bz) -> (0,bx)x(0,by)x(0,bz)
) { // ==================================================================================
  int nlevel=_NoLevels;// number of levels (0->nlevel-1) from  constructor
  // Construction vector NXYZ_val=(NX_0,NY_0,NZ_0,NX_1,NY_1,NZ_1,...NX_lev,NY_1ev,NZ_lev)
  printf("\n Generating %d level (0,...,nlevel-1) nodes for Cartesian mesh", nlevel);
  int  *NXYZ_val; NXYZ_val=new int[(nlevel)*3];
  NXYZ_val[0]=_NXNYNZ[0]+1; NXYZ_val[1]=_NXNYNZ[1]+1;     // n_points along x and y (level=0)
  NXYZ_val[2]=1; if(_dim == 3) NXYZ_val[2]=_NXNYNZ[2]+1;// n_points along z       (level=0)
  _NoNodes[0]=NXYZ_val[0]*NXYZ_val[1]*NXYZ_val[2];    // n_point totol          (level=0)

  printf("\n %d(%d) - ",_NoNodes[0], 0);
  int fcl=1;
  for(int level=1; level<nlevel; level++) {
    fcl *=2;
    NXYZ_val[level*3+0]= (NXYZ_val[0]-1) *fcl+1; NXYZ_val[level*3+1]=(NXYZ_val[1]-1)*fcl+1;
    NXYZ_val[level*3+2]=1; if(_dim == 3) NXYZ_val[level*3+2]=(NXYZ_val[2]-1) *fcl+1;
    _NoNodes[level]=NXYZ_val[level*3+0]*NXYZ_val[level*3+1]*NXYZ_val[level*3+2];//(level)
    printf(" %d(%d) - ",_NoNodes[level], level);
  }
//   int nx_pt=_NXNYNZ[0]*fcl+1; int ny_pt=_NXNYNZ[1]*fcl+1; int nz_pt=_NXNYNZ[2]*fcl+1;
//   const unsigned int n_nodes=(unsigned int)_NoNodes[nlevel-1];
  int fd=1;
  for(int ilevel=1; ilevel<nlevel; ilevel++) fd *=2;

  int nx_pt=0; int ny_pt=0; int nz_pt=0;
  int ne_x=0; int ne_y=0; int ne_z=0;
  int n_nodes=0;
  // coordinates ---------------------------------
  if(_xyz != NULL) free(_xyz);

  if(_dim == 2) {
    nx_pt=_NXNYNZ[0]*fcl+1;  ny_pt=_NXNYNZ[1]*fcl+1;
    ne_x=_NXNYNZ[0]*fd;  ne_y=_NXNYNZ[1]*fd;
    n_nodes=nx_pt*ny_pt;
    _xyz= (double *)malloc(n_nodes*_dim * sizeof(double));
    double hx=(_bLen[0]-_aLen[0])/ne_x;
    double hy=(_bLen[1]-_aLen[1])/ne_y;
    for(int i=0; i<nx_pt; i++) {
      for(int j=0; j<ny_pt; j++) {
        int indx=i+j*nx_pt;
        _xyz[indx]=i*hx;
        _xyz[indx+n_nodes]=j*hy;
      }
    }
  }
  if(_dim == 3) {
//   } else {
    nx_pt=_NXNYNZ[0]*fcl+1; ny_pt=_NXNYNZ[1]*fcl+1; nz_pt=_NXNYNZ[2]*fcl+1;
    ne_x=_NXNYNZ[0]*fd; ne_y=_NXNYNZ[1]*fd; ne_z=_NXNYNZ[2]*fd;
    n_nodes=nx_pt*ny_pt*nz_pt;
    _xyz= (double *)malloc(n_nodes*_dim * sizeof(double));
    double hx=(_bLen[0]-_aLen[0])/ne_x;
    double hy=(_bLen[1]-_aLen[1])/ne_y;
    double hz=(_bLen[2]-_aLen[2])/ne_z;
    for(int i=0;  i<nx_pt; i++) {
      for(int j=0; j<ny_pt; j++) {
        for(int k=0; k<nz_pt; k++) {
          int indx=i+j*nx_pt+k*nx_pt*ny_pt;
          _xyz[indx]=i*hx;
          _xyz[indx+n_nodes]=j*hy;
          _xyz[indx+2*n_nodes]=k*hz;
        }
      }
    }
  }

//   Elements ---------------------------------------------------------------------------
  printf("\n Generating %d levels for Cartesian elements", nlevel);
  NXYZ_val[0]=_NXNYNZ[0]; NXYZ_val[1]=_NXNYNZ[1];      // n_points along x and y (level=0)
  NXYZ_val[2]=1; if(_dim == 3) NXYZ_val[2]=_NXNYNZ[2]; // n_points along z       (level=0)
  _NoElements[0]=NXYZ_val[0]*NXYZ_val[1]*NXYZ_val[2]; // n_point totol          (level=0)
  printf("\n %d(%d) - ",_NoElements[0], 0);
  fcl=1;
  for(int level=1; level<nlevel; level++) {
    fcl *=2;
    NXYZ_val[level*3+0]= (NXYZ_val[0]) *fcl; NXYZ_val[level*3+1]=(NXYZ_val[1]) *fcl;
    NXYZ_val[level*3+2]=1; if(_dim == 3) NXYZ_val[level*3+2]=(NXYZ_val[2]) *fcl;
    _NoElements[level]=NXYZ_val[level*3+0]*NXYZ_val[level*3+1]*NXYZ_val[level*3+2];//(level)
    printf(" %d(%d) - ",_NoElements[level], level);
  }


  // cells ---------------------------------
  int nel; int np_elem=0;
//   if(_dim==3) np_elem=8;
  for(unsigned int klevel=0; klevel<_NoLevels; klevel++) {



    // connectivity
    if(_dim==2) { // element quad4

      np_elem=4;
      _offset[klevel]=new int[_NoElements[klevel]+1];
      _elem_map[klevel]=(int *)malloc(np_elem *_NoElements[klevel]* sizeof(int));
      for(unsigned int j=0; j<NXYZ_val[klevel*3+1]; j++) {
        for(unsigned int i=0; i<NXYZ_val[klevel*3+0]; i++) {
          int ielem= i+j*NXYZ_val[klevel*3+0];
          // connectivity cartesian cell	  quad4
          _elem_map[klevel][ ielem*np_elem+0]=i+j*(NXYZ_val[klevel*3+0]+1);
          _elem_map[klevel][ ielem*np_elem+1]=(i+1)+j*(NXYZ_val[klevel*3+0]+1);

          _elem_map[klevel][ ielem*np_elem+2]=(i+1)+(j+1)*(NXYZ_val[klevel*3+0]+1);
          _elem_map[klevel][ ielem*np_elem+3]=i+(j+1)*(NXYZ_val[klevel*3+0]+1);
        }
      }
    }
    if(_dim==3) { // element hex8
      np_elem=8;
      _offset[klevel]=(int *)malloc((_NoElements[klevel]+1) * sizeof(int));
      _elem_map[klevel]=(int *)malloc(np_elem *_NoElements[klevel]* sizeof(int));
      for(unsigned int k=0; k<NXYZ_val[klevel*3+2]; k++) {
        for(unsigned int j=0; j<NXYZ_val[klevel*3+1]; j++) {
          for(unsigned int i=0; i<NXYZ_val[klevel*3+0]; i++) {
            int ielem= i+j*NXYZ_val[klevel*3+0]+k*NXYZ_val[klevel*3+0]*NXYZ_val[klevel*3+1];
            // connectivity cartesian cell	  hex8
            _elem_map[klevel][ ielem*np_elem+0]=i+j*(NXYZ_val[klevel*3+0]+1)+
                                                k*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+1]=(i+1)+j*(NXYZ_val[klevel*3+0]+1)+
                                                k*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+2]=(i+1)+(j+1)*(NXYZ_val[klevel*3+0]+1)+
                                                k*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+3]=i+(j+1)*(NXYZ_val[klevel*3+0]+1)+
                                                k*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);

            _elem_map[klevel][ ielem*np_elem+0]=i+j*(NXYZ_val[klevel*3+0]+1)+
                                                (k+1)*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+1]=(i+1)+j*(NXYZ_val[klevel*3+0]+1)+
                                                (k+1)*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+2]=(i+1)+(j+1)*(NXYZ_val[klevel*3+0]+1)+
                                                (k+1)*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);
            _elem_map[klevel][ ielem*np_elem+3]=i+(j+1)*(NXYZ_val[klevel*3+0]+1)+
                                                (k+1)*(NXYZ_val[klevel*3+0]+1)*(NXYZ_val[klevel*3+1]+1);

          }
        }
      }
    }
    //offset (element nodes till to  offset)
    _offset[klevel][0]=0;
    for(int i=1; i<= _NoElements[klevel]; i++)  _offset[klevel][i]=np_elem*i;

  }
//   infile.close();
  std::cout<< "\n MGMeshC::gen_c: end  mesh generation  "  << std::endl;
  return;
}


void MGMeshC::print_hdf5(const unsigned int flag_print,const unsigned int Level) {





};



// ========================================================
/// This function prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM
void MGMeshC::print_hf5(

  const int  Level      // Level <-
)  { // ================================================

  // Level ------------------------------------------------------------------------------
  int fd=1;  for(int ilevel=1; ilevel<Level; ilevel++) fd *=2;


  // storage in hf5 (Xdmf)
  std::ostringstream namefile;  namefile  << "RESU/mesh_"<< Level-1    <<".h5";
  std::cout << namefile.str() << std::endl;
  hid_t file = H5Fcreate(namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);


  // writing _NoNodes -------------------------------------------------------------------
  hsize_t dimsf[2]; dimsf[0] = 1;  dimsf[1] = 1;
  hid_t dtsp =H5Screate_simple(2, dimsf, NULL);
  hid_t dtset=H5Dcreate(file,"/NNODES",H5T_NATIVE_INT,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT,H5P_DEFAULT
#endif
                       );
  H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, _NoNodes);
  H5Sclose(dtsp);  H5Dclose(dtset);

  // coordinates
  int nx=1; int ny=1; int nz=1;
  double hx=0.;   double hy=0.;   double hz=0.;

  if(_dim> 0) {nx=_NXNYNZ[0]*fd+1; hx=(_bLen[0]-_aLen[0])/(nx-1);}
  if(_dim> 1) {
    ny=_NXNYNZ[1]*fd+1; hy=(_bLen[1]-_aLen[1])/(ny-1);
  }
  if(_dim> 2) {nz=_NXNYNZ[2]*fd+1; hz=(_bLen[2]-_aLen[2])/(nz-1);}

  int n_nodes=nx*ny*nz;
  double *xxx=new double[n_nodes];
  double *yyy=new double[n_nodes];
  double *zzz=new double[n_nodes];
  if(_dim==2){
    for(int ix=0; ix<nx; ix++) {
      for(int jy=0; jy<ny; jy++) {
// 	 for(int kz=0; kz<nz; kz++) {
        int indx=ix+jy*nx;
        xxx[indx]=ix*hx;
        yyy[indx]=jy*hy;
        zzz[indx]=0.;//kz*hz;
//       }
      }
    }
  }
  if(_dim==3){
    for(int ix=0; ix<nx; ix++) {
      for(int jy=0; jy<ny; jy++) {
        for(int kz=0; kz<nz; kz++) {
          int indx=ix+jy*nx+kz*nx*ny;
          xxx[indx]=ix*hx;
          yyy[indx]=jy*hy;
          zzz[indx]=kz*hz;
        }
      }
    }
  }
  //  Writing  _xyz -------------------
  dimsf[0] =n_nodes ;  dimsf[1] = 1;
  // Create a group named "/MyGroup" in the file.
  hid_t   group_id = H5Gcreate(file, "NODES", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                               , H5P_DEFAULT, H5P_DEFAULT
#endif
                              );
  hid_t   group_id2 = H5Gcreate(file, "NODES/COORD", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                , H5P_DEFAULT, H5P_DEFAULT
#endif
                               );
  // coord x -> x1
  std::ostringstream Namex; Namex << "NODES/COORD/X1";
  dtsp = H5Screate_simple(2, dimsf, NULL);
//     for(int  v = 0; v < n_nodes; v++) coord[v] = _xyz[v+kc*n_nodes] + _dxdydz[v+kc*n_nodes];
  dtset=H5Dcreate(file,Namex.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT,H5P_DEFAULT
#endif
                 );
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,xxx);
  H5Sclose(dtsp); H5Dclose(dtset);

  // coord y -> x2
  std::ostringstream Namey; Namey << "NODES/COORD/X2";
  dtsp = H5Screate_simple(2, dimsf, NULL);
//     for(int  v = 0; v < n_nodes; v++) coord[v] = _xyz[v+kc*n_nodes] + _dxdydz[v+kc*n_nodes];
  dtset=H5Dcreate(file,Namey.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT,H5P_DEFAULT
#endif
                 );
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,yyy);
  H5Sclose(dtsp); H5Dclose(dtset);

  // coord z -> x3
  std::ostringstream Namez; Namez << "NODES/COORD/X3";
  dtsp = H5Screate_simple(2, dimsf, NULL);
//     for(int  v = 0; v < n_nodes; v++) coord[v] = _xyz[v+kc*n_nodes] + _dxdydz[v+kc*n_nodes];
  dtset=H5Dcreate(file,Namez.str().c_str(),H5T_NATIVE_DOUBLE,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT,H5P_DEFAULT
#endif
                 );
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,zzz);
  H5Sclose(dtsp); H5Dclose(dtset);

// Close the group
  H5Gclose(group_id2);   H5Gclose(group_id);
  delete []xxx; delete []yyy; delete []zzz;

  // connectivity ----------------------------------------------------------------------
  int  np_elem=0; int  n_elements=0;
    int * gelem_map;
  if(_dim==2) {
    np_elem=4; n_elements=(ny-1)*(nx-1);
    gelem_map=new int[n_elements*np_elem];
    for(int jy=0; jy<ny-1; jy++) {
      for(int ix=0; ix<nx-1; ix++) {
        int ielem= ix+jy*(nx-1);
// //       infile >> val;
//         // connectivity cartesian cell	  quad4
        gelem_map[ ielem*np_elem+0]=ix+jy*nx;
        gelem_map[ ielem*np_elem+1]=(ix+1)+jy*nx;
        gelem_map[ ielem*np_elem+2]=(ix+1)+(jy+1)*nx;
        gelem_map[ ielem*np_elem+3]=ix+(jy+1)*nx;
      }
    }
  }
  if(_dim==3) {
    np_elem=8; n_elements=(ny-1)*(nx-1)*(nz-1);
    gelem_map=new int[n_elements*np_elem];
    for(int kz=0; kz<nz-1; kz++) {
      for(int jy=0; jy<ny-1; jy++) {
        for(int ix=0; ix<nx-1; ix++) {
          int ielem= ix+jy*(nx-1)+kz*(nx-1)*(ny-1);
// //       infile >> val;
//         // connectivity cartesian cell	  quad4
          gelem_map[ ielem*np_elem+0]=ix+jy*nx+kz*nx*ny;
          gelem_map[ ielem*np_elem+1]=(ix+1)+jy*nx+kz*nx*ny;
          gelem_map[ ielem*np_elem+2]=(ix+1)+(jy+1)*nx+kz*nx*ny;
          gelem_map[ ielem*np_elem+3]=ix+(jy+1)*nx+kz*nx*ny;

          gelem_map[ ielem*np_elem+4]=ix+jy*nx+(kz+1)*nx*ny;
          gelem_map[ ielem*np_elem+5]=(ix+1)+jy*nx+(kz+1)*nx*ny;
          gelem_map[ ielem*np_elem+6]=(ix+1)+(jy+1)*nx+(kz+1)*nx*ny;
          gelem_map[ ielem*np_elem+7]=ix+(jy+1)*nx+(kz+1)*nx*ny;
        }
      }
    }
  }
  // Print mesh in hdf files
  std::ostringstream Name; Name << "MSHCONN";
//   hsize_t dimsf[2];
  dimsf[0] =n_elements*np_elem;  dimsf[1] = 1;
  dtsp = H5Screate_simple(2, dimsf, NULL);
  dtset = H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          ,H5P_DEFAULT,H5P_DEFAULT
#endif
                         );
  H5Dwrite(dtset,H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,gelem_map);
  H5Sclose(dtsp); H5Dclose(dtset);
  delete[]gelem_map;



  H5Fclose(file);
  return;
}

