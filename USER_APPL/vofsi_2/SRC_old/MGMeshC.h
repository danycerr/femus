#ifndef __mgmesh1_h
#define __mgmesh1_h

// C++ include files that we need
#include <fstream>
#include <iostream>
#include <string.h>
#include "stdio.h"
#include "vof_config.h"
#include "Domain_conf.h"

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>
// #include "reference_counted_object.h"
// #include "libmesh_common.h"
class MGMeshC{

 public:
  // data ------------------------------
  unsigned int _dim;
  int _NXNYNZ[3];double _aLen[3];double _bLen[3];
  unsigned  int _NoLevels;
  int *_NoNodes; 
  int *_NoElements; 
  double *_xyz;  // coordinates
  int **_offset;   // offset map 
  int **_elem_map;// node map (connectivity)
  int _iproc;// node map (connectivity)
  // MGMeshC(const unsigned int NoLevels); 
  MGMeshC(const unsigned int dim_in,const unsigned int NoLevels_in); 
  ~MGMeshC();            // destructor
 void clear();
 
 
 
int  get_nx() const{return _NXNYNZ[0];}
int  get_ny() const{return _NXNYNZ[1];}
int  get_nz() const{return _NXNYNZ[2];}
 // -------------------------------------------
void  setNXNYNZ(int NXNYNZ_in[]){
  _NXNYNZ[0]=NXNYNZ_in[0];_NXNYNZ[1]=NXNYNZ_in[1];_NXNYNZ[2]=NXNYNZ_in[2];
  return;
}
void  setalen(double aXaYaZ_in[]){
  _aLen[0]=aXaYaZ_in[0];_aLen[1]=aXaYaZ_in[1];_aLen[2]=aXaYaZ_in[2];
  return;
}
void  setblen(double bXbYbZ_in[]){
  _bLen[0]=bXbYbZ_in[0];_bLen[1]=bXbYbZ_in[1];_bLen[2]=bXbYbZ_in[2];
  return;
}
 void init(const std::string& name);
 void write(const std::string& name,int level);
 // print
 void print(const unsigned int flag_print,const unsigned int Level);
  void print_hdf5(const unsigned int flag_print,const unsigned int Level);
  // ========================================================
/// This function prints the connectivity in hdf5 format
void print_hf5(
  const int  Level      // Level <-
) ;
 void print_nodes(std::ofstream& name,const unsigned int Level,const unsigned int mode);
 void print_conn(std::ofstream& name,const unsigned int Level,const unsigned int mode);
// void gen_c(int level,int  NNodes_val[], double ab[]);
// ======================================================================================
/// This function generates a Cartesian mesh.
void gen_c(
//   int  NXYZ0[],   // number of points (each dir) (NX,NY,NZ) in vof system
//   double bxbybz[] // dimension  (bx,by,bz) -> (0,bx)x(0,by)x(0,bz)
);

 private:
 void read_c(std::ifstream& infile);
 void write_c(std::ostream& infile,int level);
 
};

#endif
