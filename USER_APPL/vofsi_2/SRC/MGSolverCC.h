#ifndef __mgsolvercc_h__
#define __mgsolvercc_h__
#include "vof_config.h"
#include "Equations_conf.h"
#include "Solverlib_conf.h"
// C++ libaries
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <float.h>
#include <stdlib.h>
#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionDD.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif
// Laspack
#include <errhandl.h>
#include <qvector.h>
#include <matrix.h>
#include <qmatrix.h>
#include <operats.h>
#include <factor.h>
#include <precond.h>
#include <rtc.h>
#include <itersolv.h>
//#include "itersolv2.h"
#include <mlsolv.h>
#include <lastypes.h>
#include <qvector.h>

// Local Includes
#include "MGMeshC.h"
// #include "MGSolver.h"
// *********************************************
class MGSol;
#define Real double

/// color function class solver
class MGSolCC {
// #ifdef TWO_PHASE
  // data -------------------------------
private:
  // Data for contraction and expansion matrix *********
  double *_tempc1;///< temp memory for expanding  CMatrix
  double *_tempc2;///< temp memory for expanding  CMatrix
  double *_tempmx1;///< temp memory for expanding  CMatrix
  double *_tempmy1;///< temp memory for expanding  CMatrix
#if DIMENSION==3
  double *_tempmz1;///< temp memory for expanding  CMatrix
#endif
  /// Matrix contraction
  void CtrRow(double tempc2[],Matrix * mtrx,unsigned int Level,unsigned int jy,unsigned int rowl);
  /// Row expansion
  void ExpRow(Matrix * mtrx,double temp[],unsigned int Level,unsigned int jy,unsigned int rowl);
  // ******************************************************  
public:
  int _name;
  int _dim;
  // cc matrix (mesh level)
  QVector _uvw;     ///< velocity field
  QVector _dsvw;     ///< displacement field
  QVector cc;      ///<  cc matrix        (QVector)
  QVector cc_old ; ///<  cc old matrix    (QVector)
  QVector _mx;     ///<  x-normal matrix  (QVector)
  QVector _my;     ///<  y-normal matrix  (QVector)
#if DIMENSION==3
  QVector  _mz;///<  normal z matrix (QVector)
#endif
  int *_invnode_dof;    ///<  map u-> cc

  // c matrix (multilevel)
  int _Nlev_cc;///< Number of c levels
  size_t *_nxyz[DIMENSION];///< c matrix dimension y
  Matrix *c1;   Matrix *c1_old;   ///<  c matrix  (CMatrix)
  Matrix *_mx1; Matrix *_mx1_old; ///<  mx matrix (CMatrix)
  Matrix *_my1; Matrix *_my1_old; ///<  my matrix (CMatrix)
#if DIMENSION==3
  Matrix *_mz1;  ///<  mz matrix (CMatrix)
  Matrix *_mz1_old;  ///<  mz matrix (CMatrix)
#endif
  // ====================================================================================
  //                               Constructor Destructor Init
  // ====================================================================================
  
  //-------------------------------------------------------------------------------------
  MGSolCC( ///< Constructor. This function set the dimensions of the multilevel cc structures
    const unsigned int NoLevels_cc_in ///< number of cc levels -> _Nlev_cc (in)
  ,int name
  );
  //-------------------------------------------------------------------------------------
  
  void init(///< This function initializes the cc color vector and the 
            ///  multilevel cc structure.
    const int ne_xyz[]
  );
  //-------------------------------------------------------------------------------------
  /// initial level function 3
   void init_level(
     const int Level,
     const int ne_xyz[]
   );
  //-------------------------------------------------------------------------------------
  /// initial dof function 4
  void init_dofCA(const unsigned int Level);
 //-------------------------------------------------------------------------------------
  // Destructor
  ~MGSolCC();     ///< destructor (-> clear) 1
  void clear();   ///< clear 2
  
  
 // field ==============================================================================
  //-------------------------------------------------------------------------------------
  // velocity field
  void GenVel(
//     const MGMeshC & mgmesh,
//     QVector &sol,
    const int Level,
    const double dt);
  //-------------------------------------------------------------------------------------
  // Color function soltion
  /// Generating  solution
  void GenSol(const int Level);
  //-------------------------------------------------------------------------------------
  void init_drop(unsigned int Level,double xc[],double r,double min=-1000.);
  /// Generating old solution
  void GenOldSol(const int Level);
  void WriteSol(const int Level,const std::string& name,const int nvars);
  void ReadSol(const int Level,const std::string& name,const int nvars);
  void read_fine(const unsigned int flag_print,const unsigned int Level);
  /// old cc color function  update
  void OldSol_update(const int Level,Matrix & sol,Matrix & sol_old);
  void CCSol_update(); ///< cc color function  update

  // Normal
  /// Normal  function
  void GenNorm(const unsigned int Level, void (MGSolCC::*pt2Member)(const unsigned int,int indxy [],double mxy[]));
  void Normal_update();///< Normal update

  // ====================================================================================
  // Multlevel solver
  /// Multilevel solver
  void MGSolve(/*MGMeshC &mgmesh,*/unsigned int nc_step,unsigned int itime,double dt);
  // Restriction - Projection Operators
  void RestSol(const int Level);///<  Restriction Operator
  void ProjSol(const int Level);///<  Projection Operator
  void ProjNormal(const int Levf); ///< Unit Normal projector
  void RestNormal(const int Levf); ///< Unit Normal restrictor

  // ---------------------------------------------
  // Computation cell quantities
  /// Computation of cell line intersect
  double get_alpha(double m1, double m2, double m3, double cc);
  double get_alpha2(double m1, double m2, double m3, double cc);

  /// Cell area computation
  double get_area(double mx,double my,double alpha,double x0,double y0);
  double get_vol(double mx,double my,double mz,double alpha,double x0,double y0,double dx, double dy);

  /// Computaton cell velocity field
  void get_vel(int ix[],int fd,double u[]);
  /// Computaton cell displacement field
  void get_disp(int ix[],int fd,double u[]);
  /// Return interface points
  void get_pts(double mx, double my, double alpha,double *pt1,double *pt2);
  /// return the surface tension
  void get_2surten(double xp[],double ff[],int ord[]);
  /// return the phase
  double get_2phase(int blck, double xp[]);

   void setFieldSource(
     double dt,
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * srcField,
    InterfaceFunctionDD * fct 
  );
         void setFieldSource_disp(
     double dt,
    int interface_name,
    int n_cmp,
    const std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField,
//     const ParaMEDMEM::MEDCouplingFieldDouble * srcField,
    InterfaceFunctionDD * fct 
  );
   void SourceVel_ext(/*const MGMeshC & mgmesh,*/double uvw_field[],  const double dt);
// ======================================================================================
//                   I/O  Print Read
// ======================================================================================
   
// --------------------------------------------------------------------------------------
    void read_fine_hdf5(
    const unsigned int flag_print,
    const unsigned int Level
  );
// --------------------------------------------------------------------------------------   
/// This function prints the time xml file (to run the single time step file cc.#.xmf)
void  print_time_cc_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,    ///<  print every print_step
    const int ndigits       ///< number of digit (namefile)
);
 
// --------------------------------------------------------------------------------------
/// This function prints the cc color file xmf for  cc (cc vector) 
/// It is called by the print fine function (print_fine_hdf5())
void print_cc_xmf(
  double time,
  std::string &flag,                    ///<    print flag (time)
  const unsigned int Level     ///< fine Level
) ;
// --------------------------------------------------------------------------------------
/// This function prints the cc color function.
void print_cc_hdf5(
  double time,
  std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level
);

// --------------------------------------------------------------------------------------
/// This function prints the time xml file (to run the single time step file ccf.#.xmf)
void  print_time_fine_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,    ///<  print every print_step
    const int ndigits       ///< number of digit (namefile)
);
// -------------------------------------------------------------------------------------- 
/// This function prints the file xmf for fine cc (c1[Level]). Only for cells 0 < C < 1.
void print_fine_xmf(
  double time,
  std::string flag,                    ///<    print flag (time)
  int n_cell_fine,             ///< number of fine cells with 0 < C < 1
  const unsigned int Level     ///< fine Level
);
// --------------------------------------------------------------------------------------
 /// This function prints the fine color function cc (c1[Level])
/// and the file xmf. Only for cells 0 < C < 1.
void print_fine_hdf5(
  double time,
 std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level 
);

#if DIMENSION==2

// --------------------------------------------------------------------------------------
/// This function prints the time xml file (to run the single time step file ccf.#.xmf)
void  print_time_interface_xmf(
  const int t_init,       ///<  intial time
  const int n_time_step,  ///<  number of time steps
  const int print_step,    ///<  print every print_step
    const int ndigits       ///< number of digit (namefile)
);
// -------------------------------------------------------------------------------------- 
/// This function prints the file xmf for fine cc (c1[Level]). Only for cells 0 < C < 1.
void print_interface_xmf(
  double time,
  std::string flag,                    ///<    print flag (time)
  int n_cell_fine,             ///< number of fine cells with 0 < C < 1
  const unsigned int Level     ///< fine Level
  
);
// --------------------------------------------------------------------------------------
 /// This function prints the fine color function cc (c1[Level])
/// and the file xmf. Only for cells 0 < C < 1.
void print_interface_hdf5(
  double time,
  std::string flag,                   ///<    print flag (time)
  const unsigned int Level    ///<    fine Level 
);
  
#endif
// --------------------------------------------------------------------------------------
/// This function prints c1 Matrix at level L in vtk format
  void print_fine(const unsigned int flag_print,const unsigned int L);
  /// This function prints the fine color file xmf. Only for cells 0 < C < 1.
void print_hf5(
  const int  Level      // Level <-
) ;
  
  // Vof Reconstruction =================================================================

  /// Young reconstruction 
  void rec_Young(const unsigned int,int indxy [],double mxy[]);
  ///  Elvira reconstruction
  void rec_elv1(const unsigned int Level,int indxy[],double mxy[]);
  ///  Alpha reconstruction
  void rec_melv1(const unsigned int Level,int ixy[],double mxy[]);
  ///  central reconstruction
  void rec_Cent(const unsigned int Level,int ixy[],double mxy[]);

  // Vof Advection ======================================================================
  
  /// Advection driver function
  void Adv(unsigned int Level,double dt, unsigned int itime=0);
  /// Split advection along the x-direction
  void lagrangeX(const unsigned int Level,const double dt, unsigned int itime=0);
  /// Split advection along the y-direction
  void lagrangeY(const unsigned int Level,const double dt, unsigned int itime=0);

#if DIMENSION==3
  /// Split advection along the z-direction
  void  lagrangeZ(const unsigned int Level,const double dt, unsigned int itime=0);
  /// Volume computation in 3D
  double get_vol3D(double m1,double m2,double m3,double alpha,double r0,
                   double dr0);
#endif
#if DIMENSION==2  
  // 2D Advection ----------------------------------------
  /// Unsplit (uncorrect) advection along the xydirection
  void AdvVelLXY(unsigned int Level, double dt);
  /// Unsplit (uncorrect) advection along the xydirection
  void AdvVelLXY1(unsigned int Level, double dt);
  /// Unsplit (uncorrect) advection along the yxdirection
  void AdvVelLYX1(unsigned int Level, double dt);
  /// Boundary conditions
  
#endif
  bool bc_outlet_read(int ix,  int jy,int indl,
    double hx, double hy,
    double cc1,
    int boundary/*,  
    double tempc1[],
    double tempmx1[],
    double tempmy1[]*/);
  
  void ic_read(
    double _tempc2[],
    unsigned int i, unsigned int j, unsigned int k,
    double hx, double hy,double hz,
    unsigned int nadd);
  
  void internal_read(
    double tempc1[],int ix,  int jy, int kz,
    int indl,
    double hx, double hy,double hz,
    double cc1);

};

#endif // end TWO_PHASE
