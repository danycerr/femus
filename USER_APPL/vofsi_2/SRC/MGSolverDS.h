#ifndef __mgsolverddds_h__
#define __mgsolverddds_h__

#include "Equations_conf.h"

// =================================
 #ifdef DS_EQUATIONS
// =================================

// classe include ---------
#include "MGSolverDA.h"
#include "UserDS.h"
// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGFEMap;




// =================================================
/// Class for mg displ solvers
class MGSolDS: public MGSolDA {

private:
  // -----------------------------------------------
    
//   int _FSI_idx;    //< Navier-Stokes flag
//    int _DS_idx;     //< energy flag
    int _FF_idx[30];    //< field equation flag
  // auxilyary element data ------------------------
  // mesh
    int _nDSdim;
  const int  _offset;  ///< offset mesh nodes (top level)
  // physics
  const double _dt;     ///< time step
  const double _rhof;   ///< fluid density
  const double _uref;   ///< ref velocity
  const double _lref;   ///< ref length
  const double _Tref;   ///< ref temp
  const double _rhos;   ///< solid density
//   int    _dir;          ///< direction

  // -------------------- class field ----------------------------
  // element boundary conditions
 int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
 int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  
  int   _bc_el[NDOF_FEM*DIMENSION];
  
  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
 double  _ub_g[3][12];                     ///< external field  (0-1-2 degree)
 double  _ub_dxg[DIMENSION*DIMENSION];     ///< external field derivative  (0-1-2 degree) 
 double  _ub_old[12*NDOF_FEM];   
  double _xx_qnds[DIMENSION*NDOF_FEM]; double _xxb_qnds[DIMENSION*NDOF_FEMB];
   double _InvJac2[DIMENSION*DIMENSION];
    double _InvJac1[DIMENSION*DIMENSION];
    
    
    double _dphijdx_g2[DIMENSION];
    double _dphijdx_g1[DIMENSION];
  double _dphiidx_g2[DIMENSION];
  double _dphiidx_g1[DIMENSION];
   double _xyz_g[DIMENSION];
//  double  _InvJac2[DIMENSION*DIMENSION]; double  _InvJac1[DIMENSION*DIMENSION];
  // =============================================================

  // -----------------------------------------------------------
public:
  // Constructor-destructor --------------------------
  /// Level constructor
  MGSolDS(MGEquationsSystem& mg_equations_map_in,
	  const int vars_in[],
            std::string eqname_in="DS",
            std::string varname_in="d"
           );

  /// Destructor (level structure)
  ~MGSolDS() {}

  
   // Setting ---------------------------------------
  /// Boundary conditions
  void bc_read(int bc_gam,int mat_gam, double xp[],int bc_Neu[], int bc_value[]);
//  void bc_intern_read(double x[],int bc_Neu[],int u[]); ///< Read ic.
  /// Initial conditions
  void ic_read(int bc_gam,int mat_gam, double xp[],double u_value[]);

  void bc_intern_read (int bc_gam,int mat_gam,double xp[], int normal[],int bc[]);
   
  // Assemblying ------------------------------------
  /// Volume Assemblying.
  void GenMatRhs(const double time,const int Level,const int mode);
  
 
//   /// Surface Assemblying
//   void GenMatRhsB(const double time, const int Level,const int mode);

  // Multigrid function -----------------------------
  /// Timestep function
  void MGTimeStep(const double time, const int /*iter*/);
      /// MG time step solver (backward Euler).
  void MGTimeStep_nl_setup(
    const double time,               // time               <-
    const int    mode                // rhs assembler flag <-   
  ); 
    int MGTimeStep_nl_iter(
    const double time,               // time               <-
    const int    mode                // rhs assembler flag <-   
  ); 
      /// MG time step solver (backward Euler).
  void MGTimeStep_nl_sol_up(
    const double time,               // time               <-
    const int    mode                // rhs assembler flag <-   
  ); 
  inline void print_ext_data(double /* table_data*/[]) { };

  void matrixrhsvol_sol_ds(DenseMatrixM &KeM, DenseVectorM &FeM,int el_ndof[], int flag_group[]);
   void matrixrhsvol_liq_ds(DenseMatrixM &KeM, DenseVectorM &FeM,int el_ndof[], int flag_group[]);
//      void set_bc_matrix ( DenseMatrixM &KeM, DenseVectorM &FeM,int dir_maxnormal,
//                          int sur_toply[], int el_ndof[],int elb_ndof[], int elb_ngauss,  double u_old[], double normal[],double Ipenalty
//                        );
     double  MGFunctional(double p1,  double & p2); 
     void set_bc_matrix(
      DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
  int dir_maxnormal,                          ///<  normal dir
  int sur_toply[],                            ///< boundary topology map
  int el_ndof[],                              ///< number of volume dofs
  int elb_ndof[],                             ///< number of boundary dofs
  int elb_ngauss,                             ///<  number of surface gaussian points
  double normal[],                            ///< normal
  const int iaxis                              ///< axisymmetric
     );
     
//    void matrixrhsbc_liq_ds(DenseMatrixM &KeM, DenseVectorM &FeM,int el_ndof[]);
};

#endif  //  #endif DS_EQUATIONS
#endif

// ====================================================================================