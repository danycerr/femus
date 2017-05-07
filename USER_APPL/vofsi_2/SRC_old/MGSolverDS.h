#ifndef __mgsolverddds_h__
#define __mgsolverddds_h__

#include "Equations_conf.h"

// =================================
 #ifdef DS_EQUATIONS
// =================================

// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGFEMap;

/// DS boundary conditions ==============================================================
enum bound_condDS {
    simm_disp=0, disp_in0=1,disp_tg0=2,wall_fix=3,
    disp_in=5,disp_tg=6,wall_disp=7,
    free_disp=10, free_dispn=12,int_disp=11,slip=13,
};


// ======================================================================================
/// Class for mg displ solvers
// ======================================================================================
class MGSolDS: public MGSolDA {

private:
  // ------------------------------------------------------------------------------------
    
  int _FSI_idx;    //< Fluid structure interaction index
  int _DS_idx;     //< Displacement index
  
  // auxilyary element data -------------------------------------------------------------
  // mesh
  int _nDSdim;          ///< Problem dimension
  const int  _offset;   ///< offset mesh nodes (top level)
  // physics
  const double _dt;     ///< time step
  const double _rhof;   ///< fluid density
  const double _uref;   ///< ref velocity
  const double _lref;   ///< ref length
  const double _Tref;   ///< ref temp
  const double _rhos;   ///< solid density
  int    _dir;          ///< direction

  // -------------------- class field ---------------------------------------------------
  // element boundary conditions
  int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
  int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  
  int   _bc_el[NDOF_FEM*DIMENSION];
  
  // ------------------ integration -----------------------
  //  fields at gaussian points
 double  _ub_g[3][12];                     ///< external field  (0-1-2 degree)
 double  _ub_dxg[DIMENSION*DIMENSION];     ///< external field derivative  (0-1-2 degree) 
 double  _ub_old[12*NDOF_FEM];             ///< old solution  (0-1-2 degree) 
 double _xx_qnds[DIMENSION*NDOF_FEM]; double _xxb_qnds[DIMENSION*NDOF_FEMB];
 double _InvJac2[DIMENSION*DIMENSION];double _InvJac1[DIMENSION*DIMENSION];
 
    
  // shape functions
  double _dphijdx_g2[DIMENSION];  ///< quadratic derivaties shape
  double _dphijdx_g1[DIMENSION];  ///< quadratic derivaties shape
  double _dphiidx_g2[DIMENSION];  ///< quadratic derivaties shape
  double _dphiidx_g1[DIMENSION];  ///< quadratic derivaties shape
  double _xyz_g[DIMENSION];
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
  /// Initial conditions
  void ic_read(int bc_gam,int mat_gam, double xp[],double u_value[]);
  /// Internal boundary conditions
  void bc_intern_read (int bc_gam,int mat_gam,double xp[], int normal[],int bc[]);
  // Assemblying ------------------------------------
  /// Volume Assemblying.
  void GenMatRhs(const double time,const int Level,const int mode);
  ///moving mesh coarse algorithm
  void MoveMesh(const int Level);

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

  void matrixrhsvol_sol_ds(DenseMatrixM &KeM, DenseVectorM &FeM,int el_ndof[]);
  void matrixrhsvol_liq_ds(DenseMatrixM &KeM, DenseVectorM &FeM,int el_ndof[], int flag_group[]);
     
  void set_bc_matrix(DenseMatrixM &KeM, DenseVectorM &FeM,
  int dir_maxnormal,int sur_toply[], 
  int el_ndof[],int elb_ndof[],  int elb_ngauss, 
  double u_old[], double normal[],double Ipenalty
  );
};

#endif  //  #endif DS_EQUATIONS
#endif

// ====================================================================================