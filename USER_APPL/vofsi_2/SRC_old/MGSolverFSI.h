#ifndef __mgsolvernfsi_h__
#define __mgsolvernfsi_h__

#include "Equations_conf.h"
// =================================
#ifdef FSI_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;


class FSI_param
{
public:
    // stabilization NS --------------------------------------------------------
    double FSI_SUPG;      // SUPG
    double FSI_UPWIND; // classical UPWIND
    double FSI_LES;         // smagorinsky-LES
    // non linear ---------------------------------------------------------------
    int    NL_ITER; // asymptotic non linear regime it
    double NL_TIME0; // initial time for  non linear regime
    int    NL_ITER0; // initial time  non linear regime it
    // time discretization
    double CRANK_NICK;
    int COMPRESSIBLE;
    

    

    void set_param() {
        // stabilization NS -----------------------------------------------------
        FSI_SUPG=1.0;       // SUPG
        FSI_UPWIND=1.0;        // classical UPWIND
        FSI_LES=.0;            // smagorinsky-LES
        // non linear ------------------------------------------------------------
        NL_ITER=0; // non linear it
        NL_TIME0=-.0001;
        NL_ITER0=0;
        // time discretization
        CRANK_NICK=1.;  // implicit (1) explicit (0) Crank-Nicolson (0.5)
        
        COMPRESSIBLE=0;
        // // 0 -> Incompressible
        // //  1-> Compressible
    }


};

// NSboundary conditions===============================================================================================
enum bound_cond {
    simm=0,velocity_in0=1,velocity_tg0=2,wall=3,
    simm1=4,velocity_in=5,velocity_tg=6,velocity=7,
    outflow=10,pressure_outlet=12,interior=11,wall_turb=13,outflow_p=14,pressure_inlet=16
};



/// Class MGSolFSI:
class MGSolFSI : public MGSolDA
{
private:
    // =====================  start data =======================================

    int  _nFSIdim;  //< dimension
    int _FSI_idx;    //< Navier-Stokes flag
    int _DS_idx;    //< Navier-Stokes flag
    int _T_idx;     //< energy flag
    int _K_idx;     //< turbulence flag

    /// a) MGSolFSI element data
    // mesh
    const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
    int   _dir;         ///< x-direction for segrgated mode

    // constant reference parameters
    const double _dt;    		///< = mg_equations_map_in.get_par("dt");
    const double _uref;  		///< = mg_equations_map_in.get_par("Uref");
    const double _lref;  		///< = mg_equations_map_in.get_par("Lref");
    double       _dirg[3];	///< =gravity
    // constant fluid properties
    const double _rhof;  		///< = mg_equations_map_in.get_par("rho0");
    const double _muf;   		///< = mg_equations_map_in.get_par("mu0");
    const double _rhos;     ///< mg_equations_map_in.get_par("rhos")
    const double _ni;  ///< mg_equations_map_in.get_par("nis")
    const double _Emod;  ///< mg_equations_map_in.get_par("Es")
    const double _hs;  ///< mg_equations_map_in.get_par("hs")


    // nondimensional numbers
    double _IRe;         		///< =Reynolds number
    double _IFr;         		///< =Froud number
    double _lambda;                     ///< =Lame' parameter
    double _mus  ;                           ///< mu_s parameter
    
    double _kappa_g[2];  		///< reference kappa
//   double _omega_g;    	///< reference omega
    double _y_dist;     		///< distance from the wall
    double _sP;         		///< turbulent tensor modulus
    double _mu_turb;    		///< effective turbulent viscosity
    double _T_nc;       		///< Reference temperature for Boussinesq term
    double _beta_mat;   		///< Thermal expansion coefficient
    double _grav[3];       	///< Gravity vector
    double _grav_mod;		///< Gravity modulus
    double  _y_bcout;
//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
    // element boundary conditions
    int   _bc_vol[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags (vol int)
    int   _bc_bd[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags  (bd int)
    int   _bc_el[NDOF_FEM* ( DIMENSION+1 )];
    // element geometry
    double _InvJac2[DIMENSION*DIMENSION];
    double _InvJac1[DIMENSION*DIMENSION];
    double _xx_qnds[NDOF_FEM*DIMENSION];
    double _xxb_qnds[NDOF_FEMB*DIMENSION];
//
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
    double  _ub_g[3][12];                  ///< external field  (0-1-2 degree)
    double  _ub_dxg[DIMENSION*DIMENSION];  ///< external field derivative  (0-1-2 degree)
    double _vof_phase;
    int _pres_order;
//   DenseMatrixM *_KeM;
//       DenseVectorM *_FeM;
    FSI_param   FSI_parameter;
    // ======================  end data ==========================================

    // ======================  start functions ===================================
public:
//  DenseMatrixM _KeM;
//      DenseVectorM _FeM;
    /// b)   Init MGSolFSI functions (constructor,destructor,external fields):
    MGSolFSI (                              ///< Constructor
        MGEquationsSystem& mg_equations_map, ///< equation map class (Mesh and parameters)
        int             nvars_in[],          ///< KLQ number of variables
        std::string     eqname_in="FS0",     ///< base name system
        std::string     varname_in="u"       ///< base name variable
    );
    ~MGSolFSI() {};                        ///< Destructor

    /// c)          Read MGSolFSI functions:
    /// This function  read bc
    void ic_read (
        int bc_gam,
        int bc_mat,
        double xp[],            // point coordinates
        int iel,
        double u_value[]        // point field values
    );

    /// This function reads ic
    void bc_read (
        int bc_gam,
        int bc_mat,
        double x[],             // point coordinates
        int bc_Neu[],           //  bc volume integral flag
        int bc_bd[]             //  bc surface integral flag
    );

    /// d)  Assemblying MGSolFSI Operators
    /// This function computes volume and surface assemblying
    void GenMatRhs (
        const double time,        // time
        const int Lev,            // Level
        const int m               // rhs assembly control
    );
// Multigrid function ------------------------------------------------------------------
    /// This function is the  time-step manager function
    void MGTimeStep (
        const double time,         ///< time
        const int /*iter*/         ///< number max of iterations
    );
    // ====================================================================================
    /// This function computes functional for optimal control
//   void MGFunctional(
//     const double time,
//     double starting_distance
//   ) {};
    // ====================================================================================
    /// This function sets the boundary condition and surface integration  in the local Matrix KeM
    void set_bc_matrix ( DenseMatrixM &KeM, DenseVectorM &FeM,int dir_maxnormal,
                         int sur_toply[], int el_ndof[],int elb_ndof[], int elb_ngauss,  double u_old[], double normal[],double Ipenalty
                       );
    // ==============================================================================================
    /// This function sets the boundary condition in the local Matrix KeM
    void matrixrhsvol_liq ( DenseMatrixM &KeM, DenseVectorM &FeM,
                        int el_ndof[],const int mode,const double eul_implicit,const double h_eff,
                        double u_old[],  double u_oold[],    double u_nl[],  double p_proj[], double dp_proj[], int flag_interface
                      );

    void matrixrhsvol_sol ( DenseMatrixM &KeM, DenseVectorM &FeM,
                            int el_ndof[],const int mode,const double eul_implicit,const double h_eff,
                            double u_old[],  double u_oold[],    double u_nl[],  double p_proj[], double dp_proj[], int flag_interface
                          );
    // ==============================================================================================
    void get_el_field_data (
        int iel, int Level,
        int el_conn [],int offset,int el_dof[],int ndof_lev,
        double u_old[],  double u_oold[],    double u_nl[],  double p_proj[], double dp_proj[]
    );

#ifdef TBK_EQUATIONS // Turbulence only (4-5) --------------------------
    double eval_var1 ( double val[] ); // computation val1=_mu_turb
    double eval_var2 ( double vel_bound[] ); // computation val2=_utau
    double eval_var3 ( double prop[] ); // computation val3=musker
//   double eval_val3(double dist, double utau, double nu);// computation val3=musker
#endif
};


#endif // endif FSI_EQUATIONS
#endif // endif _mgsolverns_h
