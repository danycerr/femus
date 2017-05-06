#ifndef __mgsolverns_h__
#define __mgsolverns_h__

#include "Equations_conf.h"
// =================================
#ifdef NS_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;



#include "UserNS.h"

// ===========================================================================
//                                Navier-Stokes equation class
//=============================================================================
/// Navier-Stokes equation  (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)  
class MGSolNS : public MGSolDA
{ 
/// Class for mg Navier-Stokes equation  with name NS_EQUATIONS. 
/// Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)  
// ===========================================================================
private:
    // =====================  start data =======================================

    int   _nNSdim;  //<dimension
    int _FF_idx[30];    //< field equation flag
    const double _dt;                  ///< =_mgutils.get_par("dt");

    /// a) MGSolNS element data
    // mesh
    const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
    int   _dir;         ///< x-direction for segrgated mode

    // constant reference parameters and  fluid properties (=_mgutils.get_par(" "))
    const double _uref /**< "Uref" */;
    const double _lref /**< "Lref" */;
    const double _rhof /**< "rhof" */;
    const double _muf  /**< "muf"  */;
    // nondimensional numbers
    double _IRe /**< Reynolds number */;
    double _IFr /**<  =Froud number */;
    double       _dirg[3]  /**< =gravity */;
    // turbulence
    double _kappa_g[2] /**< reference kappa */;
    double _y_bcout /**< wall distance */;
    double _mu_turb /**< turb eff */;
    double _sP	/**< tensor modulus*/;

//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
    int   _bc_vol[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags (vol int)
    int   _bc_bd[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags  (bd int)
    int   _bc_el[NDOF_FEM* ( DIMENSION+1 )];
    double _InvJac2[DIMENSION *DIMENSION];
    double _InvJac1[DIMENSION *DIMENSION];
    double _xx_qnds[NDOF_FEM *DIMENSION];
    double _xxb_qnds[NDOF_FEMB *DIMENSION];
    double _xyzg[DIMENSION];
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
    double  _ub_g[3][12];                  ///< external field  (0-1-2 degree)
    double  _ub_dxg[DIMENSION *DIMENSION]; ///< external field derivative  (0-1-2 degree)

    int _pres_order;
    NS_param   _NS_parameter;
    // ======================  end data ==========================================

    // ======================  start functions ===================================
public:
//  DenseMatrixM _KeM;
//      DenseVectorM _FeM;
    /// b)   Init MGSolNS functions (constructor,destructor,external fields):
    MGSolNS (                              ///< Constructor
        MGEquationsSystem &mg_equations_map, ///< equation map class (Mesh and parameters)
        int             nvars_in[],          ///< KLQ number of variables
        std::string     eqname_in="NS0",     ///< base name system
        std::string     varname_in="u"       ///< base name variable
    );
    ~MGSolNS() {};                        ///< Destructor

    /// c)          Read MGSolNS functions:
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

    /// d)  Assemblying MGSolNS Operators
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
    /// This fucntion compute functional for optimal control
//   void MGFunctional(
//     const double time,
//     double starting_distance
//   ) {};
    // ====================================================================================
//     void set_bc_matrix(
//     DenseMatrixM &KeM, DenseVectorM &FeM,
//     int dir_maxnormal,int sur_toply[],
//     int el_ndof[],int elb_ndof[], int elb_ngauss,
//     double u_old[], double normal[]
//                       );

    void set_bc_matrix (
        DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
        int dir_maxnormal,                          ///<  normal dir
        int sur_toply[],                            ///< boundary topology map
        int el_ndof[],                              ///< number of volume dofs
        int elb_ndof[],                             ///< number of boundary dofs
        int elb_ngauss,                             ///<  number of surface gaussian points
        double normal[],                            ///< normal
  const int iaxis                              ///< axisymmetric
    );

    // ==============================================================================================
    void matrixrhsvol (
        DenseMatrixM &KeM, DenseVectorM &FeM,
        int el_ndof[],
        double u_old[],  double u_oold[],   double u_nl[], double p_proj[], double dp_proj[],
        const int unsteady_flag, const int axysim,
        const double les,const int mode
    );



    // ==============================================================================================
    void get_el_field_data (
        int iel, int Level,
        int el_conn [],int offset,int el_dof[],int ndof_lev,
        double u_old[],  double u_oold[],    double u_nl[],  double p_proj[], double dp_proj[]
    );

// #ifdef TBK_EQUATIONS // Turbulence only (4-5) --------------------------
    double eval_var1 ( double val[] ); // computation val1=_mu_turb
    double eval_var2 ( double vel_bound[] ); // computation val2=_utau
    double eval_var3 ( double prop[] ); // computation val3=musker
//   double eval_val3(double dist, double utau, double nu);// computation val3=musker
// #endif
};


#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
