#ifndef __mgsolverp_h__
#define __mgsolverp_h__

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


// NSboundary conditions===============================================================================================
enum bound_cond_p {
              outflowp0=0,outflowp=4,
              vel_fix=10,interiorp=11
};





// ======================================================================================
// ====================== MGSolFSIP  =================================== =================
// ======================================================================================

#if (NS_EQUATIONS%2==0)

class MGSolNSP : public MGSolDA
{
private:
    // =====================  start data =======================================
    /// a) MGSolNS element data
    // mesh
    const int  _offset; 		///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
    int   _dir;        		///< x-direction for segrgated mode
    int   _nPdim;  //<dimension
    int _NS_idx;    //< Navier-Stokes flag
    
    // constant reference parameters
    const double _dt;    		///< =_mgutils.get_par("dt");
    const double _uref;  		///< =_mgphys.get_par("Uref");
    const double _lref;  		///< =_mgphys.get_par("Lref");
    // constant fluid properties
    const double _rhof;  		///< =_mgphys.get_par("rho0");
    const double _muf;  		///< =_mgphys.get_par("mu0");
    // nondimensional numbers
    double _IRe;         		///< =Reynolds number

//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
    int   _bc_vol[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags (vol int)
    int   _bc_bd[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags  (bd int)
//
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
    double  _ub_g[3][10];                  ///< external field  (0-1-2 degree)
    double  _ub_dxg[DIMENSION*DIMENSION];  ///< external field derivative  (0-1-2 degree)

    // ======================  end data ==========================================

    // ======================  start functions ===================================
public:



    /// b)   Init MGSolNS functions (constructor,destructor,external fields):
    MGSolNSP (                              ///< Constructor
        MGEquationsSystem& mg_equations_map, ///< equation map class (Mesh and parameters)
        int             nvars_in[],          ///< KLQ number of variables
        std::string     eqname_in="NS2P",     ///< base name system
        std::string     varname_in="p"       ///< base name variable
    );


    ~MGSolNSP() {};                        ///< Destructor

    /// c)          Read MGSolNS functions:
    /// This function  read bc
    void ic_read (
        int bc_gam,         ///< bc flag
        int bc_mat,           ///< material flag
        double xp[],          ///< xp[] is the NON-DIMENSIONAL node coordinates
        int iel,
        double u_value[]      ///< initial value (out)
    );

    /// This function reads ic
    void bc_read (
        int bc_gam,           ///< bc flag
        int bc_mat,           ///< material flag
        double xp[],          ///< xp[] is the NON-DIMENSIONAL node coordinates
        int bc_Neum[],        ///< normal
        int bc_flag[]         ///< boundary condition flag
    );

    /// d)  Assemblying MGSolNS Operators
    /// This function computes volume and surface assemblying
    void GenMatRhs (
        const double time,        // time
        const int Lev,            // Level
        const  int m              // rhs assembly control
    );
// Multigrid function ------------------------------------------------------------------
    /// This function is the  time-step manager function
    void MGTimeStep (
        const double time,         ///< time
        const int /*iter*/         ///< number max of iterations
    );
    /// MG time step solver (backward Euler).
    void MGTimeStep_nl_setup (
        const double time,               // time               <-
        const int    mode                // rhs assembler flag <-
    );
    /// MG time step solver (backward Euler).
    void MGTimeStep_nl_sol_up (
        const double time,               // time               <-
        const int    mode                // rhs assembler flag <-
    );
    /// MG time step solver (backward Euler).
    int MGTimeStep_nl_iter (
        const double time,               // time               <-
        const int    mode                // rhs assembler flag <-
    );
    // ====================================================================================
    /// This fucntion compute functional for optimal control
//   void MGFunctional(
//     const double time,
//     double starting_distance
//   ) {};

// void bc_intern_read(
//   int bc_gam,
//   int mat_gam,
//   double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
//   int bc_Neum[], // normal
//   int bc_flag[]  // boundary condition flag
// );
};
#endif // endif NS_EQUATIONS%2==0
#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
