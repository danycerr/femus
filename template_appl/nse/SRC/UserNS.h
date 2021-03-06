
#ifndef __userNS_h__
#define __userNS_h__

#include "Equations_conf.h"
// ===================================
#ifdef NS_EQUATIONS
// ==================================

/*!  \defgroup NS_param   Class Table:  Navier-Stokes equation parameters (NS_param) */    
/// \ingroup NS_param
// ============================================================================
class NS_param
{//< This class defines the physical and numerical  energy equation parameters
public:
    
    
    int AXISYM  /**< axisymmetry  */;
    // stabilization NS --------------------------------------------------------
     double SUPG/**< supg */; double UPWIND;/**< Normal Upwind */ double UPWIND2;/**< Transv Upwind */
     
    // turbulence
    double LES/**< les */; double DIST_FIX/**< distance from wall */;
    
   // non linear -------------------------------------------------------------
    int    NL_ITER;/**< NON LIN IT */    int    NL_ITER0;/**< INIT NON LIN IT */
    double NL_TIME0; // initial time for  non linear regime
    double H_EXT/**< h convective in bc  */ ;double V_EXT/**< v convective in bc  */;
    
    // time discretization
    double CRANK_NICK/**< implicit (1) explicit (0) Crank-Nicolson (0.5)*/; 
    int UNSTEADY;    /**< un(1)/steady(0) flag  */
  
     
    NS_param(){
        AXISYM=0;
         // stabilization NS -----------------------------------------------------
        SUPG=1.0/**< SUPG */; UPWIND=.0/**< Normal Upwind */; UPWIND2=.0/**< Transv Upwind */;
        
        // turbulence
        LES=.0/**< LES */;DIST_FIX=1.e-3/**< distance from wall */;
        // non linear ------------------------------------------------------------
        
        NL_ITER=0;/**< NON LIN IT */ NL_ITER0=0;/**< INIT NON LIN IT */
        NL_TIME0=-.0001;
        
        // bc  -------------------------------------------------------------
        H_EXT=1.;/**< h convective in bc  */ V_EXT=1.;/**< v convective in bc  */
        
          // time discretization
        CRANK_NICK=1.;  /**< implicit (1) explicit (0) Crank-Nicolson (0.5) */
        UNSTEADY=1;     /**< un(1)/steady(0) flag  */
       
        
    }
    inline void set_AXISYM(int val){AXISYM=val;}
    inline void set_SUPG(double    val) {SUPG=val;}
    inline void set_UPWIND(double  val) {UPWIND=val;}
    inline void set_UPWIND2(double val) {UPWIND2=val;}
    inline void set_UNSTEADY(int val){UNSTEADY=val;}
    
};




// NS boundary conditions=======================================================
/*!     \defgroup NS_Boundary_conditions     Enum Table: NS Boundary conditions  */ 
/// \ingroup NS_Boundary_conditions
    // ========================================================================
    /// Navier-Stokes boundary conditions
    enum bound_cond {
    velocity_in0=      1,     ///<  1  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$
    velocity_tg0=      2,     ///<  2  normal velocity inlet  \f$ {\bf v}.{\bf t}={ 0}{\bf t}  \f$
    wall=              3,     ///<  3  normal velocity inlet  \f$ {\bf v}={\bf 0}  \f$
    velocity_in=       5,     ///<  5  normal velocity inlet  \f$ {\bf v}.{\bf n}={\bf v}_0.{\bf n}  \f$
    velocity_tg=       6,     ///<  6  normal velocity inlet  \f$ {\bf v}.{\bf t}={\bf v}_0.{\bf n}  \f$
    velocity=          7,     ///<  7  normal velocity inlet  \f$ {\bf v}={\bf v}_0  \f$
    
    outflow=          10,     ///< 10  Neuman homogeneus         (dT.n=0)
    interior=         11,     ///<       
    pressure_outlet=  12,     ///< 12  Neuman homogeneus         (dT.n=0)
 
    wall_turb=        13,     ///< 13   Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
    outflow_p=        14,     ///< 14  Neuman nonhomogeneus      (dT.n=q_0)
    pressure_inlet=   16,     ///< 16  Neuman nonhomogeneus      (dT.n=q_0)
    
    simmx=            21,     ///< 21  simmetry      \f$ (p{\bf n}+{\bf tau} \cdot \widehat{n})\cdot \widehat{i}_x=0 \f$   
    simmy=            22,     ///< 22  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simmz=            23,     ///< 23  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simmxy=           24,     ///< 24  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simmxz=           25,     ///< 25  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simmyz=           26      ///< 26  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)     
};
    
  
  
#endif
#endif