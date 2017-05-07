
#ifndef __userDS_h__
#define __userDS_h__

#include "Equations_conf.h"
// ===================================
#ifdef DS_EQUATIONS
// ==================================

/*!  \defgroup DS_param   Class Table:  Navier-Stokes equation parameters (DS_param) */    
/// \ingroup DS_param
// ============================================================================
class DS_param
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
  int COMPRESSIBLE;
     
    DS_param(){
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
        COMPRESSIBLE=0;
        
    }
    inline void set_AXISYM(int val){AXISYM=val;}
    inline void set_SUPG(double    val) {SUPG=val;}
    inline void set_UPWIND(double  val) {UPWIND=val;}
    inline void set_UPWIND2(double val) {UPWIND2=val;}
    inline void set_UNSTEADY(int val){UNSTEADY=val;}
    
};



// // NSboundary conditions===============================================================================================
// enum bound_condDS {
//     simm_disp=0, disp_in0=1,disp_tg0=2,wall_fix=3,
//     disp_in=5,disp_tg=6,wall_disp=7,
//     free_disp=10, free_dispn=12,int_disp=11,slip=13,
// };
// NS boundary conditions=======================================================
/*!     \defgroup DS_Boundary_conditions     Enum Table: NS Boundary conditions  */ 
/// \ingroup DS_Boundary_conditions
    // ========================================================================
    /// Navier-Stokes boundary conditions
    enum bound_condDS {
//     fix_in0=      1,     ///<  1  normal disp inlet  \f$ {\bf l}.{\bf n}= 0  \f$
//     fix_tg0=      2,     ///<  2  normal disp inlet  \f$ {\bf l}.{\bf t}={ 0}{\bf t}  \f$
//     fix_disp0=  3,     ///<  3  normal disp inlet  \f$ {\bf l}={\bf 0}  \f$
//     fix_in=        5,     ///<  5  normal disp inlet  \f$ {\bf l}.{\bf n}={\bf l}_0.{\bf n}  \f$
//     fix_tg=        6,     ///<  6  normal disp inlet  \f$ {\bf l}.{\bf t}={\bf l}_0.{\bf n}  \f$
//     fix_disp=       7,     ///<  7  normal disp inlet  \f$ {\bf l}={\bf l}_0  \f$
//     free_disp= 10,     ///< 10  Neuman homogeneus         (dT.n=0)
//     inter=         11,     ///<       
//     free_disp_outlet=  12,     ///< 12  Neuman homogeneus         (dT.n=0)
//      free_wall_turb=        13,     ///< 13   Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
//      free_disp_p=        14,     ///< 14  Neuman nonhomogeneus      (dT.n=q_0)
//      free_disp_inlet=   16,     ///< 16  Neuman nonhomogeneus      (dT.n=q_0)
//     simm_dispx=            21,     ///< 21  simmetry      \f$ (p{\bf n}+{\bf tau} \cdot \widehat{n})\cdot \widehat{i}_x=0 \f$   
//     simm_dispy=            22,     ///< 22  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispz=            23,     ///< 23  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispxy=           24,     ///< 24  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispxz=           25,     ///< 25  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispyz=           26      ///< 26  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)     
   
    inter=                        11,     ///<  interiori point   
    
    free_disp=                22,     ///<  2(n) \f$  \tau_{nn}=p=0  \f$  +  2(t)  \f$ \tau_{nt}=0   \f$ 
    free_disp_outlet=    28,     ///<  2(n) \f$  \tau_{nn}=p=0  \f$  + 8(t) tg disp inlet  \f$ {\bf l}.{\bf t}= 0  \f$ 
 
    free_disp_p=            32,     ///<  3(n)    \f$ \tau_{nn}=p=p_0  \f$ +  2(t)  \f$  \tau_{nt}=0 \f$ 
    free_disp_inlet=      38,     ///<  3(n)  normal disp inlet  \f$ {\bf l}.{\bf n}= 0  \f$ + 8(t) tg disp inlet  \f$ {\bf l}.{\bf t}= 0  \f$ 
    
    fix_in0=                     88,     ///<  8(n)  normal disp inlet  \f$ {\bf l}.{\bf n}= 0  \f$ + 8(t) tg disp inlet  \f$ {\bf l}.{\bf t}= 0  \f$ 
    fix_tg0=                     88,     ///<  8(n) normal disp inlet  \f$ {\bf l}.{\bf n}= 0  \f$ +  8(t) tg disp inlet  \f$ {\bf l}.{\bf t}= 0  \f$ 
    fix_disp0=                 66,     ///<  6(n) lx=0  + 6(t) ly, lz=0  
    fix_in=                       98,     ///<  9(n)  normal disp inlet  \f$ {\bf l}.{\bf n}=  {\bf l}_0.{\bf n}  \f$ + 3(t) tg disp inlet  \f$ {\bf l}.{\bf t}= 0  \f$ 
    fix_tg=                       89,     ///<  8(n)  normal disp inlet  \f$ {\bf l}.{\bf t}=0  \f$+ 4(t) tg disp inlet  \f$ {\bf l}.{\bf t}=  {\bf l}_0.{\bf t}  \f$
     fix_disp=                  99,     ///<  9(n)  normal disp inlet  \f$ {\bf l}.{\bf n}=  {\bf l}_0.{\bf n}  \f$ +9(t) tg disp inlet  \f$ {\bf l}.{\bf t}=  {\bf l}_0.{\bf t}  \f$
    free_wall_turb=       84,     ///< 8(n)  normal disp inlet  \f$ {\bf l}.{\bf n}= 0  \f$ + 4(t)  \f$  \tau_{nt}=\alpha u_t  \f$ 
    
    simm_dispx=           82    ///< 82  simmetry   8   \f$ ({\bf l}.{\bf n}=0+  2(t)  \f$  \tau_{nt}=0 \f$    
};
    
  
  
#endif
#endif
