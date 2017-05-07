#include "Equations_conf.h"



// ======================================================================================
#ifdef DS_EQUATIONS // 3D-2D Dsiplacement  equation
// ======================================================================================

#include "MGSolverDS.h"
#include "MGGeomEl.h"

// configuration files -----------
#include "MGFE_conf.h"
#include "Printinfo_conf.h"
#include "MGEquationsSystem.h"



// local include -------------------
 #include "MGSystem.h"


// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolDS::bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
) {// ===================================
  return;
}

void MGSolDS::ic_read(
  int /*bc_gam*/,
  int /*mat_gam*/,
  double /*xp*/[],
  double u_value[]
) {// =======================================
   u_value[0] = 0.;
  return;
}

// ============================================================================
/// This function  defines the boundary conditions for the DS system:
// ============================================================================
void MGSolDS::bc_read(
  int bc_gam,
  int mat_gam,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[],        // normal
  int bc_flag[]         // boundary condition flag
) {// =========================================================================
///   enum bound_condDS -------------------------------------------------------
///   Diriclet  (point constraint < 10) ---------------------------------------
///   simm=0, disp_in0=1,disp_tg0=2,wall_fix=3, -------------------------------
///   disp_in=5,disp_tg=6,wall_disp=7, ----------------------------------------
///   Neumann (integratral constraint >9) -------------------------------------
///   free_disp=10, free_dispn=12,interior=11,slip=13, ------------------------
// ----------------------------------------------------------------------------

  bc_Neum[0]= wall_fix;bc_flag[0]=0;
  if      (bc_gam< 20)  {bc_Neum[0]= wall_fix;  }  //fixed displacement fluid boundary
  else if (bc_gam<30)   {bc_Neum[0]= free_disp; }  //free movinng surface solid
  else if (bc_gam>900)  {bc_Neum[0]= wall_fix;  }  //fixed interaface 
  if      (bc_gam== 22)  {bc_Neum[0]= wall_fix;  }  //fixed displacement fluid boundary
  return;
 }
#endif
// ============================================================================
// ============================================================================
// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 