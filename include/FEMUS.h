#ifndef __FEMUS__
#define __FEMUS__

#include <vector>
#include "Solverlib_conf.h"

#ifdef HAVE_MED
namespace ParaMEDMEM {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
class DataArrayInt;
class DataArrayDouble;
}
#endif

class MeshExtended;
class MGFemusInit;
class EquationSystemsExtendedM;
class MGSystem;
class MGUtils;
class MGGeomEl;
class MGFEMap;
class MGTimeLoop;

class BoundInterp;
#ifdef   TWO_PHASE
class MGSolCC;
#endif
#ifdef TBK_EQUATIONS
class TurbUtils;
#endif

class FEMUS {
  // ==========================================================================
  //=========   data  =========================================================
protected:

  //data communication (defined in Constructor) -------------------------------
  MGFemusInit *    _start;          // start function
  MPI_Comm         _comm;           // communicator
  bool             _local_MPI_Init; // initial mpi flag

  // param and file data (defined in init(.,.)) ------------------------------------------
  MGUtils *     _mg_utils;     // param and files

  // fem ------- (defined in init(.,.))  --------------------------------------------------
  MGGeomEl *                 _mg_geomel;    // set element type
  MGFEMap  *                 _mg_femap;     // fem map

  // data meshes --------------------------------------------------------------------------
//   int                                    /*_num_mgmesh*/;   // num of MGmeshes
  MeshExtended*                          _mg_mesh;      // FEMus-mesh
#ifdef HAVE_MED
  ParaMEDMEM::MEDCouplingUMesh *         _med_mesh;     // Med-mesh
#endif
  // system data ---------------------------------------------------------------------------
  //   _LibMeshProblem * _problem;
  EquationSystemsExtendedM*           _mg_equations_map; // system
  MGTimeLoop *                        _mg_time_loop;     // transient


  // ==========================================================================
  //=========  public functions  ============================================
public:

  // Constructor-Destructor
  FEMUS();
  FEMUS(MPI_Comm comm);
  #ifdef TBK_EQUATIONS
   void init_param(MGUtils &mgutils, TurbUtils &Parameters,int name=0 );
#endif
  void init_param(MGUtils &mgutils, int name=0);
  void init_fem(MGGeomEl & mggeomel,MGFEMap & mgfemap);
  ~FEMUS();
  void terminate();

  void setMesh();
  void setSystem(const std::vector<FIELDS> & pbName,
     int n_data_points=0,
  int n_data_cell=0
  );
  
  
  //====================================================================
/// This function gets the proc
 int  get_proc() const; 
 //====================================================================
/// This function sets up the intial set
  void solve_setup(
    int        & t_in,                 ///< initial time iteration
    double     &  time                 ///< actual time
  );
  //=============================================================================
// This function solves the problem
//   void solve();
//=============================================================================
// This function solves one step  for transient problems
  void solve_steady(
       const int & nmax_step,  ///< number max of steps
  const double & toll,  ///< tolerance
    const int  & it_step,    ///< actual num iteration
    const int  & print_step,  ///< print every
    double     &  dt,    ///< inial time step 
    const int  & eq_min=0, ///< eq min to solve -> enum  FIELDS (equations_conf.h) 
    const int  &  eq_max=30 ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  ) ; 
  
  
  
   void set_uooold(
   const int & flag,  ///<  0 xold-> x_ooold   1 x_ooold-> xold
  const double & toll,  ///< tolerance
  const double delta_t_step_in,  //   (in)
  const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int     &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
);
//=============================================================================
// This function solves one step  for transient problems
  void solve_onestep(
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double     &  time,                ///< actual time
    double     &  dt,                   ///< step time
    const int  & eq_min=0, ///< eq min to solve -> enum  FIELDS (equations_conf.h) 
    const int  &  eq_max=30 ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  ) ;
  
  //=============================================================================
// This function solves one step  for transient problems
  void solve_non_linear_onestep(
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double     &  time,                ///< actual time
    double     &  dt                   ///< step time
  ) ;
  
  void solve_control_onestep(
  const int  & t_in,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt                   ///< step time
);
  
   const  MeshExtended & get_MGMesh(){return *_mg_mesh;};
//    const EquationSystemsExtendedM& get_MGExtSystem(){return *_mg_equations_map;};
   EquationSystemsExtendedM& get_MGExtSystem(){return *_mg_equations_map;};
   
    double System_functional(
  const int  & ff,                 ///< initial time iteration
   double      parameter,             ///< functional parameter
   double     & control                   ///< step control
);
    std::string MeshName;
  inline void setMeshName(std::string meshname){MeshName = meshname;};
  inline std::string getMeshName(){return MeshName;};
  inline string GetMeshDir(){return _mg_utils->_mesh_dir;};
  
  
  #ifdef   TWO_PHASE
void set_mgcc(MGSolCC & cc);
#endif
// =================  with MED lib ===============================================
// =================  with MED lib ===============================================   
#ifdef HAVE_MED
  void setMedMesh(const std::string & dataFile);
  /// This functiongset the Med-Mesh
  const ParaMEDMEM::MEDCouplingUMesh &   getMedMesh() {
    return  *_med_mesh;
  };
  
    void init_interface(
    const int interface_name,
    const int interface_id1,
    const int interface_id2,
    int order_cmp,
    const std::string & medfile_name, // medfile name    (in)
    bool on_nodes=true,
    const int index_medmesh=0                      // med-mesh index  (in)
  ) ;
  //   init interface with two groups at ones

  void init_interface(
    const int interface_name,
    int interface_id,
    int order_cmp,
    const std::string & medfile_name, // medfile name    (in)
    const std::string & medfile_name_old, // medfile name    (in)
    const FEMUS & P_old,
    const int index_mgmesh,           // mgmesh index    (in)
    const int index_medmesh=0                      // med-mesh index  (in)
  ) ;
// // ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces function
// through the index of the interface-functions
  //! Function for creating an interface
  void init_interface(
    const int interface_name,          /**< Interface name */
    const int interface_id,            /**< Name of the group */
    int order_cmp,                     /**< Order (piecewise, linear, quadratic) */
    const std::string & medfile_name,  /**< name of the med file */
    bool on_nodes=true,                /**< Interface created on nodes */
    const int index_medmesh=0          /**< Index of the med mesh */
  ) ;
  
  //! Function for creating an interface
//   void init_interface(
//     const int interface_name,             /**< Interface name */
//     int interface_id,                     /**< Name of the group */
//     int order_cmp,                        /**< Order (piecewise, linear, quadratic) */
//     const std::string & medfile_name,     /**< medfile name    (in) */
//     const std::string & medfile_name_old, /**< medfile name    (in) */
//     const FEMUS & P_old,                  /**< Old FEMus problem          */
//     const int index_mgmesh,               /**< mgmesh index    (in) */
//     const int index_medmesh=0             /**< med-mesh index  (in) */
//   ) ;

  void init_interface(
  const int interface_name,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  bool on_nodes=true
) ;

void init_interface(
  const int interface_name,
  int order_cmp,
  ParaMEDMEM::MEDCouplingUMesh * support,
  ParaMEDMEM::DataArrayDouble* ACC,
  int interface_id = 0
);


  // =========================================================================
  /// This function sets the value from first_cmp to end_cmp of the field
  /// on the old solution x_old of the interface id
  void write_Boundary_value(
    int id_boundary_name   ,    ///< identity interface name (in)
    std::string mgsystem_name, ///< system name          (in)
    int n_cmp,             ///< from variable system (in)
    int first_cmp =0             ///< to variable system   (in)

  );


// ============================================================================
  /// This function sets the field of interface function (interface_id)
  /// with an analitic expression
  void setAnalyticSource(
    int interface_name,
    int n_cmp,
    const std::string & bcExpression      // boundary symbolic expr
  );
  /// This function sets the field of interface function (interface_id)
  void setFieldSource(
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * srcField);
  /// This function sets the field of interface function (interface_id)
  /// with an analitic expression
  void setAnalyticBoundaryValues(
    int interface_name,
    int n_cmp,
    const std::string & f);
  /// This function sets the field of interface function (interface_id)
  void setFieldBoundaryValues(
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * bcField);

   // ============================================================================  
double getAvOnBoundary_nodes(
  int                name,           // int boundary identity   (in)
  const std::string &system_name,           // system name             (in)
  int           first_cmp            // n variables       (in)
); 
// // ============================================================================
/// This function gets all the values on boundary with identity id
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(  // field (out)
    int interface_name,                  ///< boundary name (char*) (in)
    const std::string & systemName,              ///< system name           (in)
    int n_cmp,                                    ///< component             (in)
    int first_cmp=0                                     ///< component             (in)
  );
  
  
      ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnElem(
  int                  id,           ///< int boundary identity   (in)el_conn
  const std::string & system_name,           ///< system name             (in)
  int               n_cmp,           ///<  first variable       (in)
  int           first_cmp=0            ///< n variables       (in)
);
  

      
      
// ---------------------------------------------------- 
///this routin gets actual support of the interface 
 const ParaMEDMEM::MEDCouplingUMesh* getUMesh(
  int name
);  
// ---------------------------------------------------- 
///this routin gets original support of the interface 
 const ParaMEDMEM::MEDCouplingUMesh* getUMesh_orig(
  int name
); 
 
//   init interface moves the FEMUS interface
  void update_interface(
    const int interface_name, /// name of the interface
    int n_cmp, ///number of component displacement field
    const   std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField ///vector displacement field
  ) ;
  
  
    void GetInfo(
    string medfile_name,
    std::string & mesh_dir, 
    std::string & localFile, 
    std::string & filename, 
    std::vector<std::string> & meshNames, 
    std::vector<std::string> & MeshNames, 
    std::vector<std::string> & FieldNames,
    const int index_medmesh  = 0
);
#endif
// // //   end med functions

  //=============================================================================
// This function solves the problem
//   void solve();
//=============================================================================
// This function solves one step  for transient problems
  void dummy_step(
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double     &  time,                ///< actual time
    double     &  dt                   ///< step time
  ) ;
// void dummy_step(
//     const int  & t_in,                 ///< initial time iteration
//     const int  & t_step,               ///< actual time iteration
//     const int  & print_step,            ///< print every
//     double     &  time,                ///< actual time
//     double     &  dt                   ///< step time
//   ) ;
  
  
};

#endif