#ifndef __VOF__
#define __VOF__

#include <vector>
#include "Solverlib_conf.h"
#include <mpi.h>

class MGFemusInit;
class MGUtils;
class MGSolCC;
class MGMeshC;

#ifdef HAVE_MED
namespace ParaMEDMEM {
class MEDLoader;
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}
class InterfaceFunctionDD;
#endif


class VOF {
  // ==========================================================================
  //=========   data  =========================================================
protected:

  // data communication (defined in Constructor) -------------------------------
  MGFemusInit *    _start;          // start function
  MPI_Comm         _comm;           // communicator
  bool             _local_MPI_Init; // initial mpi flag

  // code data ------------------------------------------------------------------
  MGUtils *     _mg_utils;     // param and files
  MGMeshC *     _mg_mesh;      // FEMus-mesh
 
 
  MGSolCC *     _mgcc;         // system
  // interface data
#ifdef HAVE_MED
  ParaMEDMEM::MEDCouplingUMesh *    _med_mesh;     // Med-mesh
  std::map<int,InterfaceFunctionDD *> _interfaceFunMap; // MG-Med interface map
#endif

  // ==========================================================================
  //=========  public functions  ============================================
public:
  
  // Constructor-Destructor functions =================================================== 
  VOF();                          ///< Empty constructor
  VOF(MPI_Comm comm);             ///< Constructor with communicator
  ~VOF();                         ///< Destructor
  void terminate();               ///< Destructor
   
  // Set functions  =====================================================================  
  void set_param(MGUtils &mgutils);  ///< set parameters
  void setMesh();                    ///< set mesh
  void setSystem(int name );                 ///< set system
  void setFieldSource_Vinit(const int Level, const double dt);  ///< EXternal vel field
  
  // Get functions  =====================================================================
  const MGMeshC & get_MGMesh() {return *_mg_mesh;};
  MGSolCC & get_MGSystem() {return *_mgcc;};
  
// Solve functions  =====================================================================
// -------------------------------------------------------------------------------------
/// This function sets up the intial set
  void solve_setup(
    int         t_in,                 ///< initial time iteration
    double       time                 ///< actual time
  );
// --------------------------------------------------------------------------------------
// This function solves one step  for transient problems
  void solve_onestep(
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double     &  time,                ///< actual time
    double     &  dt                   ///< step time
  ) ;
  
 
// =================  with MED lib ===============================================
#ifdef HAVE_MED
/// 2. This function retrun the pointer of the InterfaceFunctionC from the id name
  InterfaceFunctionDD * get_interface_fun(
    int id                                   ///< interface identity
  ) {
    return(_interfaceFunMap.find(id)==_interfaceFunMap.end()?NULL :_interfaceFunMap[id]);
  }
// ======================================================================================
  // ======================================================================================
/// 4) This function gets the value of color function
  ParaMEDMEM::MEDCouplingFieldDouble * getValues_C_cell();

  void init_interface(
    const int interface_name,
    int interface_id,
    const std::string & medfile_name // medfile name    (in)
  ) ;
   void init_interface(
    const int interface_name,
    int interface_id_1,
    int interface_id_2,
    const std::string & medfile_name // medfile name    (in)
  ) ;
    void init_interface_oncell(
    const int interface_name,
    int interface_id,
    const std::string & medfile_name // medfile name    (in)
  ) ;
  
  /// This function sets the field of interface function (interface_id)
  void setFieldSource(
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * srcField);
  
    /// This function sets the field of interface function (interface_id)
  void setFieldSource_disp(
    int interface_name,
    int n_cmp,
    const   std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField
//     const ParaMEDMEM::MEDCouplingFieldDouble * srcField
  );

 /// This function get the mesh of a selectd interface 
 const  ParaMEDMEM::MEDCouplingUMesh* getUMesh(
  int name
);  
  
#endif // ===============================  med library ============================



};

#endif
