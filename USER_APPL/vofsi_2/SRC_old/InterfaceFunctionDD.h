#ifndef __INTERFACE_FUNCTION_VOF__
#define __INTERFACE_FUNCTION_VOF__

#include "Solverlib_conf.h"

#include <iostream>
#include <map>
#include <vector>


#ifdef HAVE_MED
namespace ParaMEDMEM {
class MEDCouplingFieldDouble;
class MEDCouplingUMesh;
}

// class MGMeshC;
//

// class MeshExtended;

class InterfaceFunctionDD {

protected:
  // data interface function ==================================================
  const ParaMEDMEM::MEDCouplingUMesh * _support_med;  ///< MEDCoupling mesh
  ParaMEDMEM::MEDCouplingFieldDouble * _field;        ///< MEDCoupling field
  int _n;          ///< number of nodes/elements in the interface
  int *_map_med;   ///< node/element map index->MGMesh  (size: _n)
  int *_map_vof;    ///< node/element map index->MEDCoupling (size: _n)

public:
  // Constructors- Destructors ==========================================================
  //
  InterfaceFunctionDD() :   ///< Constructor --------------------------------------------
    _support_med(NULL),     ///< med-sub/mesh
    _field(NULL),           ///< field
    _map_med(NULL),         ///< MED map
    _map_vof(NULL)          ///< MED vof
  {         }
  ~InterfaceFunctionDD();       ///< Destructor -----------------------------------------
  
  // set functions ======================================================================
  // ------------------------------------------------------------------------------------
  /// Setting the field (analytic)
  void set_analytic_field(
    const char *code,          ///< symbolic expression  (in)
    int nComp                  ///< n of components      (in)
  );
   // -----------------------------------------------------------------------------------
  /// Setting the field (MEDfield)
  void set_field(
    const ParaMEDMEM::MEDCouplingFieldDouble *f   ///< field  (in)
  ); 
   // -----------------------------------------------------------------------------------
//   /// Setting the field (MEDfield)
//   void set_field_source(
//     const ParaMEDMEM::MEDCouplingFieldDouble *f ///< field  (in)
//     );
 // ------------------------------------------------------------------------------------- 
 /// This function computes the mesh interface
 /// to set in the interface-function. The field in the
 /// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_nodeID(
    const ParaMEDMEM::MEDCouplingUMesh * support,///< med-mesh          (in)
    const int interface_id,                      ///< inrface identity  (in)
    const int order_cmp                          ///< order pt (1 or 2) (in)
  );
 // ------------------------------------------------------------------------------------- 
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_elemID(
    const ParaMEDMEM::MEDCouplingUMesh * support,///< med-mesh          (in)
    const int interface_id                      ///< inrface identity  (in)
  );

  // get functions ======================================================================
  // ------------------------------------------------------------------------------------
  int get_n() { return _n; } 
  // ------------------------------------------------------------------------------------
  int * get_map_vof() { if(_map_vof==NULL) return NULL; return  _map_vof; } 
  // ------------------------------------------------------------------------------------
  int * get_map_med() { if(_map_med==NULL)  return NULL; return _map_med; } 
  // ------------------------------------------------------------------------------------
  const ParaMEDMEM::MEDCouplingUMesh * getSupport() {return _support_med;}
  // ------------------------------------------------------------------------------------
  /// This function evaluates the field
  void get_val_field(
    int node_med,     ///< Med node index   (in)
    int n_cmp,        ///< number of components (vec val[]) (in)
    double val[]      ///< values (out)
  );
// ===========================================================================
  // This function prints the interface function
  void printOn(
    std::ostream & out,                     ///< ostream file
    int id                                  ///< interface name
  ) const;
};


#endif
#endif
