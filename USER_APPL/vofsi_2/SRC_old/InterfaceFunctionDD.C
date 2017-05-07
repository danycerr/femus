#include "InterfaceFunctionDD.h"


#include "Solverlib_conf.h"
#ifdef HAVE_MED

#include <iomanip>
#include <limits>
#include <boost/concept_check.hpp>

#include "VOF.h"
#include "MGMeshC.h"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"




// ======================================================
InterfaceFunctionDD::~InterfaceFunctionDD() {
  delete [] _map_vof; delete [] _map_med;
  if(_field) _field->decrRef(); if(_support_med) _support_med->decrRef();
//   delete [] _mesh_mg;
}

// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionDD::set_mesh_interface_elemID(
//   const MGMeshC & mesh_c,
//   const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const ParaMEDMEM::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id/*,*/                      ///< inrface identity  (in)
//   const int order_cmp                          ///< order pt (1 or 2) (in)
) { // ========================================================================

  // Med-mesh
   _support_med = support;                               // MED-mesh
   //   _n= _support_med->getNumberOfNodes();                 // n MED-nodes
   int n_elements_med= _support_med->getNumberOfCells(); // n MED-elements
   const ParaMEDMEM::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
   int n_nodes_med = d->getNumberOfTuples();             //  MED-nodes
   int dim_med = d->getNumberOfComponents();             //  MED dimension
   int n_nod_el=(dim_med==3)?8:4;
   
  _map_med = new int [n_elements_med];
  _map_vof =  new int [n_elements_med];
   int nxyz[3];double hxyz[3];int kxyz[3];for(int idim=0; idim<3; idim++) {
    kxyz[idim]=0;hxyz[idim]=0.;nxyz[idim]=0;
  }
  nxyz[0]=NX;  nxyz[1]=NY; /*int n_pts=1;*/
  if(dim_med==3){nxyz[2]=NZ;}
 
  for(int idim=0; idim<dim_med; idim++){hxyz[idim]=1./(nxyz[idim]);}
     
   _n=n_elements_med;
//    _n=n_nodes_med;
  double xyz_mid[3]; // MED el center
  xyz_mid[2]=0.;
  std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
  for(int idim=0; idim<dim_med; idim++)  xyz_mid[idim]=0.;
    _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    for(int inode=0; inode< n_nod_el; inode++) {
      for(int idim=0; idim<dim_med; idim++)
        xyz_mid[idim] += (d->getIJ(nodes1[inode], idim))/n_nod_el;
    } // end inode
     for(int idim=0; idim<dim_med; idim++) 
       kxyz[idim]=int((xyz_mid[idim]-0.5*hxyz[idim])/hxyz[idim]);  
     _map_med[ielem] =ielem;
     _map_vof[ielem]  = kxyz[0]+(kxyz[1])*nxyz[0]+(kxyz[2])*nxyz[0]*nxyz[1];
    
    nodes1.clear(); // clear element node vector

  }// *************************************************************************
  
  std::cout << "==================== mesh elem id  =============================== "<<  interface_id <<"\n";
  std::cout << '\n';
//   printOn(std::cout, interface_id);
  d->decrRef(); 
  return;
  
 
}



// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionDD::set_mesh_interface_nodeID(
//   const MGMeshC  &mesh_c,                ///< Femus-mesh        (in)
  const ParaMEDMEM::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id,                      ///< inrface identity  (in)
  const int order_cmp    
) { // ========================================================================
// //   // Med-mesh
   _support_med = support;                               // MED-mesh
 //   _n= _support_med->getNumberOfNodes();                 // n MED-nodes
   int n_elements_med= _support_med->getNumberOfCells(); // n MED-elements
   const ParaMEDMEM::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
   int n_nodes_med = d->getNumberOfTuples();             //  MED-nodes
   int dim_med = d->getNumberOfComponents();             //  MED dimension
//    int n_nod_el=(dim_med==3)?8:4;
   
  _map_med = new int [n_nodes_med];
  _map_vof =  new int [n_nodes_med];
   int nxyz[3];double hxyz[3];int kxyz[3];for(int idim=0; idim<3; idim++) {
    kxyz[idim]=0;hxyz[idim]=0.;nxyz[idim]=0;
  }
  nxyz[0]=NX;  nxyz[1]=NY; int n_pts=1;
  if(dim_med==3){nxyz[2]=NZ;}
 
  for(int idim=0; idim<dim_med; idim++){hxyz[idim]=1./(nxyz[idim]);n_pts *=(nxyz[idim]+1.);}
   
   
    _n=n_nodes_med;
    std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
//   for(int idim=0; idim<dim_med; idim++)  xyz_mid[ielem+idim*n_elements_med]=0.;
    _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    int n_nod_el=nodes1.size();
    for(int inode=0; inode< n_nod_el; inode++) {
      int gnode=nodes1[inode];
      for(int idim=0; idim<dim_med; idim++){
        double xyz = (d->getIJ(gnode, idim));
      kxyz[idim]=(int)((xyz+0.5*hxyz[idim])/hxyz[idim]);
      }
     _map_med[gnode] = gnode;
     _map_vof[gnode]  = kxyz[0]+kxyz[1]*(nxyz[0]+1)+kxyz[2]*(nxyz[0]+1)*(nxyz[1]+1);
    } // end inode
//      for(int idim=0; idim<dim_med; idim++) kxyz[idim]=int((xyz_mid[ielem+idim*n_elements_med])/hxyz[idim]);  
   
    
    nodes1.clear(); // clear element node vector

  }// *************************************************************************


  
  std::cout << "==================== mesh elem id  =============================== "<<  interface_id <<"\n";
  std::cout << '\n';
//   printOn(std::cout, interface_id); 
//   d->decrRef(); 
  return;
}



// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void InterfaceFunctionDD::set_analytic_field(
  const char *symbolic_eq,     // symbolic function
  int nComp             // number of componenents
) { // ========================================================================
  if(_field) _field->decrRef();
  //ON CELLS does not work (due to baryc in quad 9 or hex 27)
  ParaMEDMEM::TypeOfField type = ParaMEDMEM::ON_NODES;
  int dim=_support_med->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if(dim > 0)    vars[0] = "x";
  if(dim > 1)    vars[1] = "y";
  if(dim > 2)    vars[2] = "z";
  if(symbolic_eq==std::string("")) {
    std::cout<<
             "InterfaceFunctionDD::set_analytic_field: NULL function"; abort();
  }
  _field = _support_med->fillFromAnalytic3(type, nComp, vars, symbolic_eq);
  _field->setName(symbolic_eq);
  _field->checkCoherency();
//   std::ostream out("prova.field");
  std::cout << "InterfaceFunctionDD::set_analytic_field_interface \n";
//   printOn(std::cout);
  return;
}



// ==========================================================================
void InterfaceFunctionDD::set_field(
  const ParaMEDMEM::MEDCouplingFieldDouble *f
) {// =======================================================================
  if(_field) _field->decrRef();
//
  ParaMEDMEM::TypeOfField type = f->getTypeOfField();
  _field = ParaMEDMEM::MEDCouplingFieldDouble::New(type);
  *_field=*f;

//   _field->setName(f->getName());
  _field->checkCoherency();

  std::cout << "InterfaceFunctionDD::set_field \n";

  return;
}

// ======================================================================================
void InterfaceFunctionDD::get_val_field(
  int node_med,     ///< Med node index   (in)
  int n_cmp,        ///< number of components (vec val[]) (in)
  double val[]      ///< values (out)
) { // ==================================================================================
  if(_field == NULL) {for(int i=0; i<n_cmp; i++) val[i] =0.0;}
  else {
    for(int i = 0; i<n_cmp; i++) val[i] = _field->getIJ(node_med, i);
  }
}

// =============================================================
// This function prints the interface function
void InterfaceFunctionDD::printOn(
  std::ostream & out,                     ///< ostream file
  int id                                  ///< interface name
) const { // ===================================================

  // set up title --------------------------------------------------------
  out << std::endl
      << "================================================== \n"
      << std::endl
      << " InterfaceFunctionDD =" << id<<  std::endl << std::endl;
  // ---------------------------------------------------------------------
  // ---------------- Mapping ---------------------------------------------
  //  node map MGMesh <-> MEDCoupling ------------------------------------
  out << "node/element correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<int, int>::const_iterator itN = _nodeID.begin();
  for(int i_mg=0; i_mg < _n; i_mg++) {
    out << "  " << std::setw(4) << _map_vof[i_mg]
        << " <-> " << std::setw(4) << _map_med[i_mg] << std::endl;
  }
  out << std::endl;

//   // element map MGMesh <-> MEDCoupling ----------------------------------
//   out << "element correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<int, int>::const_iterator itE = _elemID.begin();
//   for(; itE != _elemID.end(); itE++)
//     out << "  " << std::setw(4) << itE->first
//         << " <-> " << std::setw(4) << itE->second << std::endl;
//   out << std::endl;
  // face map    MGMesh <-> MEDCoupling ----------------------------------
//   out << "face correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<std::pair<int,int>, int>::const_iterator itF = _faceID.begin();
//   for(; itF != _faceID.end(); itF++)
//     out << " (" << std::setw(4) << itF->first.first
//         << ", " << std::setw(4) << itF->first.second << ")"
//         << " <-> " << std::setw(4) << itF->second << std::endl;
//   out << std::endl;

  // ---------------------------------------------------------------------
  // ----------------  field ---------------------------------------------
  if(!_field) {
    out << "field= NULL" << std::endl;
    out << "==================================================" << std::endl;
    return;
  }
  ParaMEDMEM::TypeOfField type = _field->getTypeOfField();
  const ParaMEDMEM::DataArrayDouble * d;
  int n;
  const ParaMEDMEM::MEDCouplingUMesh * m
    = dynamic_cast<const ParaMEDMEM::MEDCouplingUMesh *>(_field->getMesh());

  out << "name : " << _field->getName() << std::endl;

    int ncf=_field->getNumberOfComponents();
    int nc;
     ParaMEDMEM::DataArrayDouble * v = _field->getArray();
  
  switch(type) {
  case ParaMEDMEM::ON_NODES:
    out << "type : ON_NODES" << std::endl;
    d = m->getCoords();
    d->incrRef();
    n = d->getNumberOfTuples();
    out << " " <<n << " nodes";
    
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;

    nc = d->getNumberOfComponents();
  for(int i=0; i<n; i++) {
    for(int j=0; j<nc; j++) out << (j == 0 ? "f( " : ", ")
                                  << std::setw(10) << std::fixed << d->getIJ(i,j);
    out << ") = ";
    for(int j=0; j<ncf; j++)  out << (j == 0 ? "( " : ", ")
                                    << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << ")  ";
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;
    break;
    
    
  case ParaMEDMEM::ON_CELLS:
    out << "type : ON_CELLS" << std::endl;
//     d = m->getBarycenterAndOwner();
//     n = d->getNumberOfTuples(); 
    nc =m->getNumberOfCells();
    out << " " << nc << " cells";
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;
  
  for(int i=0; i<nc; i++) {
     out << "f( "
                                  << std::setw(4) << std::fixed << i;
    out << ") = ";
    for(int j=0; j<ncf; j++)  out << (j == 0 ? "( " : ", ")
                                    << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << ")  ";
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;
    break;
  default :
    std::cout<< "Error type"; abort();
  }

 

//   v->decrRef();
}

#endif



