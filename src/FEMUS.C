#include <iostream>
#include <cstdlib>
#include <sstream>


// configuration files -------------------------
#include   "Printinfo_conf.h"

// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "Equations_conf.h"
#include "MGTimeLoop.h"


// class include
#include "FEMUS.h"

#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionM.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif

#include "MeshExtended.h"
#include "EquationSystemsExtendedM.h"


// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
FEMUS::FEMUS()  :
  _comm(MPI_COMM_WORLD) {  // communicator
  // Init MPI flag --------------------------------
  int flag=0;  MPI_Initialized(&flag);
  if(flag) {_local_MPI_Init = false; }
  else {_local_MPI_Init = true; }
  // femus init -----------------------------------
  int argc = 0;    char ** argv = NULL;
  _start=new  MGFemusInit(argc,argv,_comm);

  return;
}

// ============================================================================
// This function is a constructor with  communicator
FEMUS::FEMUS(
  MPI_Comm comm
):
  _comm(comm)   // communicator
//   _num_mgmesh(0)
{
  // n of femus_meshes
  // transient system
  // Init MPI flag
  int flag=0;  MPI_Initialized(&flag);
  if(flag) {_local_MPI_Init = false; }
  else {_local_MPI_Init = true; }

  // femus init
  int argc = 0;    char ** argv = NULL;
  _start= new MGFemusInit(argc,argv);

  return;
}
// =======================================================================
void FEMUS::init_param(
  MGUtils &   mgutils,
  int name
) { // ====================================================================
  _mg_utils=&mgutils;
  _mg_utils->set_name(name);
  return;
}
// ============================================================================
void FEMUS::init_fem(
  MGGeomEl & mggeomel,
  MGFEMap & mgfemap
) { // ========================================================================
// A) setting MGGeomEl
  _mg_geomel=&mggeomel;  // ***************************************************
  if(_mg_geomel == NULL) {
    std::cout<< "FEMUS::init_fem: no _mg_geomel"; abort();
  }
  /// B) setting MGFEMap (fem)
  _mg_femap=&mgfemap;  // *****************************************************
  if(_mg_femap == NULL) {
    std::cout<< "FEMUS::init_fem: no _mg_femap"; abort();
  }

  return;
}
// ============================================================================
// This function is the destructor
FEMUS::~FEMUS() {
  // ==========================================================================
  delete _start;
  delete _mg_time_loop;
//   delete _mg_equations_map;
//   delete _mg_utils;
//   delete _mg_mesh;
//   delete _mg_femap;
#ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif

}

// ============================================================================
// This function is the problem destructor
void FEMUS::terminate(
) {// =========================================================================

}

// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************
#ifdef   TWO_PHASE
void FEMUS::set_mgcc(MGSolCC & cc
  
) {
  _mg_equations_map->set_mgcc(cc);
  return;
}
#endif

// // // =============================================================================
// // // This function sets the type of problem
void FEMUS::setSystem(
  const std::vector<FIELDS> & pbName,
  int n_data_points,
  int n_data_cell
) {// ==========================================================================


  _mg_equations_map=new EquationSystemsExtendedM(*_mg_utils,*_mg_mesh,*_mg_femap,n_data_points,n_data_cell);  // MGEquationsMap class
  _mg_equations_map->read_par();
#ifdef PRINT_INFO  // ---- info ---------------
  _mg_equations_map->print_par();       // print parameters
#endif
  _mg_equations_map->init_data(0);
  _mg_equations_map->init(pbName);                              // adds the equations to the map
  _mg_equations_map->setDofBcOpIc();                            // set operators
  _mg_equations_map->set_mesh_mg(*_mg_mesh);
#ifdef HAVE_MED
  _mg_equations_map->set_mesh_med(*_med_mesh);
#endif
//   }
  if(_mg_geomel == NULL) {
    std::cout<< "FEMUS::setSystem: no _mg_equations_map"; abort();
  }

  //time loop
  _mg_time_loop=new  MGTimeLoop(*_mg_utils,*_mg_equations_map);
  if(_mg_time_loop == NULL) {
    std::cout<< "FEMUS::setSystem: no _mg_time_loop"; abort();
  }

  return;
}


// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::setMesh(
) {// ==========================================================================

  const int NoLevels= _mg_utils->get_par("nolevels");  // numb of Level
  _mg_mesh=new MeshExtended(_start->comm(), *_mg_utils,*_mg_geomel);
  // check insanity
  if(_mg_mesh == NULL) {
    std::cout<< "FEMUS::setMesh: no _mg_mesh"; abort();
  }
  if(NoLevels != _mg_mesh->_NoLevels) {
    std::cout << "Inconsistent Number of Levels between Mesh and SolBase"
              << std::endl;
    abort();
  }
  // print mesh at level NoLevels-1 (linear connectivity)
  _mg_mesh->print(NoLevels-1,0);

#ifdef HAVE_MED
  // prind mesh at level NoLevels-1 (med format)
  std::string mesh_name= _mg_utils->get_file("F_MESH_READ");//://= _mg_utils.get_file("F_MESH_READ");
  unsigned pos = mesh_name.find(".");         // position of "live" in str
  std::ostringstream name;
  if(_mg_mesh->_iproc==0) {
    name << _mg_utils->_mesh_dir <<  mesh_name.substr(0,pos) << "_fine.med" ;
    _mg_mesh->print_med(NoLevels-1,name.str().c_str());
  }
#endif

  return;
}

// // // *******************************************************************
// // // *******************************************************************
 int  FEMUS::get_proc() const{
  return _mg_mesh->_iproc;
}
// // // *******************************************************************


// // // *******************************************************************
// // // **************** Solve  *******************************************
// // // *******************************************************************
/// This function sets up the intial set
void FEMUS::solve_setup(
  int        & t_in,                 ///< initial time iteration
  double     &  time                 ///< actual time
) {
  const int restart      = _mg_utils->get_par("restart");   // restart or not
  _mg_time_loop->transient_setup(restart,t_in,time);     //  MGTimeLoop: setup
  return;
}

//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_onestep(
  const int  & t_in,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt,                   ///< step time
   const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) { // ========================================================================
  _mg_time_loop->transient_onestep(t_in,t_step,print_step,time,dt,eq_min,eq_max);    ///< step time
  return;
}
// This function solves one step  for transient problems
void  FEMUS::solve_steady(
  const int & nmax_step,  ///< number max of steps
  const double & toll,  ///< tolerance
  const int  & it_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  dt,         ///< inial time step
  const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int     &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) { // ========================================================================
  _mg_time_loop->steady(nmax_step,toll,it_step,print_step,dt,eq_min,eq_max);    ///< step time
  return;
}

// This function solves one step  for transient problems
void  FEMUS::set_uooold(
   const int & flag,  ///<  0 xold-> x_ooold   1 x_ooold-> xold
  const double & toll,  ///< tolerance
  const double delta_t_step_in,  //   (in)
  const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int     &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
){ // ========================================================================
  _mg_time_loop->set_uooold(flag, toll,delta_t_step_in,eq_min,eq_max);    ///< step time
  return;
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep(
  const int  & t_in,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt                   ///< step time
) { // ========================================================================
  _mg_time_loop->transient_control_onestep(t_in,t_step,print_step,time,dt);    ///< step time
  return;
}


double  FEMUS::System_functional(
  const int  & ff,                 ///< initial time iteration
    double      parameter,             ///< functional parameter
   double     & control                   ///< step control
) {
  return _mg_equations_map->System_functional(ff,parameter,control);
}

#ifdef HAVE_MED
void FEMUS::init_interface(
  const int interface_name,
  int interface_id_1,
  int interface_id_2,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  bool on_nodes,
  const int index_medmesh                       // med-mesh index  (in)
) {// =========================================================================
  ostringstream name_id_1; name_id_1 <<interface_id_1;
  ostringstream name_id_2; name_id_2 <<interface_id_2;
  std::vector<std::string> vG(2);
//   std::string id_name=name_id.str().c_str();
  vG[0] =name_id_1.str().c_str();
  vG[1] =name_id_2.str().c_str();

  // Reading mesh Names from med file  ------------------------------

  std::string mesh_dir=_mg_utils->_mesh_dir;
  std::string localFile=mesh_dir+medfile_name;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[index_medmesh];
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());


  int nGroups0 = GroupNames.size();
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id_1<10) {id_level=0;}
  support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
  support->zipCoords();

  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<
            interface_id_1  << "\n";
  _mg_equations_map->add_interface_fun(interface_name, interface_id_1, support,on_nodes,order_cmp);


  return;
}




// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
void FEMUS::init_interface(
  const int interface_name,
  int interface_id,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  bool on_nodes,
  const int index_medmesh                       // med-mesh index  (in)
) {// =========================================================================

  ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();
  
  std::string mesh_dir, localFile, filename;
  std::vector<std::string> meshNames, MeshNames, FieldNames;
  
  GetInfo(medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames);
  
  ParaMEDMEM::MEDCouplingFieldDouble* acc;
  ParaMEDMEM::DataArrayDouble* ACC;
      
  if(on_nodes) acc = MEDLoader::ReadField(ParaMEDMEM::ON_NODES, filename.c_str(), MeshNames[0], 0, FieldNames[0], -1,-1);  
  else acc = MEDLoader::ReadField(ParaMEDMEM::ON_CELLS, filename.c_str(), MeshNames[0], 0, FieldNames[0], -1,-1);  
  
  ACC = acc->getArray();

  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id<10) {id_level=0;}
  
  support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
  
  init_interface(interface_name, order_cmp, support, ACC, interface_id);

//   _mg_equations_map->add_interface_fun(interface_name, interface_id, support,on_nodes,order_cmp);

  return;
}

void FEMUS::init_interface(
  const int interface_name,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  bool on_nodes
) {// =========================================================================

  // Reading mesh Names from med file  ------------------------------

  std::string mesh_dir, localFile, filename;
  std::vector<std::string> meshNames, MeshNames, FieldNames;
  
  GetInfo(medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames);
  
  ParaMEDMEM::MEDCouplingFieldDouble* acc;
  ParaMEDMEM::DataArrayDouble* ACC;
  acc = MEDLoader::ReadField(ParaMEDMEM::ON_NODES, filename.c_str(),MeshNames[0],0,FieldNames[0],-1,-1);  
  ACC = acc->getArray();
  
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=0;
  
  support = MEDLoader::ReadUMeshFromFile(localFile.c_str(), meshNames[0].c_str(), id_level);
  
  init_interface(interface_name, order_cmp, support, ACC);

  return;
}

void FEMUS::init_interface(
  const int interface_name,
  int order_cmp,
  ParaMEDMEM::MEDCouplingUMesh * support,
  ParaMEDMEM::DataArrayDouble* ACC,
  int interface_id
) {
 
  std::map<int, int> MedToMg;

  int celle = support->getNumberOfCells();
  int nodes = support->getNumberOfNodes();
  int NodesPerCell = support->getNumberOfNodesInCell(0);

  ParaMEDMEM::DataArrayInt* MappaNodi = ParaMEDMEM::DataArrayInt::New();
  MappaNodi->alloc(nodes,2);
  
  std::vector<int> MgNodesIds;
  std::vector<int> MedNodesIds;
  
//   int NodesPerCell;
  for(int i =0; i<celle; i++){
    support->getNodeIdsOfCell(i,MgNodesIds);
    for(int j = 0; j<NodesPerCell; j++) MappaNodi->setIJ(MgNodesIds[i*NodesPerCell +j],0,MgNodesIds[i*NodesPerCell +j]);
  }
  
  support->zipCoords();
  
  for(int i =0; i<celle; i++){
    support->getNodeIdsOfCell(i,MedNodesIds);
    for(int j = 0; j<NodesPerCell; j++) MappaNodi->setIJ(MgNodesIds[i*NodesPerCell +j],1,MedNodesIds[i*NodesPerCell +j]);
  }
  
  for(int i=0; i<celle*NodesPerCell; i++){
    double node = ACC->getIJ(MgNodesIds[i],0);
    MedToMg[MedNodesIds[i]] = (int)ACC->getIJ(MgNodesIds[i],0);
  }
  
  sort(MedNodesIds.begin(), MedNodesIds.end());
  vector<int>::iterator it;
  it = unique(MedNodesIds.begin(), MedNodesIds.end());  

  MedNodesIds.resize(distance(MedNodesIds.begin(),it)); 

  const int NODI = MedNodesIds.size();
  int* map_med = new int[NODI];
  int* map_mg  = new int[NODI];

  for(int i = 0; i<MedNodesIds.size(); i++){
    map_med[i] = MedNodesIds[i];
    map_mg[i]  = MedToMg[MedNodesIds[i]];
  }
  
  
  MedNodesIds.clear();
  MgNodesIds.clear();
  
  if(interface_id==0) std::cout << "\n \033[1;31m Interface "<<interface_name
                                <<" support set to mesh without volume group and order "<<order_cmp 
                                << "\033[0m\n";
  else std::cout << "\n \033[1;31m Interface "<<interface_name
                 <<" support set to mesh with interface " <<interface_id 
                 <<" and order "<<order_cmp 
                 << "\033[0m\n";

  InterfaceFunctionM *fun = new InterfaceFunctionM;
  fun->set_maps(map_med,map_mg, NODI);
  fun->set_order(order_cmp);
  fun->set_mg_mesh(_mg_mesh);
  fun->set_support_med(support);
  fun->set_NumberOfNodes(NODI);
  _mg_equations_map->add_interface_fun(interface_name, fun);
  
  return;
}

void FEMUS::init_interface(
  const int interface_name,
  int interface_id,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  const std::string & medfile_name_old, // medfile name    (in)
  const FEMUS & P_old,
  const int index_mgmesh,           // mgmesh index    (in)
  const int index_medmesh                      // med-mesh index  (in)

) {// =========================================================================


  ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();

  // Reading mesh Names from med file  ------------------------------
  std::string mesh_dir=_mg_utils->_mesh_dir;
  std::string localFile=mesh_dir+medfile_name;
  std::string localFile_old=mesh_dir+medfile_name_old;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[index_medmesh];
  
  // Reading group names
  std::vector<std::string> GroupNames =  MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());

  int nGroups0 = GroupNames.size();

  // From group names to interface_mesh (id,name,support)
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1; bool on_nodes=true;
  if(interface_id<10) {id_level=0; on_nodes=false;}
  
  
  
  support= MEDLoader::ReadUMeshFromGroups(localFile_old.c_str(), meshNames[index_medmesh].c_str(), id_level, vG);
  support->zipCoords();

  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<id_name  <<"\n";
  _mg_equations_map->add_interface_fun(interface_name,interface_id, support,on_nodes,order_cmp);

  return;
}
// =========================================================================
/// This function sets the value from first_cmp to end_cmp of the field
/// on the old solution x_old of the interface id
void FEMUS::write_Boundary_value(
  int id_boundary_name ,      ///< identity interface name (in)
  std::string mgsystem_name, ///< system name          (in)
  int n_cmp,             ///< from variable system (in)
  int first_cmp               ///< to variable system   (in)

) {
  _mg_equations_map->write_Boundary_value(
    id_boundary_name,
    mgsystem_name,n_cmp,first_cmp);
  return;
}


// =============================================================================
//This function reads and sets the med-mesh
// (also the libmesh calling the other setMesh function)
void FEMUS::setMedMesh(const std::string & dataFile) {

  // Reading MED mesh ----------------------------------------------------
  // Mesh  filename
  std::string localFile;
  int l = dataFile.size();
  if(dataFile.substr(l-4) == ".med")    {localFile = dataFile;}
  else {
    std::cout<<"FEMUS::setMesh:"<<  dataFile <<"does not exist!";
  }
  std::cout << " FEMUS::setMesh: MED file "  << ": "<< localFile << std::endl;

  // Mesh names (inside MEDfile one can have different meshes
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
//   std::vector<std::string> meshNames = MEDLoader::GetMeshNames("/homesd/msandro/software/femus/USER_APPL/MESH/test1quad9_group_mat_gen.med");
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'"; abort();
  }
  // The first name is the good one
  std::string localMeshName = meshNames[0];
  _med_mesh = MEDLoader::ReadUMeshFromFile(localFile.c_str(), localMeshName.c_str(), 0);
  if(_med_mesh == NULL) {
    std::cout<<  " FEMUS::setMesh : unable to read the med-mesh'"; abort();
  }

  return;
}



// ===================================================================
void FEMUS::setAnalyticBoundaryValues(
  int interface_name,
  int n_cmp,
  const std::string & bcExpression      // boundary symbolic expr
) {
  _mg_equations_map->setBC(interface_name,n_cmp , bcExpression.c_str());
  return;
}

// ===================================================================
void FEMUS::setAnalyticSource(
//   const std::string & bcName,           // boundary name
  int interface_name,
  int n_cmp,
  const std::string & bcExpression      // boundary symbolic expr
) {
  _mg_equations_map->setBC(interface_name,n_cmp , bcExpression.c_str());
  return;
}

void FEMUS::setFieldBoundaryValues(
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * bcField) {

  if(bcField ==NULL) { return; }

  _mg_equations_map->setBC(interface_name,n_cmp, bcField);

  return;
}

void FEMUS::setFieldSource(
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * srcField) {

  if(srcField ==NULL) { return; }
  _mg_equations_map->setBC(interface_name,n_cmp, srcField);

  return;
}
// =======================================================================
double FEMUS::getAvOnBoundary_nodes(
  int                name,           // int boundary identity   (in)
  const std::string & system_name,           // system name             (in)
  int           first_cmp            // n variables       (in)
)  {
  return _mg_equations_map->getAvOnBoundary_nodes(name,system_name.c_str(),first_cmp);
}
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnElem(
  int                  id,           ///< int boundary identity   (in)el_conn
  const std::string & system_name,           ///< system name             (in)
  int               n_cmp,           ///<  first variable       (in)
  int           first_cmp            ///< n variables       (in)
)  {// ========================================================================

  return _mg_equations_map->getValuesOnElem(id,system_name.c_str(), n_cmp,first_cmp);
}

// ============================================================================
/// This function gets all the values on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary(  // field (out)
  int  interface_name,                     // boundary name (char*) (in)
  const std::string & systemName,                  // system name           (in)
  int n_cmp,                                        // component             (in)
  int first_cmp                                        // component             (in)
)  {
  return _mg_equations_map->getValuesOnBoundary_nodes(interface_name,systemName.c_str(),n_cmp,first_cmp);
}

// ============================================================================
/// This function gets the actual mesh of FEMUS problem
 const ParaMEDMEM::MEDCouplingUMesh* FEMUS::getUMesh(
  int name
){
 return _mg_equations_map->getUMeshCoupling(name);    
}
// // // ===========================================================================
/// This function gets the original mesh of FEMUS problem
 const ParaMEDMEM::MEDCouplingUMesh* FEMUS::getUMesh_orig(
  int name
){
 return _mg_equations_map->getUMeshCoupling_orig(name);    
}
// ------------------------------------------------
/// This function moves the FEMUS interface according to a given displacement field 
void FEMUS::update_interface(
  const int interface_name,
  int n_cmp,
  const   std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField
) {// =========================================================================
    ParaMEDMEM::DataArrayDouble *  mg_disp= ParaMEDMEM::DataArrayDouble::Meld(srcField[0]->getArray(),
                                                                              srcField[1]->getArray());
 if (n_cmp==3)   mg_disp= ParaMEDMEM::DataArrayDouble::Meld(mg_disp,srcField[2]->getArray());
    std::cout<< "FEMUS::UPDATE: Tuples  src_field  "<< mg_disp ->getNumberOfTuples() << 
    "  Comp    " << mg_disp ->getNumberOfComponents()
    << " and x value is " <<  mg_disp->getIJ(0,0)<< " and y value is " <<  mg_disp->getIJ(0,1)
    <<std::endl;
  std::cout << "FEMUS::UPDATE: support  set to boundary with name "<<
            interface_name  << "\n";
  
  
    ParaMEDMEM::MEDCouplingUMesh * support;    ParaMEDMEM::MEDCouplingUMesh * support_up;
//   int id_level=-1;  if(interface_id_1<10) {id_level=0;}
//   support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
  support=(getUMesh_orig(interface_name))->clone(1);
   support_up=(getUMesh(interface_name))->clone(1);
 // support->zipCoords();

   ParaMEDMEM::DataArrayDouble *coord; ParaMEDMEM::DataArrayDouble *coord_up;
  
    coord= support->getCoords();
 coord_up= support_up->getCoords();
  int npt=coord->getNumberOfTuples();
  int ncomp= coord->getNumberOfComponents();

  ParaMEDMEM::DataArrayDouble *new_cord=ParaMEDMEM::DataArrayDouble::Add(coord,mg_disp);
  
  std::cout<< "Tuples    "<< npt << "Comp    " << ncomp<<std::endl;
  support->setCoords(new_cord);
  
 InterfaceFunctionM * fct = _mg_equations_map->get_interface_fun(interface_name);
 fct->update_support(support);


  return;
}

void FEMUS::GetInfo(
    string medfile_name,
    std::string & mesh_dir, 
    std::string & localFile, 
    std::string & filename, 
    std::vector<std::string> & meshNames, 
    std::vector<std::string> & MeshNames, 
    std::vector<std::string> & FieldNames,
    const int index_medmesh   
){
  mesh_dir=_mg_utils->_mesh_dir;
  localFile=mesh_dir+medfile_name;
    
  meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
//   std::string MMMM = medfile_name;
//   int posP = MMMM.find("_");  // position of "live" in str
//   
//   filename =   mesh_dir+MMMM.substr(0,posP)  + "_MedToMg.med" ;
//   std::cout<<"\n 1:"<<localFile<<"\n 2:"<<filename<<std::endl;
  filename =localFile;
  
  FieldNames = MEDLoader::GetAllFieldNames(filename.c_str());
  MeshNames = MEDLoader::GetMeshNames(filename.c_str());
 
  return;
};

#endif
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::dummy_step(
  const int  & t_in,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt                   ///< step time
) { // ========================================================================
  _mg_time_loop->dummy_step(t_in,t_step,print_step,time,dt);    ///< step time
  return;
}


// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
