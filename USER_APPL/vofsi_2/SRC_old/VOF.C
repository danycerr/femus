// std libraries ------------------------------------------
#include <cstdlib>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::ostringstream
#include <iomanip>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSolverCC.h"

// class include
#include "VOF.h"



#ifdef HAVE_PETSCM // Petsc
#include "petsc.h" // for Petsc solver
#endif
#ifdef HAVE_MPI    // Mpi
#include <mpi.h>   //For MPI_COMM_WORLD
#endif
#ifdef HAVE_MED    // MED includes

#include "InterfaceFunctionDD.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif

// #include "MeshExtended.h"
#include "MGMeshC.h"
#include "MGSolverCC.h"
#include "vof_config.h"


// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
VOF::VOF()  :
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
VOF::VOF(
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
  _start=new  MGFemusInit(argc,argv);

  return;
}
// =======================================================================
void VOF::set_param(
  MGUtils &   mgutils
) { // ====================================================================
  _mg_utils=&mgutils;
  return;
}

// This function is the destructor
VOF::~VOF() {
  // ==========================================================================
  delete _start;
//   delete _mg_mesh;
  delete _mgcc;
#ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif

}

// ============================================================================
// This function is the problem destructor
void VOF::terminate(
) {// =========================================================================

}
// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************


// // // ****************************************************************************
// // // ****************    Set    *************************************************
// // // ================================================================================
// This function sets the type of problem
void VOF::setSystem(
//   const std::vector<FIELDS> & pbName,
//   int n_data_points,
//   int n_data_cell
  int name
) {// ===================================================================================
  int Nlev_cc=CCLEV; // number of matrix cc levels  (REFLEV)
//     int    n_steps = mgutils[0]->get_par("nsteps");
//   double      dt = _mg_utils->get_par("dt");
//     int print_step = mgutils[0]->get_par("printstep");
//     int    itime_0  = mgutils[0]->get_par("itime");
  double init_time    = 0.;
//     int nolevels = mgutils[0]->get_par("nolevels");
  int restart= _mg_utils->get_par("restart");

//    int dt=0.01; int restart=0;int init_time=0;
  // Multigrid ------------------------------------------------------------------------
//     MGSol *mgs=NULL;   MGSolT *mgsT=NULL;
  // ****************************************************

  printf("\n Vof Initializing:: with level %d \n",Nlev_cc); //-----------------------
 int ne_xyz[DIMENSION];ne_xyz[0]=NX;ne_xyz[1]=NY;
   if(DIMENSION==3) ne_xyz[2]=NZ;
  _mgcc=new MGSolCC(Nlev_cc,name); // set the multilevel cc structure dim
  _mgcc->init(/**_mg_mesh,0,*/ne_xyz);
  _mgcc->GenSol(Nlev_cc-1);
//   _mgcc=new MGSolCC(Nlev_cc); // set the multilevel cc structure dim
//   _mgcc->init(/**_mg_mesh,0,*/NX,NY,NZ);
//   _mgcc->GenSol(Nlev_cc-1);

  if(restart !=0) {
    init_time=restart;
//      _mgcc->read_fine_hdf5(restart,Nlev_cc-1);
//     mcase->read(*mgmesh1,*mgs,*mgcc,*mgsT,init_time,NoLevels-1);
  }

  return;
}
// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void VOF::setMesh(
) {// ==========================================================================
  // generation of Cartesian Mesh (0,bx)x(0,by)x(0,bz) with (NX,NY,NZ) interval ------
  fprintf(stderr," Mesh:");
  int mesh_nlevels=1;// number of mesh levels (only level=0 for the mesh)
//     MGMeshC *mgmesh1; mgmesh1=new MGMeshC(DIMENSION,mesh_nlevels);
  _mg_mesh=new MGMeshC(DIMENSION,mesh_nlevels);


  int NXYZ0[3]; double axayaz[3]; double bxbybz[3];
//     if(_dim==2){
  NXYZ0[0]=NX; NXYZ0[1]=NY; NXYZ0[2]=NZ;       //  (NX,NY,NZ)
  axayaz[0]=AX; axayaz[1]=AY; axayaz[2]=AZ; //(0,bx)x(0,by)x(0,bz)
  bxbybz[0]=BX; bxbybz[1]=BY; bxbybz[2]=BZ; //(0,bx)x(0,by)x(0,bz)

  _mg_mesh->setNXNYNZ(NXYZ0); _mg_mesh->setblen(bxbybz); _mg_mesh->setalen(axayaz);
  _mg_mesh->gen_c(); // generation of Cartesian Mesh (NX,NY)
  _mg_mesh->print_hf5(mesh_nlevels);

  return;
}

// // // *******************************************************************
// // // **************** Solve  *******************************************
// // // *******************************************************************
/// This function sets up the intial set
void VOF::solve_setup(
  int         t_in,                 ///< initial time iteration
  double       time                 ///< actual time
) {
//   const int restart      = _mg_utils->get_par("restart");   // restart or not
  int print_step = _mg_utils->get_par("printstep");
  int    n_steps = _mg_utils->get_par("nsteps");
  const int ndigits     = _mg_utils->get_par("ndigits");
  int Nlev_cc=CCLEV;

  int iproc= _mg_mesh->_iproc;
  if (iproc==0){
  _mgcc->print_time_fine_xmf(t_in,n_steps+1,print_step,ndigits);     // cell cc fine
  _mgcc->print_time_cc_xmf(t_in,n_steps+1,print_step,ndigits);       // point cc coarse
#if DIMENSION==2
  // print initial time and time series
  _mgcc->print_time_interface_xmf(t_in,n_steps+1,print_step,ndigits);// interface coarse
#endif

  if(t_in%print_step==0) {
//     std::ostringstream namefile;
//     namefile << "./RESU/cc." << std::setw(ndigits) << std::setfill('0') << t_in <<  ".h5";
    std::ostringstream flagfile;
//     namefile << "./RESU/cc." << std::setw(ndigits) << std::setfill('0') << t_in <<  ".h5";
    flagfile << std::setw(ndigits) << std::setfill('0') << t_in;
    _mgcc->print_fine_hdf5(time, flagfile.str(),Nlev_cc-1); // cell cc at t=0 and level=fine
    _mgcc->print_cc_hdf5(time, flagfile.str(),0);               // point cc at t=0 and level=0
#if DIMENSION==2
    // print initial time and time series
    _mgcc->print_interface_hdf5(time, flagfile.str(),Nlev_cc-1);// cc interface at t=0 and level=fine
#endif
  }
  }
  return;
}

//=============================================================================
// This function solves one step  for transient problems
void VOF::solve_onestep(
  const int  & /*t_in*/,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt                   ///< step time
) { // ========================================================================
  int Nlev_cc=CCLEV;
  const int ndigits     = _mg_utils->get_par("ndigits");
  std::cout << "\n\n*** Solving time step " << t_step
            << ", time = " << time << " ***" << std::endl;
  // Solving advection equation ------------------
  std::cout <<"\n Vof step: ";
  double s=1.;//  s=cos(PI*(t_step+0.5)/(double)(N_TIME_STEPS+1));

  _mgcc->MGSolve(/**_mg_mesh,*/VOF_SUBSTEP,t_step,s*dt);
  // print ---------------------------------------
    int iproc= _mg_mesh->_iproc;
  if((t_step)%print_step == 0 && iproc==0 ) {
    std::ostringstream flagfile;
//     namefile << "./RESU/cc." << std::setw(ndigits) << std::setfill('0') << t_in <<  ".h5";
    flagfile << std::setw(ndigits) << std::setfill('0') << t_step;
    _mgcc->print_fine_hdf5(time, flagfile.str(),Nlev_cc-1);
    _mgcc->print_cc_hdf5(time, flagfile.str(),0);  // print cc node vector at level=0
#if DIMENSION==2
    _mgcc->print_interface_hdf5(time, flagfile.str(),Nlev_cc-1);// print matrix cc at level=Nlev_cc-2
#endif

  }
  return;
}







#ifdef HAVE_MED
// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
void VOF::init_interface(
  const int interface_name,
  int interface_id,
  const std::string & medfile_name // medfile name    (in)
) {// =========================================================================
  
  
  // name interface (interface_id)
  std::ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();

  // Reading mesh Names from med file  ------------------------------

//   std::string mesh_dir=_mg_utils->_mesh_dir;
//   std::string localFile=mesh_dir+medfile_name;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(medfile_name.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[0];//index_medmesh=0
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(medfile_name.c_str(), localMeshName.c_str());


//   int nGroups0 = GroupNames.size();
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id<10) {id_level=0;}
  support = MEDLoader::ReadUMeshFromGroups(medfile_name.c_str(), meshNames[0].c_str(), id_level,vG);
  support->zipCoords();

  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<
            interface_id  << "\n";

  InterfaceFunctionDD *fun = new InterfaceFunctionDD;
  // set the mesh interface
//   fun->set_mesh_interface_elemID(*_mg_mesh,support,interface_id);
  fun->set_mesh_interface_nodeID(/**_mg_mesh,*/support,interface_id,2);
  // setting the fun in the map _interfaceFunMap at bd_name_id
  _interfaceFunMap[interface_name] = fun;  // added to the map

  return;
}
// // // // routine2groups
void VOF::init_interface(
  const int interface_name,
  int interface_id_1,
  int interface_id_2,
  const std::string & medfile_name // medfile name    (in)
) {// =========================================================================
  
  
  // name interface (interface_id)
  std::ostringstream name_id_1; name_id_1 <<interface_id_1;
  std::ostringstream name_id_2; name_id_2 <<interface_id_2;
  std::vector<std::string> vG(2);
//   std::string id_name=name_id.str().c_str();
  vG[0] =name_id_1.str().c_str();
  vG[1] =name_id_2.str().c_str();


  // Reading mesh Names from med file  ------------------------------

//   std::string mesh_dir=_mg_utils->_mesh_dir;
//   std::string localFile=mesh_dir+medfile_name;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(medfile_name.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[0];//index_medmesh=0
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(medfile_name.c_str(), localMeshName.c_str());


//   int nGroups0 = GroupNames.size();
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id_1<10) {id_level=0;}
  support = MEDLoader::ReadUMeshFromGroups(medfile_name.c_str(), meshNames[0].c_str(), id_level,vG);
  support->zipCoords();
   
  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<
            interface_id_1  << "\n";

  InterfaceFunctionDD *fun = new InterfaceFunctionDD;
  // set the mesh interface
//   fun->set_mesh_interface_elemID(*_mg_mesh,support,interface_id);
  fun->set_mesh_interface_nodeID(/**_mg_mesh,*/support,interface_id_1,2);
    std::cout<<"\n-----------\n"<<support->getNumberOfNodes()<<"\n";
  // setting the fun in the map _interfaceFunMap at bd_name_id
  _interfaceFunMap[interface_name] = fun;  // added to the map

  return;
}
// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
void VOF::init_interface_oncell(
  const int interface_name,          ///< interface name
  int interface_id,                  ///< volume interface (in)
  const std::string & medfile_name   ///< medfile name    (in)
) {// =========================================================================

  // name interface (interface_id)
  std::ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();

  // Reading mesh Names from med file  --------------------------------------------------
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(medfile_name.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {std::cout<<  " VOF::setMesh : no meshes in the file'";}
  std::string localMeshName = meshNames[0];//index_medmesh=0
  std::vector<std::string> GroupNames =// Reading group names
    MEDLoader::GetMeshGroupsNames(medfile_name.c_str(), localMeshName.c_str());

  // MED Mesh -> support  ---------------------------------------------------------------
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id<10) {id_level=0;}
  support = MEDLoader::ReadUMeshFromGroups(medfile_name.c_str(), meshNames[0].c_str(), id_level,vG);
  support->zipCoords();

  // Volume Interface for color function C ----------------------------------------------
  std::cout << "VOF::setInterfaces: support  set to boundary with name "<< interface_id  << "\n";
  InterfaceFunctionDD *fun = new InterfaceFunctionDD;
  // set the mesh interface (for nodes)
  fun->set_mesh_interface_nodeID(/**_mg_mesh, */support,interface_id,2);
  // store the Volume Interface in the VOF interface map (_interfaceFunMap)
  _interfaceFunMap[interface_name] = fun;  // added to the map

  return;
}


// ======================================================================================
// This function gets the  the value of the color function
//  on nodes on the domain
ParaMEDMEM::MEDCouplingFieldDouble * VOF::getValues_C_cell(
)  {// ==================================================================================

  // mesh structure ---------------------------------------------------------------------
  int nxyz[3]= {NX,NY,NZ};     // number of divisions
  int n_cell=1; int n_pts=1;
  for(int idim=0; idim<DIMENSION; idim++) {n_cell *=nxyz[idim]; n_pts *=(nxyz[idim]+1);}

  // color field from levelo 0 -> cc vector
  double *cc_tmp=new double[n_pts];
  for(int ipt=0; ipt<n_pts; ipt++)    cc_tmp[ipt]=  _mgcc->cc.Cmp[ipt+1];

  // MED field f ------------------------------------------------------------------------
  // array function to fill f
  ParaMEDMEM::DataArrayDouble *array = ParaMEDMEM::DataArrayDouble::New();
  array->alloc(n_pts,1); // array(n_of_tuples,n_of_cmp);
  for(int i=0; i < n_pts; i++) array->setIJ(i,0,cc_tmp[i]);  //  filling array C
  // field f ----------------------------------------------------------------------------
  ParaMEDMEM::MEDCouplingFieldDouble * f =
    ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(NULL);               // no set mesh -> no check please !
  f->setName("C");                // set name
  f->setArray(array);             // set array -> f

  // check and  clean -------------------------------------------------------------------
  array->decrRef();     // delete array
  delete [] cc_tmp;

  // print ------------------------------------------------------------------------------
  std::cout <<  "\n GetValuesC_nodes: " <<n_pts<< " values " << std::endl;
  return f;

}

// ============================================================================
// / This function computes the initial field velocity by using the MGSolverCC
// / class (trough InitVel)
void VOF::setFieldSource_Vinit(
  const int Level, ///< Level
  const double dt ///< time step
) { // ========================================================================
  _mgcc->GenVel(Level,dt);
  return;
}



void VOF::setFieldSource(
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * srcField) {
//
  InterfaceFunctionDD * fct = get_interface_fun(interface_name);
  double dt= _mg_utils->get_par("dt");
  _mgcc->setFieldSource(dt,interface_name,n_cmp,srcField,fct);


  return;
}
void VOF::setFieldSource_disp(
  int interface_name,
  int n_cmp,
  const std::vector<ParaMEDMEM::MEDCouplingFieldDouble*>& srcField

//   const ParaMEDMEM::MEDCouplingFieldDouble * srcField
) {
//
  InterfaceFunctionDD * fct = get_interface_fun(interface_name);
  double dt= _mg_utils->get_par("dt");
  _mgcc->setFieldSource_disp(dt,interface_name,n_cmp,srcField,fct);


  return;
}


  const ParaMEDMEM::MEDCouplingUMesh* VOF::getUMesh(
  int name
){
  InterfaceFunctionDD * fct = get_interface_fun(name);
  const  ParaMEDMEM::MEDCouplingUMesh *sourceMesh = fct->getSupport();
  
  
  return sourceMesh;
}
#endif
