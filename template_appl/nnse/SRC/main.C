// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// LibMesh library included ------------------------------
// #ifdef LM_INIT
// #include "libmesh.h" // for Libmesh library
// #endif

// solver library -------------------------------------
#include  "Solverlib_conf.h"  // Solver library options 
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
#include "MGMesh.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"
#include "FEMUS.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"

// ==========================================================================
  ParaMEDMEM::MEDCouplingFieldDouble *build_field(
   const  std::string & path_mesh,
   const  std::string & expression,
    int id_interface
  );
#endif
/// Set up
// =======================================
// Main program
// =======================================

int main(int argc, char** argv) {

  argc = argc ; argv=argv;  // no unused warning
 
  std::cout<<" ============ MGUtils ===================================== \n";
  std::cout<<" =========================================================== \n";
  // setting MGUtils -> parameters and file name ------------------------
  std::vector<MGUtils*> mgutils;
  std::string mesh_nameP[NUM_MESH];
  std::ostringstream filenameP[2];  std::ostringstream osfilenameP[2];

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    // MGUtils constructor ----------------------------------------------------
    mgutils.push_back(new MGUtils(i_mesh+1));
    // mesh names -------------------------------------------------------------
    mesh_nameP[i_mesh]= mgutils[i_mesh]->get_file("F_MESH_READ"); // name mesh
    int posP = mesh_nameP[i_mesh].find(".");  // position of "live" in str
    filenameP[i_mesh] <<   mesh_nameP[i_mesh].substr(0,posP)  << "_gen.med" ;
    osfilenameP[i_mesh]<< mgutils[i_mesh]->_mesh_dir <<filenameP[i_mesh].str();
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "
                                  << osfilenameP[i_mesh].str().c_str() <<"\n "; 
  }
  std::cout<<" ============ end loop mesh ================================ \n";
  std::cout<<" =========================================================== \n";
  
// FEM class -----------------------------------------------------------
  MGFEMap *mgfemap; mgfemap=new MGFEMap();
  // MGFEMap mg_femap;
  MGFE *dfe_q;    dfe_q=new MGFE(2,ELTYPE); dfe_q->init_qua();
  mgfemap->set_FE(dfe_q); //// initialize quadratic fem
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE); dfe_l->init_lin();
  mgfemap->set_FE(dfe_l); //initialize linear fem
  MGFE *dfe_k; dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();
  mgfemap->set_FE(dfe_k); //  initialize piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();
  
  // system problem =========================================================
  std::vector<FIELDS> myproblemP; myproblemP.resize(2);
  // Problem to solve for each mesh
   myproblemP[0]=NS_F; 
//    myproblemP[1]=FS_F; 
  myproblemP[1]=T_F;
//     myproblemP[1]=K_F;
//   myproblemP[3]=KTT_F;
  
  
//   std::vector<FIELDS> myproblemPf; myproblemP.resize(1);
//   myproblemP[0]=WF; 
  // system 1
  // MGFemusInit --------------------------------------------------------------
  FEMUS P;                                        // constructor
  P.init_param(*mgutils[0]);                      // init parameter
  P.init_fem(*mggeomel,*mgfemap);                 // init fem      
  // setting mesh -------------------------------------------------------------
//   P.setMedMesh(osfilenameP[0].str().c_str());     // set med-mesh
  P.setMesh();                                    // set MGmesh   
  // setting system -----------------------------------------------------------
  P.setSystem(myproblemP);                         // set system

  
  
// #ifdef WALL_FUNC_APP
//   FEMUS Pf;                                        // constructor
//   Pf.init_param(*mgutils[0],1);                      // init parameter
//   Pf.init_fem(*mggeomel,*mgfemap);                 // init fem      
//   // setting mesh -------------------------------------------------------------
// //   P.setMedMesh(osfilenameP[0].str().c_str());     // set med-mesh
//   Pf.setMesh();                                    // set MGmesh   
//   // setting system -----------------------------------------------------------
//   Pf.setSystem(myproblemP);                         // set system
//   
//   P.init_interface(1,1,2,filenameP[0].str().c_str());
//   P.init_interface(2,1,2,filenameP[0].str().c_str());
//   P.init_interface(3,1,2,filenameP[0].str().c_str());
//   P.init_interface(4,1,2,filenameP[0].str().c_str());
//   P.init_interface(5,1,2,filenameP[0].str().c_str());
//   P.init_interface(6,1,2,filenameP[0].str().c_str());
//   P.init_interface(7,1,2,filenameP[0].str().c_str());
//   P.init_interface(8,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(1,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(2,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(3,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(4,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(5,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(6,1,2,filenameP[0].str().c_str());
//   Pf.init_interface(7,1,2,filenameP[0].str().c_str()); 
// #if DIMENSION == 3
//   Pf.init_interface(8,1,2,filenameP[0].str().c_str());
// #endif
//   ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;
// #endif  
  // solving
  int    n_steps = mgutils[0]->get_par("nsteps");
  double      dt = mgutils[0]->get_par("dt");
  int print_step = mgutils[0]->get_par("printstep");
  int    itime_0  = mgutils[0]->get_par("itime");
  double time    = 0.;
  
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)
  
  // transient loop with  n_steps (i0time ->i0time+ n_steps)
  for(int itime=itime_0; itime<= itime_0 + n_steps; itime++) {
      P.solve_onestep(itime_0,itime,print_step,time,dt);    // solving P
      time +=dt;
  }   // end time loop
// #ifdef WALL_FUNC_APP
//     bdy = P.getValuesOnBoundary(1, "NS0X",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(1,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(1,"NS0X",1);               // write inside x_old the boundary values
//     
//     bdy = P.getValuesOnBoundary(2, "NS0Y",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(2,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(2,"NS0Y",1);               // write inside x_old the boundary values
//  
// #if DIMENSION==3    
//     bdy = P.getValuesOnBoundary(8, "NS0Z",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(8,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(8,"NS0Z",1);               // write inside x_old the boundary values
// #endif
//     
//     bdy = P.getValuesOnBoundary(3, "T",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(3,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(3,"T",1);               // write inside x_old the boundary values
//     
//     bdy = P.getValuesOnBoundary(4, "K",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(4,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(4,"K",1);               // write inside x_old the boundary values
//     
//     bdy = P.getValuesOnBoundary(5, "K2",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(5,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(5,"K2",1);               // write inside x_old the boundary values
//     
//     bdy = P.getValuesOnBoundary(6, "TK",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(6,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(6,"TK",1);               // write inside x_old the boundary values
//     
//     bdy = P.getValuesOnBoundary(7, "TK2",1);          // take field from P problem
//     Pf.setFieldBoundaryValues(7,1, bdy);                // store field inside P1 problem
//     Pf.write_Boundary_value(7,"TK2",1);               // write inside x_old the boundary values    
// //     
// //     
// Pf.solve_onestep(itime_0 + n_steps,itime_0 + n_steps + 1,1,time,dt);    // solving P
// #endif

  // end ======================================================================
  // --------------------------------------------------------------------------
  P.terminate(); 
  
// #ifdef WALL_FUNC_APP 
//   Pf.terminate();
// #endif
  // clean --------------------------------------------------------------------
  mgutils.clear();
  delete dfe_q;  
  delete dfe_l;
  delete dfe_k;  // delete   fem
  delete mggeomel; delete mgfemap;
  return 0;
}

#ifdef HAVE_MED
// ==========================================================================
  ParaMEDMEM::MEDCouplingFieldDouble *build_field(
    const std::string & path_mesh,
    const std::string & expression,
    int id_interface
  ){
  std::vector<std::string> vG(1);
  std::ostringstream s;
  s<< id_interface;
  vG[0] =s.str();
  ParaMEDMEM::MEDCouplingUMesh * support = 
  MEDLoader::ReadUMeshFromGroups(path_mesh.c_str(),"Mesh_1", 0,vG);
   
  
    ParaMEDMEM::TypeOfField type = ParaMEDMEM::ON_NODES;
  int dim=support->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if(dim > 0)    vars[0] = "x";
  if(dim > 1)    vars[1] = "y";
  if(dim > 2)    vars[2] = "z";

  const ParaMEDMEM::MEDCouplingFieldDouble * field_nodes =
    support->fillFromAnalytic3(type, 1, vars, expression.c_str());

  int n_elements_med= support->getNumberOfCells(); // n MED-elements
//   const ParaMEDMEM::DataArrayDouble * d =support->getCoords();// Med-mesh coordinates
  int n_nodes_med = field_nodes->getNumberOfTuples();             //  MED-nodes
//   int dim = field_nodes->getNumberOfComponents();             //  MED dimension

  int n_nod_el=NDOF_FEM;
  double sum;

  
 ParaMEDMEM::MEDCouplingFieldDouble *field;
 field= ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME);
 field->setMesh(support);

  ParaMEDMEM::DataArrayDouble *array;
  array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(n_elements_med,1);

  double *field_av = new double [n_elements_med]; // MED el centers
//   std::vector<std::pair<double, int > > dist_med;       // MED el distance
//   std::vector<std::pair<double, int > >::iterator itdist;

//   _elemID.clear();  // clear map elements: FEMus-mesh -> MED-mesh
  std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
    sum=0.; // zeros
    support->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    for(int inode=0; inode< n_nod_el; inode++) {
//       for(int idim=0; idim<dim; idim++)
      sum += field_nodes->getIJ(nodes1[inode],0);
    } // end inode
    nodes1.clear(); // clear element node vector
//     for(int idim=0; idim<dim; idim++)

    field_av[ielem]=sum/n_nod_el;
  }// *************************************************************************
  
   std::copy(field_av,field_av+n_elements_med,array->getPointer());

  field->setArray(array);
  field->setName("pippo");
  field->checkCoherency();
  
  return field;
  }
  
#endif

