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
#include "InterfaceProjection.h"


#include "VOF.h"
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

int main(int argc, char** argv)
{

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
    filenameP[i_mesh] <<   mesh_nameP[i_mesh].substr(0,posP)  << "_MedToMg.med" ;
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
  std::vector<FIELDS> myproblemP; myproblemP.resize(1);
  // Problem to solve for each mesh
  myproblemP[0]=FS_F;
//   myproblemP[1]=T_F;
  int group_mesh=2;
  int group_mesh2=4;
  int interface_vel=1; int interface_c=2; int interface_c2=3;

  // Navier-Stokes ----------------------------------------------------------------------
  FEMUS P1;                                            // constructor
  P1.init_param(*mgutils[0]);                          // init parameter
  P1.init_fem(*mggeomel,*mgfemap);                     // init fem
  P1.setMesh();                                        // setting mesh
  P1.setSystem(myproblemP);                            // set system: Navier-Stokes

  // NS interfaces
  P1.init_interface(interface_vel,group_mesh,group_mesh2,2,filenameP[0].str().c_str());
  P1.init_interface(111,group_mesh,group_mesh2,2,filenameP[0].str().c_str());
  P1.init_interface(interface_c2,group_mesh,group_mesh2,2,filenameP[0].str().c_str());


  // VOF --------------------------------------------------------------------------------
  VOF P2;                                        // constructor
  P2.set_param(*mgutils[0]);                     // init parameter
  P2.setMesh();                                  // set vof_mesh
  P2.setSystem(1);
  P1.set_mgcc(P2.get_MGSystem()); // new for surface tension
  // VOF interfaces
  #if DIMENSION==3
  P2.init_interface(interface_vel,group_mesh,group_mesh2,  osfilenameP[0].str().c_str());//vofsi_20x20x20_test.med
#else
  P2.init_interface(interface_vel,group_mesh,group_mesh2, osfilenameP[0].str().c_str());//vofsi_20x20x20_test.med
  P2.init_interface(112,group_mesh,group_mesh2, osfilenameP[0].str().c_str());//vofsi_20x20x20_test.med
#endif  
  const ParaMEDMEM::MEDCouplingUMesh *SourceMesh = P1.getUMesh(interface_vel);
  const ParaMEDMEM::MEDCouplingUMesh *TargetMesh = P2.getUMesh(112);  
  
  int    n_steps = mgutils[0]->get_par("nsteps");
  double      dt = mgutils[0]->get_par("dt");
  int print_step = mgutils[0]->get_par("printstep");
  int init_time=0;  double time=0.;

//     P2.setFieldSource_Vinit(0,dt);  // inital vel field from cc

  // Initial time ---------------------------------------------------------------------
  // Computation  of the initial volume
//     Matrix cc0; // cell cc
//     int Nlev_cc=CCLEV;
//     int CLevel=mgcc->_Nlev_cc-1;
//     int rlen=M_GetRowDim(&(mgcc->c1[CLevel]));
//     int clen=M_GetClmDim(&(mgcc->c1[CLevel]));
//     M_Constr(&cc0,(char *)"cc0",rlen,clen,Rowws,Normal,_LPTrue);// cc at t=0 and level=0
//     mgcc->OldSol_update(CLevel,mgcc->c1[CLevel],cc0);

//   P.solve_setup(init_time,time);
  P2.solve_setup(init_time,time);
  P1.solve_setup(init_time,time);
//   ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;
  ParaMEDMEM::MEDCouplingFieldDouble *bdy2=NULL;
  std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> disp; disp.resize(DIMENSION);
  {
    for(int i=0; i<DIMENSION; i++) disp[i]=NULL;
  }
  ParaMEDMEM::MEDCouplingFieldDouble *bdy_c2=NULL;
  bdy_c2 = P2.getValues_C_cell();///< int boundary identity   (in)
  // -------------------------------------------------------
  //  time loop
  // -------------------------------------------------------
  
  ParaMEDMEM::MEDCouplingFieldDouble *TargetField;
  
  for(int t_step=init_time+1; t_step<  n_steps+1; t_step++) {
    time += dt;
    std::cout<<"Time before Femus "<< time<<std::endl;
//     P1.setFieldSource(interface_c,1,bdy_c);
    P1.setFieldSource(interface_c2,1,bdy_c2);
    P1.solve_onestep(init_time,t_step,print_step,time,dt);    // solving P
    
    
    
    
     /// con aggiornamento spostamento mesh
    disp[0]= P1.getValuesOnBoundary(111, "SDSX",1);
    disp[1]= P1.getValuesOnBoundary(111, "SDSY",1);
    P1.update_interface(111,DIMENSION, disp);// update velocity interface
    const  ParaMEDMEM::MEDCouplingUMesh * SourceMesh_update  = P1.getUMesh(111);
    bdy2 = P1.getValuesOnBoundary(111, "FSI0",DIMENSION);    
    BoundInterp VOL_updated = BoundInterp(SourceMesh_update,TargetMesh,Volume);
    // BoundInterp VOL_updated2 = BoundInterp(SourceMesh_update,TargetMesh,Volume);
    TargetField = VOL_updated.InterpolatedField(bdy2);
    if(t_step == n_steps){
      VOL_updated.PrintMed(bdy2,"sourcefield2",0);
      VOL_updated.PrintMed(TargetField,"targetfield2",0);
    }
    P2.setFieldSource(112,DIMENSION,TargetField);
    
    
    /// senza aggiornamento spostamento mesh
//     bdy2 = P1.getValuesOnBoundary(111, "FSI0",DIMENSION);    
//     P2.setFieldSource(112,DIMENSION,bdy2);
    
//     P2.setFieldSource_disp(interface_vel,1,disp);
    P2.solve_onestep(init_time,t_step,print_step,time,dt);                   ///< step time
    bdy_c2 = P2.getValues_C_cell();///< int boundary identity   (in)
  }   // end time loop
// // --------------------------------------------------
    
// // Volume computation
// // --------------------------------------------------
// // // #ifdef TWO_PHASE
// //     double c=0.; double b=0.; double a=0.;
// //
// //     double *tempc1; double *tempc2;
// //     tempc1=new double[clen+1]; tempc2=new double[clen+1];
// //
// //     for(unsigned int jy=0; jy<rlen; jy++) {
// //       for(unsigned int ix=0; ix<=clen; ix++) {
// //         tempc1[ix]=0.; tempc2[ix]=0.;
// //       }
// //       for(unsigned int ix=0; ix<M__GetLen(&(mgcc->c1[CLevel]),jy+1); ix++)  {
// //         int pos=M__GetPos(&(mgcc->c1[CLevel]),jy+1,ix)-1;
// //         double val=M__GetVal(&(mgcc->c1[CLevel]),jy+1,ix);
// //         if(val <1.)  tempc1[pos]=val;
// //         else for(int ki=0; ki< (int) val; ki++) tempc1[pos+ki]=1.;
// //       }
// //       for(unsigned int ix=0; ix<M__GetLen(&cc0,jy+1); ix++)  {
// //         int pos=M__GetPos(&cc0,jy+1,ix)-1;
// //         double val=M__GetVal(&cc0,jy+1,ix);
// //         if(val <1.)  tempc2[pos]=val;
// //         else for(int ki=0; ki< (int) val; ki++) tempc2[pos+ki]=1.;
// //       } ~VOF();
// //   void terminate();
// //       for(unsigned int ix=0; ix<clen; ix++) {
// //         a +=fabs(tempc1[ix]-tempc2[ix]); b +=tempc1[ix];	c +=tempc2[ix];
// //       }
// //     }
// //     delete []tempc1; delete []tempc2;
// //
// //     double err_m =fabs(c-b)/fabs(c);
// //     double err_r=a/c;
// //     double err_g =a/((NX+1.)*(NX+1.)*(NX+1.));
// //     printf("\n errm %e \n errr %e \n errg %e \n ",err_m,err_r,err_g);
// #endif

// // //   clean --------------------------------
  mgutils.clear();
  myproblemP.clear();
//   P1.terminate();
//   bdy->decrRef();
//   bdy_c->decrRef();


  delete mggeomel;
  delete dfe_q;
  delete dfe_l;
  delete dfe_k;
  delete  mgfemap;
  return 0;
}

#ifdef HAVE_MED
// ==========================================================================
ParaMEDMEM::MEDCouplingFieldDouble *build_field(
  const std::string & path_mesh,
  const std::string & expression,
  int id_interface
)
{
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

