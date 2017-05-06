// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include<map>
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
#include <EquationSystemsExtendedM.h>
// #include <EquationSystemsExtendedM.h>
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
// #include "MMed.h"
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
  std::string MeshFileName[NUM_MESH];
  std::string MeshFileName_Gen[NUM_MESH];  

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    // MGUtils constructor ----------------------------------------------------
    mgutils.push_back(new MGUtils(i_mesh+1));
    // mesh names -------------------------------------------------------------
    MeshFileName[i_mesh]  = mgutils[i_mesh]->get_file("F_MESH_READ"); 
    int posP            = MeshFileName[i_mesh].find(".");  
    MeshFileName_Gen[i_mesh]    = MeshFileName[i_mesh].substr(0,posP) + "_MedToMg.med" ;    
    
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "<< mgutils[i_mesh]->_mesh_dir + MeshFileName_Gen[i_mesh]<<"\n "; 
  }
  std::cout<<" ============ end loop mesh ================================ \n";
  std::cout<<" =========================================================== \n";
  
  // FEM class -----------------------------------------------------------
  MGFEMap *mgfemap; mgfemap=new MGFEMap();
  
  // MGFEMap mg_femap;
  MGFE *dfe_q;    
  dfe_q=new MGFE(2,ELTYPE); 
  dfe_q->init_qua();
  mgfemap->set_FE(dfe_q); // initialize quadratic fem
  
  MGFE *dfe_l;  
  dfe_l=new MGFE(1,ELTYPE); 
  dfe_l->init_lin();
  mgfemap->set_FE(dfe_l); // initialize linear fem
  
  MGFE *dfe_k; 
  dfe_k=new MGFE(0,ELTYPE);  
  dfe_k->init_pie();
  mgfemap->set_FE(dfe_k); // initialize piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();
  
  // system problem =========================================================
  std::vector<FIELDS> myproblemP1; 

#ifdef NS_EQUATIONS
  myproblemP1.push_back(NS_F);
#endif
#ifdef TBK_EQUATIONS
  myproblemP1.push_back(K_F);
#endif
#ifdef TTBK_EQUATIONS
  myproblemP1.push_back(KTT_F);
#endif   
#ifdef T_EQUATIONS
  myproblemP1.push_back(T_F);
#endif
#ifdef FSI_EQUATIONS
  myproblemP1.push_back(FS_F);
#endif
#ifdef FSIA_EQUATIONS
  myproblemP1.push_back(FSA_F);
#endif
    
  
  
  
  // Object for the creation of FEMus problems interfaces
  MMed * interfaces= new BoundInterp();
  interfaces->InitFe();
  
  std::vector<FEMUS*> FEMusProblems;
  FEMusProblems.resize(NUM_MESH);
  
//   for(int i=0; i<NUM_MESH; i++){
// //     string name = std::to_string(i);
// //     FEMUS name;
//     FEMusProblems[i] = new FEMUS;
//     FEMusProblems[i] -> init_param(*mgutils[i],i);
//     FEMusProblems[i] -> init_fem(*mggeomel,*mgfemap);                 // init fem      
//     FEMusProblems[i] -> setMesh();                                    // set MGmesh   
//     FEMusProblems[i] -> setSystem(myproblemP1);                       // set system
//     FEMusProblems[i] -> setMeshName(MeshFileName[i]);
//   }
  
  FEMUS FEMus;
//     FEMus.init_param(*mgutils[0],0);
//     FEMus.init_fem(*mggeomel,*mgfemap);                 // init fem      
//     FEMus.setMesh();                                    // set MGmesh   
//     FEMus.setSystem(myproblemP1);                       // set system
//     FEMus.setMeshName(MeshFileName[0]);
  
  const int QuadInterfaceId = 22;
  const int LinInterfaceId = 11;
  
  for(int nDom=0; nDom<NUM_MESH; nDom++){
      
//     FEMUS FEMus;
    FEMus.init_param(*mgutils[nDom],nDom);
    FEMus.init_fem(*mggeomel,*mgfemap);                 // init fem      
    FEMus.setMesh();                                    // set MGmesh   
    FEMus.setSystem(myproblemP1);                       // set system
    FEMus.setMeshName(MeshFileName[nDom]); 
      
      
    std::vector<std::string> meshNames = MEDLoader::GetMeshNames(mgutils[nDom]->_mesh_dir +  MeshFileName[nDom]);
    std::vector<std::string> GroupNames = MEDLoader::GetMeshGroupsNames(mgutils[nDom]->_mesh_dir +  MeshFileName[nDom], meshNames[0]);
    int nGroups0 = GroupNames.size();
    std::cout<<"\033[1;31m Total number of system of equations: "<<myproblemP1.size()
	    <<" number of FEMus problems "<<nDom<<"\033[0m \n\n";

      interfaces->CreateInterfaces(&FEMus);
      
      std::string SourceFile = "/RESU_MED/SourceField" + std::to_string(nDom) + ".med";
      std::string QuadFieldsFile = mgutils[nDom]->_app_dir +  SourceFile;
      const  ParaMEDMEM::MEDCouplingUMesh * SourceV = MEDLoader::ReadUMeshFromFile(QuadFieldsFile, 0);  
      const  ParaMEDMEM::MEDCouplingUMesh * TargetV = FEMus.getUMesh(QuadInterfaceId);  
      MEDLoader::WriteUMesh("RESU_MED/Source.med", SourceV, 1 );
      MEDLoader::WriteUMesh("RESU_MED/Target.med", TargetV, 1 );

      interfaces->FillParameters(SourceV,TargetV,Volume);

	std::cout << "\033[1;31m ============================================================= \033[0m \n" ;
	std::cout << "\033[1;31m ********************** INTERPOLATION STEP ******************* \033[0m \n" ;  
	std::cout << "\033[1;31m ============================================================= \033[0m \n" ;  

// 	//////////////////////////////////////////////////////////////////////////////////
// 	/////				READ SOLUTION IN MED FILE
// 	//////////////////////////////////////////////////////////////////////////////////

	std::cout<<" Reading solution from fields in " <<QuadFieldsFile<<std::endl;
// 	
	std::vector<std::string> FieldNames = MEDLoader::GetAllFieldNames(QuadFieldsFile);
	std::vector<std::string> MeshNames = MEDLoader::GetMeshNames(QuadFieldsFile);
	double time = MEDLoader::GetTimeAttachedOnFieldIteration(QuadFieldsFile,FieldNames[0],-1,-1);
	std::cout<<"\n Time read from file SourceField0.med t = "<<time<<"\n";

	ParaMEDMEM::MEDCouplingFieldDouble* SourceField=NULL;
	ParaMEDMEM::MEDCouplingFieldDouble * TargetField=NULL;
	
	for(int l=0; l<FieldNames.size(); l++) {
	  std::cout<<"\033[1;36m\n Interpolation of field "<<FieldNames[l]<<"\033[0m\n";
//     
	  SourceField =  MEDLoader::ReadField(ParaMEDMEM::ON_NODES,QuadFieldsFile,MeshNames[0],0,FieldNames[l],-1,-1);      
	  int order = 2;
	  if (FieldNames[l]=="NSP"||FieldNames[l]=="NS2P") order = 1;
	                                               
	  TargetField = interfaces->InterpolatedField(SourceField,order);

          int nComp = TargetField->getNumberOfComponents();
          
          int InterfaceId = (order==1)? LinInterfaceId:QuadInterfaceId;
          
	  FEMus.setFieldBoundaryValues(InterfaceId,nComp,TargetField);
	  FEMus.write_Boundary_value(InterfaceId,FieldNames[l],nComp);
          
	  for(int i=0; i<nComp; i++){
              
            if((FieldNames[l]=="NS0" || FieldNames[l]=="FSI0" || FieldNames[l]=="FSIA0"||FieldNames[l]=="NSA0")&& i == nComp-1)  order=1; 
	    double in_norm = interfaces->Integrate(SourceField,order, 1, i, NormL2);
	    double out_norm = interfaces->Integrate(TargetField,order, 1, i, NormL2);

	    std::cout<<std::setprecision(7)<<" Initial L2 Norm of component "<< i <<":\t" << in_norm<<std::endl;
	    std::cout<<std::setprecision(7)<<" Final L2 Norm of component "<< i <<":\t\t" << out_norm<<std::endl;
	    std::cout<<std::setprecision(2)<<" Percentage of integral value variation from initial value\t \033[1;36m"
			  		   <<100.*(out_norm-in_norm)/(in_norm+1.e-20)<<" % \033[0m"<<std::endl;
	    std::cout<<"---------------------------------------------------------------------------------\n";	    
	  }
// 	  interfaces->PrintMed(SourceField,"Source_"+ FieldNames[l],1);
// 	  interfaces->PrintMed(TargetField,"Target_"+ FieldNames[l],1);
	  order=2;
	}
//     
	  SourceField->decrRef();
	  TargetField->decrRef();
          double dt=0.1;
// 	/*============================ Solution of P1 ==================================*/	 
	FEMus.dummy_step(0,1000,1,time,dt);    // solving P
	std::cout << "\033[1;31m ============================================================= \033[0m \n" ;
	std::cout << "\033[1;31m ********************** SOLUTION PRINTED ********************* \033[0m \n" ;  
	std::cout << "\033[1;31m ============================================================= \033[0m \n" ;  
// 	/*==============================================================================*/    
// 
// 	//////////////////////////////////////////////////////////////////////////////////
// 	/////				END TIME LOOP
// 	//////////////////////////////////////////////////////////////////////////////////
// 	
// 	GroupNames.clear();
	MeshNames.clear();
	FieldNames.clear();
    
//         FEMus.terminate();
//          FEMus.~FEMUS();
  }
  FEMus.terminate();
  // Wait all processes
//   MPI_Barrier( MPI_COMM_WORLD);
// 
//   
//   // Terminating FEMUS problems from last fo first
// //   for(int nDom=0; nDom<NUM_MESH; nDom++) FEMusProblems[NUM_MESH-nDom]->terminate();
// //   FEMusProblems.clear();
  std::string AppDir=mgutils[0]->_app_dir;
// 
//   
//   // clean --------------------------------------------------------------------
  mgutils.clear();
  
  delete dfe_q;  
  delete dfe_l;   
  delete dfe_k;  
  delete mggeomel; 
  delete mgfemap;
  interfaces->~MMed();
  std::cout<<"\033[1;31m \n\n\n \
  In order to use the Interpolation application you need to define \n \
   \t\"HAVE_MED\" \n \
   inside the file $"<< AppDir<<"/DATA/Solverlib_conf.h \n\
   -) Make a resu_clean \n\
   -) Set restart to 0    \n\
   -) Set the new meshes in \"param_files\"     \n\
   -) Run gencase with n processors          \n\
   -) Run Interpolator with n processors     \n\
   -) Set restart to 1000                   \n\
   -) Restart your application \033[0m \n \n\n\n";
  return 0;
}
