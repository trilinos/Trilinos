// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_config.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalData.hpp"

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"

#include "Epetra_MpiComm.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"

#include "Example_BCStrategy_Factory.hpp"
#include "Example_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Example_EquationSetFactory.hpp"

#include "AztecOO.h"

using Teuchos::RCP;
using Teuchos::rcp;


//
// This is the Mathematica code used to generate this example.
// It also generates a plot of the vector field so it is clear
// what the solution is doing.
// 
//      Needs["VectorAnalysis`"]
//
//      phi0[x_,y_]=(1-x)*(1-y)
//      phi1[x_,y_]=x*(1-y)
//      phi2[x_,y_]=x*y
//      phi3[x_,y_]=y*(1-x)
//      
//      psi0[x_,y_]={1-y,0,0}
//      psi1[x_,y_]={0,x,0}
//      psi2[x_,y_]={y,0,0}
//      psi3[x_,y_]={0,1-x,0}
//      
//      u[x_,y_]=phi2[x,y]*psi0[x,y]+phi3[x,y]*psi1[x,y]+phi0[x,y]*psi2[x,y]+phi1[x,y]*psi3[x,y]
//      f[x_,y_]=u[x,y]+Curl[Curl[u[x,y],Cartesian[x,y,z]],Cartesian[x,y,z]]
//      
//      TwoDVec[g_]={g[[1]],g[[2]]}
//      
//      DotProduct[u[0.5,0],{1,0,0}]
//      DotProduct[u[1,0.5],{0,1,0}]
//      DotProduct[u[0.5,1],{1,0,0}]
//      DotProduct[u[0,0.5],{0,1,0}]
//      
//      Out[118]= 0.
//      Out[119]= 0.
//      Out[120]= 0.
//      Out[121]= 0.
//      
//      VectorPlot[TwoDVec[u[x,y]],{x,0,1},{y,0,1}]
//      Simplify[u[x,y]]
//      Simplify[f[x,y]]
//      
//      Out[144]= {-(-1+y) y,-(-1+x) x,0}
//      Out[145]= {2+y-y^2,2+x-x^2,0}
//


void testInitialization(panzer::InputPhysicsBlock& ipb,
		       std::vector<panzer::BC>& bcs);

// calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using panzer::StrPureBasisPair;
   using panzer::StrPureBasisComp;

   Teuchos::GlobalMPISession mpiSession(&argc,&argv);
   RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setOutputToRootOnly(0);
   out.setShowProcRank(true);

   // variable declarations
   ////////////////////////////////////////////////////

   // factory definitions
   Example::EquationSetFactory eqset_factory; // where poison equation is defined
   Example::BCStrategyFactory bc_factory;    // where boundary conditions are defined 

   panzer_stk::SquareQuadMeshFactory mesh_factory;

   // other declarations
   const std::size_t workset_size = 20;

   // construction of uncommitted (no elements) mesh 
   ////////////////////////////////////////////////////////

   // set mesh factory parameters
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("X Elements",10);
   pl->set("Y Elements",10);
   mesh_factory.setParameterList(pl);

   RCP<panzer_stk::STK_Interface> mesh = mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);

   // construct input physics and physics block
   ////////////////////////////////////////////////////////

   panzer::InputPhysicsBlock ipb;
   std::vector<panzer::BC> bcs;
   std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
   {
      bool build_transient_support = false;

      testInitialization(ipb, bcs);
      
      int base_cell_dimension = mesh->getCellTopology("eblock-0_0")->getDimension();
      const panzer::CellData volume_cell_data(workset_size, base_cell_dimension,mesh->getCellTopology("eblock-0_0"));

      // GobalData sets ostream and parameter interface to physics
      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      // the physics block nows how to build and register evaluator with the field manager
      RCP<panzer::PhysicsBlock> pb 
	= rcp(new panzer::PhysicsBlock(ipb, "eblock-0_0", volume_cell_data, eqset_factory, gd, build_transient_support));

      // we can have more than one physics block, one per element block
      physicsBlocks.push_back(pb);
   }

   // finish building mesh, set required field variables and mesh bulk data
   ////////////////////////////////////////////////////////////////////////

   {
      RCP<panzer::PhysicsBlock> pb = physicsBlocks[0]; // we are assuming only one physics block

      const std::vector<StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrPureBasisPair,StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // add basis to DOF manager: block specific
      std::set<StrPureBasisPair,StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
         mesh->addSolutionField(fieldItr->first,pb->elementBlockID());

      mesh_factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   }

   // build worksets
   ////////////////////////////////////////////////////////

   Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
      = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
   Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
      = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

   // build DOF Manager and linear object factory
   /////////////////////////////////////////////////////////////
 
   // build the connection manager 
   const Teuchos::RCP<panzer::ConnManager<int,int> > 
     conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

   panzer::DOFManagerFactory<int,int> globalIndexerFactory;
   RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(MPI_COMM_WORLD,physicsBlocks,conn_manager);

   // construct some linear algebra object, build object to pass to evaluators
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
         = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm.getConst(),dofManager));

   // setup closure model
   /////////////////////////////////////////////////////////////
 
   // Add in the application specific closure model factory
   panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
   Example::ClosureModelFactory_TemplateBuilder cm_builder;
   cm_factory.buildObjects(cm_builder);

   Teuchos::ParameterList closure_models("Closure Models");
   // closure_models.sublist("solid").sublist("SOURCE_EFIELD_X").set<double>("Value",1.0); // a constant source
   // closure_models.sublist("solid").sublist("SOURCE_EFIELD_Y").set<double>("Value",2.0); // a constant source
   // closure_models.sublist("solid").sublist("SOURCE_EFIELD").set<std::string>("Scalar Names","SOURCE_EFIELD_X, SOURCE_EFIELD_Y"); // a constant source
   closure_models.sublist("solid").sublist("SOURCE_EFIELD").set<std::string>("Type","SIMPLE SOURCE"); // a constant source
      // SOURCE_EFIELD field is required by the CurlLaplacianEquationSet

   Teuchos::ParameterList user_data("User Data"); // user data can be empty here

   // setup field manager builder
   /////////////////////////////////////////////////////////////

   Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
         Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);
   fmb->setupVolumeFieldManagers(*wkstContainer,physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
   fmb->setupBCFieldManagers(*wkstContainer,bcs,physicsBlocks,eqset_factory,cm_factory,bc_factory,closure_models,
                             *linObjFactory,user_data);

   // setup assembly engine
   /////////////////////////////////////////////////////////////

   // build assembly engine: The key piece that brings together everything and 
   //                        drives and controls the assembly process. Just add
   //                        matrices and vectors
   panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
   panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,linObjFactory);
   ae_tm.buildObjects(builder);

   // assemble linear system
   /////////////////////////////////////////////////////////////

   // build linear algebra objects: Ghost is for parallel assembly, it contains
   //                               local element contributions summed, the global IDs
   //                               are not unique. The non-ghosted or "global"
   //                               container will contain the sum over all processors
   //                               of the ghosted objects. The global indices are unique.
   RCP<panzer::LinearObjContainer> ghostCont = linObjFactory->buildGhostedLinearObjContainer();
   RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
   linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                             panzer::LinearObjContainer::F |
                                             panzer::LinearObjContainer::Mat,*ghostCont);
   linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                      panzer::LinearObjContainer::F |
                                      panzer::LinearObjContainer::Mat,*container);
   ghostCont->initialize();
   container->initialize();

   panzer::AssemblyEngineInArgs input(ghostCont,container);
   input.alpha = 0;
   input.beta = 1;

   // evaluate physics: This does both the Jacobian and residual at once
   ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

   // solve linear system
   /////////////////////////////////////////////////////////////

   // convert generic linear object container to epetra container
   RCP<panzer::EpetraLinearObjContainer> ep_container 
         = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(container);

   // Setup the linear solve: notice A is used directly 
   Epetra_LinearProblem problem(&*ep_container->A,&*ep_container->x,&*ep_container->f); 

   // build the solver
   AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres); // we don't push out dirichlet conditions
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,1000);
   solver.SetAztecOption(AZ_output,10);
   solver.SetAztecOption(AZ_precond,AZ_Jacobi);

   // solve the linear system
   solver.Iterate(1000,1e-5);

   // we have now solved for the residual correction from
   // zero in the context of a Newton solve.
   //     J*e = -r = -(f - J*0) where f = J*u
   // Therefore we have  J*e=-J*u which implies e = -u
   // thus we will scale the solution vector 
   ep_container->x->Scale(-1.0);
  
   // output data (optional)
   /////////////////////////////////////////////////////////////

   // write out linear system
   if(false) {
      EpetraExt::RowMatrixToMatrixMarketFile("a_op.mm",*ep_container->A);
      EpetraExt::VectorToMatrixMarketFile("x_vec.mm",*ep_container->x);
      EpetraExt::VectorToMatrixMarketFile("b_vec.mm",*ep_container->f);
   }

   // write out solution
   if(true) {
      // redistribute solution vector to ghosted vector
      linObjFactory->globalToGhostContainer(*container,*ghostCont, panzer::EpetraLinearObjContainer::X 
                                                                 | panzer::EpetraLinearObjContainer::DxDt); 

      // get X Epetra_Vector from ghosted container
      RCP<panzer::EpetraLinearObjContainer> ep_ghostCont = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ghostCont);
      panzer_stk::write_solution_data(*rcp_dynamic_cast<panzer::DOFManager<int,int> >(dofManager),*mesh,*ep_ghostCont->x);
      mesh->writeToExodus("output.exo");
   }

   // all done!
   /////////////////////////////////////////////////////////////

   return 0;
}

void testInitialization(panzer::InputPhysicsBlock& ipb,
		       std::vector<panzer::BC>& bcs)
{
   panzer::InputEquationSet ies;
   ies.name = "CurlLaplacian";
   ies.basis = "QEdge1";
   ies.integration_order = 2;
   ies.model_id = "solid";
   ies.prefix = "";
  
   ipb.physics_block_id = "4";
   ipb.eq_sets.push_back(ies);
   
   {
      std::size_t bc_id = 0;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "EFIELD";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
  		    strategy, p);
      bcs.push_back(bc);
   }    

   {
      std::size_t bc_id = 1;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "EFIELD";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
  		    strategy, p);
      bcs.push_back(bc);
   }    

   {
      std::size_t bc_id = 2;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "EFIELD";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
  		    strategy, p);
      bcs.push_back(bc);
   }    

   {
      std::size_t bc_id = 3;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "bottom";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "EFIELD";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
  		    strategy, p);
      bcs.push_back(bc);
   }    
}
