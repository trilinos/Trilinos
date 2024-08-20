// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_PauseToAttach.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

// This unit test set tracks some particularly nasty corner cases for
// boundary conditions.

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  void testInitialzation_multiblock(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
				    std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(bc_test, basic)
  {
    using Teuchos::RCP;


    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",1);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    Teuchos::RCP<Teuchos::ParameterList> ipb =
      Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
                                 eqset_factory,
				 gd,
		   	         false,
                                 physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

   // build DOF Manager
   /////////////////////////////////////////////////////////////

   // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = eLinObjFactory;

    // setup field manager build
    /////////////////////////////////////////////////////////////

    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");

    fmb->setWorksetContainer(wkstContainer);
    fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    RCP<panzer::EpetraLinearObjContainer> eGhosted
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildGhostedLinearObjContainer());
    RCP<panzer::EpetraLinearObjContainer> eGlobal
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildLinearObjContainer());
    eLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                               panzer::EpetraLinearObjContainer::DxDt |
                                               panzer::EpetraLinearObjContainer::F |
                                               panzer::EpetraLinearObjContainer::Mat,*eGhosted);
    eLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                        panzer::EpetraLinearObjContainer::DxDt |
                                        panzer::EpetraLinearObjContainer::F |
                                        panzer::EpetraLinearObjContainer::Mat,*eGlobal);
    eGlobal->initialize();
    eGhosted->initialize();
    panzer::AssemblyEngineInArgs input(eGhosted,eGlobal);
    input.alpha = 0.0;
    input.beta = 1.0;

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
  }

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }

    // TEMPERATURE BCs
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 3.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 1;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 6.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 3;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "bottom";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 4.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }

    // ION_TEMPERATURE BCs
    {
      std::size_t bc_id = 4;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "ION_TEMPERATURE";
      std::string strategy = "Constant";
      double value = 13.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 5;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "ION_TEMPERATURE";
      std::string strategy = "Constant";
      double value = 15.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 6;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "ION_TEMPERATURE";
      std::string strategy = "Constant";
      double value = 16.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 7;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "bottom";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "ION_TEMPERATURE";
      std::string strategy = "Constant";
      double value = 14.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
  }

  // *******************************************************************************
  // *******************************************************************************
  TEUCHOS_UNIT_TEST(bc_test, multiblock)
  {
    using Teuchos::RCP;


    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",1);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    int myRank = Comm->MyPID();

    TEUCHOS_ASSERT(Comm->NumProc()==2);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation_multiblock(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
                                 eqset_factory,
				 gd,
		   	         false,
                                 physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

   // build DOF Manager
   /////////////////////////////////////////////////////////////

   // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = eLinObjFactory;

    // setup field manager build
    /////////////////////////////////////////////////////////////

    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");

    fmb->setWorksetContainer(wkstContainer);
    fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    RCP<panzer::EpetraLinearObjContainer> eGhosted
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildGhostedLinearObjContainer());
    RCP<panzer::EpetraLinearObjContainer> eGlobal
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildLinearObjContainer());
    eLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                               panzer::EpetraLinearObjContainer::DxDt |
                                               panzer::EpetraLinearObjContainer::F |
                                               panzer::EpetraLinearObjContainer::Mat,*eGhosted);
    eLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                        panzer::EpetraLinearObjContainer::DxDt |
                                        panzer::EpetraLinearObjContainer::F |
                                        panzer::EpetraLinearObjContainer::Mat,*eGlobal);
    eGhosted->initialize();
    eGlobal->initialize();
    panzer::AssemblyEngineInArgs input(eGhosted,eGlobal);
    input.alpha = 0.0;
    input.beta = 1.0;

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    const Epetra_Vector & f = *eGlobal->get_f();

    std::vector<panzer::GlobalOrdinal> GIDs;
    if(myRank==0) {
       dofManager->getElementGIDs(1,GIDs,"eblock-1_0"); // in eblock-1_0

       int gid = GIDs[3]; // top left corner at block interface
       int lid = f.Map().LID(gid);

       if(lid>=0)
          TEST_EQUALITY(f[lid],-3.0);
    }
    else if(myRank==1) {
       dofManager->getElementGIDs(0,GIDs,"eblock-0_0"); // in eblock-0_0
       int gid = GIDs[2]; // top right corner at block interface
       int lid = f.Map().LID(gid);

       if(lid>=0)
          TEST_EQUALITY(f[lid],-3.0);
    }
    else
       TEUCHOS_ASSERT(false);

    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    const Epetra_CrsMatrix & A = *eGlobal->get_A();

    if(myRank==0) {
       dofManager->getElementGIDs(1,GIDs,"eblock-1_0"); // in eblock-1_0

       panzer::GlobalOrdinal gid = GIDs[3]; // top left corner at block interface
       int lid = f.Map().LID(gid);

       if(lid>=0) {
          TEST_EQUALITY(f[lid],-3.0);

          int numEntries = 0;
          int * indices = 0;
          double * values = 0;
          A.ExtractMyRowView(lid,numEntries,values,indices);

          double sum = 0.0;
          for(int i=0;i<numEntries;i++) {
             sum += values[i];
             if(values[i]!=0.0)
             {   TEST_EQUALITY(values[i],1.0); } // diag entry should be 0
             else
             {   TEST_EQUALITY(values[i],0.0); } // all non diag entries should be 0
          }
          TEST_EQUALITY(sum,1.0); // only one entry allowed to be 1
       }
    }
    else if(myRank==1) {
       dofManager->getElementGIDs(0,GIDs,"eblock-0_0"); // in eblock-0_0
       int gid = GIDs[2]; // top right corner at block interface
       int lid = f.Map().LID(gid);

       if(lid>=0) {
          TEST_EQUALITY(f[lid],-3.0);

          int numEntries = 0;
          int * indices = 0;
          double * values = 0;
          A.ExtractMyRowView(lid,numEntries,values,indices);

          double sum = 0.0;
          for(int i=0;i<numEntries;i++) {
             sum += values[i];
             if(values[i]!=0.0)
             {   TEST_EQUALITY(values[i],1.0); } // diag entry should be 0
             else
             {   TEST_EQUALITY(values[i],0.0); } // all non diag entries should be 0
          }
          TEST_EQUALITY(sum,1.0); // only one entry allowed to be 1
       }
    }
    else {
       TEUCHOS_ASSERT(false);
    }
  }

  void testInitialzation_multiblock(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
				    std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
      p.set("Integration Order",1);
    }

    // TEMPERATURE BCs
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 3.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 1;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
		    strategy, p);
      bcs.push_back(bc);
    }
  }

} // end panzer
