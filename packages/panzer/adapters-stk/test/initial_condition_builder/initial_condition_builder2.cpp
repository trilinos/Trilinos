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
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_STKClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Panzer_InitialCondition_Builder.hpp"

#include <vector>
#include <map>
#include <string>

namespace panzer {

  void testInitialzation_blockStructure(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(initial_condition_builder2, block_structure)
  {
    using Teuchos::RCP;


    panzer_stk::STK_ExodusReaderFactory mesh_factory;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 20;

    panzer::FieldManagerBuilder fmb;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("File Name","block-decomp.exo");
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
       mesh->writeToExodus("initial_condition_builder.exo");
    }

    // setup physic blocks
    /////////////////////////////////////////////
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    {
       testInitialzation_blockStructure(ipb, bcs);

       std::map<std::string,std::string> block_ids_to_physics_ids;
       block_ids_to_physics_ids["eblock-0_0"] = "PB A";
       block_ids_to_physics_ids["eblock-1_0"] = "PB B";

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
                                  physics_blocks);
    }

    // setup worksets
    /////////////////////////////////////////////

    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physics_blocks.size();i++)
      wkstContainer->setNeeds(physics_blocks[i]->elementBlockID(),physics_blocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // get vector of element blocks
    std::vector<std::string> elementBlocks;
    mesh->getElementBlockNames(elementBlocks);

    // build volume worksets from container
    std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
    panzer::getSideWorksetsFromContainer(*wkstContainer,bcs,bc_worksets);

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager> conn_manager
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    Teuchos::RCP<const panzer::GlobalIndexerFactory > indexerFactory
          = Teuchos::rcp(new panzer::DOFManagerFactory);
    const Teuchos::RCP<panzer::GlobalIndexer> dofManager
          = indexerFactory->buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physics_blocks,conn_manager);

    // and linear object factory
    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > elof
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup field manager builder
    /////////////////////////////////////////////

    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::STKModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("SOURCE_ELECTRON_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");
    user_data.sublist("Panzer Data").set("Mesh", mesh);
    user_data.sublist("Panzer Data").set("DOF Manager", dofManager);
    user_data.sublist("Panzer Data").set("Linear Object Factory", lof);

    fmb.setWorksetContainer(wkstContainer);
    fmb.setupVolumeFieldManagers(physics_blocks,cm_factory,closure_models,*elof,user_data);
    fmb.setupBCFieldManagers(bcs,physics_blocks,*eqset_factory,cm_factory,bc_factory,closure_models,*elof,user_data);

    Teuchos::ParameterList ic_closure_models("Initial Conditions");
    ic_closure_models.sublist("eblock-0_0").sublist("TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-0_0").sublist("ELECTRON_TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-1_0").sublist("TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-1_0").sublist("ION_TEMPERATURE").set<double>("Value",3.0);

    std::map<std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > phx_ic_field_managers;
    panzer::setupInitialConditionFieldManagers(*wkstContainer,
					       physics_blocks,
					       cm_factory,
					       ic_closure_models,
					       *elof,
					       user_data,
					       true,
					       "initial_condition_test",
					       phx_ic_field_managers);


    Teuchos::RCP<panzer::LinearObjContainer> loc  = elof->buildLinearObjContainer();
    elof->initializeContainer(panzer::EpetraLinearObjContainer::X,*loc);

    Teuchos::RCP<panzer::EpetraLinearObjContainer> eloc = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
    eloc->get_x()->PutScalar(0.0);

    panzer::evaluateInitialCondition(*wkstContainer, phx_ic_field_managers, loc, *elof, 0.0);

    Teuchos::RCP<Epetra_Vector> x = eloc->get_x();
    out << x->GlobalLength() << " " << x->MyLength() << std::endl;
    for (int i=0; i < x->MyLength(); ++i)
      TEST_FLOATING_EQUALITY((*x)[i], 3.0, 1.0e-10);

  }

  void testInitialzation_blockStructure(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& /* bcs */)
  {
    // Physics block
    Teuchos::ParameterList& physics_block_a = ipb->sublist("PB A");
    {
      Teuchos::ParameterList& p = physics_block_a.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block_a.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ELECTRON_");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }

    Teuchos::ParameterList& physics_block_b = ipb->sublist("PB B");
    {
      Teuchos::ParameterList& p = physics_block_b.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block_b.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
  }

}
