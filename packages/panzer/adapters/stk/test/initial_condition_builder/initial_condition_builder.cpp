// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_GlobalData.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_STKClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Panzer_InitialCondition_Builder.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

#include <vector>
#include <map>
#include <string>

namespace panzer {

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(initial_condition_builder, exodus_restart)
  {
    using Teuchos::RCP;

    panzer_stk::SquareQuadMeshFactory mesh_factory;
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 20;

    panzer::FieldManagerBuilder<int,int> fmb;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup physic blocks
    /////////////////////////////////////////////
    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    {
       std::map<std::string,panzer::InputPhysicsBlock> 
             physics_id_to_input_physics_blocks;

       testInitialzation(ipb, bcs);

       std::map<std::string,std::string> block_ids_to_physics_ids;
       block_ids_to_physics_ids["eblock-0_0"] = "test physics";
       block_ids_to_physics_ids["eblock-1_0"] = "test physics";

       std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
       block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
       block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
    
       physics_id_to_input_physics_blocks["test physics"] = ipb;

       Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

       panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                  block_ids_to_cell_topo,
                                  physics_id_to_input_physics_blocks,
                                  2,workset_size,
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
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physics_blocks,workset_size));

    // get vector of element blocks
    std::vector<std::string> elementBlocks;
    mesh->getElementBlockNames(elementBlocks);
 
    // build volume worksets from container
    std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
    panzer::getSideWorksetsFromContainer(*wkstContainer,bcs,bc_worksets);

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager 
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    Teuchos::RCP<const panzer::UniqueGlobalIndexerFactory<int,int,int,int> > indexerFactory
          = Teuchos::rcp(new panzer::DOFManagerFactory<int,int>);
    const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
          = indexerFactory->buildUniqueGlobalIndexer(MPI_COMM_WORLD,physics_blocks,conn_manager);

    // and linear object factory
    Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > elof 
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(comm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup field manager builder
    /////////////////////////////////////////////
      
    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
    user_app::STKModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");
    user_data.sublist("Panzer Data").set("Mesh", mesh);
    user_data.sublist("Panzer Data").set("DOF Manager", dofManager);
    user_data.sublist("Panzer Data").set("Linear Object Factory", lof);

    fmb.setupVolumeFieldManagers(*wkstContainer,physics_blocks,cm_factory,closure_models,*elof,user_data);
    fmb.setupBCFieldManagers(*wkstContainer,bcs,physics_blocks,eqset_factory,cm_factory,bc_factory,closure_models,*elof,user_data);

    Teuchos::ParameterList ic_closure_models("Initial Conditions");
    ic_closure_models.sublist("eblock-0_0").sublist("TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-0_0").sublist("ION_TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-1_0").sublist("TEMPERATURE").set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-1_0").sublist("ION_TEMPERATURE").set<double>("Value",3.0);    

    std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > phx_ic_field_managers;
    panzer::setupInitialConditionFieldManagers(*wkstContainer,
					       physics_blocks,
					       cm_factory,
					       ic_closure_models,
					       *elof,
					       user_data,
					       true,
					       "initial_condition_test",
					       phx_ic_field_managers);

    
    Teuchos::RCP<panzer::LinearObjContainer> loc = elof->buildLinearObjContainer();
    elof->initializeContainer(panzer::EpetraLinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::EpetraLinearObjContainer> eloc = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
    eloc->get_x()->PutScalar(0.0);
    panzer::evaluateInitialCondition(fmb.getWorksets(), phx_ic_field_managers, loc, 0.0);
    
    Teuchos::RCP<Epetra_Vector> x = eloc->get_x();
    for (int i=0; i < x->MyLength(); ++i)
      TEST_FLOATING_EQUALITY((*x)[i], 3.0, 1.0e-10);

  }

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = "solid";
    ies_1.prefix = "";

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "ion solid";
    ies_2.prefix = "ION_";

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);


    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
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
    {
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
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

}
