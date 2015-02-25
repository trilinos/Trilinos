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

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  void getNodeIds(stk_classic::mesh::EntityRank nodeRank,const stk_classic::mesh::Entity * element,
		  std::vector<stk_classic::mesh::EntityId> & nodeIds);

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);


  TEUCHOS_UNIT_TEST(workset_builder, stk_edge)
  {
    PHX::KokkosDeviceSession session;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);  // in each block
    pl->set("Y Elements",2);  // in each block

    // This mesh will have element IDs
    //    4 5 6 7
    //    0 1 2 3

    int myRank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    panzer_stk_classic::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    mesh->writeToExodus("test.exo");

    std::vector<std::string> element_blocks;
    mesh->getElementBlockNames(element_blocks);
    const std::size_t workset_size = 20;

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      const int default_integration_order = 1;

      Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

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

    {
      std::string sideset = "vertical_0";
      Teuchos::RCP<std::vector<panzer::Workset> > worksets = panzer_stk_classic::buildWorksets(*mesh,
                                          *(panzer::findPhysicsBlock(element_blocks[0],physicsBlocks)),
                                          *(panzer::findPhysicsBlock(element_blocks[1],physicsBlocks)),
                                          sideset);
     
      if(myRank==0) {
        TEST_EQUALITY((*worksets).size(),0); // no elements on this processor
      }
      else {
        TEST_EQUALITY((*worksets).size(),1);
        TEST_EQUALITY((*worksets)[0].num_cells,2);
        TEST_EQUALITY((*worksets)[0].subcell_dim,1);
        TEST_EQUALITY((*worksets)[0].details.size(),2);
  
        // this is identical to details[0]
        TEST_EQUALITY((*worksets)[0].subcell_index, 1);
        TEST_EQUALITY((*worksets)[0].block_id, "eblock-0_0");
        TEST_EQUALITY((*worksets)[0].cell_local_ids.size(),2);
        TEST_EQUALITY((*worksets)[0].cell_local_ids[0],0);
        TEST_EQUALITY((*worksets)[0].cell_local_ids[1],2);
        TEST_EQUALITY((*worksets)[0].ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].bases.size(),2);

        TEST_EQUALITY((*worksets)[0].details[0]->subcell_index, 1);
        TEST_EQUALITY((*worksets)[0].details[0]->block_id, "eblock-0_0");
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids.size(),2);
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids[0],0);
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids[1],2);
        TEST_EQUALITY((*worksets)[0].details[0]->ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].details[0]->int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].details[0]->basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].details[0]->bases.size(),2);
  
        TEST_EQUALITY((*worksets)[0].details[1]->subcell_index, 3);
        TEST_EQUALITY((*worksets)[0].details[1]->block_id, "eblock-1_0");
        TEST_EQUALITY((*worksets)[0].details[1]->cell_local_ids[0],5);
        TEST_EQUALITY((*worksets)[0].details[1]->cell_local_ids[1],7);
        TEST_EQUALITY((*worksets)[0].details[1]->ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].details[1]->int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].details[1]->basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].details[1]->bases.size(),2);
      }
    }

    {
      std::string sideset = "vertical_0";
      Teuchos::RCP<std::vector<panzer::Workset> > worksets = panzer_stk_classic::buildWorksets(*mesh,
                                          *(panzer::findPhysicsBlock(element_blocks[1],physicsBlocks)),
                                          *(panzer::findPhysicsBlock(element_blocks[0],physicsBlocks)),
                                          sideset);
     
      if(myRank==1) {
        TEST_EQUALITY((*worksets).size(),0); // no elements on this processor
      }
      else {
        TEST_EQUALITY((*worksets).size(),1);
        TEST_EQUALITY((*worksets)[0].num_cells,2);
        TEST_EQUALITY((*worksets)[0].subcell_dim,1);
        TEST_EQUALITY((*worksets)[0].details.size(),2);
  
        // this is identical to details[0]
        TEST_EQUALITY((*worksets)[0].subcell_index, 3);
        TEST_EQUALITY((*worksets)[0].block_id, "eblock-1_0");
        TEST_EQUALITY((*worksets)[0].cell_local_ids.size(),2);
        TEST_EQUALITY((*worksets)[0].cell_local_ids[0],1);
        TEST_EQUALITY((*worksets)[0].cell_local_ids[1],3);
        TEST_EQUALITY((*worksets)[0].ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].bases.size(),2);

        TEST_EQUALITY((*worksets)[0].details[0]->subcell_index, 3);
        TEST_EQUALITY((*worksets)[0].details[0]->block_id, "eblock-1_0");
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids.size(),2);
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids[0],1);
        TEST_EQUALITY((*worksets)[0].details[0]->cell_local_ids[1],3);
        TEST_EQUALITY((*worksets)[0].details[0]->ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].details[0]->int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].details[0]->basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].details[0]->bases.size(),2);
  
        TEST_EQUALITY((*worksets)[0].details[1]->subcell_index, 1);
        TEST_EQUALITY((*worksets)[0].details[1]->block_id, "eblock-0_0");
        TEST_EQUALITY((*worksets)[0].details[1]->cell_local_ids[0],4);
        TEST_EQUALITY((*worksets)[0].details[1]->cell_local_ids[1],6);
        TEST_EQUALITY((*worksets)[0].details[1]->ir_degrees->size(),1);
        TEST_EQUALITY((*worksets)[0].details[1]->int_rules.size(),1);
        TEST_EQUALITY((*worksets)[0].details[1]->basis_names->size(),2);
        TEST_EQUALITY((*worksets)[0].details[1]->bases.size(),2);
      }
    }
    
  }

  void getNodeIds(stk_classic::mesh::EntityRank nodeRank,const stk_classic::mesh::Entity * element,
		  std::vector<stk_classic::mesh::EntityId> & nodeIds)
  {
    stk_classic::mesh::PairIterRelation nodeRel = element->relations(nodeRank);
    
    stk_classic::mesh::PairIterRelation::iterator itr;
    for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
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
      p.set("Basis Order",2);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }

    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "UX";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }    
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "UX";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }   
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "UX";
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
