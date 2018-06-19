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
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs, const int integration_order);

  void testIpMatch(const panzer::WorksetDetails& d0, const panzer::WorksetDetails& d1,
                   const index_t num_cells, Teuchos::FancyOStream& out, bool& success);

  TEUCHOS_UNIT_TEST(workset_builder, stk_edge)
  {

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

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    mesh->writeToExodus("test.exo");

    std::vector<std::string> element_blocks;
    mesh->getElementBlockNames(element_blocks);
    const std::size_t workset_size = 20;

    // Use a high integration order to test ordering of the fields in
    // IntegrationValues2.
    const int default_integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs, default_integration_order);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
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
      Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[0],physicsBlocks);
      Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[1],physicsBlocks);
      Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets = panzer_stk::buildBCWorksets(
          *mesh, pb_a->getWorksetNeeds(),pb_a->elementBlockID(),
                 pb_b->getWorksetNeeds(),pb_b->elementBlockID(), sideset);
     
      if(myRank==0) {
        TEST_EQUALITY(worksets->size(),0); // no elements on this processor
      }
      else {
        TEST_EQUALITY(worksets->size(),1);
        panzer::Workset& workset = worksets->begin()->second;

        TEST_EQUALITY((*worksets).size(),1);
        TEST_EQUALITY(workset.num_cells,2);
        TEST_EQUALITY(workset.subcell_dim,1);
        TEST_EQUALITY(workset.numDetails(),2);
  
        // this is identical to workset(0)
        TEST_EQUALITY(workset.subcell_index, 1);
        TEST_EQUALITY(workset.block_id, "eblock-0_0");
        TEST_EQUALITY(workset.cell_local_ids.size(),2);
        TEST_EQUALITY(workset.cell_local_ids[0],0);
        TEST_EQUALITY(workset.cell_local_ids[1],2);
        TEST_EQUALITY(workset.ir_degrees->size(),1);
        TEST_EQUALITY(workset.int_rules.size(),1);
        TEST_EQUALITY(workset.basis_names->size(),2);
        TEST_EQUALITY(workset.bases.size(),2);

        TEST_EQUALITY(workset(0).subcell_index, 1);
        TEST_EQUALITY(workset(0).block_id, "eblock-0_0");
        TEST_EQUALITY(workset(0).cell_local_ids.size(),2);
        TEST_EQUALITY(workset(0).cell_local_ids[0],0);
        TEST_EQUALITY(workset(0).cell_local_ids[1],2);
        TEST_EQUALITY(workset(0).ir_degrees->size(),1);
        TEST_EQUALITY(workset(0).int_rules.size(),1);
        TEST_EQUALITY(workset(0).basis_names->size(),2);
        TEST_EQUALITY(workset(0).bases.size(),2);
  
        TEST_EQUALITY(workset(1).subcell_index, 3);
        TEST_EQUALITY(workset(1).block_id, "eblock-1_0");
        TEST_EQUALITY(workset(1).cell_local_ids[0],5);
        TEST_EQUALITY(workset(1).cell_local_ids[1],7);
        TEST_EQUALITY(workset(1).ir_degrees->size(),1);
        TEST_EQUALITY(workset(1).int_rules.size(),1);
        TEST_EQUALITY(workset(1).basis_names->size(),2);
        TEST_EQUALITY(workset(1).bases.size(),2);

        testIpMatch(workset(0), workset(1), workset.num_cells, out, success);
      }
    }

    {
      std::string sideset = "vertical_0";
      Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[1],physicsBlocks);
      Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[0],physicsBlocks);
      Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets = panzer_stk::buildBCWorksets(
          *mesh, pb_a->getWorksetNeeds(),pb_a->elementBlockID(),
                 pb_b->getWorksetNeeds(),pb_b->elementBlockID(), sideset);
     
      if(myRank==1) {
        TEST_EQUALITY(worksets->size(),0); // no elements on this processor
      }
      else {
        TEST_EQUALITY(worksets->size(),1);
        panzer::Workset& workset = worksets->begin()->second;

        TEST_EQUALITY(workset.num_cells,2);
        TEST_EQUALITY(workset.subcell_dim,1);
        TEST_EQUALITY(workset.numDetails(),2);
  
        // this is identical to details[0]
        TEST_EQUALITY(workset.subcell_index, 3);
        TEST_EQUALITY(workset.block_id, "eblock-1_0");
        TEST_EQUALITY(workset.cell_local_ids.size(),2);
        TEST_EQUALITY(workset.cell_local_ids[0],1);
        TEST_EQUALITY(workset.cell_local_ids[1],3);
        TEST_EQUALITY(workset.ir_degrees->size(),1);
        TEST_EQUALITY(workset.int_rules.size(),1);
        TEST_EQUALITY(workset.basis_names->size(),2);
        TEST_EQUALITY(workset.bases.size(),2);

        TEST_EQUALITY(workset(0).subcell_index, 3);
        TEST_EQUALITY(workset(0).block_id, "eblock-1_0");
        TEST_EQUALITY(workset(0).cell_local_ids.size(),2);
        TEST_EQUALITY(workset(0).cell_local_ids[0],1);
        TEST_EQUALITY(workset(0).cell_local_ids[1],3);
        TEST_EQUALITY(workset(0).ir_degrees->size(),1);
        TEST_EQUALITY(workset(0).int_rules.size(),1);
        TEST_EQUALITY(workset(0).basis_names->size(),2);
        TEST_EQUALITY(workset(0).bases.size(),2);
  
        TEST_EQUALITY(workset(1).subcell_index, 1);
        TEST_EQUALITY(workset(1).block_id, "eblock-0_0");
        TEST_EQUALITY(workset(1).cell_local_ids[0],4);
        TEST_EQUALITY(workset(1).cell_local_ids[1],6);
        TEST_EQUALITY(workset(1).ir_degrees->size(),1);
        TEST_EQUALITY(workset(1).int_rules.size(),1);
        TEST_EQUALITY(workset(1).basis_names->size(),2);
        TEST_EQUALITY(workset(1).bases.size(),2);

        testIpMatch(workset(0), workset(1), workset.num_cells, out, success);
      }
    }
    
  }

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs, const int integration_order)
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
      p.set("Integration Order",integration_order);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",integration_order);
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

  void testIpMatch(const panzer::WorksetDetails& d0, const panzer::WorksetDetails& d1,
                   const index_t num_cells, Teuchos::FancyOStream& out, bool& success)
  {
#define TED01(m) TEST_EQUALITY(d0.m, d1.m)
    TED01(int_rules.size());
    for (std::size_t iri = 0; iri < d0.int_rules.size(); ++iri) {
      const std::size_t num_ip = d0.int_rules[iri]->cub_points.extent(0),
        num_dim = d0.int_rules[iri]->cub_points.extent(1);
      for (index_t cell = 0; cell < num_cells; ++cell)
        for (std::size_t ip = 0; ip < num_ip; ++ip)
          for (std::size_t dim = 0; dim < num_dim; ++dim)
            TED01(int_rules[iri]->ip_coordinates(cell, ip, dim));
    }
#undef TED01
  }

}
