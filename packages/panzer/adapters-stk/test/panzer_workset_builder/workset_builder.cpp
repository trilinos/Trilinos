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
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);


  TEUCHOS_UNIT_TEST(workset_builder, volume)
  {

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);  // in each block
    pl->set("Y Elements",2);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    if(mesh->isWritable())
      mesh->writeToExodus("blocked_mesh.exo");

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

    std::vector< Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

    for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); ++i) {

      std::vector<std::size_t> local_cell_ids;
      Kokkos::DynRankView<double,PHX::Device> cell_node_coordinates;

      panzer_stk::workset_utils::getIdsAndNodes(*mesh, element_blocks[i], local_cell_ids,
				cell_node_coordinates);

      Teuchos::RCP<shards::CellTopology> topo
         = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      Teuchos::RCP<const panzer::PhysicsBlock> pb = panzer::findPhysicsBlock(element_blocks[i],physicsBlocks);
      worksets.push_back(panzer::buildWorksets(pb->getWorksetNeeds(),pb->elementBlockID(),
					       local_cell_ids,
					       cell_node_coordinates));

      auto& cur_workset = (*worksets[i])[0];
      auto workset_cell_node_coordinates_view = cur_workset.cell_node_coordinates.get_view();
      auto workset_cell_node_coordinates_h = Kokkos::create_mirror_view(workset_cell_node_coordinates_view);
      Kokkos::deep_copy(workset_cell_node_coordinates_h, workset_cell_node_coordinates_view);

      auto cell_node_coordinates_h = Kokkos::create_mirror_view(cell_node_coordinates);
      Kokkos::deep_copy(cell_node_coordinates_h, cell_node_coordinates);

      TEST_EQUALITY(workset_cell_node_coordinates_h(0,0,0), cell_node_coordinates_h(0,0,0));
      TEST_EQUALITY(workset_cell_node_coordinates_h(2,3,1), cell_node_coordinates_h(2,3,1));

      TEST_ASSERT(cur_workset.cell_local_ids==local_cell_ids);

      workset_cell_node_coordinates_view = cur_workset(0).cell_node_coordinates.get_view();
      Kokkos::deep_copy(workset_cell_node_coordinates_h, workset_cell_node_coordinates_view);

      TEST_EQUALITY(workset_cell_node_coordinates_h(0,0,0), cell_node_coordinates_h(0,0,0));
      TEST_EQUALITY(workset_cell_node_coordinates_h(2,3,1), cell_node_coordinates_h(2,3,1));
    }


    TEST_EQUALITY(worksets.size(), 2);
    TEST_EQUALITY(worksets[0]->size(), 1);
    TEST_EQUALITY(worksets[1]->size(), 1);

    TEST_EQUALITY((*worksets[0])[0].num_cells, 4);
    TEST_EQUALITY((*worksets[1])[0].num_cells, 4);

    TEST_EQUALITY((*worksets[0])[0].block_id, element_blocks[0]);
    TEST_EQUALITY((*worksets[1])[0].block_id, element_blocks[1]);

  }

  TEUCHOS_UNIT_TEST(workset_builder, edge)
  {

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);  // in each block
    pl->set("Y Elements",2);  // in each block

    // This mesh will have element IDs
    //    4 5 6 7
    //    0 1 2 3

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

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

    Teuchos::RCP<std::vector<panzer::Workset> > worksets;

    {
      std::vector<std::size_t> local_cell_ids_a, local_cell_ids_b;
      std::vector<std::size_t> local_side_ids_a, local_side_ids_b;

      local_cell_ids_a.push_back(1);
      local_cell_ids_a.push_back(5);
      local_cell_ids_b.push_back(2);
      local_cell_ids_b.push_back(6);

      local_side_ids_a.push_back(1);
      local_side_ids_a.push_back(1);
      local_side_ids_b.push_back(3);
      local_side_ids_b.push_back(3);

      Kokkos::DynRankView<double,PHX::Device> cell_node_coordinates_a, cell_node_coordinates_b;
      mesh->getElementNodes(local_cell_ids_a,cell_node_coordinates_a);
      mesh->getElementNodes(local_cell_ids_b,cell_node_coordinates_b);

      Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[0],physicsBlocks);
      Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[1],physicsBlocks);
      worksets = panzer::buildEdgeWorksets(
                                 pb_a->getWorksetNeeds(),pb_a->elementBlockID(),
 	                               local_cell_ids_a,
				                         local_side_ids_a,
				                         cell_node_coordinates_a,
                                 pb_b->getWorksetNeeds(),pb_b->elementBlockID(),
			                           local_cell_ids_b,
			                           local_side_ids_b,
			                           cell_node_coordinates_b);


      TEST_EQUALITY((*worksets).size(),1);
      TEST_EQUALITY((*worksets)[0].num_cells,2);
      TEST_EQUALITY((*worksets)[0].subcell_dim,1);

      auto cell_node_coordinates_0_view = (*worksets)[0].cell_node_coordinates.get_view();
      auto cell_node_coordinates_0_h = Kokkos::create_mirror_view(cell_node_coordinates_0_view);
      Kokkos::deep_copy(cell_node_coordinates_0_h, cell_node_coordinates_0_view);

      auto cell_node_coordinates_a_h = Kokkos::create_mirror_view(cell_node_coordinates_a);
      Kokkos::deep_copy(cell_node_coordinates_a_h, cell_node_coordinates_a);

      // this is identical to (*worksets)[0](0)
      TEST_EQUALITY(cell_node_coordinates_0_h(0,0,0), cell_node_coordinates_a_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_0_h(1,3,1), cell_node_coordinates_a_h(1,3,1));
      TEST_EQUALITY((*worksets)[0].subcell_index, 1);
      TEST_EQUALITY((*worksets)[0].block_id, "eblock-0_0");
      TEST_EQUALITY((*worksets)[0].cell_local_ids.size(),2);
      TEST_EQUALITY((*worksets)[0].cell_local_ids[0],1);
      TEST_EQUALITY((*worksets)[0].cell_local_ids[1],5);
      TEST_EQUALITY((*worksets)[0].ir_degrees->size(),1);
      TEST_EQUALITY((*worksets)[0].int_rules.size(),1);
      TEST_EQUALITY((*worksets)[0].basis_names->size(),2);
      TEST_EQUALITY((*worksets)[0].bases.size(),2);

      cell_node_coordinates_0_view = (*worksets)[0](0).cell_node_coordinates.get_view();
      cell_node_coordinates_0_h = Kokkos::create_mirror_view(cell_node_coordinates_0_view);
      Kokkos::deep_copy(cell_node_coordinates_0_h, cell_node_coordinates_0_view);

      TEST_EQUALITY(cell_node_coordinates_0_h(0,0,0), cell_node_coordinates_a_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_0_h(1,3,1), cell_node_coordinates_a_h(1,3,1));
      TEST_EQUALITY((*worksets)[0](0).subcell_index, 1);
      TEST_EQUALITY((*worksets)[0](0).block_id, "eblock-0_0");
      TEST_EQUALITY((*worksets)[0](0).cell_local_ids.size(),2);
      TEST_EQUALITY((*worksets)[0](0).cell_local_ids[0],1);
      TEST_EQUALITY((*worksets)[0](0).cell_local_ids[1],5);
      TEST_EQUALITY((*worksets)[0](0).ir_degrees->size(),1);
      TEST_EQUALITY((*worksets)[0](0).int_rules.size(),1);
      TEST_EQUALITY((*worksets)[0](0).basis_names->size(),2);
      TEST_EQUALITY((*worksets)[0](0).bases.size(),2);

      auto cell_node_coordinates_1_view = (*worksets)[0](1).cell_node_coordinates.get_view();
      auto cell_node_coordinates_1_h = Kokkos::create_mirror_view(cell_node_coordinates_1_view);
      Kokkos::deep_copy(cell_node_coordinates_1_h, cell_node_coordinates_1_view);

      auto cell_node_coordinates_b_h = Kokkos::create_mirror_view(cell_node_coordinates_b);
      Kokkos::deep_copy(cell_node_coordinates_b_h, cell_node_coordinates_b);

      TEST_EQUALITY(cell_node_coordinates_1_h(0,0,0), cell_node_coordinates_b_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_1_h(1,3,1), cell_node_coordinates_b_h(1,3,1));
      TEST_EQUALITY((*worksets)[0](1).subcell_index, 3);
      TEST_EQUALITY((*worksets)[0](1).block_id, "eblock-1_0");
      TEST_EQUALITY((*worksets)[0](1).cell_local_ids[0],2);
      TEST_EQUALITY((*worksets)[0](1).cell_local_ids[1],6);
      TEST_EQUALITY((*worksets)[0](1).ir_degrees->size(),1);
      TEST_EQUALITY((*worksets)[0](1).int_rules.size(),1);
      TEST_EQUALITY((*worksets)[0](1).basis_names->size(),2);
      TEST_EQUALITY((*worksets)[0](1).bases.size(),2);
    }

  }

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

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

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
      Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[0],physicsBlocks);
      Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[1],physicsBlocks);
      Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets = panzer_stk::buildBCWorksets(
          *mesh, pb_a->getWorksetNeeds(),pb_a->elementBlockID(),
                 pb_b->getWorksetNeeds(),pb_b->elementBlockID(), sideset);

      std::vector<std::size_t> local_cell_ids_a, local_cell_ids_b;
      std::vector<std::size_t> local_side_ids_a, local_side_ids_b;

      local_cell_ids_a.push_back(1);
      local_cell_ids_a.push_back(5);
      local_cell_ids_b.push_back(2);
      local_cell_ids_b.push_back(6);

      local_side_ids_a.push_back(1);
      local_side_ids_a.push_back(1);
      local_side_ids_b.push_back(3);
      local_side_ids_b.push_back(3);

      Kokkos::DynRankView<double,PHX::Device> cell_node_coordinates_a, cell_node_coordinates_b;
      mesh->getElementNodes(local_cell_ids_a,cell_node_coordinates_a);
      mesh->getElementNodes(local_cell_ids_b,cell_node_coordinates_b);

      TEST_EQUALITY(worksets->size(),1);
      panzer::Workset& workset = worksets->begin()->second;

      TEST_EQUALITY(workset.num_cells,2);
      TEST_EQUALITY(workset.subcell_dim,1);
      TEST_EQUALITY(workset.numDetails(),2);

      auto cell_node_coordinates_0_view = workset.cell_node_coordinates.get_view();
      auto cell_node_coordinates_0_h = Kokkos::create_mirror_view(cell_node_coordinates_0_view);
      Kokkos::deep_copy(cell_node_coordinates_0_h, cell_node_coordinates_0_view);

      auto cell_node_coordinates_a_h = Kokkos::create_mirror_view(cell_node_coordinates_a);
      Kokkos::deep_copy(cell_node_coordinates_a_h, cell_node_coordinates_a);

      // this is identical to workset(0)
      TEST_EQUALITY(cell_node_coordinates_0_h(0,0,0), cell_node_coordinates_a_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_0_h(1,3,1), cell_node_coordinates_a_h(1,3,1));
      TEST_EQUALITY(workset.subcell_index, 1);
      TEST_EQUALITY(workset.block_id, "eblock-0_0");
      TEST_EQUALITY(workset.cell_local_ids.size(),2);
      TEST_EQUALITY(workset.cell_local_ids[0],1);
      TEST_EQUALITY(workset.cell_local_ids[1],5);
      TEST_EQUALITY(workset.ir_degrees->size(),1);
      TEST_EQUALITY(workset.int_rules.size(),1);
      TEST_EQUALITY(workset.basis_names->size(),2);
      TEST_EQUALITY(workset.bases.size(),2);

      cell_node_coordinates_0_view = workset(0).cell_node_coordinates.get_view();
      cell_node_coordinates_0_h = Kokkos::create_mirror_view(cell_node_coordinates_0_view);
      Kokkos::deep_copy(cell_node_coordinates_0_h, cell_node_coordinates_0_view);

      TEST_EQUALITY(cell_node_coordinates_0_h(0,0,0), cell_node_coordinates_a_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_0_h(1,3,1), cell_node_coordinates_a_h(1,3,1));
      TEST_EQUALITY(workset(0).subcell_index, 1);
      TEST_EQUALITY(workset(0).block_id, "eblock-0_0");
      TEST_EQUALITY(workset(0).cell_local_ids.size(),2);
      TEST_EQUALITY(workset(0).cell_local_ids[0],1);
      TEST_EQUALITY(workset(0).cell_local_ids[1],5);
      TEST_EQUALITY(workset(0).ir_degrees->size(),1);
      TEST_EQUALITY(workset(0).int_rules.size(),1);
      TEST_EQUALITY(workset(0).basis_names->size(),2);
      TEST_EQUALITY(workset(0).bases.size(),2);

      auto cell_node_coordinates_1_view = workset(1).cell_node_coordinates.get_view();
      auto cell_node_coordinates_1_h = Kokkos::create_mirror_view(cell_node_coordinates_1_view);
      Kokkos::deep_copy(cell_node_coordinates_1_h, cell_node_coordinates_1_view);

      auto cell_node_coordinates_b_h = Kokkos::create_mirror_view(cell_node_coordinates_b);
      Kokkos::deep_copy(cell_node_coordinates_b_h, cell_node_coordinates_b);

      TEST_EQUALITY(cell_node_coordinates_1_h(0,0,0), cell_node_coordinates_b_h(0,0,0));
      TEST_EQUALITY(cell_node_coordinates_1_h(1,3,1), cell_node_coordinates_b_h(1,3,1));
      TEST_EQUALITY(workset(1).subcell_index, 3);
      TEST_EQUALITY(workset(1).block_id, "eblock-1_0");
      TEST_EQUALITY(workset(1).cell_local_ids[0],2);
      TEST_EQUALITY(workset(1).cell_local_ids[1],6);
      TEST_EQUALITY(workset(1).ir_degrees->size(),1);
      TEST_EQUALITY(workset(1).int_rules.size(),1);
      TEST_EQUALITY(workset(1).basis_names->size(),2);
      TEST_EQUALITY(workset(1).bases.size(),2);
    }

  }

  TEUCHOS_UNIT_TEST(workset_builder, sidesets)
  {
    using Teuchos::RCP;


    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);

    Teuchos::RCP<shards::CellTopology> topo
       = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    unsigned dim = mesh->getDimension();

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      const int workset_size = 20;
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

    std::vector<Teuchos::RCP<std::map<unsigned,panzer::Workset> > >
      bc_worksets;

    for (const auto& bc : bcs) {
      std::vector<stk::mesh::Entity> sideEntities;
      mesh->getMySides(bc.sidesetID(),bc.elementBlockID(),sideEntities);

      std::vector<stk::mesh::Entity> elements;
      std::vector<std::size_t> local_cell_ids;
      std::vector<std::size_t> local_side_ids;
      panzer_stk::workset_utils::getSideElements(*mesh, bc.elementBlockID(),
		      sideEntities,local_side_ids,elements);

      Kokkos::DynRankView<double,PHX::Device> nodes =
	      Kokkos::createDynRankView(nodes,"nodes",elements.size(),4,dim);
      auto nodes_h = Kokkos::create_mirror_view(nodes);

      // loop over elements of this block
      for(std::size_t elm=0;elm<elements.size();++elm) {
	      std::vector<stk::mesh::EntityId> elem_nodes;
	      stk::mesh::Entity element = elements[elm];

	      local_cell_ids.push_back(mesh->elementLocalId(element));
        mesh->getNodeIdsForElement(element,elem_nodes);

	      TEUCHOS_ASSERT(elem_nodes.size()==4);

	      for(std::size_t v=0;v<elem_nodes.size();++v) {
	        const double * coord = mesh->getNodeCoordinates(elem_nodes[v]);

	        for(unsigned d=0;d<dim;++d)
	          nodes_h(elm,v,d) = coord[d];
	      }
      }

      Kokkos::deep_copy(nodes, nodes_h);

      const auto pb = panzer::findPhysicsBlock(bc.elementBlockID(),physicsBlocks);
      const auto workset = buildBCWorkset(pb->getWorksetNeeds(),bc.elementBlockID(), local_cell_ids, local_side_ids, nodes);

      bc_worksets.push_back(workset);
    }


    TEST_EQUALITY(bc_worksets[0]->size(), 1);
    TEST_EQUALITY(bc_worksets[1]->size(), 1);
    TEST_EQUALITY(bc_worksets[2]->size(), 1);

    std::map<unsigned,panzer::Workset>& workset_bc0 = *bc_worksets[0];
    TEST_EQUALITY(workset_bc0[3].num_cells, 4);
    TEST_EQUALITY(workset_bc0[3].block_id, "eblock-0_0");
    TEST_EQUALITY(workset_bc0[3].numDetails(), 1);
    std::map<unsigned,panzer::Workset>& workset_bc1 = *bc_worksets[1];
    TEST_EQUALITY(workset_bc1[1].num_cells, 4);
    TEST_EQUALITY(workset_bc1[1].block_id, "eblock-1_0");
    TEST_EQUALITY(workset_bc1[1].numDetails(), 1);
    std::map<unsigned,panzer::Workset>& workset_bc2 = *bc_worksets[2];
    TEST_EQUALITY(workset_bc2[2].num_cells, 6);
    TEST_EQUALITY(workset_bc2[2].block_id, "eblock-1_0");
    TEST_EQUALITY(workset_bc2[2].numDetails(),1);

    // for debugging purposes
    out << "\nWORKSEST_0 - Side 3: ";
    for (std::size_t i=0; i < 4; i++ )
       out << mesh->elementGlobalId(workset_bc0[3].cell_local_ids[i]) << " ";
    out << std::endl;

    // these two loops make sure all the global IDs on the boundary are included
    // as local IDs in the workset

    for (std::size_t i=1; i < 38; i+=12 )
      TEST_ASSERT(std::find(workset_bc0[3].cell_local_ids.begin(),
			    workset_bc0[3].cell_local_ids.end(), mesh->elementLocalId(i)) !=
		  workset_bc0[3].cell_local_ids.end());

    for (std::size_t i=43; i < 49; i++)
      TEST_ASSERT(std::find(workset_bc2[2].cell_local_ids.begin(),
			    workset_bc2[2].cell_local_ids.end(), mesh->elementLocalId(i)) !=
		  workset_bc2[2].cell_local_ids.end());

    std::size_t cell_index =
      std::distance(workset_bc0[3].cell_local_ids.begin(),
		    std::find(workset_bc0[3].cell_local_ids.begin(),
			      workset_bc0[3].cell_local_ids.end(), 0)
		    );

    auto workset_bc0_3_cell_node_coordinates_v = workset_bc0[3].cell_node_coordinates.get_view();
    auto workset_bc0_3_cell_node_coordinates_h = Kokkos::create_mirror_view(workset_bc0_3_cell_node_coordinates_v);
    Kokkos::deep_copy(workset_bc0_3_cell_node_coordinates_h, workset_bc0_3_cell_node_coordinates_v);

    TEST_ASSERT(workset_bc0_3_cell_node_coordinates_h(cell_index,0,0) == 0.0);

    cell_index =
    std::distance(workset_bc2[2].cell_local_ids.begin(),
		  std::find(workset_bc2[2].cell_local_ids.begin(),
			    workset_bc2[2].cell_local_ids.end(), 47)
		  );

    auto workset_bc2_2_cell_node_coordinates_v = workset_bc2[2].cell_node_coordinates.get_view();
    auto workset_bc2_2_cell_node_coordinates_h = Kokkos::create_mirror_view(workset_bc2_2_cell_node_coordinates_v);
    Kokkos::deep_copy(workset_bc2_2_cell_node_coordinates_h, workset_bc2_2_cell_node_coordinates_v);

    TEST_EQUALITY(workset_bc2_2_cell_node_coordinates_h(cell_index,2,0), 1.0);
  }

  TEUCHOS_UNIT_TEST(workset_builder, side_element_cascade)
  {
    using Teuchos::RCP;


    // excercise subcell entities capability
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",1);
      pl->set("Y Blocks",1);
      pl->set("Z Blocks",1);
      pl->set("X Elements",6);
      pl->set("Y Elements",4);
      pl->set("Z Elements",2);

      panzer_stk::CubeTetMeshFactory factory;
      factory.setParameterList(pl);
      RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

      std::vector<stk::mesh::Entity> sideEntities;
      mesh->getMySides("left","eblock-0_0_0",sideEntities);

      std::vector<std::vector<stk::mesh::Entity> > subcells;
      panzer_stk::workset_utils::getSubcellEntities(*mesh,sideEntities,subcells);

      TEST_EQUALITY(subcells.size(),2);
      TEST_EQUALITY(subcells[0].size(),15);
      TEST_EQUALITY(subcells[1].size(),30);
    }

    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",1);
      pl->set("Y Blocks",1);
      pl->set("Z Blocks",1);
      pl->set("X Elements",1);
      pl->set("Y Elements",1);
      pl->set("Z Elements",1);

      panzer_stk::CubeTetMeshFactory factory;
      factory.setParameterList(pl);
      RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

      std::vector<stk::mesh::Entity> sideEntities;
      mesh->getMySides("left","eblock-0_0_0",sideEntities);

      std::vector<std::size_t> localSubcellDim,localSubcellIds;
      std::vector<stk::mesh::Entity> elements;

      // TOUCHING TABLE:
      // the following elements touch the side
      //    1 2 3 4 5 6 9 10 11 12       // element global IDs
      //    N E N E F F N  E  E  N       // face (F), edge (E), node (N)

      panzer_stk::workset_utils::getSideElementCascade(*mesh,"eblock-0_0_0",sideEntities,
                                                       localSubcellDim,localSubcellIds,elements);

      TEST_EQUALITY(elements.size(),30);
      TEST_EQUALITY(elements.size(),localSubcellDim.size());
      TEST_EQUALITY(elements.size(),localSubcellIds.size());

      // if you use the "TOUCHING TABLE" you can determine that
      // 2 elements touch on the face, 10 elements touch on the edge, and 18 elements touch on the node
      // of course the elements are repeated for multiple edges and nodes that touch. For instance an
      // edge touches the side, this implies that the nodes also touch that side, an element containing that
      // edge will then be included once for the edge, and twice for each node contained in that edge.

      {
        bool nodes = true; for(int i= 0;i<18;i++) nodes &= (localSubcellDim[i]==0); TEST_ASSERT(nodes);
        bool edges = true; for(int i=18;i<28;i++) edges &= (localSubcellDim[i]==1); TEST_ASSERT(edges);
        bool faces = true; for(int i=28;i<30;i++) edges &= (localSubcellDim[i]==2); TEST_ASSERT(faces);
      }

      // check that each element is assigned the correct dimension
      {
        std::set<stk::mesh::EntityId> nodeE, edgeE, faceE;
        nodeE.insert(1); nodeE.insert(2); nodeE.insert(3); nodeE.insert(4); nodeE.insert(5);
        nodeE.insert(6); nodeE.insert(9); nodeE.insert(10); nodeE.insert(11); nodeE.insert(12);

        edgeE.insert(2); edgeE.insert(4); edgeE.insert(10); edgeE.insert(11); edgeE.insert(5); edgeE.insert(6);

        faceE.insert(5); faceE.insert(6);

        bool nodes = true; for(int i= 0;i<18;i++) nodes &= (nodeE.find(mesh->elementGlobalId(elements[i]))!=nodeE.end()); TEST_ASSERT(nodes);
        bool edges = true; for(int i=18;i<28;i++) edges &= (edgeE.find(mesh->elementGlobalId(elements[i]))!=edgeE.end()); TEST_ASSERT(edges);
        bool faces = true; for(int i=28;i<30;i++) faces &= (faceE.find(mesh->elementGlobalId(elements[i]))!=faceE.end()); TEST_ASSERT(faces);
      }
    }
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
