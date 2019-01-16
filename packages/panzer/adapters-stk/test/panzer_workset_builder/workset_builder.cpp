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
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_WorksetUtilities.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace {

const panzer::Workset &
getWorksetWithSide(const std::vector<panzer::Workset> & worksets,
                   const unsigned int subcell_index)
{
  unsigned int count=0;
  unsigned int idx;
  unsigned int i=0;
  for(const auto & workset : worksets){
    if(workset.getSubcellIndex() == subcell_index){
      ++count;
      idx = i;
    }
    ++i;
  }

  TEUCHOS_ASSERT(count == 1);

  return worksets[idx];

}

}


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

    auto mesh_info_rcp = panzer_stk::generateLocalMeshInfo(*mesh);
    auto & mesh_info = *mesh_info_rcp;

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

      worksets.push_back(buildWorksets(mesh_info, WorksetDescriptor(element_blocks[i])));

      std::vector<std::size_t> local_cell_ids;
      Kokkos::DynRankView<double,PHX::Device> local_cell_vertices;
      panzer_stk::workset_utils::getIdsAndVertices(*mesh, element_blocks[i], local_cell_ids, local_cell_vertices);

      const auto & workset_cell_vertices = (*worksets[i])[0].getCellVertices();
      const auto & other_workset_cell_vertices = (*worksets[i])[0](0).getCellVertices();

      TEST_EQUALITY(workset_cell_vertices(0,0,0), local_cell_vertices(0,0,0));
      TEST_EQUALITY(workset_cell_vertices(2,3,1), local_cell_vertices(2,3,1));

      TEST_EQUALITY(other_workset_cell_vertices(0,0,0), local_cell_vertices(0,0,0));
      TEST_EQUALITY(other_workset_cell_vertices(2,3,1), local_cell_vertices(2,3,1));

      for(unsigned int j=0; j<local_cell_ids.size(); ++j)
        TEST_EQUALITY((*worksets[i])[0].getLocalCellIDs()(j), local_cell_ids[j]);

    }
    

    TEST_EQUALITY(worksets.size(), 2);
    TEST_EQUALITY(worksets[0]->size(), 1);
    TEST_EQUALITY(worksets[1]->size(), 1);

    TEST_EQUALITY((*worksets[0])[0].numCells(), 4);
    TEST_EQUALITY((*worksets[1])[0].numCells(), 4);
    
    TEST_EQUALITY((*worksets[0])[0].getElementBlock(), element_blocks[0]);
    TEST_EQUALITY((*worksets[1])[0].getElementBlock(), element_blocks[1]);
    
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

      Kokkos::DynRankView<double,PHX::Device> cell_vertex_coordinates_a, cell_vertex_coordinates_b;
      mesh->getElementVertices(local_cell_ids_a,cell_vertex_coordinates_a);
      mesh->getElementVertices(local_cell_ids_b,cell_vertex_coordinates_b);

      Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[0],physicsBlocks);
      Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[1],physicsBlocks);
      worksets = panzer::buildGroupedSubcellWorksets(pb_a->elementBlockID(), "", pb_a->cellData().getCellTopology(),
                                                    local_cell_ids_a,
                                                    local_side_ids_a,
                                                    cell_vertex_coordinates_a,
                                                    pb_b->elementBlockID(), "", pb_b->cellData().getCellTopology(),
                                                    local_cell_ids_b,
                                                    local_side_ids_b,
                                                    cell_vertex_coordinates_b,
                                                    workset_size);

     
      TEST_EQUALITY((*worksets).size(),1);
      TEST_EQUALITY((*worksets)[0].numCells(),2);
      TEST_EQUALITY((*worksets)[0].getSubcellDimension(),1);
      TEST_EQUALITY((*worksets)[0].size(),2);

      // this is identical to (*worksets)[0](0)
      TEST_EQUALITY((*worksets)[0].getCellVertices()(0,0,0), cell_vertex_coordinates_a(0,0,0));
      TEST_EQUALITY((*worksets)[0].getCellVertices()(1,3,1), cell_vertex_coordinates_a(1,3,1));
      TEST_EQUALITY((*worksets)[0].getSubcellIndex(), 1);
      TEST_EQUALITY((*worksets)[0].getElementBlock(), "eblock-0_0");
      TEST_EQUALITY((*worksets)[0].getLocalCellIDs().size(),2);
      TEST_EQUALITY((*worksets)[0].getLocalCellIDs()(0),1);
      TEST_EQUALITY((*worksets)[0].getLocalCellIDs()(1),5);
      
      TEST_EQUALITY((*worksets)[0](0).getCellVertices()(0,0,0), cell_vertex_coordinates_a(0,0,0));
      TEST_EQUALITY((*worksets)[0](0).getCellVertices()(1,3,1), cell_vertex_coordinates_a(1,3,1));
      TEST_EQUALITY((*worksets)[0](0).getSubcellIndex(), 1);
      TEST_EQUALITY((*worksets)[0](0).getElementBlock(), "eblock-0_0");
      TEST_EQUALITY((*worksets)[0](0).getLocalCellIDs().size(),2);
      TEST_EQUALITY((*worksets)[0](0).getLocalCellIDs()(0),1);
      TEST_EQUALITY((*worksets)[0](0).getLocalCellIDs()(1),5);

      TEST_EQUALITY((*worksets)[0](1).getCellVertices()(0,0,0), cell_vertex_coordinates_b(0,0,0));
      TEST_EQUALITY((*worksets)[0](1).getCellVertices()(1,3,1), cell_vertex_coordinates_b(1,3,1));
      TEST_EQUALITY((*worksets)[0](1).getSubcellIndex(), 3);
      TEST_EQUALITY((*worksets)[0](1).getElementBlock(), "eblock-1_0");
      TEST_EQUALITY((*worksets)[0](1).getLocalCellIDs()(0),2);
      TEST_EQUALITY((*worksets)[0](1).getLocalCellIDs()(1),6);

      // Make sure accessing outside of the available workests fails
      TEST_THROW((*worksets)[0](2), std::logic_error);
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

    auto mesh_info_rcp = panzer_stk::generateLocalMeshInfo(*mesh);
    auto & mesh_info = *mesh_info_rcp;

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

      auto worksets = panzer::buildWorksets(mesh_info, WorksetDescriptor(pb_a->elementBlockID(), pb_b->elementBlockID(), sideset, sideset));

      const std::vector<std::size_t> local_cell_ids_a {5,1}, local_cell_ids_b {6,2};

      Kokkos::DynRankView<double,PHX::Device> cell_vertex_coordinates_a, cell_vertex_coordinates_b;
      mesh->getElementVertices(local_cell_ids_a,cell_vertex_coordinates_a);
      mesh->getElementVertices(local_cell_ids_b,cell_vertex_coordinates_b);

      TEST_EQUALITY(worksets->size(),1);
      panzer::Workset& workset = *worksets->begin();
     
      TEST_EQUALITY(workset.numCells(),2);
      TEST_EQUALITY(workset.getSubcellDimension(),1);
      TEST_EQUALITY(workset.size(),2);

      // this is identical to workset(0)
      TEST_EQUALITY(workset.getCellVertices()(0,0,0), cell_vertex_coordinates_a(0,0,0));
      TEST_EQUALITY(workset.getCellVertices()(1,3,1), cell_vertex_coordinates_a(1,3,1));
      TEST_EQUALITY(workset.getSubcellIndex(), 1);
      TEST_EQUALITY(workset.getElementBlock(), "eblock-0_0");
      TEST_EQUALITY(workset.getLocalCellIDs().size(),2);
      TEST_EQUALITY(workset.getLocalCellIDs()(0),(int) local_cell_ids_a[0]);
      TEST_EQUALITY(workset.getLocalCellIDs()(1),(int) local_cell_ids_a[1]);

      TEST_EQUALITY(workset(0).getCellVertices()(0,0,0), cell_vertex_coordinates_a(0,0,0));
      TEST_EQUALITY(workset(0).getCellVertices()(1,3,1), cell_vertex_coordinates_a(1,3,1));
      TEST_EQUALITY(workset(0).getSubcellIndex(), 1);
      TEST_EQUALITY(workset(0).getElementBlock(), "eblock-0_0");
      TEST_EQUALITY(workset(0).getLocalCellIDs().size(),2);
      TEST_EQUALITY(workset(0).getLocalCellIDs()(0),(int) local_cell_ids_a[0]);
      TEST_EQUALITY(workset(0).getLocalCellIDs()(1),(int) local_cell_ids_a[1]);
      
      TEST_EQUALITY(workset(1).getCellVertices()(0,0,0), cell_vertex_coordinates_b(0,0,0));
      TEST_EQUALITY(workset(1).getCellVertices()(1,3,1), cell_vertex_coordinates_b(1,3,1));
      TEST_EQUALITY(workset(1).getSubcellIndex(), 3);
      TEST_EQUALITY(workset(1).getElementBlock(), "eblock-1_0");
      TEST_EQUALITY(workset(1).getLocalCellIDs().size(),2);
      TEST_EQUALITY(workset(1).getLocalCellIDs()(0),(int) local_cell_ids_b[0]);
      TEST_EQUALITY(workset(1).getLocalCellIDs()(1),(int) local_cell_ids_b[1]);
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

    std::vector<Teuchos::RCP<std::vector<panzer::Workset> > > bc_worksets;
    
    for (const auto & bc : bcs) {
      
      std::vector<stk::mesh::Entity> sideEntities; 
      mesh->getMySides(bc.sidesetID(),bc.elementBlockID(),sideEntities);

      std::vector<stk::mesh::Entity> elements;
      std::vector<std::size_t> local_cell_ids;
      std::vector<std::size_t> local_side_ids;
      panzer_stk::workset_utils::getSideElements(*mesh, bc.elementBlockID(),
                                                 sideEntities,local_side_ids,elements);

      Kokkos::DynRankView<double,PHX::Device> vertices = 
          Kokkos::createDynRankView(vertices,"vertices",elements.size(),4,dim);

      // loop over elements of this block
      for(std::size_t elm=0;elm<elements.size();++elm) {
        std::vector<stk::mesh::EntityId> nodes;
        stk::mesh::Entity element = elements[elm];

        local_cell_ids.push_back(mesh->elementLocalId(element));
        mesh->getNodeIdsForElement(element,nodes);

        TEUCHOS_ASSERT(nodes.size()==4);

        for(std::size_t v=0;v<nodes.size();++v) {
          const double * coord = mesh->getNodeCoordinates(nodes[v]);

          for(unsigned d=0;d<dim;++d)
            vertices(elm,v,d) = coord[d];
        }
      }

      Teuchos::RCP<const panzer::PhysicsBlock> pb = panzer::findPhysicsBlock(bc.elementBlockID(),physicsBlocks);
      auto worksets = buildGroupedSubcellWorksets(bc.elementBlockID(), bc.sidesetID(), pb->cellData().getCellTopology(),
                                                  local_cell_ids, local_side_ids, vertices,true);

      bc_worksets.push_back(worksets);
    }
    
    
    TEST_EQUALITY(bc_worksets[0]->size(), 1);
    TEST_EQUALITY(bc_worksets[1]->size(), 1);
    TEST_EQUALITY(bc_worksets[2]->size(), 1);

    const panzer::Workset& workset_bc0_3 = getWorksetWithSide(*bc_worksets[0],3);
    TEST_EQUALITY(workset_bc0_3.numCells(), 4);
    TEST_EQUALITY(workset_bc0_3.getElementBlock(), "eblock-0_0");
    TEST_EQUALITY(workset_bc0_3.size(), 1);
    const panzer::Workset& workset_bc1_1 = getWorksetWithSide(*bc_worksets[1],1);
    TEST_EQUALITY(workset_bc1_1.numCells(), 4);
    TEST_EQUALITY(workset_bc1_1.getElementBlock(), "eblock-1_0");
    TEST_EQUALITY(workset_bc1_1.size(), 1);
    const panzer::Workset& workset_bc2_2 = getWorksetWithSide(*bc_worksets[2],2);
    TEST_EQUALITY(workset_bc2_2.numCells(), 6);
    TEST_EQUALITY(workset_bc2_2.getElementBlock(), "eblock-1_0");
    TEST_EQUALITY(workset_bc2_2.size(),1);

    // for debugging purposes
    out << "\nWORKSEST_0 - Side 3: ";
    for (std::size_t i=0; i < 4; i++ )
       out << mesh->elementGlobalId(workset_bc0_3.getLocalCellIDs()(i)) << " ";
    out << std::endl;

    // these two loops make sure all the global IDs on the boundary are included
    // as local IDs in the workset

    for (std::size_t i=1; i < 38; i+=12 ){
      bool exists = false;
      const auto local_id = mesh->elementLocalId(i);
      const auto & local_cell_ids = workset_bc0_3.getLocalCellIDs();
      for(size_t j=0; j<local_cell_ids.size(); ++j){
        if(local_cell_ids(j) == local_id){
          exists = true;
          break;
        }
      }
      TEST_ASSERT(exists);
    }

    for (std::size_t i=43; i < 49; i++){
      bool exists = false;
      const auto local_id = mesh->elementLocalId(i);
      const auto & local_cell_ids = workset_bc2_2.getLocalCellIDs();
      for(size_t j=0; j<local_cell_ids.size(); ++j){
        if(local_cell_ids(j) == local_id){
          exists = true;
          break;
        }
      }
      TEST_ASSERT(exists);
    }

    {
      bool exists = false;
      const auto & local_cell_ids = workset_bc0_3.getLocalCellIDs();
      for(size_t j=0; j<local_cell_ids.size(); ++j){
        if(local_cell_ids(j) == 0){
          TEST_ASSERT(workset_bc0_3.getCellVertices()(j,0,0) == 0.0);
          exists = true;
          break;
        }
      }
      TEUCHOS_ASSERT(exists);
    }

    {
      bool exists = false;
      const auto & local_cell_ids = workset_bc2_2.getLocalCellIDs();
      for(size_t j=0; j<local_cell_ids.size(); ++j){
        if(local_cell_ids(j) == 47){
          TEST_ASSERT(workset_bc2_2.getCellVertices()(j,2,0) == 1.0);
          exists = true;
          break;
        }
      }
      TEUCHOS_ASSERT(exists);
    }

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
