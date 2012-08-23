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
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"



namespace panzer {

  void getNodeIds(stk::mesh::EntityRank nodeRank,const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds);

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
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
    int base_cell_dimension = 2;

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);


    std::vector< Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

    for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); 
	 ++i) {

      std::vector<std::size_t> local_cell_ids;
      Intrepid::FieldContainer<double> cell_vertex_coordinates;

      panzer_stk::workset_utils::getIdsAndVertices(*mesh, element_blocks[i], local_cell_ids, 
				cell_vertex_coordinates);

      Teuchos::RCP<shards::CellTopology> topo
         = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));


      worksets.push_back(panzer::buildWorksets(element_blocks[i],topo,
					       local_cell_ids,
					       cell_vertex_coordinates,
					       ipb,
					       workset_size,
					       base_cell_dimension));
    
      TEST_EQUALITY((*worksets[i])[0].cell_vertex_coordinates(0,0,0), cell_vertex_coordinates(0,0,0));
      TEST_EQUALITY((*worksets[i])[0].cell_vertex_coordinates(2,3,1), cell_vertex_coordinates(2,3,1));

      TEST_ASSERT((*worksets[i])[0].cell_local_ids == local_cell_ids);
    }
    

    TEST_EQUALITY(worksets.size(), 2);
    TEST_EQUALITY(worksets[0]->size(), 1);
    TEST_EQUALITY(worksets[1]->size(), 1);

    TEST_EQUALITY((*worksets[0])[0].num_cells, 4);
    TEST_EQUALITY((*worksets[1])[0].num_cells, 4);
    
    TEST_EQUALITY((*worksets[0])[0].block_id, element_blocks[0]);
    TEST_EQUALITY((*worksets[1])[0].block_id, element_blocks[1]);
    
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

    int base_cell_dimension = 2;

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    std::vector<Teuchos::RCP<std::map<unsigned,panzer::Workset> > > 
      bc_worksets;
    
    for (std::vector<panzer::BC>::const_iterator bc = bcs.begin();
	 bc != bcs.end(); ++bc) {
      
      std::vector<stk::mesh::Entity*> sideEntities; 
      mesh->getMySides(bc->sidesetID(),bc->elementBlockID(),sideEntities);
   
      
      std::vector<stk::mesh::Entity*> elements;
      std::vector<std::size_t> local_cell_ids;
      std::vector<std::size_t> local_side_ids;
      panzer_stk::workset_utils::getSideElements(*mesh, bc->elementBlockID(),
		      sideEntities,local_side_ids,elements);

      Intrepid::FieldContainer<double> vertices;
      vertices.resize(elements.size(),4,dim);  
      
      // loop over elements of this block
      for(std::size_t elm=0;elm<elements.size();++elm) {
	std::vector<stk::mesh::EntityId> nodes;
	stk::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh->elementLocalId(element));
	getNodeIds(mesh->getNodeRank(),element,nodes);
	
	TEUCHOS_ASSERT(nodes.size()==4);
	
	for(std::size_t v=0;v<nodes.size();++v) {
	  const double * coord = mesh->getNodeCoordinates(nodes[v]);
          
	  for(unsigned d=0;d<dim;++d) 
	    vertices(elm,v,d) = coord[d]; 
	}
      }
      
      Teuchos::RCP<std::map<unsigned,panzer::Workset> > workset = 
	buildBCWorkset(*bc,topo, local_cell_ids, local_side_ids,
		       vertices, ipb, base_cell_dimension);
      
      bc_worksets.push_back(workset);
    }
    
    
    TEST_EQUALITY(bc_worksets[0]->size(), 1);
    TEST_EQUALITY(bc_worksets[1]->size(), 1);
    TEST_EQUALITY(bc_worksets[2]->size(), 1);

    std::map<unsigned,panzer::Workset>& workset_bc0 = *bc_worksets[0];
    TEST_EQUALITY(workset_bc0[3].num_cells, 4);
    TEST_EQUALITY(workset_bc0[3].block_id, "eblock-0_0");
    std::map<unsigned,panzer::Workset>& workset_bc1 = *bc_worksets[1];
    TEST_EQUALITY(workset_bc1[1].num_cells, 4);
    TEST_EQUALITY(workset_bc1[1].block_id, "eblock-1_0");
    std::map<unsigned,panzer::Workset>& workset_bc2 = *bc_worksets[2];
    TEST_EQUALITY(workset_bc2[2].num_cells, 6);
    TEST_EQUALITY(workset_bc2[2].block_id, "eblock-1_0");

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

    TEST_ASSERT(workset_bc0[3].cell_vertex_coordinates(cell_index,0,0)
		== 0.0);

    cell_index =
    std::distance(workset_bc2[2].cell_local_ids.begin(), 
		  std::find(workset_bc2[2].cell_local_ids.begin(),
			    workset_bc2[2].cell_local_ids.end(), 47)
		  );

    TEST_EQUALITY(workset_bc2[2].cell_vertex_coordinates(cell_index,2,0),
		  1.0);

  }


  void getNodeIds(stk::mesh::EntityRank nodeRank,const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds)
  {
    stk::mesh::PairIterRelation nodeRel = element->relations(nodeRank);
    
    stk::mesh::PairIterRelation::iterator itr;
    for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
  }
  
  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Momentum";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = "fluid";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    panzer::InputEquationSet ies_2;
    ies_2.name = "Continuity";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "fluid";
    ies_2.prefix = "ION_";
    ies_2.params.set<int>("junk", 1);

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);


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
