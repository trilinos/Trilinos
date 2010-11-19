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

namespace panzer {

  template<typename Array>
  void getIdsAndVertices(const panzer_stk::STK_Interface& si,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 Array& vertices);

  void getNodeIds(const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds);

  TEUCHOS_UNIT_TEST(workset_builder, basic)
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

    panzer::InputEquationSet ies_1;
    ies_1.name = "Momentum";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = 6;
    ies_1.model_factory = "rf";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    panzer::InputEquationSet ies_2;
    ies_2.name = "Continuity";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = 6;
    ies_2.model_factory = "rf";
    ies_2.prefix = "ION_";
    ies_2.params.set<int>("junk", 1);

    panzer::InputPhysicsBlock ipb;
    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);

    std::vector< Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

    for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); 
	 ++i) {

      std::vector<std::size_t> local_cell_ids;
      Intrepid::FieldContainer<double> cell_vertex_coordinates;

      panzer::getIdsAndVertices(*mesh, element_blocks[i], local_cell_ids, 
				cell_vertex_coordinates);

      worksets.push_back(panzer::buildWorksets(element_blocks[i],
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
  
  template<typename Array>
  void getIdsAndVertices(const panzer_stk::STK_Interface& si,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 Array& vertices) {
    
    std::vector<stk::mesh::Entity*> elements;
    si.getMyElements(blockId,elements);
    
    unsigned dim = si.getDimension();
    
    vertices.resize(elements.size(),4,dim);
    
    // loop over elements of this block
    for(std::size_t elm=0;elm<elements.size();++elm) {
      std::vector<stk::mesh::EntityId> nodes;
      stk::mesh::Entity * element = elements[elm];
      
      localIds.push_back(si.elementLocalId(element));
      getNodeIds(element,nodes);
      
      TEUCHOS_ASSERT(nodes.size()==4);
      
      for(std::size_t v=0;v<nodes.size();++v) {
	const double * coord = si.getNodeCoordinates(nodes[v]);
	
	for(unsigned d=0;d<dim;++d) 
	  vertices(elm,v,d) = coord[d]; 
      }
    }
  }

  void getNodeIds(const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds)
  {
    stk::mesh::PairIterRelation nodeRel = element->relations(stk::mesh::Node);
    
    stk::mesh::PairIterRelation::iterator itr;
    for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
  }

}

