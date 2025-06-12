#include "gtest/gtest.h"
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <string>

namespace
{

TEST(TestGraph, update_graph_edges)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  //Build this mesh (2D triangles):
  //
  //   P0      P1
  //
  //              6* 
  //             / | 
  //            /E4| 
  //     3*    3*  | 
  //    / |    | \ |
  //   /E1|    |E3\| 
  // 1*   |    |  *5
  //  | \ |    |  /
  //  |E2\|    | /  
  //  |   *2  2*    
  //  |  /   
  //  | /   
  //  4*     
  std::string meshDesc = "0,1,TRI_3_2D, 1,2,3\n"
                         "0,2,TRI_3_2D, 1,4,2\n"
                         "1,3,TRI_3_2D, 2,5,3\n"
                         "1,4,TRI_3_2D, 3,5,6";
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                          .set_spatial_dimension(2)
                                          .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                          .create();
  std::vector<double> coords = {0,0, 1,-1, 1,1,
                                0,-2, 2,0, 2,2};

  stk::unit_test_util::setup_text_mesh(*bulkPtr, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  ASSERT_TRUE(bulkPtr->has_face_adjacent_element_graph());

  //The initial mesh should have 1 graph-edge on the proc boundary:
  {
    const stk::mesh::ElemElemGraph& elemGraph = bulkPtr->get_face_adjacent_element_graph();
    EXPECT_EQ(1u, elemGraph.get_parallel_graph().get_parallel_graph_info().size());
  }

  bulkPtr->modification_begin();

  if (bulkPtr->parallel_rank() == 0) {
    //delete element 1 and shared-node 3
    stk::mesh::Entity elem1 = bulkPtr->get_entity(stk::topology::ELEM_RANK,1);
    stk::mesh::Entity node3 = bulkPtr->get_entity(stk::topology::NODE_RANK,3);
    EXPECT_TRUE(bulkPtr->destroy_entity(elem1));
    EXPECT_TRUE(bulkPtr->destroy_entity(node3));
  }
  else {
    //disconnect node 5 from element 3 and then reconnect them.
    //this is just to make it modified, to trigger a call to
    //the elem-graph method 'update_graph_edges'
    stk::mesh::Entity elem3 = bulkPtr->get_entity(stk::topology::ELEM_RANK,3);
    stk::mesh::Entity node5 = bulkPtr->get_entity(stk::topology::NODE_RANK,5);
    stk::mesh::ConnectivityOrdinal ord = 1;
    bulkPtr->destroy_relation(elem3, node5, ord);
    bulkPtr->declare_relation(elem3, node5, ord);
  }

  //when elem-graph 'update_graph_edges' communicates the shared side-nodes
  //of modified elements, proc 0 will recv nodes 2 and 3 for the formerly-shared
  //side between elements 1 and 3. But node 3 is no longer valid on proc 0 and
  //this was resulting in a throw in an NGS improve_mesh scenario. The code now
  //correctly handles this by ignoring the no-longer-shared side on proc 0.
  EXPECT_NO_THROW(bulkPtr->modification_end());

  //Neither proc should have any parallel graph info (there should be no
  //graph-edges at the proc boundary).
  {
    const stk::mesh::ElemElemGraph& elemGraph = bulkPtr->get_face_adjacent_element_graph();
    EXPECT_TRUE(elemGraph.get_parallel_graph().get_parallel_graph_info().empty());
  }
}

}
