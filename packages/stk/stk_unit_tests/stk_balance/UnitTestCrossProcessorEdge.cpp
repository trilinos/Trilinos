#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>

#include "stk_tools/mesh_tools/CustomAura.hpp"

namespace
{

////////////////////////////////////////////////////////////////////////////////////////////

void test_if_graph_has_cross_processor_edge(stk::mesh::BulkData &bulk);
void create_edges_for_mesh_and_check_for_cross_processor_edge(stk::mesh::BulkData &bulk);
void test_cross_processor_edge_is_in_graph(const stk::mesh::BulkData &bulk, const std::vector<stk::balance::GraphEdge> &graphEdges);
std::vector<stk::balance::GraphEdge> make_graph_edges(stk::mesh::BulkData &bulk);
std::vector<stk::balance::GraphEdge> get_graph_edges_using_graph_settings(stk::mesh::BulkData &bulk, stk::balance::GraphCreationSettings graphSettings);
typedef std::pair<stk::mesh::Entity, stk::mesh::EntityId> Edge;
Edge get_cross_processor_edge_given_proc_id(const stk::mesh::BulkData &bulk);
bool is_edge_in_graph(const Edge &edge, const std::vector<stk::balance::GraphEdge> &graphEdges);
bool are_edges_equal(const Edge &edge, const stk::balance::GraphEdge &graphEdge);

////////////////////////////////////////////////////////////////////////////////////////////

class GraphCrossProc : public stk::unit_test_util::MeshFixture {};

TEST_F(GraphCrossProc, checkEdgeWithAura)
{
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
  test_if_graph_has_cross_processor_edge(get_bulk());
}

TEST_F(GraphCrossProc, checkEdgeWithNoAura)
{
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
  test_if_graph_has_cross_processor_edge(get_bulk());
}

void test_if_graph_has_cross_processor_edge(stk::mesh::BulkData &bulk)
{
  if(bulk.parallel_size() == 2)
    create_edges_for_mesh_and_check_for_cross_processor_edge(bulk);
}

void create_edges_for_mesh_and_check_for_cross_processor_edge(stk::mesh::BulkData &bulk)
{
  std::vector<stk::balance::GraphEdge> graphEdges = make_graph_edges(bulk);
  test_cross_processor_edge_is_in_graph(bulk, graphEdges);
}

std::vector<stk::balance::GraphEdge> make_graph_edges(stk::mesh::BulkData &bulk)
{
  stk::balance::GraphCreationSettings graphSettings;
  return get_graph_edges_using_graph_settings(bulk, graphSettings);
}

std::vector<stk::balance::GraphEdge> get_graph_edges_using_graph_settings(stk::mesh::BulkData &bulk, stk::balance::GraphCreationSettings graphSettings)
{
  stk::mesh::Ghosting *customAura = nullptr;
  if(!bulk.is_automatic_aura_on())
  {
    bulk.modification_begin();
    customAura = &bulk.create_ghosting("customAura");
    bulk.modification_end();

    stk::mesh::EntityProcVec entitiesToGhost ;
    stk::tools::fill_list_of_entities_to_send_for_aura_like_ghosting(bulk, bulk.mesh_meta_data().globally_shared_part(), entitiesToGhost);
    bulk.batch_add_to_ghosting(*customAura , entitiesToGhost);
  }

  std::vector<stk::balance::GraphEdge> graphEdges;
  Zoltan2ParallelGraph graphData;

  stk::mesh::impl::LocalIdMapper localIds(bulk, stk::topology::ELEM_RANK);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(bulk, counts);
  graphData.set_num_global_elements(counts[stk::topology::ELEM_RANK]);
  graphData.set_spatial_dim(bulk.mesh_meta_data().spatial_dimension());

  graphData.createGraphEdgesUsingNodeConnectivity(bulk,
                                                  bulk.mesh_meta_data().universal_part(),
                                                  graphSettings,
                                                  graphData.get_num_global_elements(),
                                                  graphEdges,
                                                  localIds);

  if(!bulk.is_automatic_aura_on())
  {
    bulk.modification_begin();
    bulk.destroy_ghosting(*customAura);
    bulk.modification_end();
  }
  return graphEdges;
}

void test_cross_processor_edge_is_in_graph(const stk::mesh::BulkData &bulk, const std::vector<stk::balance::GraphEdge> &graphEdges)
{
  Edge edge = get_cross_processor_edge_given_proc_id(bulk);
  ASSERT_TRUE(bulk.is_valid(edge.first));
  EXPECT_TRUE(is_edge_in_graph(edge, graphEdges));
}

bool is_edge_in_graph(const Edge &edge, const std::vector<stk::balance::GraphEdge> &graphEdges)
{
  for(const stk::balance::GraphEdge& graphEdge : graphEdges)
    if(are_edges_equal(edge, graphEdge))
      return true;
  return false;
}

bool are_edges_equal(const Edge &edge, const stk::balance::GraphEdge &graphEdge)
{
  return (graphEdge.vertex1() == edge.first && graphEdge.vertex2_id() == edge.second);
}

Edge get_cross_processor_edge_given_proc_id(const stk::mesh::BulkData &bulk)
{
  std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId>> elementIds = {{2, 3}, {3, 2}};
  stk::mesh::Entity local_element = bulk.get_entity(stk::topology::ELEM_RANK, elementIds[bulk.parallel_rank()].first);
  stk::mesh::EntityId remote_element_id = elementIds[bulk.parallel_rank()].second;
  return Edge(local_element, remote_element_id);
}

}
