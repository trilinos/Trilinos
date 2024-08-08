#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_util/parallel/DistributedIndex.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/MeshDiagnostics.hpp"

namespace
{

class MeshChecker : public stk::unit_test_util::MeshFixture
{
public:
};

void test_non_unique_key_results(const stk::mesh::BulkData& bulkData, const std::vector<std::vector<stk::mesh::EntityKeyProc>>& gold_values)
{
  std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(bulkData);
  EXPECT_EQ(gold_values[bulkData.parallel_rank()].size(), badKeyProcs.size());
  EXPECT_EQ(gold_values[bulkData.parallel_rank()], badKeyProcs);
  EXPECT_THROW(stk::mesh::throw_if_any_proc_has_false(bulkData.parallel(), badKeyProcs.empty()), std::logic_error);
}

void test_keys_are_unique(const stk::mesh::BulkData& bulkData)
{
  std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(bulkData);
  EXPECT_TRUE(badKeyProcs.empty());
  EXPECT_NO_THROW(stk::mesh::throw_if_any_proc_has_false(bulkData.parallel(), badKeyProcs.empty()));
}

TEST_F(MeshChecker, check_throw_on_non_unique_ids)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0}},
    };

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().modification_begin();
    get_bulk().declare_node(1);
    get_bulk().modification_end();

    test_non_unique_key_results(get_bulk(), gold_values);
  }
}

TEST_F(MeshChecker, check_no_throw_on_shared_nodes)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().modification_begin();
    stk::mesh::Entity node1 = get_bulk().declare_node(1);
    get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
    get_bulk().modification_end();

    test_keys_are_unique(get_bulk());
  }
}

TEST_F(MeshChecker, check_throw_on_partial_shared_nodes)
{
  if(stk::parallel_machine_size(get_comm())==3)
  {
    std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
       { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
    };

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().modification_begin();
    stk::mesh::Entity node1 = get_bulk().declare_node(1);
    if(get_bulk().parallel_rank()<2)
    {
      get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
    }
    get_bulk().modification_end();

    test_non_unique_key_results(get_bulk(), gold_values);
  }
}

TEST_F(MeshChecker, check_throw_on_partial_shared_nodes_on_four_procs)
{
  if(stk::parallel_machine_size(get_comm())==4)
  {
    std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2},
       { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 3}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2},
       { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 3}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
       { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
      {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
       { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
    };

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().modification_begin();
    stk::mesh::Entity node1 = get_bulk().declare_node(1);
    if(get_bulk().parallel_rank()<2)
    {
      get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
    }
    else
    {
      const int other_proc = get_bulk().parallel_rank()==2 ? 3 : 2;
      get_bulk().add_node_sharing(node1, other_proc);
    }
    get_bulk().modification_end();

    test_non_unique_key_results(get_bulk(), gold_values);
  }
}

TEST_F(MeshChecker, check_throw_on_non_unique_ids_for_elements)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
      {{ stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 1}},
      {{ stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 0}},
    };

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().modification_begin();
    stk::mesh::EntityId nodeId = static_cast<stk::mesh::EntityId>(get_bulk().parallel_rank()+1u);
    stk::mesh::declare_element(get_bulk(),get_meta().get_topology_root_part(stk::topology::PARTICLE),1u,stk::mesh::EntityIdVector{nodeId});
    get_bulk().modification_end();

    test_non_unique_key_results(get_bulk(), gold_values);
  }
}

}


