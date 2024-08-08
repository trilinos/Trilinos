#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for get_cell_topology, MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include <stk_io/WriteMesh.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_mesh/base/GetEntities.hpp>


namespace
{

class TwoHexWithInternalSideset : public stk::unit_test_util::MeshFixture
{
protected:
  void create_mesh_with_internal_side()
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    get_bulk().set_use_entity_ids_for_resolving_sharing(true);

    get_bulk().modification_begin();
    stk::mesh::Entity element1 = get_bulk().get_entity(stk::topology::ELEM_RANK,1);
    stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEM_RANK,2);

    if(get_bulk().is_valid(element1) && get_bulk().bucket(element1).owned())
      get_bulk().declare_element_side(element1, 5, stk::mesh::ConstPartVector{});

    if(get_bulk().is_valid(element2) && get_bulk().bucket(element2).owned())
      get_bulk().declare_element_side(element2, 4, stk::mesh::ConstPartVector{});

    get_bulk().modification_end();
    verify_side_is_internal_not_shared();
  }

  void verify_side_is_internal_not_shared()
  {
    stk::mesh::EntityVector sides;
    stk::mesh::Selector shared_owned = get_meta().locally_owned_part() | get_meta().globally_shared_part();
    stk::mesh::get_selected_entities(shared_owned, get_bulk().buckets(get_meta().side_rank()), sides);

    std::vector<size_t> numGoldSidesPerProc = { 1, 1 };
    ASSERT_EQ(numGoldSidesPerProc[get_bulk().parallel_rank()], sides.size());

    unsigned num_elements = get_bulk().num_elements(sides[0]);
    ASSERT_EQ(1u, num_elements);
  }

  void test_throw_on_add_node_sharing(bool use_ids_for_sharing)
  {
    get_bulk().set_use_entity_ids_for_resolving_sharing(use_ids_for_sharing);
    stk::mesh::EntityId uniqueNodeId = 100;
    get_bulk().modification_begin();
    stk::mesh::Entity node = get_bulk().declare_node(uniqueNodeId, stk::mesh::ConstPartVector{});
    int other_proc = 1 - get_bulk().parallel_rank();
    get_bulk().add_node_sharing(node, other_proc);
    get_bulk().modification_end();
  }
};

TEST_F(TwoHexWithInternalSideset, DISABLED_addNodeSharingClassic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    create_mesh_with_internal_side();
    bool use_ids_for_resolving_sharing = false;
    EXPECT_THROW(test_throw_on_add_node_sharing(use_ids_for_resolving_sharing), std::runtime_error);
  }
}

TEST_F(TwoHexWithInternalSideset, addNodeSharingCurrent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    create_mesh_with_internal_side();
    bool use_ids_for_resolving_sharing = true;
    EXPECT_NO_THROW(test_throw_on_add_node_sharing(use_ids_for_resolving_sharing));
  }
}

}
