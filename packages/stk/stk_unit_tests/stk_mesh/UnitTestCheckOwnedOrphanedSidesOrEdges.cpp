#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/MeshDiagnostics.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace
{

class MeshCheckerOwnedOrphans : public stk::unit_test_util::MeshFixture
{
protected:
  MeshCheckerOwnedOrphans()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

void fill_mesh_with_orphaned_owned_sides(stk::mesh::BulkData& bulkData)
{
  bulkData.modification_begin();
  if(bulkData.parallel_rank()==1)
  {
    stk::mesh::EntityVector nodes(4);
    for(stk::mesh::EntityId id=1;id<=4;++id)
      nodes[id-1] = bulkData.declare_node(id);

    stk::mesh::Entity face = bulkData.declare_solo_side(1, {&bulkData.mesh_meta_data().get_topology_root_part(stk::topology::QUADRILATERAL_4)});
    for(size_t i=0;i<nodes.size();++i)
      bulkData.declare_relation(face, nodes[i], i);
  }
  bulkData.modification_end();
}

TEST_F(MeshCheckerOwnedOrphans, check_mesh_without_orphaned_owned_sides)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<stk::mesh::Entity> orphanedOwnedSides = stk::mesh::get_orphaned_owned_sides(get_bulk());
    EXPECT_TRUE(orphanedOwnedSides.empty());
  }
}

TEST_F(MeshCheckerOwnedOrphans, check_mesh_with_orphaned_owned_sides)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<stk::mesh::EntityKey> > badSidesPerProc = {
      {},
      {stk::mesh::EntityKey(stk::topology::FACE_RANK, 1)}
    };

    fill_mesh_with_orphaned_owned_sides(get_bulk());
    std::vector<stk::mesh::Entity> orphanedOwnedSides = stk::mesh::get_orphaned_owned_sides(get_bulk());
    std::vector<stk::mesh::EntityKey> orphanedOwnedSidesKeys;
    for(stk::mesh::Entity entity : orphanedOwnedSides)
      orphanedOwnedSidesKeys.push_back(get_bulk().entity_key(entity));

    EXPECT_EQ(badSidesPerProc[get_bulk().parallel_rank()].size(), orphanedOwnedSides.size());
    EXPECT_EQ(badSidesPerProc[get_bulk().parallel_rank()], orphanedOwnedSidesKeys);
    EXPECT_FALSE(stk::is_true_on_all_procs(get_bulk().parallel(), orphanedOwnedSides.empty()));
    const std::vector<std::string> & orphanedOwnedMessages = get_messages_for_orphaned_owned_sides(get_bulk(), orphanedOwnedSides);
    for (const std::string & errorMessage : orphanedOwnedMessages) {
      std::cerr << errorMessage;
    }
  }
}

class MeshCheckerWithElements : public stk::unit_test_util::MeshFixture
{
protected:
  void delete_element5_on_proc0()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank()==0)
    {
      stk::mesh::EntityId elementId = 5;
      stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elementId);
      get_bulk().destroy_entity(element);
    }
    get_bulk().modification_end();
  }

  void make_an_orphaned_face_on_proc0()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank()==0)
    {
      stk::mesh::EntityIdVector ids{18, 19, 22, 23};
      stk::mesh::Entity face = get_bulk().declare_solo_side(1, {&get_bulk().mesh_meta_data().get_topology_root_part(stk::topology::QUADRILATERAL_4)});
      for(size_t i=0;i<ids.size();++i)
      {
        stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, ids[i]);
        get_bulk().declare_relation(face, node, i);
      }
    }
    get_bulk().modification_end();
  }

  void test_if_orphaned_faces_are_found()
  {
    std::vector<std::vector<stk::mesh::EntityKey> > badSidesPerProc = {
      {stk::mesh::EntityKey(stk::topology::FACE_RANK, 1)},
      {},
    };

    std::vector<stk::mesh::Entity> orphanedOwnedSides = get_orphaned_owned_sides(get_bulk());
    std::vector<stk::mesh::EntityKey> orphanedOwnedSidesKeys;
    for(stk::mesh::Entity entity : orphanedOwnedSides)
      orphanedOwnedSidesKeys.push_back(get_bulk().entity_key(entity));
    EXPECT_EQ(badSidesPerProc[get_bulk().parallel_rank()], orphanedOwnedSidesKeys);
  }
};

TEST_F(MeshCheckerWithElements, withoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh("generated:3x1x4", stk::mesh::BulkData::NO_AUTO_AURA);
    delete_element5_on_proc0();
    make_an_orphaned_face_on_proc0();
    test_if_orphaned_faces_are_found();
  }
}

TEST_F(MeshCheckerWithElements, withAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh("generated:3x1x4", stk::mesh::BulkData::AUTO_AURA);
    delete_element5_on_proc0();
    make_an_orphaned_face_on_proc0();
    test_if_orphaned_faces_are_found();
  }
}

}
