#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "../FaceCreatorFixture.hpp"
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/MeshDiagnostics.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
namespace
{

class UnitTestRuleThreeViolation : public stk::unit_test_util::MeshFixture
{
protected:
  stk::mesh::EntityVector get_shared_nodes_of_element(stk::mesh::Entity elem)
  {
    stk::mesh::EntityVector sharedNodes;

    const stk::mesh::Entity *nodes = get_bulk().begin_nodes(elem);
    unsigned numNodes = get_bulk().num_nodes(elem);
    for(unsigned i = 0; i<numNodes; ++i)
    {
      if(get_bulk().bucket(nodes[i]).shared())
        sharedNodes.push_back(nodes[i]);
    }

    return sharedNodes;
  }

  stk::mesh::EntityVector get_shared_nodes_of_element5_on_p1()
  {
    stk::mesh::EntityVector sharedNodes;

    if(get_proc_rank() == 1)
    {
      stk::mesh::Entity elem5 = get_bulk().get_entity(stk::topology::ELEM_RANK, 5u);
      sharedNodes = get_shared_nodes_of_element(elem5);
    }
    return sharedNodes;
  }

  void delete_element1_on_p0()
  {
    if(get_proc_rank() == 0)
    {
      stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1u);

      bool isDestroyed = get_bulk().destroy_entity(elem1);
      EXPECT_TRUE(isDestroyed);
    }
  }

  stk::mesh::EntityVector delete_element2_on_p0_and_get_shared_nodes()
  {
    stk::mesh::EntityVector sharedNodes;
    if(get_proc_rank() == 0)
    {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2u);
      sharedNodes = get_shared_nodes_of_element(elem2);

      bool isDestroyed = get_bulk().destroy_entity(elem2);
      EXPECT_TRUE(isDestroyed);
    }

    return sharedNodes;
  }

  int get_proc_rank()
  {
    return get_bulk().parallel_rank();
  }

  void create_orphaned_side_on_p0()
  {
    if (get_proc_rank() == 0)
    {
      stk::topology quad4Topology = stk::topology::QUAD_4;
      stk::mesh::Part &quad4Part = get_meta().get_topology_root_part(quad4Topology);
      std::vector<stk::mesh::Entity> nodes(quad4Topology.num_nodes());
      for(unsigned i=0; i<quad4Topology.num_nodes(); ++i) {
        stk::mesh::EntityId nodeId = i+1;
        nodes[i] = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      }
      stk::mesh::Entity face1 = get_bulk().declare_solo_side(1, {&quad4Part});
      stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
      for(unsigned i=0; i<quad4Topology.num_nodes(); ++i) {
        get_bulk().declare_relation(face1, nodes[i], i, perm);
      }
    }
  }

  void create_face_with_shared_nodes(const stk::mesh::EntityVector &sharedNodes)
  {
    get_bulk().modification_begin();

    stk::mesh::Part &quad4Part = get_meta().get_topology_root_part(stk::topology::QUAD_4);

    if (get_proc_rank() == 0) {
      stk::mesh::Entity face = get_bulk().declare_solo_side({&quad4Part});
      stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation> (0);
      for(unsigned i = 0; i<sharedNodes.size(); ++i)
        get_bulk().declare_relation(face, sharedNodes[i], i, perm);
    }
    else if (get_proc_rank() == 1) {
      stk::mesh::Entity elem5 = get_bulk().get_entity(stk::topology::ELEM_RANK, 5u);
      get_bulk().declare_element_side(elem5, 4, stk::mesh::ConstPartVector{&quad4Part});
    }

    get_bulk().modification_end();
  }

  void run(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:3x1x2", auraOption);

    get_bulk().modification_begin();
    stk::mesh::EntityVector sharedNodes = delete_element2_on_p0_and_get_shared_nodes();
    get_bulk().modification_end();

    if(get_proc_rank() == 1)
      sharedNodes = get_shared_nodes_of_element5_on_p1();

    create_face_with_shared_nodes(sharedNodes);
    std::vector<stk::mesh::Entity> orphanedSides = stk::mesh::get_orphaned_sides_with_attached_element_on_different_proc(get_bulk());

    bool faceLeftBehind = !orphanedSides.empty();
    bool faceLeftBehindOnAnyProc = stk::is_true_on_any_proc(get_bulk().parallel(), faceLeftBehind);
    ASSERT_TRUE(faceLeftBehindOnAnyProc);
  }

  void runTrueOrphanedSide(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x2", auraOption);

    get_bulk().modification_begin();
    delete_element1_on_p0();
    get_bulk().modification_end();

    get_bulk().modification_begin();
    create_orphaned_side_on_p0();
    get_bulk().modification_end();

    std::vector<stk::mesh::Entity> orphanedSides = stk::mesh::get_orphaned_sides_with_attached_element_on_different_proc(get_bulk());

    bool faceLeftBehind = !orphanedSides.empty();
    bool faceLeftBehindOnAnyProc = stk::is_true_on_any_proc(get_bulk().parallel(), faceLeftBehind);
    ASSERT_FALSE(faceLeftBehindOnAnyProc);
  }
};

///////////////////////////////////////////////////////////////////

TEST_F(UnitTestRuleThreeViolation, keyholeWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    run(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(UnitTestRuleThreeViolation, keyholeWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    run(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(UnitTestRuleThreeViolation, trueOrphanedSideWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    runTrueOrphanedSide(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(UnitTestRuleThreeViolation, trueOrphanedSideWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    runTrueOrphanedSide(stk::mesh::BulkData::NO_AUTO_AURA);
}

}
