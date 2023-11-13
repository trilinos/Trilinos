#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "../FaceCreatorFixture.hpp"
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace
{

class UnitTestFaceSharingUsingGraph : public FaceCreatorFixture
{
protected:
  virtual void test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face()
  {
    only_proc_0_makes_a_face();
    test_that_num_sides_is_expected_value(1);
    unsigned num_faces_this_proc = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(get_meta().side_rank()));
    test_that_num_sides_is_correct(num_faces_this_proc);
  }

  virtual void create_faces_only_one_proc(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank()==0)
    {
      create_face_per_proc(element, nodes_of_face);
    }
    test_that_num_sides_is_expected_value(1);

    get_bulk().modification_end();
  }

  int get_proc_rank()
  {
    return get_bulk().parallel_rank();
  }

  void test_that_num_sides_is_correct(unsigned num_faces_this_proc)
  {
    std::vector<unsigned> gold_num_faces_per_proc={1,0};
    EXPECT_EQ(gold_num_faces_per_proc[get_proc_rank()], num_faces_this_proc);
  }

  void run(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x2", auraOption);
    test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();

    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(get_meta().globally_shared_part(), get_bulk().buckets(get_meta().side_rank()), sides);
    EXPECT_EQ(4u, stk::mesh::count_selected_entities(get_meta().globally_shared_part(), get_bulk().buckets(stk::topology::NODE_RANK)));
    ASSERT_EQ(1u, sides.size());
    ASSERT_EQ(4u, get_bulk().num_nodes(sides[0]));
    const stk::mesh::Entity *nodes = get_bulk().begin_nodes(sides[0]);
    EXPECT_EQ(5u, get_bulk().identifier(nodes[0]));
    EXPECT_EQ(6u, get_bulk().identifier(nodes[1]));
    EXPECT_EQ(8u, get_bulk().identifier(nodes[2]));
    EXPECT_EQ(7u, get_bulk().identifier(nodes[3]));
  }
};

///////////////////////////////////////////////////////////////////

TEST_F(UnitTestFaceSharingUsingGraph, twoHexesTwoProcsCreateOneFaceWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    run(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(UnitTestFaceSharingUsingGraph, twoHexesTwoProcsCreateOneFaceWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
    run(stk::mesh::BulkData::NO_AUTO_AURA);
}

}
