#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_util/parallel/Parallel.hpp>

/* Ticket summary:
 * 13227: can we detect a doubly sided sideset? App might want api is_doubly_sided_sideset(sideset id);
 *      Need definition of doubly_sideset sideset: if any (element1_side1, element2_sid2) is an edge in the graph
 *              and from the same sideset, then it is a doubly sided sideset.
 */


class Ticket13227 : public stk::unit_test_util::MeshFixture
{
};

TEST_F(Ticket13227, DISABLED_doubly_sided_sideset)
{
  if(stk::parallel_machine_size(get_comm())==1)
  {
    setup_mesh("meshed_enclosure_double_sideset.g", stk::mesh::BulkData::NO_AUTO_AURA);
    // EXPECT_TRUE(sideset_2 is doubly_sided);
    // EXPECT_TRUE(sideset_1 is singly_sided);
  }
}

TEST_F(Ticket13227, DISABLED_singly_sided_sideset)
{
  if(stk::parallel_machine_size(get_comm())==1)
  {
    setup_mesh("meshed_enclosure_single_sideset.g", stk::mesh::BulkData::NO_AUTO_AURA);
    // EXPECT_TRUE(sideset_2 is singly_sided);
    // EXPECT_TRUE(sideset_1 is singly_sided);
  }
}
