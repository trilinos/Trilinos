#include "UnitTestCEO4ElemRotate.hpp"
#include "UnitTestCEOCommonUtils.hpp"
#include "stk_mesh/base/Field.hpp"      // for Field
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE

namespace CEOUtils
{

void fill_mesh_4elem_4proc_check_p0(stk::unit_test_util::BulkDataTester& mesh,
                                    stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                    stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                    stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                    stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, shared_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void fill_mesh_4elem_4proc_check_p1(stk::unit_test_util::BulkDataTester& mesh,
                                    stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                    stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                    stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                    stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 0, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, owned_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void fill_mesh_4elem_4proc_check_p2(stk::unit_test_util::BulkDataTester& mesh,
                                    stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                    stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                    stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                    stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, owned_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void fill_mesh_4elem_4proc_check_p3(stk::unit_test_util::BulkDataTester& mesh,
                                    stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                    stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                    stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                    stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

//////////////////////////////////// 4Elem4ProcRotate //////////////////////////////////////

void fillMeshfor4Elem4ProcRotateAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
  //
  //     id/proc                id/proc
  //      7/3---8/2---9/2        7/2---8/1---9/1
  //       |     |     |          |     |     |
  //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
  //       |     |     |          |     |     |
  //      4/0---5/0---6/1  -->   4/3---5/3---6/0
  //       |     |     |          |     |     |
  //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
  //       |     |     |          |     |     |
  //      1/0---2/0---3/1        1/3---2/3---3/0

  const int p_rank = mesh.parallel_rank();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_4 = meta.declare_part_with_topology("block_4", stk::topology::QUAD_4_2D);
  meta.commit();

  // 1 elem-id for each proc
  stk::mesh::EntityId proc_elemIDs[] = {1, 2, 3, 4};

  // list of node-ids for each element
  stk::mesh::EntityIdVector elem_nodeIDs[] {
    {1, 2, 5, 4},
    {2, 3, 6, 5},
    {5, 6, 9, 8},
    {4, 5, 8, 7}
  };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  int shared_nodeIDs_and_procs[][3] = {
    {0, 2, 1}, {0, 5, 1}, {0, 5, 2}, {0, 5, 3}, {0, 4, 3}, // proc 0
    {1, 2, 0}, {1, 6, 2}, {1, 5, 0}, {1, 5, 2}, {1, 5, 3}, // proc 1
    {2, 5, 0}, {2, 5, 1}, {2, 5, 3}, {2, 6, 1}, {2, 8, 3}, // proc 2
    {3, 4, 0}, {3, 5, 0}, {3, 5, 1}, {3, 5, 2}, {3, 8, 2}  // proc 3
  };
  int numSharedNodeTriples = 20;

  mesh.modification_begin();

  stk::mesh::EntityId elemID = proc_elemIDs[p_rank];
  if (0 == p_rank) {
    stk::mesh::declare_element(mesh, block_1, elemID, elem_nodeIDs[p_rank]);
  }
  else if (1 == p_rank) {
    stk::mesh::declare_element(mesh, block_2, elemID, elem_nodeIDs[p_rank]);
  }
  else if (2 == p_rank) {
    stk::mesh::declare_element(mesh, block_3, elemID, elem_nodeIDs[p_rank]);
  }
  else if (3 == p_rank) {
    stk::mesh::declare_element(mesh, block_4, elemID, elem_nodeIDs[p_rank]);
  }

  for (int proc = 0; proc < numSharedNodeTriples; ++proc) {
    if (p_rank == shared_nodeIDs_and_procs[proc][0]) {
      int nodeID = shared_nodeIDs_and_procs[proc][1];
      int sharingProc = shared_nodeIDs_and_procs[proc][2];
      stk::mesh::Entity node = mesh.get_entity(NODE_RANK, nodeID);
      mesh.add_node_sharing(node, sharingProc);
    }
  }

  mesh.modification_end();

  stk::mesh::Part * block1_part = meta.get_part("block_1");
  stk::mesh::Part * block2_part = meta.get_part("block_2");
  stk::mesh::Part * block3_part = meta.get_part("block_3");
  stk::mesh::Part * block4_part = meta.get_part("block_4");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  // Check the initial state

  //these are the same everywhere, everything is valid, always
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));

  if (p_rank == 0) {
    fill_mesh_4elem_4proc_check_p0(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                   block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 1) {
    fill_mesh_4elem_4proc_check_p1(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                   block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 2) {
    fill_mesh_4elem_4proc_check_p2(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                   block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 3) {
    fill_mesh_4elem_4proc_check_p3(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                   block1_part, block2_part, block3_part, block4_part);
  }
}

void check_state_after_CEO_4elem_4proc_rotate_p0(stk::unit_test_util::BulkDataTester& mesh,
                                                 stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                 stk::mesh::Part* shared_part, stk::mesh::Part* topo_part,
                                                 stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                 stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1));  //shared with 3 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 2));  //connected to elem 1 & 2 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));  //ghosted to 1,2,3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));  //related to 1,2,3,4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED)); //shared with 1 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO)); //ghosted to 2&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, block2_part, block3_part)); //shared_part, universal_part, owned_part, topo_part, block2_part, block3_part)); at end
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));  //2 and 3 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));

}

void check_state_after_CEO_4elem_4proc_rotate_p1(stk::unit_test_util::BulkDataTester& mesh,
                                                 stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                 stk::mesh::Part* shared_part, stk::mesh::Part* topo_part,
                                                 stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                 stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO)); //0,2&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));  //1,2,3&4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2)); // 0 at end, not 2
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3)); //2&3 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED)); //2 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO)); //0&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, topo_part, block3_part, block4_part));//add shared_part at end
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3)); //3&4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO)); //0,2&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void check_state_after_CEO_4elem_4proc_rotate_p2(stk::unit_test_util::BulkDataTester& mesh,
                                                 stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                 stk::mesh::Part* shared_part, stk::mesh::Part* topo_part,
                                                 stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                 stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO)); //0,1&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));  //3 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, topo_part, block1_part, block4_part));  //shared_part as well at end
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 4)); //1 &4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 4)); //1,2,3&4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO)); //0,1&3 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3));  //1 at end, not 3
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 4)); //3&4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void check_state_after_CEO_4elem_4proc_rotate_p3(stk::unit_test_util::BulkDataTester& mesh,
                                                 stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                 stk::mesh::Part* shared_part, stk::mesh::Part* topo_part,
                                                 stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                 stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO)); //0,1&2 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO)); //0,1&2 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED)); //0 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO)); //1&2 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, block1_part, block2_part));  //add shares_part by end
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));  //1 & 2 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0));  //shared to 2 at end
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));  //0&1 at end
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1)); // 1 & 4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1));  // 1,2,3 & 4 at end

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void checkStatesAfterCEO_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
  // check intermediate state after change_entity_owner() but before internal_modification_end()
  //     id/proc                id/proc
  //      7/3---8/2---9/2        7/2---8/1---9/1
  //       |     |     |          |     |     |
  //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
  //       |     |     |          |     |     |
  //      4/0---5/0---6/1  -->   4/3---5/3---6/0
  //       |     |     |          |     |     |
  //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
  //       |     |     |          |     |     |
  //      1/0---2/0---3/1        1/3---2/3---3/0
  //
  //NOTES: sharing for nodes that have changed owner has not been deleted, still in previous state
  //       ghosting deleted entirely
  //       shared_part updated to match the commmap/list version of sharing
  //       aura_part never present (no ghosting)

  const int p_rank = mesh.parallel_rank();

  stk::mesh::Part * block1_part = meta.get_part("block_1");
  stk::mesh::Part * block2_part = meta.get_part("block_2");
  stk::mesh::Part * block3_part = meta.get_part("block_3");
  stk::mesh::Part * block4_part = meta.get_part("block_4");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  //  Part * aura_part      = &meta.aura_part(); // gone

  if (p_rank == 0) {
    check_state_after_CEO_4elem_4proc_rotate_p0(mesh, universal_part, owned_part, shared_part, topo_part,
                                                block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 1) {
    check_state_after_CEO_4elem_4proc_rotate_p1(mesh, universal_part, owned_part, shared_part, topo_part,
                                                block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 2) {
    check_state_after_CEO_4elem_4proc_rotate_p2(mesh, universal_part, owned_part, shared_part, topo_part,
                                                block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 3) {
    check_state_after_CEO_4elem_4proc_rotate_p3(mesh, universal_part, owned_part, shared_part, topo_part,
                                                block1_part, block2_part, block3_part, block4_part);
  }
}

void check_state_after_CEOME_4elem_4proc_rotate_p0(stk::unit_test_util::BulkDataTester& mesh,
                                                   stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                   stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                                   stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                   stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, owned_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void check_state_after_CEOME_4elem_4proc_rotate_p1(stk::unit_test_util::BulkDataTester& mesh,
                                                   stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                   stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                                   stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                   stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, owned_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void check_state_after_CEOME_4elem_4proc_rotate_p2(stk::unit_test_util::BulkDataTester& mesh,
                                                   stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                   stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                                   stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                   stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void check_state_after_CEOME_4elem_4proc_rotate_p3(stk::unit_test_util::BulkDataTester& mesh,
                                                   stk::mesh::Part* universal_part, stk::mesh::Part* owned_part,
                                                   stk::mesh::Part* shared_part, stk::mesh::Part* aura_part, stk::mesh::Part* topo_part,
                                                   stk::mesh::Part* block1_part, stk::mesh::Part* block2_part,
                                                   stk::mesh::Part* block3_part, stk::mesh::Part* block4_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 0, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 0, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 0, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
}

void checkStatesAfterCEOME_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta)
{

  //     id/proc                id/proc
  //      7/3---8/2---9/2        7/2---8/1---9/1
  //       |     |     |          |     |     |
  //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
  //       |     |     |          |     |     |
  //      4/0---5/0---6/1  -->   4/3---5/3---6/0
  //       |     |     |          |     |     |
  //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
  //       |     |     |          |     |     |
  //      1/0---2/0---3/1        1/3---2/3---3/0

  const int p_rank = mesh.parallel_rank();

  stk::mesh::Part * block1_part = meta.get_part("block_1");
  stk::mesh::Part * block2_part = meta.get_part("block_2");
  stk::mesh::Part * block3_part = meta.get_part("block_3");
  stk::mesh::Part * block4_part = meta.get_part("block_4");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  //these are the same everywhere, everything is valid, always
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));

  if (p_rank == 0) {
    check_state_after_CEOME_4elem_4proc_rotate_p0(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                                  block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 1) {
    check_state_after_CEOME_4elem_4proc_rotate_p1(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                                  block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 2) {
    check_state_after_CEOME_4elem_4proc_rotate_p2(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                                  block1_part, block2_part, block3_part, block4_part);
  }
  else if (p_rank == 3) {
    check_state_after_CEOME_4elem_4proc_rotate_p3(mesh, universal_part, owned_part, shared_part, aura_part, topo_part,
                                                  block1_part, block2_part, block3_part, block4_part);
  }
}

}
