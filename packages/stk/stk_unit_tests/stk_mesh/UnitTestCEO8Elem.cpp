#include "UnitTestCEO8Elem.hpp"
#include "UnitTestCEOCommonUtils.hpp"
#include "Setup8Quad4ProcMesh.hpp"
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE

namespace CEOUtils
{

void check_state_proc_0(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part,
                        Part* owned_part,
                        stk::mesh::Part* topo_part,
                        stk::mesh::Part* block_1,
                        Part* aura_part,
                        Part* shared_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 1, 2, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 2, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
}

void check_state_proc_1(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part,
                        Part* aura_part,
                        stk::mesh::Part* topo_part,
                        stk::mesh::Part* block_1,
                        Part* owned_part,
                        Part* shared_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 1, 2, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
}

void check_state_proc_2(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part,
                        Part* aura_part,
                        stk::mesh::Part* topo_part,
                        stk::mesh::Part* block_1,
                        Part* owned_part,
                        Part* shared_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 2, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3, 4, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void check_state_proc_3(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part,
                        Part* aura_part,
                        stk::mesh::Part* topo_part,
                        stk::mesh::Part* block_1,
                        Part* owned_part,
                        Part* shared_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3, 4, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void fillMeshfor8Elem4ProcMoveTopAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
  //
  //     id/proc                           id/proc
  //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
  //
  // This test moves ownership of elements 6 and 7 (as well as their locally-owned
  // nodes) to procs 3 and 0, respectively.
  //

  const int p_rank = mesh.parallel_rank();

  //we don't need this because the setup8Quad4ProcMesh2D declares this part already as block_1
  //stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);


  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  setup8Quad4ProcMesh2D(mesh);

  stk::mesh::Part * block_1 = meta.get_part("block_1");

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  // Check the initial state
  if(p_rank == 0)
  {
    check_state_proc_0(mesh, universal_part, owned_part, topo_part, block_1, aura_part, shared_part);
  }
  else if(p_rank == 1)
  {
    check_state_proc_1(mesh, universal_part, aura_part, topo_part, block_1, owned_part, shared_part);
  }
  else if(p_rank == 2)
  {
    check_state_proc_2(mesh, universal_part, aura_part, topo_part, block_1, owned_part, shared_part);
  }
  else if(p_rank == 3)
  {
    check_state_proc_3(mesh, universal_part, aura_part, topo_part, block_1, owned_part, shared_part);
  }
}




//////////////////////////////////// 8Elem4ProcMoveTop //////////////////////////////////////

void check_state_proc_0(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part, Part* owned_part, Part* shared_part,
                        stk::mesh::Part* topo_part, stk::mesh::Part* block_1)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 2, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1 /*, 2 */));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 1 /* , 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 5 /* 1, 2, 5, 6*/ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED /* STATE_SHARED, 1, 2, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, /* shared_part, */  topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, /* 2, 3, 6, */  7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED /* STATE_SHARED, 2, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_TO, 1 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 7 /* 3, 4, 7, 8 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 1 /* 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5 /* , 6 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED /* STATE_SHARED, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, /* 6, */  7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED /* STATE_SHARED, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_TO, 1, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7 /* , 8 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void check_state_proc_1(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part, Part* owned_part, Part* shared_part,
                        stk::mesh::Part* topo_part, stk::mesh::Part* block_1)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, /* 1, */ 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 3 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2 /* , 3 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0 /* , 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2 /* 1, 2, 5, 6 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 2 /* 0, 2, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2 /* , 3, 6, 7 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 2 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
}

void check_state_proc_2(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part, Part* owned_part, Part* shared_part,
                        stk::mesh::Part* topo_part, stk::mesh::Part* block_1)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, /* 2, */  3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3 /*, 4 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 1 /* 0, 1, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 3 /* 2, 3, 6, 7 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, /* 0, */ 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3 /* , 4, 7, 8 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void check_state_proc_3(stk::unit_test_util::BulkDataTester& mesh,
                        Part* universal_part, Part* owned_part, Part* shared_part,
                        stk::mesh::Part* topo_part, stk::mesh::Part* block_1)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /** STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM  /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM /** STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, /* 3, */ 4 ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED /* STATE_SHARED, 0, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, /* 1, 2, 5, */ 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED /* STATE_SHARED, 0, 1, 2 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, /* 2, 3, */  6 /*, 7 */));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, /* 0, */  2   ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, /*3,*/ 4, /*7,*/ 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED /* STATE_SHARED, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, /* 5, */ 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED /* STATE_SHARED, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6 /*, 7*/ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 2 /* 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, /* 7, */ 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}


void checkStatesAfterCEO_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh)
{
  //
  //     id/proc                           id/proc
  //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
  //
  // This test moves ownership of elements 6 and 7 (as well as their locally-owned
  // nodes) to procs 3 and 0, respectively.
  //

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  // Part * aura_part      = &meta.aura_part();

  int p_rank = mesh.parallel_rank();

  // Check the final state
  if (p_rank == 0) {
    check_state_proc_0(mesh, universal_part, owned_part, shared_part, topo_part, block_1);
  }
  else if (p_rank == 1) {
    check_state_proc_1(mesh, universal_part, owned_part, shared_part, topo_part, block_1);
  }
  else if (p_rank == 2) {
    check_state_proc_2(mesh, universal_part, owned_part, shared_part, topo_part, block_1);
  }
  else if (p_rank == 3) {
    check_state_proc_3(mesh, universal_part, owned_part, shared_part, topo_part, block_1);
  }
}


void p0_elem_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_TO, 1, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));
}

void p0_node_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_TO, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 1, 3    ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 1, 2, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 2, 3    ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_TO, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_TO, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void p1_elem_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
}

void p1_node_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0, 3   ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 2, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
}

void p2_elem_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));
}

void p2_node_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_TO, 0, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 1, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 0, 3   ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void p3_elem_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_TO, 0, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));
}

void p3_node_checks_for_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester& mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * block_1 = meta.get_part("block_1");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0, 1   ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 1, 2));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 0, 2   ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_TO, 1, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
}

void checkStatesAfterCEOME_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh)
{
  //
  //     id/proc                           id/proc
  //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
  //
  // This test moves ownership of elements 6 and 7 (as well as their locally-owned
  // nodes) to procs 3 and 0, respectively.
  //

  int p_rank = mesh.parallel_rank();

  // Check the final state
  if (p_rank == 0) {
    p0_elem_checks_for_8Elem4ProcMoveTop(mesh);
    p0_node_checks_for_8Elem4ProcMoveTop(mesh);
  }
  else if (p_rank == 1) {
    p1_elem_checks_for_8Elem4ProcMoveTop(mesh);
    p1_node_checks_for_8Elem4ProcMoveTop(mesh);
  }
  else if (p_rank == 2) {
    p2_elem_checks_for_8Elem4ProcMoveTop(mesh);
    p2_node_checks_for_8Elem4ProcMoveTop(mesh);
  }
  else if (p_rank == 3) {
    p3_elem_checks_for_8Elem4ProcMoveTop(mesh);
    p3_node_checks_for_8Elem4ProcMoveTop(mesh);
  }
}


}
