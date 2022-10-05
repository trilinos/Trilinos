#include "UnitTestCEO4ElemEdge.hpp"
#include "UnitTestCEOCommonUtils.hpp"
#include "stk_mesh/base/Field.hpp"      // for Field
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE

namespace CEOUtils
{

void fill_mesh_4elem_4proc_check_p0(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                    Part* universal_part, Part* owned_part, Part* elem_part, Part* elem_topo_part,
                                    Part* aura_part, Part* node_part, Part* node_topo_part,
                                    Part* shared_part, Part* elem_block)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED,  0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED,  1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 2 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 2 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void fill_mesh_4elem_4proc_check_p1(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh, Part* universal_part,
                                    Part* aura_part, Part* elem_part, Part* elem_topo_part,
                                    Part* owned_part, Part* elem_block, Part* edge_part,
                                    Part* edge_topo_part, Part* node_part, Part* node_topo_part,
                                    Part* shared_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 2  ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_GHOSTED_FROM, 2));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void fill_mesh_4elem_4proc_check_p2(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                    Part* universal_part, Part* aura_part, Part* elem_part, Part* elem_topo_part,
                                    Part* owned_part, Part* elem_block, Part* shared_part, Part* edge_part,
                                    Part* edge_topo_part, Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 2   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 3  ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, owned_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

void fill_mesh_4elem_4proc_check_p3(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                    Part* universal_part, Part* aura_part, Part* elem_part, Part* elem_topo_part,
                                    Part* elem_block, Part* owned_part, Part* shared_part, Part* edge_part,
                                    Part* edge_topo_part, Part* node_part, Part* node_topo_part)
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
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 2   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 2  ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 2));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

//////////////////////////////////// 4Elem4ProcEdge //////////////////////////////////////

void fillMeshfor4Elem4ProcEdgeAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data,
                                      EntityKey &elem_key_chg_own)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    ||     |          |      |     |    ||     |
  //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
  //        |      |     |    ||     |          |      |     |    ||     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
  //  this edge moves to p0 --^
  //  element 3 moves to proc 0.
  //  nodes 7&8 move to proc 0.
  //  proc 2 forgets everything.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  int p_rank = mesh.parallel_rank();

  Part* elem_part = &meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  Part* edge_part = &meta_data.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  Part* node_part = &meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  Part* elem_block = &meta_data.declare_part("elem_block", stk::topology::ELEMENT_RANK);
  typedef stk::mesh::Field<double> rank_zero_field;
  rank_zero_field  & f0 = meta_data.declare_field<double>( NODE_RANK, "elem_field" );
  stk::mesh::put_field_on_mesh( f0 , *elem_block , nullptr);
  meta_data.commit();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  elem_key_chg_own= EntityKey(ELEM_RANK, 3 /*id*/);
  EntityKey edge_key_chg_own(EDGE_RANK, 33);  // Third side of elem 3
  EntityKey node_A_key_chg_own(NODE_RANK, 7 /*id*/);
  EntityKey node_B_key_chg_own(NODE_RANK, 8 /*id*/);
  EntityKey node_C_key(NODE_RANK, 5 /*id*/);
  EntityKey node_D_key(NODE_RANK, 6 /*id*/);

  // Create element
  stk::mesh::EntityId elemId = p_rank + 1;
  stk::mesh::Entity elem = mesh.declare_element(elemId, stk::mesh::ConstPartVector{elem_part});

  // If it is 2nd to last element, it is the one changing
  if(p_rank == 2)
  {
    EXPECT_TRUE(elem_key_chg_own == mesh.entity_key(elem));
  }

  // Create nodes
  stk::mesh::EntityVector nodes;
  if(p_rank == 0)
  {
    nodes.push_back(mesh.declare_node(1, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(2, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(3, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(4, stk::mesh::ConstPartVector{node_part}));
  }
  else if(p_rank == 1)
  {
    nodes.push_back(mesh.declare_node(3, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(4, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(5, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(6, stk::mesh::ConstPartVector{node_part}));
  }
  else if(p_rank == 2)
  {
    PartVector add_parts;
    PartVector remove_parts;
    add_parts.push_back(elem_block);
    mesh.change_entity_parts(elem,add_parts,remove_parts);
    nodes.push_back(mesh.declare_node(5, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(6, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(7, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(8, stk::mesh::ConstPartVector{node_part}));
  }
  else
  { // p_rank == 3
    nodes.push_back(mesh.declare_node(7, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(8, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(9, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(10,stk::mesh::ConstPartVector{node_part}));
  }

  // Add element relations to nodes
  unsigned rel_id = 0;
  for(stk::mesh::EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
  {
    mesh.declare_relation(elem, *itr, rel_id);
  }

  if(p_rank == 0)
  {
    mesh.add_node_sharing(nodes[2], 1);
    mesh.add_node_sharing(nodes[3], 1);
  }
  else if((p_rank == 1) || (p_rank == 2))
  {
    mesh.add_node_sharing(nodes[0], p_rank - 1);
    mesh.add_node_sharing(nodes[1], p_rank - 1);
    mesh.add_node_sharing(nodes[2], p_rank + 1);
    mesh.add_node_sharing(nodes[3], p_rank + 1);
  }
  else
  { // p_rank ==3
    mesh.add_node_sharing(nodes[0], p_rank - 1);
    mesh.add_node_sharing(nodes[1], p_rank - 1);
  }

  mesh.modification_end();

  // Create edge on p2 (it will be shared to p3 and aura'd to p1)

  mesh.initialize_face_adjacent_element_graph();

  mesh.modification_begin();
  if (p_rank == 2) {
    mesh.declare_element_side(elem, 2, stk::mesh::ConstPartVector{edge_part});
  }
  mesh.modification_end();

  Part* elem_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  Part* edge_topo_part = &meta_data.get_topology_root_part(stk::topology::LINE_2);
  Part* node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part* universal_part = &meta_data.universal_part();
  Part* owned_part     = &meta_data.locally_owned_part();
  Part* shared_part    = &meta_data.globally_shared_part();
  Part* aura_part      = &meta_data.aura_part();

  //test pre-conditions
  if(p_rank == 0)
  {
    fill_mesh_4elem_4proc_check_p0(edge_key_chg_own, mesh, universal_part,
                                   owned_part, elem_part, elem_topo_part, aura_part, node_part,
                                   node_topo_part, shared_part, elem_block);
  }
  else if(p_rank == 1)
  {
    fill_mesh_4elem_4proc_check_p1(edge_key_chg_own, mesh, universal_part,
                                   aura_part, elem_part, elem_topo_part, owned_part, elem_block,
                                   edge_part, edge_topo_part, node_part, node_topo_part, shared_part);
  }
  else if(p_rank == 2)
  {
    fill_mesh_4elem_4proc_check_p2(edge_key_chg_own, mesh, universal_part,
                                   aura_part, elem_part, elem_topo_part, owned_part, elem_block,
                                   shared_part, edge_part, edge_topo_part, node_part, node_topo_part);
  }
  else if(p_rank == 3)
  {
    fill_mesh_4elem_4proc_check_p3(edge_key_chg_own, mesh, universal_part,
                                   aura_part, elem_part, elem_topo_part, elem_block, owned_part,
                                   shared_part, edge_part, edge_topo_part, node_part, node_topo_part);
  }
}


void check_state_after_CEO_4elem_4proc_edge_p0(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                               Part* universal_part, Part* owned_part, Part* shared_part, Part* aura_part,
                                               Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                               Part* edge_part, Part* edge_topo_part,
                                               Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_block, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, owned_part, /* shared_part, */ elem_block, elem_part, elem_topo_part, edge_part, edge_topo_part));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3 /*, 4 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED            /* STATE_SHARED, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1       /* STATE_NOT_GHOSTED_FROM */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, /* shared_part, */ aura_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED           /* STATE_SHARED, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1      /* STATE_NOT_GHOSTED_FROM */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, /* shared_part, */ aura_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, /* shared_part,*/ elem_block, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3 /*, 4 */ ));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, /* shared_part, */ elem_block, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3 /*, 4 */));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM    /* STATE_GHOSTED_FROM, 3 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

void check_state_after_CEO_4elem_4proc_edge_p1(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                               Part* universal_part, Part* owned_part, Part* shared_part, Part* aura_part,
                                               Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                               Part* /*edge_part*/, Part* /*edge_topo_part*/,
                                               Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));
  stk::mesh::Entity entity = mesh.get_entity(EntityKey(ELEM_RANK,1));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Unchanged);

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));
  entity = mesh.get_entity(EntityKey(ELEM_RANK,2));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Modified);

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  //  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));
  entity = mesh.get_entity(EntityKey(NODE_RANK,1));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Unchanged);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));
  entity = mesh.get_entity(EntityKey(NODE_RANK,2));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Unchanged);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));
  entity = mesh.get_entity(EntityKey(NODE_RANK,3));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Unchanged);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));
  entity = mesh.get_entity(EntityKey(NODE_RANK,4));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Unchanged);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 2 /* 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0 /* 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2 /* , 3 */ ));
  entity = mesh.get_entity(EntityKey(NODE_RANK,5));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Modified);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2 /* 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0/*, 3 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2 /* , 3 */));
  entity = mesh.get_entity(EntityKey(NODE_RANK,6));
  EXPECT_TRUE(mesh.state(entity) == stk::mesh::Modified);

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM /*STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
  //                          edge_part, edge_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
  //              edge_part, edge_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void check_state_after_CEO_4elem_4proc_edge_p2(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                               Part* /*universal_part*/, Part* /*owned_part*/, Part* /*shared_part*/, Part* /*aura_part*/,
                                               Part* /*elem_part*/, Part* /*elem_topo_part*/, Part* /*elem_block*/,
                                               Part* /*edge_part*/, Part* /*edge_topo_part*/,
                                               Part* /*node_part*/, Part* /*node_topo_part*/)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 3 /*STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void check_state_after_CEO_4elem_4proc_edge_p3(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                               Part* universal_part, Part* owned_part, Part* shared_part, Part* /*aura_part*/,
                                               Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                               Part* edge_part, Part* edge_topo_part,
                                               Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 0 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 2 /*, 0 */ ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, shared_part, elem_block, elem_part, elem_topo_part, edge_part, edge_topo_part));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, /* 3, */ 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  //  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  //  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  //  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 2 /* 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_block, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, /* 3, */ 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 /* 0 */ ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_block, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, /* 3, */ 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0 */ ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0 */));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

void checkStatesAfterCEO_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    ||     |          |      |     |    ||     |
  //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
  //        |      |     |    ||     |          |      |     |    ||     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
  //  this edge moves to p0 --^
  //  element 3 moves to proc 0.
  //  nodes 7&8 move to proc 0.
  //  proc 2 forgets everything.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  int p_rank = mesh.parallel_rank();

  const stk::mesh::MetaData &meta_data = mesh.mesh_meta_data();

  Part* elem_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  Part* edge_topo_part = &meta_data.get_topology_root_part(stk::topology::LINE_2);
  Part* node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part* elem_part =      meta_data.get_part("elem_part");
  Part* elem_block =     meta_data.get_part("elem_block");
  Part* edge_part =      meta_data.get_part("edge_part");
  Part* node_part =      meta_data.get_part("node_part");
  Part* universal_part = &meta_data.universal_part();
  Part* owned_part     = &meta_data.locally_owned_part();
  Part* shared_part    = &meta_data.globally_shared_part();
  Part* aura_part      = &meta_data.aura_part();
  EntityKey edge_key_chg_own(EDGE_RANK, 33);  // Third side of elem 3

  //test post condition
  if (p_rank == 0) {
    check_state_after_CEO_4elem_4proc_edge_p0(edge_key_chg_own, mesh,
                                              universal_part, owned_part, shared_part, aura_part,
                                              elem_part, elem_topo_part, elem_block,
                                              edge_part, edge_topo_part,
                                              node_part, node_topo_part);
  }
  else if (p_rank == 1) {
    check_state_after_CEO_4elem_4proc_edge_p1(edge_key_chg_own, mesh,
                                              universal_part, owned_part, shared_part, aura_part,
                                              elem_part, elem_topo_part, elem_block,
                                              edge_part, edge_topo_part,
                                              node_part, node_topo_part);
  }
  else if (p_rank == 2) {
    check_state_after_CEO_4elem_4proc_edge_p2(edge_key_chg_own, mesh,
                                              universal_part, owned_part, shared_part, aura_part,
                                              elem_part, elem_topo_part, elem_block,
                                              edge_part, edge_topo_part,
                                              node_part, node_topo_part);
  }
  else if (p_rank == 3) {
    check_state_after_CEO_4elem_4proc_edge_p3(edge_key_chg_own, mesh,
                                              universal_part, owned_part, shared_part, aura_part,
                                              elem_part, elem_topo_part, elem_block,
                                              edge_part, edge_topo_part,
                                              node_part, node_topo_part);
  }
}


void check_state_after_CEOME_4elem_4proc_edge_p0(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                                 Part* universal_part, Part* owned_part, Part* shared_part, Part* aura_part,
                                                 Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                                 Part* edge_part, Part* edge_topo_part,
                                                 Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 3  ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_GHOSTED_TO, 1));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, owned_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 1 ));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

void check_state_after_CEOME_4elem_4proc_edge_p1(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                                 Part* universal_part, Part* owned_part, Part* shared_part, Part* aura_part,
                                                 Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                                 Part* edge_part, Part* edge_topo_part,
                                                 Part* node_part, Part* node_topo_part)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 3));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void check_state_after_CEOME_4elem_4proc_edge_p2(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                                 Part* /*universal_part*/, Part* /*owned_part*/, Part* /*shared_part*/, Part* /*aura_part*/,
                                                 Part* /*elem_part*/, Part* /*elem_topo_part*/, Part* /*elem_block*/,
                                                 Part* /*edge_part*/, Part* /*edge_topo_part*/,
                                                 Part* /*node_part*/, Part* /*node_topo_part*/)
{
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
}

void check_state_after_CEOME_4elem_4proc_edge_p3(EntityKey edge_key_chg_own, stk::unit_test_util::BulkDataTester& mesh,
                                                 Part* universal_part, Part* owned_part, Part* shared_part, Part* aura_part,
                                                 Part* elem_part, Part* elem_topo_part, Part* elem_block,
                                                 Part* edge_part, Part* edge_topo_part,
                                                 Part* node_part, Part* node_topo_part)
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
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0));
  EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
  EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, edge_key_chg_own.id()));

  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_VALID));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_OWNED, 0   ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_SHARED, 0  ));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, edge_key_chg_own       , STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, edge_key_chg_own       , universal_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , NODE_RANK, 7, 8));
  EXPECT_TRUE(check_relns(mesh, edge_key_chg_own       , ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 0 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_part, elem_topo_part,
                          edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, edge_key_chg_own.id()));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 0));
  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
}

void checkStatesAfterCEOME_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    ||     |          |      |     |    ||     |
  //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
  //        |      |     |    ||     |          |      |     |    ||     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
  //  this edge moves to p0 --^
  //  element 3 moves to proc 0.
  //  nodes 7&8 move to proc 0.
  //  proc 2 forgets everything.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  int p_rank = mesh.parallel_rank();

  const stk::mesh::MetaData &meta_data = mesh.mesh_meta_data();

  Part* elem_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  Part* edge_topo_part = &meta_data.get_topology_root_part(stk::topology::LINE_2);
  Part* node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part* elem_part =      meta_data.get_part("elem_part");
  Part* edge_part =      meta_data.get_part("edge_part");
  Part* node_part =      meta_data.get_part("node_part");
  Part* elem_block =     meta_data.get_part("elem_block");
  Part* universal_part = &meta_data.universal_part();
  Part* owned_part     = &meta_data.locally_owned_part();
  Part* shared_part    = &meta_data.globally_shared_part();
  Part* aura_part      = &meta_data.aura_part();
  EntityKey edge_key_chg_own(EDGE_RANK, 33);  // Third side of elem 3

  //test post condition
  if (p_rank == 0) {
    check_state_after_CEOME_4elem_4proc_edge_p0(edge_key_chg_own, mesh,
                                                universal_part, owned_part, shared_part, aura_part,
                                                elem_part, elem_topo_part, elem_block,
                                                edge_part, edge_topo_part,
                                                node_part, node_topo_part);
  }
  else if (p_rank == 1) {
    check_state_after_CEOME_4elem_4proc_edge_p1(edge_key_chg_own, mesh,
                                                universal_part, owned_part, shared_part, aura_part,
                                                elem_part, elem_topo_part, elem_block,
                                                edge_part, edge_topo_part,
                                                node_part, node_topo_part);
  }
  else if (p_rank == 2) {
    check_state_after_CEOME_4elem_4proc_edge_p2(edge_key_chg_own, mesh,
                                                universal_part, owned_part, shared_part, aura_part,
                                                elem_part, elem_topo_part, elem_block,
                                                edge_part, edge_topo_part,
                                                node_part, node_topo_part);
  }
  else if (p_rank == 3) {
    check_state_after_CEOME_4elem_4proc_edge_p3(edge_key_chg_own, mesh,
                                                universal_part, owned_part, shared_part, aura_part,
                                                elem_part, elem_topo_part, elem_block,
                                                edge_part, edge_topo_part,
                                                node_part, node_topo_part);
  }
}

}
