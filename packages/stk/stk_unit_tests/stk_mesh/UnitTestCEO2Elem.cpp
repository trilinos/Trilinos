#include "UnitTestCEO2Elem.hpp"
#include "UnitTestCEOCommonUtils.hpp"
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE

namespace CEOUtils
{

//////////////////////////////////// 2Elem2ProcMove //////////////////////////////////////
void fillMeshfor2Elem2ProcMoveAndTest(stk::unit_test_util::BulkDataTester& bulk, stk::mesh::MetaData &meta, std::vector<stk::mesh::Entity>& elems)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  stk::mesh::EntityId element_ids[2] = {1, 2};
  stk::mesh::EntityIdVector elem_node_ids[] {
    {1, 2, 3, 4},
    {4, 3, 6, 5}
  };

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0

  bulk.modification_begin();
  int p_rank = bulk.parallel_rank();
  if(p_rank == 0)
  {
    elems.push_back(stk::mesh::declare_element(bulk, *elem_part, element_ids[0], elem_node_ids[0]));
    elems.push_back(stk::mesh::declare_element(bulk, *elem_part, element_ids[1], elem_node_ids[1]));
  }
  bulk.modification_end();

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  //    Part * shared_part    = &meta.globally_shared_part();
  //    Part * aura_part      = &meta.aura_part();

  // Check initial state
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  }
}

void checkStatesAfterCEO_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  // Part * shared_part    = &meta.globally_shared_part();
  // Part * aura_part      = &meta.aura_part();

  int p_rank = bulk.parallel_rank();
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  }
  else if(p_rank == 1)
  {
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

// valid vs not valid, owned.... shared vs not shared, ghosted_from vs not ghosted_from
// ghosted_to vs not_ghosted_to

void checkStatesAfterCEOME_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  int p_rank = bulk.parallel_rank();
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0 ));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0 ));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

//////////////////////////////////// 2Elem2ProcFlip //////////////////////////////////////

void fillMeshfor2Elem2ProcFlipAndTest(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::EntityId element_ids[2] = {1, 2};
  stk::mesh::EntityIdVector elem_node_ids[] {
    {1, 2, 3, 4},
    {4, 3, 6, 5}
  };

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  int p_rank = mesh.parallel_rank();

  mesh.modification_begin();
  if(p_rank == 0)
  {
    elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[0], elem_node_ids[0]));
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 1);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 1);
  }
  else if(p_rank == 1)
  {
    elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[1], elem_node_ids[1]));
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 0);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 0);
  }

  mesh.modification_end();

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

void checkStatesAfterCEOME_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  int p_rank = mesh.parallel_rank();

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

void checkStatesAfterCEO_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  //    Part * aura_part      = &meta.aura_part();

  int p_rank = mesh.parallel_rank();

  if(p_rank == 0)
  {
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    //        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 1 ));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    //        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 1 ));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    //        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 1 ));
    //        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    //        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    //        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));  // Was (1, 2)

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));  // Was (1, 2)

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));  // Was (1, 2)

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));  // Was (1, 2)

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  }
}

//these tests are for turning regenerate_aura off in various places

void checkStatesAfterCEOME_2Elem2ProcMove_no_ghost(stk::unit_test_util::BulkDataTester &bulk)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1

  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  //Part * aura_part      = &meta.aura_part(); //no aura, no aura part

  int p_rank = bulk.parallel_rank();
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1)); //used to be 1 & 2

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));  //used to be 1 & 2

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    //    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    //    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    //    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 2)); //used to be 1 & 2

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
    /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));//used to be 1 & 2

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

void fillMeshfor2Elem2ProcFlipAndTest_no_ghost(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::EntityId element_ids[2] = {1, 2};
  stk::mesh::EntityIdVector elem_node_ids[] {
    {1, 2, 3, 4},
    {4, 3, 6, 5}
  };

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  int p_rank = mesh.parallel_rank();

  mesh.modification_begin();
  if(p_rank == 0)
  {
    elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[0], elem_node_ids[0]));
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 1);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 1);
  }
  else if(p_rank == 1)
  {
    elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[1], elem_node_ids[1]));
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 0);
    mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 0);
  }

  mesh.modification_end();  //call IME through the tester to pass regenerate_aura = false

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  //Part * aura_part      = &meta.aura_part(); //no aura

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, topo_part, elem_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, elem_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

void checkStatesAfterCEOME_2Elem2ProcFlip_no_ghost(stk::unit_test_util::BulkDataTester& mesh)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);

  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  //Part * aura_part      = &meta.aura_part(); //no ghosts

  int p_rank = mesh.parallel_rank();

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    //      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    //      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
    //      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
  }
}

}
