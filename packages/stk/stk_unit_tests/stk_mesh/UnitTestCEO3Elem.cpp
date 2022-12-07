#include <UnitTestCEO3Elem.hpp>
#include "UnitTestCEOCommonUtils.hpp"
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE

namespace CEOUtils
{
//////////////////////////////////// 3Elem2ProcMoveRight //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveRightAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector& elements)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

  stk::mesh::Part* elem_part = &meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part* node_part = &meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  int p_rank = mesh.parallel_rank();
  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  if(p_rank == 0)
  {
    nodes.push_back(mesh.declare_node(1, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(2, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(3, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(4, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(5, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(6, stk::mesh::ConstPartVector{node_part}));
    elements.push_back(mesh.declare_element(1, stk::mesh::ConstPartVector{elem_part}));
    elements.push_back(mesh.declare_element(2, stk::mesh::ConstPartVector{elem_part}));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.declare_relation(elements[1], nodes[2], 0);
    mesh.declare_relation(elements[1], nodes[3], 1);
    mesh.declare_relation(elements[1], nodes[4], 2);
    mesh.declare_relation(elements[1], nodes[5], 3);

    mesh.add_node_sharing(nodes[4], 1);
    mesh.add_node_sharing(nodes[5], 1);
  }
  else if(p_rank == 1)
  {
    nodes.push_back(mesh.declare_node(5, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(6, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(7, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(8, stk::mesh::ConstPartVector{node_part}));
    elements.push_back(mesh.declare_element(3, stk::mesh::ConstPartVector{elem_part}));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.add_node_sharing(nodes[0], 0);
    mesh.add_node_sharing(nodes[1], 0);
  }

  mesh.modification_end();

  stk::mesh::Part * quad_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta_data.universal_part();
  Part * owned_part     = &meta_data.locally_owned_part();
  Part * shared_part    = &meta_data.globally_shared_part();
  Part * aura_part      = &meta_data.aura_part();

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

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
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

void checkStatesAfterCEO_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

  stk::mesh::MetaData & meta_data = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta_data.get_part("elem_part");
  stk::mesh::Part * node_part = meta_data.get_part("node_part");
  stk::mesh::Part * quad_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta_data.universal_part();
  Part * owned_part     = &meta_data.locally_owned_part();
  Part * shared_part    = &meta_data.globally_shared_part();

  int p_rank = mesh.parallel_rank();
  if(p_rank == 0)
  {
    // /**/ is the delta from after modification_end to before modification_end

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  }
  else if(p_rank == 1)
  {

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

void checkStatesAfterCEOME_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

  stk::mesh::MetaData & meta_data = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta_data.get_part("elem_part");
  stk::mesh::Part * node_part = meta_data.get_part("node_part");
  stk::mesh::Part * quad_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta_data.universal_part();
  Part * owned_part     = &meta_data.locally_owned_part();
  Part * shared_part    = &meta_data.globally_shared_part();
  Part * aura_part      = &meta_data.aura_part();

  int p_rank = mesh.parallel_rank();
  if(p_rank == 0)
  {

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  }
  else if(p_rank == 1)
  {

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

//////////////////////////////////// 3Elem2ProcMoveLeft //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveLeftAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector &elements)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

  int p_rank = mesh.parallel_rank();

  stk::mesh::Part* elem_part = & meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part* node_part = & meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  if(p_rank == 0)
  {
    nodes.push_back(mesh.declare_node(1, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(2, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(3, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(4, stk::mesh::ConstPartVector{node_part}));
    elements.push_back(mesh.declare_element(1, stk::mesh::ConstPartVector{elem_part}));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.add_node_sharing(nodes[2], 1);
    mesh.add_node_sharing(nodes[3], 1);
  }
  else if(p_rank == 1)
  {
    nodes.push_back(mesh.declare_node(3, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(4, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(5, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(6, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(7, stk::mesh::ConstPartVector{node_part}));
    nodes.push_back(mesh.declare_node(8, stk::mesh::ConstPartVector{node_part}));
    elements.push_back(mesh.declare_element(2, stk::mesh::ConstPartVector{elem_part}));
    elements.push_back(mesh.declare_element(3, stk::mesh::ConstPartVector{elem_part}));

    mesh.declare_relation(elements[0], nodes[0], 0);
    mesh.declare_relation(elements[0], nodes[1], 1);
    mesh.declare_relation(elements[0], nodes[2], 2);
    mesh.declare_relation(elements[0], nodes[3], 3);

    mesh.declare_relation(elements[1], nodes[2], 0);
    mesh.declare_relation(elements[1], nodes[3], 1);
    mesh.declare_relation(elements[1], nodes[4], 2);
    mesh.declare_relation(elements[1], nodes[5], 3);

    mesh.add_node_sharing(nodes[0], 0);
    mesh.add_node_sharing(nodes[1], 0);
  }

  mesh.modification_end();

  stk::mesh::Part * quad_topo_part = &meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta_data.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta_data.universal_part();
  Part * owned_part     = &meta_data.locally_owned_part();
  Part * shared_part    = &meta_data.globally_shared_part();
  Part * aura_part      = &meta_data.aura_part();

  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

void checkStatesAfterCEO_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * node_part = meta.get_part("node_part");
  stk::mesh::Part * quad_topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  // Part * aura_part      = &meta.aura_part();

  int p_rank = mesh.parallel_rank();
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM)); ///////
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM)); ///////
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
    /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

void checkStatesAfterCEOME_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part * elem_part = meta.get_part("elem_part");
  stk::mesh::Part * node_part = meta.get_part("node_part");
  stk::mesh::Part * quad_topo_part = &meta.get_topology_root_part(stk::topology::QUAD_4_2D);
  stk::mesh::Part * node_topo_part = &meta.get_topology_root_part(stk::topology::NODE);
  Part * universal_part = &meta.universal_part();
  Part * owned_part     = &meta.locally_owned_part();
  Part * shared_part    = &meta.globally_shared_part();
  Part * aura_part      = &meta.aura_part();

  int p_rank = mesh.parallel_rank();
  if(p_rank == 0)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
  else if(p_rank == 1)
  {
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

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
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0));
    EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
    EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
  }
}

}


