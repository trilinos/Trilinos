#include <iostream>
#include <sstream>

#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace {

template <stk::mesh::ConnectivityType Type>
void do_test(stk::mesh::impl::BucketConnectivity<stk::topology::NODE_RANK, Type>& connectivity,
             stk::mesh::impl::BucketConnectivity<stk::topology::NODE_RANK, Type>& connectivity2,
             unsigned num_from_entities,
             unsigned num_node_connectivity)
{
  using namespace stk::mesh;

    // Start entities at 1 b/c still discussing whether Entity will be POD.

  const unsigned from_entities_start_num = 1;
  std::vector<Entity> from_entities(num_from_entities);
  for (unsigned i = 0; i < num_from_entities; ++i)
  {
    from_entities[i].set_local_offset(i + from_entities_start_num);
    connectivity.add_entity();
  }

  const unsigned num_to_entities = num_from_entities * num_node_connectivity;
  const unsigned to_entities_start_num = 100;
  std::vector<Entity> to_entities(num_to_entities);
  std::vector<Node> nodes(num_to_entities);
  for (unsigned i = 0; i < num_to_entities; ++i)
  {
    to_entities[i].set_local_offset(i + to_entities_start_num);
    nodes[i] = static_cast<Node>(i);
  }

  //assign initial connectivities
  {
    Permutation dummy_perm = INVALID_PERMUTATION;
    for (unsigned i = 0; i < num_from_entities; ++i)
    {
      for (unsigned ord = 0; ord < num_node_connectivity; ++ord)
      {
        connectivity.add_connectivity(i, to_entities[i * num_node_connectivity + ord],
                                      static_cast<ConnectivityOrdinal>(ord), dummy_perm);
      }
    }
  }

  STKUNIT_EXPECT_EQ(num_node_connectivity * num_from_entities,
                    connectivity.end(num_from_entities - 1) - connectivity.begin(0));

  {
    unsigned ilc = 0;
    for (unsigned i = 0; i < num_from_entities; ++i)
    {
      const Entity * b = connectivity.begin_entities(i);
      const Entity * e = connectivity.end_entities(i);
      const ConnectivityOrdinal * cb = connectivity.begin_ordinals(i);
      EXPECT_EQ(e - b, num_node_connectivity);
      for (unsigned ord = 0; b<e; ++b, ++ord, ++cb, ++ilc)
      {
        Entity to = to_entities[i * num_node_connectivity + ord];
        EXPECT_EQ(*b,  to);
        EXPECT_EQ(*cb, ord);
      }
    }
    STKUNIT_EXPECT_EQ(ilc, num_from_entities * num_node_connectivity);
  }

  //reverse connectivities by swapping first with last

  {
    for (unsigned i = 0, e = num_from_entities - 1; i < e; ++i, --e) {
      connectivity.swap(i, connectivity, e);
    }
  }

  {
    unsigned ilc = 0;
    for (unsigned i = 0; i < num_from_entities; ++i)
    {
      unsigned to_idx = (num_from_entities - i - 1) * num_node_connectivity;
      const Entity * b = connectivity.begin_entities(i);
      const Entity * e = connectivity.end_entities(i);
      const ConnectivityOrdinal * cb = connectivity.begin_ordinals(i);
      EXPECT_EQ(e - b, num_node_connectivity);
      for (unsigned ord = 0; b<e; ++b, ++ord, ++cb, ++ilc)
      {
        Entity to = to_entities[to_idx++];
        EXPECT_EQ(*b, to);
        EXPECT_EQ(*cb, ord);
      }
    }
    STKUNIT_EXPECT_EQ(ilc, num_from_entities * num_node_connectivity);
  }

  //move entities to new instance

  for (unsigned i = 0; i < num_from_entities; ++i)
  {
    connectivity.move_entity(connectivity2);
  }

  {
    unsigned ilc = 0;
    for (unsigned i = 0; i < num_from_entities; ++i)
    {
      const Entity * b = connectivity2.begin_entities(i);
      const Entity * e = connectivity2.end_entities(i);
      const ConnectivityOrdinal * cb = connectivity2.begin_ordinals(i);
      EXPECT_EQ(e - b, num_node_connectivity);
      for (unsigned ord = 0; b<e; ++b, ++ord, ++cb, ++ilc)
      {
        Entity to = to_entities[i * num_node_connectivity + ord];
        EXPECT_EQ(*b, to);
        EXPECT_EQ(*cb, ord);
      }
    }
    STKUNIT_EXPECT_EQ(ilc, num_from_entities * num_node_connectivity);
  }

  //move back to original partition

  for (unsigned i = 0; i < num_from_entities; ++i)
  {
    connectivity2.move_entity(connectivity);
  }

  {
    unsigned ilc = 0;
    for (unsigned i = 0; i < num_from_entities; ++i)
    {
      unsigned to_idx = (num_from_entities - i - 1) * num_node_connectivity;
      const Entity * b = connectivity.begin_entities(i);
      const Entity * e = connectivity.end_entities(i);
      const ConnectivityOrdinal * cb = connectivity.begin_ordinals(i);
      EXPECT_EQ(e - b, num_node_connectivity);
      for (unsigned ord = 0; b<e; ++b, ++ord, ++cb, ++ilc)
      {
        Entity to = to_entities[to_idx++];
        EXPECT_EQ(*b, to);
        EXPECT_EQ(*cb, ord);
      }
    }
    STKUNIT_EXPECT_EQ(ilc, num_from_entities * num_node_connectivity);
  }
}

}

// Ported from stk_samba unit tests.
TEST(BucketConnectivity, Fixed)
{
  using namespace stk::mesh;

  typedef impl::BucketConnectivity<stk::topology::NODE_RANK, FIXED_CONNECTIVITY> fixed_node_type;

  const unsigned num_from_entities = 3;
  const unsigned num_node_connectivity = 4;

  fixed_node_type nodes_rel(num_node_connectivity);
  fixed_node_type nodes_rel2(num_node_connectivity);

  do_test(nodes_rel, nodes_rel2, num_from_entities, num_node_connectivity);
}

TEST(BucketConnectivity, Dynamic)
{
  using namespace stk::mesh;

  typedef impl::BucketConnectivity<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY> dynamic_node_type;

  const unsigned num_from_entities = 3;
  const unsigned num_node_connectivity = 4;

  dynamic_node_type nodes_rel(stk::topology::ELEMENT_RANK, 0);
  dynamic_node_type nodes_rel2(stk::topology::ELEMENT_RANK, 0);

  do_test(nodes_rel, nodes_rel2, num_from_entities, num_node_connectivity);
}

