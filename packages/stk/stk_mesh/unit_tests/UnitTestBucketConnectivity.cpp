#include <stddef.h>                     // for NULL
#include <stdint.h>                     // for uint64_t
#include <algorithm>                    // for copy
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper, TEST

// CRW: this should be in BucketConnectivity.hpp, but circular dependency for now
#include "stk_mesh/base/BulkData.hpp"

#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class BulkData; } }

using namespace stk::mesh;

namespace {

typedef impl::BucketConnectivity<stk::topology::NODE_RANK, FIXED_CONNECTIVITY> fixed_conn;
typedef impl::BucketConnectivity<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY> dynamic_conn;

void check_uninit_conn_size(fixed_conn& conn, unsigned num_conn, unsigned ordinal)
{
  STKUNIT_EXPECT_EQ(conn.num_connectivity(ordinal), num_conn);
}

void check_uninit_conn_size(dynamic_conn& conn, unsigned num_conn, unsigned ordinal)
{
  STKUNIT_EXPECT_EQ(conn.num_connectivity(ordinal), 0u);
}

void check_even_conn_removed(fixed_conn& conn, unsigned num_conn, unsigned ordinal)
{
  STKUNIT_EXPECT_EQ(conn.num_connectivity(ordinal), num_conn);

  Entity const* targets = conn.begin(ordinal);
  for (unsigned i = 0; i < num_conn; ++i) {
    Entity e_to = {ordinal * num_conn + i + 1};
    if ( (i % 2) == 0 ) {
      STKUNIT_EXPECT_EQ(targets[i], Entity());
    }
    else {
      STKUNIT_EXPECT_EQ(targets[i], e_to);
    }
  }
}

void check_even_conn_removed(dynamic_conn& conn, unsigned num_conn, unsigned ordinal)
{
  STKUNIT_EXPECT_EQ(conn.num_connectivity(ordinal), num_conn / 2);

  Entity const* targets = conn.begin(ordinal);
  ConnectivityOrdinal const* ordinals = conn.begin_ordinals(ordinal);
  for (unsigned i = 0; i < num_conn / 2; ++i) {
    Entity e_to = {ordinal * num_conn + ((2*i) + 1) + 1};
    STKUNIT_EXPECT_EQ(targets[i], e_to);
    STKUNIT_EXPECT_EQ(ordinals[i], static_cast<ConnectivityOrdinal>((2*i) + 1));
  }
}

template <typename Connectivity>
void test_simple_add(Connectivity& connectivity, unsigned num_entities_to_add, unsigned num_to_add)
{
  // Populate connectivity all at once for each entity

  STKUNIT_EXPECT_EQ(connectivity.size(), 0u);

  for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
    connectivity.add_entity();

    STKUNIT_EXPECT_EQ(connectivity.size(), ord + 1);
    check_uninit_conn_size(connectivity, num_to_add, ord);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity e_to = {ord * num_to_add + i + 1};
      connectivity.add_connectivity(ord, e_to, static_cast<ConnectivityOrdinal>(i));
    }

    STKUNIT_EXPECT_EQ(connectivity.num_connectivity(ord), num_to_add);

    Entity const* begin = connectivity.begin(ord);
    Entity const* end   = connectivity.end(ord);
    ConnectivityOrdinal const* begin_ord = connectivity.begin_ordinals(ord);
    ConnectivityOrdinal const* end_ord   = connectivity.end_ordinals(ord);

    STKUNIT_EXPECT_EQ(end - begin, num_to_add);
    STKUNIT_EXPECT_EQ(end_ord - begin_ord, num_to_add);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity expected_to = {ord * num_to_add + i + 1};
      STKUNIT_EXPECT_EQ(expected_to, begin[i]);
      STKUNIT_EXPECT_EQ(static_cast<ConnectivityOrdinal>(i), begin_ord[i]);
    }
  }
}

template <typename Connectivity>
void test_complex_add(Connectivity& connectivity, unsigned num_entities_to_add, unsigned num_to_add)
{
  // Populate connectivity one at a time for each entity

  STKUNIT_EXPECT_EQ(connectivity.size(), 0u);

  for (uint64_t i = 0; i < num_to_add; ++i) {
    for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
      if (i == 0) {
        connectivity.add_entity();
      }

      if (i == 0) {
        STKUNIT_EXPECT_EQ(connectivity.size(), ord + 1);
      }
      else {
        STKUNIT_EXPECT_EQ(connectivity.size(), num_entities_to_add);
      }

      Entity e_to = {ord * num_to_add + i + 1};
      connectivity.add_connectivity(ord, e_to, static_cast<ConnectivityOrdinal>(i));
    }
  }

  for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
    STKUNIT_EXPECT_EQ(connectivity.num_connectivity(ord), num_to_add);

    Entity const* begin = connectivity.begin(ord);
    Entity const* end   = connectivity.end(ord);
    ConnectivityOrdinal const* begin_ord = connectivity.begin_ordinals(ord);
    ConnectivityOrdinal const* end_ord   = connectivity.end_ordinals(ord);

    STKUNIT_EXPECT_EQ(end - begin, num_to_add);
    STKUNIT_EXPECT_EQ(end_ord - begin_ord, num_to_add);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity expected_to = {ord * num_to_add + i + 1};
      STKUNIT_EXPECT_EQ(expected_to, begin[i]);
      STKUNIT_EXPECT_EQ(static_cast<ConnectivityOrdinal>(i), begin_ord[i]);
    }
  }
}

template <typename Connectivity>
void test_remove(Connectivity& connectivity, unsigned num_entities, unsigned num_to_add)
{
  test_simple_add(connectivity, num_entities, num_to_add);

  unsigned ord_to_remove_from = num_entities / 2;

  for (uint64_t i = 0; i < num_to_add; ++i) {
    Entity e_to = {ord_to_remove_from * num_to_add + i + 1};
    if ( (i % 2) == 0 ) {
      bool rv = connectivity.remove_connectivity(ord_to_remove_from, e_to, static_cast<ConnectivityOrdinal>(i));
      STKUNIT_EXPECT_TRUE(rv);
    }
  }

  check_even_conn_removed(connectivity, num_to_add, ord_to_remove_from);
}

template <typename Connectivity>
void test_inter_conn_copy(Connectivity& connectivity, unsigned num_entities, unsigned num_to_add)
{
  // TODO
}

template <typename Connectivity>
void test_intra_conn_copy(Connectivity& connectivity, unsigned num_entities, unsigned num_to_add)
{
  // TODO
}

template <typename Connectivity>
void test_mod_end(Connectivity& connectivity, unsigned num_entities, unsigned num_to_add)
{
  // TODO
}

}

TEST(BucketConnectivity, fixed_simple_add)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_simple_add(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_simple_add)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_simple_add(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, fixed_complex_add)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_complex_add(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_complex_add)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_complex_add(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, fixed_remove)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_remove(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_remove)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_remove(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, fixed_intra_conn_copy)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_intra_conn_copy(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_intra_conn_copy)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_intra_conn_copy(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, fixed_inter_conn_copy)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_inter_conn_copy(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_inter_conn_copy)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_inter_conn_copy(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, fixed_mod_end)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_mod_end(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnectivity, dynamic_mod_end)
{
  const unsigned num_to_add = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  dynamic_conn conn(stk::topology::ELEMENT_RANK, bulk);

  test_mod_end(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}
