// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stddef.h>                     // for NULL
#include <stdint.h>                     // for uint64_t
#include <algorithm>                    // for copy
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper, TEST

// CRW: this should be in BucketConnectivity.hpp, but circular dependency for now
#include "stk_mesh/base/BulkData.hpp"

#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_mesh/baseImpl/BucketConnDynamic.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class BulkData; } }

using namespace stk::mesh;

namespace {

typedef impl::BucketConnectivity<stk::topology::NODE_RANK, FIXED_CONNECTIVITY> fixed_conn;

void check_uninit_conn_size(fixed_conn& conn, unsigned num_conn, unsigned ordinal)
{
  EXPECT_EQ(conn.num_connectivity(ordinal), num_conn);
}

void check_even_conn_removed(fixed_conn& conn, unsigned num_conn, unsigned ordinal)
{
  EXPECT_EQ(conn.num_connectivity(ordinal), num_conn);

  Entity const* targets = conn.begin(ordinal);
  for (unsigned i = 0; i < num_conn; ++i) {
    Entity e_to(ordinal * num_conn + i + 1);
    if ( (i % 2) == 0 ) {
      EXPECT_EQ(targets[i], Entity());
    }
    else {
      EXPECT_EQ(targets[i], e_to);
    }
  }
}

template <typename Connectivity>
void test_simple_add(Connectivity& connectivity, unsigned num_entities_to_add, unsigned num_to_add)
{
  // Populate connectivity all at once for each entity

  EXPECT_EQ(connectivity.size(), 0u);

  for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
    connectivity.add_entity();

    EXPECT_EQ(connectivity.size(), ord + 1);
    check_uninit_conn_size(connectivity, num_to_add, ord);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity e_to(ord * num_to_add + i + 1);
      connectivity.add_connectivity(ord, e_to, static_cast<ConnectivityOrdinal>(i));
    }

    EXPECT_EQ(connectivity.num_connectivity(ord), num_to_add);

    Entity const* begin = connectivity.begin(ord);
    Entity const* end   = connectivity.end(ord);
    ConnectivityOrdinal const* begin_ord = connectivity.begin_ordinals(ord);
    ConnectivityOrdinal const* end_ord   = connectivity.end_ordinals(ord);

    EXPECT_EQ(end - begin, num_to_add);
    EXPECT_EQ(end_ord - begin_ord, num_to_add);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity expected_to(ord * num_to_add + i + 1);
      EXPECT_EQ(expected_to, begin[i]);
      EXPECT_EQ(static_cast<ConnectivityOrdinal>(i), begin_ord[i]);
    }
  }
}

template <typename Connectivity>
void test_complex_add(Connectivity& connectivity, unsigned num_entities_to_add, unsigned num_to_add)
{
  // Populate connectivity one at a time for each entity

  EXPECT_EQ(connectivity.size(), 0u);

  for (uint64_t i = 0; i < num_to_add; ++i) {
    for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
      if (i == 0) {
        connectivity.add_entity();
      }

      if (i == 0) {
        EXPECT_EQ(connectivity.size(), ord + 1);
      }
      else {
        EXPECT_EQ(connectivity.size(), num_entities_to_add);
      }

      Entity e_to(ord * num_to_add + i + 1);
      connectivity.add_connectivity(ord, e_to, static_cast<ConnectivityOrdinal>(i));
    }
  }

  for (unsigned ord = 0; ord < num_entities_to_add; ++ord) {
    EXPECT_EQ(connectivity.num_connectivity(ord), num_to_add);

    Entity const* begin = connectivity.begin(ord);
    Entity const* end   = connectivity.end(ord);
    ConnectivityOrdinal const* begin_ord = connectivity.begin_ordinals(ord);
    ConnectivityOrdinal const* end_ord   = connectivity.end_ordinals(ord);

    EXPECT_EQ(end - begin, num_to_add);
    EXPECT_EQ(end_ord - begin_ord, num_to_add);

    for (uint64_t i = 0; i < num_to_add; ++i) {
      Entity expected_to(ord * num_to_add + i + 1);
      EXPECT_EQ(expected_to, begin[i]);
      EXPECT_EQ(static_cast<ConnectivityOrdinal>(i), begin_ord[i]);
    }
  }
}

template <typename Connectivity>
void test_remove(Connectivity& connectivity, unsigned num_entities, unsigned num_to_add)
{
  test_simple_add(connectivity, num_entities, num_to_add);

  unsigned ord_to_remove_from = num_entities / 2;

  for (uint64_t i = 0; i < num_to_add; ++i) {
    Entity e_to(ord_to_remove_from * num_to_add + i + 1);
    if ( (i % 2) == 0 ) {
      bool rv = connectivity.remove_connectivity(ord_to_remove_from, e_to, static_cast<ConnectivityOrdinal>(i));
      EXPECT_TRUE(rv);
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

TEST(BucketConnectivity, fixed_complex_add)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

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

TEST(BucketConnectivity, fixed_intra_conn_copy)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

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

TEST(BucketConnectivity, fixed_mod_end)
{
  const unsigned num_to_add   = 8;
  const unsigned num_entities = 100;
  BulkData * bulk = NULL;
  fixed_conn conn(num_to_add);

  test_mod_end(conn, num_entities, num_to_add);
  conn.end_modification(bulk);
}

TEST(BucketConnDynamic, basic)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);
  unsigned maxOrdinal = bucketCapacity-1;
  conn.grow_if_necessary(maxOrdinal);
  for(unsigned i=0; i<bucketCapacity; ++i) {
    EXPECT_EQ(0u,  conn.num_connectivity(i));
    EXPECT_EQ(conn.begin(i), conn.end(i));
    EXPECT_EQ(conn.begin_ordinals(i), conn.end_ordinals(i));
    EXPECT_EQ(conn.begin_permutations(i), conn.end_permutations(i));
  }

#ifndef NDEBUG
  EXPECT_ANY_THROW(conn.num_connectivity(bucketCapacity));
  EXPECT_ANY_THROW(conn.begin(bucketCapacity));
  EXPECT_ANY_THROW(conn.begin_ordinals(bucketCapacity));
  EXPECT_ANY_THROW(conn.begin_permutations(bucketCapacity));
#endif
}

TEST(BucketConnDynamic, addConnectivity)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity(1);
  stk::mesh::ConnectivityOrdinal ordinal = 3;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity, ordinal));
  EXPECT_EQ(1u, conn.num_connectivity(0));
  ASSERT_EQ(1u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity, conn.begin(0)[0]);
  EXPECT_EQ(ordinal, conn.begin_ordinals(0)[0]);
}

TEST(BucketConnDynamic, addTwoSeparateConnectivities)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1);
  stk::mesh::ConnectivityOrdinal ordinal3 = 3;
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal3));

  stk::mesh::Entity entity11(11);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1;
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal1));

  EXPECT_EQ(1u, conn.num_connectivity(0));
  ASSERT_EQ(1u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(0)[0]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  ASSERT_EQ(1u, std::distance(conn.begin(1), conn.end(1)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(1), conn.end_ordinals(1)));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(1)[0]);
}

TEST(BucketConnDynamic, addTwoSeparateConnectivitiesInReverseOrder)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity11(11);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1;
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal1));

  stk::mesh::Entity entity1(1);
  stk::mesh::ConnectivityOrdinal ordinal3 = 3;
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal3));

  EXPECT_EQ(1u, conn.num_connectivity(0));
  ASSERT_EQ(1u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(0)[0]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  ASSERT_EQ(1u, std::distance(conn.begin(1), conn.end(1)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(1), conn.end_ordinals(1)));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(1)[0]);
}

TEST(BucketConnDynamic, addConnectivitiesOutOfOrderForSingleBktOrdinal)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::EntityVector entities = {Entity(5), Entity(9), Entity(7)};
  std::vector<stk::mesh::ConnectivityOrdinal> ordinals = {2, 7, 5};

  for(unsigned i=0; i<entities.size(); ++i) {
    EXPECT_TRUE(conn.add_connectivity(0, entities[i], ordinals[i]));
  }

  EXPECT_EQ(3u, conn.num_connectivity(0));
  ASSERT_EQ(3u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(3u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  const stk::mesh::Entity* connEntities = conn.begin(0);
  const stk::mesh::ConnectivityOrdinal* connOrdinals = conn.begin_ordinals(0);
  std::sort(entities.begin(), entities.end());
  std::sort(ordinals.begin(), ordinals.end());
  for(unsigned i=0; i<3u; ++i) {
    EXPECT_EQ(entities[i], connEntities[i]);
    EXPECT_EQ(ordinals[i], connOrdinals[i]);
  }
}

TEST(BucketConnDynamic, addConnectivitiesOutOfOrderRepeatedOrdinalsForSingleBktOrdinal)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::EntityVector entities = {
      stk::mesh::Entity(1), stk::mesh::Entity(8), stk::mesh::Entity(5), stk::mesh::Entity(9),
      stk::mesh::Entity(6), stk::mesh::Entity(3), stk::mesh::Entity(7), stk::mesh::Entity(11)
                                     };
  std::vector<stk::mesh::ConnectivityOrdinal> ordinals = {
      2, 1, 2, 2,
      1, 2, 1, 1
      };

  for(unsigned i=0; i<entities.size(); ++i) {
    EXPECT_TRUE(conn.add_connectivity(0, entities[i], ordinals[i]));
  }

  EXPECT_EQ(8u, conn.num_connectivity(0));
  ASSERT_EQ(8u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(8u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  const stk::mesh::Entity* connEntities = conn.begin(0);
  const stk::mesh::ConnectivityOrdinal* connOrdinals = conn.begin_ordinals(0);
  entities = {
      stk::mesh::Entity(6), stk::mesh::Entity(7), stk::mesh::Entity(8), stk::mesh::Entity(11),
      stk::mesh::Entity(1), stk::mesh::Entity(3), stk::mesh::Entity(5), stk::mesh::Entity(9)
             };
  std::sort(ordinals.begin(), ordinals.end());
  for(unsigned i=0; i<6u; ++i) {
    EXPECT_EQ(entities[i], connEntities[i]);
    EXPECT_EQ(ordinals[i], connOrdinals[i]);
  }
}

TEST(BucketConnDynamic, addDuplicateConnectivity_noOp)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::EntityVector entities = {Entity(5), Entity(9), Entity(7)};
  std::vector<stk::mesh::ConnectivityOrdinal> ordinals = {2, 7, 5};

  for(unsigned i=0; i<entities.size(); ++i) {
    EXPECT_TRUE(conn.add_connectivity(0, entities[i], ordinals[i]));
  }

  EXPECT_FALSE(conn.add_connectivity(0, entities[1], ordinals[1]));

  EXPECT_EQ(3u, conn.num_connectivity(0));
  ASSERT_EQ(3u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(3u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  const stk::mesh::Entity* connEntities = conn.begin(0);
  const stk::mesh::ConnectivityOrdinal* connOrdinals = conn.begin_ordinals(0);
  std::sort(entities.begin(), entities.end());
  std::sort(ordinals.begin(), ordinals.end());
  for(unsigned i=0; i<3u; ++i) {
    EXPECT_EQ(entities[i], connEntities[i]);
    EXPECT_EQ(ordinals[i], connOrdinals[i]);
  }
}

TEST(BucketConnDynamic, addTwoSeparateConnectivitiesThenAnotherForFirstBktOrdinal)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1);
  stk::mesh::ConnectivityOrdinal ordinal3 = 3;
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal3));

  stk::mesh::Entity entity11(11);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1;
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal1));

  stk::mesh::Entity entity2(2);
  stk::mesh::ConnectivityOrdinal ordinal4 = 4;
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal4));

  EXPECT_EQ(2u, conn.num_connectivity(0));
  ASSERT_EQ(2u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(2u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(entity2, conn.begin(0)[1]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(ordinal4, conn.begin_ordinals(0)[1]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  ASSERT_EQ(1u, std::distance(conn.begin(1), conn.end(1)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(1), conn.end_ordinals(1)));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(1)[0]);

  EXPECT_EQ(3u, conn.total_num_connectivity());
  EXPECT_EQ(1u, conn.num_unused_entries());

  conn.compress_connectivity();

  EXPECT_EQ(3u, conn.total_num_connectivity());
  EXPECT_EQ(0u, conn.num_unused_entries());
}

TEST(BucketConnDynamic, addFiveSeparateConnectivitiesThenAnotherForFirstThreeToTestCompress)
{
  constexpr unsigned bucketCapacity = 15;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1);
  stk::mesh::ConnectivityOrdinal ordinal4 = 4;
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal4));

  stk::mesh::Entity entity2(2);
  stk::mesh::ConnectivityOrdinal ordinal3 = 3;
  EXPECT_TRUE(conn.add_connectivity(1, entity2, ordinal3));

  stk::mesh::Entity entity21(21);
  stk::mesh::ConnectivityOrdinal ordinal31 = 31;
  EXPECT_TRUE(conn.add_connectivity(1, entity21, ordinal31));

  stk::mesh::Entity entity3(3);
  stk::mesh::ConnectivityOrdinal ordinal2 = 2;
  EXPECT_TRUE(conn.add_connectivity(2, entity3, ordinal2));

  stk::mesh::Entity entity4(4);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1;
  EXPECT_TRUE(conn.add_connectivity(3, entity4, ordinal1));

  stk::mesh::Entity entity5(5);
  stk::mesh::ConnectivityOrdinal ordinal0 = 0;
  EXPECT_TRUE(conn.add_connectivity(4, entity5, ordinal0));

  stk::mesh::Entity entity22(22);
  stk::mesh::ConnectivityOrdinal ordinal32 = 32;
  EXPECT_TRUE(conn.add_connectivity(1, entity22, ordinal32));

  stk::mesh::Entity entity31(31);
  stk::mesh::ConnectivityOrdinal ordinal21 = 21;
  EXPECT_TRUE(conn.add_connectivity(2, entity31, ordinal21));

  stk::mesh::Entity entity51(51);
  stk::mesh::ConnectivityOrdinal ordinal01 = 100;
  EXPECT_TRUE(conn.add_connectivity(4, entity51, ordinal01));

  EXPECT_EQ(1u, conn.num_connectivity(0));
  ASSERT_EQ(1u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(1u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal4, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(1)[0]);

  EXPECT_EQ(3u, conn.num_connectivity(1));
  ASSERT_EQ(3u, std::distance(conn.begin(1), conn.end(1)));
  ASSERT_EQ(3u, std::distance(conn.begin_ordinals(1), conn.end_ordinals(1)));
  EXPECT_EQ(entity2, conn.begin(1)[0]);
  EXPECT_EQ(entity21, conn.begin(1)[1]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(1)[0]);

  EXPECT_EQ(9u, conn.total_num_connectivity());
  EXPECT_EQ(4u, conn.num_unused_entries());

  conn.compress_connectivity();

  EXPECT_EQ(9u, conn.total_num_connectivity());
  EXPECT_EQ(0u, conn.num_unused_entries());

  EXPECT_EQ(3u, conn.num_connectivity(1));
  EXPECT_EQ(entity2, conn.begin(1)[0]);
  EXPECT_EQ(entity21, conn.begin(1)[1]);
  EXPECT_EQ(entity22, conn.begin(1)[2]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(1)[0]);
  EXPECT_EQ(ordinal31, conn.begin_ordinals(1)[1]);
  EXPECT_EQ(ordinal32, conn.begin_ordinals(1)[2]);

  EXPECT_EQ(2u, conn.num_connectivity(4));
  EXPECT_EQ(entity5, conn.begin(4)[0]);
  EXPECT_EQ(entity51, conn.begin(4)[1]);
  EXPECT_EQ(ordinal0, conn.begin_ordinals(4)[0]);
  EXPECT_EQ(ordinal01, conn.begin_ordinals(4)[1]);
}

TEST(BucketConnDynamic, addConnectivity_copyFirstToEnd_compress)
{
  constexpr unsigned bucketCapacity = 3;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1);
  stk::mesh::ConnectivityOrdinal ordinal4 = 4;
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal4));

  stk::mesh::Entity entity2(2);
  stk::mesh::ConnectivityOrdinal ordinal3 = 3;
  EXPECT_TRUE(conn.add_connectivity(1, entity2, ordinal3));

  stk::mesh::Entity entity3(3);
  stk::mesh::ConnectivityOrdinal ordinal2 = 2;
  EXPECT_TRUE(conn.add_connectivity(2, entity3, ordinal2));

  stk::mesh::Entity entity4(4);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1;
  EXPECT_TRUE(conn.add_connectivity(0, entity4, ordinal1));

  EXPECT_EQ(2u, conn.num_connectivity(0));
  ASSERT_EQ(2u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(2u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity4, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(ordinal4, conn.begin_ordinals(0)[1]);
  EXPECT_EQ(entity2, conn.begin(1)[0]);

  EXPECT_EQ(4u, conn.total_num_connectivity());
  EXPECT_EQ(1u, conn.num_unused_entries());

  conn.compress_connectivity();

  EXPECT_EQ(4u, conn.total_num_connectivity());
  EXPECT_EQ(0u, conn.num_unused_entries());

  EXPECT_EQ(2u, conn.num_connectivity(0));
  ASSERT_EQ(2u, std::distance(conn.begin(0), conn.end(0)));
  ASSERT_EQ(2u, std::distance(conn.begin_ordinals(0), conn.end_ordinals(0)));
  EXPECT_EQ(entity4, conn.begin(0)[0]);
  EXPECT_EQ(entity1, conn.begin(0)[1]);
  EXPECT_EQ(entity2, conn.begin(1)[0]);
  EXPECT_EQ(entity3, conn.begin(2)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(ordinal4, conn.begin_ordinals(0)[1]);
  EXPECT_EQ(ordinal3, conn.begin_ordinals(1)[0]);
  EXPECT_EQ(ordinal2, conn.begin_ordinals(2)[0]);
}

TEST(BucketConnDynamic, removeTheOnlyConnectivity)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2));
  EXPECT_TRUE(conn.add_connectivity(2, entity12, ordinal8));

  EXPECT_EQ(2u, conn.num_connectivity(0));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(entity2, conn.begin(0)[1]);
  EXPECT_EQ(ordinal2, conn.begin_ordinals(0)[1]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal7, conn.begin_ordinals(1)[0]);

  EXPECT_EQ(1u, conn.num_connectivity(2));
  EXPECT_EQ(entity12, conn.begin(2)[0]);
  EXPECT_EQ(ordinal8, conn.begin_ordinals(2)[0]);
}

TEST(BucketConnDynamic, addConnectivityInWeirdOrder)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2));
  EXPECT_TRUE(conn.add_connectivity(2, entity12, ordinal8));

  EXPECT_EQ(2u, conn.num_connectivity(0));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(entity2, conn.begin(0)[1]);
  EXPECT_EQ(ordinal2, conn.begin_ordinals(0)[1]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal7, conn.begin_ordinals(1)[0]);

  EXPECT_EQ(1u, conn.num_connectivity(2));
  EXPECT_EQ(entity12, conn.begin(2)[0]);
  EXPECT_EQ(ordinal8, conn.begin_ordinals(2)[0]);
}

TEST(BucketConnDynamic, addConnectivityInWeirdOrder_withPermutations)
{
  constexpr unsigned bucketCapacity = 10;
  const bool hasPermutations = true;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity, hasPermutations);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  stk::mesh::Permutation perm1 = static_cast<stk::mesh::Permutation>(1);
  stk::mesh::Permutation perm2 = static_cast<stk::mesh::Permutation>(2);
  stk::mesh::Permutation perm7 = static_cast<stk::mesh::Permutation>(7);
  stk::mesh::Permutation perm8 = static_cast<stk::mesh::Permutation>(8);
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1, perm1));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7, perm7));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2, perm2));
  EXPECT_TRUE(conn.add_connectivity(2, entity12, ordinal8, perm8));

  EXPECT_EQ(2u, conn.num_connectivity(0));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
  EXPECT_EQ(perm1, conn.begin_permutations(0)[0]);
  EXPECT_EQ(entity2, conn.begin(0)[1]);
  EXPECT_EQ(ordinal2, conn.begin_ordinals(0)[1]);
  EXPECT_EQ(perm2, conn.begin_permutations(0)[1]);

  EXPECT_EQ(1u, conn.num_connectivity(1));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal7, conn.begin_ordinals(1)[0]);
  EXPECT_EQ(perm7, conn.begin_permutations(1)[0]);

  EXPECT_EQ(1u, conn.num_connectivity(2));
  EXPECT_EQ(entity12, conn.begin(2)[0]);
  EXPECT_EQ(ordinal8, conn.begin_ordinals(2)[0]);
  EXPECT_EQ(perm8, conn.begin_permutations(2)[0]);
}

TEST(BucketConnDynamic, addAndRemoveMultipleConnectivity)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2));
  EXPECT_TRUE(conn.add_connectivity(1, entity12, ordinal8));

  EXPECT_TRUE(conn.remove_connectivity(0, entity2, ordinal2));
  EXPECT_EQ(1u, conn.num_connectivity(0));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);

  EXPECT_TRUE(conn.remove_connectivity(1, entity11, ordinal7));
  EXPECT_EQ(1u, conn.num_connectivity(1));
  EXPECT_EQ(entity12, conn.begin(1)[0]);
  EXPECT_EQ(ordinal8, conn.begin_ordinals(1)[0]);

  EXPECT_TRUE(conn.remove_connectivity(1, entity12, ordinal8));
  EXPECT_EQ(0u, conn.num_connectivity(1));
  EXPECT_EQ(conn.begin(1), conn.end(1));
  EXPECT_EQ(conn.begin_ordinals(1), conn.begin_ordinals(1));

  EXPECT_FALSE(conn.remove_connectivity(0, entity2, ordinal2));
  EXPECT_EQ(1u, conn.num_connectivity(0));
  EXPECT_EQ(entity1, conn.begin(0)[0]);
  EXPECT_EQ(ordinal1, conn.begin_ordinals(0)[0]);
}

TEST(BucketConnDynamic, removeAllConnectivityForBktOrdinal)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2));
  EXPECT_TRUE(conn.add_connectivity(1, entity12, ordinal8));

  EXPECT_TRUE(conn.remove_connectivity(0));
  EXPECT_EQ(0u, conn.num_connectivity(0));
  EXPECT_EQ(conn.begin(0), conn.end(0));
  EXPECT_EQ(conn.begin_ordinals(0), conn.begin_ordinals(0));

  EXPECT_EQ(2u, conn.num_connectivity(1));
  EXPECT_EQ(entity11, conn.begin(1)[0]);
  EXPECT_EQ(ordinal7, conn.begin_ordinals(1)[0]);
  EXPECT_EQ(entity12, conn.begin(1)[1]);
  EXPECT_EQ(ordinal8, conn.begin_ordinals(1)[1]);
}

TEST(BucketConnDynamic, makeHoleThenFillHoleWithoutRaisingCapacity)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  stk::mesh::Entity entity1(1), entity2(2), entity11(11), entity12(12);
  stk::mesh::Entity entity21(21), entity22(22);
  stk::mesh::ConnectivityOrdinal ordinal1 = 1, ordinal2 = 2, ordinal7 = 7, ordinal8 = 8;
  stk::mesh::ConnectivityOrdinal ordinal14 = 14, ordinal15 = 15;
  
  EXPECT_TRUE(conn.add_connectivity(0, entity1, ordinal1));
  EXPECT_TRUE(conn.add_connectivity(0, entity2, ordinal2));
  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(1, entity12, ordinal8));
  EXPECT_TRUE(conn.add_connectivity(2, entity21, ordinal14));
  EXPECT_TRUE(conn.add_connectivity(2, entity22, ordinal15));

  size_t totalCapacity = conn.total_capacity();
  size_t totalNumConnectivity = conn.total_num_connectivity();
  size_t numUnused = conn.num_unused_entries();

  EXPECT_TRUE(conn.remove_connectivity(1));
  EXPECT_EQ(0u, conn.num_connectivity(1));
  EXPECT_EQ(conn.begin(1), conn.end(1));
  EXPECT_EQ(conn.begin_ordinals(1), conn.begin_ordinals(1));

  EXPECT_TRUE(conn.add_connectivity(1, entity11, ordinal7));
  EXPECT_TRUE(conn.add_connectivity(1, entity12, ordinal8));

  EXPECT_TRUE(totalCapacity == conn.total_capacity());
  EXPECT_TRUE(totalNumConnectivity == conn.total_num_connectivity());
  EXPECT_TRUE(numUnused == conn.num_unused_entries());
}

TEST(BucketConnDynamic, replaceConnectivity_sameLength)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  std::vector<std::vector<stk::mesh::Entity>> entities = {
                           {stk::mesh::Entity(1),  stk::mesh::Entity(2)},
                           {stk::mesh::Entity(11), stk::mesh::Entity(12)},
                           {stk::mesh::Entity(21), stk::mesh::Entity(22)}
                                                         };
  std::vector<std::vector<stk::mesh::ConnectivityOrdinal>> ordinals = {
                                                               {1, 2},
                                                               {7, 8},
                                                               {14, 15}
                                                                 };
 
  for(unsigned bktOrdinal=0; bktOrdinal<entities.size(); ++bktOrdinal) {
    for(unsigned j=0; j<entities[bktOrdinal].size(); ++j) {
      EXPECT_TRUE(conn.add_connectivity(bktOrdinal, entities[bktOrdinal][j], ordinals[bktOrdinal][j]));
    }
  }

  std::vector<stk::mesh::Entity> newEntities = {Entity(31), Entity(32)};
  std::vector<stk::mesh::ConnectivityOrdinal> newOrdinals = {24, 25};
  const stk::mesh::Permutation* perms = nullptr;

  EXPECT_TRUE(conn.replace_connectivity(1, newEntities.size(), newEntities.data(), newOrdinals.data(), perms));
  EXPECT_EQ(2u, conn.num_connectivity(1));
  EXPECT_EQ(newEntities[0], conn.begin(1)[0]);
  EXPECT_EQ(newEntities[1], conn.begin(1)[1]);
  EXPECT_EQ(newOrdinals[0], conn.begin_ordinals(1)[0]);
  EXPECT_EQ(newOrdinals[1], conn.begin_ordinals(1)[1]);
}

TEST(BucketConnDynamic, replaceConnectivity_longer)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  std::vector<std::vector<stk::mesh::Entity>> entities = {
                           {stk::mesh::Entity(1),  stk::mesh::Entity(2)},
                           {stk::mesh::Entity(11), stk::mesh::Entity(12)},
                           {stk::mesh::Entity(21), stk::mesh::Entity(22)}
                                                         };
  std::vector<std::vector<stk::mesh::ConnectivityOrdinal>> ordinals = {
                                                               {1, 2},
                                                               {7, 8},
                                                               {14, 15}
                                                                 };
 
  for(unsigned bktOrdinal=0; bktOrdinal<entities.size(); ++bktOrdinal) {
    for(unsigned j=0; j<entities[bktOrdinal].size(); ++j) {
      EXPECT_TRUE(conn.add_connectivity(bktOrdinal, entities[bktOrdinal][j], ordinals[bktOrdinal][j]));
    }
  }

  std::vector<stk::mesh::Entity> newEntities = {Entity(31), Entity(32), Entity(33)};
  std::vector<stk::mesh::ConnectivityOrdinal> newOrdinals = {24, 25, 26};
  const stk::mesh::Permutation* perms = nullptr;

  EXPECT_TRUE(conn.replace_connectivity(1, newEntities.size(), newEntities.data(), newOrdinals.data(), perms));
  EXPECT_EQ(3u, conn.num_connectivity(1));
  EXPECT_EQ(newEntities[0], conn.begin(1)[0]);
  EXPECT_EQ(newEntities[1], conn.begin(1)[1]);
  EXPECT_EQ(newEntities[2], conn.begin(1)[2]);
  EXPECT_EQ(newOrdinals[0], conn.begin_ordinals(1)[0]);
  EXPECT_EQ(newOrdinals[1], conn.begin_ordinals(1)[1]);
  EXPECT_EQ(newOrdinals[2], conn.begin_ordinals(1)[2]);
}

TEST(BucketConnDynamic, replaceConnectivity_shorter)
{
  constexpr unsigned bucketCapacity = 10;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  std::vector<std::vector<stk::mesh::Entity>> entities = {
                           {stk::mesh::Entity(1),  stk::mesh::Entity(2)},
                           {stk::mesh::Entity(11), stk::mesh::Entity(12)},
                           {stk::mesh::Entity(21), stk::mesh::Entity(22)}
                                                         };
  std::vector<std::vector<stk::mesh::ConnectivityOrdinal>> ordinals = {
                                                               {1, 2},
                                                               {7, 8},
                                                               {14, 15}
                                                                 };
 
  for(unsigned bktOrdinal=0; bktOrdinal<entities.size(); ++bktOrdinal) {
    for(unsigned j=0; j<entities[bktOrdinal].size(); ++j) {
      EXPECT_TRUE(conn.add_connectivity(bktOrdinal, entities[bktOrdinal][j], ordinals[bktOrdinal][j]));
    }
  }

  std::vector<stk::mesh::Entity> newEntities = {Entity(31)};
  std::vector<stk::mesh::ConnectivityOrdinal> newOrdinals = {24};
  const stk::mesh::Permutation* perms = nullptr;

  EXPECT_TRUE(conn.replace_connectivity(1, newEntities.size(), newEntities.data(), newOrdinals.data(), perms));
  EXPECT_EQ(1u, conn.num_connectivity(1));
  EXPECT_EQ(newEntities[0], conn.begin(1)[0]);
  EXPECT_EQ(newOrdinals[0], conn.begin_ordinals(1)[0]);
}

void add_100k_connectivities()
{
  constexpr unsigned bucketCapacity = 50000;
  stk::mesh::impl::BucketConnDynamic conn(bucketCapacity);

  for(unsigned i=0; i<bucketCapacity; ++i) {
    conn.add_connectivity(i, stk::mesh::Entity(1), stk::mesh::ConnectivityOrdinal(0));
    conn.add_connectivity(i, stk::mesh::Entity(2), stk::mesh::ConnectivityOrdinal(1));
  }
}

TEST(BucketConnDynamic, test_65k_limit)
{
#ifdef STK_16BIT_UPWARDCONN_INDEX_TYPE
  EXPECT_ANY_THROW(add_100k_connectivities());
#else
  EXPECT_NO_THROW(add_100k_connectivities());
#endif
}

