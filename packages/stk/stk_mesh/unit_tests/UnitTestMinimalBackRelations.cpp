/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

namespace {

using namespace stk::mesh;

fixtures::HexFixture* set_up_mesh(ConnectivityMap const& conn_map)
{
  const unsigned NX = 9;
  const unsigned NY = 9;
  const unsigned NZ = 9;

  fixtures::HexFixture* rv = new fixtures::HexFixture(MPI_COMM_WORLD, NX, NY, NZ, &conn_map);
  rv->m_meta.commit();
  rv->generate_mesh();

  stk::mesh::skin_mesh(rv->m_bulk_data);
  stk::mesh::create_edges(rv->m_bulk_data);

  return rv;
}

struct SortComparator
{
  bool operator()(std::pair<Entity, ConnectivityOrdinal> const& lhs, std::pair<Entity, ConnectivityOrdinal> const& rhs) const
  {
    return impl::HigherConnectivityCompare()(lhs.first, lhs.second, rhs.first, rhs.second);
  }
};

void check_equiv_conn(Bucket const& bucket_full_conn, Bucket const& bucket_min_conn, size_t ord, EntityRank rank)
{
  BulkData& mesh_min_conn  = bucket_min_conn.mesh();

  ThrowRequire(!mesh_min_conn.connectivity_map().valid(bucket_min_conn.entity_rank(), rank));

  EntityVector temp_entities;
  std::vector<ConnectivityOrdinal> temp_ordinals;
  Entity const* rel_entities_min;
  ConnectivityOrdinal const* rel_ordinals_min;
  size_t num_min_upward = get_connectivity(mesh_min_conn,
                                           bucket_min_conn[ord],
                                           rank,
                                           temp_entities,
                                           temp_ordinals);
  rel_entities_min = &*temp_entities.begin();
  rel_ordinals_min = &*temp_ordinals.begin();

  STKUNIT_ASSERT_EQ(bucket_full_conn.num_connectivity(ord, rank), num_min_upward);

  Entity const* rel_entities_full              = bucket_full_conn.begin(ord, rank);
  ConnectivityOrdinal const* rel_ordinals_full = bucket_full_conn.begin_ordinals(ord, rank);

  // NOTE: computed back-connectivity may not be returned in the same order as it would
  //       be if it were stored, so we have to sort
  std::vector< std::pair<Entity, ConnectivityOrdinal> > temp;
  for (size_t i = 0; i < num_min_upward; ++i) {
    temp.push_back( std::make_pair( rel_entities_min[i], rel_ordinals_min[i] ) );
  }
  std::sort(temp.begin(), temp.end(), SortComparator());

  for (size_t i = 0; i < num_min_upward; ++i) {
    temp_entities[i] = temp[i].first;
    temp_ordinals[i] = temp[i].second;
  }

  for (size_t i = 0; i < num_min_upward; ++i) {
    STKUNIT_EXPECT_EQ( rel_entities_min[i], rel_entities_full[i] );
    STKUNIT_EXPECT_EQ( rel_ordinals_min[i], rel_ordinals_full[i] );
  }
}

STKUNIT_UNIT_TEST( UnitTestMinimalBackRelation, simpleHex )
{
  fixtures::HexFixture* fixture_with_full_conn = set_up_mesh(ConnectivityMap::classic_stk_mesh());
  fixtures::HexFixture* fixture_with_min_conn  = set_up_mesh(ConnectivityMap::minimal_back_relations_map());

  BulkData& mesh_full_conn = fixture_with_full_conn->m_bulk_data;
  BulkData& mesh_min_conn  = fixture_with_min_conn->m_bulk_data;

  {
    for (stk::topology::rank_t rank = stk::topology::NODE_RANK; rank < stk::topology::ELEMENT_RANK; ++rank) {
      const BucketVector & buckets_full_conn = mesh_full_conn.buckets(rank);
      const BucketVector & buckets_min_conn  = mesh_min_conn.buckets(rank);
      STKUNIT_ASSERT_EQ(buckets_full_conn.size(), buckets_min_conn.size());

      for (size_t ib=0, endb=buckets_full_conn.size(); ib < endb; ++ib) {
        const Bucket & bucket_full_conn = *buckets_full_conn[ib];
        const Bucket & bucket_min_conn  = *buckets_min_conn[ib];
        STKUNIT_ASSERT_EQ(bucket_full_conn.size(), bucket_min_conn.size());

        for (size_t ord=0, end=bucket_full_conn.size(); ord<end; ++ord) {
          if ( rank > stk::topology::NODE_RANK ) {
            STKUNIT_EXPECT_EQ(bucket_min_conn.num_elements(ord),0u); // no stored back-rels to elements except for nodes
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::ELEMENT_RANK);
          }
          if ( rank < stk::topology::FACE_RANK) {
            STKUNIT_EXPECT_EQ(bucket_min_conn.num_faces(ord),0u);    // no stored back-rels to faces
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::FACE_RANK);
          }
          if ( rank < stk::topology::EDGE_RANK) {
            STKUNIT_EXPECT_EQ(bucket_min_conn.num_edges(ord),0u);    // no stored back-rels to edges
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::EDGE_RANK);
          }

          // Check that all downward rels are the same
          for (EntityRank irank = stk::topology::NODE_RANK; irank < rank; ++irank) {
            STKUNIT_EXPECT_EQ(bucket_full_conn.num_connectivity(ord, irank), bucket_min_conn.num_connectivity(ord, irank));
          }
        }
      }
    }
  }

  delete fixture_with_full_conn;
  delete fixture_with_min_conn;
}

} // namespace
