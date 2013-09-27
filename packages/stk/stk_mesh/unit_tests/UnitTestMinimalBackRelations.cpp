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


STKUNIT_UNIT_TEST( UnitTestMinimalBackRelation, simpleHex )
{
  // TODO:  Make this test pass in parallel
  const unsigned NX = 9;
  const unsigned NY = 9;
  const unsigned NZ = 9;

  ConnectivityMap connectivity_map = ConnectivityMap::minimal_back_relations_map();
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD,NX,NY,NZ,&connectivity_map);
  fixture.m_meta.commit();
  fixture.generate_mesh();

  stk::mesh::skin_mesh(fixture.m_bulk_data, stk::topology::ELEMENT_RANK, NULL);
  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    for (stk::topology::rank_t rank = stk::topology::NODE_RANK; rank < stk::topology::ELEMENT_RANK; ++rank) {
      const BucketVector & buckets = fixture.m_bulk_data.buckets(rank);
      for (size_t ib=0, endb=buckets.size(); ib < endb; ++ib) {
        const Bucket & b = *buckets[ib];
        for (size_t ord=0, end=b.size(); ord<end; ++ord) {
          if ( rank > stk::topology::NODE_RANK) STKUNIT_EXPECT_EQ(b.num_elements(ord),0u);
          if ( rank <= stk::topology::FACE_RANK) STKUNIT_EXPECT_EQ(b.num_faces(ord),0u);
          if ( rank <= stk::topology::EDGE_RANK) STKUNIT_EXPECT_EQ(b.num_edges(ord),0u);
        }
      }
    }
  }
}

} // namespace
