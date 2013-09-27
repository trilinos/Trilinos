
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>

#include <iomanip>
#include <algorithm>
#include <sstream>

using stk::mesh::MetaData;

STKUNIT_UNIT_TEST ( UnitTestCreateEdges, Quad_2x1 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 2, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 6u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 0u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 2u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 6u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 7u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 2u ); // elements
  }
}

STKUNIT_UNIT_TEST ( UnitTestCreateEdges, Quad_3x1 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 3, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 8u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 0u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 3u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 8u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 10u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 3u ); // elements
  }
}

STKUNIT_UNIT_TEST ( UnitTestCreateEdges, Quad_2x2 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 2, 2);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  stk::mesh::skin_mesh(fixture.m_bulk_data, stk::topology::ELEMENT_RANK, NULL);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 9u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 8u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 4u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 9u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 12u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 4u ); // elements
  }

  //shouldn't do anything
  stk::mesh::create_edges(fixture.m_bulk_data);
}

STKUNIT_UNIT_TEST ( UnitTestCreateEdges, Hex_2x1x1 )
{
  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 1, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 12u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 0u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 2u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 12u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 20u ); // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u ); // faces
    STKUNIT_EXPECT_EQ( counts[3] , 2u ); // elements
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateEdges , Hex_3x1x1 )
{
  using namespace stk::mesh;

  fixtures::HexFixture fixture(MPI_COMM_WORLD, 3, 1, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 16u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u );  // faces
    STKUNIT_EXPECT_EQ( counts[3] , 3u ); // elements
  }

  create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 16u );
    STKUNIT_EXPECT_EQ( counts[1] , 28u );
    STKUNIT_EXPECT_EQ( counts[2] , 0u );
    STKUNIT_EXPECT_EQ( counts[3] , 3u );
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateEdges , testCreateEdges3x3x3 )
{
  const size_t elem_rank = MetaData::ELEMENT_RANK;
  const size_t face_rank = MetaData::FACE_RANK;
  const size_t edge_rank = MetaData::EDGE_RANK;
  const size_t node_rank = MetaData::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[node_rank] , 64u ); // nodes
    STKUNIT_EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[face_rank] , 0u );  // faces
    STKUNIT_EXPECT_EQ( counts[elem_rank] , 27u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( 64u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 144u, counts[edge_rank] );  // edges
    STKUNIT_EXPECT_EQ( 0u, counts[face_rank] );  // faces
    STKUNIT_EXPECT_EQ( 27u, counts[elem_rank]  ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( 0u, b.num_faces(i) );
      STKUNIT_EXPECT_EQ( 12u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( 8u,  b.num_nodes(i) );
    }
  }

  stk::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = face_buckets.begin();
       b_itr != face_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( 4u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( 4u, b.num_nodes(i) );
    }
  }

  stk::mesh::BucketVector  edge_buckets = fixture.m_bulk_data.buckets(edge_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = edge_buckets.begin();
       b_itr != edge_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( unsigned edge_ordinal = 0; edge_ordinal< b.size(); ++edge_ordinal) {
      STKUNIT_EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateEdges , testSkinAndCreateEdges3x3x3 )
{
  const size_t elem_rank = MetaData::ELEMENT_RANK;
  const size_t face_rank = MetaData::FACE_RANK;
  const size_t edge_rank = MetaData::EDGE_RANK;
  const size_t node_rank = MetaData::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[node_rank] , 64u ); // nodes
    STKUNIT_EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[face_rank] , 0u );  // faces
    STKUNIT_EXPECT_EQ( counts[elem_rank] , 27u ); // elements
  }

  stk::mesh::skin_mesh(fixture.m_bulk_data, MetaData::ELEMENT_RANK, NULL);
  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( 64u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 144u, counts[edge_rank] );  // edges
    STKUNIT_EXPECT_EQ( 54u, counts[face_rank] );  // faces
    STKUNIT_EXPECT_EQ( 27u, counts[elem_rank]  ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned elem_ordinal = i;
      STKUNIT_EXPECT_EQ( 12u, b.num_edges(elem_ordinal) );
      STKUNIT_EXPECT_EQ( 8u,  b.num_nodes(elem_ordinal) );
    }
  }

  stk::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = face_buckets.begin();
       b_itr != face_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned face_ordinal = i;
      STKUNIT_EXPECT_EQ( 4u, b.num_edges(face_ordinal) );
      STKUNIT_EXPECT_EQ( 4u, b.num_nodes(face_ordinal) );
    }
  }

  stk::mesh::BucketVector  edge_buckets = fixture.m_bulk_data.buckets(edge_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = edge_buckets.begin();
       b_itr != edge_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned edge_ordinal = i;
      STKUNIT_EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateEdges , testCreateEdges3x3 )
{
  const size_t elem_rank = MetaData::ELEMENT_RANK;
  const size_t edge_rank = MetaData::EDGE_RANK;
  const size_t node_rank = MetaData::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;

  stk::mesh::fixtures::QuadFixture fixture(MPI_COMM_WORLD, NX, NY);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[node_rank] , 16u ); // nodes
    STKUNIT_EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[elem_rank] , 9u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( 16u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 24u, counts[edge_rank] );  // edges
    STKUNIT_EXPECT_EQ( 9u, counts[elem_rank]  ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned elem_ordinal = i;
      STKUNIT_EXPECT_EQ( 4u, b.num_edges(elem_ordinal) );
      STKUNIT_EXPECT_EQ( 4u, b.num_nodes(elem_ordinal) );
    }
  }

  stk::mesh::BucketVector  edge_buckets = fixture.m_bulk_data.buckets(edge_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = edge_buckets.begin();
       b_itr != edge_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned edge_ordinal = i;
      STKUNIT_EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}
