
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>

#include <iomanip>
#include <algorithm>
#include <sstream>

STKUNIT_UNIT_TEST( UnitTestStkMeshSkinning , testCreateAdjacentEntities3x1x1 )
{
  const size_t NX = 3;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk_classic::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_fem_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 16u ); // nodes
    STKUNIT_EXPECT_EQ( counts[1] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[2] , 0u );  // faces
    STKUNIT_EXPECT_EQ( counts[3] , 3u ); // elements
  }

  stk_classic::mesh::PartVector empty_add_parts;

  stk_classic::mesh::create_adjacent_entities(fixture.m_bulk_data, empty_add_parts);

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[0] , 16u );
    STKUNIT_EXPECT_EQ( counts[1] , 28u );
    STKUNIT_EXPECT_EQ( counts[2] , 16u );
    STKUNIT_EXPECT_EQ( counts[3] , 3u );
  }
}

STKUNIT_UNIT_TEST( UnitTestStkMeshSkinning , testCreateAdjacentEntities3x3x3 )
{
  const size_t elem_rank = 3;
  const size_t face_rank = 2;
  const size_t edge_rank = 1;
  const size_t node_rank = 0;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk_classic::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_fem_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[node_rank] , 64u ); // nodes
    STKUNIT_EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[face_rank] , 0u );  // faces
    STKUNIT_EXPECT_EQ( counts[elem_rank] , 27u ); // elements
  }

  stk_classic::mesh::PartVector empty_add_parts;

  stk_classic::mesh::create_adjacent_entities(fixture.m_bulk_data, empty_add_parts);

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( 64u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 144u, counts[edge_rank] );  // edges
    STKUNIT_EXPECT_EQ( 108u, counts[face_rank] );  // faces
    STKUNIT_EXPECT_EQ( 27u, counts[elem_rank]  ); // elements
  }

  stk_classic::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk_classic::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk_classic::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      stk_classic::mesh::Entity & elem = b[i];

      STKUNIT_EXPECT_EQ( 6u, elem.relations(face_rank).size() );
      STKUNIT_EXPECT_EQ( 12u, elem.relations(edge_rank).size() );
      STKUNIT_EXPECT_EQ( 8u,  elem.relations(node_rank).size() );

    }
  }

  stk_classic::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk_classic::mesh::BucketVector::iterator b_itr = face_buckets.begin();
       b_itr != face_buckets.end();
       ++b_itr
      )
  {
    stk_classic::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      stk_classic::mesh::Entity & face = b[i];
      STKUNIT_EXPECT_EQ( 4u,face.relations(edge_rank).size());
      STKUNIT_EXPECT_EQ( 4u, face.relations(node_rank).size() );
    }
  }

  stk_classic::mesh::BucketVector  edge_buckets = fixture.m_bulk_data.buckets(edge_rank);
  for ( stk_classic::mesh::BucketVector::iterator b_itr = edge_buckets.begin();
       b_itr != edge_buckets.end();
       ++b_itr
      )
  {
    stk_classic::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      stk_classic::mesh::Entity & edge = b[i];
      STKUNIT_EXPECT_EQ( 2u, edge.relations(node_rank).size() );
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestStkMeshSkinning , testCreateAdjacentEntities3x3 )
{
  const size_t elem_rank = 2;
  const size_t edge_rank = 1;
  const size_t node_rank = 0;

  const size_t NX = 3;
  const size_t NY = 3;

  stk_classic::mesh::fixtures::QuadFixture fixture(MPI_COMM_WORLD, NX, NY);

  fixture.m_fem_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( counts[node_rank] , 16u ); // nodes
    STKUNIT_EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    STKUNIT_EXPECT_EQ( counts[elem_rank] , 9u ); // elements
  }

  stk_classic::mesh::PartVector empty_add_parts;

  stk_classic::mesh::create_adjacent_entities(fixture.m_bulk_data, empty_add_parts);

  {
    std::vector<size_t> counts ;
    stk_classic::mesh::fem::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( 16u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 24u, counts[edge_rank] );  // edges
    STKUNIT_EXPECT_EQ( 9u, counts[elem_rank]  ); // elements
  }

  stk_classic::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk_classic::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk_classic::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      stk_classic::mesh::Entity & elem = b[i];

      STKUNIT_EXPECT_EQ( 4u, elem.relations(edge_rank).size() );
      STKUNIT_EXPECT_EQ( 4u,  elem.relations(node_rank).size() );

    }
  }

  stk_classic::mesh::BucketVector  edge_buckets = fixture.m_bulk_data.buckets(edge_rank);
  for ( stk_classic::mesh::BucketVector::iterator b_itr = edge_buckets.begin();
       b_itr != edge_buckets.end();
       ++b_itr
      )
  {
    stk_classic::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      stk_classic::mesh::Entity & edge = b[i];
      STKUNIT_EXPECT_EQ( 2u, edge.relations(node_rank).size() );
    }
  }
}
