#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateEdges.hpp>  // for create_edges
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc

using stk::mesh::MetaData;

TEST ( UnitTestCreateEdges, Quad_2x1 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 2, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 6u ); // nodes
    EXPECT_EQ( counts[1] , 0u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 2u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 6u ); // nodes
    EXPECT_EQ( counts[1] , 7u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 2u ); // elements
  }
}

TEST ( UnitTestCreateEdges, Quad_3x1 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 3, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 8u ); // nodes
    EXPECT_EQ( counts[1] , 0u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 3u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 8u ); // nodes
    EXPECT_EQ( counts[1] , 10u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 3u ); // elements
  }
}

TEST ( UnitTestCreateEdges, Quad_2x2 )
{
  stk::mesh::fixtures::QuadFixture fixture( MPI_COMM_WORLD, 2, 2);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  stk::mesh::skin_mesh(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 9u ); // nodes
    EXPECT_EQ( counts[1] , 8u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 4u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 9u ); // nodes
    EXPECT_EQ( counts[1] , 12u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 4u ); // elements
  }

  //shouldn't do anything
  stk::mesh::create_edges(fixture.m_bulk_data);
}

TEST ( UnitTestCreateEdges, Hex_2x1x1 )
{
  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 1, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 12u ); // nodes
    EXPECT_EQ( counts[1] , 0u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 2u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 12u ); // nodes
    EXPECT_EQ( counts[1] , 20u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 2u ); // elements
  }
}

TEST( UnitTestCreateEdges , Hex_3x1x1 )
{
  using namespace stk::mesh;

  fixtures::HexFixture fixture(MPI_COMM_WORLD, 3, 1, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 16u ); // nodes
    EXPECT_EQ( counts[1] , 0u );  // edges
    EXPECT_EQ( counts[2] , 0u );  // faces
    EXPECT_EQ( counts[3] , 3u ); // elements
  }

  create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 16u );
    EXPECT_EQ( counts[1] , 28u );
    EXPECT_EQ( counts[2] , 0u );
    EXPECT_EQ( counts[3] , 3u );
  }
}

TEST( UnitTestCreateEdges , testCreateEdges3x3x3 )
{
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
  const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[node_rank] , 64u ); // nodes
    EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    EXPECT_EQ( counts[face_rank] , 0u );  // faces
    EXPECT_EQ( counts[elem_rank] , 27u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( 64u, counts[node_rank] ); // nodes
    EXPECT_EQ( 144u, counts[edge_rank] );  // edges
    EXPECT_EQ( 0u, counts[face_rank] );  // faces
    EXPECT_EQ( 27u, counts[elem_rank]  ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      EXPECT_EQ( 0u, b.num_faces(i) );
      EXPECT_EQ( 12u, b.num_edges(i) );
      EXPECT_EQ( 8u,  b.num_nodes(i) );
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
      EXPECT_EQ( 4u, b.num_edges(i) );
      EXPECT_EQ( 4u, b.num_nodes(i) );
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
      EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}

TEST( UnitTestCreateEdges , testSkinAndCreateEdges3x3x3 )
{
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
  const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[node_rank] , 64u ); // nodes
    EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    EXPECT_EQ( counts[face_rank] , 0u );  // faces
    EXPECT_EQ( counts[elem_rank] , 27u ); // elements
  }

  stk::mesh::skin_mesh(fixture.m_bulk_data);
  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( 64u, counts[node_rank] ); // nodes
    EXPECT_EQ( 144u, counts[edge_rank] );  // edges
    EXPECT_EQ( 54u, counts[face_rank] );  // faces
    EXPECT_EQ( 27u, counts[elem_rank]  ); // elements
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
      EXPECT_EQ( 12u, b.num_edges(elem_ordinal) );
      EXPECT_EQ( 8u,  b.num_nodes(elem_ordinal) );
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
      EXPECT_EQ( 4u, b.num_edges(face_ordinal) );
      EXPECT_EQ( 4u, b.num_nodes(face_ordinal) );
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
      EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}

TEST( UnitTestCreateEdges , testCreateEdges3x3 )
{
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;

  stk::mesh::fixtures::QuadFixture fixture(MPI_COMM_WORLD, NX, NY);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[node_rank] , 16u ); // nodes
    EXPECT_EQ( counts[edge_rank] , 0u );  // edges
    EXPECT_EQ( counts[elem_rank] , 9u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( 16u, counts[node_rank] ); // nodes
    EXPECT_EQ( 24u, counts[edge_rank] );  // edges
    EXPECT_EQ( 9u, counts[elem_rank]  ); // elements
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
      EXPECT_EQ( 4u, b.num_edges(elem_ordinal) );
      EXPECT_EQ( 4u, b.num_nodes(elem_ordinal) );
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
      EXPECT_EQ( 2u, b.num_nodes(edge_ordinal) );
    }
  }
}
