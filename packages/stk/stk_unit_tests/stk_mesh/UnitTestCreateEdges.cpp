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

#include "Setup2Block2HexMesh.hpp"  // for setup2Block2HexMesh
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FEMHelpers.hpp"  // for declare_element
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator!, Selector
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateEdges.hpp>  // for create_edges
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/StkMeshFromGeneratedMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <string>                       // for string
#include <vector>                       // for vector, vector<>::iterator

using stk::mesh::MetaData;
using stk::unit_test_util::build_mesh;

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

TEST ( UnitTestCreateEdges, Quad9_2x1 )
{
  stk::mesh::fixtures::Quad9Fixture fixture( MPI_COMM_WORLD, 2, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] ,15u ); // nodes
    EXPECT_EQ( counts[1] , 0u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 2u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] ,15u ); // nodes
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


TEST ( UnitTestCreateEdges, Quad9_3x1 )
{
  stk::mesh::fixtures::Quad9Fixture fixture( MPI_COMM_WORLD, 3, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] ,21u ); // nodes
    EXPECT_EQ( counts[1] , 0u ); // edges
    EXPECT_EQ( counts[2] , 0u ); // faces
    EXPECT_EQ( counts[3] , 3u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 21u ); // nodes
    EXPECT_EQ( counts[1] , 10u ); // edges
    EXPECT_EQ( counts[2] ,  0u ); // faces
    EXPECT_EQ( counts[3] ,  3u ); // elements
  }
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
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, 3, 1, 1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( counts[0] , 16u ); // nodes
    EXPECT_EQ( counts[1] , 0u );  // edges
    EXPECT_EQ( counts[2] , 0u );  // faces
    EXPECT_EQ( counts[3] , 3u ); // elements
  }

  stk::mesh::create_edges(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

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

TEST( UnitTestCreateEdges , TwoBlockTwoHexTwoProc )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;

  setup2Block2HexMesh(bulk);

  stk::mesh::create_edges(bulk, *meta.get_part("block_1"));

  unsigned num_elems = stk::mesh::count_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::ELEM_RANK));
  unsigned num_edges = stk::mesh::count_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::EDGE_RANK));
  unsigned expected_num_elems = 2;//1 owned, 1 ghost on each proc
  unsigned expected_num_edges = 12;//edges only on the block_1 elem
  EXPECT_EQ(expected_num_elems, num_elems);
  EXPECT_EQ(expected_num_edges, num_edges);
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

TEST( UnitTestCreateEdges , hex1x1x4 )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:1x1x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();

    const stk::mesh::BucketVector &buckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);

    EXPECT_EQ(0u, buckets.size());

    stk::mesh::create_edges(stkMeshBulkData);

    EXPECT_NE(0u, buckets.size());

    stkMeshBulkData.modification_begin();
    stkMeshBulkData.modification_end();
  }
}

TEST( UnitTestCreateEdges, hybrid_HexPyrTet )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.1
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.1  |
  //       |   |          |   |          |   | <-- (Undrawable transition pyramid between node (5,6,7,8,9)
  //       |   |   1.0    |   |          |   |      and 4 tets contained in volume on the right)
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size != 2)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;
  const int p_rank = mesh.parallel_rank();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
  stk::mesh::Part * pyrPart = &meta.declare_part_with_topology("pyr_part", stk::topology::PYRAMID_5);
  stk::mesh::Part * tetPart = &meta.declare_part_with_topology("tet_part", stk::topology::TET_4);
  meta.commit();

  const size_t numHex = 1;
  stk::mesh::EntityIdVector hexNodeIDs[] {
    { 1, 2, 3, 4, 5, 6, 7, 8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  const size_t numPyr = 1;
  stk::mesh::EntityIdVector pyrNodeIDs[] {
    { 5, 6, 7, 8, 9 }
  };
  stk::mesh::EntityId pyrElemIDs[] = { 2 };

  const size_t numTet = 4;
  stk::mesh::EntityIdVector tetNodeIDs[] {
    { 7, 8, 9, 12 },
    { 6, 9, 10, 7 },
    { 7, 9, 10, 12 },
    { 7, 12, 10, 11 }
  };
  stk::mesh::EntityId tetElemIDs[] = { 3, 4, 5, 6 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  int shared_nodeIDs_and_procs[][3] =
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };
  int numSharedNodeTriples = 8;

  mesh.modification_begin();

  if (p_rank == 0) {
    for (size_t i = 0; i < numHex; ++i) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  else {
    for (size_t i = 0; i < numPyr; ++i) {
      stk::mesh::declare_element(mesh, *pyrPart, pyrElemIDs[i], pyrNodeIDs[i]);
    }
    for (size_t i = 0; i < numTet; ++i) {
      stk::mesh::declare_element(mesh, *tetPart, tetElemIDs[i], tetNodeIDs[i]);
    }
  }

  for (int nodeIdx = 0; nodeIdx < numSharedNodeTriples; ++nodeIdx) {
    if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
      stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
      int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
      mesh.add_node_sharing(node, sharingProc);
    }
  }

  mesh.modification_end();

  stk::mesh::create_edges(mesh, *meta.get_part("pyr_part"));

  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(mesh, counts);

    EXPECT_EQ( counts[stk::topology::NODE_RANK], 12u ); // nodes
    EXPECT_EQ( counts[stk::topology::EDGE_RANK], 8u );  // edges
    EXPECT_EQ( counts[stk::topology::ELEM_RANK], 6u );  // elements
  }
}

TEST ( UnitTestCreateEdges, Hex_2x1x1_select_out_a_face )
{
  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 1, 1);

  stk::mesh::Part & facePart = fixture.m_meta.declare_part_with_topology("face_part_to_exclude", stk::topology::QUADRILATERAL_4, true);
  stk::mesh::PartVector facePartVector(1, &facePart);
  fixture.m_meta.commit();
  fixture.generate_mesh();
  stk::mesh::BulkData & mesh = fixture.m_bulk_data;
  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(mesh , counts);

    EXPECT_EQ( counts[stk::topology::NODE_RANK] , 12u );
    EXPECT_EQ( counts[stk::topology::EDGE_RANK] , 0u );
    EXPECT_EQ( counts[stk::topology::FACE_RANK] , 0u );
    EXPECT_EQ( counts[stk::topology::ELEM_RANK] , 2u );
  }
  stk::mesh::create_faces(mesh);
  const stk::mesh::BucketVector & buckets = mesh.buckets(stk::topology::FACE_RANK);
  mesh.modification_begin();
  for (unsigned bucketCount = 0 ; bucketCount < buckets.size() ; ++bucketCount) {
    stk::mesh::Bucket & bucket = *buckets[bucketCount];
    for (unsigned count = 0; count < bucket.size(); ++count) {
      stk::mesh::Entity face = bucket[count];
      if (mesh.bucket(face).owned()) {
        mesh.change_entity_parts(face, facePartVector);
      }
    }
  }
  mesh.modification_end();

  stk::mesh::Selector selectOutFaces(!facePart);

  stk::mesh::create_edges(mesh, selectOutFaces);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( mesh, counts);

    EXPECT_EQ( counts[stk::topology::NODE_RANK] , 12u );
    EXPECT_EQ( counts[stk::topology::EDGE_RANK] , 20u );
    EXPECT_EQ( counts[stk::topology::FACE_RANK] , 11u );
    EXPECT_EQ( counts[stk::topology::ELEM_RANK] , 2u );
  }
}

TEST( UnitTestCreateEdges, particle )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if (p_size != 1) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * particlePart = &meta.declare_part_with_topology("particle_part", stk::topology::PARTICLE);
  meta.commit();

  stk::mesh::EntityIdVector particleNodeIDs { 1 };
  stk::mesh::EntityId particleElemID( 1 );

  mesh.modification_begin();
  stk::mesh::declare_element(mesh, *particlePart, particleElemID, particleNodeIDs);
  mesh.modification_end();

  stk::mesh::create_edges(mesh);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  EXPECT_EQ( 1u, counts[stk::topology::NODE_RANK] );
  EXPECT_EQ( 0u, counts[stk::topology::EDGE_RANK] );  // Shouldn't create any; particles don't have edges
  EXPECT_EQ( 0u, counts[stk::topology::FACE_RANK] );
  EXPECT_EQ( 1u, counts[stk::topology::ELEM_RANK] );
}

