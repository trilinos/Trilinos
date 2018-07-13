// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_unit_tests/stk_mesh_fixtures/HexFixture.hpp>  // for HexFixture
#include <vector>                       // for vector, vector<>::iterator
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc

using stk::mesh::MetaData;

namespace {
  const size_t nodes_per_hex = 8;
  const size_t faces_per_hex = 6;
  const size_t nodes_per_quad= 4;

  size_t exp_hex_face_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_face = 3 * nx * ny * nz;
    exp_face += ny * nz + nz * nx + nx * ny;
    return exp_face;
  }

  size_t exp_node_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_node = (nx+1) * (ny+1) * (nz+1);
    return exp_node;
  }

  size_t exp_hex_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_elem = nx * ny * nz;
    return exp_elem;
  }

  TEST( UnitTestCreateFaces , testSkinAndCreateFaces3x3x3 )
  {
    const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
    const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
    const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
    const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

    const size_t NX = 3;
    const size_t NY = 3;
    const size_t NZ = 3;

    stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

    fixture.m_meta.commit();
    fixture.generate_mesh();

    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
      EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
      EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
      EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
    }

    stk::mesh::skin_mesh(fixture.m_bulk_data);
    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( 64u, counts[node_rank] ); // nodes
      EXPECT_EQ( 0u,  counts[edge_rank] );  // edges
      EXPECT_EQ( 54u, counts[face_rank] );  // faces
      EXPECT_EQ( 27u, counts[elem_rank] ); // elements
    }

    stk::mesh::create_faces(fixture.m_bulk_data);

    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
      EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
      EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
      EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
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
        EXPECT_EQ( faces_per_hex,  b.num_faces(elem_ordinal) );
        EXPECT_EQ( nodes_per_hex,  b.num_nodes(elem_ordinal) );
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
        EXPECT_EQ( 0u, b.num_edges(face_ordinal) );
        EXPECT_EQ( nodes_per_quad, b.num_nodes(face_ordinal) );
      }
    }

  }


  TEST( UnitTestCreateFaces , testCreateFacesThenSkin3x3x3 )
  {
    const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
    const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
    const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
    const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

    const size_t NX = 3;
    const size_t NY = 3;
    const size_t NZ = 3;

    stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

    fixture.m_meta.commit();
    fixture.generate_mesh();

    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
      EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
      EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
      EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
    }

    stk::mesh::create_faces(fixture.m_bulk_data);

    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
      EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
      EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
      EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
    }

    stk::mesh::skin_mesh(fixture.m_bulk_data);

    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

      EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
      EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
      EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
      EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
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
        EXPECT_EQ( faces_per_hex,  b.num_faces(elem_ordinal) );
        EXPECT_EQ( nodes_per_hex,  b.num_nodes(elem_ordinal) );
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
        EXPECT_EQ( 0u, b.num_edges(face_ordinal) );
        EXPECT_EQ( nodes_per_quad, b.num_nodes(face_ordinal) );
      }
    }

  }

} //end empty namespace

