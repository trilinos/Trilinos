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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size
#include <vector>                       // for vector, vector<>::iterator
#include "mpi.h"                        // for MPI_COMM_WORLD

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/GearsFixture.hpp"  // for GearsFixture, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/TetFixture.hpp"  // for TetFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/degenerate_mesh.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"

using stk::mesh::MetaData;
using stk::unit_test_util::build_mesh;

namespace {
const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

const size_t nodes_per_hex = 8;
const size_t faces_per_hex = 6;
const size_t nodes_per_quad= 4;

const size_t nodes_per_tet = 4;
const size_t faces_per_tet = 4;
const size_t nodes_per_tri = 3;

size_t exp_hex_face_count(size_t nx, size_t ny, size_t nz)
{
  size_t exp_face = 3 * nx * ny * nz;
  exp_face += ny * nz + nz * nx + nx * ny;
  return exp_face;
}

size_t exp_tet_face_count(size_t nx, size_t ny, size_t nz)
{
  size_t exp_face = 12 * nx * ny * nz;
  exp_face += 2* (ny * nz + nz * nx + nx * ny);
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

size_t exp_tet_count(size_t nx, size_t ny, size_t nz)
{
  size_t exp_elem = 6 * nx * ny * nz;
  return exp_elem;
}


TEST ( UnitTestCreateFaces, Hex_2x1x1 )
{
  const size_t NX = 2;
  const size_t NY = 1;
  const size_t NZ = 1;

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

}

TEST ( UnitTestCreateFaces, Tet_2x1x1 )
{
  const size_t NX = 2;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

}

TEST( UnitTestCreateFaces , Hex_3x1x1 )
{
  const size_t NX = 3;
  const size_t NY = 1;
  const size_t NZ = 1;

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
}

TEST( UnitTestCreateFaces , Tet_3x1x1 )
{
  const size_t NX = 3;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }
}

TEST( UnitTestCreateFaces , testCreateFaces3x3x3 )
{
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

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
        b_itr != elem_buckets.end();
        ++b_itr
        )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      EXPECT_EQ( faces_per_hex, b.num_faces(i) );
      EXPECT_EQ( 0u, b.num_edges(i) );
      EXPECT_EQ( nodes_per_hex,  b.num_nodes(i) );
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
      EXPECT_EQ( 0u, b.num_edges(i) );
      EXPECT_EQ( nodes_per_quad, b.num_nodes(i) );
    }
  }

  unsigned num_ghosted_faces = 0;

  stk::mesh::Selector ghosted_part = fixture.m_meta.universal_part() &
      !(fixture.m_meta.locally_owned_part() | fixture.m_meta.globally_shared_part());

  stk::mesh::BucketVector  ghosted_elem_buckets = fixture.m_bulk_data.get_buckets(elem_rank,ghosted_part);

  for ( stk::mesh::BucketVector::iterator b_itr = ghosted_elem_buckets.begin();
        b_itr != ghosted_elem_buckets.end();
        ++b_itr
        )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      const stk::mesh::Entity & elem = b[i];

      const unsigned num_faces = fixture.m_bulk_data.num_faces(elem);
      EXPECT_EQ(faces_per_hex, num_faces);

      const stk::mesh::Entity * faces = fixture.m_bulk_data.begin_faces(elem);
      unsigned elem_num_ghosted_faces = 0;
      unsigned elem_num_shared_faces = 0;

      const int my_proc_id = fixture.m_bulk_data.parallel_rank();

      for (unsigned k = 0; k < num_faces; ++k)
      {
        const stk::mesh::Entity & face = faces[k];

        const bool is_ghosted = fixture.m_bulk_data.in_receive_ghost(fixture.m_bulk_data.entity_key(face));
        if (is_ghosted)
        {
          elem_num_ghosted_faces += 1;
          const int face_owner = fixture.m_bulk_data.parallel_owner_rank(face);
          EXPECT_NE(face_owner, my_proc_id);
        }

        const bool is_shared = fixture.m_bulk_data.in_shared(fixture.m_bulk_data.entity_key(face));
        if (is_shared) elem_num_shared_faces += 1;

        const bool is_ghosted_or_shared = is_ghosted || is_shared;
        ASSERT_TRUE(is_ghosted_or_shared);

      }

      EXPECT_LE(1u, elem_num_ghosted_faces);
      num_ghosted_faces += elem_num_ghosted_faces;

      const unsigned shared_and_ghosted = elem_num_shared_faces + elem_num_ghosted_faces;
      EXPECT_EQ(faces_per_hex, shared_and_ghosted);
    }
  }

  if (fixture.m_bulk_data.parallel_size() > 1)
  {
    EXPECT_LE(1u, num_ghosted_faces);
  }
}

TEST( UnitTestCreateFaces , testCreateTetFaces3x3x3 )
{
  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
        b_itr != elem_buckets.end();
        ++b_itr
        )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      EXPECT_EQ( faces_per_tet, b.num_faces(i) );
      EXPECT_EQ( 0u, b.num_edges(i) );
      EXPECT_EQ( nodes_per_tet,  b.num_nodes(i) );
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
      EXPECT_EQ( 0u, b.num_edges(i) );
      EXPECT_EQ( nodes_per_tri, b.num_nodes(i) );
    }
  }

  unsigned num_ghosted_faces = 0;

  stk::mesh::Selector ghosted_part = fixture.m_meta.universal_part() &
      !(fixture.m_meta.locally_owned_part() | fixture.m_meta.globally_shared_part());

  stk::mesh::BucketVector  ghosted_elem_buckets = fixture.m_bulk_data.get_buckets(elem_rank,ghosted_part);

  for ( stk::mesh::BucketVector::iterator b_itr = ghosted_elem_buckets.begin();
        b_itr != ghosted_elem_buckets.end();
        ++b_itr
        )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      const stk::mesh::Entity & elem = b[i];

      const unsigned num_faces = fixture.m_bulk_data.num_faces(elem);
      EXPECT_EQ(faces_per_tet, num_faces);

      const stk::mesh::Entity * faces = fixture.m_bulk_data.begin_faces(elem);
      unsigned elem_num_ghosted_faces = 0;
      unsigned elem_num_shared_faces = 0;

      const int my_proc_id = fixture.m_bulk_data.parallel_rank();

      for (unsigned k = 0; k < num_faces; ++k)
      {
        const stk::mesh::Entity & face = faces[k];

        const bool is_ghosted = fixture.m_bulk_data.in_receive_ghost(fixture.m_bulk_data.entity_key(face));
        if (is_ghosted)
        {
          elem_num_ghosted_faces += 1;
          const int face_owner = fixture.m_bulk_data.parallel_owner_rank(face);
          EXPECT_NE(face_owner, my_proc_id);
        }

        const bool is_shared = fixture.m_bulk_data.in_shared(fixture.m_bulk_data.entity_key(face));
        if (is_shared) elem_num_shared_faces += 1;

        const bool is_ghosted_or_shared = is_ghosted || is_shared;
        ASSERT_TRUE(is_ghosted_or_shared);

      }

      EXPECT_LE(1u, elem_num_ghosted_faces);
      num_ghosted_faces += elem_num_ghosted_faces;

      const unsigned shared_and_ghosted = elem_num_shared_faces + elem_num_ghosted_faces;
      EXPECT_EQ(faces_per_tet, shared_and_ghosted);
    }
  }

  if (fixture.m_bulk_data.parallel_size() > 1)
  {
    EXPECT_LE(1u, num_ghosted_faces);
  }
}

TEST ( UnitTestCreateFaces, Gears )
{
  stk::mesh::fixtures::GearsFixture fixture( MPI_COMM_WORLD, 1,
                                                            stk::mesh::fixtures::GearParams(0.1, 0.4, 1.0, -0.4, 0.4));

  fixture.meta_data.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.bulk_data , counts);

    EXPECT_EQ( 4340u,  counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    EXPECT_EQ( 3348u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.bulk_data , counts);

    EXPECT_EQ( 4340u, counts[node_rank] ); // nodes
    EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    EXPECT_EQ( 10974u, counts[face_rank] ); // faces
    EXPECT_EQ( 3348u, counts[elem_rank] ); // elements
  }

}

void heterogeneous_create_faces_test(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numprocs > 1)
    return;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD, autoAuraOption);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;
  stk::mesh::fixtures::VectorFieldType & node_coord =
      meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh( node_coord , meta_data.universal_part() , 3, nullptr);

  stk::mesh::fixtures::heterogeneous_mesh_meta_data( meta_data , node_coord );
  meta_data.commit();

  stk::mesh::fixtures::heterogeneous_mesh_bulk_data( bulk_data , node_coord );

  /*
   * Three hexes, three wedges, three tets, two pyramids,
   * three quad shells, and three triangle shells.
   *
   *  Z = 0 plane:
   *
   *    Y
   *    ^   9      10
   *    !   *-------*
   *    !  / \     / \
   *    ! /   \   /   \
   *     /     \ /     \
   *    *-------*-------*-------*
   *   5|      6|      7|      8|
   *    |       |       |       |
   *    |       |       |       |
   *    *-------*-------*-------*    ----> X
   *    1       2       3       4
   *
   *  Z = -1 plane:
   *
   *    Y
   *    ^  19      20
   *    !   *-------*
   *    !  / \     / \
   *    ! /   \   /   \
   *     /     \ /     \
   *    *-------*-------*-------*
   *  15|     16|     17|     18|
   *    |       |       |       |
   *    |       |       |       |
   *    *-------*-------*-------*    ----> X
   *   11      12      13      14
   *
   *
   *  Last node (#21) at Z = -2, translated from node #16
   */

  /*
   * Total elements = 21 = 3 + 3 + 3 + 2 + 3 + 3
   * Total faces    = 45 = 6 front + 6 back + 9 perimeter + 6 internal + 7 from tet + 5 from pyr + 6 from extra (2sided) shells
   */
  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    EXPECT_EQ( 21u, counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,  counts[edge_rank] ); // edges
    EXPECT_EQ( 0u,  counts[face_rank] ); // faces
    EXPECT_EQ( 17u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    EXPECT_EQ( 21u, counts[node_rank] ); // nodes
    EXPECT_EQ( 0u,  counts[edge_rank] ); // edges
    EXPECT_EQ( 45u, counts[face_rank] ); // faces
    EXPECT_EQ( 17u, counts[elem_rank] ); // elements
  }
}

TEST ( UnitTestCreateFaces, HeterogeneousWithAura )
{
  heterogeneous_create_faces_test(stk::mesh::BulkData::AUTO_AURA );
}

TEST ( UnitTestCreateFaces, HeterogeneousNoAura )
{
  heterogeneous_create_faces_test(stk::mesh::BulkData::NO_AUTO_AURA );
}

TEST ( UnitTestCreateFaces, Degenerate )
{
  int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numprocs > 1)
    return;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;
  stk::mesh::fixtures::VectorFieldType & node_coord =
      meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh( node_coord , meta_data.universal_part() , 3, nullptr);

  stk::mesh::fixtures::degenerate_mesh_meta_data( meta_data , node_coord );
  meta_data.commit();

  stk::mesh::fixtures::degenerate_mesh_bulk_data( bulk_data , node_coord );

  /*
   *  Z = 0 plane:
   *
   *    Y
   *    ^
   *    !
   *
   *   4*       *5
   *    |\     /|
   *    | \   / |
   *    |  \ /  |
   *    *---*---* ----> X
   *    1   2   3
   *
   *  Z = -1 plane:
   *
   *    Y
   *    ^
   *    !
   *
   *   9*       *10
   *    |\     /|
   *    | \   / |
   *    |  \ /  |
   *    *---*---* ----> X
   *    6   7   8
   *
   * Total elements = 2 hexes
   * Total faces    = 11 = 2 front + 2 back + 6 perimeter + 1 internal degenerate
   */
  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    EXPECT_EQ(10u, counts[node_rank] ); // nodes
    EXPECT_EQ( 0u, counts[edge_rank] ); // edges
    EXPECT_EQ( 0u, counts[face_rank] ); // faces
    EXPECT_EQ( 2u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    EXPECT_EQ(10u, counts[node_rank] ); // nodes
    EXPECT_EQ( 0u, counts[edge_rank] ); // edges
    EXPECT_EQ(11u, counts[face_rank] ); // faces
    EXPECT_EQ( 2u, counts[elem_rank] ); // elements
  }
}

} //end empty namespace

