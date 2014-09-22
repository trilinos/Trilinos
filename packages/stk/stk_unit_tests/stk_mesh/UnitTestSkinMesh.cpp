// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_sum
#include <gtest/gtest.h>
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_mesh/base/CreateEdges.hpp"

TEST( SkinMesh, SimpleHex)
{
  const unsigned X = 5, Y = 5, Z = 5;
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, X, Y, Z);

  stk::mesh::EntityRank side_rank = fixture.m_meta.side_rank();

  stk::mesh::Part & skin_part = fixture.m_meta.declare_part("SkinPart", side_rank);
  stk::mesh::Part & skin_part_2 = fixture.m_meta.declare_part("SkinPart_2", side_rank);
  stk::mesh::Part & locally_owned = fixture.m_meta.locally_owned_part();

  fixture.m_meta.commit();

  fixture.generate_mesh();

  stk::mesh::BulkData & mesh = fixture.m_bulk_data;

  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

  stk::mesh::create_edges(mesh, fixture.m_meta.universal_part());

  std::cout<<"created "<<stk::mesh::count_selected_entities(fixture.m_meta.universal_part(), mesh.buckets(stk::topology::EDGE_RANK)) << " edges."<<std::endl;

  // skin the mesh
  {
    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(mesh, add_parts);
  }

  {
    size_t local_counts[2] = {}, global_counts[2] = {};
    local_counts[0] = stk::mesh::count_selected_entities( skin_part & locally_owned, mesh.buckets(stk::topology::NODE_RANK));
    local_counts[1] = stk::mesh::count_selected_entities( skin_part & locally_owned, mesh.buckets(side_rank));

    stk::all_reduce_sum( mesh.parallel(), local_counts, global_counts, 2);

    EXPECT_EQ( 152u, global_counts[0] );
    EXPECT_EQ( 150u, global_counts[1] );
  }

  // trying skinning the mesh again but put skin into part 2
  // skin the mesh
  {
    stk::mesh::PartVector add_parts(1,&skin_part_2);
    stk::mesh::skin_mesh(mesh, add_parts);
  }

  {
    size_t local_counts[2] = {}, global_counts[2] = {};
    local_counts[0] = stk::mesh::count_selected_entities( skin_part_2 & locally_owned, mesh.buckets(stk::topology::NODE_RANK));
    local_counts[1] = stk::mesh::count_selected_entities( skin_part_2 & locally_owned, mesh.buckets(side_rank));

    stk::all_reduce_sum( mesh.parallel(), local_counts, global_counts, 2);

    EXPECT_EQ( 152u, global_counts[0] );
    EXPECT_EQ( 150u, global_counts[1] );
  }

}

TEST( SkinMesh, SimpleQuad)
{
  const unsigned X = 5, Y = 5;
  stk::mesh::fixtures::QuadFixture fixture(MPI_COMM_WORLD, X, Y);

  stk::mesh::EntityRank side_rank = fixture.m_meta.side_rank();

  stk::mesh::Part & skin_part = fixture.m_meta.declare_part("SkinPart", side_rank);
  stk::mesh::Part & skin_part_2 = fixture.m_meta.declare_part("SkinPart_2", side_rank);
  stk::mesh::Part & locally_owned = fixture.m_meta.locally_owned_part();

  fixture.m_meta.commit();

  fixture.generate_mesh();

  stk::mesh::BulkData & mesh = fixture.m_bulk_data;

  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

  // skin the mesh
  {
    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(mesh, add_parts);
  }

  {
    size_t local_counts[2] = {}, global_counts[2] = {};
    local_counts[0] = stk::mesh::count_selected_entities( skin_part & locally_owned, mesh.buckets(stk::topology::NODE_RANK));
    local_counts[1] = stk::mesh::count_selected_entities( skin_part & locally_owned, mesh.buckets(side_rank));

    stk::all_reduce_sum( mesh.parallel(), local_counts, global_counts, 2);

    EXPECT_EQ( 20u, global_counts[0] );
    EXPECT_EQ( 20u, global_counts[1] );
  }

  // trying skinning the mesh again but put skin into part 2
  // skin the mesh
  {
    stk::mesh::PartVector add_parts(1,&skin_part_2);
    stk::mesh::skin_mesh(mesh, add_parts);
  }

  {
    size_t local_counts[2] = {}, global_counts[2] = {};
    local_counts[0] = stk::mesh::count_selected_entities( skin_part_2 & locally_owned, mesh.buckets(stk::topology::NODE_RANK));
    local_counts[1] = stk::mesh::count_selected_entities( skin_part_2 & locally_owned, mesh.buckets(side_rank));

    stk::all_reduce_sum( mesh.parallel(), local_counts, global_counts, 2);

    EXPECT_EQ( 20u, global_counts[0] );
    EXPECT_EQ( 20u, global_counts[1] );
  }

}
