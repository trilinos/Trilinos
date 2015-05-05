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
#include "stk_unit_test_utils/ioUtils.hpp"

void test_skin_mesh_with_hexes(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  const int spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::EntityRank side_rank = meta.side_rank();

  stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
  stk::mesh::Part & skin_part_2 = meta.declare_part("SkinPart_2", side_rank);
  stk::mesh::Part & locally_owned = meta.locally_owned_part();

  stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD, autoAuraOption);
  stk::unit_test_util::fill_mesh_using_stk_io("generated:5x5x5", mesh, MPI_COMM_WORLD);

  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

  stk::mesh::create_edges(mesh, meta.universal_part());

  std::cout<<"created "<<stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::EDGE_RANK)) << " edges."<<std::endl;

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

TEST( SkinMesh, SimpleHexWithAura)
{
    test_skin_mesh_with_hexes(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, SimpleHexWithoutAura)
{
    test_skin_mesh_with_hexes(stk::mesh::BulkData::NO_AUTO_AURA);
}

void move_element2_into_part(stk::mesh::BulkData& mesh, stk::mesh::EntityId element_id, stk::mesh::Part& part)
{
    stk::mesh::EntityVector entities;
    std::vector<stk::mesh::PartVector> add_parts_per_entity;
    std::vector<stk::mesh::PartVector> remove_parts_per_entity;

    stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEM_RANK, element_id);
    if ( mesh.is_valid(element) && mesh.bucket(element).owned() )
    {
        entities.push_back(element);
        stk::mesh::PartVector add_parts;
        stk::mesh::PartVector rm_parts;
        add_parts.push_back(&part);
        add_parts_per_entity.push_back(add_parts);
        remove_parts_per_entity.push_back(rm_parts);
    }

    mesh.batch_change_entity_parts(entities, add_parts_per_entity, remove_parts_per_entity);
}

void test_2_hex_2_block(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);

        stk::mesh::EntityRank side_rank = meta.side_rank();

        stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
        stk::mesh::Part & block_2 = meta.declare_part("block_2", stk::topology::FACE_RANK);

        stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD, autoAuraOption);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x2", mesh, MPI_COMM_WORLD);

        ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
        ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

        stk::mesh::EntityId element_id = 2;
        move_element2_into_part(mesh, element_id, block_2);

        stk::mesh::PartVector skin_parts;
        skin_parts.push_back(&skin_part);

        stk::mesh::skin_mesh(mesh, block_2, skin_parts);

        stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEM_RANK, element_id);

        if (mesh.is_valid(element2))
        {
            unsigned num_faces = mesh.num_faces(element2);
            EXPECT_EQ(5u, num_faces);
        }

        stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);

        // with correct face connection behavior, shouldn't this be 1 for num_faces?
        if (mesh.is_valid(element1))
        {
            unsigned num_faces = mesh.num_faces(element1);
            EXPECT_EQ(0u, num_faces);
        }
    }
}

TEST( SkinMesh, test_2_hex_2_block_with_aura)
{
    test_2_hex_2_block(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, test_2_hex_2_block_without_aura)
{
    test_2_hex_2_block(stk::mesh::BulkData::NO_AUTO_AURA);
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
