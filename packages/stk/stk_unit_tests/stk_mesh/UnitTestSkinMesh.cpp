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
// //     * Redistributions in binary form must reproduce the above
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

#include <array>
#include <algorithm>
#include <stk_mesh/base/Comm.hpp>
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/TetFixture.hpp>  // for QuadFixture
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_sum
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

void put_all_faces_in_io_part(stk::mesh::BulkData &mesh, stk::mesh::Selector locally_owned, stk::mesh::Part& face_output_part)
{
    mesh.modification_begin();
    stk::mesh::PartVector part_vector;
    part_vector.push_back(&face_output_part);
    stk::mesh::EntityVector stk_faces;
    stk::mesh::get_entities(mesh, stk::topology::FACE_RANK, stk_faces);
    for(size_t count = 0; count < stk_faces.size(); ++count)
    {
        if(locally_owned(mesh.bucket(stk_faces[count])))
        {
            mesh.change_entity_parts(stk_faces[count], part_vector);
        }
    }
    mesh.modification_end();
}

void put_all_elements_in_io_part(stk::mesh::BulkData &mesh, stk::mesh::Selector locally_owned, stk::mesh::Part& element_output_part)
{
    mesh.modification_begin();
    stk::mesh::PartVector part_vector;
    part_vector.push_back(&element_output_part);
    stk::mesh::EntityVector stk_elements;
    stk::mesh::get_entities(mesh, stk::topology::ELEMENT_RANK, stk_elements);
    for(size_t count = 0; count < stk_elements.size(); ++count)
    {
        if(locally_owned(mesh.bucket(stk_elements[count])))
        {
            mesh.change_entity_parts(stk_elements[count], part_vector);
        }
    }
    mesh.modification_end();
}

void writeStkDebuggingFile(stk::io::StkMeshIoBroker &stkMeshIoBroker, stk::mesh::BulkData &mesh, const std::string &output_name)
{
    stk::mesh::Part & face_output_part = mesh.mesh_meta_data().declare_part_with_topology("output_face_name", stk::topology::TRI_3);
    stk::io::put_io_part_attribute(face_output_part);
    put_all_faces_in_io_part(mesh, mesh.mesh_meta_data().locally_owned_part(), face_output_part);

    size_t resultFileIndex = stkMeshIoBroker.create_output_mesh(output_name, stk::io::WRITE_RESULTS);
    stkMeshIoBroker.write_output_mesh(resultFileIndex);
}

template <std::size_t SIZE>
bool check_if_one_owned_face_with_these_nodes_exists(const std::array <uint64_t, SIZE> &nodes, stk::mesh::BulkData &mesh)
{
  stk::mesh::EntityVector face_vector;
  unsigned numNodes = nodes.size();
  stk::mesh::EntityVector node_vector(nodes.size());

  for (unsigned i = 0; i < numNodes; ++i)
  {
    node_vector[i] = mesh.get_entity(stk::topology::NODE_RANK, nodes[i]);
  }

  stk::mesh::impl::find_faces_these_nodes_have_in_common(mesh, numNodes, &node_vector[0], face_vector);

  // remove all non-locally owned entities
  face_vector.erase(
      std::remove_if(face_vector.begin(),
                     face_vector.end(),
                     [&](stk::mesh::Entity e){return !mesh.bucket(e).owned();}),
      face_vector.end());

  return face_vector.size() == 1;
}

void test_skin_mesh_with_hexes(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
     //  ID.proc
    //
    //          4.0------------8.0-----------12.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      3.0------------7.0-----------11.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |   |          |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);
  if(p_size > 2)
  {
    return;
  }

  const int spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::EntityRank side_rank = meta.side_rank();

  stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
  stk::mesh::Part & locally_owned = meta.locally_owned_part();

  stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD, autoAuraOption);
  // node ordering for generated mesh is different than hand-crafted meshes for other unit tests below
  stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x2", mesh, MPI_COMM_WORLD);
  const int p_rank = mesh.parallel_rank();

  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

  // skin the mesh
  stk::mesh::PartVector add_parts(1,&skin_part);
  stk::mesh::skin_mesh(mesh, add_parts);

  stk::mesh::Selector skin = skin_part & locally_owned;
  std::vector<size_t> counts(4);
  counts[0] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::NODE_RANK));
  counts[1] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::EDGE_RANK));
  counts[2] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::FACE_RANK));
  counts[3] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::ELEM_RANK));

//  stk::io::StkMeshIoBroker stkio(pm);
//  stkio.set_bulk_data(mesh);
//  writeStkDebuggingFile(stkio, mesh, "hex1x1x2.exo");

  if (p_rank == 0)
  {
    // check number of entities in skin part
    EXPECT_EQ( p_size == 2 ? 8u : 12u, counts[stk::topology::NODE_RANK] );  // nodes
    EXPECT_EQ(                     0u, counts[stk::topology::EDGE_RANK] );  // edges
    EXPECT_EQ( p_size == 2 ? 5u : 10u, counts[stk::topology::FACE_RANK] );  // face
    EXPECT_EQ(                     0u, counts[stk::topology::ELEM_RANK] );  // elements
    // check boundary faces are created
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 2, 4, 3 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 7, 3 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 3, 7, 8, 4 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 6, 2 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 2, 6, 8, 4 }}, mesh));

    if (p_size == 1)
    {
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 11, 7 }},   mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 7, 11, 12, 8 }},  mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 9, 10, 12, 11 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 10, 6 }},   mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 6, 10, 12, 8 }},  mesh));
    }
  }
  else if (p_rank == 1)
  {
    // check number of entities in skin part
    EXPECT_EQ( 4u, counts[stk::topology::NODE_RANK] );  // nodes
    EXPECT_EQ( 0u, counts[stk::topology::EDGE_RANK] );  // edges
    EXPECT_EQ( 5u, counts[stk::topology::FACE_RANK] );  // face
    EXPECT_EQ( 0u, counts[stk::topology::ELEM_RANK] );  // elements
    // check boundary faces are created
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 11, 7 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 7, 11, 12, 8 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 9, 10, 12, 11 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 10, 6 }}, mesh));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 6, 10, 12, 8 }}, mesh));
  }
}

TEST( SkinMesh, SkinHexWithAura )
{
  test_skin_mesh_with_hexes(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, SkinHexWithoutAura )
{
  test_skin_mesh_with_hexes(stk::mesh::BulkData::NO_AUTO_AURA);
}

void test_skin_mesh_with_tets(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);
  if(p_size > 2)
  {
    return;
  }

  const size_t NX = 1;
  const size_t NY = 1;
  const size_t NZ = 2;

  // fixture generates six tets from NXxNYxNZ hex
  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ, autoAuraOption);

  fixture.m_meta.commit();
  fixture.generate_mesh();
  const int p_rank = fixture.m_bulk_data.parallel_rank();

  stk::mesh::EntityRank side_rank = fixture.m_meta.side_rank();

  stk::mesh::Part & skin_part = fixture.m_meta.declare_part("SkinPart", side_rank);
  stk::mesh::Part & locally_owned = fixture.m_meta.locally_owned_part();

  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, fixture.m_bulk_data.buckets(stk::topology::NODE_RANK)) );
  ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, fixture.m_bulk_data.buckets(side_rank)) );

  // skin the mesh
  {
    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(fixture.m_bulk_data, add_parts);
  }
  stk::mesh::Selector skin = skin_part & locally_owned;
  std::vector<size_t> counts(4);
  counts[0] = stk::mesh::count_selected_entities( skin, fixture.m_bulk_data.buckets(stk::topology::NODE_RANK));
  counts[1] = stk::mesh::count_selected_entities( skin, fixture.m_bulk_data.buckets(stk::topology::EDGE_RANK));
  counts[2] = stk::mesh::count_selected_entities( skin, fixture.m_bulk_data.buckets(stk::topology::FACE_RANK));
  counts[3] = stk::mesh::count_selected_entities( skin, fixture.m_bulk_data.buckets(stk::topology::ELEM_RANK));

  if (0 == p_rank)
  {
    // check number of entities in skin part
    EXPECT_EQ( 1 == p_size ? 12u :  8u, counts[stk::topology::NODE_RANK] ); // nodes
    EXPECT_EQ(                      0u, counts[stk::topology::EDGE_RANK] );  // edges
    EXPECT_EQ( 1 == p_size ? 20u : 10u, counts[stk::topology::FACE_RANK] );  // face
    EXPECT_EQ(                      0u, counts[stk::topology::ELEM_RANK] );  // elements
    // check boundary faces are created
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 3, 4, 8 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 3, 7, 8 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 5, 7 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 3, 7 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 5, 6 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 2, 6 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 3, 4 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 2, 4 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 6, 8 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 4, 8 }}, fixture.m_bulk_data));

    if (1 == p_size) {
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 7, 12 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 7, 11 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 11, 12 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 11, 9 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 11, 12 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 10, 12 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 10 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 6, 10 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 12 }}, fixture.m_bulk_data));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 8, 12 }}, fixture.m_bulk_data));
    }
  }
  else if (1 == p_rank)
  {
    // check number of entities in skin part
    EXPECT_EQ( 4u, counts[stk::topology::NODE_RANK] ); // nodes
    EXPECT_EQ( 0u, counts[stk::topology::EDGE_RANK] );  // edges
    EXPECT_EQ( 10u, counts[stk::topology::FACE_RANK] );  // face
    EXPECT_EQ( 0u, counts[stk::topology::ELEM_RANK] );  // elements
    // check boundary faces are created
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 7, 12 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 7, 11 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 11, 12 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 11, 9 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 11, 12 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 10, 12 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 10 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 6, 10 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 12 }}, fixture.m_bulk_data));
    EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 8, 12 }}, fixture.m_bulk_data));
  }
}

TEST( SkinMesh, SkinTetWithAura )
{
  test_skin_mesh_with_tets(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, SkinTetWithoutAura )
{
  test_skin_mesh_with_tets(stk::mesh::BulkData::NO_AUTO_AURA);
}

void test_skin_mesh_with_wedge(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption, bool addShells)
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

    if(p_size > 2)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm, autoAuraOption);
    const int p_rank = mesh.parallel_rank();

    stk::mesh::EntityRank side_rank = meta.side_rank();

    stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
    stk::mesh::Part & locally_owned = meta.locally_owned_part();

    stk::mesh::Part * wedgePart = &meta.declare_part_with_topology("wedge_part", stk::topology::WEDGE_6);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUADRILATERAL_4);
    meta.commit();

    const size_t numWedges = 4;
    stk::mesh::EntityId elementIDsToProc[][2] =
    {
        { 1, 0 },  // proc 0
        { 2, 0 },  // proc 0
        { 3, 1 },  // proc 1
        { 4, 1 }   // proc 1
    };

    stk::mesh::EntityIdVector wedgeNodeIDs[] {
        { 1, 5,  2, 4,  8,  3 },
        { 2, 5,  6, 3,  8,  7 },
        { 5, 9,  6, 8, 12,  7 },
        { 6, 9, 10, 7, 12, 11 }
    };

    const size_t numShells = 1;
    stk::mesh::EntityId shellIDsToProc[][2] =
    {
        { 5, 0 }  // proc 0
    };

    stk::mesh::EntityIdVector shellNodeIDs[] {
        { 1, 2, 3, 4 }
    };

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

    for (size_t i = 0; i < numWedges; ++i) {
      if ( (1 == p_size) || (elementIDsToProc[i][1] == static_cast<unsigned>(p_rank)) )
      {
        stk::mesh::declare_element(mesh, *wedgePart, elementIDsToProc[i][0], wedgeNodeIDs[i]);
      }
    }

    if (addShells)
    {
      for (size_t i = 0; i < numShells; ++i) {
        if ( (1 == p_size) || (shellIDsToProc[i][1] == static_cast<unsigned>(p_rank)) )
        {
          stk::mesh::declare_element(mesh, *shellPart, shellIDsToProc[i][0], shellNodeIDs[i]);
        }
      }
    }

    if (p_size > 1)
    {
      for (int nodeIdx = 0; nodeIdx < numSharedNodeTriples; ++nodeIdx) {
          if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
              stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
              int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
              stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
              mesh.add_node_sharing(node, sharingProc);
          }
      }
    }

    mesh.modification_end();

    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(mesh, add_parts);

    stk::mesh::Selector skin = skin_part & locally_owned;

    std::vector<size_t> counts(4);
    counts[0] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::NODE_RANK));
    counts[1] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::EDGE_RANK));
    counts[2] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::FACE_RANK));
    counts[3] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::ELEM_RANK));

    if (p_rank == 0)
    {
      // check number of entities
      EXPECT_EQ( 1 == p_size ? 12u : 8u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ(                     0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 1 == p_size ? 14u : 7u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ(                     0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 8, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 2, 3, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 2, 6, 7, 3 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 2, 5 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 6, 2 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 4, 8, 3 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 7, 3 }}, mesh));

      if (1 == p_size) {
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 12, 8 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 9, 10, 11, 12 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 6, 10, 11, 7 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 12, 7 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 11, 7, 12 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 6 }}, mesh));
        EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 6, 9 }}, mesh));
      }
    }
    else if (p_rank == 1)
    {
      // check number of entities in skin part
      EXPECT_EQ( 4u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ( 0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 7u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ( 0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 12, 8 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 9, 10, 11, 12 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 6, 10, 11, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 12, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 11, 7, 12 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 6 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 6, 9 }}, mesh));
    }
}

TEST( SkinMesh, SkinWedgeWithAura )
{
  test_skin_mesh_with_wedge(stk::mesh::BulkData::AUTO_AURA, false);
}

TEST( SkinMesh, SkinWedgeWithoutAura )
{
  test_skin_mesh_with_wedge(stk::mesh::BulkData::NO_AUTO_AURA, false);
}

TEST( SkinMesh, SkinWedgeWithAuraWithShell )
{
  test_skin_mesh_with_wedge(stk::mesh::BulkData::AUTO_AURA, true);
}

TEST( SkinMesh, SkinWedgeWithoutAuraWithShell )
{
  test_skin_mesh_with_wedge(stk::mesh::BulkData::NO_AUTO_AURA, true);
}

void test_skin_mesh_with_pyramid(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
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

    if(p_size > 2)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm, autoAuraOption);
    const int p_rank = mesh.parallel_rank();

    stk::mesh::EntityRank side_rank = meta.side_rank();

    stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
    stk::mesh::Part & locally_owned = meta.locally_owned_part();

    stk::mesh::Part * pyramidPart = &meta.declare_part_with_topology("pyramid_part", stk::topology::PYRAMID_5);
    meta.commit();

    const size_t numPyramids = 6;
    stk::mesh::EntityId elementIDsToProc[][2] =
    {
        { 1, 0 },  // proc 0
        { 2, 0 },  // proc 0
        { 3, 0 },  // proc 0
        { 4, 1 },  // proc 1
        { 5, 1 },  // proc 1
        { 6, 1 }   // proc 1
    };

    stk::mesh::EntityIdVector pyramidNodeIDs[] {
        { 1, 4,  8,  5,  2 },
        { 5, 8,  7,  6,  2 },
        { 3, 7,  8,  4,  2 },
        { 5, 8, 12,  9, 10 },
        { 8, 7, 11, 12, 10 },
        { 5, 6,  7,  8, 10 }
    };

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

    for (size_t i = 0; i < numPyramids; ++i) {
      if ( (1 == p_size) || (elementIDsToProc[i][1] == static_cast<unsigned>(p_rank)) )
      {
        stk::mesh::Entity element = stk::mesh::declare_element(mesh, *pyramidPart, elementIDsToProc[i][0], pyramidNodeIDs[i]);
        ASSERT_TRUE(mesh.is_valid(element));
      }
    }

    if (p_size > 1)
    {
      for (int nodeIdx = 0; nodeIdx < numSharedNodeTriples; ++nodeIdx) {
          if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
              stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
              int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
              stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
              mesh.add_node_sharing(node, sharingProc);
          }
      }
    }

    mesh.modification_end();

    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(mesh, add_parts);

    stk::mesh::Selector skin = skin_part & locally_owned;

    std::vector<size_t> counts(4);
    counts[0] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::NODE_RANK));
    counts[1] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::EDGE_RANK));
    counts[2] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::FACE_RANK));
    counts[3] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::ELEM_RANK));

    if (p_rank == 0)
    {
      // check number of entities
      EXPECT_EQ( 1 == p_size ? 12u : 8u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ(                     0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 1 == p_size ? 16u : 8u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ(                     0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 8, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 4, 8, 7, 3 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 6, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 7, 3 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 2, 5 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 5, 6 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 1, 2, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 2, 3, 4 }}, mesh));

      if (1 == p_size) {
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 12, 8 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 8, 12, 11, 7 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 7 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 7, 11 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 10 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 10, 6 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 10, 12 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 11, 12 }}, mesh));
      }
    }
    else if (p_rank == 1)
    {
      // check number of entities in skin part
      EXPECT_EQ( 4u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ( 0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 8u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ( 0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 5, 9, 12, 8 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 8, 12, 11, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 7, 11 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 10 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 10, 6 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 10, 12 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 11, 12 }}, mesh));
    }
}

TEST( SkinMesh, SkinPyramidWithAura )
{
  test_skin_mesh_with_pyramid(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, SkinPyramidWithoutAura )
{
  test_skin_mesh_with_pyramid(stk::mesh::BulkData::NO_AUTO_AURA);
}

void test_skin_hybrid_mesh(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
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

    if(p_size > 2)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm, autoAuraOption);
    const int p_rank = mesh.parallel_rank();

    stk::mesh::EntityRank side_rank = meta.side_rank();

    stk::mesh::Part & skin_part = meta.declare_part("SkinPart", side_rank);
    stk::mesh::Part & locally_owned = meta.locally_owned_part();

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

    if (0 == p_rank) {
        for (size_t i = 0; i < numHex; ++i) {
          stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    if ( (1 == p_rank) || (1 == p_size) )  { // setup the pyramids/tets for either np 2 or serial
        for (size_t i = 0; i < numPyr; ++i) {
          stk::mesh::declare_element(mesh, *pyrPart, pyrElemIDs[i], pyrNodeIDs[i]);
        }
        for (size_t i = 0; i < numTet; ++i) {
          stk::mesh::declare_element(mesh, *tetPart, tetElemIDs[i], tetNodeIDs[i]);
        }
    }

    if (p_size > 1)
    {
      for (int nodeIdx = 0; nodeIdx < numSharedNodeTriples; ++nodeIdx) {
          if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
              stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
              int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
              stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
              mesh.add_node_sharing(node, sharingProc);
          }
      }
    }

    mesh.modification_end();

    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(stk::topology::NODE_RANK)) );
    ASSERT_EQ( 0u, stk::mesh::count_selected_entities( skin_part, mesh.buckets(side_rank)) );

    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(mesh, add_parts);

    stk::mesh::Selector skin = skin_part & locally_owned;

    std::vector<size_t> counts(4);
    counts[0] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::NODE_RANK));
    counts[1] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::EDGE_RANK));
    counts[2] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::FACE_RANK));
    counts[3] = stk::mesh::count_selected_entities( skin, mesh.buckets(stk::topology::ELEM_RANK));

    if (p_rank == 0)
    {
      // check number of entities
      EXPECT_EQ( 1 == p_size ? 12u : 8u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ(                     0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 1 == p_size ? 15u : 5u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ(                     0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 2, 3, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 8, 4 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 4, 8, 7, 3 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 1, 5, 6, 2 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 4 > {{ 2, 6, 7, 3 }}, mesh));

      if (1 == p_size)
      {
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 8 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 8, 12 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 7, 12 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 12, 11 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 12, 10 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 12, 11 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 6 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 9, 10 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 7 }}, mesh));
          EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 10, 11 }}, mesh));
      }
    }
    else if (p_rank == 1)
    {
      // check number of entities in skin part
      EXPECT_EQ(  4u, counts[stk::topology::NODE_RANK] ); // nodes
      EXPECT_EQ(  0u, counts[stk::topology::EDGE_RANK] );  // edges
      EXPECT_EQ( 10u, counts[stk::topology::FACE_RANK] );  // face
      EXPECT_EQ(  0u, counts[stk::topology::ELEM_RANK] );  // elements
      // check boundary faces are created
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 8 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 8, 12 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 8, 7, 12 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 12, 11 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 9, 12, 10 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 10, 12, 11 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 5, 9, 6 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 9, 10 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 6, 10, 7 }}, mesh));
      EXPECT_TRUE(check_if_one_owned_face_with_these_nodes_exists(std::array< uint64_t, 3 > {{ 7, 10, 11 }}, mesh));
    }
}

TEST( SkinMesh, SkinHybridMeshWithAura )
{
  test_skin_hybrid_mesh(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, SkinHybridMeshWithoutAura )
{
  test_skin_hybrid_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
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

void test_2_hex_2_block_with_second_selector(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
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

        stk::mesh::Selector onlyConsiderTheseElementsWhenFindingBoundary(!stk::mesh::Selector(*meta.get_part("block_1")));
        stk::mesh::skin_mesh(mesh, block_2, skin_parts, &onlyConsiderTheseElementsWhenFindingBoundary);

        stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEM_RANK, element_id);

        if (mesh.is_valid(element2))
        {
            unsigned num_faces = mesh.num_faces(element2);
            EXPECT_EQ(6u, num_faces);
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

TEST( SkinMesh, test_2_hex_2_block_with_second_selector_with_aura)
{
    test_2_hex_2_block_with_second_selector(stk::mesh::BulkData::AUTO_AURA);
}

TEST( SkinMesh, test_2_hex_2_block_with_second_selector_without_aura)
{
    test_2_hex_2_block_with_second_selector(stk::mesh::BulkData::NO_AUTO_AURA);
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
