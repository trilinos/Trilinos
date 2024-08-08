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

#include <stddef.h>                     // for NULL
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element, etc
#include <stk_mesh/base/MetaData.hpp>   // for get_cell_topology, MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/MeshCommVerify.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include <stk_unit_test_utils/BulkDataTester.hpp>

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct EntitySideComponent; } }

using stk::ParallelMachine;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;
using stk::mesh::EntitySideComponent;
using stk::mesh::PairIterRelation;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::unit_test_util::build_mesh;

class TopologyHelpersTestingFixture
{
public:
  TopologyHelpersTestingFixture(ParallelMachine pm);
  ~TopologyHelpersTestingFixture() {}

  const int spatial_dimension;
  std::shared_ptr<BulkData> bulkPtr;
  BulkData& bulk;
  MetaData& meta;
  const EntityRank element_rank;
  const EntityRank side_rank;
  Part & generic_element_part;
  Part & element_tet_part;
  Part & element_wedge_part;
  Part & element_hex_part;
  Part & generic_face_part;
  Part & another_generic_face_part;
  Part & face_tri_part;
  Part & face_quad_part;
  Part & another_generic_element_part;

  EntityId nextEntityId()
  { return psize*(++entity_id)+prank; }

  Entity create_entity( EntityRank rank, Part& part_membership)
  {
    PartVector part_intersection;
    part_intersection.push_back ( &part_membership );
    if (rank != meta.side_rank()) {
      return bulk.declare_entity(rank, nextEntityId(), part_intersection);
    }
    return Entity();
  }

private:
  EntityId entity_id;
  const int psize;
  const int prank;
};

TopologyHelpersTestingFixture::TopologyHelpersTestingFixture(ParallelMachine pm)
  : spatial_dimension( 3 ),
    bulkPtr( build_mesh(3, pm, stk::mesh::BulkData::AUTO_AURA) ),
    bulk( *bulkPtr ),
    meta( bulkPtr->mesh_meta_data() ),
    element_rank( stk::topology::ELEMENT_RANK ),
    side_rank( meta.side_rank()),
    generic_element_part( meta.declare_part("another part", element_rank ) ),
    element_tet_part( meta.declare_part_with_topology( "block_left_1", stk::topology::TET_4 ) ),
    element_wedge_part( meta.declare_part_with_topology( "block_left_2", stk::topology::WEDGE_15 ) ),
    element_hex_part( meta.declare_part_with_topology( "block_left_3", stk::topology::HEX_8 ) ),
    generic_face_part( meta.declare_part_with_topology( "A_1", stk::topology::QUAD_4 ) ),
    another_generic_face_part( meta.declare_part("A_2", side_rank ) ),
    face_tri_part( meta.declare_part_with_topology("A_3_0", stk::topology::TRI_3)),
    face_quad_part( meta.declare_part("A_3", side_rank ) ),
    another_generic_element_part( meta.declare_part("B_3", element_rank ) ),
    entity_id(0u),
    psize(bulk.parallel_size()),
    prank(bulk.parallel_rank())
{
  meta.commit();
}

namespace {

TEST( testTopologyHelpers, get_cell_topology_based_on_part)
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  if (fix.bulk.parallel_size() != 1) {
    return;
  }

  fix.bulk.modification_begin();
  Entity elem1  = fix.create_entity( fix.element_rank, fix.element_hex_part );
  std::vector<Entity> elem_node(8);
  for (int i = 0; i < 8; ++i) {
    elem_node[i] = fix.bulk.declare_node(100 + i);
    fix.bulk.declare_relation(elem1, elem_node[i], i);
  }

  Entity side1 = fix.bulk.declare_element_side(elem1, 0, stk::mesh::ConstPartVector{&fix.generic_face_part});

  PartVector tmp(1);
  tmp[0] = & fix.face_quad_part;
  fix.bulk.change_entity_parts ( side1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(side1).topology(), stk::topology::QUAD_4 );
  fix.bulk.change_entity_parts ( side1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(side1).topology(), stk::topology::QUAD_4 );
  tmp[0] = & fix.another_generic_face_part;
  fix.bulk.change_entity_parts ( side1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(side1).topology(), stk::topology::QUAD_4 );
  ASSERT_NE( fix.bulk.bucket(side1).topology(), stk::topology::WEDGE_15 );

  fix.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_no_topology )
{
  // Coverage for declare_element_side - TopologyHelpers.cpp - "Cannot discern element topology"
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();
  //fix.bulk.modification_end();


  {
    stk::mesh::EntityIdVector elem_node {1, 2, 3, 4};
    //fix.bulk.modification_begin();
    // Cannot declare an element without a topology defined
    ASSERT_THROW(
          stk::mesh::declare_element(fix.bulk, fix.generic_element_part, fix.nextEntityId(), elem_node),
          std::runtime_error
          );
  }
}

TEST( testTopologyHelpers, declare_element_side_wrong_bulk_data)
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "BulkData for 'elem' and 'side' are different"
  TopologyHelpersTestingFixture fix1(MPI_COMM_WORLD);

  fix1.bulk.modification_begin();

  TopologyHelpersTestingFixture fix2(MPI_COMM_WORLD);
  fix2.bulk.modification_begin();
  fix2.bulk.modification_end();
}


TEST( testTopologyHelpers, declare_element_side_full )
{
  // Go all way the through declare_element_side - use new element
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();

  stk::mesh::EntityIdVector elem_node {1, 2, 3, 4};

  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );

  const EntityId zero_side_count = 0;
  Entity face2 = fix.bulk.declare_element_side(element, zero_side_count, stk::mesh::PartVector{&fix.face_tri_part});
  fix.bulk.modification_end();

  const stk::mesh::ConnectedEntities rel2_nodes = fix.bulk.get_connected_entities(face2, stk::topology::NODE_RANK);
  ASSERT_TRUE(rel2_nodes.size() != 0u);

  ASSERT_TRUE( true );  // This test is to check compilation.
}

TEST( testTopologyHelpers, element_side_polarity_valid )
{
  // Coverage of element_side_polarity in TopologyHelpers.cpp 168-181 and 200-215
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  stk::mesh::EntityIdVector elem_node {1, 2, 3, 4};

  fix.bulk.modification_begin();
  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  const EntityId zero_side_count = 0;
  Entity face2 = fix.bulk.declare_element_side(element, zero_side_count, stk::mesh::PartVector{&fix.face_tri_part});
  fix.bulk.modification_end();

  const int local_side_id = 0;
  ASSERT_TRUE( stk::mesh::element_side_polarity(fix.bulk, element, face2, local_side_id) );

}

TEST( testTopologyHelpers, element_side_polarity_invalid_1 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  stk::mesh::EntityIdVector elem_node {1, 2, 3, 4};

  // Coverage of element_side_polarity in TopologyHelpers.cpp
  {
    fix.bulk.modification_begin();
    Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
    const EntityId zero_side_count = 0;
    Entity face = fix.bulk.declare_element_side(element, zero_side_count, stk::mesh::PartVector{&fix.face_tri_part});
    fix.bulk.modification_end();

    const unsigned invalid_local_side_id = static_cast<unsigned>(-1);
    // Hits "Unsuported local_side_id" error condition:
    ASSERT_THROW(
          stk::mesh::element_side_polarity(fix.bulk, element, face, invalid_local_side_id),
          std::runtime_error
          );
  }
}

TEST( testTopologyHelpers, element_side_polarity_invalid_2 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  stk::mesh::EntityIdVector elem_node { 1, 2, 3, 4 };

  // Coverage of element_side_polarity in TopologyHelpers.cpp
  fix.bulk.modification_begin();

  stk::mesh::ConstPartVector part_intersection;
  part_intersection.push_back ( &fix.generic_element_part);
  Entity element = fix.bulk.declare_element(fix.nextEntityId(), part_intersection);
  ASSERT_TRUE( fix.bulk.bucket(element).topology() == stk::topology::INVALID_TOPOLOGY );

  Entity element_with_top = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  ASSERT_TRUE( fix.bulk.bucket(element_with_top).topology() != stk::topology::INVALID_TOPOLOGY );

  const EntityId zero_side_count = 0;
  Entity face_with_top = fix.bulk.declare_element_side(element_with_top, zero_side_count, stk::mesh::PartVector{&fix.face_tri_part});
  const int valid_local_side_id = 0;
  ASSERT_THROW(
        element_side_polarity(fix.bulk, element, face_with_top, valid_local_side_id),
        std::runtime_error
        );

  //modification_end is not called due to difference in expected behavior for release and debug builds - debug should throw, release should not
  //difference occurs within check_for_connected_nodes method
  //ASSERT_THROW(fix.bulk.modification_end(), std::logic_error);

}

TEST(stk_topology_permutations, lexicographical_smallest_permutation_preserve_polarity)
{
  {
    stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;
    unsigned shell_node_ids[3] = {10, 8, 12};
    {
      unsigned triangle_node_ids[3] = {12, 10, 8};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 2;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned triangle_node_ids[3] = {10, 8, 12};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 2;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned triangle_node_ids[3] = {8, 12, 10};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 2;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned triangle_node_ids[3] = {12, 8, 10};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 5;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned triangle_node_ids[3] = {10, 12, 8};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 5;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned triangle_node_ids[3] = {8, 10, 12};

      unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 5;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
  }
}

TEST(stk_topology_permutations, quad_lexicographical_smallest_permutation_preserve_polarity)
{
  {
    stk::topology quad_shell = stk::topology::SHELL_QUAD_4;
    unsigned shell_node_ids[4] = {1, 2, 3, 4};
    {
      unsigned quad_node_ids[4] = {1, 2, 3, 4};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 0;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {4, 1, 2, 3};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 0;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {3, 4, 1, 2};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 0;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {2, 3, 4, 1};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 0;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {1, 4, 3, 2};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 4;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {2, 1, 4, 3};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 4;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {3, 2, 1, 4};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 4;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {4, 3, 2, 1};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 4;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
    {
      unsigned quad_node_ids[4] = {4, 2, 3, 1};

      unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
      unsigned gold_lexicographical_smallest_permutation_index = 8;
      // driven by vertices, NOT mid-edge nodes
      EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
  }
}

int check_permutation_given(stk::mesh::BulkData& mesh, stk::mesh::Entity elem, unsigned face_ord, stk::mesh::Permutation claimed_permutation, stk::mesh::Entity elem_face)
{
  stk::mesh::Entity face_nodes_buff[100], mapped_face_nodes_buff[100];

  stk::topology elem_topo = mesh.bucket(elem).topology();
  const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);

  elem_topo.face_nodes(elem_nodes, face_ord, face_nodes_buff);
  stk::topology face_topo = elem_topo.face_topology(face_ord);

  stk::mesh::ConnectivityOrdinal const *face_ordinals_for_elem = mesh.begin_face_ordinals(elem);

  const unsigned num_faces = mesh.num_faces(elem);
  EXPECT_GT(num_faces, 0u);
  unsigned face_idx;
  for (face_idx = 0; face_idx < num_faces; ++face_idx)
  {
    if (face_ord == face_ordinals_for_elem[face_idx])
      break;
  }
  EXPECT_LT(face_idx, num_faces);

  stk::mesh::Permutation const *face_perms_for_elem = mesh.begin_face_permutations(elem);
  stk::mesh::Permutation face_perm = face_perms_for_elem[face_idx];
  face_topo.permutation_nodes(face_nodes_buff, face_perm, mapped_face_nodes_buff);

  // Another way to get to face's nodes.
  const stk::mesh::ConnectedEntities face_nodes = mesh.get_connected_entities(elem_face, stk::topology::NODE_RANK);
  STK_ThrowAssert(face_nodes.size() == face_topo.num_nodes());

  int innermost_hits = 0;
  for (unsigned ni = 0; ni < face_topo.num_nodes(); ++ni)
  {
    ++innermost_hits;
    EXPECT_EQ(face_nodes[ni], mapped_face_nodes_buff[ni]);
  }

  // Indeed, find_permutation computes what was stored!
  stk::mesh::Permutation perm = stk::mesh::find_permutation(mesh, elem_topo, elem_nodes, face_topo, face_nodes.data(), face_ord);
  EXPECT_EQ(perm, claimed_permutation);
  return innermost_hits;
}

TEST (stkTopologyFunctions, use_permutations_Hex_2x1x1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    const size_t NX = 2;
    const size_t NY = 1;
    const size_t NZ = 1;

    stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

    fixture.m_meta.commit();
    fixture.generate_mesh();

    stk::mesh::BulkData &mesh = fixture.m_bulk_data;
    stk::mesh::create_faces(mesh);

    const std::vector<stk::mesh::Bucket *> &elem_buckets =
        mesh.get_buckets(stk::topology::ELEMENT_RANK, fixture.m_meta.universal_part());

    int innermost_hits = 0;
    int num_elems      = 0;

    for (unsigned bkt_idx = 0; bkt_idx < elem_buckets.size(); ++bkt_idx)
    {
      const stk::mesh::Bucket &elem_bkt = *elem_buckets[bkt_idx];
      stk::topology elem_topo           = elem_bkt.topology();

      unsigned num_faces = elem_topo.num_faces();

      for (unsigned elem_idx = 0; elem_idx < elem_bkt.size(); ++elem_idx)
      {
        ++num_elems;

        stk::mesh::Permutation const *elem_face_perms = elem_bkt.begin_face_permutations(elem_idx);
        stk::mesh::Entity const *elem_faces           = elem_bkt.begin_faces(elem_idx);

        for (unsigned face_ord = 0; face_ord < num_faces; ++face_ord)
        {
          innermost_hits += check_permutation_given(mesh, elem_bkt[elem_idx], face_ord, elem_face_perms[face_ord], elem_faces[face_ord]);
        }
      }
    }
    std::ostringstream msg_buff;
    msg_buff << "P" << mesh.parallel_rank() << ": knows " << num_elems << " elements" << std::endl;
    std::cout << msg_buff.str();
    EXPECT_EQ(innermost_hits, num_elems * 24);

    if (num_elems > 0)
    {
      EXPECT_EQ(num_elems, 2);
    }
  }
}

void test_side_creation(unsigned *gold_side_ids,unsigned local_side_id)
{
  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  std::string name = "generated:1x1x1";
  stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

  stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, 1);

  mesh.modification_begin();
  stk::mesh::Part &quad4_part = mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
  stk::mesh::Entity side = mesh.declare_element_side(elem, local_side_id, stk::mesh::PartVector{&quad4_part});
  mesh.modification_end();

  stk::mesh::Permutation identity_permutation = static_cast<stk::mesh::Permutation>(0);
  check_permutation_given(mesh, elem, local_side_id, identity_permutation, side);

  //    std::string out_filename = "testStk.exo";
  //    size_t resultFileIndex = stkMeshIoBroker.create_output_mesh(out_filename, stk::io::WRITE_RESULTS);
  //    stkMeshIoBroker.write_output_mesh(resultFileIndex);

  const stk::mesh::Entity *side_nodes = mesh.begin_nodes(side);
  unsigned num_nodes = mesh.num_nodes(side);
  ASSERT_EQ(4u, num_nodes);

  for(unsigned i = 0; i < num_nodes; ++i)
  {
    EXPECT_EQ(gold_side_ids[i], mesh.identifier(side_nodes[i]));
  }
}

void test_side_creation_with_permutation(unsigned *gold_side_ids,unsigned local_side_id, stk::mesh::Permutation perm)
{
  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  std::string name = "generated:1x1x1";
  stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

  stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, 1);

  stk::mesh::PartVector parts;
  parts.push_back(&mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4));
  unsigned gold_num_nodes = 4;
  stk::mesh::EntityVector nodes(gold_num_nodes);
  for(unsigned i = 0; i < gold_num_nodes; ++i) {
    stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_ids[i]);
    nodes[i] = node;
  }

  mesh.modification_begin();
  stk::mesh::Entity side = mesh.declare_element_side(elem, local_side_id, parts);
  mesh.modification_end();

  check_permutation_given(mesh, elem, local_side_id, perm, side);

  const stk::mesh::Entity *side_nodes = mesh.begin_nodes(side);
  unsigned num_nodes = mesh.num_nodes(side);
  ASSERT_EQ(gold_num_nodes, num_nodes);

  for(unsigned i = 0; i < num_nodes; ++i)
  {
    //unsigned ordinal = nodes.size() - i - 1;
    EXPECT_EQ(gold_side_ids[i], mesh.identifier(side_nodes[i]));
  }
}

TEST(stkTopologyFunctions, check_permutation_with_FEM_Helper)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    unsigned gold_side_ids[6][4] = {
      {1,2,6,5},
      {2,4,8,6},
      {4,3,7,8},
      {1,5,7,3},
      {1,3,4,2},
      {5,6,8,7}
    };

    for(unsigned local_side_id = 0; local_side_id < 6; ++local_side_id)
    {
      test_side_creation(gold_side_ids[local_side_id],local_side_id);
    }
  }
}

TEST(stkTopologyFunctions, check_permutation_without_FEM_Helper)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    unsigned gold_side_ids[6][4] = {
      {1,2,6,5},
      {2,4,8,6},
      {4,3,7,8},
      {1,5,7,3},
      {1,3,4,2},
      {5,6,8,7}
    };
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
    unsigned local_side_id = 0;
    test_side_creation_with_permutation(gold_side_ids[local_side_id],local_side_id, perm);
  }
}

TEST(stkTopologyFunctions, check_permutations_for_Hex_1x1x1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    unsigned gold_side_ids[4] = {5,6,8,7};
    unsigned local_side_id = 5; // nodes 2,4,8,5 are nodes of face 'local_side_id' of a 1x1x1 hex element (generated)

    unsigned gold_permutations[8][4] = {
      {5,6,8,7},
      {7,5,6,8},
      {8,7,5,6},
      {6,8,7,5},
      {5,7,8,6},
      {7,8,6,5},
      {8,6,5,7},
      {6,5,7,8}
    };

    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    std::string name = "generated:1x1x1";
    stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, 1);

    stk::mesh::PartVector parts;
    parts.push_back(&mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4));
    unsigned gold_num_nodes = 4;
    stk::mesh::EntityVector nodes(gold_num_nodes);
    for(unsigned i = 0; i < gold_num_nodes; ++i) {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_ids[i]);
      nodes[i] = node;
    }

    mesh.modification_begin();
    stk::mesh::Entity side = mesh.declare_element_side(elem, local_side_id, parts);
    mesh.modification_end();

    stk::mesh::EntityVector mapped_face_nodes(gold_num_nodes);

    stk::topology elem_topo = mesh.bucket(elem).topology();
    const stk::mesh::Entity* face_nodes = mesh.begin_nodes(side);

    stk::topology face_topo = elem_topo.face_topology(local_side_id);
    unsigned num_permutations = face_topo.num_permutations();
    for(unsigned i = 0; i < num_permutations; ++i) {
      stk::mesh::Permutation local_perm = static_cast<stk::mesh::Permutation>(i);
      face_topo.permutation_nodes(face_nodes, local_perm, &mapped_face_nodes[0]);

      for (unsigned j = 0; j < gold_num_nodes; ++j) {
        EXPECT_EQ(gold_permutations[i][j], mesh.identifier(mapped_face_nodes[j])) << " for permutation " << i;
      }
    }
  }
}

void pack_downward_relations_for_entity(stk::mesh::BulkData& mesh, stk::mesh::Entity some_entity, std::vector<stk::mesh::Relation>& recv_relations1)
{
  unsigned bucket_ordinal = mesh.bucket_ordinal(some_entity);
  const stk::mesh::Bucket& bucket = mesh.bucket(some_entity);
  for(EntityRank irank=stk::topology::BEGIN_RANK; irank<mesh.entity_rank(some_entity);++irank)
  {
    Entity const *rels_itr = bucket.begin(bucket_ordinal, irank);
    Entity const *rels_end = bucket.end(bucket_ordinal, irank);
    stk::mesh::ConnectivityOrdinal const *ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);

    for(;rels_itr!=rels_end;++rels_itr,++ords_itr)
    {
      recv_relations1.push_back(stk::mesh::Relation(*rels_itr, mesh.entity_rank(*rels_itr), *ords_itr));
    }

  }
}

TEST(stkTopologyFunctions, permutation_consistency_check_3d)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    unsigned gold_side_ids[4] = {2,4,8,6};
    unsigned local_side_id = 1; // nodes 2,4,8,5 are nodes of face 'local_side_id' of a 1x1x1 hex element (generated)

    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    std::string name = "generated:1x1x1";
    stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, 1);

    stk::mesh::PartVector parts;
    parts.push_back(&mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4));
    unsigned gold_num_nodes = 4;
    stk::mesh::EntityVector nodes(gold_num_nodes);
    for(unsigned i = 0; i < gold_num_nodes; ++i) {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_ids[i]);
      nodes[i] = node;
    }

    mesh.modification_begin();
    stk::mesh::Entity side = mesh.declare_element_side(elem, local_side_id, parts);
    mesh.modification_end();

    bool rel_bad = false;

    std::vector<stk::mesh::Relation> recv_relations;
    pack_downward_relations_for_entity(mesh, side, recv_relations);
    stk::mesh::impl::unpack_not_owned_verify_compare_closure_relations(mesh, side, recv_relations, rel_bad);
    EXPECT_FALSE(rel_bad);

    std::vector<stk::mesh::Relation> relations;
    pack_downward_relations_for_entity(mesh, elem, relations);
    stk::mesh::impl::unpack_not_owned_verify_compare_closure_relations(mesh, elem, relations, rel_bad);
    EXPECT_FALSE(rel_bad);
  }
}

TEST(stkTopologyFunctions, check_permutation_consistency_parallel)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    unsigned gold_side_ids[4] = {5,6,8,7};

    stk::mesh::MetaData meta(3);
    stk::unit_test_util::BulkDataFaceSharingTester mesh(meta, MPI_COMM_WORLD);

    const std::string generatedMeshSpec = "generated:1x1x2";
    stk::io::fill_mesh(generatedMeshSpec, mesh);

    unsigned elem_id = 0;
    unsigned local_side_id = 0;

    if (mesh.parallel_rank()==0)
    {
      local_side_id=5;
      elem_id = 1;
    }
    else
    {
      local_side_id=4;
      elem_id = 2;
    }

    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);
    EXPECT_TRUE(mesh.bucket(elem).owned());

    stk::mesh::PartVector parts;
    parts.push_back(&mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4));
    unsigned gold_num_nodes = 4;
    stk::mesh::EntityVector nodes(gold_num_nodes);
    for(unsigned i = 0; i < gold_num_nodes; ++i) {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_ids[i]);
      nodes[i] = node;
    }

    mesh.modification_begin();
    mesh.declare_element_side(elem, local_side_id, parts);
    mesh.modification_end();

    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(mesh, mesh_counts);
    size_t numFacesGlobal = 1u;
    EXPECT_EQ(numFacesGlobal, mesh_counts[stk::topology::FACE_RANK]);
  }
}

TEST(stkTopologyFunctions, permutation_consistency_check_2d)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(2, MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;

    mesh.modification_begin();
    stk::mesh::Entity Quad9 = mesh.declare_element(1, stk::mesh::ConstPartVector{&meta.get_topology_root_part(stk::topology::QUAD_9_2D)});
    unsigned node_ids[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    stk::mesh::EntityVector nodes(9);
    for(unsigned i=0;i<9;++i)
    {
      nodes[i] = mesh.declare_node(node_ids[i]);
    }

    for(size_t i=0;i<nodes.size();++i)
    {
      mesh.declare_relation(Quad9, nodes[i], i);
    }

    mesh.modification_end();

    mesh.modification_begin();

    stk::mesh::EntityVector sides(4);
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
    for(size_t i=0;i<sides.size();++i)
    {
      sides[i] = mesh.declare_element_side(Quad9, i, stk::mesh::ConstPartVector{&meta.get_topology_root_part(stk::topology::LINE_3)});
    }

    EXPECT_NO_THROW(mesh.modification_end());

    EXPECT_TRUE(stk::mesh::check_permutation(mesh, Quad9, sides[0], 0, perm)) << "for side 0";
    EXPECT_TRUE(stk::mesh::check_permutation(mesh, Quad9, sides[1], 1, perm)) << "for side 1";
    EXPECT_TRUE(stk::mesh::check_permutation(mesh, Quad9, sides[2], 2, perm)) << "for side 2";
    EXPECT_TRUE(stk::mesh::check_permutation(mesh, Quad9, sides[3], 3, perm)) << "for side 3";

    EXPECT_TRUE(stk::mesh::impl::check_permutations_on_all(mesh));
  }
}

struct SuperTopologySideData
{
  stk::mesh::EntityIdVector elemIDsPerProc;
  std::vector<stk::mesh::EntityIdVector> nodeIDsPerProc;
  stk::mesh::EntityIdVector sharedNodeIds;
  stk::mesh::EntityId sharedFaceId;
  std::vector<stk::mesh::ConnectivityOrdinal> ordinalPerProc;
  std::vector<stk::mesh::Permutation> permPerProc;
};

class SuperTopologies : public stk::unit_test_util::MeshFixture
{
protected:
  SuperTopologies(int dim)
    : stk::unit_test_util::MeshFixture(dim)
  {
  }

  void expect_mesh_correct(stk::topology superElem, stk::topology superSide, const SuperTopologySideData &s)
  {
    expect_entity_has_topology(stk::topology::ELEM_RANK, s.elemIDsPerProc[0], superElem);
    expect_entity_has_topology(stk::topology::ELEM_RANK, s.elemIDsPerProc[1], superElem);
    expect_entity_has_topology(get_meta().side_rank(), s.sharedFaceId, superSide);
    expect_elem_to_side_connections(superSide, s.elemIDsPerProc[0], s.ordinalPerProc[0], s.permPerProc[0]);
    expect_elem_to_side_connections(superSide, s.elemIDsPerProc[1], s.ordinalPerProc[1], s.permPerProc[1]);
  }

  void expect_entity_has_topology(stk::mesh::EntityRank rank, stk::mesh::EntityId id, stk::topology topo)
  {
    stk::mesh::Entity entity = get_bulk().get_entity(rank, id);
    EXPECT_EQ(topo, get_bulk().bucket(entity).topology());
    EXPECT_EQ(topo.num_nodes(), get_bulk().num_nodes(entity));
  }

  void expect_elem_to_side_connections(stk::topology superSide,
                                       stk::mesh::EntityId elemId,
                                       stk::mesh::ConnectivityOrdinal expectedOrdinal,
                                       stk::mesh::Permutation expectedPerm)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    unsigned numSides = get_bulk().num_sides(elem);
    EXPECT_EQ(1u, numSides);
    const stk::mesh::Entity side = *get_bulk().begin(elem, get_meta().side_rank());
    EXPECT_EQ(superSide, get_bulk().bucket(side).topology());
    EXPECT_EQ(superSide.num_nodes(), get_bulk().num_nodes(side));
    const stk::mesh::ConnectivityOrdinal elemSide = *get_bulk().begin_ordinals(elem, get_meta().side_rank());
    EXPECT_EQ(expectedOrdinal, elemSide);
    const stk::mesh::Permutation perm = *get_bulk().begin_permutations(elem, get_meta().side_rank());
    EXPECT_EQ(expectedPerm, perm);
  }
  void create_mesh(stk::topology superElem, stk::topology superSide, const SuperTopologySideData &s)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    int procId = get_bulk().parallel_rank();
    stk::mesh::PartVector parts = {&get_meta().declare_part_with_topology("superElempart",superElem)};
    stk::mesh::PartVector sideParts = {&get_meta().declare_part_with_topology("superSidepart",superSide)};

    get_bulk().modification_begin();
    stk::mesh::declare_element(get_bulk(), parts, s.elemIDsPerProc[procId], s.nodeIDsPerProc[procId]);
    add_shared_nodes(s.sharedNodeIds);
    get_bulk().modification_end();
    create_sides(superSide, sideParts, s);
  }

  void create_sides(stk::topology superSide, const stk::mesh::PartVector& parts, const SuperTopologySideData &s)
  {
    get_bulk().modification_begin();
    stk::mesh::Entity side = get_bulk().declare_solo_side(s.sharedFaceId, parts);
    for(unsigned i=0; i<s.sharedNodeIds.size(); ++i) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, s.sharedNodeIds[i]);
      get_bulk().declare_relation(side, node, i);
    }

    int procId = get_bulk().parallel_rank();
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, s.elemIDsPerProc[procId]);
    get_bulk().declare_relation(elem, side, s.ordinalPerProc[procId], s.permPerProc[procId]);
    get_bulk().modification_end();
  }

  void add_shared_nodes(const stk::mesh::EntityIdVector &sharedNodeIds)
  {
    int otherProc = 1-get_bulk().parallel_rank();
    for(stk::mesh::EntityId nodeId : sharedNodeIds) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      get_bulk().add_node_sharing(node, otherProc);
    }
  }

};

class SuperTopologies3d : public SuperTopologies
{
protected:
  SuperTopologies3d() : SuperTopologies(3)
  {
    superTopo3d.elemIDsPerProc = {1, 2};
    superTopo3d.nodeIDsPerProc = { {1, 2, 3, 4, 5, 6, 7, 8}, {5, 6, 7, 8, 9, 10, 11, 12} };
    superTopo3d.sharedNodeIds = {5, 6, 7, 8};
    superTopo3d.sharedFaceId = 1;
    superTopo3d.ordinalPerProc = {stk::mesh::ConnectivityOrdinal(4), stk::mesh::ConnectivityOrdinal(5)};
    superTopo3d.permPerProc = {stk::mesh::Permutation(0), stk::mesh::Permutation(4)};
  }
  struct SuperTopologySideData superTopo3d;
};

TEST_F(SuperTopologies3d, twoElemsTwoProcs)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    stk::topology super8 = stk::create_superelement_topology(8);
    stk::topology superface4 = stk::create_superface_topology(4);
    create_mesh(super8, superface4, superTopo3d);
    expect_mesh_correct(super8, superface4, superTopo3d);
  }
}

TEST_F(SuperTopologies3d, twoElemsTwoProcsElemGraph)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    stk::topology super8 = stk::create_superelement_topology(8);
    stk::topology superface4 = stk::create_superface_topology(4);
    create_mesh(super8, superface4, superTopo3d);
    stk::mesh::ElemElemGraph elemGraph(get_bulk());
    EXPECT_EQ(0u, elemGraph.num_edges());
    EXPECT_EQ(0u, elemGraph.num_parallel_edges());
  }
}

class SuperTopologies2d : public SuperTopologies
{
protected:
  SuperTopologies2d() : SuperTopologies(2)
  {
    superTopo2d.elemIDsPerProc = {1, 2};
    superTopo2d.nodeIDsPerProc = { {1, 2, 3, 4}, {4, 3, 5, 6} };
    superTopo2d.sharedNodeIds = {3, 4};
    superTopo2d.sharedFaceId = 1;
    superTopo2d.ordinalPerProc = {stk::mesh::ConnectivityOrdinal(3), stk::mesh::ConnectivityOrdinal(0)};
    superTopo2d.permPerProc = {stk::mesh::Permutation(0), stk::mesh::Permutation(1)};
  }
  struct SuperTopologySideData superTopo2d;
};

TEST_F(SuperTopologies2d, twoElemsTwoProcs)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    stk::topology superElem4 = stk::create_superelement_topology(4);
    stk::topology superEdge2 = stk::create_superedge_topology(2);
    create_mesh(superElem4, superEdge2, superTopo2d);
    expect_mesh_correct(superElem4, superEdge2, superTopo2d);
  }
}

}

