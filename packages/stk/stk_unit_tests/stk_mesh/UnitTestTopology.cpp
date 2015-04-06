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

#include <stddef.h>                     // for NULL
#include <Shards_BasicTopologies.hpp>   // for getCellTopologyData, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element, etc
#include <stk_mesh/base/MetaData.hpp>   // for get_cell_topology, MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include "Shards_CellTopologyData.h"    // for CellTopologyData
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct EntitySideComponent; } }

namespace stk {
namespace mesh {
void unpack_not_owned_verify_compare_closure_relations( const BulkData &             mesh,
                                                        Entity                       entity,
                                                        std::vector<stk::mesh::Relation> const& recv_relations,
                                                        bool&                        bad_rel);
}
}

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

class TopologyHelpersTestingFixture
{
 public:
  TopologyHelpersTestingFixture(ParallelMachine pm);
  ~TopologyHelpersTestingFixture() {}

  const int spatial_dimension;
  MetaData meta;
  BulkData bulk;
  const EntityRank element_rank;
  const EntityRank side_rank;
  Part & generic_element_part;
  Part & element_tet_part;
  Part & element_wedge_part;
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
    return bulk.declare_entity(rank, nextEntityId(), part_intersection);
  }

 private:
  EntityId entity_id;
  const int psize;
  const int prank;
};

TopologyHelpersTestingFixture::TopologyHelpersTestingFixture(ParallelMachine pm)
  : spatial_dimension( 3 )
  , meta( spatial_dimension )
  , bulk( meta, pm, 100 )
  , element_rank( stk::topology::ELEMENT_RANK )
  , side_rank( meta.side_rank())
  , generic_element_part( meta.declare_part("another part", element_rank ) )
  , element_tet_part( meta.declare_part_with_topology( "block_left_1", stk::topology::TET_4 ) )
  , element_wedge_part( meta.declare_part_with_topology( "block_left_2", stk::topology::WEDGE_15 ) )
  , generic_face_part( meta.declare_part_with_topology( "A_1", stk::topology::QUAD_4 ) )
  , another_generic_face_part( meta.declare_part("A_2", side_rank ) )
  , face_tri_part( meta.declare_part_with_topology("A_3_0", stk::topology::TRI_3))
  , face_quad_part( meta.declare_part("A_3", side_rank ) )
  , another_generic_element_part( meta.declare_part("B_3", element_rank ) )
  , entity_id(0u)
  , psize(bulk.parallel_size())
  , prank(bulk.parallel_rank())
{
  meta.commit();
}

namespace {

TEST( testTopologyHelpers, get_cell_topology_based_on_part)
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();
  Entity elem1  = fix.create_entity( fix.side_rank, fix.generic_face_part );

  std::vector<Entity> elem_node(4);
  for (int i = 0; i < 4; ++i) {
    elem_node[i] = fix.bulk.declare_entity(stk::topology::NODE_RANK, 100 + i);
    fix.bulk.declare_relation(elem1, elem_node[i], i);
  }

  PartVector tmp(1);
  tmp[0] = & fix.face_quad_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(elem1).topology(), stk::topology::QUAD_4 );
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(elem1).topology(), stk::topology::QUAD_4 );
  tmp[0] = & fix.another_generic_face_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( fix.bulk.bucket(elem1).topology(), stk::topology::QUAD_4 );
  ASSERT_NE( fix.bulk.bucket(elem1).topology(), stk::topology::WEDGE_15 );

  fix.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_no_topology )
{
  // Coverage for declare_element_side - TopologyHelpers.cpp - "Cannot discern element topology"
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();
  Entity elem4  = fix.create_entity( fix.element_rank , fix.generic_element_part );
  ASSERT_THROW(
    stk::mesh::declare_element_side( fix.bulk, fix.element_rank, elem4, fix.nextEntityId(), &fix.element_wedge_part ),
    std::runtime_error
      );
  //fix.bulk.modification_end();


  {
    EntityId elem_node[4];
    elem_node[0] = 1;
    elem_node[1] = 2;
    elem_node[2] = 3;
    elem_node[3] = 4;
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

TEST( testTopologyHelpers, declare_element_side_no_topology_2 )
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "No element topology found and cell side id exceeds..."
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;
  Entity element  = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node);
  stk::topology elem_top = fix.bulk.bucket(element).topology();
  const EntityId nSideCount = elem_top.num_sides() + 10 ;
  ASSERT_THROW(
    stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, nSideCount, &fix.element_tet_part ),
    std::runtime_error
      );
  fix.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_full )
{
  // Go all way the through declare_element_side - use new element
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );

  const EntityId zero_side_count = 0;
  Entity face2 = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                  &fix.face_tri_part);
  fix.bulk.modification_end();

  stk::mesh::Entity const *rel2_nodes = fix.bulk.begin_nodes(face2);
  ASSERT_TRUE(rel2_nodes != 0);

  ASSERT_TRUE( true );  // This test is to check compilation.
}

TEST( testTopologyHelpers, element_side_polarity_valid )
{
  // Coverage of element_side_polarity in TopologyHelpers.cpp 168-181 and 200-215
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  fix.bulk.modification_begin();
  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  const EntityId zero_side_count = 0;
  Entity face2 = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                  &fix.face_tri_part);
  fix.bulk.modification_end();

  const int local_side_id = 0;
  ASSERT_TRUE( fix.bulk.element_side_polarity( element, face2, local_side_id) );

}

TEST( testTopologyHelpers, element_side_polarity_invalid_1 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  // Coverage of element_side_polarity in TopologyHelpers.cpp
  {
    fix.bulk.modification_begin();
    Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
    const EntityId zero_side_count = 0;
    Entity face = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                   &fix.face_tri_part);
    fix.bulk.modification_end();

    const unsigned invalid_local_side_id = static_cast<unsigned>(-1);
    // Hits "Unsuported local_side_id" error condition:
    ASSERT_THROW(
        fix.bulk.element_side_polarity( element, face, invalid_local_side_id),
        std::runtime_error
        );
  }
}

TEST( testTopologyHelpers, element_side_polarity_invalid_2 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  // Coverage of element_side_polarity in TopologyHelpers.cpp
  fix.bulk.modification_begin();

  PartVector part_intersection;
  part_intersection.push_back ( &fix.generic_element_part);
  Entity element = fix.bulk.declare_entity(fix.element_rank, fix.nextEntityId(), part_intersection);
  ASSERT_TRUE( fix.bulk.bucket(element).topology() == stk::topology::INVALID_TOPOLOGY );

  Entity element_with_top = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  ASSERT_TRUE( fix.bulk.bucket(element_with_top).topology() != stk::topology::INVALID_TOPOLOGY );

  const EntityId zero_side_count = 0;
  Entity face_with_top = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element_with_top, zero_side_count,
                                                          &fix.face_tri_part);
  const int valid_local_side_id = 0;
  ASSERT_THROW(
      fix.bulk.element_side_polarity( element, face_with_top, valid_local_side_id),
      std::runtime_error
      );

  //modification_end is not called due to difference in expected behavior for release and debug builds - debug should throw, release should not
  //difference occurs within check_for_connected_nodes method
  //ASSERT_THROW(fix.bulk.modification_end(), std::logic_error);

  // Hits "Element has no defined topology" error condition:
  //ASSERT_TRUE( stk::mesh::get_cell_topology( fix.bulk.bucket(element) ).getCellTopologyData() == NULL );



}

TEST(stk_topology_permutations, lexicographical_smallest_permutation_preserve_polarity)
{
    {
        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;
        unsigned shell_node_ids[3] = {10, 8, 12};
        {
            unsigned triangle_node_ids[3] = {12, 10, 8};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 2;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned triangle_node_ids[3] = {10, 8, 12};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 2;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned triangle_node_ids[3] = {8, 12, 10};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 2;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned triangle_node_ids[3] = {12, 8, 10};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 5;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned triangle_node_ids[3] = {10, 12, 8};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 5;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned triangle_node_ids[3] = {8, 10, 12};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity(triangle_node_ids, shell_node_ids);
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

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 0;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {4, 1, 2, 3};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 0;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {3, 4, 1, 2};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 0;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {2, 3, 4, 1};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 0;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {1, 4, 3, 2};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 4;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {2, 1, 4, 3};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 4;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {3, 2, 1, 4};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 4;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {4, 3, 2, 1};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 4;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
        {
            unsigned quad_node_ids[4] = {4, 2, 3, 1};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity(quad_node_ids, shell_node_ids);
            unsigned gold_lexicographical_smallest_permutation_index = 8;
            // driven by vertices, NOT mid-edge nodes
            EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
        }
    }
}

int check_permutation_given(stk::mesh::BulkData& mesh, stk::mesh::Entity elem, unsigned face_ord, stk::mesh::Permutation elem_perm, stk::mesh::Entity elem_face)
{
    stk::mesh::Entity face_nodes_buff[100], mapped_face_nodes_buff[100];

    stk::topology elem_topo = mesh.bucket(elem).topology();
    const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);

    elem_topo.face_nodes(elem_nodes, face_ord, face_nodes_buff);
    stk::topology face_topo = elem_topo.face_topology(face_ord);
    face_topo.permutation_nodes(face_nodes_buff, elem_perm, mapped_face_nodes_buff);

    // Another way to get to face's nodes.
    stk::mesh::Entity const *face_nodes = mesh.begin_nodes(elem_face);

    int innermost_hits = 0;
    for (unsigned ni = 0; ni < face_topo.num_nodes(); ++ni)
    {
        ++innermost_hits;
        EXPECT_EQ(face_nodes[ni], mapped_face_nodes_buff[ni]);
    }

    // Indeed, find_permutation computes what was stored!
    stk::mesh::Permutation perm = mesh.find_permutation(elem_topo, elem_nodes, face_topo, face_nodes, face_ord);
    EXPECT_EQ(perm, elem_perm);
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
    unsigned global_side_id = 1;
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    std::string name = "generated:1x1x1";
    stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, 1);

    mesh.modification_begin();
    stk::mesh::Part &quad4_part = mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
    stk::mesh::Entity side = stk::mesh::declare_element_side(mesh, global_side_id, elem, local_side_id, &quad4_part);
    mesh.modification_end();

    stk::mesh::Permutation elem_perm = static_cast<stk::mesh::Permutation>(0);
    check_permutation_given(mesh, elem, local_side_id, elem_perm, side);

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
    unsigned global_side_id = 1;

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
    stk::mesh::Entity side = mesh.declare_entity(stk::topology::FACE_RANK, global_side_id, parts);
    for(unsigned i = 0; i < nodes.size(); ++i) {
        unsigned ordinal = nodes.size() - i - 1;
        mesh.declare_relation(side, nodes[i], ordinal);
    }
    mesh.declare_relation(elem, side, local_side_id, perm);
    mesh.modification_end();

    check_permutation_given(mesh, elem, local_side_id, perm, side);

    const stk::mesh::Entity *side_nodes = mesh.begin_nodes(side);
    unsigned num_nodes = mesh.num_nodes(side);
    ASSERT_EQ(gold_num_nodes, num_nodes);

    for(unsigned i = 0; i < num_nodes; ++i)
    {
        unsigned ordinal = nodes.size() - i - 1;
        EXPECT_EQ(gold_side_ids[ordinal], mesh.identifier(side_nodes[i]));
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
        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(5);
        unsigned local_side_id = 0;
        test_side_creation_with_permutation(gold_side_ids[local_side_id],local_side_id, perm);
    }
}

TEST(stkTopologyFunctions, check_permutations_for_Hex_1x1x1)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        unsigned global_side_id = 1;
        unsigned gold_side_ids[4] = {5,6,8,7};
        unsigned local_side_id = 5; // nodes 2,4,8,5 are nodes of face 'local_side_id' of a 1x1x1 hex element (generated)

        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
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
        stk::mesh::Entity side = mesh.declare_entity(stk::topology::FACE_RANK, global_side_id, parts);
        for(unsigned i = 0; i < nodes.size(); ++i) {
            //unsigned ordinal = nodes.size() - i - 1;
            mesh.declare_relation(side, nodes[i], i);
        }
        mesh.declare_relation(elem, side, local_side_id, perm);
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
        unsigned global_side_id = 1;
        unsigned gold_side_ids[4] = {2,4,8,6};
        unsigned local_side_id = 1; // nodes 2,4,8,5 are nodes of face 'local_side_id' of a 1x1x1 hex element (generated)

        unsigned perm_value = 0;
        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(perm_value);

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
        stk::mesh::Entity side = mesh.declare_entity(stk::topology::FACE_RANK, global_side_id, parts);
        for(unsigned i = 0; i < nodes.size(); ++i) {
            //unsigned ordinal = nodes.size() - i - 1;
            mesh.declare_relation(side, nodes[i], i);
        }
        mesh.declare_relation(elem, side, local_side_id, perm);
        mesh.modification_end();

        bool rel_bad = false;

        std::vector<stk::mesh::Relation> recv_relations;
        pack_downward_relations_for_entity(mesh, side, recv_relations);
        unpack_not_owned_verify_compare_closure_relations(mesh, side, recv_relations, rel_bad);
        EXPECT_FALSE(rel_bad);

        std::vector<stk::mesh::Relation> relations;
        pack_downward_relations_for_entity(mesh, elem, relations);
        unpack_not_owned_verify_compare_closure_relations(mesh, elem, relations, rel_bad);
        EXPECT_FALSE(rel_bad);
    }
}

TEST(stkTopologyFunctions, check_permutation_consistency_parallel)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        unsigned global_side_id = 1;
        unsigned gold_side_ids[4] = {5,6,8,7};

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        std::string name = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        unsigned elem_id = 0;
        unsigned local_side_id = 0;
        unsigned perm_value = 0;

        if (mesh.parallel_rank()==0)
        {
            local_side_id=5;
            elem_id = 1;
            perm_value = 0;
        }
        else
        {
            local_side_id=4;
            elem_id = 2;
            perm_value = 4;
        }

        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(perm_value);

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
        stk::mesh::Entity side = mesh.declare_entity(stk::topology::FACE_RANK, global_side_id, parts);
        for(unsigned i = 0; i < nodes.size(); ++i) {
            mesh.declare_relation(side, nodes[i], i);
        }
        mesh.declare_relation(elem, side, local_side_id, perm);
        EXPECT_NO_THROW(mesh.modification_end());

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
        stk::mesh::MetaData meta(2);
        stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);

        mesh.modification_begin();
        stk::mesh::Entity Quad9 = mesh.declare_entity(stk::topology::ELEM_RANK, 1, meta.get_topology_root_part(stk::topology::QUAD_9_2D));
        unsigned node_ids[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        stk::mesh::EntityVector nodes(9);
        for(unsigned i=0;i<9;++i)
        {
            nodes[i] = mesh.declare_entity(stk::topology::NODE_RANK, node_ids[i]);
        }

        for(size_t i=0;i<nodes.size();++i)
        {
            mesh.declare_relation(Quad9, nodes[i], i);
        }

        mesh.modification_end();

        unsigned gold_side_node_ids[4][3] = {
                {1, 2, 5},
                {2, 3, 6},
                {3, 4, 7},
                {4, 1, 8}
        };

        mesh.modification_begin();

        unsigned global_side_id[] = {1, 2, 3, 4};
        stk::mesh::EntityVector sides(4);
        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(1);
        for(size_t i=0;i<sides.size();++i)
        {
            sides[i] = mesh.declare_entity(stk::topology::EDGE_RANK, global_side_id[i], meta.get_topology_root_part(stk::topology::LINE_3));
            mesh.declare_relation(Quad9, sides[i], i, perm);
            mesh.declare_relation(sides[i], nodes[gold_side_node_ids[i][0]-1], 1);
            mesh.declare_relation(sides[i], nodes[gold_side_node_ids[i][1]-1], 0);
            mesh.declare_relation(sides[i], nodes[gold_side_node_ids[i][2]-1], 2);
        }

        EXPECT_NO_THROW(mesh.modification_end());

        EXPECT_TRUE(mesh.check_permutation(Quad9, sides[0], 0, perm)) << "for side 0";
        EXPECT_TRUE(mesh.check_permutation(Quad9, sides[1], 1, perm)) << "for side 1";
        EXPECT_TRUE(mesh.check_permutation(Quad9, sides[2], 2, perm)) << "for side 2";
        EXPECT_TRUE(mesh.check_permutation(Quad9, sides[3], 3, perm)) << "for side 3";

        EXPECT_TRUE(stk::mesh::impl::check_permutations_on_all(mesh));
    }
}

std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> get_ordinal_and_permutation(stk::mesh::BulkData& mesh, stk::mesh::Entity element, stk::mesh::EntityRank to_rank, stk::mesh::EntityVector &nodes_of_sub_rank)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = std::make_pair(stk::mesh::INVALID_CONNECTIVITY_ORDINAL,
            stk::mesh::INVALID_PERMUTATION);

    const Entity* elemNodes = mesh.begin_nodes(element);
    stk::topology elemTopology = mesh.bucket(element).topology();
    unsigned num_entities_of_sub_topology = elemTopology.num_sub_topology(to_rank);
    unsigned max_nodes_possible = 100;
    stk::mesh::EntityVector nodes_of_sub_topology(max_nodes_possible);

    for (unsigned i=0;i<num_entities_of_sub_topology;++i)
    {
        stk::topology sub_topology = elemTopology.sub_topology(to_rank, i);
        unsigned num_nodes = sub_topology.num_nodes();
        ThrowRequireMsg(num_nodes<=max_nodes_possible, "Program error. Exceeded expected array dimensions. Contact sierra-help for support.");
        nodes_of_sub_topology.resize(num_nodes);
        elemTopology.sub_topology_nodes(elemNodes, to_rank, i, nodes_of_sub_topology.begin());
        std::pair<bool, unsigned> result = sub_topology.equivalent(nodes_of_sub_rank, nodes_of_sub_topology);
        if (result.first == true)
        {
            ordinalAndPermutation.first = static_cast<stk::mesh::ConnectivityOrdinal>(i);
            ordinalAndPermutation.second = static_cast<stk::mesh::Permutation>(result.second);
        }
    }

    return ordinalAndPermutation;
}

stk::mesh::Entity declare_element_to_sub_topology_with_nodes(stk::mesh::BulkData &mesh, stk::mesh::Entity elem, stk::mesh::EntityVector &side_nodes,
        stk::mesh::EntityId global_side_id, stk::mesh::EntityRank to_rank, stk::mesh::Part &part)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(mesh, elem, to_rank, side_nodes);

    stk::mesh::Entity side = mesh.declare_entity(to_rank, global_side_id, part);
    for (unsigned i=0;i<side_nodes.size();++i)
    {
        mesh.declare_relation(side, side_nodes[i], i);
    }

    mesh.declare_relation(elem, side, ordinalAndPermutation.first, ordinalAndPermutation.second);
    return side;
}

TEST(FEMHelper, get_ordinal_and_permutation)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        unsigned gold_side_node_ids[4] = {5,6,8,7};
        unsigned gold_num_nodes = 4;

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        std::string name = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        unsigned elem_id = 0;
        unsigned gold_local_side_id = 0;
        unsigned perm_value = 0;

        if (mesh.parallel_rank()==0)
        {
            gold_local_side_id=5;
            elem_id = 1;
            perm_value = 0;
        }
        else
        {
            gold_local_side_id=4;
            elem_id = 2;
            perm_value = 4;
        }

        stk::mesh::Permutation gold_permutation = static_cast<stk::mesh::Permutation>(perm_value);

        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);
        EXPECT_TRUE(mesh.bucket(elem).owned());

        stk::mesh::EntityVector side_nodes(gold_num_nodes);
        for(unsigned i = 0; i < gold_num_nodes; ++i)
        {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i]);
            side_nodes[i] = node;
        }

        stk::mesh::EntityRank to_rank = stk::topology::FACE_RANK;
        std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(mesh, elem, to_rank, side_nodes);

        ASSERT_TRUE(ordinalAndPermutation.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL);
        ASSERT_TRUE(ordinalAndPermutation.second != stk::mesh::INVALID_PERMUTATION);

        EXPECT_EQ(gold_local_side_id, ordinalAndPermutation.first);
        EXPECT_EQ(gold_permutation, ordinalAndPermutation.second);
    }
}

TEST(stkTopologyFunctions, check_permutation_consistency_using_FEMHelper_parallel)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::EntityId global_side_id = 1;
        unsigned gold_side_node_ids[4] = {5,6,8,7};
        unsigned gold_num_nodes = 4;

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        std::string name = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        unsigned elem_id = 0;

        if (mesh.parallel_rank()==0)
        {
            elem_id = 1;
        }
        else
        {
            elem_id = 2;
        }

        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);
        EXPECT_TRUE(mesh.bucket(elem).owned());

        stk::mesh::EntityVector side_nodes(gold_num_nodes);
        for(unsigned i = 0; i < gold_num_nodes; ++i) {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i]);
            side_nodes[i] = node;
        }

        stk::mesh::Part &part = mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
        mesh.modification_begin();
        stk::mesh::Entity side = declare_element_to_sub_topology_with_nodes(mesh, elem, side_nodes, global_side_id, stk::topology::FACE_RANK, part);
        EXPECT_NO_THROW(mesh.modification_end());

        std::vector<size_t> mesh_counts;
        stk::mesh::comm_mesh_counts(mesh, mesh_counts);
        size_t numFacesGlobal = 1u;
        EXPECT_EQ(numFacesGlobal, mesh_counts[stk::topology::FACE_RANK]);
        EXPECT_TRUE(mesh.is_valid(side));
    }
}

}

