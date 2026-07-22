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

#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector, vector<>::iterator

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp"  // for RingFixture
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace mesh { class Ghosting; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }

using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityVector;
using stk::mesh::Part;
using stk::mesh::Relation;
using stk::mesh::Selector;
using stk::mesh::EntityId;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Ghosting;
using stk::mesh::fixtures::BoxFixture;
using stk::mesh::fixtures::HexFixture;
using stk::mesh::fixtures::RingFixture;
using stk::unit_test_util::build_mesh;

namespace {

TEST(UnitTestingOfRelation, testRelationCoverage)
{
  // Test some relation error cases

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier ( MPI_COMM_WORLD );

  // Just use any random fixture for convenience
  HexFixture fixture(pm, 3 /*x*/, 3 /*y*/, 3 /*z*/);
  MetaData& meta  = fixture.m_meta;
  BulkData& bulk  = fixture.m_bulk_data;

  meta.commit();

  fixture.generate_mesh();

  Entity node0 = (*bulk.buckets(stk::topology::NODE_RANK)[0])[0];

  bulk.modification_begin();
  std::vector<Part*> empty_parts;
  stk::mesh::EntityId  new_id = bulk.parallel_rank() + 1;
  Entity edge = bulk.declare_edge(new_id, empty_parts);

  // Cannot declare a back relation
  ASSERT_THROW ( bulk.declare_relation ( node0 , edge , 0 /*rel ord*/) , std::runtime_error );
}

TEST(UnitTestingOfRelation, testRelationNoGhosting)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  const unsigned nPerProc = 10;

  RingFixture mesh( pm , nPerProc , false /* No element parts */, stk::mesh::BulkData::NO_AUTO_AURA );
  mesh.m_meta_data.commit();

  BulkData& ring_bulk = mesh.m_bulk_data;

  ring_bulk.modification_begin();
  mesh.generate_mesh( );
  ASSERT_TRUE(ring_bulk.modification_end());

  mesh.fixup_node_ownership();

  // This process' first element in the loop
  // if a parallel mesh has a shared node

  Entity elementnew = ring_bulk.get_entity( stk::topology::ELEMENT_RANK , mesh.m_element_ids[ nPerProc * p_rank ] );

  ring_bulk.modification_begin();
  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    if ( p != p_rank ) {
      ASSERT_EQ( ring_bulk.in_shared( ring_bulk.entity_key(elementnew) , p ), false );
      ASSERT_EQ( ring_bulk.in_send_ghost( ring_bulk.entity_key(elementnew) , p ), false );
    }
  }

  elementnew = ring_bulk.get_entity( stk::topology::ELEMENT_RANK , mesh.m_element_ids[ nPerProc * p_rank ] );
  ASSERT_EQ( ring_bulk.in_send_ghost( ring_bulk.entity_key(elementnew) , p_rank+100 ), false );

  Entity node = ring_bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nPerProc * p_rank ] );
  ASSERT_EQ( ring_bulk.in_shared( ring_bulk.entity_key(node) , p_rank+100 ), false );
}

TEST(UnitTestingOfRelation, testRelationWithGhosting)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  const unsigned nPerProc = 10;

  if ( 1 < p_size ) { // With ghosting
    RingFixture mesh( pm , nPerProc , false /* No element parts */ );
    mesh.m_meta_data.commit();
    BulkData& bulk = mesh.m_bulk_data;

    bulk.modification_begin();
    mesh.generate_mesh();
    ASSERT_TRUE(bulk.modification_end());

    mesh.fixup_node_ownership();

    const unsigned nNotOwned = nPerProc * p_rank ;

    // The not-owned shared entity:
    Entity node = bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nNotOwned ] );

    bulk.modification_begin();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      if ( p != p_rank ) {
        ASSERT_EQ( bulk.in_send_ghost( bulk.entity_key(node) , p ), false );
      }
    }

    //not owned and not shared
    Entity node2 = bulk.get_entity( stk::topology::NODE_RANK , mesh.m_node_ids[ nPerProc * p_rank ] );

    ASSERT_EQ( bulk.in_shared( bulk.entity_key(node2) , p_rank+100 ), false );
    ASSERT_EQ( bulk.in_send_ghost( bulk.entity_key(node) , p_rank+100 ), false );
  }
}

TEST(UnitTestingOfRelation, testDegenerateRelation)
{
  // Test that, if you set up degenerate relations, only of the relations
  // is deleted when you destroy one of the degenerate relations.
  // BulkData::destroy_relation has been changed to take a relation-id so
  // that it can work this way.
  //
  // To test this, we set up an element that has several relations
  // to the same node and then delete them one by one.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dim, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;
  meta_data.commit();
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create element
  Entity elem = mesh.declare_element(p_rank+1 /*elem_id*/, empty_parts);

  // Create node
  Entity node = mesh.declare_node(p_rank+1 /*node_id*/, empty_parts);

  // Add degenerate relations
  const unsigned nodes_per_elem = 4;
  for (unsigned i = 0; i < nodes_per_elem; ++i) {
    mesh.declare_relation( elem, node, i );
  }

  // Elem should have nodes-per-elem relations
  ASSERT_EQ( nodes_per_elem, static_cast<unsigned>(mesh.num_nodes(elem)) );

  // Destroy relation one-by-one, always checking that appropriate number
  // of relations remain.
  for (unsigned i = 0; i < nodes_per_elem; ++i) {
    mesh.destroy_relation( elem, node, i );
    ASSERT_EQ( nodes_per_elem - (i+1), static_cast<unsigned>(mesh.num_nodes(elem)) );
  }
}

TEST(UnitTestingOfRelation, testRelationExtraRank)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned spatial_dim = 3;
  const unsigned num_ext_ranks = 1;
  const char * ext_rank_names[num_ext_ranks] = { "EXT_RANK_0" };
  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  for (unsigned i = 0; i < num_ext_ranks; ++i) {
    entity_rank_names.push_back(ext_rank_names[i]);
  }

  stk::mesh::MeshBuilder builder(pm);
  builder.set_spatial_dimension(spatial_dim);
  builder.set_entity_rank_names(entity_rank_names);
  std::shared_ptr<BulkData> bulkPtr = builder.create();

  MetaData& meta_data = bulkPtr->mesh_meta_data();
  Part &hex8_part = meta_data.declare_part_with_topology("Hex8_Part", stk::topology::HEX_8);
  meta_data.commit();
  BulkData& mesh = *bulkPtr;

  std::vector<EntityRank> ext_ranks;
  for (unsigned i = 0;  i < num_ext_ranks; ++i) {
    ext_ranks.push_back(meta_data.entity_rank(ext_rank_names[i]));
  }

  const unsigned p_rank = mesh.parallel_rank();
  stk::mesh::PartVector empty_parts;
  unsigned new_ent_id = (p_rank << 16) + 1;

  // baseline
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&hex8_part);
  mesh.modification_begin();
  Entity elem = mesh.declare_element(new_ent_id++, elem_parts);
  Entity node;
  for (unsigned i = 0; i < 8; ++i)
  {
    node = mesh.declare_node(new_ent_id++, empty_parts);
    mesh.declare_relation(elem, node, i);
  }
  mesh.modification_end();

  ASSERT_FALSE(mesh.num_connectivity(elem, ext_ranks[0]));

  stk::mesh::ConnectivityOrdinal ord_count = static_cast<stk::mesh::ConnectivityOrdinal>(0);

  // simplest test for get_others_offset_range(..) and other_entities_have_single_rank(..)
  mesh.modification_begin();
  Entity ere_0_0 = mesh.declare_entity(ext_ranks[0], new_ent_id++, empty_parts);
  Entity ere_0_1 = mesh.declare_entity(ext_ranks[0], new_ent_id++, empty_parts);
  mesh.declare_relation(ere_0_0, node, ++ord_count);
  mesh.declare_relation(ere_0_1, node, ++ord_count);
  mesh.declare_relation(ere_0_0, elem, ++ord_count);
  mesh.declare_relation(ere_0_1, elem, ++ord_count);
  mesh.modification_end();
  {
    EXPECT_EQ(mesh.num_connectivity(node, ext_ranks[0]), 2u);
    EXPECT_EQ(mesh.num_connectivity(elem, ext_ranks[0]), 2u);
  }
}

TEST(UnitTestingOfRelation, testDoubleDeclareOfRelation)
{
  // It should be legal to declare the same relation between shared
  // entities on two procs.
  //
  // 1---3---5
  // | 1 | 2 |
  // 2---4---6
  //
  // To test this, we use the mesh above, with elem 1 going on rank 0 and
  // elem 2 going on rank 1. Nodes 3,4 are shared along with the edge between
  // nodes 3 and 4. On both procs we declare relations from the shared edge
  // to the shared nodes on both procs.
  //
  // TODO: If we change how declare_relation works, not requiring all
  // sharers to declare the same relations, but instead allowing just
  // the owner to declare relations, that should be tested here.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dim, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;
  Part &quad4_part = meta_data.declare_part_with_topology("quad4_part", stk::topology::QUAD_4_2D);
  Part &line2_part = meta_data.declare_part_with_topology("LINE2_part", stk::topology::LINE_2);
  meta_data.commit();
  unsigned p_rank = mesh.parallel_rank();
  unsigned p_size = mesh.parallel_size();

  if (p_size != 2) {
    return;
  }

  Entity edge = Entity();
  EntityVector nodes;
  const unsigned nodes_per_elem = 4, nodes_per_side = 2;

  // Begin modification cycle so we can create the entities and relations
  mesh.modification_begin();

  // entities with EDGE_RANK, FACE_RANK, or ELEMENT_RANK need topology before
  // modification_end() is called.
  stk::mesh::PartVector elems_parts, sides_parts, empty_parts;
  elems_parts.push_back(&quad4_part);
  sides_parts.push_back(&line2_part);

  // Create element
  Entity elem = mesh.declare_element(p_rank+1 /*elem_id*/, elems_parts);

  // Create nodes
  const unsigned starting_node_id = p_rank * nodes_per_side + 1;
  for (unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id) {
    nodes.push_back(mesh.declare_node(id, empty_parts));
  }

  // Add relations to nodes
  unsigned rel_id = 0;
  for (EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id) {
    mesh.declare_relation( elem, *itr, rel_id );
  }

  stk::mesh::EntityVector side_nodes(2);
  side_nodes[0] = mesh.get_entity(stk::topology::NODE_RANK, 3);
  side_nodes[1] = mesh.get_entity(stk::topology::NODE_RANK, 4);
  int otherProc = (1 - p_rank);
  mesh.add_node_sharing(side_nodes[0], otherProc);
  mesh.add_node_sharing(side_nodes[1], otherProc);

  // Create side (edge)
  unsigned local_side_id = 2;
  if (p_rank == 1)
  {
    local_side_id = 0;
  }
  edge = mesh.declare_element_side(elem, local_side_id, sides_parts);

  stk::topology elem_top = mesh.bucket(elem).topology();
  stk::mesh::Permutation perm1 = stk::mesh::find_permutation(mesh,
      elem_top, nodes.data(), elem_top.side_topology(local_side_id), side_nodes.data(), local_side_id);
  ASSERT_TRUE(perm1 != stk::mesh::Permutation::INVALID_PERMUTATION);

  // Set up duplicate relations from edge to nodes
  rel_id = 0;
  const unsigned starting_node_idx = (1 - p_rank) * nodes_per_side;
  for (unsigned node_idx = starting_node_idx;
       node_idx < starting_node_idx + nodes_per_side;
       ++node_idx, ++rel_id) {
    mesh.declare_relation( edge, nodes[node_idx], rel_id );
  }
  mesh.modification_end();

}

}
