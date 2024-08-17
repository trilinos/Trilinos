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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <map>                          // for map, map<>::mapped_type
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CreateEdges.hpp>  // for create_edges
#include <stk_topology/topology.hpp>    // for topology, etc
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD

#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityRank
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;

TEST( UnitTestGloballyShared, DISABLED_keyhole_3x1 )
{
  // layout:
  // [ e_1, e_2, e_3 ] elements
  // [ p_0, p_1, p_0 ] processors
  //
  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 3;
  const unsigned NY = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,5 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(3);
    parallel_distribution[0] = element_ids;
    element_ids.clear();
    element_ids.push_back(2);
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::QuadFixture qf(MPI_COMM_WORLD,NX,NY);
  qf.m_meta.commit();
  if (p_rank <= 1) {
    qf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    qf.generate_mesh( empty_vector ) ;
  }

  BulkData & mesh = qf.m_bulk_data;

  stk::mesh::create_edges(mesh);

  // Quad edge ordinals:
  //            2
  //          -----
  //         |     |
  //       3 |     | 1
  //         |     |
  //          -----
  //            0

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  Entity element_1 = mesh.get_entity(element_rank, 1);
  Entity element_2 = mesh.get_entity(element_rank, 2);
  Entity element_3 = mesh.get_entity(element_rank, 3);
  if (p_rank == 0) {
    ASSERT_TRUE( mesh.is_valid(element_1) );
    ASSERT_TRUE( mesh.is_valid(element_2) );
    ASSERT_TRUE( mesh.is_valid(element_3) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(element_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(element_2) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(element_3) );
    // Verify global sharing of edges on element_1 and element_3
    // element_1:  edge_1 should be globally shared
    // element_3:  edge_3 should be globally shared
    stk::mesh::Entity const* element_1_edge_relations = mesh.begin_edges(element_1);
    const int num_element_1_edges = mesh.num_edges(element_1);
    ASSERT_EQ( 4, num_element_1_edges );
    Entity element_1_edge_1 = element_1_edge_relations[1];
    ASSERT_TRUE( mesh.is_valid(element_1_edge_1) );
    EXPECT_TRUE( mesh.in_shared(mesh.entity_key(element_1_edge_1)) );

    stk::mesh::Entity const* element_3_edge_relations = mesh.begin_edges(element_3);
    const int num_element_3_edges = mesh.num_edges(element_3);
    ASSERT_EQ( 4, num_element_3_edges );
    Entity element_3_edge_3 = element_3_edge_relations[3];
    ASSERT_TRUE( mesh.is_valid(element_3_edge_3) );
    EXPECT_TRUE( mesh.in_shared(mesh.entity_key(element_3_edge_3)) );

    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_1)) );
    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_2)) );
    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_3)) );
  }
  else if (p_rank == 1) {
    ASSERT_TRUE( mesh.is_valid(element_1) );
    ASSERT_TRUE( mesh.is_valid(element_2) );
    ASSERT_TRUE( mesh.is_valid(element_3) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(element_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(element_2) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(element_3) );
    // Verify global sharing of edges on element_2
    // element_2:  edge_0 and edge_2 should _not be_ globally shared
    //             edge_1 and edge_3 should _be_ globally shared
    stk::mesh::Entity const* element_2_edge_relations = mesh.begin_edges(element_2);
    const int num_element_2_edges = mesh.num_edges(element_2);
    ASSERT_EQ( 4, num_element_2_edges );

    Entity element_2_edge_0 = element_2_edge_relations[0];
    ASSERT_TRUE( mesh.is_valid(element_2_edge_0) );

    Entity element_2_edge_2 = element_2_edge_relations[2];
    ASSERT_TRUE( mesh.is_valid(element_2_edge_2) );

    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_2_edge_0)) );
    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_2_edge_2)) );

    Entity element_2_edge_1 = element_2_edge_relations[1];
    ASSERT_TRUE( mesh.is_valid(element_2_edge_1) );
    Entity element_2_edge_3 = element_2_edge_relations[3];
    ASSERT_TRUE( mesh.is_valid(element_2_edge_3) );

    EXPECT_TRUE( mesh.in_shared(mesh.entity_key(element_2_edge_1)) );
    EXPECT_TRUE( mesh.in_shared(mesh.entity_key(element_2_edge_3)) );

    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_1)) );
    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_2)) );
    EXPECT_FALSE( mesh.in_shared(mesh.entity_key(element_3)) );
  }
  else {
    EXPECT_TRUE( !mesh.is_valid(element_1) );
    EXPECT_TRUE( !mesh.is_valid(element_2) );
    EXPECT_TRUE( !mesh.is_valid(element_3) );
  }
}

