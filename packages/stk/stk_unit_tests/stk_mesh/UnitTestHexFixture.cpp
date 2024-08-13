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
#include <map>                          // for map, map<>::mapped_type
#include <ostream>                      // for basic_ostream::operator<<
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD

#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc

namespace {

using stk::mesh::MetaData;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;

TEST( UnitTestHexFixture, elem_ids_1d_x )
{
  // Test 3x1x1 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(1,0,0), 2u );
  EXPECT_EQ( hf.elem_id(2,0,0), 3u );
}

TEST( UnitTestHexFixture, elem_ids_3d_x )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(1,0,0), 2u );
  EXPECT_EQ( hf.elem_id(2,0,0), 3u );
}

TEST( UnitTestHexFixture, elem_ids_1d_y )
{
  // Test 1x3x1 HexFixture structure
  const unsigned NX = 1;
  const unsigned NY = 3;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(0,1,0), 2u );
  EXPECT_EQ( hf.elem_id(0,2,0), 3u );
}

TEST( UnitTestHexFixture, elem_ids_3d_y )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(0,1,0), 4u );
  EXPECT_EQ( hf.elem_id(0,2,0), 7u );
}

TEST( UnitTestHexFixture, elem_ids_1d_z )
{
  // Test 1x1x3 HexFixture structure
  const unsigned NX = 1;
  const unsigned NY = 1;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(0,0,1), 2u );
  EXPECT_EQ( hf.elem_id(0,0,2), 3u );
}

TEST( UnitTestHexFixture, elem_ids_3d_z )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(0,0,1), 10u );
  EXPECT_EQ( hf.elem_id(0,0,2), 19u );
}

TEST( UnitTestHexFixture, elem_ids_3d_diag )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.generate_mesh();
  EXPECT_EQ( hf.elem_id(0,0,0), 1u );
  EXPECT_EQ( hf.elem_id(1,1,1), 14u );
  EXPECT_EQ( hf.elem_id(2,2,2), 27u );
}

TEST( UnitTestHexFixture, trivial_parallel_2 )
{
  // Test a customized element distribution with one element on proc 0 and 1
  // and none on the other procs

  const int p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 2;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. One element will go on rank 0, the other on rank 1, any
  // other ranks get nothing.
  std::map<int,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    parallel_distribution[0] = element_ids;
    element_ids[0] = 2;
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  if (p_rank <= 1) {
    hf.fill_node_map(parallel_distribution);
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify element_id 1 is owned by proc 0
  // Verify element_id 2 is owned by proc 1
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  Entity entity_1 = mesh.get_entity(element_rank, 1);
  Entity entity_2 = mesh.get_entity(element_rank ,2);
  if (p_rank <= 1) {
    ASSERT_TRUE( mesh.is_valid(entity_1) );
    ASSERT_TRUE( mesh.is_valid(entity_2) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_2) );
  }
  else {
    EXPECT_TRUE( !mesh.is_valid(entity_1) );
    EXPECT_TRUE( !mesh.is_valid(entity_2) );
  }
}

TEST( UnitTestHexFixture, disjoint_parallel_psizex1x1 )
{
  // Test a customized element distribution with one element on each proc

  const int p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  const unsigned NX = p_size;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign each processor an element such that rank+1 = elem_id
  std::map<int,std::vector<EntityId> > parallel_distribution;
  for (int p=0 ; p < p_size ; ++p) {
    std::vector< EntityId> element_ids;
    element_ids.push_back(p+1); // element id's start at 1
    parallel_distribution[p] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  hf.fill_node_map(parallel_distribution);
  hf.generate_mesh(parallel_distribution[p_rank]);
  stk::mesh::BulkData & mesh = hf.m_bulk_data;
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;

  // We should always know about, and own, the element assigned to us
  Entity my_entity    = mesh.get_entity(element_rank, p_rank + 1);
  ASSERT_TRUE( mesh.is_valid(my_entity) );
  EXPECT_EQ( p_rank, mesh.parallel_owner_rank(my_entity) );

  // If applicable, we know about the element on adjacent lower rank
  if (p_rank > 0) {
    Entity prior_entity = mesh.get_entity(element_rank, p_rank);
    ASSERT_TRUE( mesh.is_valid(prior_entity) );
    EXPECT_EQ( p_rank - 1, mesh.parallel_owner_rank(prior_entity) );
  }

  // If applicable, we know about the element on adjacent higher rank
  if (p_rank < p_size - 1) {
    Entity next_entity   = mesh.get_entity(element_rank, p_rank + 2);
    ASSERT_TRUE( mesh.is_valid(next_entity) );
    EXPECT_EQ( p_rank + 1, mesh.parallel_owner_rank(next_entity) );
  }
}

TEST( UnitTestHexFixture, disjoint_parallel_4x2x1 )
{
  // layout: 4x2x1 hex mesh
  // elements:
  // [ e_1, e_2, e_3, e_4 ]
  // [ e_5, e_6, e_7, e_8 ]
  // processors:
  // [ p_0, p_1, p_1, p_1 ]
  // [ p_1, p_0, p_1, p_1 ]
  //

  const int p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 4;
  const unsigned NY = 2;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,6 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<int,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(6);
    parallel_distribution[0] = element_ids; // proc 0
    element_ids.clear();
    element_ids.push_back(2);
    element_ids.push_back(3);
    element_ids.push_back(4);
    element_ids.push_back(5);
    element_ids.push_back(7);
    element_ids.push_back(8);
    parallel_distribution[1] = element_ids; // proc 1
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  if (p_rank <= 1) {
    hf.fill_node_map(parallel_distribution);
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  Entity entity_1 = mesh.get_entity(element_rank, 1);
  Entity entity_2 = mesh.get_entity(element_rank, 2);
  Entity entity_3 = mesh.get_entity(element_rank, 3);
  Entity entity_4 = mesh.get_entity(element_rank, 4);
  Entity entity_5 = mesh.get_entity(element_rank, 5);
  Entity entity_6 = mesh.get_entity(element_rank, 6);
  Entity entity_7 = mesh.get_entity(element_rank, 7);
  Entity entity_8 = mesh.get_entity(element_rank, 8);
  if (p_rank == 0) {
    ASSERT_TRUE( mesh.is_valid(entity_1) );
    ASSERT_TRUE( mesh.is_valid(entity_2) );
    ASSERT_TRUE( mesh.is_valid(entity_3) );
    ASSERT_TRUE( !mesh.is_valid(entity_4) );
    ASSERT_TRUE( mesh.is_valid(entity_5) );
    ASSERT_TRUE( mesh.is_valid(entity_6) );
    ASSERT_TRUE( mesh.is_valid(entity_7) );
    ASSERT_TRUE( !mesh.is_valid(entity_8) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_2) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_3) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_5) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_6) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_7) );
  }
  else if (p_rank == 1) {
    ASSERT_TRUE( mesh.is_valid(entity_1) );
    ASSERT_TRUE( mesh.is_valid(entity_2) );
    ASSERT_TRUE( mesh.is_valid(entity_3) );
    ASSERT_TRUE( mesh.is_valid(entity_4) );
    ASSERT_TRUE( mesh.is_valid(entity_5) );
    ASSERT_TRUE( mesh.is_valid(entity_6) );
    ASSERT_TRUE( mesh.is_valid(entity_7) );
    ASSERT_TRUE( mesh.is_valid(entity_8) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_2) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_3) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_4) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_5) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_6) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_7) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_8) );
  }
  else {
    EXPECT_TRUE( !mesh.is_valid(entity_1) );
    EXPECT_TRUE( !mesh.is_valid(entity_2) );
    EXPECT_TRUE( !mesh.is_valid(entity_3) );
    EXPECT_TRUE( !mesh.is_valid(entity_4) );
    EXPECT_TRUE( !mesh.is_valid(entity_5) );
    EXPECT_TRUE( !mesh.is_valid(entity_6) );
    EXPECT_TRUE( !mesh.is_valid(entity_7) );
    EXPECT_TRUE( !mesh.is_valid(entity_8) );
  }
}

TEST( UnitTestHexFixture, disjoint_parallel_5x1x1 )
{
  // layout:
  // [ e_1, e_2, e_3, e_4, e_5 ] elements
  // [ p_0, p_1, p_1, p_1, p_0 ] processors
  //
  const int p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 5;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,5 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<int,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(5);
    parallel_distribution[0] = element_ids;
    element_ids.clear();
    element_ids.push_back(2);
    element_ids.push_back(3);
    element_ids.push_back(4);
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_meta.commit();
  if (p_rank <= 1) {
    hf.fill_node_map(parallel_distribution);
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  Entity entity_1 = mesh.get_entity(element_rank, 1);
  Entity entity_2 = mesh.get_entity(element_rank, 2);
  Entity entity_3 = mesh.get_entity(element_rank, 3);
  Entity entity_4 = mesh.get_entity(element_rank, 4);
  Entity entity_5 = mesh.get_entity(element_rank, 5);
  if (p_rank == 0) {
    ASSERT_TRUE( mesh.is_valid(entity_1) );
    ASSERT_TRUE( mesh.is_valid(entity_2) );
    ASSERT_TRUE( !mesh.is_valid(entity_3) );
    ASSERT_TRUE( mesh.is_valid(entity_4) );
    ASSERT_TRUE( mesh.is_valid(entity_5) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_2) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_4) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_5) );
  }
  else if (p_rank == 1) {
    ASSERT_TRUE( mesh.is_valid(entity_1) );
    ASSERT_TRUE( mesh.is_valid(entity_2) );
    ASSERT_TRUE( mesh.is_valid(entity_3) );
    ASSERT_TRUE( mesh.is_valid(entity_4) );
    ASSERT_TRUE( mesh.is_valid(entity_5) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_1) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_2) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_3) );
    EXPECT_EQ( 1, mesh.parallel_owner_rank(entity_4) );
    EXPECT_EQ( 0, mesh.parallel_owner_rank(entity_5) );
  }
  else {
    EXPECT_TRUE( !mesh.is_valid(entity_1) );
    EXPECT_TRUE( !mesh.is_valid(entity_2) );
    EXPECT_TRUE( !mesh.is_valid(entity_3) );
    EXPECT_TRUE( !mesh.is_valid(entity_4) );
    EXPECT_TRUE( !mesh.is_valid(entity_5) );
  }
}

} // end namespace
