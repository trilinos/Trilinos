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

#include "mpi.h"                        // for MPI_Barrier, MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp"  // for RingFixture
#include <gtest/gtest.h>                // for ASSERT_TRUE, AssertHelper, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <vector>                       // for vector
namespace stk { namespace mesh { class Selector; } }

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::fixtures::RingFixture;

//----------------------------------------------------------------------

TEST(UnitTestingOfBulkData, testChangeParts_ringmesh)
{
  // This unit test tests part operations and verifies operations
  // by looking at bucket supersets.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const unsigned nPerProc   = 10;
  const int p_rank     = stk::parallel_machine_rank( pm );
  const int p_size     = stk::parallel_machine_size( pm );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;

  // Create the ring mesh

  RingFixture ring_mesh( pm , nPerProc , true /* generate parts */, stk::mesh::BulkData::NO_AUTO_AURA);
  ring_mesh.m_meta_data.commit();
  BulkData& bulk = ring_mesh.m_bulk_data;

  bulk.modification_begin();
  ring_mesh.generate_mesh( );
  ASSERT_TRUE(bulk.modification_end());

  ring_mesh.fixup_node_ownership();

  Part & part_owns = ring_mesh.m_meta_data.locally_owned_part();
  Part & part_univ = ring_mesh.m_meta_data.universal_part();

  // Check that local elements are in the expected parts. Note that the
  // RingMesh puts each element in its own part.
  for ( unsigned i = 0 ; i < nLocalElement ; ++i ) {
    const unsigned n = i + nPerProc * p_rank ;
    Entity const element = bulk.get_entity( stk::topology::ELEMENT_RANK /*entity rank*/,
                                            ring_mesh.m_element_ids[n] );
    ASSERT_TRUE( bulk.is_valid(element) );
    ASSERT_TRUE( bulk.bucket(element).member( part_univ ) );
    ASSERT_TRUE( bulk.bucket(element).member( part_owns ) );
    ASSERT_TRUE( bulk.bucket(element).member( * ring_mesh.m_element_parts[ n % ring_mesh.m_element_parts.size() ] ) );
  }

  // Check that local nodes are in the expected parts. Note that the relations
  // that nodes have to elements should cause induced membership of the node
  // in the parts of both elements it touches.
  for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
    const unsigned n = ( i + nPerProc * p_rank ) % ring_mesh.m_node_ids.size();
    const unsigned e0 = n ;
    const unsigned e1 = ( n + ring_mesh.m_element_ids.size() - 1 ) % ring_mesh.m_element_ids.size();
    const unsigned ns = ring_mesh.m_element_parts.size();
    const unsigned n0 = e0 % ns ;
    const unsigned n1 = e1 % ns ;
    Part * const epart_0 = ring_mesh.m_element_parts[ n0 < n1 ? n0 : n1 ];
    Part * const epart_1 = ring_mesh.m_element_parts[ n0 < n1 ? n1 : n0 ];

    Entity const node = bulk.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[n] );
    ASSERT_TRUE( bulk.is_valid(node) );
    if ( bulk.parallel_owner_rank(node) == p_rank ) {
      ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
      ASSERT_TRUE( bulk.bucket(node).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(node).member( *epart_0 ) );
      ASSERT_TRUE( bulk.bucket(node).member( *epart_1 ) );
    }
    else {
      ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
      ASSERT_TRUE( ! bulk.bucket(node).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(node).member( * epart_0 ) );
      ASSERT_TRUE( bulk.bucket(node).member( * epart_1 ) );
    }
  }

  bulk.modification_begin();

  // On rank 0, change all locally owned elements to the extra-part then check
  // for correct part membership
  if ( 0 == p_rank ) {
    for ( unsigned i = 0 ; i < nLocalElement ; ++i ) {
      const unsigned n = i + nPerProc * p_rank ;

      PartVector add(1); add[0] = & ring_mesh.m_element_part_extra ;
      PartVector rem(1); rem[0] = ring_mesh.m_element_parts[ n % ring_mesh.m_element_parts.size() ];

      Entity const element = bulk.get_entity( stk::topology::ELEMENT_RANK , ring_mesh.m_element_ids[n] );
      bulk.change_entity_parts( element , add , rem );
      ASSERT_TRUE( bulk.bucket(element).member( part_univ ) );
      ASSERT_TRUE( bulk.bucket(element).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(element).member(ring_mesh.m_element_part_extra ) );
    }
  }

  bulk.modification_end();

  // Modification end has been called, check that the part changes made
  // in the previous step are reflected across the other procs.
  for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
    const unsigned n = ( i + nPerProc * p_rank ) % ring_mesh.m_node_ids.size();
    const unsigned e0 = n ;
    const unsigned e1 = ( n + ring_mesh.m_element_ids.size() - 1 ) % ring_mesh.m_element_ids.size();
    const unsigned ns = ring_mesh.m_element_parts.size();
    const unsigned n0 = e0 % ns ;
    const unsigned n1 = e1 % ns ;
    Part * ep_0 = e0 < nLocalElement ? & ring_mesh.m_element_part_extra : ring_mesh.m_element_parts[n0] ;
    Part * ep_1 = e1 < nLocalElement ? & ring_mesh.m_element_part_extra : ring_mesh.m_element_parts[n1] ;

    Part * epart_0 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_0 : ep_1 ;
    Part * epart_1 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_1 : ep_0 ;

    Entity const node = bulk.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[n] );
    ASSERT_TRUE( bulk.is_valid(node) );
    if ( bulk.parallel_owner_rank(node) == p_rank ) {
      ASSERT_TRUE( bulk.bucket(node).member( part_owns ) );
    }
    else {
      ASSERT_TRUE( ! bulk.bucket(node).member( part_owns ) );
    }

    ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
    ASSERT_TRUE( bulk.bucket(node).member( *epart_0 ) );
    ASSERT_TRUE( bulk.bucket(node).member( *epart_1 ) );
  }
}
