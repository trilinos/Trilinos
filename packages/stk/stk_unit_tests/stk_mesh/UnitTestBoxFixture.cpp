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

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator|
#include <gtest/gtest.h>
#include <vector>                       // for vector, etc

#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp"  // for BoxFixture::BOX, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_util/util/PairIter.hpp"   // for PairIter

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::BoxFixture;

namespace {

const EntityRank NODE_RANK = stk::topology::NODE_RANK;

class ExposePartition : public BoxFixture
{
public:
  static void expose_box_partition( int ip , int up , int axis ,
                                    const BOX box ,
                                    BOX p_box[] )
  {
    box_partition(ip, up, axis, box, p_box);
  }
};

}

TEST( UnitTestBoxFixture, verifyBoxFixture )
{
  // A unit test to verify the correctness of the BoxFixture fixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create the box fixture we'll be testing

  // box specifications
  const BoxFixture::BOX root_box = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  BoxFixture::BOX local_box = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  BoxFixture fixture(pm);
  MetaData& meta = fixture.fem_meta();
  stk::unit_test_util::BulkDataTester& bulk = fixture.bulk_data();

  const EntityRank element_rank = stk::topology::ELEMENT_RANK;

  const int p_rank = bulk.parallel_rank();
  const int p_size = bulk.parallel_size();

  BoxFixture::BOX * const p_box = new BoxFixture::BOX[ p_size ];
  ExposePartition::expose_box_partition( 0, p_size, 2, root_box, &p_box[0] );

  meta.commit();

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );

  const unsigned nx = local_box[0][1] - local_box[0][0] ;
  const unsigned ny = local_box[1][1] - local_box[1][0] ;
  const unsigned nz = local_box[2][1] - local_box[2][0] ;

  const unsigned e_local = nx * ny * nz ;
  const unsigned n_local = ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 );

  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;

  std::vector<size_t> local_count ;

  // Verify that the correct entities are on this process

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
    for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
      for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {

        const EntityId n0= 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
        const EntityId n1= 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
        const EntityId n2= 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
        const EntityId n3= 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
        const EntityId n4= 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
        const EntityId n5= 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
        const EntityId n6= 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
        const EntityId n7= 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

        const EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;

        std::vector<Entity> nodes(8);
        nodes[0] = bulk.get_entity( NODE_RANK, n0 );
        nodes[1] = bulk.get_entity( NODE_RANK, n1 );
        nodes[2] = bulk.get_entity( NODE_RANK, n2 );
        nodes[3] = bulk.get_entity( NODE_RANK, n3 );
        nodes[4] = bulk.get_entity( NODE_RANK, n4 );
        nodes[5] = bulk.get_entity( NODE_RANK, n5 );
        nodes[6] = bulk.get_entity( NODE_RANK, n6 );
        nodes[7] = bulk.get_entity( NODE_RANK, n7 );

        Entity elem = bulk.get_entity( element_rank, elem_id );

        std::vector<Entity> elems ;
        stk::mesh::get_entities_through_relations(bulk, nodes, element_rank, elems);
        ASSERT_EQ( elems.size() , size_t(1) );
        ASSERT_TRUE( elems[0] == elem );
      }
    }
  }

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used = select_owned |
      meta.globally_shared_part();
  Selector select_all( meta.universal_part() );

  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( e_local , local_count[3] );
  ASSERT_EQ( 0u , local_count[2] );
  ASSERT_EQ( 0u , local_count[1] );
  ASSERT_EQ( n_local , local_count[0] );

  ASSERT_TRUE(bulk.modification_end());

  // Verify declarations and sharing

  stk::mesh::count_entities( select_used , bulk , local_count );
  ASSERT_EQ( local_count[3] , e_local );
  ASSERT_EQ( local_count[2] , 0u );
  ASSERT_EQ( local_count[1] , 0u );
  ASSERT_EQ( local_count[0] , n_local );

  for ( int k = local_box[2][0] ; k <= local_box[2][1] ; ++k ) {
    for ( int j = local_box[1][0] ; j <= local_box[1][1] ; ++j ) {
      for ( int i = local_box[0][0] ; i <= local_box[0][1] ; ++i ) {
        EntityRank node_type = stk::topology::NODE_RANK;
        EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
        Entity const node = bulk.get_entity( node_type , node_id );
        ASSERT_TRUE( bulk.is_valid(node) );
        // Shared if on a processor boundary.
        const bool shared =
            ( k == local_box[2][0] && k != root_box[2][0] ) ||
            ( k == local_box[2][1] && k != root_box[2][1] ) ||
            ( j == local_box[1][0] && j != root_box[1][0] ) ||
            ( j == local_box[1][1] && j != root_box[1][1] ) ||
            ( i == local_box[0][0] && i != root_box[0][0] ) ||
            ( i == local_box[0][1] && i != root_box[0][1] );
        if (bulk.parallel_size() > 1) {
          std::vector<int> shared_procs;
          bulk.comm_shared_procs(bulk.entity_key(node),shared_procs);
          ASSERT_EQ( shared , ! shared_procs.empty() );
        }
      }
    }
  }

  size_t count_shared_node_pairs = 0 ;
  for ( int p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
    for ( int k = p_box[p][2][0] ; k <= p_box[p][2][1] ; ++k )
      if ( local_box[2][0] <= k && k <= local_box[2][1] ) {

        for ( int j = p_box[p][1][0] ; j <= p_box[p][1][1] ; ++j )
          if ( local_box[1][0] <= j && j <= local_box[1][1] ) {

            for ( int i = p_box[p][0][0] ; i <= p_box[p][0][1] ; ++i )
              if ( local_box[0][0] <= i && i <= local_box[0][1] ) {

                EntityRank node_type = stk::topology::NODE_RANK;
                EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
                Entity const node = bulk.get_entity( node_type , node_id );
                ASSERT_TRUE( bulk.is_valid(node) );
                // Must be shared with 'p'
                std::vector<int> shared_procs;
                bulk.comm_shared_procs(bulk.entity_key(node),shared_procs);
                std::vector<int>::const_iterator it=std::find(shared_procs.begin(),shared_procs.end(),p);
                ASSERT_TRUE( it != shared_procs.end() );

                ++count_shared_node_pairs ;
              }
          }
      }
  }

  size_t count_shared_entities = 0 ;
  for (stk::mesh::EntityCommListInfoVector::const_iterator
       i = bulk.my_internal_comm_list().begin() ;
       i != bulk.my_internal_comm_list().end() ;
       ++i) {
    std::vector<int> shared_procs;
    bulk.comm_shared_procs(i->key,shared_procs);
    count_shared_entities += shared_procs.size();
  }
  ASSERT_EQ( count_shared_entities , count_shared_node_pairs );

  delete [] p_box;
}
