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

#include <stddef.h>                     // for size_t, NULL
#include <algorithm>                    // for sort
#include <stdexcept>                    // for logic_error
#include <stk_mesh/base/CreateEdges.hpp>  // for create_edges
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <gtest/gtest.h>
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_FALSE, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, get_connectivity
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequire

namespace {

using namespace stk::mesh;

fixtures::HexFixture* set_up_mesh(ConnectivityMap const& conn_map)
{
  const unsigned NX = 9;
  const unsigned NY = 9;
  const unsigned NZ = 9;

  fixtures::HexFixture* rv = new fixtures::HexFixture(MPI_COMM_WORLD, NX, NY, NZ, &conn_map);
  rv->m_meta.commit();
  rv->generate_mesh();

  stk::mesh::skin_mesh(rv->m_bulk_data);
  stk::mesh::create_edges(rv->m_bulk_data);

  return rv;
}

struct SortComparator
{
  bool operator()(std::pair<Entity, ConnectivityOrdinal> const& lhs, std::pair<Entity, ConnectivityOrdinal> const& rhs) const
  {
    return impl::HigherConnectivityCompare()(lhs.first, lhs.second, rhs.first, rhs.second);
  }
};

void check_equiv_conn(Bucket const& bucket_full_conn, Bucket const& bucket_min_conn, size_t ord, EntityRank rank)
{
  BulkData& mesh_min_conn  = bucket_min_conn.mesh();

  ThrowRequire(!mesh_min_conn.connectivity_map().valid(bucket_min_conn.entity_rank(), rank));

  EntityVector temp_entities;
  std::vector<ConnectivityOrdinal> temp_ordinals;
  Entity const* rel_entities_min = NULL;
  ConnectivityOrdinal const* rel_ordinals_min = NULL;
  size_t num_min_upward = get_connectivity(mesh_min_conn,
                                           bucket_min_conn[ord],
                                           rank,
                                           temp_entities,
                                           temp_ordinals);
  rel_entities_min = &*temp_entities.begin();
  rel_ordinals_min = &*temp_ordinals.begin();

  ASSERT_EQ(bucket_full_conn.num_connectivity(ord, rank), num_min_upward);

  Entity const* rel_entities_full              = bucket_full_conn.begin(ord, rank);
  ConnectivityOrdinal const* rel_ordinals_full = bucket_full_conn.begin_ordinals(ord, rank);

  // NOTE: computed back-connectivity may not be returned in the same order as it would
  //       be if it were stored, so we have to sort
  std::vector< std::pair<Entity, ConnectivityOrdinal> > temp;
  for (size_t i = 0; i < num_min_upward; ++i) {
    temp.push_back( std::make_pair( rel_entities_min[i], rel_ordinals_min[i] ) );
  }
  std::sort(temp.begin(), temp.end(), SortComparator());

  for (size_t i = 0; i < num_min_upward; ++i) {
    temp_entities[i] = temp[i].first;
    temp_ordinals[i] = temp[i].second;
  }

  for (size_t i = 0; i < num_min_upward; ++i) {
    EXPECT_EQ( rel_entities_min[i], rel_entities_full[i] );
    EXPECT_EQ( rel_ordinals_min[i], rel_ordinals_full[i] );
  }
}

TEST( UnitTestMinimalBackRelation, simpleHex )
{
  fixtures::HexFixture* fixture_with_full_conn = set_up_mesh(ConnectivityMap::classic_stk_mesh());
  fixtures::HexFixture* fixture_with_min_conn  = set_up_mesh(ConnectivityMap::minimal_upward_connectivity_map());

  BulkData& mesh_full_conn = fixture_with_full_conn->m_bulk_data;
  BulkData& mesh_min_conn  = fixture_with_min_conn->m_bulk_data;

  {
    for (EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEMENT_RANK; ++rank) {
      const BucketVector & buckets_full_conn = mesh_full_conn.buckets(rank);
      const BucketVector & buckets_min_conn  = mesh_min_conn.buckets(rank);
      ASSERT_EQ(buckets_full_conn.size(), buckets_min_conn.size());

      for (size_t ib=0, endb=buckets_full_conn.size(); ib < endb; ++ib) {
        const Bucket & bucket_full_conn = *buckets_full_conn[ib];
        const Bucket & bucket_min_conn  = *buckets_min_conn[ib];
        ASSERT_EQ(bucket_full_conn.size(), bucket_min_conn.size());

        for (size_t ord=0, end=bucket_full_conn.size(); ord<end; ++ord) {
          if ( rank > stk::topology::NODE_RANK ) {
            EXPECT_EQ(bucket_min_conn.num_elements(ord),0u); // no stored back-rels to elements except for nodes
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::ELEMENT_RANK);
          }
          if ( rank < stk::topology::FACE_RANK) {
            EXPECT_EQ(bucket_min_conn.num_faces(ord),0u);    // no stored back-rels to faces
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::FACE_RANK);
          }
          if ( rank < stk::topology::EDGE_RANK) {
            EXPECT_EQ(bucket_min_conn.num_edges(ord),0u);    // no stored back-rels to edges
            check_equiv_conn(bucket_full_conn, bucket_min_conn, ord, stk::topology::EDGE_RANK);
          }

          // Check that all downward rels are the same
          for (EntityRank irank = stk::topology::NODE_RANK; irank < rank; ++irank) {
            EXPECT_EQ(bucket_full_conn.num_connectivity(ord, irank), bucket_min_conn.num_connectivity(ord, irank));
          }
        }
      }
    }
  }

  delete fixture_with_full_conn;
  delete fixture_with_min_conn;
}

TEST( UnitTestNoUpwardConnectivity, simpleTri )
{
   const unsigned spatialDimension = 2;
   stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());
   stk::mesh::Part &tri3_part = metaData.declare_part_with_topology("tri3_part", stk::topology::TRI_3_2D);
   stk::mesh::Part &line2_part = metaData.declare_part_with_topology("line2_part", stk::topology::LINE_2);
   metaData.commit();

   stk::mesh::ConnectivityMap custom_connectivity = stk::mesh::ConnectivityMap::none();
   //Now set which connectivities we want enabled:
   custom_connectivity(stk::topology::ELEM_RANK, stk::topology::NODE_RANK) = stk::mesh::ConnectivityMap::fixed();
   custom_connectivity(stk::topology::ELEM_RANK, stk::topology::EDGE_RANK) = stk::mesh::ConnectivityMap::fixed();
   custom_connectivity(stk::topology::EDGE_RANK, stk::topology::NODE_RANK) = stk::mesh::ConnectivityMap::fixed();

   //Now verify that node->edge, node->elem and edge->elem connections are disabled, but
   //elem->node and edge->node connections are allowed:
   EXPECT_FALSE(custom_connectivity.valid(stk::topology::NODE_RANK, stk::topology::ELEM_RANK));
   EXPECT_FALSE(custom_connectivity.valid(stk::topology::NODE_RANK, stk::topology::EDGE_RANK));
   EXPECT_FALSE(custom_connectivity.valid(stk::topology::EDGE_RANK, stk::topology::ELEM_RANK));
   EXPECT_TRUE(custom_connectivity.valid(stk::topology::ELEM_RANK, stk::topology::NODE_RANK));
   EXPECT_TRUE(custom_connectivity.valid(stk::topology::EDGE_RANK, stk::topology::NODE_RANK));

#ifdef SIERRA_MIGRATION
   bool add_fmwk_data = false;
   stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD, add_fmwk_data, &custom_connectivity);
#else
   stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD, &custom_connectivity);
#endif
   if (mesh.parallel_size() > 1) {
     return;//this test can't run in parallel
   }

   mesh.modification_begin();

   //set up 1 element (3-node triangle) with elem->node and edge->node connections
   stk::mesh::PartVector elem_parts, side_parts;
   elem_parts.push_back(&tri3_part);
   side_parts.push_back(&line2_part);
   stk::mesh::EntityId elemId = 1;
   stk::mesh::EntityId elemNodeIds[] = {1, 2, 3};
   stk::mesh::EntityId elemEdgeIds[] = {6, 7, 8};
   stk::mesh::Entity elemNodes[3];
   stk::mesh::Entity elemEdges[3];
   stk::mesh::Entity elem = mesh.declare_entity(stk::topology::ELEM_RANK, elemId, elem_parts);
   elemNodes[0] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[0]);
   elemNodes[1] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[1]);
   elemNodes[2] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[2]);

   elemEdges[0] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[0], side_parts);
   elemEdges[1] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[1], side_parts);
   elemEdges[2] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[2], side_parts);

   //downward element -> node connectivity
   mesh.declare_relation(elem, elemNodes[0], 0);
   mesh.declare_relation(elem, elemNodes[1], 1);
   mesh.declare_relation(elem, elemNodes[2], 2);

   //downward edge -> node connectivity
   mesh.declare_relation(elemEdges[0], elemNodes[0], 0); mesh.declare_relation(elemEdges[0], elemNodes[1], 1);
   mesh.declare_relation(elemEdges[1], elemNodes[1], 0); mesh.declare_relation(elemEdges[1], elemNodes[2], 1);
   mesh.declare_relation(elemEdges[2], elemNodes[2], 0); mesh.declare_relation(elemEdges[2], elemNodes[0], 1);
   mesh.modification_end();

   //verify that get_connectivity throws (rather than infinitely recursing) since node->elem connections are disabled.
   stk::mesh::EntityVector connected_edges;
   EXPECT_THROW(stk::mesh::get_connectivity(mesh, elemNodes[0], stk::topology::EDGE_RANK, connected_edges), std::logic_error);
}

} // namespace
