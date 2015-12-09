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
   stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA, add_fmwk_data, &custom_connectivity);
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
