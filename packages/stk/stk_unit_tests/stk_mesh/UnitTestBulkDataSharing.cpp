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

#include <gtest/gtest.h>
#include <vector>                         // for vector
#include <stk_mesh/base/BulkData.hpp>     // for BulkData
#include <stk_util/parallel/Parallel.hpp> // for ParallelMachine
#include "stk_mesh/base/Selector.hpp"     // for Selector
#include "stk_mesh/base/Entity.hpp"       // for Entity
#include "stk_mesh/base/MetaData.hpp"     // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"        // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"      // for topology, etc
#include "stk_mesh/base/GetEntities.hpp"  // for count_entities

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;

TEST(UnitTestingOfBulkData, node_sharing)
{
  //
  // Test to make sure user-provided node sharing information is propagated
  // correctly.  This test manually builds up a 2D mesh of 2x2 quad elements,
  // and only runs on 4 processors (one element per processor).  This creates
  // the special case of a central node that has 4-way sharing between all
  // processors.
  //
  //              O-----------O-----------O
  //              |10         |20         |31
  //              |           |           |
  //          p0  |    100    |    101    |  p1
  //              |           |           |
  //              |           |           |
  //              O-----------O-----------O
  //              |40         |50         |61
  //              |           |           |
  //          p2  |    102    |    103    |  p3
  //              |           |           |
  //              |           |           |
  //              O-----------O-----------O
  //               72          82          93
  //
  // Numbering Scheme:
  //  Element ID: 103                      Node ID: 61
  //                ^                               ^^
  //                |                               ||
  //                \--- Owning proc (3)            |\--- Owning proc (1)
  //                                                \---- Node number (6)
  //
  // Connectivity information for each processor:
  //   {elem, node0, node1, node2, node3}
  EntityId connInfo[][5] = { {100, 10, 40, 50, 20},    // p0
                             {101, 20, 50, 61, 31},    // p1
                             {102, 40, 72, 82, 50},    // p2
                             {103, 50, 82, 93, 61} };  // p3

  // Sharing information for each created node for each processor.  Data
  // is ordered in pairs as:
  //   {node, sharingProc, node, sharingProc, ...}
  const unsigned numSharings = 5u;
  EntityId sharingInfo[][10] = { {20, 1, 40, 2, 50, 1, 50, 2, 50, 3},    // p0
                                 {20, 0, 61, 3, 50, 0, 50, 2, 50, 3},    // p1
                                 {40, 0, 82, 3, 50, 0, 50, 1, 50, 3},    // p2
                                 {61, 1, 82, 2, 50, 0, 50, 1, 50, 2} };  // p3

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  MetaData meta_data(spatial_dim, entity_rank_names);
  stk::mesh::Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Only run in 4-processor configuration
  if (p_size != 4) return;


  mesh.modification_begin();

  // Create the single element and four nodes for this processor.  Nodes
  // on the boundary with other processors will be shared.
  //
  Entity createdElem = mesh.declare_entity(stk::topology::ELEM_RANK, connInfo[p_rank][0], elem_part);

  // Create all 4 nodes for this element
  EntityVector createdNodes;
  createdNodes.push_back( mesh.declare_entity(stk::topology::NODE_RANK, connInfo[p_rank][1], node_part) );
  createdNodes.push_back( mesh.declare_entity(stk::topology::NODE_RANK, connInfo[p_rank][2], node_part) );
  createdNodes.push_back( mesh.declare_entity(stk::topology::NODE_RANK, connInfo[p_rank][3], node_part) );
  createdNodes.push_back( mesh.declare_entity(stk::topology::NODE_RANK, connInfo[p_rank][4], node_part) );

  // Add relations to nodes
  mesh.declare_relation( createdElem, createdNodes[0], 0 );
  mesh.declare_relation( createdElem, createdNodes[1], 1 );
  mesh.declare_relation( createdElem, createdNodes[2], 2 );
  mesh.declare_relation( createdElem, createdNodes[3], 3 );

  // Declare all processors that share any of our nodes
  for (unsigned i = 0; i < numSharings; ++i)
  {
    Entity node = mesh.get_entity(stk::topology::NODE_RANK, sharingInfo[p_rank][i*2]);
    int sharingProc = sharingInfo[p_rank][i*2+1];
    mesh.add_node_sharing(node, sharingProc);
  }

  mesh.modification_end();


  // Make sure we know about all nodes and elements (including aura, which
  // includes *all* entities in our small mesh)
  //
  std::vector<unsigned> countsAll;
  count_entities(meta_data.universal_part(), mesh, countsAll);

  EXPECT_EQ( 9u, countsAll[stk::topology::NODE_RANK] );
  EXPECT_EQ( 4u, countsAll[stk::topology::ELEM_RANK] );

  // Count how many entities each proc owns, which will be different for each
  // proc (because of the lower parallel rank owning shared nodes on the
  // bondaries).  Check values in the processor-specific sections below.
  //
  std::vector<unsigned> countsOwned;
  count_entities(meta_data.locally_owned_part(), mesh, countsOwned);

  std::vector<int> sharingProcs;

  if (p_rank == 0)
  {
    EXPECT_EQ( 4u, countsOwned[stk::topology::NODE_RANK] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[0]), sharingProcs);
    EXPECT_TRUE( sharingProcs.empty() );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[1]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 2, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[2]), sharingProcs);
    EXPECT_EQ( 3u, sharingProcs.size() );
    EXPECT_EQ( 1, sharingProcs[0] );
    EXPECT_EQ( 2, sharingProcs[1] );
    EXPECT_EQ( 3, sharingProcs[2] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[3]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 1, sharingProcs[0] );
  }
  else if (p_rank == 1)
  {
    EXPECT_EQ( 2u, countsOwned[stk::topology::NODE_RANK] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[0]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 0, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[1]), sharingProcs);
    EXPECT_EQ( 3u, sharingProcs.size() );
    EXPECT_EQ( 0, sharingProcs[0] );
    EXPECT_EQ( 2, sharingProcs[1] );
    EXPECT_EQ( 3, sharingProcs[2] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[2]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 3, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[3]), sharingProcs);
    EXPECT_TRUE( sharingProcs.empty() );
  }
  else if (p_rank == 2)
  {
    EXPECT_EQ( 2u, countsOwned[stk::topology::NODE_RANK] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[0]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 0, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[1]), sharingProcs);
    EXPECT_TRUE( sharingProcs.empty() );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[2]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 3, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[3]), sharingProcs);
    EXPECT_EQ( 3u, sharingProcs.size() );
    EXPECT_EQ( 0, sharingProcs[0] );
    EXPECT_EQ( 1, sharingProcs[1] );
    EXPECT_EQ( 3, sharingProcs[2] );
  }
  else if (p_rank == 3)
  {
    EXPECT_EQ( 1u, countsOwned[stk::topology::NODE_RANK] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[0]), sharingProcs);
    EXPECT_EQ( 3u, sharingProcs.size());
    EXPECT_EQ( 0, sharingProcs[0] );
    EXPECT_EQ( 1, sharingProcs[1] );
    EXPECT_EQ( 2, sharingProcs[2] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[1]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 2, sharingProcs[0] );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[2]), sharingProcs);
    EXPECT_TRUE( sharingProcs.empty() );

    mesh.comm_shared_procs(mesh.entity_key(createdNodes[3]), sharingProcs);
    EXPECT_EQ( 1u, sharingProcs.size() );
    EXPECT_EQ( 1, sharingProcs[0] );
  }


}

