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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>                             // for AssertHelper, etc
#include <stddef.h>                                  // for size_t
#include <iostream>                                  // for operator<<, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/parallel/Parallel.hpp>            // for ParallelMachine
#include <string>                                    // for string, etc
#include <vector>                                    // for vector
#include "mpi.h"                                     // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Entity.hpp"                  // for Entity, etc
#include "stk_mesh/base/EntityKey.hpp"               // for operator<<, etc
#include "stk_mesh/base/GetEntities.hpp"             // for count_entities
#include "stk_mesh/base/MetaData.hpp"                // for MetaData, etc
#include "stk_mesh/base/Types.hpp"                   // for EntityVector, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_io/FillMesh.hpp"

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::MeshBuilder;
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
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MeshBuilder builder(pm);
  builder.set_spatial_dimension(spatial_dim);
  std::shared_ptr<BulkData> meshPtr = builder.create();
  BulkData& mesh = *meshPtr;
  MetaData& meta_data = mesh.mesh_meta_data();
  stk::mesh::Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Only run in 4-processor configuration
  if (p_size != 4) return;


  mesh.modification_begin();

  // Create the single element and four nodes for this processor.  Nodes
  // on the boundary with other processors will be shared.
  //
  Entity createdElem = mesh.declare_element(connInfo[p_rank][0], stk::mesh::ConstPartVector{&elem_part});

  // Create all 4 nodes for this element
  EntityVector createdNodes;
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][1], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][2], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][3], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][4], stk::mesh::ConstPartVector{&node_part}) );

  // Add relations to nodes
  mesh.declare_relation( createdElem, createdNodes[0], 0 );
  mesh.declare_relation( createdElem, createdNodes[1], 1 );
  mesh.declare_relation( createdElem, createdNodes[2], 2 );
  mesh.declare_relation( createdElem, createdNodes[3], 3 );

  // Declare all processors that share any of our nodes
  std::ostringstream oss;
  for (unsigned i = 0; i < numSharings; ++i)
  {
    Entity node = mesh.get_entity(stk::topology::NODE_RANK, sharingInfo[p_rank][i*2]);
    oss << "P" << p_rank << ": entity " << node << " = " <<  mesh.entity_key(node) << std::endl;
    int sharingProc = sharingInfo[p_rank][i*2+1];
    mesh.add_node_sharing(node, sharingProc);
  }
  std::cout << oss.str() << std::flush;

  mesh.modification_end();

  // Make sure we know about all nodes and elements (including aura, which
  // includes *all* entities in our small mesh)
  //
  std::vector<size_t> countsAll;
  count_entities(meta_data.universal_part(), mesh, countsAll);

  EXPECT_EQ( 9u, countsAll[stk::topology::NODE_RANK] );
  EXPECT_EQ( 4u, countsAll[stk::topology::ELEM_RANK] );

  // Count how many entities each proc owns, which will be different for each
  // proc (because of the lower parallel rank owning shared nodes on the
  // bondaries).  Check values in the processor-specific sections below.
  //
  std::vector<size_t> countsOwned;
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

TEST(UnitTestingOfBulkData, sharedProcsIntersection)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {
    return;
  }

  unsigned spatialDim = 3;
  MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x4", bulk);

  stk::mesh::Entity sharedNode9 = bulk.get_entity(stk::topology::NODE_RANK, 9);
  stk::mesh::Entity sharedNode10 = bulk.get_entity(stk::topology::NODE_RANK, 10);
  stk::mesh::EntityId unsharedNodeId = bulk.parallel_rank()==0 ? 1 : 13;
  stk::mesh::Entity unsharedNode = bulk.get_entity(stk::topology::NODE_RANK, unsharedNodeId);

  EXPECT_TRUE(bulk.bucket(sharedNode9).shared());
  EXPECT_TRUE(bulk.bucket(sharedNode10).shared());
  EXPECT_FALSE(bulk.bucket(unsharedNode).shared());
  EXPECT_FALSE(bulk.bucket(unsharedNode).in_aura());

  stk::mesh::EntityVector nodes;
  std::vector<int> sharingProcs;

  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_TRUE(sharingProcs.empty());

  nodes = {sharedNode9};
  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_EQ(1u, sharingProcs.size());
  const int otherProc = 1 - bulk.parallel_rank();
  EXPECT_EQ(otherProc, sharingProcs[0]);

  nodes = {unsharedNode};
  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_TRUE(sharingProcs.empty());

  nodes = {unsharedNode, sharedNode9};
  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_TRUE(sharingProcs.empty());

  nodes = {sharedNode9, unsharedNode};
  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_TRUE(sharingProcs.empty());

  nodes = {sharedNode9, sharedNode10};
  bulk.shared_procs_intersection(nodes, sharingProcs);
  EXPECT_EQ(1u, sharingProcs.size());
  EXPECT_EQ(otherProc, sharingProcs[0]);
}

TEST(UnitTestingOfBulkData, node_sharing_with_dangling_nodes)
{
  //
  // Test to make sure user-provided node sharing information is propagated
  // correctly.  This test manually builds up a 2D mesh of 2x2 quad elements,
  // and only runs on 4 processors (one element per processor).  This creates
  // the special case of a central node that has 4-way sharing between all
  // processors.
  //
  // Test is same as node_sharing, but creates dangling nodes emulating
  // what happens in Percept adaptivity - extra nodes 100 and 113 created
  //
  //
  //              O-----------O-----------O
  //              |10         |20         |31
  //              |           |           |
  //          p0  |    100    |    101    |  p1
  //              |           |           |
  //              |           |           |
  //              O---100-----O----113----O
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

  const unsigned numSharings1 = 1u;
  EntityId sharingInfo1[][10] = { {100, 2} ,    // p0
                                  {113, 3} ,    // p1
                                  {100, 0} ,    // p2
                                  {113, 1}  };  // p3

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MeshBuilder builder(pm);
  builder.set_spatial_dimension(spatial_dim);
  std::shared_ptr<BulkData> meshPtr = builder.create();
  BulkData& mesh = *meshPtr;
  MetaData& meta_data = mesh.mesh_meta_data();
  stk::mesh::Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
  meta_data.commit();

  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  // Only run in 4-processor configuration
  if (p_size != 4) return;


  mesh.modification_begin();

  // Create the single element and four nodes for this processor.  Nodes
  // on the boundary with other processors will be shared.
  //
  Entity createdElem = mesh.declare_element(connInfo[p_rank][0], stk::mesh::ConstPartVector{&elem_part});

  // Create all 4 nodes for this element
  EntityVector createdNodes;
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][1], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][2], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][3], stk::mesh::ConstPartVector{&node_part}) );
  createdNodes.push_back( mesh.declare_node(connInfo[p_rank][4], stk::mesh::ConstPartVector{&node_part}) );


  bool doNew = true;
  if (doNew)
  {
    mesh.declare_node(sharingInfo1[p_rank][0], stk::mesh::ConstPartVector{&node_part}) ;
    std::cout << "P[" << p_rank << "] declared node= " << sharingInfo1[p_rank][0] << std::endl;
  }

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

  if (doNew)
  {
    for (unsigned i = 0; i < numSharings1; ++i)
    {
      int sharingProc = sharingInfo1[p_rank][i*2+1];
      Entity node = mesh.get_entity(stk::topology::NODE_RANK, sharingInfo1[p_rank][i*2]);
      std::cout << "P[" << p_rank << "] add_node_sharing node= " << sharingInfo1[p_rank][i*2] << " sharingProc= " << sharingProc <<  " is_valid= " << mesh.is_valid(node) << std::endl;
      mesh.add_node_sharing(node, sharingProc);
    }
  }

  EXPECT_NO_THROW(mesh.modification_end());

  // Make sure we know about all nodes and elements (including aura, which
  // includes *all* entities in our small mesh)
  //
  std::vector<size_t> countsAll;
  count_entities(meta_data.universal_part(), mesh, countsAll);

  EXPECT_EQ( 10u, countsAll[stk::topology::NODE_RANK] ); //
  EXPECT_EQ( 4u,  countsAll[stk::topology::ELEM_RANK] );

  // Count how many entities each proc owns, which will be different for each
  // proc (because of the lower parallel rank owning shared nodes on the
  // bondaries).  Check values in the processor-specific sections below.
  //
  std::vector<size_t> countsOwned;
  count_entities(meta_data.locally_owned_part(), mesh, countsOwned);

  std::vector<int> sharingProcs;

  if (p_rank == 0)
  {
    EXPECT_EQ( 5u, countsOwned[stk::topology::NODE_RANK] ); //

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
    EXPECT_EQ( 3u, countsOwned[stk::topology::NODE_RANK] ); //

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

void move_elem(stk::mesh::BulkData& bulk, stk::mesh::EntityId id,
               stk::mesh::Part& addPart, stk::mesh::Part& rmPart)
{
  bulk.modification_begin();
  stk::mesh::PartVector addParts = {&addPart};
  stk::mesh::PartVector rmParts = {&rmPart};

  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, id);
  if (bulk.is_valid(elem) && bulk.bucket(elem).owned()) {
    bulk.change_entity_parts(elem, addParts, rmParts);
  }
  bulk.modification_end();
}

void create_faces(stk::mesh::BulkData& bulk)
{
  bulk.modification_begin();

  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
    stk::mesh::Entity elem3 = bulk.get_entity(stk::topology::ELEM_RANK, 3);
    ASSERT_TRUE(bulk.is_valid(elem1) && bulk.is_valid(elem2) && bulk.is_valid(elem3));
    stk::mesh::ConnectivityOrdinal sideOrdinal = 5;
    bulk.declare_element_side(elem1, sideOrdinal, stk::mesh::ConstPartVector{});
    bulk.declare_element_side(elem2, sideOrdinal, stk::mesh::ConstPartVector{});
    bulk.declare_element_side(elem3, sideOrdinal, stk::mesh::ConstPartVector{});
  }

  bulk.modification_end();
}

void detach_face(stk::mesh::BulkData& bulk, stk::mesh::EntityId elemId)
{
  bulk.modification_begin();
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
  if (bulk.is_valid(elem) && bulk.bucket(elem).owned()) {
    unsigned numFaces = bulk.num_faces(elem);
    ASSERT_EQ(1u, numFaces);
    const stk::mesh::Entity* faces = bulk.begin_faces(elem);
    stk::mesh::Entity face = faces[0];
    const stk::mesh::ConnectivityOrdinal* ords = bulk.begin_face_ordinals(elem);
    stk::mesh::ConnectivityOrdinal ord = ords[0];
    bulk.destroy_relation(elem, face, ord);
  }
  bulk.modification_end();
}

void assert_shared_faces(const stk::mesh::BulkData& bulk, size_t expectedNumSharedFaces)
{
  stk::mesh::Selector shared = bulk.mesh_meta_data().globally_shared_part();
  ASSERT_TRUE(expectedNumSharedFaces == stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, shared));
}

TEST(UnitTestBulkData, resolveSharedModifiedFaceNextToUnmodifiedFaces)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_spatial_dimension(3);
  std::shared_ptr<BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::io::fill_mesh("generated:3x1x2", bulk);
  stk::mesh::Part& block1 = *meta.get_part("block_1");

  move_elem(bulk, 2, block2, block1);

  create_faces(bulk);
  assert_shared_faces(bulk, 3);

  detach_face(bulk, 2);
  assert_shared_faces(bulk, 3);
}

