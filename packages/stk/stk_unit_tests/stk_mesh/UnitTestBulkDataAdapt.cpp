// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_mesh/fixtures/FixtureNodeSharing.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string, char_traits
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { struct EntityKey; } }


using stk::mesh::Part;
using stk::mesh::Bucket;
using stk::mesh::PairIterRelation;
using stk::mesh::PairIterEntityComm;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::BaseEntityRank;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::BucketVector;
using stk::mesh::fixtures::BoxFixture;

namespace {
const EntityRank NODE_RANK = stk::topology::NODE_RANK;
} // empty namespace


void printBuckets(std::ostringstream& msg, BulkData& mesh)
{
  const BucketVector & buckets = mesh.buckets(NODE_RANK);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      const Bucket& bucket = *buckets[i];
      msg << " bucket[" << i << "] = ";
      size_t bucket_size = bucket.size();
      for (unsigned ie=0; ie < bucket_size; ie++)
        {
          msg << mesh.identifier(bucket[ie]) << ", ";
        }
    }
}

static void checkBuckets( BulkData& mesh)
{
  const BucketVector & buckets = mesh.buckets(NODE_RANK);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      Bucket* bucket = buckets[i];
      ASSERT_TRUE(bucket->assert_correct());
    }
}

TEST(UnitTestingOfBulkData, test_other_ghosting_2)
{
  //
  // testing if modification flags propagate properly for ghosted entities
  //
  // To test this, we focus on a single node shared on 2 procs, ghosted on others
  //

  /**
   * 1D Mesh (node,owner)--[elem,owner]---(...)
   *
   * <---(70,0)--[500,0]--(41,1)--[301,1]---(42,2)---[402,2]---(70,0)--->
   *
   * <---(50,0)--[100,0]--(21,1)--[201,1]---(32,2)---[302,2]---(50,0)--->
   *
   */

  // elem, node0, node1, owner
  EntityId elems_0[][4] = { {100, 21, 50, 0}, {201, 21, 32, 1}, {302, 32, 50, 2},
                            {500, 41, 70, 0}, {301, 41, 42, 1}, {402, 42, 70, 2}  };
  // node, owner
  EntityId nodes_0[][2] = { {21,1}, {50,0}, {32, 2}, {41, 1}, {42, 1}, {70, 0} };

  unsigned nelems = sizeof(elems_0)/4/sizeof(EntityId);
  unsigned nnodes = sizeof(nodes_0)/2/sizeof(EntityId);

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 1;

  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  entity_rank_names.push_back("FAMILY_TREE");

  MetaData meta_data(spatial_dim, entity_rank_names);
 Part & elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::LINE_2_1D);
 Part & node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);


  meta_data.commit();
  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if (p_size != 3) return;

  // Build map for node sharing
  stk::mesh::fixtures::NodeToProcsMMap nodes_to_procs;
  {
    for (unsigned ielem=0; ielem < nelems; ielem++) {
      int e_owner = static_cast<int>(elems_0[ielem][3]);
      stk::mesh::fixtures::AddToNodeProcsMMap(nodes_to_procs, elems_0[ielem][2], e_owner);
      stk::mesh::fixtures::AddToNodeProcsMMap(nodes_to_procs, elems_0[ielem][1], e_owner);
    }
  }

  //
  // Begin modification cycle so we can create the entities and relations
  //

  // Create elements
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  Entity elem = Entity();

  mesh.modification_begin();

  for (unsigned ielem=0; ielem < nelems; ielem++)
    {
      if (static_cast<int>(elems_0[ielem][3]) == p_rank)
        {
          elem = mesh.declare_entity(elem_rank, elems_0[ielem][0], elem_part);

          EntityVector nodes;
          // Create node on all procs
          nodes.push_back( mesh.declare_entity(NODE_RANK, elems_0[ielem][2], node_part) );
          nodes.push_back( mesh.declare_entity(NODE_RANK, elems_0[ielem][1], node_part) );

          // Add relations to nodes
          mesh.declare_relation( elem, nodes[0], 0 );
          mesh.declare_relation( elem, nodes[1], 1 );

          // Node sharing
          stk::mesh::fixtures::DoAddNodeSharings(mesh, nodes_to_procs, mesh.identifier(nodes[0]), nodes[0]);
          stk::mesh::fixtures::DoAddNodeSharings(mesh, nodes_to_procs, mesh.identifier(nodes[1]), nodes[1]);
        }
    }

  mesh.modification_end();

  Entity node1 = Entity();

  // change node owners
  mesh.modification_begin();

  std::vector<EntityProc> change;

  for (unsigned inode=0; inode < nnodes; inode++)
    {
      node1 = mesh.get_entity(stk::topology::NODE_RANK, nodes_0[inode][0]);
      if (mesh.is_valid(node1) && mesh.parallel_owner_rank(node1) == p_rank)
        {
          int dest = nodes_0[inode][1];
          EntityProc eproc(node1, dest);
          change.push_back(eproc);
        }
    }

  mesh.change_entity_owner( change );

  mesh.modification_end();

  checkBuckets(mesh);

  MPI_Barrier(MPI_COMM_WORLD);


  // attempt to delete a node and its elems but on a ghosted proc
  mesh.modification_begin();

  if (p_rank == 2)
    {
      node1 = mesh.get_entity(stk::topology::NODE_RANK, 21);
      Entity elem1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 201);
      Entity elem2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 100);

      bool did_it_elem = mesh.destroy_entity(elem1);
      did_it_elem = did_it_elem & mesh.destroy_entity(elem2);
      ASSERT_TRUE(did_it_elem);
      bool did_it = mesh.destroy_entity(node1);
      ASSERT_TRUE(did_it);
    }

  mesh.modification_end();

  checkBuckets(mesh);

  // this node should no longer exist anywhere
  node1 = mesh.get_entity(stk::topology::NODE_RANK, 21);

  // uncomment to force failure of test
  // ASSERT_TRUE(node1 == 0);

}

