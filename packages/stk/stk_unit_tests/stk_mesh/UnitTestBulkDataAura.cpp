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

#include <gtest/gtest.h>                // for ASSERT_TRUE, AssertHelper, etc
#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for string, allocator, etc
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_Barrier, MPI_COMM_WORLD, etc

#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp"

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::MeshBuilder;
using stk::mesh::PartVector;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityVector;

TEST(UnitTestingOfBulkData, aura1DRing_RestoreDeletedAuraEntity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

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

  MeshBuilder builder(pm);
  builder.set_spatial_dimension(spatial_dim);
  builder.set_entity_rank_names(entity_rank_names);
  std::shared_ptr<BulkData> bulkPtr = builder.create();
  BulkData& mesh = *bulkPtr;
  MetaData& meta_data = mesh.mesh_meta_data();
  Part & elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::LINE_2_1D);
  Part & node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);

  meta_data.commit();

  int p_rank = mesh.parallel_rank();

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
  Entity elem = Entity();

  mesh.modification_begin();

  for (unsigned ielem=0; ielem < nelems; ielem++)
  {
    if (static_cast<int>(elems_0[ielem][3]) == p_rank)
    {
      elem = mesh.declare_element(elems_0[ielem][0], stk::mesh::ConstPartVector{&elem_part});

      EntityVector nodes;
      // Create node on all procs
      nodes.push_back( mesh.declare_node(elems_0[ielem][2], stk::mesh::ConstPartVector{&node_part}) );
      nodes.push_back( mesh.declare_node(elems_0[ielem][1], stk::mesh::ConstPartVector{&node_part}) );

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

  // The node that we deleted on proc 2 was an aura node. so modification_end
  // should have restored the node (and also restored the aura elements we deleted).
  // The way this mesh is arranged, every proc has every element.
  node1 = mesh.get_entity(stk::topology::NODE_RANK, 21);
  ASSERT_TRUE(mesh.is_valid(node1));
  elem = mesh.get_entity(stk::topology::ELEM_RANK, 201);
  ASSERT_TRUE(mesh.is_valid(elem));
  elem = mesh.get_entity(stk::topology::ELEM_RANK, 100);
  ASSERT_TRUE(mesh.is_valid(elem));
}

