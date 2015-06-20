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
#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc

inline
void setup8Quad4ProcMesh2D(stk::mesh::BulkData& bulk)
{
    ThrowRequire(bulk.parallel_size() == 4);
//
//     id/proc
//     11/0--12/0--13/1--14/2--15/3
//       |     |     |     |     |
//       | 5/0 | 6/1 | 7/2 | 8/3 |
//       |     |     |     |     |
//      6/0---7/0---8/1---9/2--10/3
//       |     |     |     |     |
//       | 1/0 | 2/1 | 3/2 | 4/3 |
//       |     |     |     |     |
//      1/0---2/0---3/1---4/2---5/3

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  meta.commit();

  bulk.modification_begin();

  //2 elem-ids for each proc
  stk::mesh::EntityId proc_elemIDs[][2] = {{1, 5}, {2, 6}, {3, 7}, {4, 8}};

  //list of node-ids for each element
  stk::mesh::EntityIdVector elem_nodeIDs[] {
    {1, 2, 7, 6},
    {2, 3, 8, 7},
    {3, 4, 9, 8},
    {4, 5, 10, 9},
    {6, 7, 12, 11},
    {7, 8, 13, 12},
    {8, 9, 14, 13},
    {9, 10, 15, 14}
  };

  //list of triples: (owner-proc, shared-nodeID, sharing-proc)
  int shared_nodeIDs_and_procs[][3] = {
  {0, 2, 1}, {0, 7, 1}, {0, 12, 1},                                   //proc 0
  {1, 2, 0}, {1, 7, 0}, {1, 12, 0}, {1, 3, 2}, {1, 8, 2}, {1, 13, 2}, //proc 1
  {2, 3, 1}, {2, 8, 1}, {2, 13, 1}, {2, 4, 3}, {2, 9, 3}, {2, 14, 3}, //proc 2
  {3, 4, 2}, {3, 9, 2}, {3, 14, 2}                                    //proc 3
  };

  stk::mesh::EntityId elemID = proc_elemIDs[bulk.parallel_rank()][0];
  int elemIdx = elemID - 1;
  stk::mesh::declare_element(bulk, block_1, elemID, elem_nodeIDs[elemIdx]);

  elemID = proc_elemIDs[bulk.parallel_rank()][1];
  elemIdx = elemID - 1;
  stk::mesh::declare_element(bulk, block_1, elemID, elem_nodeIDs[elemIdx]);

  int numSharedNodeTriples = 18;

  for(int i=0; i<numSharedNodeTriples; ++i) {
     if (bulk.parallel_rank() == shared_nodeIDs_and_procs[i][0]) {
         int nodeID = shared_nodeIDs_and_procs[i][1];
         int sharingProc = shared_nodeIDs_and_procs[i][2];
         stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeID);
         bulk.add_node_sharing(node, sharingProc);
     }
  }

  bulk.modification_end();

  std::vector<unsigned> counts(meta.entity_rank_count());
  stk::mesh::Selector owned_or_shared = meta.locally_owned_part() | meta.globally_shared_part();
  stk::mesh::count_entities(owned_or_shared, bulk, counts);

  unsigned numNodes = counts[stk::topology::NODE_RANK];
  unsigned numElems = counts[stk::topology::ELEM_RANK];

  unsigned expectedNumNodes = 6;
  unsigned expectedNumElems = 2;
  EXPECT_EQ(expectedNumNodes, numNodes);
  EXPECT_EQ(expectedNumElems, numElems);
}

