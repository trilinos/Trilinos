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
void setup2Block2HexMesh(stk::mesh::BulkData& bulk)
{
//
//   proc 0      proc 1
//             |
//    block_1  |  block_2
//             |
//      8----7 |    7----12
//     /    /| |   /    / |
//    5----6 3 |  6----11 10
//    | 1  |/  |  | 2  | /
//    1----2   |  2----9
//             |
//             |
//             |
//
//shared nodes 2, 3, 6, 7
//

  if (bulk.parallel_size() > 2) {
    return;
  }

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::topology hex = stk::topology::HEX_8;
  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", hex);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", hex);
  meta.commit();

  bulk.modification_begin();

  stk::mesh::EntityId elem1_nodes[] = {1, 2, 3, 4, 5, 6, 7, 8};
  stk::mesh::EntityId elem2_nodes[] = {2, 9, 10, 3, 6, 11, 12, 7};

  stk::mesh::EntityId elemId = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::declare_element(bulk, block_1, elemId, elem1_nodes);
  }
  if (bulk.parallel_rank() == 1 || bulk.parallel_size() == 1) {
    elemId = 2;
    stk::mesh::declare_element(bulk, block_2, elemId, elem2_nodes);
  }
  if(bulk.parallel_rank() == 0 && bulk.parallel_size() == 2)
  {
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 2), 1);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 3), 1);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 6), 1);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 7), 1);
  }
  if(bulk.parallel_rank() == 1 && bulk.parallel_size() == 2)
  {
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 2), 0);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 3), 0);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 6), 0);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK , 7), 0);
  }

  bulk.modification_end();
}

