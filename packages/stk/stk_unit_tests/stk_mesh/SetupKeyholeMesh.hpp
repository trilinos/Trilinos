#ifndef setupkeyholemeshhpp
#define setupkeyholemeshhpp
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
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
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
#include <stk_unit_test_utils/ioUtils.hpp>


namespace stk { namespace mesh { namespace unit_test {

inline
void setupKeyholeMesh2D_case1(stk::mesh::BulkData& bulk)
{
  //
  //   proc 0      proc 1
  //            |
  //            |  block_2 block_3
  //            |
  //  block_1   |  10---9  9----12
  //            |  | 3  |  |  4  |
  //    4----3  |  3----8  8----11
  //    | 1  |  |
  //    1----2  |  2----7
  //            |  | 2  |
  //            |  5----6
  //            |
  //
  //shared nodes 2 and 3 should be members of block_1 and block_2 on both procs
  //nodes 8 and 9 are ghosts on proc 0, and should be members of block_2 and block_3
  //
  //if edges are added, the edge between nodes 2 and 3 should be a member of block_1 not block_2.
  //
  //also, the edge between nodes 8 and 9 should be a member of block_2 and block_3 on both procs.

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::QUAD_4_2D);
  meta.commit();

  bulk.modification_begin();

  stk::mesh::EntityIdVector elem1_nodes {1, 2, 3, 4};
  stk::mesh::EntityIdVector elem2_nodes {5, 6, 7, 2};
  stk::mesh::EntityIdVector elem3_nodes {3, 8, 9, 10};
  stk::mesh::EntityIdVector elem4_nodes {8, 11, 12, 9};

  stk::mesh::EntityId elemId = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::declare_element(bulk, block_1, elemId, elem1_nodes);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    bulk.add_node_sharing(node2, 1);
    bulk.add_node_sharing(node3, 1);
  }
  else if (bulk.parallel_rank() == 1) {
    elemId = 2;
    stk::mesh::declare_element(bulk, block_2, elemId, elem2_nodes);
    elemId = 3;
    stk::mesh::declare_element(bulk, block_2, elemId, elem3_nodes);
    elemId = 4;
    stk::mesh::declare_element(bulk, block_3, elemId, elem4_nodes);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    bulk.add_node_sharing(node2, 0);
    bulk.add_node_sharing(node3, 0);
  }

  bulk.modification_end();
}

inline
void setupKeyholeMesh2D_case2(stk::mesh::BulkData& bulk)
{
  //
  //   proc 0      proc 1
  //            |
  //            | block_2 block_3
  //            |
  // block_1    |         12---11
  //            |         | 4  |
  //    4----3  | 3----6  6----10
  //    | 1  |  | |  2 |
  //    1----2  | 2----5  5----9
  //            |         | 3  |
  //            |         7----8
  //            |
  //
  //nodes 5 and 6 are ghosts (aura) on proc 0,
  //and should be members of block_2 and block_3 on proc 0
  //if edges are added, the edge between nodes 5 and 6 should
  //be a member of block_2 not block_3.
  //

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::QUAD_4_2D);
  meta.commit();

  bulk.modification_begin();

  stk::mesh::EntityIdVector elem1_nodes {1, 2, 3, 4};
  stk::mesh::EntityIdVector elem2_nodes {2, 5, 6, 3};
  stk::mesh::EntityIdVector elem3_nodes {7, 8, 9, 5};
  stk::mesh::EntityIdVector elem4_nodes {6, 10, 11, 12};

  stk::mesh::EntityId elemId = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::declare_element(bulk, block_1, elemId, elem1_nodes);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    bulk.add_node_sharing(node2, 1);
    bulk.add_node_sharing(node3, 1);
  }
  else if (bulk.parallel_rank() == 1) {
    elemId = 2;
    stk::mesh::declare_element(bulk, block_2, elemId, elem2_nodes);
    elemId = 3;
    stk::mesh::declare_element(bulk, block_3, elemId, elem3_nodes);
    elemId = 4;
    stk::mesh::declare_element(bulk, block_3, elemId, elem4_nodes);
    stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    bulk.add_node_sharing(node2, 0);
    bulk.add_node_sharing(node3, 0);
  }

  bulk.modification_end();
}


// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  5/1  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|
inline
void setupKeyholeMesh3D_case1(stk::mesh::BulkData& bulk)
{
  STK_ThrowRequire(bulk.parallel_size() == 3);
  stk::io::fill_mesh("generated:3x1x3", bulk);

  stk::mesh::EntityProcVec elementProcChanges;
  if (bulk.parallel_rank() == 1) {
    elementProcChanges.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK,4u),2));
    elementProcChanges.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK,6u),2));
  }
  bulk.change_entity_owner(elementProcChanges);
}

// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  n/a  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|
// The element in the middle has been deleted

inline
void setupKeyholeMesh3D_case2(stk::mesh::BulkData& bulk)
{
  STK_ThrowRequire(bulk.parallel_size() == 3);
  stk::io::fill_mesh("generated:3x1x3", bulk);

  stk::mesh::EntityProcVec elementProcChanges;
  if (bulk.parallel_rank() == 1) {
    elementProcChanges.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK,4),2));
    elementProcChanges.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK,6),2));
  }
  bulk.change_entity_owner(elementProcChanges);
  bulk.modification_begin();
  if (bulk.parallel_rank() == 1) {
    stk::mesh::Entity local_element5 = bulk.get_entity(stk::topology::ELEM_RANK,5);
    const bool delete_success = bulk.destroy_entity(local_element5);
    STK_ThrowRequire(delete_success);
  }
  bulk.modification_end();
}

} } } // namespace stk mesh unit_test

#endif
