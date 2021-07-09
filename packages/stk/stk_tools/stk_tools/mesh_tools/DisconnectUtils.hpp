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

#ifndef DISCONNECT_UTILS_H
#define DISCONNECT_UTILS_H

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_tools/mesh_tools/DisconnectTypes.hpp"

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace tools {
namespace impl {

struct PartPairLess {

  inline bool operator()( const BlockPair& lhs, const BlockPair& rhs ) const
  {
    if(lhs.first->mesh_meta_data_ordinal() < rhs.first->mesh_meta_data_ordinal()) {
      return true;
    }
    else if(lhs.first->mesh_meta_data_ordinal() == rhs.first->mesh_meta_data_ordinal()) {
      return lhs.second->mesh_meta_data_ordinal() < rhs.second->mesh_meta_data_ordinal();
    }
    return false;
  }
};

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part);

stk::mesh::Part* get_block_part_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

int get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

void fill_block_membership(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::PartVector& members);

stk::tools::BlockPair get_block_pair(stk::mesh::Part* block1, stk::mesh::Part* block2);

void fill_ordered_block_pairs(stk::mesh::PartVector& allBlocksInMesh, BlockPairVector& orderedBlockPairsInMesh);

void insert_block_pair(stk::mesh::Part* block1, stk::mesh::Part* block2,
                       std::vector<stk::tools::BlockPair>& blockPairs);

void populate_blocks_to_reconnect(const stk::mesh::BulkData& bulk, const BlockPairVector& orderedBlockPairsInMesh,
                                  const BlockPairVector& blockPairsToDisconnect,
                                  BlockPairVector& blockPairsToReconnect);

stk::tools::BlockPairVector get_local_reconnect_list(const stk::mesh::BulkData& bulk, const stk::tools::BlockPairVector& disconnectList);

void get_all_blocks_in_mesh(const stk::mesh::BulkData & bulk, stk::mesh::PartVector & blocksInMesh);

}}}

#endif
