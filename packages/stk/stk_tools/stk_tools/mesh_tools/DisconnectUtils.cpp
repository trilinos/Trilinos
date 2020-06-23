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

#include "stk_tools/mesh_tools/DisconnectUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>
#include <map>

namespace stk {
namespace tools {
namespace impl {

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  const bool isElementPart = (part.primary_entity_rank() == stk::topology::ELEM_RANK);
  const bool isIoPart      = stk::io::has_io_part_attribute(part);
  return (isElementPart && isIoPart);
}

stk::mesh::Part* get_block_part_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element)
{
  const stk::mesh::PartVector & elementParts = bulk.bucket(element).supersets();
  for (stk::mesh::Part * part : elementParts) {
    if (is_block(bulk, *part)) {
      return part;
    }
  }
  return nullptr;
}

int get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element)
{
  const stk::mesh::PartVector & elementParts = bulk.bucket(element).supersets();
  for (stk::mesh::Part * part : elementParts) {
    if (is_block(bulk, *part)) {
      return part->mesh_meta_data_ordinal();
    }
  }
  return -1;
}

void fill_block_membership(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::PartVector& members)
{
  const unsigned numElems = bulk.num_elements(node);
  const stk::mesh::Entity* elements = bulk.begin_elements(node);

  for(unsigned i=0; i<numElems; ++i) {
    stk::mesh::Part* block = get_block_part_for_element(bulk, elements[i]);
    ThrowRequire(block != nullptr);
    stk::util::insert_keep_sorted_and_unique(block, members, stk::mesh::PartLess());
  }
}

BlockPair get_block_pair(stk::mesh::Part* block1, stk::mesh::Part* block2)
{
  ThrowRequire(nullptr != block1 && nullptr != block2);
  ThrowRequire(block1 != block2);

  if(block2->mesh_meta_data_ordinal() > block1->mesh_meta_data_ordinal()) {
    return BlockPair(block1, block2);
  }
  return BlockPair(block2, block1);
}

void insert_block_pair(stk::mesh::Part* block1, stk::mesh::Part* block2,
                       BlockPairVector& blockPairs)
{
  BlockPair blockPair = get_block_pair(block1, block2);
  stk::util::insert_keep_sorted_and_unique(blockPair, blockPairs, PartPairLess());
}

void populate_blocks_to_reconnect(const stk::mesh::BulkData& bulk, const BlockPairVector& orderedBlockPairsInMesh,
                                  const BlockPairVector& blockPairsToDisconnect,
                                  BlockPairVector& blockPairsToReconnect)
{
  BlockPairVector orderedBlockPairsInMeshCopy;
  for(const BlockPair& blockPair : orderedBlockPairsInMesh) {
    orderedBlockPairsInMeshCopy.push_back(blockPair);
  }

  for(const BlockPair& connectPair : blockPairsToDisconnect) {
    BlockPair parts = get_block_pair(connectPair.first, connectPair.second);

    auto it = std::find(orderedBlockPairsInMeshCopy.begin(), orderedBlockPairsInMeshCopy.end(), parts);
    if(it != orderedBlockPairsInMeshCopy.end()) {
      orderedBlockPairsInMeshCopy.erase(it);
    }
  }

  for(const BlockPair blockPair : orderedBlockPairsInMeshCopy) {
    blockPairsToReconnect.push_back(blockPair);
  }
}

void fill_ordered_block_pairs(stk::mesh::PartVector& blockParts, BlockPairVector& orderedBlockPairs)
{
  if(blockParts.empty()) { return; }

  for(unsigned i = 0; i < blockParts.size() - 1; i++) {
    for(unsigned j = i+1; j < blockParts.size(); j++) {
      insert_block_pair(blockParts[i], blockParts[j], orderedBlockPairs);
    }
  }
}

stk::tools::BlockPairVector get_local_reconnect_list(const stk::mesh::BulkData& bulk, const stk::tools::BlockPairVector& disconnectList)
{
  stk::tools::BlockPairVector reconnectList;

  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<stk::mesh::Entity> sideNodes;
  stk::mesh::PartVector partMembers;

  for(const stk::tools::BlockPair& blockPair : disconnectList) {
    const stk::mesh::Part & firstBlock  = *blockPair.first;
    const stk::mesh::Part & secondBlock = *blockPair.second;
    ThrowAssert(secondBlock.mesh_meta_data_ordinal() > firstBlock.mesh_meta_data_ordinal());

    stk::mesh::Selector boundaryBetweenBlocks = firstBlock & secondBlock & (meta.locally_owned_part() | meta.globally_shared_part());

    const stk::mesh::BucketVector & nodesOnBoundaryBetweenBlocks = bulk.get_buckets(stk::topology::NODE_RANK, boundaryBetweenBlocks);
    for (const stk::mesh::Bucket * bucket : nodesOnBoundaryBetweenBlocks) {
      for (const stk::mesh::Entity node : *bucket) {
        stk::tools::impl::fill_block_membership(bulk, node, partMembers);
        stk::tools::impl::fill_ordered_block_pairs(partMembers, reconnectList);
        partMembers.clear();
      }
    }
  }

  for(const stk::tools::BlockPair& disconnectPair : disconnectList) {
    auto it = std::find_if(reconnectList.begin(), reconnectList.end(),
        [disconnectPair](const stk::tools::BlockPair& lhs)
                        {
                         return (lhs.first->mesh_meta_data_ordinal() == disconnectPair.first->mesh_meta_data_ordinal() &&
                                 lhs.second->mesh_meta_data_ordinal() == disconnectPair.second->mesh_meta_data_ordinal());
                        }
    );

    if(it != reconnectList.end()) {
      reconnectList.erase(it);
    }
  }

  return reconnectList;
}

void get_all_blocks_in_mesh(const stk::mesh::BulkData & bulk, stk::mesh::PartVector & blocksInMesh)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const stk::mesh::PartVector & allParts = meta.get_parts();
  for (stk::mesh::Part * part : allParts) {
    if (is_block(bulk, *part)) {
      blocksInMesh.push_back(part);
    }
  }
}
}}}
