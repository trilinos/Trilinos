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

#include "stk_tools/mesh_tools/DisconnectBlocks.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshUtils.hpp"
#include "stk_tools/mesh_tools/DetectHingesImpl.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_util/environment/WallTime.hpp"

namespace stk
{
namespace tools
{

namespace impl
{
void disconnect_and_reconnect_blocks(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, BlockPairVector& blocksToReconnect, LinkInfo& info)
{
  info.setupTime = stk::wall_time();

  disconnect_block_pairs(bulk, blocksToDisconnect, info);

  info.disconnectTime = stk::wall_time();

  reconnect_block_pairs(bulk, blocksToReconnect, info);

  info.reconnectTime = stk::wall_time();
}

bool has_nodes_in_part(const stk::mesh::BulkData& bulk, const stk::mesh::Part* part)
{
  if(part == nullptr) { return false; }

  unsigned localCount = stk::mesh::count_selected_entities(*part, bulk.buckets(stk::topology::NODE_RANK));
  unsigned globalCount;

  MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED, MPI_SUM, bulk.parallel());

  return (globalCount > 0);
}

void disconnect_user_blocks_locally(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, LinkInfo& info)
{
  info.startTime = stk::wall_time();

  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::BlockPairVector sortedBlocksToDisconnect;

  for(const BlockPair& blockPair : blocksToDisconnect) {
    if(!has_nodes_in_part(bulk, blockPair.first) || !has_nodes_in_part(bulk, blockPair.second)) {
      continue;
    }
    stk::tools::impl::insert_block_pair(blockPair.first, blockPair.second, sortedBlocksToDisconnect);
  }
  blockPairsToReconnect = get_local_reconnect_list(bulk, sortedBlocksToDisconnect);

  disconnect_and_reconnect_blocks(bulk, sortedBlocksToDisconnect, blockPairsToReconnect, info);
}

void snip_hinges(stk::mesh::BulkData& bulk, impl::HingeNodeVector& preservedHingeNodes, const BlockPairVector& blocksToDisconnect, LinkInfo& info)
{
  stk::mesh::EntityVector affectedNodes = get_affected_nodes(bulk, blocksToDisconnect);

  snip_all_hinges_for_input_nodes(bulk, affectedNodes, preservedHingeNodes);

  info.snipTime = stk::wall_time();
}

void populate_hinge_node_list(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, impl::HingeNodeVector& preservedHingeNodes)
{
  stk::mesh::EntityVector commonNodes = get_affected_nodes(bulk, blocksToDisconnect);
  HingeNodeVector commonHingeNodes = impl::get_hinge_nodes(bulk, commonNodes);

  for(auto hingeNode : commonHingeNodes) {
    preservedHingeNodes.push_back(hingeNode);
  }
}

void populate_block_pairs(stk::mesh::BulkData& bulk, const BlockNamePairVector& blockNamesToDisconnect, BlockPairVector& blocksToDisconnect)
{
  for(const BlockNamePair& blockNamePair : blockNamesToDisconnect) {
    stk::mesh::Part* block1 = bulk.mesh_meta_data().get_part(blockNamePair.first);
    stk::mesh::Part* block2 = bulk.mesh_meta_data().get_part(blockNamePair.second);
    blocksToDisconnect.push_back(get_block_pair(block1, block2));
  }
}

void print_timings(stk::mesh::BulkData& bulk, LinkInfo& info)
{
  if (bulk.parallel_rank() == 0) {
    std::cout << "Setup time = " << (info.setupTime - info.startTime) << " s" << std::endl;
    std::cout << "Disconnect time = " << (info.disconnectTime - info.setupTime) << " s" << std::endl;
    std::cout << "Reconnect time = " << (info.reconnectTime - info.disconnectTime) << " s" << std::endl;
    std::cout << "Hinge snip time = " << (info.snipTime - info.reconnectTime) << " s" << std::endl;
    std::cout << "Overall runtime = " << (info.snipTime - info.startTime) << " s" << std::endl;
  }
}
}

void disconnect_all_blocks(stk::mesh::BulkData & bulk, bool preserveOrphans)
{
  impl::LinkInfo info;
  info.preserveOrphans = preserveOrphans;
  disconnect_all_blocks(bulk, info, preserveOrphans);
}

void disconnect_all_blocks(stk::mesh::BulkData & bulk, impl::LinkInfo& info, bool preserveOrphans)
{
  std::vector<BlockPair> blockPairsToDisconnect = impl::get_block_pairs_to_disconnect(bulk);

  impl::disconnect_block_pairs(bulk, blockPairsToDisconnect, info);
}

void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect,
                            SnipOption snipOption)
{
  impl::LinkInfo info;
  info.preserveOrphans = (snipOption == PRESERVE_INITIAL_HINGES) ? true : false;

  impl::HingeNodeVector preservedHingeNodes;
  if(info.preserveOrphans) {
    impl::populate_hinge_node_list(bulk, blocksToDisconnect, preservedHingeNodes);
  }

  impl::disconnect_user_blocks_locally(bulk, blocksToDisconnect, info);

  impl::snip_hinges(bulk, preservedHingeNodes, blocksToDisconnect, info);

  // impl::print_timings(bulk, info);
}

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Sep 2024
void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockPairVector& blockPairsToDisconnect,
                            DisconnectBlocksOption options)
{
  disconnect_user_blocks(bulk, blockPairsToDisconnect, options.snipOption);
}
#endif

}
}
