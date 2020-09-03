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
#include "stk_tools/mesh_tools/CustomAura.hpp"
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

void disconnect_user_blocks_locally(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, LinkInfo& info)
{
  stk::tools::BlockPairVector blockPairsToReconnect;
  blockPairsToReconnect = get_local_reconnect_list(bulk, blocksToDisconnect);

  disconnect_and_reconnect_blocks(bulk, blocksToDisconnect, blockPairsToReconnect, info);
}

void disconnect_user_blocks_globally(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect,
                                         LinkInfo& info)
{
  info.startTime = stk::wall_time();

  stk::mesh::PartVector allBlocksInMesh;
  BlockPairVector orderedBlockPairsInMesh;
  BlockPairVector blockPairsToReconnect;

  get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);
  populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blocksToDisconnect, blockPairsToReconnect);

  disconnect_and_reconnect_blocks(bulk, orderedBlockPairsInMesh, blockPairsToReconnect, info);
}

void snip_hinges_locally(stk::mesh::BulkData& bulk, impl::HingeNodeVector& preservedHingeNodes, LinkInfo& info)
{
  stk::mesh::EntityVector affectedNodes = extract_nodes(bulk, info);

  snip_all_hinges_for_input_nodes(bulk, affectedNodes, preservedHingeNodes);

  info.snipTime = stk::wall_time();
}

void snip_hinges_globally(stk::mesh::BulkData& bulk, LinkInfo& info)
{
  snip_all_hinges_between_blocks(bulk);

  info.snipTime = stk::wall_time();
}

void populate_hinge_node_list(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, impl::HingeNodeVector& preservedHingeNodes)
{
  stk::mesh::Selector selector;
  for(auto blockPair : blocksToDisconnect) {
    selector |= *blockPair.first & *blockPair.second;
  }

  selector &= bulk.mesh_meta_data().locally_owned_part();

  stk::mesh::EntityVector commonNodes;
  stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), commonNodes);

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
                            DisconnectOption disconnectOption, SnipOption snipOption)
{
  impl::LinkInfo info;
  info.preserveOrphans = true;

  impl::HingeNodeVector preservedHingeNodes;
  if(snipOption == SNIP_LOCAL) {
    impl::populate_hinge_node_list(bulk, blocksToDisconnect, preservedHingeNodes);
  }

  if(disconnectOption == DISCONNECT_GLOBAL) {
    impl::disconnect_user_blocks_globally(bulk, blocksToDisconnect, info);
  } else {
    impl::disconnect_user_blocks_locally(bulk, blocksToDisconnect, info);
  }

  if(snipOption == SNIP_GLOBAL) {
    impl::snip_hinges_globally(bulk, info);
  } else {
    impl::snip_hinges_locally(bulk, preservedHingeNodes, info);
  }

  impl::print_timings(bulk, info);
}

void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockNamePairVector& blockNamesToDisconnect,
                            DisconnectOption disconnectOption, SnipOption snipOption)
{
  BlockPairVector blocksToDisconnect;
  impl::populate_block_pairs(bulk, blockNamesToDisconnect, blocksToDisconnect);

  disconnect_user_blocks(bulk, blocksToDisconnect, disconnectOption, snipOption);
}

}
}