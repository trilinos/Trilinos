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
#include "stk_tools/mesh_tools/CustomAura.hpp"
#include "stk_tools/mesh_tools/DetectHingesImpl.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_util/environment/WallTime.hpp"

namespace stk
{
namespace tools
{

void disconnect_all_blocks(stk::mesh::BulkData & bulk, bool preserveOrphans)
{
  if(bulk.parallel_rank() == 0) {
    std::cout << "Constructing block pairs for disconnect" << std::endl;
  }
  std::vector<BlockPair> blockPairsToDisconnect = impl::get_block_pairs_to_disconnect(bulk);

  impl::LinkInfo info;
  info.preserveOrphans = preserveOrphans;

  impl::disconnect_block_pairs(bulk, blockPairsToDisconnect, info);
}

void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect, int debugLevel)
{
  impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = debugLevel;
  stk::mesh::PartVector allBlocksInMesh;
  double startTime = stk::wall_time();

#ifndef PRINT_DEBUG
  info.debugLevel = 0;
#endif

  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);

  stk::tools::BlockPairVector orderedBlockPairsInMesh;
  stk::tools::impl::fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);
  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::impl::populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blocksToDisconnect, blockPairsToReconnect);

  double setupTime = stk::wall_time();

  stk::tools::impl::disconnect_block_pairs(bulk, orderedBlockPairsInMesh, info);

  double disconnectTime = stk::wall_time();

  stk::tools::impl::reconnect_block_pairs(bulk, blockPairsToReconnect, info);

  double reconnectTime = stk::wall_time();

  stk::tools::impl::snip_all_hinges_between_blocks(bulk, info.debugLevel > 0);

  double snipTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << "Setup time = " << (setupTime - startTime) << " s" << std::endl;
    std::cout << "Disconnect time = " << (disconnectTime - setupTime) << " s" << std::endl;
    std::cout << "Reconnect time = " << (reconnectTime- disconnectTime) << " s" << std::endl;
    std::cout << "Hinge snip time = " << (snipTime - reconnectTime) << " s" << std::endl;
    std::cout << "Overall runtime = " << (snipTime - startTime) << " s" << std::endl;
  }
}

void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockNamePairVector& blockNamesToDisconnect, int debugLevel)
{
  BlockPairVector blocksToDisconnect;
  for(const BlockNamePair& blockNamePair : blockNamesToDisconnect) {
    stk::mesh::Part* block1 = bulk.mesh_meta_data().get_part(blockNamePair.first);
    stk::mesh::Part* block2 = bulk.mesh_meta_data().get_part(blockNamePair.second);
    blocksToDisconnect.push_back(impl::get_block_pair(block1, block2));
  }

  disconnect_user_blocks(bulk, blocksToDisconnect, debugLevel);
}
}

}
