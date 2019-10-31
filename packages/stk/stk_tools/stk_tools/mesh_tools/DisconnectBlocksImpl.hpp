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
#ifndef _DisconnectBlocksImpl_hpp_
#define _DisconnectBlocksImpl_hpp_

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_tools/mesh_tools/DisconnectGroup.hpp"
#include "stk_tools/mesh_tools/DisconnectUtils.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include <utility>
#include <vector>
#include <map>

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace tools {
namespace impl {

struct NodeMapKey
{
  NodeMapKey(stk::mesh::Entity _parentNode, const DisconnectGroup& _disconnectedGroup)
    : parentNode(_parentNode),
      disconnectedGroup(_disconnectedGroup) {}
  ~NodeMapKey() = default;

  stk::mesh::Entity parentNode;
  DisconnectGroup disconnectedGroup;
};

struct NodeMapValue
{
  NodeMapValue()
    : oldNodeId(stk::mesh::InvalidEntityId),
      newNodeId(stk::mesh::InvalidEntityId),
      reconnectNodeId(stk::mesh::InvalidEntityId)
      {}
  NodeMapValue(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
    : oldNodeId(bulk.identifier(node)),
      newNodeId(stk::mesh::InvalidEntityId),
      reconnectNodeId(stk::mesh::InvalidEntityId)
  {
    fill_block_membership(bulk, node, oldBlockMembership);
  }

  ~NodeMapValue() = default;

  stk::mesh::EntityId oldNodeId;
  stk::mesh::EntityId newNodeId;
  stk::mesh::EntityId reconnectNodeId;

  std::vector<int> sharingProcs;
  stk::mesh::PartVector oldBlockMembership;
};

class NodeMapLess {
public:
  NodeMapLess() = default;
  inline bool operator()(const NodeMapKey & lhs, const NodeMapKey & rhs) const
  {
    if (lhs.parentNode != rhs.parentNode) {
      return (lhs.parentNode < rhs.parentNode);
    }
    return (lhs.disconnectedGroup < rhs.disconnectedGroup);
  }
};

using BlockPairType = std::pair<stk::mesh::Part*, stk::mesh::Part*>;
using SideSetType = std::vector<stk::mesh::SideSetEntry>;
using NodeMapType = std::map<NodeMapKey, NodeMapValue, NodeMapLess>;
using DisconnectGroupVector = std::vector<DisconnectGroup>;
using PreservedSharingInfo = std::map<stk::mesh::EntityId, std::vector<int>>;

struct ReconnectNodeInfo {
  stk::mesh::EntityId reconnectNodeId = stk::mesh::InvalidEntityId;
  std::vector<int> reconnectProcs;
  stk::mesh::EntityVector relatedNodes;
};

typedef std::map<stk::mesh::EntityId, ReconnectNodeInfo> ReconnectMap;

struct LinkInfo
{
  PreservedSharingInfo sharedInfo;
  NodeMapType clonedNodeMap;
  NodeMapType preservedNodeMap;
  bool preserveOrphans = false;
  int debugLevel = 0;
  std::ostringstream os;
  ReconnectMap reconnectMap;

  void flush() {
    if(debugLevel > 1) {
      std::cerr << os.str();
    }
    os.str("");
    os.clear();
  }
};


void update_node_id(stk::mesh::EntityId newNodeId, int proc,
                    LinkInfo& info, const DisconnectGroup& group);

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part);

stk::mesh::Part* get_block_part_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

unsigned get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

void add_to_sharing_lookup(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, PreservedSharingInfo& info);

void add_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                             const BlockPairType & blockPair,
                             impl::LinkInfo& info);

void create_new_duplicate_node_IDs(stk::mesh::BulkData & bulk, LinkInfo& info);

void communicate_shared_node_information(stk::mesh::BulkData & bulk, LinkInfo& info);

void get_all_blocks_in_mesh(const stk::mesh::BulkData & bulk, stk::mesh::PartVector & blocksInMesh);

std::vector<BlockPairType> get_block_pairs_to_disconnect(const stk::mesh::BulkData & bulk);

void disconnect_elements(stk::mesh::BulkData& bulk, const DisconnectGroup& group, LinkInfo& info);

void disconnect_elements(stk::mesh::BulkData & bulk, const BlockPairType & blockPair, LinkInfo& info);

void reconnect_elements(stk::mesh::BulkData& bulk, const BlockPairType & blockPair, const DisconnectGroup& group, LinkInfo& info);

void reconnect_block_pair(stk::mesh::BulkData& bulk, const BlockPairType & blockPair, LinkInfo& info);

const std::vector<int>& find_preserved_sharing_data(stk::mesh::EntityId oldNodeId, const PreservedSharingInfo& info);

void restore_node_sharing(stk::mesh::BulkData& bulk, stk::mesh::Entity node, LinkInfo& info);

void sanitize_node_map(NodeMapType& nodeMap, std::ostringstream& os);

void determine_reconnect_node_id(stk::mesh::BulkData& bulk, const std::vector<impl::BlockPairType>& blockPairsToReconnect,
                                 LinkInfo& info);

void fix_indirect_node_sharing(stk::mesh::BulkData& bulk, const std::vector<BlockPairType>& blockPairsToReconnect, LinkInfo& info);

} } }

#endif
