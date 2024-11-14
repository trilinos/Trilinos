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
#include "stk_tools/mesh_tools/DisconnectTypes.hpp"
#include "stk_tools/mesh_tools/DisconnectUtils.hpp"
#include "stk_tools/mesh_tools/ConvexGroup.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include <map>
#include <utility>
#include <vector>

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
      reconnectNodeId(stk::mesh::InvalidEntityId),
      reconnectGroupId(std::numeric_limits<unsigned>::max())
  {}
  NodeMapValue(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
    : boundaryNode(node),
      oldNodeId(bulk.identifier(node)),
      newNodeId(stk::mesh::InvalidEntityId),
      reconnectNodeId(stk::mesh::InvalidEntityId),
      reconnectGroupId(std::numeric_limits<unsigned>::max())
  {
    fill_block_membership(bulk, node, oldBlockMembership);
  }

  ~NodeMapValue() = default;

  stk::mesh::Entity boundaryNode;
  stk::mesh::EntityId oldNodeId;
  stk::mesh::EntityId newNodeId;
  stk::mesh::EntityId reconnectNodeId;

  std::vector<int> sharingProcs;
  stk::mesh::PartVector oldBlockMembership;

  unsigned reconnectGroupId;
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


using SideSetType = std::vector<stk::mesh::SideSetEntry>;
using NodeMapType = std::map<NodeMapKey, NodeMapValue, NodeMapLess>;
using DisconnectGroupVector = std::vector<DisconnectGroup>;
using PreservedSharingInfo = std::map<stk::mesh::EntityId, std::vector<int>>;
using NodeMapIterator = NodeMapType::iterator;
using ReconnectMapKey = std::pair<stk::mesh::EntityId, unsigned>;
using InternalSidePerPair = std::vector<std::pair<stk::mesh::Entity, BlockPair>>;

struct ReconnectNodeInfo {
  stk::mesh::EntityId reconnectNodeId = stk::mesh::InvalidEntityId;
  std::vector<int> reconnectProcs;
  stk::mesh::EntityVector relatedNodes;
  stk::mesh::PartVector reconnectParts;
  //  ConvexGroup<BlockPair,BlockPairIdGetter,stk::mesh::PartLess> group;
};

typedef std::map<ReconnectMapKey, ReconnectNodeInfo> ReconnectMap;

class NullStream : public std::ostream {
  class NullBuffer : public std::streambuf {
  public:
    int overflow( int c ) { return c; }
  } m_nb;
public:
  NullStream() : std::ostream( &m_nb ) {}
};

enum class InternalFaceDisconnectState {KEEP, DESTROY, NOOP};

inline std::ostream &operator<<(std::ostream &out, const InternalFaceDisconnectState &t)
{
  switch (t) {
  case InternalFaceDisconnectState::KEEP: return out << "KEEP"; break;
  case InternalFaceDisconnectState::DESTROY: return out << "DESTROY"; break;
  case InternalFaceDisconnectState::NOOP: return out << "NOOP"; break;
  default: return out << "INVALID"; break;
  }
  return out << "INVALID[" << (unsigned)t << "]";
}

struct InternalFaceInfo
{
  BlockPair blockPair;
  stk::mesh::Entity elemInFirstBlock;
  stk::mesh::ConnectivityOrdinal faceOrdinalInFirstBlock;
  stk::mesh::Entity elemInSecondBlock;
  stk::mesh::ConnectivityOrdinal faceOrdinalInSecondBlock;
  stk::mesh::Entity internalFace;
  std::pair<bool, bool> isInBlockPairSideset;
  std::map<std::pair<stk::mesh::Entity,stk::mesh::ConnectivityOrdinal>, stk::mesh::Entity> nodeMap;

  InternalFaceInfo(const BlockPair& blockPair_,
                   const stk::mesh::Entity elemInFirstBlock_,
                   const stk::mesh::ConnectivityOrdinal faceOrdinalInFirstBlock_,
                   const stk::mesh::Entity elemInSecondBlock_,
                   const stk::mesh::ConnectivityOrdinal faceOrdinalInSecondBlock_,
                   const stk::mesh::Entity internalFace_,
                   const std::pair<bool, bool>& isInBlockPairSideset_)
  : blockPair(blockPair_)
  , elemInFirstBlock(elemInFirstBlock_)
  , faceOrdinalInFirstBlock(faceOrdinalInFirstBlock_)
  , elemInSecondBlock(elemInSecondBlock_)
  , faceOrdinalInSecondBlock(faceOrdinalInSecondBlock_)
  , internalFace(internalFace_)
  , isInBlockPairSideset(isInBlockPairSideset_) {}

  operator stk::mesh::Entity() const { return internalFace; }

  bool operator<(const InternalFaceInfo &rhs) const
  {
    return internalFace < rhs.internalFace;
  };

  bool operator<(const stk::mesh::Entity &rhs) const
  {
    return internalFace < rhs;
  };

  bool operator==(const InternalFaceInfo &rhs) const
  {
    return internalFace == rhs.internalFace;
  };

  bool operator==(const stk::mesh::Entity &rhs) const
  {
    return internalFace == rhs;
  };

  bool operator!=(const InternalFaceInfo &rhs) const
  {
    return internalFace != rhs.internalFace;
  };

  bool operator!=(const stk::mesh::Entity &rhs) const
  {
    return internalFace != rhs;
  };
};

struct LinkInfo
{
  using DisconnectInternalFaces = std::vector<InternalFaceInfo>;

  PreservedSharingInfo sharedInfo;
  NodeMapType clonedNodeMap;
  NodeMapType originalNodeMap;
  bool preserveOrphans = false;
  int debugLevel = 0;
  std::string debugString = "";
  std::ostringstream os;
  NullStream ns;
  ReconnectMap reconnectMap;
  DisconnectInternalFaces internalSides;

  double startTime;
  double setupTime;
  double disconnectTime;
  double reconnectTime;
  double snipTime;

  void flush(std::ostream& stream) {
    stream << os.str();
    os.str("");
    os.clear();
  }

  void flush() {
    flush(std::cerr);
  }

  std::ostream& print_debug_msg(int userDebugLevel, int mpiRank, bool prefixMsg = true) {
    if(userDebugLevel <= debugLevel) {
      if(prefixMsg) {
        os << "P" << mpiRank << ": ";
      }
      return os;
    } else {
      return ns;
    }
  }

  std::ostream& print_debug_msg_p0(int userDebugLevel, int mpiRank, bool prefixMsg = true) {
    if(mpiRank == 0) {
      return print_debug_msg(userDebugLevel, prefixMsg);
    }
    return ns;
  }
};

void clean_up_aura(stk::mesh::BulkData& bulk, LinkInfo& info);

void update_node_id(stk::mesh::EntityId newNodeId, int proc,
                    LinkInfo& info, const DisconnectGroup& group);

bool should_be_reconnected(const DisconnectGroup& disconnectedGroup, const NodeMapValue& nodeMapValue, const stk::mesh::Part& blockToReconnect,
                           const stk::mesh::Part& srcBlock, LinkInfo& info);

bool can_be_reconnected(const DisconnectGroup& disconnectedGroup, const NodeMapValue& nodeMapValue, const BlockPair& blockPair, LinkInfo& info);

bool can_be_reconnected(const DisconnectGroup& disconnectedGroup, const NodeMapValue& nodeMapValue, const BlockPair& blockPair, stk::mesh::Entity currentEntity, LinkInfo& info);

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part);

stk::mesh::Part* get_block_part_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

int get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

void add_to_sharing_lookup(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, PreservedSharingInfo& info);

void add_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                             const BlockPair & blockPair,
                             LinkInfo& info);

void add_sidesets_to_disconnect(stk::mesh::BulkData const& bulk,
                                BlockPair const& blockPair,
                                LinkInfo& info);

void create_new_duplicate_node_IDs(stk::mesh::BulkData & bulk, LinkInfo& info);

void communicate_shared_node_information(stk::mesh::BulkData & bulk, LinkInfo& info);

void get_all_blocks_in_mesh(const stk::mesh::BulkData & bulk, stk::mesh::PartVector & blocksInMesh);

std::vector<BlockPair> get_block_pairs_to_disconnect(const stk::mesh::BulkData & bulk);

void disconnect_elements(stk::mesh::BulkData& bulk, const NodeMapKey& key, NodeMapValue& value, LinkInfo& info);

void disconnect_elements(stk::mesh::BulkData & bulk, const BlockPair & blockPair, LinkInfo& info);

void reconnect_elements(stk::mesh::BulkData& bulk, const BlockPair & blockPair, const NodeMapKey& key, const NodeMapValue& value, LinkInfo& info);

void reconnect_block_pair(stk::mesh::BulkData& bulk, const BlockPair & blockPair, LinkInfo& info);

const std::vector<int>& find_preserved_sharing_data(stk::mesh::EntityId oldNodeId, const PreservedSharingInfo& info);

void restore_node_sharing(stk::mesh::BulkData& bulk, const NodeMapValue& value, stk::mesh::Entity node, LinkInfo& info);

void sanitize_node_map(NodeMapType& nodeMap, LinkInfo& os);

stk::mesh::EntityVector extract_nodes(const stk::mesh::BulkData& bulk, LinkInfo& info);

stk::mesh::EntityVector get_affected_nodes(const stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect);

void disconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToDisconnect,
                            LinkInfo& info);
void reconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToDisconnect,
                           LinkInfo& info);
} } }

#endif
