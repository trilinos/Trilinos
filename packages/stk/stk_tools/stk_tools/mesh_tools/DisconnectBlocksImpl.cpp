#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocks.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>
#include <map>
#include <stk_util/environment/RuntimeWarning.hpp>

namespace stk {
namespace tools {
namespace impl {


void add_to_sharing_lookup(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, PreservedSharingInfo& info)
{
  stk::mesh::EntityId id = bulk.identifier(node);

  auto iter = info.find(id);
  if(iter == info.end()) {
    std::vector<int> commSharedProcs;
    bulk.comm_shared_procs(bulk.entity_key(node), commSharedProcs);
    info.insert(std::make_pair(id, commSharedProcs));
  }
}

void create_new_node_map_entry(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node, const BlockPair& blockPair, LinkInfo& info)
{
  const stk::mesh::Part & secondBlock = *blockPair.second;
  const stk::mesh::Entity * elems = bulk.begin_elements(node);
  const unsigned numElems = bulk.num_elements(node);
  bool needToCloneNode = false;
  for (unsigned i = 0; i < numElems; ++i) {
    const stk::mesh::Bucket & elemBucket = bulk.bucket(elems[i]);
    if (elemBucket.member(secondBlock) && elemBucket.owned()) {
      needToCloneNode = true;
      break;
    }
  }
#ifdef PRINT_DEBUG
  const stk::mesh::Part & firstBlock  = *blockPair.first;
  info.print_debug_msg(1) << "blockpairs " << firstBlock.name() << " " << secondBlock.name()
              << " needToClone: " << needToCloneNode << " node id: " << bulk.identifier(node) << std::endl;
#endif
  add_to_sharing_lookup(bulk, node, info.sharedInfo);
  if (needToCloneNode && (bulk.state(node) != stk::mesh::Created)) {
    DisconnectGroup group(bulk, &secondBlock, node);
    info.clonedNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
#ifdef PRINT_DEBUG
    info.print_debug_msg(1) << "cloning node " << bulk.identifier(node) << " group id " <<
        group.id() << " part name: " << secondBlock.name() << std::endl;
#endif
  } else {
    DisconnectGroup group(bulk, &secondBlock, node);
    info.preservedNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
#ifdef PRINT_DEBUG
    info.print_debug_msg(1) << "skipping node " << bulk.identifier(node) << " group id " <<
        group.id() << " part name: " << secondBlock.name() << std::endl;
#endif
  }
}

bool has_reconnect_node_id(LinkInfo& info, stk::mesh::EntityId referenceId)
{
  auto reconnectMapIter = info.reconnectMap.find(referenceId);
  return (reconnectMapIter != info.reconnectMap.end());
}

void set_reconnect_node_id(LinkInfo& info, stk::mesh::EntityId referenceId, stk::mesh::EntityId reconnectId)
{
   auto reconnectMapIter = info.reconnectMap.find(referenceId);
   if(reconnectMapIter != info.reconnectMap.end()) {
     reconnectMapIter->second.reconnectNodeId = reconnectId;
   }
}

std::vector<int> get_reconnect_procs(LinkInfo& info, stk::mesh::EntityId referenceId)
{
  auto reconnectMapIter = info.reconnectMap.find(referenceId);
  ThrowRequire(reconnectMapIter != info.reconnectMap.end());
  return reconnectMapIter->second.reconnectProcs;
}

stk::mesh::EntityId get_reconnect_node_id(LinkInfo& info, stk::mesh::EntityId referenceId)
{
  auto reconnectMapIter = info.reconnectMap.find(referenceId);
  ThrowRequire(reconnectMapIter != info.reconnectMap.end());
  return reconnectMapIter->second.reconnectNodeId;
}

void add_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                             const BlockPair & blockPair,
                             LinkInfo& info)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<stk::mesh::Entity> sideNodes;

  const stk::mesh::Part & firstBlock  = *blockPair.first;
  const stk::mesh::Part & secondBlock = *blockPair.second;
  ThrowAssert(secondBlock.mesh_meta_data_ordinal() > firstBlock.mesh_meta_data_ordinal());

  stk::mesh::Selector boundaryBetweenBlocks = firstBlock & secondBlock & (meta.locally_owned_part() | meta.globally_shared_part());

  const stk::mesh::BucketVector & nodesOnBoundaryBetweenBlocks = bulk.get_buckets(stk::topology::NODE_RANK, boundaryBetweenBlocks);
  for (const stk::mesh::Bucket * bucket : nodesOnBoundaryBetweenBlocks) {
    for (const stk::mesh::Entity node : *bucket) {
      create_new_node_map_entry(bulk, node, blockPair, info);
    }
  }
}

void create_new_duplicate_node_IDs(stk::mesh::BulkData & bulk, LinkInfo& info)
{
  std::vector<stk::mesh::EntityId> newNodeIDs;
  bulk.generate_new_ids(stk::topology::NODE_RANK, info.clonedNodeMap.size(), newNodeIDs);

  size_t newNodeIdx = 0;
  for (auto & nodeMapEntry : info.clonedNodeMap) {
    nodeMapEntry.second.newNodeId = newNodeIDs[newNodeIdx++];
#ifdef PRINT_DEBUG
    info.print_debug_msg(2) << "locally assigning new id " << newNodeIDs[newNodeIdx-1] << " to node " <<
        bulk.identifier(nodeMapEntry.first.parentNode) << std::endl;
#endif
  }
}

void pack_shared_node_information(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, LinkInfo& info)
{
  for (const auto & nodeMapEntry : info.clonedNodeMap) {
    const stk::mesh::Entity node = nodeMapEntry.first.parentNode;
    if (!bulk.bucket(node).shared()) continue;
    const DisconnectGroup & group = nodeMapEntry.first.disconnectedGroup;
    ThrowRequire(node == group.get_node());

    const stk::mesh::EntityId newNodeId = nodeMapEntry.second.newNodeId;
    bool needToCommunicate = group.needs_to_communicate();

    if (needToCommunicate) {
      std::vector<int> sharingProcs;
      bulk.comm_shared_procs(bulk.entity_key(node), sharingProcs);
      for (const int proc : sharingProcs) {
        stk::CommBuffer& procBuff = commSparse.send_buffer(proc);
        group.pack_group_info(procBuff, newNodeId, proc, info.print_debug_msg(3,false));
      }
    }
  }
}

void unpack_shared_node_information(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, LinkInfo& info)
{
  for (int proc = 0; proc < bulk.parallel_size(); ++proc) {
    while (commSparse.recv_buffer(proc).remaining()) {
      stk::CommBuffer & procBuff = commSparse.recv_buffer(proc);
      stk::mesh::EntityId newNodeId;
      DisconnectGroup group(bulk);
      group.unpack_group_info(procBuff, newNodeId, proc, info.print_debug_msg(3,false));
      update_node_id(newNodeId, proc, info, group);
    }
  }
}

void communicate_shared_node_information(stk::mesh::BulkData & bulk, LinkInfo& info)
{
  stk::CommSparse commSparse(bulk.parallel());

  for(int phase = 0; phase < 2; ++phase) {
    pack_shared_node_information(bulk, commSparse, info);

    if (phase == 0) {
      commSparse.allocate_buffers();
    } else {
      commSparse.communicate();
    }
  }
  unpack_shared_node_information(bulk, commSparse, info);
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

std::vector<BlockPair> get_block_pairs_to_disconnect(const stk::mesh::BulkData & bulk)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<BlockPair> blockPairsToDisconnect;
  stk::mesh::PartVector allBlocksInMesh;
  get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  if (allBlocksInMesh.size() < 2u) return blockPairsToDisconnect;

  for (size_t firstBlockIdx = 0; firstBlockIdx < allBlocksInMesh.size()-1; ++firstBlockIdx) {
    for (size_t secondBlockIdx = firstBlockIdx+1; secondBlockIdx < allBlocksInMesh.size(); ++secondBlockIdx) {
      stk::mesh::Part & firstBlock  = *allBlocksInMesh[firstBlockIdx];
      stk::mesh::Part & secondBlock = *allBlocksInMesh[secondBlockIdx];
      stk::mesh::Selector boundaryBetweenBlocks = firstBlock & secondBlock & (meta.locally_owned_part() | meta.globally_shared_part());
      const stk::mesh::BucketVector & nodesOnBoundaryBetweenBlocks = bulk.get_buckets(stk::topology::NODE_RANK, boundaryBetweenBlocks);
      if (!nodesOnBoundaryBetweenBlocks.empty()) {
        blockPairsToDisconnect.emplace_back(&firstBlock, &secondBlock);
      }
    }
  }
  return blockPairsToDisconnect;
}

void update_disconnected_element_relation(stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::Entity newNode,
                                          stk::mesh::Entity elem, LinkInfo& info)
{
  unsigned numNodes = bulk.num_connectivity(elem, stk::topology::NODE_RANK);
  const stk::mesh::Entity * elemNodes = bulk.begin(elem, stk::topology::NODE_RANK);
  stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(elem, stk::topology::NODE_RANK);
  for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
    if (elemNodes[iNode] == node) {
      bulk.destroy_relation(elem, node, nodeOrdinals[iNode]);
      bulk.declare_relation(elem, newNode, nodeOrdinals[iNode]);
    }
  }

  if (bulk.num_elements(node) == 0 && !info.preserveOrphans) {
    bulk.destroy_entity(node);
  }
}

void disconnect_elements(stk::mesh::BulkData& bulk, const NodeMapKey& key, NodeMapValue& value, LinkInfo& info)
{
  const DisconnectGroup& group = key.disconnectedGroup;
  if (group.is_active()) {
    const stk::mesh::Entity node = key.parentNode;

    for (stk::mesh::Entity elem : group.get_entities()) {
      stk::mesh::EntityId newNodeId = value.newNodeId;
      stk::mesh::Entity newNode = bulk.declare_node(newNodeId);
      value.boundaryNode = newNode;
      bulk.copy_entity_fields(node, newNode);

      for (int sharingProc : value.sharingProcs) {
        bulk.add_node_sharing(newNode, sharingProc);
      }
      update_disconnected_element_relation(bulk, node, newNode, elem, info);
    }
    group.set_active(false);
  }
}

void disconnect_elements(stk::mesh::BulkData & bulk, const BlockPair & blockPair, LinkInfo& info)
{
  const stk::mesh::Part & blockToDisconnect = *blockPair.second;

  for (auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    const DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    const stk::mesh::ConstPartVector& parts = disconnectedGroup.get_parts();

    if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == blockToDisconnect.mesh_meta_data_ordinal()) {
      disconnect_elements(bulk, nodeMapEntryIt->first, nodeMapEntryIt->second, info);
    }
  }
}

std::vector<int> get_node_sharing_for_restoration(stk::mesh::BulkData& bulk, const stk::mesh::Part* blockPart,
                                                  const DisconnectGroup& group, stk::mesh::Entity destNode, LinkInfo& info)
{
  std::vector<int> sharingProcs;

  ThrowRequire(bulk.is_valid(destNode));
  std::vector<int> procs = group.get_sharing_procs(*blockPart);

#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << "printing sharing procs for part: " << blockPart->name() << std::endl;
#endif
  for(int proc : procs) {
    if(proc == bulk.parallel_rank()) { continue; }
    sharingProcs.push_back(proc);
#ifdef PRINT_DEBUG
    info.print_debug_msg(3,false) << "\tP" << proc;
#endif
  }
#ifdef PRINT_DEBUG
  info.print_debug_msg(3,false) << std::endl;
#endif

  return sharingProcs;
}

void restore_node_sharing(stk::mesh::BulkData& bulk, const stk::mesh::EntityId referenceId,
                          stk::mesh::Entity destNode, LinkInfo& info)
{
  std::vector<int> commonProcs = get_reconnect_procs(info, referenceId);

  for(int proc : commonProcs) {
    bulk.add_node_sharing(destNode, proc);
#ifdef PRINT_DEBUG
    stk::mesh::EntityId nodeId = bulk.identifier(destNode);
    info.print_debug_msg(3) << "sharing node "  << nodeId << " to proc " << proc << std::endl;
#endif
  }
}

void update_reconnected_element_relation(stk::mesh::BulkData& bulk, const DisconnectGroup& group, stk::mesh::Entity newNode, stk::mesh::Entity reconnectNode)
{
  for (stk::mesh::Entity elem : group.get_entities()) {
    bulk.copy_entity_fields(newNode, reconnectNode);

    unsigned numNodes = bulk.num_connectivity(elem, stk::topology::NODE_RANK);
    const stk::mesh::Entity * elemNodes = bulk.begin(elem, stk::topology::NODE_RANK);
    stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(elem, stk::topology::NODE_RANK);
    for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
      if (elemNodes[iNode] == newNode) {
        bulk.destroy_relation(elem, newNode, nodeOrdinals[iNode]);
        bulk.declare_relation(elem, reconnectNode, nodeOrdinals[iNode]);
      }
    }
  }
}

void reconnect_elements(stk::mesh::BulkData& bulk, const BlockPair & blockPair, const NodeMapKey& key, const NodeMapValue& value,
                        LinkInfo& info)
{
  auto& disconnectedGroup = key.disconnectedGroup;

  stk::mesh::EntityId reconnectId = get_reconnect_node_id(info, value.oldNodeId);
  stk::mesh::Entity reconnectNode = bulk.get_entity(stk::topology::NODE_RANK, reconnectId);

  if(!bulk.is_valid(reconnectNode)) {
#ifdef PRINT_DEBUG
    info.print_debug_msg(2) << "declaring reconnect node " << reconnectId << std::endl;
#endif
    reconnectNode = bulk.declare_node(reconnectId);
  }

  stk::mesh::EntityId newNodeId = value.newNodeId;
  stk::mesh::Entity newNode = bulk.get_entity(stk::topology::NODE_RANK, newNodeId);
  ThrowRequire(bulk.is_valid(newNode));

  restore_node_sharing(bulk, value.oldNodeId, reconnectNode, info);

#ifdef PRINT_DEBUG
  info.print_debug_msg(2) << "Restoring node " << newNodeId << " to " << reconnectId  << std::endl;
#endif
  update_reconnected_element_relation(bulk, disconnectedGroup, newNode, reconnectNode);
}

struct BlockPairComplement {
  stk::mesh::Part* commonPart;
  BlockPair blockPair;

  BlockPairComplement() :
    commonPart(nullptr),
    blockPair(BlockPair(nullptr, nullptr)) {}
};

BlockPairComplement get_block_pair_complement(const BlockPair blockPair1, const BlockPair blockPair2)
{
  BlockPairComplement complement;

  if(blockPair1.first == blockPair2.first) {
    complement.commonPart = blockPair1.first;
    complement.blockPair.first = blockPair1.second;
    complement.blockPair.second = blockPair2.second;
  }
  else if(blockPair1.first == blockPair2.second) {
    complement.commonPart = blockPair1.first;
    complement.blockPair.first = blockPair1.second;
    complement.blockPair.second = blockPair2.first;
  }
  else if(blockPair1.second == blockPair2.first) {
    complement.commonPart = blockPair1.second;
    complement.blockPair.first = blockPair1.first;
    complement.blockPair.second = blockPair2.second;
  }
  else if(blockPair1.second == blockPair2.second) {
    complement.commonPart = blockPair1.second;
    complement.blockPair.first = blockPair1.first;
    complement.blockPair.second = blockPair2.first;
  }

  return complement;
}

const std::vector<int>& find_preserved_sharing_data(stk::mesh::EntityId oldNodeId, const PreservedSharingInfo& info)
{
  auto iter = info.find(oldNodeId);
  ThrowRequire(iter != info.end());
  return iter->second;
}

bool is_in_node_map(NodeMapType& nodeMap, const DisconnectGroup& group)
{
  for(auto it = nodeMap.begin(); it != nodeMap.end(); ++it){
    if(it->first.disconnectedGroup.get_node() == group.get_node()) {
      return true;
    }
  }
  return false;
}

bool should_be_reconnected(NodeMapType::iterator nodeMapEntryIt, const stk::mesh::Part& blockToReconnect,
                           const stk::mesh::Part& srcBlock, LinkInfo& info)
{
  const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
  const stk::mesh::PartVector& blockMembership = nodeMapEntryIt->second.oldBlockMembership;
  const stk::mesh::ConstPartVector& parts = disconnectedGroup.get_parts();
  bool isOriginalMember = false;

#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << "parts size: " << parts.size() << " : name "
            << parts[0]->name() << " group id: " << disconnectedGroup.id() << std::endl;
#endif

  if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == blockToReconnect.mesh_meta_data_ordinal()) {
    isOriginalMember = std::binary_search(blockMembership.begin(), blockMembership.end(), &srcBlock, stk::mesh::PartLess());
#ifdef PRINT_DEBUG
    const stk::mesh::BulkData& bulk = disconnectedGroup.get_bulk();
    info.print_debug_msg(3) << "Checking to see if part " << srcBlock.name() << " is original member of old node " << nodeMapEntryIt->second.oldNodeId << " : isOriginalMember = "<< isOriginalMember << std::endl;

    for(stk::mesh::Part* part : blockMembership) {
      info.print_debug_msg(4, false) << "\t\tP" << bulk.parallel_rank() << ": " << part->name() << std::endl;
    }
#endif
  }
  return isOriginalMember;
}

bool should_update_sharing_info(NodeMapType& nodeMap, const DisconnectGroup& group,
                                const stk::mesh::Part& blockToReconnect, LinkInfo& info)
{
#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << "Checking to see if to update sharing info for old node " << group.get_node_id() << std::endl;
#endif
  for(auto it = nodeMap.begin(); it != nodeMap.end(); ++it){
    if(it->first.disconnectedGroup.get_node() == group.get_node()) {
      const stk::mesh::ConstPartVector& parts = it->first.disconnectedGroup.get_parts();
      if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() != blockToReconnect.mesh_meta_data_ordinal()) {
        const stk::mesh::Part* block1 = parts[0];
        const stk::mesh::Part* block2 = &blockToReconnect;

        if(blockToReconnect.mesh_meta_data_ordinal() < parts[0]->mesh_meta_data_ordinal()) {
          block1 = &blockToReconnect;
          block2 = parts[0];
        }

        if(should_be_reconnected(it, *block2, *block1, info)) {
#ifdef PRINT_DEBUG
          info.print_debug_msg(4) << "returning true" << std::endl;
#endif
          return true;
        }
      }
    }
  }
#ifdef PRINT_DEBUG
  info.print_debug_msg(4) << "returning false" << std::endl;
#endif
  return false;
}

void sanitize_node_map(NodeMapType& nodeMap, LinkInfo& info)
{
  for (NodeMapType::iterator nodeMapEntryIt = nodeMap.begin(); nodeMapEntryIt != nodeMap.end();) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if(!disconnectedGroup.get_bulk().is_valid(disconnectedGroup.get_node())) {
#ifdef PRINT_DEBUG
      info.print_debug_msg(2)  << "sanitizing node id " << nodeMapEntryIt->second.oldNodeId << std::endl;
#endif
      nodeMapEntryIt = nodeMap.erase(nodeMapEntryIt);
    } else {
      nodeMapEntryIt++;
    }
  }
}

void insert_uniquely_reconnect_info(const int myRank, const std::vector<int>& procs, ReconnectNodeInfo& reconnectInfo)
{
  for(unsigned i = 0; i < procs.size(); i++) {
    if(myRank == procs[i]) { continue; }
    stk::util::insert_keep_sorted_and_unique(procs[i], reconnectInfo.reconnectProcs);
  }
}

void pack_reconnect_node_information(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, LinkInfo& info)
{
  for (const auto & reconnectMapEntry : info.reconnectMap) {
    for (const int proc : reconnectMapEntry.second.reconnectProcs) {
      stk::CommBuffer& procBuff = commSparse.send_buffer(proc);
      procBuff.pack<stk::mesh::EntityId>(reconnectMapEntry.first);
      procBuff.pack<stk::mesh::EntityId>(reconnectMapEntry.second.reconnectNodeId);
#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << "passing info about ref node: " << reconnectMapEntry.first << " to proc " << proc << std::endl;
#endif
    }
  }
}

void unpack_reconnect_node_information(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, LinkInfo& info)
{
  for (int proc = 0; proc < bulk.parallel_size(); ++proc) {
    while (commSparse.recv_buffer(proc).remaining()) {
      stk::CommBuffer & procBuff = commSparse.recv_buffer(proc);
      stk::mesh::EntityId referenceNodeId;
      stk::mesh::EntityId reconnectNodeId;

      procBuff.unpack<stk::mesh::EntityId>(referenceNodeId);
      procBuff.unpack<stk::mesh::EntityId>(reconnectNodeId);

#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << "recved info about ref node: " << referenceNodeId << " from proc " << proc << std::endl;
#endif
      auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
      ThrowRequire(reconnectMapIter != info.reconnectMap.end());

#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << "comparing stored reconnectNodeId: " << reconnectMapIter->second.reconnectNodeId
          << " and recv reconnectNodeId: " << reconnectNodeId << std::endl;
#endif
      reconnectMapIter->second.reconnectNodeId = std::min(reconnectMapIter->second.reconnectNodeId, reconnectNodeId);
    }
  }
}

void communicate_reconnect_node_information(stk::mesh::BulkData & bulk, LinkInfo& info)
{
  stk::CommSparse commSparse(bulk.parallel());

  for(int phase = 0; phase < 2; ++phase) {
    pack_reconnect_node_information(bulk, commSparse, info);

    if (phase == 0) {
      commSparse.allocate_buffers();
    }
    else {
      commSparse.communicate();
    }
  }
  unpack_reconnect_node_information(bulk, commSparse, info);

#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << "reconnect node info " << std::endl;

  for(auto iter : info.reconnectMap) {
    info.print_debug_msg(3,false) << "\tRef: " << iter.first << " reconnectId: "<< iter.second.reconnectNodeId << std::endl;
    for(stk::mesh::Entity node : iter.second.relatedNodes) {
      info.print_debug_msg(3,false) << "\t\tNodes: " << bulk.identifier(node);
    }
    info.print_debug_msg(3,false) << std::endl;
  }
  info.print_debug_msg(3,false) << std::endl;
#endif
}

void update_reconnect_procs_for_block_pair(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& transitiveBlockList,
                                           NodeMapType::iterator it, stk::mesh::Entity currentEntity, LinkInfo& info)
{
  const DisconnectGroup& disconnectedGroup = it->first.disconnectedGroup;
  stk::mesh::Entity boundaryNode = it->second.boundaryNode;

#ifdef PRINT_DEBUG
  info.print_debug_msg(1) << "checking ref id for " << it->second.oldNodeId << std::endl;
  info.print_debug_msg(1) << "update transitive reconnect id. first: " << transitiveBlockList[0]->name()
                      << " second: " << transitiveBlockList[1]->name() << " third: " << transitiveBlockList[2]->name()
                      << " boundary node id: " << bulk.identifier(boundaryNode)
                      << " currentEntity id: " << bulk.identifier(currentEntity) << std::endl;
#endif
  auto reconnectMapIter = info.reconnectMap.find(it->second.oldNodeId);
  if(reconnectMapIter == info.reconnectMap.end()) { return; }
  ReconnectNodeInfo* reconnectInfo = &reconnectMapIter->second;
  std::vector<int> sharingProcs;

  for(const stk::mesh::Part* blockPart : transitiveBlockList) {
    sharingProcs = get_node_sharing_for_restoration(bulk, blockPart, disconnectedGroup, boundaryNode, info);

    for(int proc : sharingProcs) {
      stk::util::insert_keep_sorted_and_unique(proc, reconnectInfo->reconnectProcs);
    }
  }
}

template <typename Functor>
void traverse_transitive_relations(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                   NodeMapType::iterator nodeMapEntryIt, stk::mesh::Entity currentEntity,
                                   LinkInfo& info, Functor func)
{
  if(blockPairsToReconnect.empty()) { return; }
  if(!bulk.is_valid(currentEntity)) { return; }

  BlockPairVector blockPairs;
  stk::mesh::PartLess partLess;
  const stk::mesh::PartVector& oldBlockMembership = nodeMapEntryIt->second.oldBlockMembership;
  std::vector<stk::mesh::Part*> transitiveTriplet(3);

#ifdef PRINT_DEBUG
  stk::mesh::EntityId refNodeId = nodeMapEntryIt->second.oldNodeId;
  info.print_debug_msg(3) << "fill_transitive_reconnect_node_info: ref node id: " << refNodeId << " current entity: " << bulk.identifier(currentEntity) << std::endl;
#endif

  for(const BlockPair& blockPair : blockPairsToReconnect) {
    bool foundFirst = std::binary_search(oldBlockMembership.begin(), oldBlockMembership.end(), blockPair.first, partLess);
    bool foundSecond = std::binary_search(oldBlockMembership.begin(), oldBlockMembership.end(), blockPair.second, partLess);

    if(foundFirst && foundSecond) {
      blockPairs.push_back(blockPair);
    }
  }
  if(!blockPairs.empty()) {
    for(unsigned i = 0; i < blockPairs.size()-1; i++) {
      for(unsigned j = i+1; j < blockPairs.size(); j++) {
        BlockPairComplement complement = get_block_pair_complement(blockPairs[i], blockPairs[j]);

        if(complement.commonPart != nullptr) {
          transitiveTriplet[0] = complement.blockPair.first;
          transitiveTriplet[1] = complement.blockPair.second;
          transitiveTriplet[2] = complement.commonPart;
#ifdef PRINT_DEBUG
          info.print_debug_msg(3) << "Considering triplet: {" << transitiveTriplet[0]->name() << "," << transitiveTriplet[1]->name() << "," << transitiveTriplet[2]->name() << "}" << std::endl;
#endif
          func(transitiveTriplet, nodeMapEntryIt, currentEntity);
        }
      }
    }
  }
}

template <typename Functor>
void update_procs_for_transitive_reconnection(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                              NodeMapType::iterator it, stk::mesh::Entity currentEntity, LinkInfo& info, Functor transitiveFunc)
{
  traverse_transitive_relations(bulk, blockPairsToReconnect, it, currentEntity, info, transitiveFunc);
}

void initialize_reconnect_node_id_info(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
{
  for(auto nodeMapEntryIt = info.preservedNodeMap.begin(); nodeMapEntryIt != info.preservedNodeMap.end(); ++nodeMapEntryIt) {
    stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
    auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
    if(reconnectMapIter == info.reconnectMap.end()) {
      ReconnectNodeInfo newReconnectInfo;
      newReconnectInfo.reconnectNodeId = referenceNodeId;
      info.reconnectMap.insert(std::make_pair(referenceNodeId, newReconnectInfo));
    }
  }

  for(auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    stk::mesh::EntityId id = stk::mesh::InvalidEntityId;
    stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
    stk::mesh::Entity currentEntity = nodeMapEntryIt->second.boundaryNode;
    stk::mesh::Entity referenceNode = bulk.get_entity(stk::topology::NODE_RANK, referenceNodeId);

    for(const BlockPair& blockPair : blockPairsToReconnect) {
      const stk::mesh::Part* otherBlock = nullptr;

      if(bulk.bucket(currentEntity).member(*blockPair.first)) {
        otherBlock = blockPair.second;
      }
      else if(bulk.bucket(currentEntity).member(*blockPair.second)) {
        otherBlock = blockPair.first;
      }
      if(otherBlock != nullptr) {

        if(bulk.is_valid(referenceNode)) {
          if(bulk.bucket(referenceNode).member(*otherBlock)) {
            id = std::min(id, referenceNodeId);
          }
        }
      }
    }
    auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
    if(reconnectMapIter != info.reconnectMap.end()) {
      reconnectMapIter->second.reconnectNodeId = std::min(reconnectMapIter->second.reconnectNodeId, id);
    } else {
      ReconnectNodeInfo newReconnectInfo;
      newReconnectInfo.reconnectNodeId = id;
      info.reconnectMap.insert(std::make_pair(referenceNodeId, newReconnectInfo));
    }
  }
}

void fill_direct_reconnect_node_info(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                     NodeMapType::iterator nodeMapEntryIt, stk::mesh::Entity currentEntity,
                                     LinkInfo& info)
{
  for(BlockPair blockPair : blockPairsToReconnect) {
    const stk::mesh::Part & srcBlock = *blockPair.first;
    const stk::mesh::Part & blockToReconnect = *blockPair.second;

    bool isInEitherBlockPair = bulk.bucket(currentEntity).member(srcBlock) || bulk.bucket(currentEntity).member(blockToReconnect);
    bool shouldBeReconnected = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, info) ||
                               should_be_reconnected(nodeMapEntryIt, srcBlock, blockToReconnect, info);

    if(isInEitherBlockPair && shouldBeReconnected) {
      const DisconnectGroup& disconnectGroup = nodeMapEntryIt->first.disconnectedGroup;
      stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
      ReconnectNodeInfo* reconnectInfo = nullptr;
#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << "ref Id: " << referenceNodeId << " currentEntityId: " << bulk.identifier(currentEntity) << std::endl;
#endif
      auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
      ThrowRequire(reconnectMapIter != info.reconnectMap.end());
      reconnectInfo = &reconnectMapIter->second;
#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << "found ref Id with reconnect node id: " << reconnectInfo->reconnectNodeId << std::endl;
#endif

      reconnectInfo->reconnectNodeId = std::min(bulk.identifier(currentEntity), reconnectInfo->reconnectNodeId);
#ifdef PRINT_DEBUG
      info.print_debug_msg(3) << " setting reconnectNodeId to " << reconnectInfo->reconnectNodeId << " for current entity " << bulk.identifier(currentEntity) << std::endl;;
#endif
      reconnectInfo->relatedNodes.push_back(currentEntity);

      std::vector<int> srcOwners = disconnectGroup.get_sharing_procs(srcBlock);
      std::vector<int> destOwners = disconnectGroup.get_sharing_procs(blockToReconnect);

      insert_uniquely_reconnect_info(bulk.parallel_rank(), srcOwners, *reconnectInfo);
      insert_uniquely_reconnect_info(bulk.parallel_rank(), destOwners, *reconnectInfo);
    }
  }
}

template <typename Functor>
void fill_reconnect_node_info(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                              NodeMapType::iterator it, stk::mesh::Entity currentEntity,
                              LinkInfo& info, Functor transitiveFunc)
{
  fill_direct_reconnect_node_info(bulk, blockPairsToReconnect, it, currentEntity, info);


#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << " using currentEntity with id: " << bulk.identifier(currentEntity) << " for fill_transitive" << std::endl;
#endif
  traverse_transitive_relations(bulk, blockPairsToReconnect, it, currentEntity, info, transitiveFunc);
#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << " using node with id: " << bulk.identifier(nodeMapEntryIt->first.disconnectedGroup.get_node()) << " for fill_transitive" << std::endl;
#endif
  traverse_transitive_relations(bulk, blockPairsToReconnect, it, it->first.disconnectedGroup.get_node(), info, transitiveFunc);
}

void determine_local_reconnect_node_id(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                       LinkInfo& info)
{
  auto fill_reconnect_node_info_func =
     [&](const stk::mesh::PartVector& transitiveBlockList_, NodeMapType::iterator it, stk::mesh::Entity currentEntity_)
     {
       auto reconnectMapIter = info.reconnectMap.find(it->second.oldNodeId);
       if(reconnectMapIter == info.reconnectMap.end()) { return; }
       ReconnectNodeInfo* reconnectInfo = &reconnectMapIter->second;

#ifdef PRINT_DEBUG
       info_.print_debug_msg(3) << "setting transitive reconnect node id: original: " << reconnectInfo->reconnectNodeId
           << " is now set to " << std::min(bulk.identifier(currentEntity_), reconnectInfo->reconnectNodeId) << std::endl;
#endif
       reconnectInfo->reconnectNodeId = std::min(bulk.identifier(currentEntity_), reconnectInfo->reconnectNodeId);
     };

  initialize_reconnect_node_id_info(bulk, blockPairsToReconnect, info);

#ifdef PRINT_DEBUG
  info.print_debug_msg(4) << "checking cloned map for determining direct node reconnect" << std::endl;
#endif
  for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, it->second.boundaryNode, info, fill_reconnect_node_info_func);
  }
#ifdef PRINT_DEBUG
  info.print_debug_msg(4) << "checking preserved map for determining direct node reconnect" << std::endl;
#endif
  for(auto it = info.preservedNodeMap.begin(); it != info.preservedNodeMap.end(); ++it) {
    const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
    stk::mesh::Entity currNode = disconnectGroup.get_node();
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, currNode, info, fill_reconnect_node_info_func);
  }
}

void update_reconnect_node_sharing(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                   LinkInfo& info)
{
  auto update_sharing_info_func =
      [&](const stk::mesh::PartVector& transitiveBlockList, NodeMapType::iterator it, stk::mesh::Entity currentEntity_)
      {
         const DisconnectGroup& disconnectedGroup = it->first.disconnectedGroup;
         stk::mesh::Entity boundaryNode = it->second.boundaryNode;

#ifdef PRINT_DEBUG
         info.print_debug_msg(1) << "checking ref id for " << it->second.oldNodeId << std::endl;
         info.print_debug_msg(1) << "update transitive reconnect id. first: " << transitiveBlockList[0]->name()
                          << " second: " << transitiveBlockList[1]->name() << " third: " << transitiveBlockList[2]->name()
                          << " boundary node id: " << bulk.identifier(boundaryNode)
                          << " currentEntity id: " << bulk.identifier(currentEntity_) << std::endl;
#endif
         auto reconnectMapIter = info.reconnectMap.find(it->second.oldNodeId);
         if(reconnectMapIter == info.reconnectMap.end()) { return; }
         ReconnectNodeInfo* reconnectInfo = &reconnectMapIter->second;
         std::vector<int> sharingProcs;

         for(const stk::mesh::Part* blockPart : transitiveBlockList) {
           sharingProcs = get_node_sharing_for_restoration(bulk, blockPart, disconnectedGroup, boundaryNode, info);

           for(int proc : sharingProcs) {
             stk::util::insert_keep_sorted_and_unique(proc, reconnectInfo->reconnectProcs);
           }
         }
      };
#ifdef PRINT_DEBUG
  info.print_debug_msg(4) << "checking cloned map for determining transitive node reconnect" << std::endl;
#endif
  for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
    update_procs_for_transitive_reconnection(bulk, blockPairsToReconnect, it, it->second.boundaryNode, info, update_sharing_info_func);
  }
#ifdef PRINT_DEBUG
  info.print_debug_msg(4) << "checking preserved map for determining transitive node reconnect" << std::endl;
#endif
  for(auto it = info.preservedNodeMap.begin(); it != info.preservedNodeMap.end(); ++it) {
    const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
    stk::mesh::Entity currNode = disconnectGroup.get_node();
    update_procs_for_transitive_reconnection(bulk, blockPairsToReconnect, it, currNode, info, update_sharing_info_func);
  }
}

const DisconnectGroup* get_group_for_block(NodeMapType& nodeMap, const stk::mesh::Part& blockPart, stk::mesh::Entity node)
{
  const DisconnectGroup* group = nullptr;

  for (NodeMapType::iterator nodeMapEntryIt = nodeMap.begin(); nodeMapEntryIt != nodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectGroup = nodeMapEntryIt->first.disconnectedGroup;
    const stk::mesh::ConstPartVector& parts = disconnectGroup.get_parts();

    if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == blockPart.mesh_meta_data_ordinal() &&
        node == disconnectGroup.get_node()) {
      group = &disconnectGroup;
    }
  }
  return group;
}

void reconnect_block_pair(stk::mesh::BulkData& bulk, const BlockPair & blockPair, impl::LinkInfo& info)
{
  const stk::mesh::Part & srcBlock = *blockPair.first;
  const stk::mesh::Part & blockToReconnect = *blockPair.second;

#ifdef PRINT_DEBUG
  info.print_debug_msg(2) << "considering block pair: " << blockPair.first->name() << " and " <<
      blockPair.second->name() << std::endl;

  info.print_debug_msg(2) << "Checking cloned list" << std::endl;
#endif
  for (NodeMapType::iterator nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    bool shouldReconnect = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, info) ||
                           should_be_reconnected(nodeMapEntryIt, srcBlock, blockToReconnect, info);

    if(shouldReconnect) {
      impl::reconnect_elements(bulk, blockPair, nodeMapEntryIt->first, nodeMapEntryIt->second, info);
    }
  }

#ifdef PRINT_DEBUG
  info.print_debug_msg(2) << "Checking preserved list" << std::endl;
#endif
  for (NodeMapType::iterator nodeMapEntryIt = info.preservedNodeMap.begin(); nodeMapEntryIt != info.preservedNodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    stk::mesh::Entity preservedNode = disconnectedGroup.get_node();

    ThrowRequire(bulk.is_valid(preservedNode));
    bool isInEitherBlockPair = bulk.bucket(preservedNode).member(srcBlock) || bulk.bucket(preservedNode).member(blockToReconnect);
    bool shouldReconnect = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, info);

    if(isInEitherBlockPair && shouldReconnect) {
#ifdef PRINT_DEBUG
      info.print_debug_msg(2) << "restoring node sharing" << std::endl;
#endif
      impl::restore_node_sharing(bulk, nodeMapEntryIt->second.oldNodeId, disconnectedGroup.get_node(), info);
    }
  }
}

void update_node_id(stk::mesh::EntityId newNodeId, int proc, LinkInfo& info, const DisconnectGroup& group) {
  ThrowRequire(!group.get_parts().empty());
  NodeMapKey nodeMapKey(group.get_node(), group);
  const auto & nodeMapIt = info.clonedNodeMap.find(nodeMapKey);

#ifdef PRINT_DEBUG
  info.print_debug_msg(3) << " update_node_id called" << std::endl;
#endif

  if (nodeMapIt != info.clonedNodeMap.end()) {
    NodeMapValue & newNodeData = nodeMapIt->second;
    newNodeData.newNodeId = std::min(newNodeData.newNodeId, newNodeId);
    newNodeData.oldNodeId = group.get_node_id();
    newNodeData.sharingProcs.push_back(proc);
#ifdef PRINT_DEBUG
    info.print_debug_msg(3) << "globally assigning id " << newNodeData.newNodeId << " to " << newNodeData.oldNodeId << std::endl;
#endif
  }
}

void clean_up_aura(stk::mesh::BulkData& bulk, LinkInfo& info)
{
  stk::mesh::EntityVector allNodes;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), allNodes);

  for(stk::mesh::Entity node : allNodes) {
    unsigned numElems = bulk.num_connectivity(node, stk::topology::ELEMENT_RANK);
    const stk::mesh::Entity* elems = bulk.begin(node, stk::topology::ELEMENT_RANK);
    int numLocallyOwnedElems = 0;

    for(unsigned i = 0; i < numElems; i++) {
      stk::mesh::Entity elem = elems[i];
      if(bulk.bucket(elem).owned()) {
        numLocallyOwnedElems++;
      }
    }

    if(numLocallyOwnedElems == 0) {
      for(unsigned j = 0; j < numElems; j++) {
        stk::mesh::Entity elem = elems[j];
#ifdef PRINT_DEBUG
        info.print_debug_msg(2) << "Destroying element: " << bulk.identifier(elem) << std::endl;;
#endif
        bulk.destroy_entity(elem);
      }
#ifdef PRINT_DEBUG
      info.print_debug_msg(2) << "Destroying node: " << bulk.identifier(node) << std::endl;
#endif
      bulk.destroy_entity(node);
    }
  }
}

void disconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToDisconnect,
                            LinkInfo& info)
{
  if(blockPairsToDisconnect.empty()) {
    stk::RuntimeWarningP0() << "No block pairs to disconnect" << std::endl;
  }

  bulk.modification_begin();

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    std::cout << "Adding nodes for disconnect" << std::endl;
  }

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
#ifdef PRINT_DEBUG
    info.print_debug_msg(1) << "First " << blockPairsToDisconnect[i].first->name() << " second " <<
         blockPairsToDisconnect[i].second->name() << std::endl;
#endif
    add_nodes_to_disconnect(bulk, blockPairsToDisconnect[i], info);
  }

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Creating new duplicate node IDs" << std::endl;
  }
  create_new_duplicate_node_IDs(bulk, info);

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Communicating shared node info" << std::endl;
  }
  communicate_shared_node_information(bulk, info);

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Disconnecting elements" << std::endl;
  }
  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    disconnect_elements(bulk, blockPairsToDisconnect[i], info);
  }

  clean_up_aura(bulk, info);

  info.flush();

  bulk.modification_end();

  if (bulk.has_face_adjacent_element_graph()) {
    bulk.delete_face_adjacent_element_graph();
    bulk.initialize_face_adjacent_element_graph();
  }
}

void reconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                           LinkInfo& info)
{
  bulk.modification_begin();

  sanitize_node_map(info.preservedNodeMap, info);

  info.flush();

  determine_local_reconnect_node_id(bulk, blockPairsToReconnect, info);

  info.flush();

  update_reconnect_node_sharing(bulk, blockPairsToReconnect, info);

  info.flush();

  communicate_reconnect_node_information(bulk, info);

  info.flush();

  for(const BlockPair & blockPair : blockPairsToReconnect) {
    reconnect_block_pair(bulk, blockPair, info);
  }
  info.flush();

  clean_up_aura(bulk, info);

  bulk.modification_end();
}
}
}
}
