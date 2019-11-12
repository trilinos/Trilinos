#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>
#include <map>

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

void remove_from_node_map(NodeMapType& nodeMap, DisconnectGroup& group, std::ostringstream& os)
{
  for(auto it = nodeMap.begin(); it != nodeMap.end(); ) {
    if(it->first.disconnectedGroup.get_node() == group.get_node()) {
      nodeMap.erase(it++);
      os << "P" << group.get_bulk().parallel_rank() << ": " << "removing node " << group.get_bulk().identifier(group.get_node()) << std::endl;
    } else {
      ++it;
    }
  }
}

void create_new_node_map_entry(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node, const BlockPairType& blockPair, LinkInfo& info)
{
  const stk::mesh::Part & firstBlock  = *blockPair.first;
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
  info.os << "P" << bulk.parallel_rank() << ": blockpairs " << firstBlock.name() << " " << secondBlock.name()
              << " needToClone: " << needToCloneNode << " node id: " << bulk.identifier(node) << std::endl;
  add_to_sharing_lookup(bulk, node, info.sharedInfo);
  if (needToCloneNode && (bulk.state(node) != stk::mesh::Created)) {
    DisconnectGroup group(bulk, &secondBlock, node);
    info.clonedNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
    info.os << "P" << bulk.parallel_rank() << ": " << "cloning node " << bulk.identifier(node) << " group id " <<
        group.id() << " part name: " << secondBlock.name() << std::endl;
  } else {
    DisconnectGroup group(bulk, &secondBlock, node);
    info.preservedNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
    info.os << "P" << bulk.parallel_rank() << ": " << "skipping node " << bulk.identifier(node) << " group id " <<
        group.id() << " part name: " << secondBlock.name() << std::endl;
  }
}

void add_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                             const BlockPairType & blockPair,
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
    info.os << "P" << bulk.parallel_rank() << ": locally assigning new id " << newNodeIDs[newNodeIdx-1] << " to node " <<
        bulk.identifier(nodeMapEntry.first.parentNode) << std::endl;
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
        group.pack_group_info(procBuff, newNodeId, proc, info.os);
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
      group.unpack_group_info(procBuff, newNodeId, proc, info.os);
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

std::vector<BlockPairType> get_block_pairs_to_disconnect(const stk::mesh::BulkData & bulk)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<BlockPairType> blockPairsToDisconnect;
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

void disconnect_elements(stk::mesh::BulkData& bulk, const DisconnectGroup& group, LinkInfo& info)
{
  for (auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    auto& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if (disconnectedGroup == group && disconnectedGroup.is_active()) {
      const stk::mesh::Entity node = nodeMapEntryIt->first.parentNode;

      for (stk::mesh::Entity elem : group.get_entities()) {
        stk::mesh::EntityId newNodeId = nodeMapEntryIt->second.newNodeId;
        stk::mesh::Entity newNode = bulk.declare_node(newNodeId);
        bulk.copy_entity_fields(node, newNode);

        for (int sharingProc : nodeMapEntryIt->second.sharingProcs) {
          bulk.add_node_sharing(newNode, sharingProc);
        }
        update_disconnected_element_relation(bulk, node, newNode, elem, info);
      }
      disconnectedGroup.set_active(false);
    }
  }
}

void disconnect_elements(stk::mesh::BulkData & bulk,
                         const BlockPairType & blockPair,
                         LinkInfo& info)
{
  const stk::mesh::Part & blockToDisconnect = *blockPair.second;

  for (auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    const DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    const stk::mesh::ConstPartVector& parts = disconnectedGroup.get_parts();

    if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == blockToDisconnect.mesh_meta_data_ordinal()) {
      disconnect_elements(bulk, disconnectedGroup, info);
    }
  }
}


void restore_node_sharing(stk::mesh::BulkData& bulk, const BlockPairType & blockPair,
                          const DisconnectGroup& group, stk::mesh::Entity destNode, LinkInfo& info)
{
  ThrowRequire(bulk.is_valid(destNode));
  stk::mesh::EntityId nodeId = bulk.identifier(destNode);
  std::vector<int> srcOwners = group.get_sharing_procs(*blockPair.first);
  std::vector<int> destOwners = group.get_sharing_procs(*blockPair.second);
  std::vector<int> commonOwners;

  info.os << "P" << bulk.parallel_rank() << ": printing owners for srcOwners: " << std::endl;
  for(int proc : srcOwners) {
    info.os << "\tP" << proc;
  }
  info.os << std::endl;
  info.os << "P" << bulk.parallel_rank() << ": printing owners for destOwners: " << std::endl;
  for(int proc : destOwners) {
    info.os << "\tP" << proc;
  }
  info.os << std::endl;

  for(int owner : destOwners) {
    if(owner == bulk.parallel_rank()) { continue; }
    bulk.add_node_sharing(destNode, owner);
    info.os << "P" << bulk.parallel_rank() << ": sharing node "  << nodeId << " to proc " << owner << std::endl;
  }
  for(int owner : srcOwners) {
    if(owner == bulk.parallel_rank()) { continue; }
    bulk.add_node_sharing(destNode, owner);
    info.os << "P" << bulk.parallel_rank() << ": sharing node "  << nodeId << " to proc " << owner << std::endl;
  }
}

stk::mesh::EntityId get_reconnect_node_id(LinkInfo& info, stk::mesh::EntityId referenceId)
{
   auto reconnectMapIter = info.reconnectMap.find(referenceId);
   ThrowRequire(reconnectMapIter != info.reconnectMap.end());
   return reconnectMapIter->second.reconnectNodeId;
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

void reconnect_elements(stk::mesh::BulkData& bulk, const BlockPairType & blockPair,
                        const DisconnectGroup& group, LinkInfo& info)
{
  for (auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    auto& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if (disconnectedGroup == group) {
      stk::mesh::EntityId reconnectId = get_reconnect_node_id(info, nodeMapEntryIt->second.oldNodeId);
      stk::mesh::Entity reconnectNode = bulk.get_entity(stk::topology::NODE_RANK, reconnectId);

      if(!bulk.is_valid(reconnectNode)) {
        info.os << "P" << bulk.parallel_rank() << ": Redeclaring old node " << reconnectId << std::endl;
        reconnectNode = bulk.declare_node(reconnectId);
      }

      stk::mesh::EntityId newNodeId = nodeMapEntryIt->second.newNodeId;
      stk::mesh::Entity newNode = bulk.get_entity(stk::topology::NODE_RANK, newNodeId);
      ThrowRequire(bulk.is_valid(newNode));

      restore_node_sharing(bulk, blockPair, group, reconnectNode, info);

      info.os << "P" << bulk.parallel_rank() << ": Restoring node " << newNodeId << " to " << reconnectId  << std::endl;
      update_reconnected_element_relation(bulk, group, newNode, reconnectNode);
    }
  }
}

struct BlockPairComplement {
  stk::mesh::Part* commonPart;
  BlockPairType blockPair;

  BlockPairComplement() :
    commonPart(nullptr),
    blockPair(BlockPairType(nullptr, nullptr)) {}
};

BlockPairComplement get_block_pair_complement(const BlockPairType blockPair1, const BlockPairType blockPair2)
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

void fix_indirect_node_sharing_for_block_pair(stk::mesh::BulkData& bulk, NodeMapType& nodeMap,
                                              const BlockPairComplement& complement, LinkInfo& info)
{
  for(NodeMapType::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it) {

    const DisconnectGroup& disconnectedGroup = it->first.disconnectedGroup;
    const stk::mesh::PartVector& oldBlockMembership = it->second.oldBlockMembership;
    bool foundFirstBlock = std::binary_search(oldBlockMembership.begin(), oldBlockMembership.end(), complement.blockPair.first, stk::mesh::PartLess());
    bool foundSecondBlock = std::binary_search(oldBlockMembership.begin(), oldBlockMembership.end(), complement.blockPair.second, stk::mesh::PartLess());
    bool foundCommonBlock = std::binary_search(oldBlockMembership.begin(), oldBlockMembership.end(), complement.commonPart, stk::mesh::PartLess());

    if(foundCommonBlock && foundFirstBlock && foundSecondBlock) {
      stk::mesh::EntityId oldNodeId = it->second.oldNodeId;
      stk::mesh::Entity oldNode = bulk.get_entity(stk::topology::NODE_RANK, oldNodeId);
      ThrowRequire(bulk.is_valid(oldNode));
      info.os << "P" << bulk.parallel_rank() << ": Fixing indirect node sharing for node " << oldNodeId
          << " between " << complement.blockPair.first->name() << " and " << complement.blockPair.second->name() << std::endl;
      restore_node_sharing(bulk, complement.blockPair, disconnectedGroup, oldNode, info);
    }
  }
}

void fix_indirect_node_sharing(stk::mesh::BulkData& bulk, const std::vector<BlockPairType>& blockPairsToReconnect, LinkInfo& info) {

  if(blockPairsToReconnect.empty()) { return; }

  for(unsigned i = 0; i < blockPairsToReconnect.size()-1; i++) {
    for(unsigned j = i+1; j < blockPairsToReconnect.size(); j++) {
      BlockPairComplement complement = get_block_pair_complement(blockPairsToReconnect[i], blockPairsToReconnect[j]);

      if(complement.commonPart != nullptr) {
        fix_indirect_node_sharing_for_block_pair(bulk, info.clonedNodeMap, complement, info);
        fix_indirect_node_sharing_for_block_pair(bulk, info.preservedNodeMap, complement, info);
      }
    }
  }
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
                           const stk::mesh::Part& srcBlock, std::ostringstream& os)
{
  const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
  const stk::mesh::PartVector& blockMembership = nodeMapEntryIt->second.oldBlockMembership;
  const stk::mesh::BulkData& bulk = disconnectedGroup.get_bulk();
  const stk::mesh::ConstPartVector& parts = disconnectedGroup.get_parts();
  bool isOriginalMember = false;

  os << "P" << bulk.parallel_rank() << ": parts size: " << parts.size() << " : name "
            << parts[0]->name() << " group id: " << disconnectedGroup.id() << std::endl;

  if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == blockToReconnect.mesh_meta_data_ordinal()) {
    isOriginalMember = std::binary_search(blockMembership.begin(), blockMembership.end(), &srcBlock, stk::mesh::PartLess());

    os << "P" << bulk.parallel_rank() << ": Checking to see if part " << srcBlock.name() << " is original member of old node " << nodeMapEntryIt->second.oldNodeId << " : isOriginalMember = "<< isOriginalMember << std::endl;
    for(stk::mesh::Part* part : blockMembership) {
      os << "\t\tP" << bulk.parallel_rank() << ": " << part->name() << std::endl;
    }
  }
  return isOriginalMember;
}

bool should_update_sharing_info(NodeMapType& nodeMap, const DisconnectGroup& group,
                                const stk::mesh::Part& blockToReconnect, std::ostringstream& os)
{
  os << "P" << group.get_bulk().parallel_rank() << ": Checking to see if to update sharing info for old node " << group.get_node_id() << std::endl;
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

        if(should_be_reconnected(it, *block2, *block1, os)) {
          os << "P" << group.get_bulk().parallel_rank() << " returning true" << std::endl;
          return true;
        }
      }
    }
  }
  os << "P" << group.get_bulk().parallel_rank() << " returning false" << std::endl;
  return false;
}

void sanitize_node_map(NodeMapType& nodeMap, std::ostringstream& os)
{
  for (NodeMapType::iterator nodeMapEntryIt = nodeMap.begin(); nodeMapEntryIt != nodeMap.end();) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if(!disconnectedGroup.get_bulk().is_valid(disconnectedGroup.get_node())) {
      os << "P" << disconnectedGroup.get_bulk().parallel_rank() << " sanitizing node id "
          << nodeMapEntryIt->second.oldNodeId << std::endl;
      nodeMapEntryIt = nodeMap.erase(nodeMapEntryIt);
    } else {
      nodeMapEntryIt++;
    }
  }
}

void insert_uniquely_reconnect_info(const std::vector<int>& procs, ReconnectNodeInfo& reconnectInfo)
{
  for(unsigned i = 0; i < procs.size(); i++) {
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
      info.os << "P" << bulk.parallel_rank() << ": passing info about ref node: " << reconnectMapEntry.first << " to proc " << proc << std::endl;
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

      info.os << "P" << bulk.parallel_rank() << ": recved info about ref node: " << referenceNodeId << " from proc " << proc << std::endl;
      auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
      ThrowRequire(reconnectMapIter != info.reconnectMap.end());

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
}

ReconnectNodeInfo get_initialized_node_info_for_block_pair(stk::mesh::BulkData& bulk, const BlockPairType& blockPair,
                                                           stk::mesh::Entity currentEntity, stk::mesh::EntityId referenceNodeId,
                                                           LinkInfo& info)
{
  ReconnectNodeInfo reconnectInfo;
  const stk::mesh::Part* otherBlock = nullptr;

  if(bulk.bucket(currentEntity).member(*blockPair.first)) {
    otherBlock = blockPair.second;
  }
  else if(bulk.bucket(currentEntity).member(*blockPair.second)) {
    otherBlock = blockPair.first;
  }

  if(otherBlock != nullptr) {
    stk::mesh::Entity referenceNode = bulk.get_entity(stk::topology::NODE_RANK, referenceNodeId);

    if(bulk.is_valid(referenceNode)) {
      if(bulk.bucket(referenceNode).member(*otherBlock)) {
        info.os << "P" << bulk.parallel_rank() << " Initialization1: entityid: " << bulk.identifier(currentEntity) << " reconnectNodeId " << reconnectInfo.reconnectNodeId << " set to ref id " << referenceNodeId << std::endl;
        reconnectInfo.reconnectNodeId = referenceNodeId;
      }
    }

    for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
      const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
      const stk::mesh::ConstPartVector& parts = disconnectGroup.get_parts();

      if (parts.size() == 1u && parts[0]->mesh_meta_data_ordinal() == otherBlock->mesh_meta_data_ordinal() &&
          it->second.oldNodeId == referenceNodeId) {
        info.os << "P" << bulk.parallel_rank() << " Initialization2: entityid: " << bulk.identifier(currentEntity) << " reconnectNodeId " << reconnectInfo.reconnectNodeId << " set to "
            << std::min(reconnectInfo.reconnectNodeId, it->second.newNodeId) << std::endl;
        reconnectInfo.reconnectNodeId = std::min(reconnectInfo.reconnectNodeId, it->second.newNodeId);
      }
    }
  }

  return reconnectInfo;
}

void fill_reconnect_node_info(stk::mesh::BulkData& bulk, const std::vector<BlockPairType>& blockPairsToReconnect,
                              NodeMapType::iterator nodeMapEntryIt, stk::mesh::Entity currentEntity,
                              LinkInfo& info)
{
  for(BlockPairType blockPair : blockPairsToReconnect) {
    const stk::mesh::Part & srcBlock = *blockPair.first;
    const stk::mesh::Part & blockToReconnect = *blockPair.second;

    bool isInEitherBlockPair = bulk.bucket(currentEntity).member(srcBlock) || bulk.bucket(currentEntity).member(blockToReconnect);
    bool shouldBeReconnected = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, info.os) ||
                               should_be_reconnected(nodeMapEntryIt, srcBlock, blockToReconnect, info.os);

    if(isInEitherBlockPair && shouldBeReconnected) {
      const DisconnectGroup& disconnectGroup = nodeMapEntryIt->first.disconnectedGroup;
      stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
      ReconnectNodeInfo* reconnectInfo = nullptr;
      info.os << "P" << bulk.parallel_rank() << ": ref Id: " << referenceNodeId << " currentEntityId: " << bulk.identifier(currentEntity) << std::endl;
      auto reconnectMapIter = info.reconnectMap.find(referenceNodeId);
      if(reconnectMapIter != info.reconnectMap.end()) {
        reconnectInfo = &reconnectMapIter->second;
        info.os << "P" << bulk.parallel_rank() << ": found ref Id with reconnect node id: " << reconnectInfo->reconnectNodeId << std::endl;
      } else {
        info.os << "P" << bulk.parallel_rank() << ": did not find ref Id "  << std::endl;
        ReconnectNodeInfo newReconnectInfo = get_initialized_node_info_for_block_pair(bulk, blockPair, currentEntity, referenceNodeId, info);
        auto entryStatus = info.reconnectMap.insert(std::make_pair(referenceNodeId, newReconnectInfo));
        reconnectInfo = &entryStatus.first->second;
      }

      reconnectInfo->reconnectNodeId = std::min(bulk.identifier(currentEntity), reconnectInfo->reconnectNodeId);
      reconnectInfo->relatedNodes.push_back(currentEntity);

      std::vector<int> srcOwners = disconnectGroup.get_sharing_procs(srcBlock);
      std::vector<int> destOwners = disconnectGroup.get_sharing_procs(blockToReconnect);

      insert_uniquely_reconnect_info(srcOwners, *reconnectInfo);
      insert_uniquely_reconnect_info(destOwners, *reconnectInfo);
    }
  }

}

void determine_reconnect_node_id(stk::mesh::BulkData& bulk, const std::vector<BlockPairType>& blockPairsToReconnect,
                                 LinkInfo& info)
{
  info.os << "P" << bulk.parallel_rank() << ": checking cloned map for determining node reconnect" << std::endl;
  for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
    stk::mesh::EntityId currId = it->second.newNodeId;
    stk::mesh::Entity currNode = bulk.get_entity(stk::topology::NODE_RANK, currId);
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, currNode, info);
  }
  info.os << "P" << bulk.parallel_rank() << ": checking preserved map for determining node reconnect" << std::endl;
  for(auto it = info.preservedNodeMap.begin(); it != info.preservedNodeMap.end(); ++it) {
    const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
    stk::mesh::Entity currNode = disconnectGroup.get_node();
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, currNode, info);
  }

  info.flush();

  communicate_reconnect_node_information(bulk, info);

  info.os << "P" << bulk.parallel_rank() << ": reconnect node info " << std::endl;

  for(auto iter : info.reconnectMap) {
    info.os << "\tRef: " << iter.first << " reconnectId: "<< iter.second.reconnectNodeId << std::endl;
    for(stk::mesh::Entity node : iter.second.relatedNodes) {
      info.os << "\t\tNodes: " << bulk.identifier(node);
    }
    info.os << std::endl;
  }
  info.os << std::endl;
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

void reconnect_block_pair(stk::mesh::BulkData& bulk, const impl::BlockPairType & blockPair, impl::LinkInfo& info)
{
  std::ostringstream &os = info.os;

  const stk::mesh::Part & srcBlock = *blockPair.first;
  const stk::mesh::Part & blockToReconnect = *blockPair.second;

  os << "P" << bulk.parallel_rank() << ": considering block pair: " << blockPair.first->name() << " and " <<
      blockPair.second->name() << std::endl;

  os << "P" << bulk.parallel_rank() << ": Checking cloned list" << std::endl;
  for (NodeMapType::iterator nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    bool shouldReconnect = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, os) ||
                           should_be_reconnected(nodeMapEntryIt, srcBlock, blockToReconnect, os);

    if(shouldReconnect) {
      impl::reconnect_elements(bulk, blockPair, disconnectedGroup, info);
    }
  }

  os << "P" << bulk.parallel_rank() << ": Checking preserved list" << std::endl;
  for (NodeMapType::iterator nodeMapEntryIt = info.preservedNodeMap.begin(); nodeMapEntryIt != info.preservedNodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    stk::mesh::Entity preservedNode = disconnectedGroup.get_node();

    ThrowRequire(bulk.is_valid(preservedNode));
    bool isInEitherBlockPair = bulk.bucket(preservedNode).member(srcBlock) || bulk.bucket(preservedNode).member(blockToReconnect);
    bool shouldReconnect = should_be_reconnected(nodeMapEntryIt, blockToReconnect, srcBlock, os);

    if(isInEitherBlockPair && shouldReconnect) {
      os << "P" << bulk.parallel_rank() << ": restoring node sharing" << std::endl;
      impl::restore_node_sharing(bulk, blockPair, disconnectedGroup, disconnectedGroup.get_node(), info);
    }
  }
}

void update_node_id(stk::mesh::EntityId newNodeId, int proc, LinkInfo& info, const DisconnectGroup& group) {
  ThrowRequire(!group.get_parts().empty());
  NodeMapKey nodeMapKey(group.get_node(), group);
  const auto & nodeMapIt = info.clonedNodeMap.find(nodeMapKey);

  info.os << "P" << group.get_bulk().parallel_rank() << " update_node_id called" << std::endl;

  if (nodeMapIt != info.clonedNodeMap.end()) {
    NodeMapValue & newNodeData = nodeMapIt->second;
    newNodeData.newNodeId = std::min(newNodeData.newNodeId, newNodeId);
    newNodeData.oldNodeId = group.get_node_id();
    newNodeData.sharingProcs.push_back(proc);
    info.os << "P" << group.get_bulk().parallel_rank() << ": globally assigning id " << newNodeData.newNodeId << " to " << newNodeData.oldNodeId << std::endl;
  }
}

}
}
}
