#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/DestroyRelations.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/SideSetUtil.hpp"
#include "stk_mesh/base/SkinMeshUtil.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocks.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/environment/RuntimeWarning.hpp"
#include <algorithm>
#include <map>

#include <string>
#include <sstream>

namespace stk {
namespace tools {
namespace impl {

class ReconnectGroup {
public:
  ReconnectGroup(const DisconnectGroup& disconnectGroup, const NodeMapValue& nodeMapValue,
                 const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
  {
    set_block_pairs(disconnectGroup, nodeMapValue, blockPairsToReconnect, info);
    set_parts();
    set_id();
  }

  ReconnectGroup(const ReconnectGroup& group)
  {
    groupBlockPairs = group.groupBlockPairs;
    groupBlocks = group.groupBlocks;
    groupId = group.groupId;
  }

  unsigned get_id() const { return groupId; }
  const stk::mesh::PartVector& get_blocks() const { return groupBlocks; }
  const std::vector<BlockPair>& get_block_pairs() const { return groupBlockPairs; }

private:
  std::vector<BlockPair> groupBlockPairs;
  stk::mesh::PartVector groupBlocks;
  unsigned groupId = std::numeric_limits<unsigned>::max();

  ReconnectGroup() {}

  void set_id()
  {
    groupId = std::numeric_limits<unsigned>::max();

    for(const BlockPair& blockPair : groupBlockPairs) {
      groupId = std::min(groupId, blockPair.first->mesh_meta_data_ordinal());
      groupId = std::min(groupId, blockPair.second->mesh_meta_data_ordinal());
    }
  }

  void set_parts()
  {
    stk::mesh::PartLess compare;

    for(const BlockPair& blockPair : groupBlockPairs) {
      stk::util::insert_keep_sorted_and_unique(blockPair.first, groupBlocks, compare);
      stk::util::insert_keep_sorted_and_unique(blockPair.second, groupBlocks, compare);
    }
  }

  void set_block_pairs(const DisconnectGroup& disconnectGroup, const NodeMapValue& nodeMapValue,
                       const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
  {
    for(const BlockPair& blockPair : blockPairsToReconnect) {
      bool canBeReconnected = can_be_reconnected(disconnectGroup, nodeMapValue, blockPair, info);

      if(canBeReconnected) {
        groupBlockPairs.push_back(blockPair);
      }
    }
  }
};

ReconnectGroup get_group_for_reconnect_node(const DisconnectGroup& disconnectGroup, const NodeMapValue& nodeMapValue,
                                            const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
{
  ReconnectGroup reconnectGroup(disconnectGroup, nodeMapValue, blockPairsToReconnect, info);

  return reconnectGroup;
}

bool is_in_disconnected_block(const stk::mesh::BulkData& bulk, const stk::mesh::Entity elem, const BlockPair& blockPair)
{
  return bulk.bucket(elem).member(*blockPair.second);
}

bool is_in_disconnected_block(const stk::mesh::Bucket & elemBucket, const BlockPair& blockPair)
{
  return elemBucket.member(*blockPair.second);
}

bool is_in_preserved_block(const stk::mesh::BulkData& bulk, const stk::mesh::Entity elem, const BlockPair& blockPair)
{
  return bulk.bucket(elem).member(*blockPair.first);
}

bool is_in_preserved_block(const stk::mesh::Bucket & elemBucket, const BlockPair& blockPair)
{
  return elemBucket.member(*blockPair.first);
}

void insert_parts_uniquely(const stk::mesh::PartVector& fromParts, stk::mesh::PartVector& toParts)
{
  stk::mesh::PartLess compare;

  for(stk::mesh::Part* part : fromParts) {
    stk::util::insert_keep_sorted_and_unique(part, toParts, compare);
  }
}

void add_to_sharing_lookup(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, PreservedSharingInfo& info)
{
  stk::mesh::EntityId id = bulk.identifier(node);

  auto iter = info.find(id);
  if(iter == info.end()) {
    std::vector<int> commSharedProcs;
    bulk.comm_shared_procs(node, commSharedProcs);
    info.insert(std::make_pair(id, commSharedProcs));
  }
}

void create_new_node_map_entry(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node, const BlockPair& blockPair, LinkInfo& info)
{
  const stk::mesh::Entity * elems = bulk.begin_elements(node);
  const unsigned numElems = bulk.num_elements(node);
  bool needToCloneNode = false;
  for (unsigned i = 0; i < numElems; ++i) {
    const stk::mesh::Bucket & elemBucket = bulk.bucket(elems[i]);
    if (is_in_disconnected_block(elemBucket, blockPair) && elemBucket.owned()) {
      needToCloneNode = true;
      break;
    }
  }

  add_to_sharing_lookup(bulk, node, info.sharedInfo);
  if (needToCloneNode && (bulk.state(node) != stk::mesh::Created)) {
    DisconnectGroup group(bulk, blockPair, node);
    info.clonedNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
  } else {
    DisconnectGroup group(bulk, blockPair, node);
    info.originalNodeMap[NodeMapKey(node, group)] = NodeMapValue(bulk, node);
  }
}

std::vector<int> get_reconnect_procs(LinkInfo& info, const ReconnectMapKey& mapKey)
{
  auto reconnectMapIter = info.reconnectMap.find(mapKey);
  STK_ThrowRequire(reconnectMapIter != info.reconnectMap.end());
  return reconnectMapIter->second.reconnectProcs;
}

stk::mesh::EntityId get_reconnect_node_id(LinkInfo& info, const ReconnectMapKey& mapKey)
{
  auto reconnectMapIter = info.reconnectMap.find(mapKey);
  STK_ThrowRequire(reconnectMapIter != info.reconnectMap.end());
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
  STK_ThrowAssert(secondBlock.mesh_meta_data_ordinal() > firstBlock.mesh_meta_data_ordinal());

  stk::mesh::Selector boundaryBetweenBlocks = firstBlock & secondBlock & (meta.locally_owned_part() | meta.globally_shared_part());

  const stk::mesh::BucketVector & nodesOnBoundaryBetweenBlocks = bulk.get_buckets(stk::topology::NODE_RANK, boundaryBetweenBlocks);
  for (const stk::mesh::Bucket * bucket : nodesOnBoundaryBetweenBlocks) {
    for (const stk::mesh::Entity node : *bucket) {
      create_new_node_map_entry(bulk, node, blockPair, info);
    }
  }
}

void add_faces_to_disconnect(stk::mesh::BulkData& bulk,
                             BlockPair const& blockPair,
                             LinkInfo& info)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const stk::mesh::Part& firstBlock  = *blockPair.first;
  const stk::mesh::Part& secondBlock = *blockPair.second;

  auto firstBlockSurfaces = meta.get_surfaces_touched_by_block(blockPair.first);
  auto secondBlockSurfaces = meta.get_surfaces_touched_by_block(blockPair.second);

  auto check_is_in_blocks = [](auto& inBlock, auto& surfaces, auto partOrdinal) {
    for (auto surface : surfaces) {
      if (partOrdinal == surface->mesh_meta_data_ordinal()) {
        inBlock = true;
        break;
      }
    }};

  stk::mesh::Selector blockBoundary = firstBlock & secondBlock &
                                      (meta.locally_owned_part() | meta.globally_shared_part());
  const stk::mesh::BucketVector & onBlockBoundaryFaces = bulk.get_buckets(meta.side_rank(), blockBoundary);
  for (auto bucket : onBlockBoundaryFaces) {
    stk::mesh::PartVector bucketSidesetParts;

    const stk::mesh::PartVector& superset = bucket->supersets();
    for (stk::mesh::Part* part : superset) {
      if (stk::mesh::is_side_set(*part)) {
        bucketSidesetParts.push_back(part);
      }
    }

    for (auto face : *bucket) {
      auto elems = bulk.begin_elements(face);
      auto numElems = bulk.num_elements(face);
      STK_ThrowRequire(numElems == 2u);

      std::pair<bool, bool> isInBlocks = {false, false};
      for (auto part : bucketSidesetParts) {
        auto ordinal = part->mesh_meta_data_ordinal();

        check_is_in_blocks(isInBlocks.first, firstBlockSurfaces, ordinal);
        check_is_in_blocks(isInBlocks.second, secondBlockSurfaces, ordinal);
      }

      auto ordinals = bulk.begin_ordinals(face, stk::topology::ELEM_RANK);

      for (unsigned i = 0; i < numElems; ++i) {
        if (is_in_disconnected_block(bulk, elems[i], blockPair)) {
          InternalFaceInfo faceInfo(blockPair,
                                    elems[1u-i], ordinals[1u-i],
                                    elems[i], ordinals[i], face,
                                    isInBlocks);
          stk::util::insert_keep_sorted_and_unique(faceInfo, info.internalSides);
        }
      }
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
  }
}

void pack_shared_node_information(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, LinkInfo& info)
{
  for (const auto & nodeMapEntry : info.clonedNodeMap) {
    const stk::mesh::Entity node = nodeMapEntry.first.parentNode;
    if (!bulk.bucket(node).shared()) continue;
    const DisconnectGroup & group = nodeMapEntry.first.disconnectedGroup;
    STK_ThrowRequire(node == group.get_node());

    const stk::mesh::EntityId newNodeId = nodeMapEntry.second.newNodeId;
    bool needToCommunicate = group.needs_to_communicate();

    if (needToCommunicate) {
      std::vector<int> sharingProcs;
      bulk.comm_shared_procs(node, sharingProcs);
      for (const int proc : sharingProcs) {
        stk::CommBuffer& procBuff = commSparse.send_buffer(proc);
        group.pack_group_info(procBuff, newNodeId, proc);
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
      group.unpack_group_info(procBuff, newNodeId, proc);
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

bool update_disconnected_entity_relation(stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::Entity newNode,
                                         stk::mesh::Entity entity)
{
  bool updatedNode = false;
  unsigned numNodes = bulk.num_connectivity(entity, stk::topology::NODE_RANK);
  const stk::mesh::Entity * entityNodes = bulk.begin(entity, stk::topology::NODE_RANK);
  stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(entity, stk::topology::NODE_RANK);
  for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
    if (entityNodes[iNode] == node) {
      stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrdinals[iNode];
      bulk.destroy_relation(entity, node, nodeOrd);
      bulk.declare_relation(entity, newNode, nodeOrd);
      updatedNode = true;
    }
  }
  return updatedNode;
}

void update_internal_face_node_map(stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::Entity newNode,
                                   InternalFaceInfo& faceInfo)
{
  stk::mesh::Entity face = faceInfo.internalFace;
  unsigned numNodes = bulk.num_connectivity(face, stk::topology::NODE_RANK);
  const stk::mesh::Entity * entityNodes = bulk.begin(face, stk::topology::NODE_RANK);
  stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(face, stk::topology::NODE_RANK);
  for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
    if (entityNodes[iNode] == node) {
      faceInfo.nodeMap[std::make_pair(node,nodeOrdinals[iNode])] = newNode;
    }
  }
}

void update_internal_face_entity_relation(stk::mesh::BulkData& bulk, InternalFaceInfo& faceInfo)
{
  stk::mesh::Entity face = faceInfo.internalFace;

  for (auto iter = faceInfo.nodeMap.begin(); iter != faceInfo.nodeMap.end(); ++iter) {
    const auto& key = iter->first;
    const auto& newNode = iter->second;

    stk::mesh::Entity node = key.first;
    stk::mesh::ConnectivityOrdinal ordinal = key.second;

    bulk.destroy_relation(face, node, ordinal);
    bulk.declare_relation(face, newNode, ordinal);
  }
}

void disconnect_sub_rank(stk::mesh::BulkData& bulk, LinkInfo& info, stk::mesh::EntityRank rank,
                         stk::mesh::Entity newNode, const stk::mesh::Entity oldNode,
                         stk::mesh::Entity elem, stk::tools::BlockPair& blockPair)
{
  if (blockPair.is_valid() && is_in_preserved_block(bulk, elem, blockPair)) { return; }

  auto sideRank = bulk.mesh_meta_data().side_rank();
  const stk::mesh::Entity* subEntities = bulk.begin(elem, rank);
  unsigned numSubEntities = bulk.num_connectivity(elem, rank);
  for(unsigned i = 0; i < numSubEntities; i++) {
    auto entity = subEntities[i];

    if(sideRank == rank) {
      auto iter = std::lower_bound(info.internalSides.begin(), info.internalSides.end(), entity);
      if (!(iter == info.internalSides.end()) && !(entity < *iter)) {
        update_internal_face_node_map(bulk, oldNode, newNode, *iter);
        continue;
      }
    }

    update_disconnected_entity_relation(bulk, oldNode, newNode, entity);
  }
}

void disconnect_sub_ranks(stk::mesh::BulkData& bulk, LinkInfo& info,
                          stk::mesh::Entity newNode, const stk::mesh::Entity oldNode,
                          stk::mesh::Entity elem, stk::tools::BlockPair& blockPair)
{
  disconnect_sub_rank(bulk, info, bulk.mesh_meta_data().side_rank(), newNode, oldNode, elem, blockPair);

  if(bulk.mesh_meta_data().side_rank() != stk::topology::EDGE_RANK) {
    disconnect_sub_rank(bulk, info, stk::topology::EDGE_RANK, newNode, oldNode, elem, blockPair);
  }
}

void disconnect_elements(stk::mesh::BulkData& bulk, const NodeMapKey& key, NodeMapValue& value, LinkInfo& info)
{
  const DisconnectGroup& group = key.disconnectedGroup;
  auto blockPair = group.get_block_pair();

  if (group.is_active()) {
    const stk::mesh::Entity node = key.parentNode;
    stk::mesh::EntityId newNodeId = value.newNodeId;
    stk::mesh::Entity newNode = bulk.declare_node(newNodeId);
    value.boundaryNode = newNode;
    bulk.copy_entity_fields(node, newNode);

    for (stk::mesh::Entity elem : group.get_entities()) {
      for (int sharingProc : value.sharingProcs) {
        bulk.add_node_sharing(newNode, sharingProc);
      }
      bool updatedNode = update_disconnected_entity_relation(bulk, node, newNode, elem);

      if(updatedNode) {
        disconnect_sub_ranks(bulk, info, newNode, node, elem, blockPair);
      }

      if (bulk.num_elements(node) == 0 && !info.preserveOrphans) {
        bulk.destroy_entity(node);
      }
    }

    group.set_active(false);
  }
}

void disconnect_elements(stk::mesh::BulkData & bulk, const BlockPair & blockPair, LinkInfo& info)
{
  for (auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    const DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if (disconnectedGroup.has_block_pair() && disconnectedGroup.get_block_pair() == blockPair) {
      disconnect_elements(bulk, nodeMapEntryIt->first, nodeMapEntryIt->second, info);
    }
  }
}

void destroy_internal_face(stk::mesh::BulkData& bulk, const InternalFaceInfo& faceInfo)
{
  bulk.destroy_relation(faceInfo.elemInSecondBlock, faceInfo.internalFace, faceInfo.faceOrdinalInSecondBlock);
  bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);
  bulk.destroy_entity(faceInfo.internalFace, false);
}

void duplicate_and_destroy_internal_face(stk::mesh::BulkData& bulk,
                                         const InternalFaceInfo& faceInfo,
                                         const InternalFaceDisconnectState firstBlockState,
                                         const InternalFaceDisconnectState secondBlockState,
                                         stk::mesh::ConstPartVector& firstBlockSurfaces,
                                         stk::mesh::ConstPartVector& secondBlockSurfaces)
{
  bulk.destroy_relation(faceInfo.elemInSecondBlock, faceInfo.internalFace, faceInfo.faceOrdinalInSecondBlock);
  bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);

  if(secondBlockState == InternalFaceDisconnectState::KEEP) {
    stk::mesh::Entity newFace = bulk.declare_element_side(faceInfo.elemInSecondBlock,
                                                          faceInfo.faceOrdinalInSecondBlock,
                                                          stk::mesh::PartVector{faceInfo.blockPair.second});
    STK_ThrowRequire(newFace != faceInfo.internalFace);

    bulk.copy_entity_fields(faceInfo.internalFace, newFace);
    bulk.change_entity_parts(newFace, secondBlockSurfaces, firstBlockSurfaces);
  }

  if(firstBlockState == InternalFaceDisconnectState::KEEP) {
    stk::mesh::Entity newFace = bulk.declare_element_side(faceInfo.elemInFirstBlock,
                                                          faceInfo.faceOrdinalInFirstBlock,
                                                          stk::mesh::PartVector{faceInfo.blockPair.first});
    STK_ThrowRequire(newFace != faceInfo.internalFace);

    bulk.copy_entity_fields(faceInfo.internalFace, newFace);
    bulk.change_entity_parts(newFace, firstBlockSurfaces, secondBlockSurfaces);
  }

  bulk.destroy_entity(faceInfo.internalFace, false);
}

void connect_element_side_to_internal_face(stk::mesh::BulkData& bulk,
                                           stk::mesh::Entity elem,
                                           stk::mesh::ConnectivityOrdinal sideOrdinal,
                                           stk::mesh::Entity face,
                                           stk::mesh::ConstPartVector& addSurfaces,
                                           stk::mesh::ConstPartVector& removeSurfaces)
{
  stk::mesh::destroy_relations(bulk, face, stk::topology::ELEM_RANK);
  stk::mesh::destroy_relations(bulk, face, stk::topology::NODE_RANK);

  auto elemTopology = bulk.bucket(elem).topology();
  auto sideTopology = elemTopology.side_topology(sideOrdinal);

  stk::mesh::impl::connect_element_to_entity(bulk, elem, face, sideOrdinal, addSurfaces, sideTopology);
  bulk.change_entity_parts(face, addSurfaces, removeSurfaces);
}

void duplicate_non_owned_internal_face(stk::mesh::BulkData& bulk,
                                       const InternalFaceInfo& faceInfo,
                                       stk::mesh::ConstPartVector& firstBlockSurfaces,
                                       stk::mesh::ConstPartVector& secondBlockSurfaces)
{
  bulk.destroy_relation(faceInfo.elemInSecondBlock, faceInfo.internalFace, faceInfo.faceOrdinalInSecondBlock);
  bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);

  stk::mesh::Entity tempSoloSide = bulk.declare_solo_side(stk::mesh::PartVector{faceInfo.blockPair.second});
  STK_ThrowRequire(tempSoloSide != faceInfo.internalFace);
  bulk.copy_entity_fields(faceInfo.internalFace, tempSoloSide);

  bulk.destroy_entity(faceInfo.internalFace, false);

  stk::mesh::Entity newFaceOnFirstBlock = bulk.declare_element_side(faceInfo.elemInFirstBlock,
                                                                    faceInfo.faceOrdinalInFirstBlock,
                                                                    stk::mesh::PartVector{faceInfo.blockPair.first});
  STK_ThrowRequire(newFaceOnFirstBlock != tempSoloSide);

  bulk.copy_entity_fields(tempSoloSide, newFaceOnFirstBlock);
  bulk.change_entity_parts(newFaceOnFirstBlock, firstBlockSurfaces, secondBlockSurfaces);

  bulk.destroy_entity(tempSoloSide, false);
}

void detach_internal_face(stk::mesh::BulkData& bulk,
                          const InternalFaceInfo& faceInfo,
                          const InternalFaceDisconnectState firstBlockState,
                          const InternalFaceDisconnectState secondBlockState,
                          stk::mesh::ConstPartVector& firstBlockSurfaces,
                          stk::mesh::ConstPartVector& secondBlockSurfaces)
{
  stk::mesh::EntityId faceIdFromSecondBlock = stk::mesh::impl::side_id_formula(bulk.identifier(faceInfo.elemInSecondBlock),
                                                                               faceInfo.faceOrdinalInSecondBlock);

  if(bulk.identifier(faceInfo.internalFace) == faceIdFromSecondBlock) {
    if(secondBlockState == InternalFaceDisconnectState::KEEP) {
      bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);
    }

    if(firstBlockState == InternalFaceDisconnectState::KEEP) {
      duplicate_and_destroy_internal_face(bulk, faceInfo, firstBlockState, InternalFaceDisconnectState::NOOP,
                                          firstBlockSurfaces, secondBlockSurfaces);
    }
  } else {
    if(secondBlockState == InternalFaceDisconnectState::KEEP) {
      duplicate_and_destroy_internal_face(bulk, faceInfo, InternalFaceDisconnectState::NOOP, secondBlockState,
                                          firstBlockSurfaces, secondBlockSurfaces);
    }

    if(firstBlockState == InternalFaceDisconnectState::KEEP) {
      if(bulk.bucket(faceInfo.internalFace).owned()) {
        connect_element_side_to_internal_face(bulk, faceInfo.elemInFirstBlock, faceInfo.faceOrdinalInFirstBlock,
                                              faceInfo.internalFace, firstBlockSurfaces, secondBlockSurfaces);
      } else {
        duplicate_non_owned_internal_face(bulk, faceInfo, firstBlockSurfaces, secondBlockSurfaces);
      }
    }
  }
}

void detach_and_duplicate_internal_face(stk::mesh::BulkData& bulk,
                                        const InternalFaceInfo& faceInfo,
                                        const InternalFaceDisconnectState firstBlockState,
                                        const InternalFaceDisconnectState secondBlockState,
                                        stk::mesh::ConstPartVector& firstBlockSurfaces,
                                        stk::mesh::ConstPartVector& secondBlockSurfaces)
{
  stk::mesh::EntityId faceIdFromSecondBlock = stk::mesh::impl::side_id_formula(bulk.identifier(faceInfo.elemInSecondBlock),
                                                                               faceInfo.faceOrdinalInSecondBlock);

  if(bulk.identifier(faceInfo.internalFace) == faceIdFromSecondBlock) {
    bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);

    stk::mesh::Entity newFace = bulk.declare_element_side(faceInfo.elemInFirstBlock,
                                                          faceInfo.faceOrdinalInFirstBlock,
                                                          stk::mesh::PartVector{faceInfo.blockPair.first});
    STK_ThrowRequire(newFace != faceInfo.internalFace);

    bulk.copy_entity_fields(faceInfo.internalFace, newFace);
    bulk.change_entity_parts(faceInfo.internalFace, secondBlockSurfaces, firstBlockSurfaces);
    bulk.change_entity_parts(newFace, firstBlockSurfaces, secondBlockSurfaces);
  } else {
    bulk.destroy_relation(faceInfo.elemInSecondBlock, faceInfo.internalFace, faceInfo.faceOrdinalInSecondBlock);
    bulk.destroy_relation(faceInfo.elemInFirstBlock, faceInfo.internalFace, faceInfo.faceOrdinalInFirstBlock);

    stk::mesh::Entity newFaceOnSecondBlock = bulk.declare_element_side(faceInfo.elemInSecondBlock,
                                                                       faceInfo.faceOrdinalInSecondBlock,
                                                                       stk::mesh::PartVector{faceInfo.blockPair.second});
    STK_ThrowRequire(newFaceOnSecondBlock != faceInfo.internalFace);
    bulk.copy_entity_fields(faceInfo.internalFace, newFaceOnSecondBlock);
    bulk.change_entity_parts(newFaceOnSecondBlock, secondBlockSurfaces, firstBlockSurfaces);

    bulk.destroy_entity(faceInfo.internalFace, false);

    stk::mesh::Entity newFaceOnFirstBlock = bulk.declare_element_side(faceInfo.elemInFirstBlock,
                                                                      faceInfo.faceOrdinalInFirstBlock,
                                                                      stk::mesh::PartVector{faceInfo.blockPair.first});
    STK_ThrowRequire(newFaceOnFirstBlock != newFaceOnSecondBlock);

    bulk.copy_entity_fields(newFaceOnSecondBlock, newFaceOnFirstBlock);
    bulk.change_entity_parts(newFaceOnFirstBlock, firstBlockSurfaces, secondBlockSurfaces);
  }
}

InternalFaceDisconnectState get_internal_face_disconnect_state(stk::mesh::BulkData& bulk,
                                                               stk::mesh::Entity elemInBlock,
                                                               bool isSidesetAttachedToBlock)
{
  InternalFaceDisconnectState blockState = InternalFaceDisconnectState::NOOP;

  if(isSidesetAttachedToBlock) {
    if(bulk.is_valid(elemInBlock) && bulk.bucket(elemInBlock).owned()) {
      blockState = InternalFaceDisconnectState::KEEP;
    } else {
      blockState = InternalFaceDisconnectState::DESTROY;
    }
  }

  return blockState;
}

std::vector<const stk::mesh::Part*>
get_surface_parts_for_internal_face_block(stk::mesh::BulkData& bulk, stk::mesh::Part* block, stk::mesh::Entity face)
{
  auto& meta = bulk.mesh_meta_data();
  auto surfaces = meta.get_surfaces_touched_by_block(block);

  const stk::mesh::Bucket& faceBucket = bulk.bucket(face);

  auto is_not_selected = [&](const stk::mesh::Part* surface) {
    return !faceBucket.member(*surface) ;
  };
  auto end = std::remove_if(surfaces.begin(), surfaces.end(), is_not_selected);
  surfaces.erase(end, surfaces.end());

  return surfaces;
}

void disconnect_internal_face(stk::mesh::BulkData& bulk, const InternalFaceInfo& faceInfo)
{
  auto firstBlockSurfaces = get_surface_parts_for_internal_face_block(bulk, faceInfo.blockPair.first, faceInfo.internalFace);
  firstBlockSurfaces.push_back(faceInfo.blockPair.first);

  auto secondBlockSurfaces = get_surface_parts_for_internal_face_block(bulk, faceInfo.blockPair.second, faceInfo.internalFace);
  secondBlockSurfaces.push_back(faceInfo.blockPair.second);

  InternalFaceDisconnectState firstBlockState = get_internal_face_disconnect_state(bulk,
                                                                                   faceInfo.elemInFirstBlock,
                                                                                   faceInfo.isInBlockPairSideset.first);
  InternalFaceDisconnectState secondBlockState = get_internal_face_disconnect_state(bulk,
                                                                                    faceInfo.elemInSecondBlock,
                                                                                    faceInfo.isInBlockPairSideset.second);

  bool wedgedSidesetAttachedToBothBlocksWithOneLocalElement =
      ( firstBlockState == InternalFaceDisconnectState::KEEP && secondBlockState == InternalFaceDisconnectState::DESTROY) ||
      (secondBlockState == InternalFaceDisconnectState::KEEP &&  firstBlockState == InternalFaceDisconnectState::DESTROY);

  if(wedgedSidesetAttachedToBothBlocksWithOneLocalElement) {
      detach_internal_face(bulk, faceInfo, firstBlockState, secondBlockState, firstBlockSurfaces, secondBlockSurfaces);
  }

  bool wedgedSidesetAttachedToOneBlockWithLocalAndRemoteElement =
      ( firstBlockState == InternalFaceDisconnectState::NOOP && secondBlockState == InternalFaceDisconnectState::DESTROY) ||
      (secondBlockState == InternalFaceDisconnectState::NOOP &&  firstBlockState == InternalFaceDisconnectState::DESTROY);

  if(wedgedSidesetAttachedToOneBlockWithLocalAndRemoteElement) {
    destroy_internal_face(bulk, faceInfo);
  }

  bool wedgedSidesetAttachedToOneBlockWithLocalElementInSidesetBlock =
      ( firstBlockState == InternalFaceDisconnectState::NOOP && secondBlockState == InternalFaceDisconnectState::KEEP) ||
      (secondBlockState == InternalFaceDisconnectState::NOOP &&  firstBlockState == InternalFaceDisconnectState::KEEP);

  if(wedgedSidesetAttachedToOneBlockWithLocalElementInSidesetBlock) {
    detach_internal_face(bulk, faceInfo, firstBlockState, secondBlockState, firstBlockSurfaces, secondBlockSurfaces);
  }

  bool wedgedSidesetAttachedToBothBlocksWithBothLocalElements =
      (firstBlockState == InternalFaceDisconnectState::KEEP && secondBlockState == InternalFaceDisconnectState::KEEP);

  if(wedgedSidesetAttachedToBothBlocksWithBothLocalElements) {
    detach_and_duplicate_internal_face(bulk, faceInfo, firstBlockState, secondBlockState, firstBlockSurfaces, secondBlockSurfaces);
  }
}

void disconnect_faces(stk::mesh::BulkData& bulk, LinkInfo& info)
{
  for (auto faceInfo : info.internalSides) {
    update_internal_face_entity_relation(bulk, faceInfo);
  }

  if (bulk.has_face_adjacent_element_graph()) {
    bulk.get_face_adjacent_element_graph().fill_from_mesh();
  }

  for (auto faceInfo : info.internalSides) {
    disconnect_internal_face(bulk, faceInfo);
  }
}

std::vector<int> get_node_sharing_for_restoration(stk::mesh::BulkData& bulk, const stk::mesh::Part* blockPart,
                                                  const DisconnectGroup& group, stk::mesh::Entity destNode, LinkInfo& info)
{
  std::vector<int> sharingProcs;

  STK_ThrowRequire(bulk.is_valid(destNode));
  std::vector<int> procs = group.get_sharing_procs(*blockPart);

  for(int proc : procs) {
    if(proc == bulk.parallel_rank()) { continue; }
    sharingProcs.push_back(proc);
  }

  return sharingProcs;
}

void restore_node_sharing(stk::mesh::BulkData& bulk, const NodeMapValue& nodeMapValue,
                          stk::mesh::Entity destNode, LinkInfo& info)
{
  ReconnectMapKey mapKey = std::make_pair(nodeMapValue.oldNodeId, nodeMapValue.reconnectGroupId);
  std::vector<int> commonProcs = get_reconnect_procs(info, mapKey);

  for(int proc : commonProcs) {
    bulk.add_node_sharing(destNode, proc);
  }
}

void update_reconnected_entity_relation(stk::mesh::BulkData& bulk, stk::mesh::Entity newNode, stk::mesh::Entity reconnectNode, stk::mesh::Entity entity)
{
  bulk.copy_entity_fields(newNode, reconnectNode);

  unsigned numNodes = bulk.num_connectivity(entity, stk::topology::NODE_RANK);
  const stk::mesh::Entity * elemNodes = bulk.begin(entity, stk::topology::NODE_RANK);
  stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(entity, stk::topology::NODE_RANK);
  for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
    if (elemNodes[iNode] == newNode) {
      stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrdinals[iNode];
      bulk.destroy_relation(entity, newNode, nodeOrd);
      bulk.declare_relation(entity, reconnectNode, nodeOrd);
    }
  }
}

void reconnect_sub_rank(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::Entity newNode,
                        const stk::mesh::Entity reconnectNode, stk::mesh::Entity elem)
{
  const stk::mesh::Entity* subEntities = bulk.begin(elem, rank);
  unsigned numSubEntities = bulk.num_connectivity(elem, rank);
  for(unsigned i = 0; i < numSubEntities; i++) {
    update_reconnected_entity_relation(bulk, newNode, reconnectNode, subEntities[i]);
  }
}

void reconnect_sub_ranks(stk::mesh::BulkData& bulk, stk::mesh::Entity newNode,
                         const stk::mesh::Entity reconnectNode, stk::mesh::Entity elem)
{
  reconnect_sub_rank(bulk, bulk.mesh_meta_data().side_rank(), newNode, reconnectNode, elem);

  if(bulk.mesh_meta_data().side_rank() != stk::topology::EDGE_RANK) {
    reconnect_sub_rank(bulk, stk::topology::EDGE_RANK, newNode, reconnectNode, elem);
  }
}

void reconnect_elements(stk::mesh::BulkData& bulk, const BlockPair & blockPair, const NodeMapKey& key, const NodeMapValue& value,
                        LinkInfo& info)
{
  auto& disconnectedGroup = key.disconnectedGroup;

  ReconnectMapKey mapKey = std::make_pair(value.oldNodeId, value.reconnectGroupId);
  stk::mesh::EntityId reconnectId = get_reconnect_node_id(info, mapKey);
  stk::mesh::Entity reconnectNode = bulk.get_entity(stk::topology::NODE_RANK, reconnectId);

  if(!bulk.is_valid(reconnectNode)) {
    reconnectNode = bulk.declare_node(reconnectId);
  }

  stk::mesh::EntityId newNodeId = value.newNodeId;
  stk::mesh::Entity newNode = bulk.get_entity(stk::topology::NODE_RANK, newNodeId);
  STK_ThrowRequire(bulk.is_valid(newNode));

  restore_node_sharing(bulk, value, reconnectNode, info);

  for (stk::mesh::Entity elem : disconnectedGroup.get_entities()) {
    update_reconnected_entity_relation(bulk, newNode, reconnectNode, elem);
    reconnect_sub_ranks(bulk, newNode, reconnectNode, elem);
  }
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
  STK_ThrowRequire(iter != info.end());
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


bool can_be_reconnected(const DisconnectGroup& disconnectedGroup, const NodeMapValue& nodeMapValue, const BlockPair& blockPair, LinkInfo& info)
{
  return can_be_reconnected(disconnectedGroup, nodeMapValue, blockPair, nodeMapValue.boundaryNode, info);
}

bool can_be_reconnected(const DisconnectGroup& disconnectedGroup, const NodeMapValue& nodeMapValue, const BlockPair& blockPair, stk::mesh::Entity currentEntity, LinkInfo& info)
{
  const stk::mesh::BulkData& bulk = disconnectedGroup.get_bulk();

  if(!bulk.is_valid(currentEntity)) { return false; }

  const stk::mesh::PartVector& blockMembership = nodeMapValue.oldBlockMembership;
  bool isOriginalMember = false;
  bool isInOneOfBlockPair = is_in_preserved_block(bulk,currentEntity,blockPair) != is_in_disconnected_block(bulk,currentEntity,blockPair);

  if (disconnectedGroup.has_block_pair() && isInOneOfBlockPair) {
    isOriginalMember = std::binary_search(blockMembership.begin(), blockMembership.end(), blockPair.first, stk::mesh::PartLess()) &&
        std::binary_search(blockMembership.begin(), blockMembership.end(), blockPair.second, stk::mesh::PartLess());
  }

  return isOriginalMember;
}

void sanitize_node_map(NodeMapType& nodeMap, LinkInfo& info)
{
  for (NodeMapType::iterator nodeMapEntryIt = nodeMap.begin(); nodeMapEntryIt != nodeMap.end();) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    if(!disconnectedGroup.get_bulk().is_valid(disconnectedGroup.get_node())) {
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
      procBuff.pack<stk::mesh::EntityId>(reconnectMapEntry.first.first);
      procBuff.pack<unsigned>(reconnectMapEntry.first.second);
      procBuff.pack<stk::mesh::EntityId>(reconnectMapEntry.second.reconnectNodeId);
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
      unsigned groupId;

      procBuff.unpack<stk::mesh::EntityId>(referenceNodeId);
      procBuff.unpack<unsigned>(groupId);
      procBuff.unpack<stk::mesh::EntityId>(reconnectNodeId);

      ReconnectMapKey mapKey = std::make_pair(referenceNodeId, groupId);
      auto reconnectMapIter = info.reconnectMap.find(mapKey);

      STK_ThrowRequire(reconnectMapIter != info.reconnectMap.end());

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

  info.flush();
  unpack_reconnect_node_information(bulk, commSparse, info);
}

template <typename Functor>
void traverse_transitive_relations(stk::mesh::BulkData& bulk, const BlockPairVector& blockPairsToReconnect,
                                   NodeMapType::iterator nodeMapEntryIt, stk::mesh::Entity currentEntity,
                                   LinkInfo& info, Functor func)
{
  if(blockPairsToReconnect.empty()) { return; }
  if(!bulk.is_valid(currentEntity)) { return; }

  BlockPairVector blockPairs;
  stk::mesh::PartLess partLess;
  const stk::mesh::PartVector& oldBlockMembership = nodeMapEntryIt->second.oldBlockMembership;
  std::vector<stk::mesh::Part*> transitiveTriplet(3);

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

          bool isTransitive = func(transitiveTriplet, nodeMapEntryIt, currentEntity);

          if(isTransitive) {
            ReconnectMapKey mapKey = std::make_pair(nodeMapEntryIt->second.oldNodeId, nodeMapEntryIt->second.reconnectGroupId);
            auto reconnectMapIter = info.reconnectMap.find(mapKey);
            STK_ThrowAssert(reconnectMapIter != info.reconnectMap.end());
            insert_parts_uniquely(transitiveTriplet, reconnectMapIter->second.reconnectParts);
          }
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

void initialize_reconnect_node_id_info_for_node_entry(stk::mesh::BulkData& bulk, NodeMapType::iterator nodeMapEntryIt,
                                                      const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
{
  stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
  stk::mesh::Entity referenceNode = bulk.get_entity(stk::topology::NODE_RANK, referenceNodeId);
  stk::mesh::Entity currentEntity = nodeMapEntryIt->second.boundaryNode;

  stk::mesh::EntityId currentNodeId = bulk.identifier(currentEntity);
  ReconnectGroup reconnectGroup = get_group_for_reconnect_node(nodeMapEntryIt->first.disconnectedGroup, nodeMapEntryIt->second, blockPairsToReconnect, info);
  unsigned groupId = reconnectGroup.get_id();
  nodeMapEntryIt->second.reconnectGroupId = groupId;

  ReconnectMapKey mapKey = std::make_pair(referenceNodeId, groupId);
  auto reconnectMapIter = info.reconnectMap.find(mapKey);

  if(bulk.is_valid(referenceNode) && bulk.bucket(referenceNode).member_any(reconnectGroup.get_blocks())) {
    currentNodeId = referenceNodeId;
  }

  ReconnectNodeInfo* reconnectInfo = nullptr;

  if(reconnectMapIter != info.reconnectMap.end()) {
    reconnectInfo = &reconnectMapIter->second;
    reconnectMapIter->second.reconnectNodeId = std::min(reconnectMapIter->second.reconnectNodeId, currentNodeId);
  } else {
    ReconnectNodeInfo newReconnectInfo;
    newReconnectInfo.reconnectNodeId = currentNodeId;
    auto value = info.reconnectMap.insert(std::make_pair(mapKey, newReconnectInfo));
    reconnectInfo = &value.first->second;
  }
  insert_parts_uniquely(reconnectGroup.get_blocks(), reconnectInfo->reconnectParts);
}

void initialize_reconnect_node_id_info(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect, LinkInfo& info)
{
  for(auto nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    initialize_reconnect_node_id_info_for_node_entry(bulk, nodeMapEntryIt, blockPairsToReconnect, info);
  }

  for(auto nodeMapEntryIt = info.originalNodeMap.begin(); nodeMapEntryIt != info.originalNodeMap.end(); ++nodeMapEntryIt) {
    initialize_reconnect_node_id_info_for_node_entry(bulk, nodeMapEntryIt, blockPairsToReconnect, info);
  }
}

void fill_direct_reconnect_node_info(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                     NodeMapType::iterator nodeMapEntryIt, stk::mesh::Entity currentEntity,
                                     LinkInfo& info)
{
  const DisconnectGroup& disconnectGroup = nodeMapEntryIt->first.disconnectedGroup;

  for(const BlockPair& blockPair : blockPairsToReconnect) {
    const stk::mesh::Part & srcBlock = *blockPair.first;
    const stk::mesh::Part & blockToReconnect = *blockPair.second;

    if(can_be_reconnected(disconnectGroup, nodeMapEntryIt->second, blockPair, info)) {
      stk::mesh::EntityId referenceNodeId = nodeMapEntryIt->second.oldNodeId;
      ReconnectNodeInfo* reconnectInfo = nullptr;
      unsigned groupId = nodeMapEntryIt->second.reconnectGroupId;
      STK_ThrowRequire(groupId != std::numeric_limits<unsigned>::max());
      ReconnectMapKey mapKey = std::make_pair(referenceNodeId, groupId);

      auto reconnectMapIter = info.reconnectMap.find(mapKey);
      STK_ThrowRequire(reconnectMapIter != info.reconnectMap.end());
      reconnectInfo = &reconnectMapIter->second;
      reconnectInfo->reconnectNodeId = std::min(bulk.identifier(currentEntity), reconnectInfo->reconnectNodeId);
      stk::util::insert_keep_sorted_and_unique(currentEntity, reconnectInfo->relatedNodes);

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

  traverse_transitive_relations(bulk, blockPairsToReconnect, it, currentEntity, info, transitiveFunc);
  traverse_transitive_relations(bulk, blockPairsToReconnect, it, it->first.disconnectedGroup.get_node(), info, transitiveFunc);
}

unsigned get_group_id_for_reconnect_parts(const stk::mesh::PartVector& reconnectParts)
{
  unsigned groupId = std::numeric_limits<unsigned>::max();

  for(stk::mesh::Part* part : reconnectParts) {
    groupId = std::min(groupId, part->mesh_meta_data_ordinal());
  }
  return groupId;
}

typedef std::map<ReconnectMapKey,unsigned> MergeGroupsMap;

void update_group_id_in_node_map(NodeMapType& nodeMap, const MergeGroupsMap& mergeGroupsMap, LinkInfo& info)
{
  for(auto it = nodeMap.begin(); it != nodeMap.end(); ++it) {
    unsigned groupId = it->second.reconnectGroupId;
    ReconnectMapKey mapKey = std::make_pair(it->second.oldNodeId, groupId);
    auto mergeGroupIt = mergeGroupsMap.find(mapKey);
    if(mergeGroupIt != mergeGroupsMap.end()) {
      it->second.reconnectGroupId = mergeGroupIt->second;
    }
  }
}

void merge_reconnect_groups(const stk::mesh::BulkData& bulk, LinkInfo& info)
{
  MergeGroupsMap mergeGroupsMap;

  for(auto iter = info.reconnectMap.begin(); iter != info.reconnectMap.end(); ) {
    unsigned groupId = iter->first.second;
    unsigned mergedGroupId = get_group_id_for_reconnect_parts(iter->second.reconnectParts);
    if(groupId != mergedGroupId) {
      ReconnectMapKey mapKey = std::make_pair(iter->first.first, mergedGroupId);
      auto srcGroupIter = info.reconnectMap.find(mapKey);
      if(srcGroupIter == info.reconnectMap.end()) {
        info.reconnectMap.insert(std::make_pair(mapKey, iter->second));
        iter = info.reconnectMap.begin();
        continue;
      }
      STK_ThrowRequire(srcGroupIter->second.reconnectParts == iter->second.reconnectParts);

      mergeGroupsMap.insert(std::make_pair(iter->first, mergedGroupId));
      iter = info.reconnectMap.erase(iter);
    }
    else {
      ++iter;
    }
  }

  update_group_id_in_node_map(info.clonedNodeMap, mergeGroupsMap, info);
  update_group_id_in_node_map(info.originalNodeMap, mergeGroupsMap, info);
}

void determine_local_reconnect_node_id(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                       LinkInfo& info)
{
  auto fill_reconnect_node_info_func =
      [&](const stk::mesh::PartVector& transitiveBlockList_, NodeMapType::iterator it, stk::mesh::Entity currentEntity_)
  {
    unsigned groupId = it->second.reconnectGroupId;

    ReconnectMapKey mapKey = std::make_pair(it->second.oldNodeId, groupId);
    auto reconnectMapIter = info.reconnectMap.find(mapKey);
    if(reconnectMapIter == info.reconnectMap.end()) { return false; }
    ReconnectNodeInfo* reconnectInfo = &reconnectMapIter->second;

    reconnectInfo->reconnectNodeId = std::min(bulk.identifier(currentEntity_), reconnectInfo->reconnectNodeId);
    return true;
  };

  initialize_reconnect_node_id_info(bulk, blockPairsToReconnect, info);

  for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, it->second.boundaryNode, info, fill_reconnect_node_info_func);
  }
  for(auto it = info.originalNodeMap.begin(); it != info.originalNodeMap.end(); ++it) {
    const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
    stk::mesh::Entity currNode = disconnectGroup.get_node();
    fill_reconnect_node_info(bulk, blockPairsToReconnect, it, currNode, info, fill_reconnect_node_info_func);
  }
  merge_reconnect_groups(bulk, info);
}

void update_reconnect_node_sharing(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                                   LinkInfo& info)
{
  auto update_sharing_info_func =
      [&](const stk::mesh::PartVector& transitiveBlockList, NodeMapType::iterator it, stk::mesh::Entity currentEntity_)
  {
    const DisconnectGroup& disconnectedGroup = it->first.disconnectedGroup;
    stk::mesh::Entity boundaryNode = it->second.boundaryNode;
    unsigned groupId = it->second.reconnectGroupId;

    ReconnectMapKey mapKey = std::make_pair(it->second.oldNodeId, groupId);
    auto reconnectMapIter = info.reconnectMap.find(mapKey);
    if(reconnectMapIter == info.reconnectMap.end()) { return false; }
    ReconnectNodeInfo* reconnectInfo = &reconnectMapIter->second;
    std::vector<int> sharingProcs;

    for(const stk::mesh::Part* blockPart : transitiveBlockList) {
      sharingProcs = get_node_sharing_for_restoration(bulk, blockPart, disconnectedGroup, boundaryNode, info);

      for(int proc : sharingProcs) {
        stk::util::insert_keep_sorted_and_unique(proc, reconnectInfo->reconnectProcs);
      }
    }

    return false;
  };

  for(auto it = info.clonedNodeMap.begin(); it != info.clonedNodeMap.end(); ++it) {
    update_procs_for_transitive_reconnection(bulk, blockPairsToReconnect, it, it->second.boundaryNode, info, update_sharing_info_func);
  }

  for(auto it = info.originalNodeMap.begin(); it != info.originalNodeMap.end(); ++it) {
    const DisconnectGroup& disconnectGroup = it->first.disconnectedGroup;
    stk::mesh::Entity currNode = disconnectGroup.get_node();
    update_procs_for_transitive_reconnection(bulk, blockPairsToReconnect, it, currNode, info, update_sharing_info_func);
  }
}

void reconnect_block_pair(stk::mesh::BulkData& bulk, const BlockPair & blockPair, impl::LinkInfo& info)
{
  const stk::mesh::Part & srcBlock = *blockPair.first;
  const stk::mesh::Part & blockToReconnect = *blockPair.second;

  for (NodeMapType::iterator nodeMapEntryIt = info.clonedNodeMap.begin(); nodeMapEntryIt != info.clonedNodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;

    bool shouldReconnect = can_be_reconnected(disconnectedGroup, nodeMapEntryIt->second, blockPair, info);
    if(shouldReconnect && bulk.is_valid(nodeMapEntryIt->second.boundaryNode)) {
      impl::reconnect_elements(bulk, blockPair, nodeMapEntryIt->first, nodeMapEntryIt->second, info);
    }
  }

  for (NodeMapType::iterator nodeMapEntryIt = info.originalNodeMap.begin(); nodeMapEntryIt != info.originalNodeMap.end(); ++nodeMapEntryIt) {
    const impl::DisconnectGroup& disconnectedGroup = nodeMapEntryIt->first.disconnectedGroup;
    stk::mesh::Entity originalNode = disconnectedGroup.get_node();

    STK_ThrowRequire(bulk.is_valid(originalNode));
    bool isInEitherBlockPair = bulk.bucket(originalNode).member(srcBlock) || bulk.bucket(originalNode).member(blockToReconnect);
    bool shouldReconnect = can_be_reconnected(disconnectedGroup, nodeMapEntryIt->second, blockPair, info);
    if(isInEitherBlockPair && shouldReconnect) {
      impl::restore_node_sharing(bulk, nodeMapEntryIt->second, disconnectedGroup.get_node(), info);
    }
  }
}

void update_node_id(stk::mesh::EntityId newNodeId, int proc, LinkInfo& info, const DisconnectGroup& group) {
  STK_ThrowRequire(!group.get_parts().empty());
  NodeMapKey nodeMapKey(group.get_node(), group);
  const auto & nodeMapIt = info.clonedNodeMap.find(nodeMapKey);

  if (nodeMapIt != info.clonedNodeMap.end()) {
    NodeMapValue & newNodeData = nodeMapIt->second;
    newNodeData.newNodeId = std::min(newNodeData.newNodeId, newNodeId);
    newNodeData.oldNodeId = group.get_node_id();
    newNodeData.sharingProcs.push_back(proc);
  }
}

void clean_up_aura(stk::mesh::BulkData& bulk, LinkInfo& info)
{
  stk::mesh::EntityVector allNodes;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, bulk.mesh_meta_data().locally_owned_part(), allNodes);

  for(stk::mesh::Entity node : allNodes) {
    const int numElems = bulk.num_connectivity(node, stk::topology::ELEMENT_RANK);
    const stk::mesh::Entity* elems = bulk.begin(node, stk::topology::ELEMENT_RANK);
    int numDestroyedElems = 0;

    for(int i = numElems-1; i >= 0; --i) {
      stk::mesh::Entity elem = elems[i];
      if(!bulk.bucket(elem).owned()) {
        bulk.destroy_entity(elem);
        ++numDestroyedElems;
      }
    }

    if(numDestroyedElems == numElems) {
      bulk.destroy_entity(node);
    }
  }
}

void disconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToDisconnect,
                            LinkInfo& info)
{
  bulk.modification_begin();

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    add_nodes_to_disconnect(bulk, blockPairsToDisconnect[i], info);
    add_faces_to_disconnect(bulk, blockPairsToDisconnect[i], info);
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

  disconnect_faces(bulk, info);

  clean_up_aura(bulk, info);

  bulk.modification_end();

  if (bulk.has_face_adjacent_element_graph()) {
    bulk.delete_face_adjacent_element_graph();
    bulk.initialize_face_adjacent_element_graph();
  }
}

stk::mesh::EntityVector extract_nodes(const stk::mesh::BulkData& bulk, LinkInfo& info)
{
  stk::mesh::EntityVector nodes;

  for(auto mapPair : info.clonedNodeMap) {
    if(!bulk.is_valid(mapPair.second.boundaryNode)) { continue; }
    stk::util::insert_keep_sorted_and_unique(mapPair.second.boundaryNode, nodes);
  }

  for(auto mapPair : info.originalNodeMap) {
    if(!bulk.is_valid(mapPair.second.boundaryNode)) { continue; }
    stk::util::insert_keep_sorted_and_unique(mapPair.second.boundaryNode, nodes);
  }

  return nodes;
}

stk::mesh::PartVector get_block_parts(const BlockPairVector& blocksToDisconnect)
{
  stk::mesh::PartVector blockParts;

  for(const BlockPair& blockPair : blocksToDisconnect) {
    stk::util::insert_keep_sorted_and_unique(blockPair.first, blockParts);
    stk::util::insert_keep_sorted_and_unique(blockPair.second, blockParts);
  }

  return blockParts;
}

stk::mesh::EntityVector get_affected_nodes(const stk::mesh::BulkData& bulk, const BlockPairVector& blocksToDisconnect)
{
  stk::mesh::EntityVector nodes;
  stk::mesh::PartVector blockParts = get_block_parts(blocksToDisconnect);
  stk::mesh::Selector partSelector = stk::mesh::selectUnion(blockParts);

  stk::mesh::get_selected_entities(partSelector, bulk.buckets(stk::topology::NODE_RANK), nodes);

  return nodes;
}

void reconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<BlockPair>& blockPairsToReconnect,
                           LinkInfo& info)
{
  bulk.modification_begin();

  sanitize_node_map(info.originalNodeMap, info);

  determine_local_reconnect_node_id(bulk, blockPairsToReconnect, info);

  update_reconnect_node_sharing(bulk, blockPairsToReconnect, info);

  communicate_reconnect_node_information(bulk, info);

  for(const BlockPair & blockPair : blockPairsToReconnect) {
    reconnect_block_pair(bulk, blockPair, info);
  }

  clean_up_aura(bulk, info);

  bulk.modification_end();
}
}
}
}
