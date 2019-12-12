#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include <map>
#include <algorithm>

namespace stk {
namespace tools {
namespace impl {

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  const bool isElementPart = (part.primary_entity_rank() == stk::topology::ELEM_RANK);
  const bool isIoPart      = stk::io::has_io_part_attribute(part);
  return (isElementPart && isIoPart);
}

unsigned get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element)
{
  const stk::mesh::PartVector & elementParts = bulk.bucket(element).supersets();
  for (stk::mesh::Part * part : elementParts) {
    if (is_block(bulk, *part)) {
      return part->mesh_meta_data_ordinal();
    }
  }
  return -1;
}

void add_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                             const BlockPairType & blockPair,
                             impl::NodeMapType & nodeMap)
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
      if (needToCloneNode && (bulk.state(node) != stk::mesh::Created)) {
        nodeMap[NodeMapKey(node, secondBlock)] = NodeMapValue();
      }
    }
  }
}

void create_new_duplicate_node_IDs(stk::mesh::BulkData & bulk, NodeMapType & nodeMap)
{
  std::vector<stk::mesh::EntityId> newNodeIDs;
  bulk.generate_new_ids(stk::topology::NODE_RANK, nodeMap.size(), newNodeIDs);

  size_t newNodeIdx = 0;
  for (auto & nodeMapEntry : nodeMap) {
    nodeMapEntry.second.newNodeId = newNodeIDs[newNodeIdx++];
  }
}

void communicate_shared_node_information(stk::mesh::BulkData & bulk, NodeMapType & nodeMap)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::CommSparse commSparse(bulk.parallel());

  for(int phase = 0; phase < 2; ++phase) {

    for (const auto & nodeMapEntry : nodeMap) {
      const stk::mesh::Entity node = nodeMapEntry.first.parentNode;
      if (!bulk.bucket(node).shared()) continue;

      const stk::mesh::Part & blockToDisconnect = nodeMapEntry.first.disconnectedBlock;
      const stk::mesh::EntityId newNodeId = nodeMapEntry.second.newNodeId;

      const stk::mesh::Entity * elems = bulk.begin_elements(node);
      const unsigned numElems = bulk.num_elements(node);
      bool needToCommunicate = false;
      for (unsigned i = 0; i < numElems; ++i) {
        const stk::mesh::Bucket & elemBucket = bulk.bucket(elems[i]);
        if (elemBucket.member(blockToDisconnect) && elemBucket.owned()) {
          needToCommunicate = true;
          break;
        }
      }

      if (needToCommunicate) {
        std::vector<int> sharingProcs;
        bulk.comm_shared_procs(bulk.entity_key(node), sharingProcs);
        for (const int proc : sharingProcs) {
          stk::CommBuffer& procBuff = commSparse.send_buffer(proc);
          procBuff.pack<stk::mesh::EntityId>(bulk.identifier(node));
          procBuff.pack<stk::mesh::PartOrdinal>(blockToDisconnect.mesh_meta_data_ordinal());
          procBuff.pack<stk::mesh::EntityId>(newNodeId);
        }
      }
    }

    if (phase == 0) {
      commSparse.allocate_buffers();
    }
    else {
      commSparse.communicate();
    }
  }


  for (int proc = 0; proc < bulk.parallel_size(); ++proc) {
    while (commSparse.recv_buffer(proc).remaining()) {
      stk::CommBuffer & procBuff = commSparse.recv_buffer(proc);
      stk::mesh::EntityId parentNodeId;
      stk::mesh::PartOrdinal blockToDisconnectOrdinal;
      stk::mesh::EntityId newNodeId;
      procBuff.unpack<stk::mesh::EntityId>(parentNodeId);
      procBuff.unpack<stk::mesh::PartOrdinal>(blockToDisconnectOrdinal);
      procBuff.unpack<stk::mesh::EntityId>(newNodeId);

      NodeMapKey nodeMapKey(bulk.get_entity(stk::topology::NODE_RANK, parentNodeId), meta.get_part(blockToDisconnectOrdinal));
      const auto & nodeMapIt = nodeMap.find(nodeMapKey);
      if (nodeMapIt != nodeMap.end()) {
        NodeMapValue & newNodeData = nodeMapIt->second;
        newNodeData.newNodeId = std::min(newNodeData.newNodeId, newNodeId);
        newNodeData.sharingProcs.push_back(proc);
      }
    }
  }
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

  return std::move(blockPairsToDisconnect);
}

void disconnect_elements(stk::mesh::BulkData & bulk,
                         const BlockPairType & blockPair,
                         NodeMapType & nodeMap)
{
  const stk::mesh::Part & blockToDisconnect = *blockPair.second;

  for (auto nodeMapEntryIt = nodeMap.begin(); nodeMapEntryIt != nodeMap.end(); ) {
    if (nodeMapEntryIt->first.disconnectedBlock == blockToDisconnect) {
      const stk::mesh::Entity node = nodeMapEntryIt->first.parentNode;
      const std::vector<stk::mesh::Entity> connectedElems(bulk.begin_elements(node), bulk.end_elements(node));

      for (stk::mesh::Entity elem : connectedElems) {
        if (bulk.bucket(elem).member(blockToDisconnect)) {
          stk::mesh::EntityId newNodeId = nodeMapEntryIt->second.newNodeId;
          stk::mesh::Entity newNode = bulk.declare_node(newNodeId);
          bulk.copy_entity_fields(node, newNode);

          for (int sharingProc : nodeMapEntryIt->second.sharingProcs) {
            bulk.add_node_sharing(newNode, sharingProc);
          }

          unsigned numNodes = bulk.num_connectivity(elem, stk::topology::NODE_RANK);
          const stk::mesh::Entity * elemNodes = bulk.begin(elem, stk::topology::NODE_RANK);
          stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(elem, stk::topology::NODE_RANK);
          for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
            if (elemNodes[iNode] == node) {
              bulk.destroy_relation(elem, node, nodeOrdinals[iNode]);
              bulk.declare_relation(elem, newNode, nodeOrdinals[iNode]);
            }
          }

          if (bulk.num_elements(node) == 0) {
            bulk.destroy_entity(node);
          }
        }
      }
      nodeMapEntryIt = nodeMap.erase(nodeMapEntryIt);
    }
    else {
      ++nodeMapEntryIt;
    }
  }
}

}
}
}
