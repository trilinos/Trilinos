#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_mesh/base/SkinMeshUtil.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_io/IossBridge.hpp"
#include <map>

namespace stk {
namespace tools {
namespace impl {

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  const bool isElementPart = (part.primary_entity_rank() == stk::topology::ELEM_RANK);
  const bool isIoPart      = stk::io::has_io_part_attribute(part);
  const bool hasPartId     = (part.id() > 0);
  return (isElementPart && isIoPart && hasPartId);
}

int64_t get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element)
{
  const stk::mesh::PartVector & elementParts = bulk.bucket(element).supersets();
  for (stk::mesh::Part * part : elementParts) {
    if (is_block(bulk, *part)) {
      return part->id();
    }
  }
  return -1;
}

void get_nodes_for_element_side(const stk::mesh::BulkData & bulk,
                                stk::mesh::Entity element,
                                stk::mesh::ConnectivityOrdinal sideOrdinal,
                                std::vector<stk::mesh::Entity> & sideNodes)
{
  const stk::mesh::Entity * elemNodes = bulk.begin_nodes(element);
  const stk::topology elemTopology = bulk.bucket(element).topology();
  const stk::topology sideTopology = elemTopology.side_topology(sideOrdinal);
  sideNodes.resize(sideTopology.num_nodes());
  elemTopology.side_nodes(elemNodes, sideOrdinal, sideNodes.data());
}

void get_node_ordinals_for_element_side(const stk::mesh::BulkData & bulk,
                                        stk::mesh::Entity element,
                                        stk::mesh::ConnectivityOrdinal sideOrdinal,
                                        std::vector<stk::mesh::ConnectivityOrdinal> & sideNodeOrdinals)
{
  const stk::topology elemTopology = bulk.bucket(element).topology();
  const stk::topology sideTopology = elemTopology.side_topology(sideOrdinal);
  sideNodeOrdinals.resize(sideTopology.num_nodes());
  elemTopology.side_node_ordinals(sideOrdinal, sideNodeOrdinals.data());
}

void disconnect_side_nodes(stk::mesh::BulkData & bulk,
                           const NodeMapType nodeMap,
                           stk::mesh::Entity elem,
                           const std::vector<stk::mesh::Entity> & sideNodes,
                           const std::vector<stk::mesh::ConnectivityOrdinal> & sideNodeOrdinals)
{
  for (size_t i = 0; i < sideNodes.size(); ++i) {
    const stk::mesh::Entity node = sideNodes[i];
    const stk::mesh::ConnectivityOrdinal nodeOrdinal = sideNodeOrdinals[i];
    const auto nodeIt = nodeMap.find(node);
    if (nodeIt != nodeMap.end()) {
      stk::mesh::Entity newNode = nodeIt->second;
//      std::cout << "  Swapping " << bulk.entity_key(node) << " to " << bulk.entity_key(newNode) << std::endl;
      bulk.destroy_relation(elem, node, nodeOrdinal);
      bulk.declare_relation(elem, newNode, nodeOrdinal);
    }
    else {
//      std::cout << "  Failed to query new node mapped from " << bulk.entity_key(node) << std::endl;
    }
  }
}

NodeMapType get_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                                    const BlockPairType & blockPair,
                                    const std::vector<stk::mesh::SideSetEntry> & sidesToDisconnect)
{
  const int64_t blockIdToDisconnect = (blockPair.first->id() > blockPair.second->id()) ? blockPair.first->id() : blockPair.second->id();
  std::vector<stk::mesh::Entity> sideNodes;
  NodeMapType nodeMap;

  for (const stk::mesh::SideSetEntry & sideSetEntry : sidesToDisconnect) {
    if (get_block_id_for_element(bulk, sideSetEntry.element) == blockIdToDisconnect) {
      get_nodes_for_element_side(bulk, sideSetEntry.element, sideSetEntry.side, sideNodes);
      for (const stk::mesh::Entity node : sideNodes) {
        if (bulk.state(node) != stk::mesh::Created) {
          nodeMap[node] = stk::mesh::Entity();
        }
      }
    }
  }

  return std::move(nodeMap);
}

void create_new_duplicate_nodes(stk::mesh::BulkData & bulk, NodeMapType & nodeMap)
{
  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<size_t> requests(meta.entity_rank_count(), 0);
  requests[stk::topology::NODE_RANK] = nodeMap.size();

  stk::mesh::EntityVector newNodes;
  bulk.generate_new_entities(requests, newNodes);

  size_t newNodeIdx = 0;
  for (std::pair<const stk::mesh::Entity, stk::mesh::Entity> & nodeMapEntry : nodeMap) {
    nodeMapEntry.second = newNodes[newNodeIdx++];
//    std::cout << "Creating node " << bulk.entity_key(nodeMapEntry.second) << " for old node " << bulk.entity_key(nodeMapEntry.first) << " (state=" << bulk.state(nodeMapEntry.first) << ")" << std::endl;

    bulk.copy_entity_fields(nodeMapEntry.first, nodeMapEntry.second);
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
  std::vector<BlockPairType> blockPairsToDisconnect;
  stk::mesh::PartVector allBlocksInMesh;
  get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  if (allBlocksInMesh.size() < 2u) return blockPairsToDisconnect;

  for (size_t firstBlockIdx = 0; firstBlockIdx < allBlocksInMesh.size()-1; ++firstBlockIdx) {
    for (size_t secondBlockIdx = firstBlockIdx+1; secondBlockIdx < allBlocksInMesh.size(); ++secondBlockIdx) {
      stk::mesh::Part & firstBlock  = *allBlocksInMesh[firstBlockIdx];
      stk::mesh::Part & secondBlock = *allBlocksInMesh[secondBlockIdx];
      stk::mesh::Selector boundaryBetweenBlocks = firstBlock & secondBlock;
      const stk::mesh::BucketVector & nodesOnBoundaryBetweenBlocks = bulk.get_buckets(stk::topology::NODE_RANK, boundaryBetweenBlocks);
      if (!nodesOnBoundaryBetweenBlocks.empty()) {
        blockPairsToDisconnect.emplace_back(&firstBlock, &secondBlock);
      }
    }
  }

  return std::move(blockPairsToDisconnect);
}

std::vector<SideSetType> get_sidesets_to_disconnect(stk::mesh::BulkData & bulk,
                                                    const std::vector<BlockPairType> & blockPairsToDisconnect)
{
  std::vector<SideSetType> sideSetsToDisconnect;
  sideSetsToDisconnect.resize(blockPairsToDisconnect.size());
  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    stk::mesh::Selector blockPair = *blockPairsToDisconnect[i].first | *blockPairsToDisconnect[i].second;
    sideSetsToDisconnect[i] = stk::mesh::SkinMeshUtil::get_interior_sideset(bulk, blockPair);
  }
  return std::move(sideSetsToDisconnect);
}

void disconnect_block_pair_boundary(stk::mesh::BulkData & bulk,
                                    const BlockPairType & blockPair,
                                    const SideSetType & sideSet,
                                    const NodeMapType & nodeMap)
{
  std::vector<stk::mesh::Entity> otherSideNodes;
  std::vector<stk::mesh::ConnectivityOrdinal> otherSideNodeOrdinals;
  const stk::mesh::ElemElemGraph & graph = bulk.get_face_adjacent_element_graph();

  for (const stk::mesh::SideSetEntry & sideSetEntry : sideSet) {
    const stk::mesh::Entity elem = sideSetEntry.element;
    const stk::mesh::ConnectivityOrdinal elemSide = sideSetEntry.side;
    const stk::mesh::GraphEdgesForElement & graphEdges = graph.get_edges_for_element(graph.get_local_element_id(elem));

    for (const stk::mesh::GraphEdge & edge : graphEdges) {
      if (edge.side1() == elemSide) {
        stk::mesh::Entity otherElem = graph.get_entity_from_local_id(edge.elem2());
        const int64_t elemBlockId      = get_block_id_for_element(bulk, elem);
        const int64_t otherElemBlockId = get_block_id_for_element(bulk, otherElem);

        if (elemBlockId < otherElemBlockId) {
//          std::cout << "Disconnecting " << bulk.entity_key(otherElem) << " from " << bulk.entity_key(elem) << std::endl;
          stk::mesh::ConnectivityOrdinal otherElemSide = edge.side2();
          get_nodes_for_element_side(bulk, otherElem, otherElemSide, otherSideNodes);
          get_node_ordinals_for_element_side(bulk, otherElem, otherElemSide, otherSideNodeOrdinals);

          disconnect_side_nodes(bulk, nodeMap, otherElem, otherSideNodes, otherSideNodeOrdinals);
          break;
        }
      }
    }
  }
}

void disconnect_face_adjacent_elements(stk::mesh::BulkData & bulk,
                                       const BlockPairType & blockPair,
                                       const NodeMapType & nodeMap)
{
  const stk::mesh::Part * firstBlock  = blockPair.first;
  const stk::mesh::Part * secondBlock = blockPair.second;
  if (firstBlock->id() > secondBlock->id()) {
    firstBlock = blockPair.second;
    secondBlock = blockPair.first;
  }

  for (const auto & nodeMapEntry : nodeMap) {
    const stk::mesh::Entity node = nodeMapEntry.first;
    const std::vector<stk::mesh::Entity> connectedElems(bulk.begin_elements(node), bulk.end_elements(node));

    for (stk::mesh::Entity elem : connectedElems) {
//      std::cout << "Considering " << bulk.entity_key(elem) << " for disconnect" << std::endl;
      if (bulk.bucket(elem).member(*secondBlock)) {
        const auto nodeIt = nodeMap.find(node);
//        std::cout << "  Considering node " << bulk.entity_key(node) << " for disconnect" << std::endl;
        if (nodeIt != nodeMap.end()) {
          stk::mesh::Entity newNode = nodeIt->second;
//          std::cout << "  Swapping " << bulk.entity_key(node) << " to " << bulk.entity_key(newNode) << " on " << bulk.entity_key(elem) << std::endl;

          unsigned numNodes = bulk.num_connectivity(elem, stk::topology::NODE_RANK);
          const stk::mesh::Entity * elemNodes = bulk.begin(elem, stk::topology::NODE_RANK);
          stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(elem, stk::topology::NODE_RANK);
          for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
            if (elemNodes[iNode] == node) {
              bulk.destroy_relation(elem, node, nodeOrdinals[iNode]);
              bulk.declare_relation(elem, newNode, nodeOrdinals[iNode]);
            }
          }
        }
      }
    }
  }
}

NodeMapType get_non_face_adjacent_nodes(stk::mesh::BulkData & bulk,
                                        const std::vector<BlockPairType> & blockPairsToDisconnect)
{
  stk::mesh::Selector selector;
  for (const BlockPairType & blockPair : blockPairsToDisconnect) {
    selector |= (*blockPair.first & *blockPair.second);
  }

  stk::mesh::EntityVector nodesToDisconnect;
  stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), nodesToDisconnect);

  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  std::vector<size_t> requests(meta.entity_rank_count(), 0);
  requests[stk::topology::NODE_RANK] = nodesToDisconnect.size();

  stk::mesh::EntityVector newNodes;
  bulk.generate_new_entities(requests, newNodes);

  NodeMapType nonFaceAdjacentNodeMap;
  for (unsigned i = 0; i < nodesToDisconnect.size(); ++i) {
    nonFaceAdjacentNodeMap[nodesToDisconnect[i]] = newNodes[i];
    bulk.copy_entity_fields(nodesToDisconnect[i], newNodes[i]);
  }

  return std::move(nonFaceAdjacentNodeMap);
}

void disconnect_non_face_adjacent_elements(stk::mesh::BulkData & bulk,
                                           const BlockPairType & blockPair,
                                           const NodeMapType & nodesToDisconnect)
{
  const stk::mesh::Part * firstBlock  = blockPair.first;
  const stk::mesh::Part * secondBlock = blockPair.second;
  if (firstBlock->id() > secondBlock->id()) {
    firstBlock = blockPair.second;
    secondBlock = blockPair.first;
  }

  stk::mesh::EntityVector firstBlockElems;
  stk::mesh::EntityVector secondBlockElems;

  for (auto nodePair : nodesToDisconnect) {
    firstBlockElems.clear();
    secondBlockElems.clear();
    const stk::mesh::Entity nodeToDisconnect = nodePair.first;
    const stk::mesh::Entity newNode = nodePair.second;
    const stk::mesh::Entity * elems = bulk.begin_elements(nodeToDisconnect);
    const unsigned numElems = bulk.num_elements(nodeToDisconnect);
    for (unsigned i = 0; i < numElems; ++i) {
      if (bulk.bucket(elems[i]).member(*firstBlock)) {
        firstBlockElems.push_back(elems[i]);
      }
      else if (bulk.bucket(elems[i]).member(*secondBlock)) {
        secondBlockElems.push_back(elems[i]);
      }
    }

//    std::cout << "Considering disconnection of " << bulk.entity_key(nodeToDisconnect) << " from:" << std::endl;
//    std::cout << "    " << firstBlock->name() << std::endl;
//    for (stk::mesh::Entity elem : firstBlockElems) {
//       std::cout << "        " << bulk.entity_key(elem) << std::endl;
//    }
//    std::cout << "    " << secondBlock->name() << std::endl;
//    for (stk::mesh::Entity elem : secondBlockElems) {
//       std::cout << "        " << bulk.entity_key(elem) << std::endl;
//    }

    if (firstBlockElems.empty() || secondBlockElems.empty()) continue;

//    for (stk::mesh::Entity firstElem : firstBlockElems) {
//      for (stk::mesh::Entity secondElem : secondBlockElems) {
//        std::vector<stk::mesh::Entity> firstElemNodes(bulk.begin_nodes(firstElem), bulk.end_nodes(firstElem));
//        std::vector<stk::mesh::Entity> secondElemNodes(bulk.begin_nodes(secondElem), bulk.end_nodes(secondElem));
//        std::sort(firstElemNodes.begin(), firstElemNodes.end());
//        std::sort(secondElemNodes.begin(), secondElemNodes.end());
//
//        std::vector<stk::mesh::Entity> commonNodes(std::max(firstElemNodes.size(), secondElemNodes.size()));
//        auto it = std::set_intersection(firstElemNodes.begin(), firstElemNodes.end(),
//                                        secondElemNodes.begin(), secondElemNodes.end(),
//                                        commonNodes.begin());
//        commonNodes.resize(it-commonNodes.begin());
//      }
//    }

    for (stk::mesh::Entity elem : secondBlockElems) {
      unsigned numNodes = bulk.num_connectivity(elem, stk::topology::NODE_RANK);
      const stk::mesh::Entity * elemNodes = bulk.begin(elem, stk::topology::NODE_RANK);
      stk::mesh::ConnectivityOrdinal const * nodeOrdinals = bulk.begin_ordinals(elem, stk::topology::NODE_RANK);
      for (unsigned iNode = 0; iNode < numNodes; ++iNode) {
        if (elemNodes[iNode] == nodeToDisconnect) {
          bulk.destroy_relation(elem, nodeToDisconnect, nodeOrdinals[iNode]);
          bulk.declare_relation(elem, newNode, nodeOrdinals[iNode]);
        }
      }
    }
  }
}

}
}
}
