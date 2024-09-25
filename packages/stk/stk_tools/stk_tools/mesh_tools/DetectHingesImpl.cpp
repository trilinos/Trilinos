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

#include "stk_tools/mesh_tools/DetectHingesImpl.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_util/util/GraphCycleDetector.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include <vector>

namespace stk {
namespace tools {
namespace impl {

void print_node_count(stk::mesh::BulkData& bulk, const std::string str)
{
  const unsigned nodeCount = stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part());

  std::cout << str << std::endl;
  std::cout << "p:" << bulk.parallel_rank() << " num nodes: " << nodeCount << std::endl;
}

bool are_equal(const stk::mesh::Entity* nodesPtr, const stk::mesh::EntityVector& nodesVec)
{
  for(stk::mesh::Entity node : nodesVec) {
    if (node != *nodesPtr) {
      return false;
    }
    ++nodesPtr;
  }
  return true;
}

class SideFinder
{
public:
  SideFinder(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
    : m_bulk(bulkData), m_elem(element)
  {
    const unsigned numNodes = m_bulk.num_nodes(m_elem);
    if (numNodes > MAX_NUM_NODES) {
      m_heapScratchSpace.resize(numNodes);
      m_scratchSpace = m_heapScratchSpace.data();
    }
    else {
      m_scratchSpace = m_stackScratchSpace;
    }
  }

  virtual ~SideFinder() {}

  bool put_nodes_in_side_order(stk::mesh::EntityVector& nodes)
  {
    stk::topology elementTopo = m_bulk.bucket(m_elem).topology();
    stk::mesh::EntityRank subRank = m_bulk.mesh_meta_data().side_rank();

    const unsigned numSubTopo = elementTopo.num_sub_topology(subRank);
    const unsigned numNodes = nodes.size();
    const stk::mesh::Entity* elemNodes = m_bulk.begin_nodes(m_elem);

    for(unsigned i = 0; i < numSubTopo; i++) {
      stk::topology subTopo = elementTopo.sub_topology(subRank, i);
      if(numNodes != subTopo.num_nodes()) {
        continue;
      }

      elementTopo.sub_topology_nodes(elemNodes, subRank, i, m_scratchSpace);
      std::sort(m_scratchSpace, m_scratchSpace+subTopo.num_nodes());

      if(are_equal(m_scratchSpace, nodes)) {
        elementTopo.sub_topology_nodes(elemNodes, subRank, i, nodes.data());
        return true;
      }
    }
    return false;
  }

private:
  const stk::mesh::BulkData& m_bulk;
  stk::mesh::Entity m_elem;
  static constexpr unsigned MAX_NUM_NODES = 32;
  stk::mesh::Entity m_stackScratchSpace[MAX_NUM_NODES];
  stk::mesh::EntityVector m_heapScratchSpace;
  stk::mesh::Entity* m_scratchSpace;
};

std::pair<stk::mesh::EntityVector,bool> get_pairwise_common_nodes(const stk::mesh::BulkData& bulk,
                                                                  stk::mesh::Entity elem1,
                                                                  stk::mesh::Entity elem2)
{
  std::pair<stk::mesh::EntityVector,bool> result;
  stk::mesh::EntityVector& commonNodes = result.first;
  commonNodes.assign(bulk.begin_nodes(elem1), bulk.end_nodes(elem1));

  const stk::mesh::Entity* elem2Nodes = bulk.begin_nodes(elem2);
  const unsigned numElem2Nodes = bulk.num_nodes(elem2);
  unsigned numIntersect = 0;
  const unsigned numElem1Nodes = commonNodes.size();
  for(unsigned n=0; n<numElem1Nodes; ++n) {
    for(unsigned m=0; m<numElem2Nodes; ++m) {
      if (commonNodes[n] == elem2Nodes[m]) {
        if (n > numIntersect) {
          commonNodes[numIntersect] = commonNodes[n];
        }
        ++numIntersect;
        break;
      }
    }
  }
  commonNodes.resize(numIntersect);
  std::sort(commonNodes.begin(), commonNodes.end());

  SideFinder sideFinder(bulk, elem1);
  bool commonNodesAreSideNodes = sideFinder.put_nodes_in_side_order(commonNodes);
  result.second = commonNodesAreSideNodes;

  return result;
}

stk::mesh::EntityVector get_common_elements(const stk::mesh::BulkData& bulk, stk::mesh::Entity node1, stk::mesh::Entity node2)
{
  stk::mesh::EntityVector commonElements;
  stk::mesh::EntityVector node1Elems(bulk.begin_elements(node1), bulk.begin_elements(node1)+bulk.num_elements(node1));
  stk::mesh::EntityVector node2Elems(bulk.begin_elements(node2), bulk.begin_elements(node2)+bulk.num_elements(node2));

  std::sort(node1Elems.begin(), node1Elems.end());
  std::sort(node2Elems.begin(), node2Elems.end());
  std::set_intersection(node1Elems.begin(), node1Elems.end(),
                        node2Elems.begin(), node2Elems.end(),
                        std::back_inserter(commonElements));

  return commonElements;
}


void populate_pairwise_side_info(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem1,
                                 stk::mesh::Entity elem2, PairwiseSideInfoVector& infoVec)
{
  infoVec.emplace_back(bulk, elem1, elem2);

  if(infoVec.back().get_common_nodes().empty()) {
    infoVec.pop_back();
  }
}

void fill_common_nodes_for_connected_elems(const stk::mesh::BulkData& bulk, stk::mesh::Entity node,
                                           PairwiseSideInfoVector& infoVec)
{
  unsigned numConnectedElem = bulk.num_elements(node);
  if(numConnectedElem == 0) { return; }
  const stk::mesh::Entity* elems = bulk.begin_elements(node);

  for(unsigned i = 0; i < numConnectedElem-1; i++) {
    for(unsigned j = i+1; j < numConnectedElem; j++) {
      populate_pairwise_side_info(bulk, elems[i], elems[j], infoVec);
    }
  }
}

bool nodes_are_elem_side_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem,
                               const stk::mesh::EntityVector& nodesToCheck)
{
  stk::topology elementTopo = bulk.bucket(elem).topology();

  const stk::mesh::EntityRank sideRank = bulk.mesh_meta_data().side_rank();
  const unsigned numSubTopo = elementTopo.num_sub_topology(sideRank);

  for(unsigned i = 0; i < numSubTopo; i++) {
    const stk::topology subTopo = elementTopo.sub_topology(sideRank, i);
    if(nodesToCheck.size() != subTopo.num_nodes()) {
      continue;
    }

    if(stk::mesh::is_side_equivalent(bulk, elem, i, nodesToCheck.data())) {
      return true;
    }
  }
  return false;
}

bool common_nodes_are_part_of_a_side (const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElements, const stk::mesh::EntityVector& commonNodes)
{
  for(stk::mesh::Entity elem : commonElements) {
    if (nodes_are_elem_side_nodes(bulk, elem, commonNodes)) {
      return true;
    }
  }
  return false;
}

bool common_nodes_are_part_of_an_edge(const stk::mesh::BulkData& bulk, stk::mesh::Entity node1, stk::mesh::Entity node2)
{
  stk::mesh::EntityVector commonElements = get_common_elements(bulk, node1, node2);

  if(commonElements.size() > 0) {

    stk::mesh::Entity elem = commonElements[0];
    stk::topology elementTopo = bulk.bucket(elem).topology();

    unsigned numSubTopo = elementTopo.num_sub_topology(stk::topology::EDGE_RANK);

    for(unsigned i = 0; i < numSubTopo; i++) {
      stk::topology subTopo = elementTopo.sub_topology(stk::topology::EDGE_RANK, i);
      if(2 != subTopo.num_nodes()) {
        continue;
      }

      if(stk::mesh::is_edge_equivalent(bulk, elem, i, stk::mesh::EntityVector{node1,node2}.data())) {
        return true;
      }
    }
  }

  return false;
}

unsigned get_side_count(const PairwiseSideInfoVector& infoVec)
{
  unsigned sideCount = 0;

  for(const PairwiseSideInfo& info : infoVec) {
    if(info.is_adjacent()) {
      sideCount++;
    }
  }

  return sideCount;
}

PairwiseSideInfoVector get_hinge_info_vec(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  PairwiseSideInfoVector infoVec;

  fill_common_nodes_for_connected_elems(bulk, node, infoVec);

  return infoVec;
}

HingeNode convert_to_hinge_node (const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  PairwiseSideInfoVector infoVec = get_hinge_info_vec(bulk, node);

  HingeGroupVector groupVec;
  insert_into_group(infoVec, groupVec);

  if(groupVec.size() >= 2) {
    return HingeNode(node, infoVec);
  }
  return HingeNode();
}

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& nodes)
{
  HingeNodeVector hingeNodes;

  for(stk::mesh::Entity node : nodes) {
    HingeNode hingeNode = convert_to_hinge_node(bulk, node);
    if(hingeNode.is_a_hinge()) {
      hingeNode.set_is_owned( hinge_node_is_locally_owned(bulk, hingeNode) );
      hingeNodes.push_back(hingeNode);
    }
  }
  return hingeNodes;
}

stk::mesh::EntityVector get_mesh_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocks)
{
  stk::mesh::EntityVector nodes;
  stk::mesh::Selector selector;

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::PartVector parts;
  for(const auto& block : blocks) {
    stk::mesh::Part* part = meta.get_part(block);
    if(nullptr != part &&
       part->primary_entity_rank() == stk::topology::ELEM_RANK &&
       part->id() != stk::mesh::Part::INVALID_ID) {
      parts.push_back(part);
    }

  }
  selector = stk::mesh::selectUnion(parts);

  if(parts.empty()) {
    selector = meta.universal_part();
  }
  selector &= (meta.locally_owned_part() | meta.globally_shared_part());

  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, selector, nodes);

  return nodes;
}

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk)
{
  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);

  HingeNodeVector hingeNodes = get_hinge_nodes(bulk, nodes);

  return hingeNodes;
}

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect)
{
  stk::mesh::EntityVector nodes = get_mesh_nodes(bulk, blocksToDetect);

  HingeNodeVector hingeNodes = get_hinge_nodes(bulk, nodes);

  return hingeNodes;
}

void fill_hinge_edges_for_hinge_node(const stk::mesh::BulkData& bulk, const HingeNodeVector& hingeNodes, const HingeNode& hingeNode, HingeEdgeVector& hingeEdges)
{
  stk::mesh::Entity currentNode;
  stk::mesh::Entity otherEdgeNode;

  PairwiseSideInfoVector infoVec = hingeNode.get_info();
  currentNode = hingeNode.get_node();

  for(PairwiseSideInfo info : infoVec) {
    const stk::mesh::EntityVector& commonNodes = info.get_common_nodes();

    if(commonNodes.size() == 2) {
      otherEdgeNode = (commonNodes[0] == currentNode) ? commonNodes[1] : commonNodes[0];
      if(currentNode < otherEdgeNode) {
        if(common_nodes_are_part_of_an_edge(bulk, currentNode, otherEdgeNode)) {
          auto iter = std::find(hingeNodes.begin(), hingeNodes.end(), otherEdgeNode);
          if(iter != hingeNodes.end()) {
            stk::util::insert_keep_sorted_and_unique( {hingeNode, *iter}, hingeEdges);
          }
        }
      }
    }
  }
}

HingeEdgeVector get_hinge_edges(const stk::mesh::BulkData& bulk, const HingeNodeVector& hingeNodes)
{
  HingeEdgeVector hingeEdges;

  for(const HingeNode& hingeNode : hingeNodes) {
    fill_hinge_edges_for_hinge_node(bulk, hingeNodes, hingeNode, hingeEdges);
  }

  return hingeEdges;
}

void remove_entity_from_list(stk::mesh::EntityVector& entityVec, stk::mesh::Entity entity)
{
  auto it = std::find(entityVec.begin(), entityVec.end(), entity);
  if(it != entityVec.end()) {
    entityVec.erase(it);
  }
}

void prune_hinge_edge_node(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElems, const HingeNode& hingeNode, HingeNodeVector& hingeNodes)
{
  stk::mesh::Entity node = hingeNode.get_node();
  stk::mesh::EntityVector nodeElems(bulk.begin_elements(node), bulk.begin_elements(node)+bulk.num_elements(node));
  std::sort(nodeElems.begin(), nodeElems.end());

  for(const PairwiseSideInfo& info : hingeNode.get_info()) {
    if(info.is_adjacent()) {
      remove_entity_from_list(nodeElems, info.get_element1());
      remove_entity_from_list(nodeElems, info.get_element2());
    }
  }

  for(stk::mesh::Entity elem : commonElems) {
    remove_entity_from_list(nodeElems, elem);
  }

  if(nodeElems.size() == 0) {
    auto it = std::find(hingeNodes.begin(), hingeNodes.end(), hingeNode);
    if(it != hingeNodes.end()) {
      hingeNodes.erase(it);
    }
  }
}

void prune_hinge_edge(const stk::mesh::BulkData& bulk, const HingeEdge& hingeEdge, HingeNodeVector& hingeNodes)
{
  stk::mesh::Entity node1 = hingeEdge.first.get_node();
  stk::mesh::Entity node2 = hingeEdge.second.get_node();
  stk::mesh::EntityVector commonElems = get_common_elements(bulk, node1, node2);

  prune_hinge_edge_node(bulk, commonElems, hingeEdge.first, hingeNodes);
  prune_hinge_edge_node(bulk, commonElems, hingeEdge.second, hingeNodes);
}

void prune_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, const HingeEdgeVector& hingeEdges)
{
  std::sort(hingeNodes.begin(), hingeNodes.end());

  for(const HingeEdge& edge : hingeEdges) {
    prune_hinge_edge(bulk, edge, hingeNodes);
  }
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes)
{
  std::vector<std::string> blocksToDetect;
  fill_mesh_hinges( bulk,  blocksToDetect, hingeNodes);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes)
{
  hingeNodes = get_hinge_nodes(bulk, blocksToDetect);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges)
{
  std::vector<std::string> blocksToDetect;
  fill_mesh_hinges( bulk,  blocksToDetect, hingeNodes, hingeEdges);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges)
{
  hingeNodes = get_hinge_nodes(bulk, blocksToDetect);

  if(hingeNodes.size() != 0) {
    hingeEdges = get_hinge_edges(bulk, hingeNodes);
  }

  prune_hinge_nodes(bulk, hingeNodes, hingeEdges);
}

bool hinge_node_is_locally_owned(const stk::mesh::BulkData& bulk, const HingeNode& node)
{
  return bulk.bucket(node.get_node()).owned();
}

bool hinge_edge_is_locally_owned(const stk::mesh::BulkData& bulk, const HingeEdge& edge)
{
  bool ownedFirstNode = bulk.bucket(edge.first.get_node()).owned();
  bool ownedSecondNode = bulk.bucket(edge.second.get_node()).owned();
  int otherNodeOwnerId = -1;

  if(ownedFirstNode || ownedSecondNode) {
    if(ownedFirstNode) {
      otherNodeOwnerId = bulk.parallel_owner_rank(edge.second.get_node());
    } else {
      otherNodeOwnerId = bulk.parallel_owner_rank(edge.first.get_node());
    }

    if(otherNodeOwnerId >= bulk.parallel_rank()) {
      return true;
    }
  }
  return false;
}

std::pair<unsigned, unsigned> get_hinge_count(const stk::mesh::BulkData& bulk)
{
  HingeNodeVector hingeNodes;
  HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  unsigned localHingeCount[2] = {0, 0};
  unsigned globalHingeCount[2];

  for(const HingeNode& node : hingeNodes) {
    if(node.is_owned()) {
      localHingeCount[0]++;
    }
  }
  for(const HingeEdge& edge : hingeEdges) {
    if(hinge_edge_is_locally_owned(bulk, edge)) {
      localHingeCount[1]++;
    }
  }

  stk::all_reduce_sum(bulk.parallel(), localHingeCount, globalHingeCount, 2);

  return std::make_pair(globalHingeCount[0], globalHingeCount[1]);
}

// Convex groupings
void merge_groups(HingeGroupVector& groupVec, int idx1, int idx2, const stk::mesh::EntityLess& compareLess)
{
  if(idx1 == idx2 || idx1 >= (int)groupVec.size() || idx2 >= (int)groupVec.size() ||
     idx1 < 0 || idx2 < 0) { return; }

  auto it1 = groupVec.begin()+idx1;
  auto it2 = groupVec.begin()+idx2;

  for(stk::mesh::Entity elem : *it2) {
    stk::util::insert_keep_sorted_and_unique( elem, *it1, compareLess );
  }
  groupVec.erase(it2);
}

int find_element_in_groups(const HingeGroupVector& groupVec, stk::mesh::Entity elem)
{
  for(unsigned i = 0; i < groupVec.size(); i++) {
    if(std::find(groupVec[i].begin(), groupVec[i].end(), elem) != groupVec[i].end()) {
      return i;
    }
  }
  return -1;
}

void populate_group(const PairwiseSideInfo& info, const int elem1Idx, const int elem2Idx, HingeGroupVector& groupVec)
{
  if (info.is_adjacent()) {
    stk::mesh::EntityLess compareLess(info.get_bulk());

    if (elem1Idx < 0 && elem2Idx < 0) {
      stk::util::insert_keep_sorted_and_unique({ info.get_element1(), info.get_element2() }, groupVec);
    } else if (elem1Idx >= 0 && elem2Idx >= 0) {
      merge_groups(groupVec, elem1Idx, elem2Idx, compareLess);
    } else if (elem1Idx < 0 && elem2Idx >= 0) {
      stk::util::insert_keep_sorted_and_unique(info.get_element1(), groupVec[elem2Idx], compareLess);
    } else {
      stk::util::insert_keep_sorted_and_unique(info.get_element2(), groupVec[elem1Idx], compareLess);
    }
  } else {
    stk::mesh::EntityLess compareLess(info.get_bulk());
    if (elem1Idx < 0) {
      stk::util::insert_keep_sorted_and_unique( { info.get_element1() }, groupVec);
    }
    if (elem2Idx < 0) {
      stk::util::insert_keep_sorted_and_unique( { info.get_element2() }, groupVec);
    }
  }
}

void insert_into_group(const PairwiseSideInfoVector& infoVec, HingeGroupVector& groupVec)
{
  for(const PairwiseSideInfo& info : infoVec) {
    int elem1Idx = find_element_in_groups(groupVec, info.get_element1());
    int elem2Idx = find_element_in_groups(groupVec, info.get_element2());

    populate_group(info, elem1Idx, elem2Idx, groupVec);
  }
}

void insert_into_group(const PairwiseSideInfoVector& node1InfoVec, const PairwiseSideInfoVector& node2InfoVec,
                       const stk::mesh::EntityVector& commonElem, HingeGroupVector& groupVec)
{
  for(const PairwiseSideInfo& info : node1InfoVec) {
    int elem1Idx = find_element_in_groups(groupVec, info.get_element1());
    int elem2Idx = find_element_in_groups(groupVec, info.get_element2());

    bool elem1IsCommon = std::find(commonElem.begin(), commonElem.end(), info.get_element1()) != commonElem.end();
    bool elem2IsCommon = std::find(commonElem.begin(), commonElem.end(), info.get_element2()) != commonElem.end();

    if(elem1IsCommon && elem2IsCommon) {
      populate_group(info, elem1Idx, elem2Idx, groupVec);
    }
  }
}


HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  HingeGroupVector groupVec;
  if(!bulk.is_valid(node)) {
    return groupVec;
  }
  HingeNode hingeNode = convert_to_hinge_node(bulk, node);
  const PairwiseSideInfoVector& infoVec = hingeNode.get_info();

  if(hingeNode.is_a_hinge()) {
    insert_into_group(infoVec, groupVec);
    STK_ThrowAssert(groupVec.size() != 0);
  } else {
    stk::mesh::EntityVector entityVec(bulk.begin_elements(node), bulk.begin_elements(node)+bulk.num_elements(node));
    if(entityVec.size() != 0) {
      groupVec.push_back(entityVec);
    }
  }

  return groupVec;
}

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeNode& node1, const HingeNode& node2)
{
  HingeGroupVector groupVec;
  if(!bulk.is_valid(node1.get_node()) || !bulk.is_valid(node2.get_node())) {
    return groupVec;
  }
  const PairwiseSideInfoVector& infoVec1 = node1.get_info();
  const PairwiseSideInfoVector& infoVec2 = node2.get_info();
  stk::mesh::EntityVector commonElem = get_common_elements(bulk, node1.get_node(), node2.get_node());
  insert_into_group(infoVec1, infoVec2, commonElem, groupVec);

  return groupVec;
}

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeNode& node)
{
  return get_convex_groupings(bulk, node.get_node());
}

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeEdge& edge)
{
  return get_convex_groupings(bulk, edge.first, edge.second);
}

stk::mesh::ConstPartVector get_blocks_for_hinge_group(const stk::mesh::BulkData& bulk, const HingeGroup& group)
{
  stk::mesh::ConstPartVector blocks;

  for(stk::mesh::Entity elem : group) {
    const stk::mesh::Part* part = get_block_part_for_element(bulk, elem);
    STK_ThrowRequire(part != nullptr);
    blocks.push_back(part);
  }

  stk::util::sort_and_unique(blocks, stk::mesh::PartLess());

  return blocks;
}

std::string indent(unsigned n)
{
  std::string indentation;

  for(unsigned i=0; i<n; ++i) {
    indentation += "    ";
  }
  return indentation;
}

void print_pairwise_side_info(const stk::mesh::BulkData& bulk, const PairwiseSideInfoVector& infoVec, unsigned indentLevel, std::ostringstream& os)
{
  for(const PairwiseSideInfo& info : infoVec) {
    os << indent(indentLevel)
       << "Element pair: {" << bulk.identifier(info.get_element1()) << ", " << bulk.identifier(info.get_element2()) << "}  "
       << "Is adjacent: " << info.is_adjacent() << std::endl;
  }
}

void print_hinge_node_info(const stk::mesh::BulkData& bulk, const HingeNode& hingeNode, unsigned indentLevel, std::ostringstream& os)
{
  os << indent(indentLevel+0) << "Hinge node global id: " << bulk.identifier(hingeNode.get_node()) << std::endl;
  os << indent(indentLevel+1) << "Is a hinge : " << hingeNode.is_a_hinge() << std::endl;
  os << indent(indentLevel+1) << "Is owned: " << hingeNode.is_owned() << std::endl;
  os << indent(indentLevel+1) << "Pairwise side info:" << std::endl;
  print_pairwise_side_info(bulk, hingeNode.get_info(), indentLevel+2, os);
}

void print_hinge_group_info(const stk::mesh::BulkData& bulk, const HingeGroupVector& hingeGroups, unsigned indentLevel, std::ostringstream& os)
{
  for(unsigned i=0; i< hingeGroups.size(); ++i) {
    const HingeGroup& group = hingeGroups[i];
    os << indent(indentLevel+0) << "Hinge group (" << i+1 << ")" << std::endl;
    for(stk::mesh::Entity entity : group) {
      os << indent(indentLevel+1) << bulk.entity_key(entity) << " : " << bulk.bucket(entity).topology() << std::endl;
    }
  }
}

void print_disconnect_group_info(const DisconnectGroup& group, unsigned indentLevel, std::ostringstream& os)
{
  os << indent(indentLevel+0) << "Disconnect group id: " << group.id() << std::endl;
  os << indent(indentLevel+1) << "Node: " << group.get_node_id() << std::endl;
  os << indent(indentLevel+1) << "Parts list:" << std::endl;
  for(const stk::mesh::Part* part : group.get_parts()) {
    os << indent(indentLevel+2) << part->name() << std::endl;
  }
  os << indent(indentLevel+1) << "Entity list:" << std::endl;
  for(stk::mesh::Entity entity : group.get_entities()) {
    os << indent(indentLevel+2) << group.get_bulk().entity_key(entity) << std::endl;
  }
}

void snip_all_hinges(stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes)
{
  HingeGroupVector hingeGroups;
  std::vector<stk::mesh::EntityId> newNodeIdVec;
  std::pair<stk::mesh::Part*, stk::mesh::PartVector> hingeBlocks;

  LinkInfo info;

  bulk.modification_begin();

  for(HingeNode& hinge : hingeNodes) {
    hingeGroups = get_convex_groupings(bulk, hinge);

    for (size_t firstGroupIdx = 0; firstGroupIdx < hingeGroups.size()-1; ++firstGroupIdx) {
      for (size_t secondGroupIdx = firstGroupIdx+1; secondGroupIdx < hingeGroups.size(); ++secondGroupIdx) {
        HingeGroup * firstGroup  = &hingeGroups[firstGroupIdx];
        HingeGroup * secondGroup = &hingeGroups[secondGroupIdx];

        if(bulk.identifier((*firstGroup)[0]) > bulk.identifier((*secondGroup)[0]))
        {
          firstGroup  = &hingeGroups[secondGroupIdx];
          secondGroup = &hingeGroups[firstGroupIdx];
        }

        stk::mesh::ConstPartVector blocksForSecondGroup = get_blocks_for_hinge_group(bulk, *secondGroup);
        DisconnectGroup group(bulk, blocksForSecondGroup, hinge.get_node());

        add_to_sharing_lookup(bulk, hinge.get_node(), info.sharedInfo);

        if(group.needs_to_communicate()) {
          NodeMapKey key(hinge.get_node(), group);
          info.clonedNodeMap[key] = NodeMapValue(bulk, hinge.get_node());
        } else {
          NodeMapKey key(hinge.get_node(), group);
          info.originalNodeMap[key] = NodeMapValue(bulk, hinge.get_node());
        }
      }
    }
  }

  create_new_duplicate_node_IDs(bulk, info);

  communicate_shared_node_information(bulk, info);

  for(auto mapEntry = info.clonedNodeMap.begin(); mapEntry != info.clonedNodeMap.end(); ++mapEntry) {
    const NodeMapKey& key = mapEntry->first;
    disconnect_elements(bulk, key, mapEntry->second, info);
  }
  info.flush(std::cout);
  bulk.modification_end();
}

HingeNodeVector get_cyclic_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes)
{
  HingeEdgeVector hingeEdges = get_hinge_edges(bulk, hingeNodes);
  stk::mesh::impl::LocalIdMapper localIdMapper;
  localIdMapper.set_size(bulk);

  for(unsigned i  = 0; i < hingeNodes.size(); i++) {
    localIdMapper.add_new_entity_with_local_id(hingeNodes[i].get_node(), i);
  }

  GraphCycleDetector hingeGraph(hingeNodes.size());

  for(auto hingeEdge : hingeEdges) {
    unsigned hingeNode1Id = localIdMapper.entity_to_local(hingeEdge.first.get_node());
    unsigned hingeNode2Id = localIdMapper.entity_to_local(hingeEdge.second.get_node());
    hingeGraph.add_edge(hingeNode1Id, hingeNode2Id);
  }

  HingeNodeVector cyclicHingeNodes;
  const std::vector<unsigned>& nodesInCycles = hingeGraph.get_nodes_in_cycles();

  for(unsigned nodeId : nodesInCycles) {
    cyclicHingeNodes.push_back(hingeNodes[nodeId]);
  }

  return cyclicHingeNodes;
}

void prune_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, const HingeNodeVector& prunedHingeNodes)
{
  for(auto node : prunedHingeNodes) {
    auto it = std::find(hingeNodes.begin(), hingeNodes.end(), node);

    if(it != hingeNodes.end()) {
      hingeNodes.erase(it);
    }
  }
}

void prune_cyclic_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes)
{
  HingeNodeVector cyclicHingeNodes = get_cyclic_hinge_nodes(bulk, hingeNodes);
  prune_hinge_nodes(bulk, hingeNodes, cyclicHingeNodes);
}

void snip_all_hinges_for_input_nodes(stk::mesh::BulkData& bulk, const stk::mesh::EntityVector nodes)
{
  snip_all_hinges_for_input_nodes(bulk, nodes, HingeNodeVector{});
}

void snip_all_hinges_for_input_nodes(stk::mesh::BulkData& bulk, const stk::mesh::EntityVector nodes,
                                     const HingeNodeVector& preservedHingeNodes)
{
  HingeNodeVector hingeNodes = get_hinge_nodes(bulk, nodes);
  // prune_cyclic_hinge_nodes(bulk, hingeNodes);
  prune_hinge_nodes(bulk, hingeNodes, preservedHingeNodes);
  snip_all_hinges(bulk, hingeNodes);
}

void snip_all_hinges_between_blocks(stk::mesh::BulkData& bulk)
{
  HingeNodeVector hingeNodes = get_hinge_nodes(bulk);
  // prune_cyclic_hinge_nodes(bulk, hingeNodes);
  snip_all_hinges(bulk, hingeNodes);
}

} } }
