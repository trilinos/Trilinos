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
#include "stk_mesh/base/BulkData.hpp"
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
  stk::mesh::EntityVector nodes;
  bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), nodes);

  std::cout << str << std::endl;
  std::cout << "p:" << bulk.parallel_rank() << " node vec size: " << nodes.size() << std::endl;
}

void fill_pairwise_common_nodes_by_face(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem, stk::mesh::EntityVector& commonNodes)
{
  if(bulk.mesh_meta_data().spatial_dimension() != 3) {
    return;
  }

  stk::mesh::EntityVector elemSubTopoNodes(bulk.num_nodes(elem));
  stk::mesh::EntityVector sortedElemSubTopoNodes(bulk.num_nodes(elem));
  stk::topology elementTopo = bulk.bucket(elem).topology();
  stk::mesh::EntityRank subRank = bulk.mesh_meta_data().side_rank();

  unsigned numSubTopo = elementTopo.num_sub_topology(subRank);

  for(unsigned i = 0; i < numSubTopo; i++) {
    stk::topology subTopo = elementTopo.sub_topology(subRank, i);
    if(commonNodes.size() != subTopo.num_nodes()) {
      continue;
    }

    elementTopo.sub_topology_nodes(bulk.begin_nodes(elem), subRank, i, elemSubTopoNodes.begin());
    elemSubTopoNodes.resize(commonNodes.size());
    sortedElemSubTopoNodes = elemSubTopoNodes;
    std::sort(sortedElemSubTopoNodes.begin(), sortedElemSubTopoNodes.end());

    if(sortedElemSubTopoNodes == commonNodes) {
      elemSubTopoNodes.swap(commonNodes);
      return;
    }
  }
}

stk::mesh::EntityVector get_pairwise_common_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem1,
                                         stk::mesh::Entity elem2)
{
  stk::mesh::EntityVector commonNodes;
  stk::mesh::EntityVector sortedElem1Nodes(bulk.begin_nodes(elem1), bulk.begin_nodes(elem1)+bulk.num_nodes(elem1));
  stk::mesh::EntityVector sortedElem2Nodes(bulk.begin_nodes(elem2), bulk.begin_nodes(elem2)+bulk.num_nodes(elem2));
  std::sort(sortedElem1Nodes.begin(), sortedElem1Nodes.end());
  std::sort(sortedElem2Nodes.begin(), sortedElem2Nodes.end());
  std::set_intersection(sortedElem1Nodes.begin(), sortedElem1Nodes.end(),
                        sortedElem2Nodes.begin(), sortedElem2Nodes.end(),
                        std::back_inserter(commonNodes));

  fill_pairwise_common_nodes_by_face(bulk, elem1, commonNodes);
  return commonNodes;
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
    PairwiseSideInfo info(bulk, elem1, elem2);

    if(info.get_common_nodes().size() > 0) {
      infoVec.push_back(info);
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

bool common_nodes_are_part_of_a_side (const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElements, const stk::mesh::EntityVector commonNodes)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank sideRank = meta.side_rank();
  for(stk::mesh::Entity elem : commonElements) {
    stk::topology elementTopo = bulk.bucket(elem).topology();

    unsigned numSubTopo = elementTopo.num_sub_topology(sideRank);

    for(unsigned i = 0; i < numSubTopo; i++) {
      stk::topology subTopo = elementTopo.sub_topology(sideRank, i);
      if(commonNodes.size() != subTopo.num_nodes()) {
        continue;
      }

      if(stk::mesh::is_side_equivalent(bulk, elem, i, commonNodes.data())) {
        return true;
      }
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

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk)
{
  stk::mesh::EntityVector nodes;
  bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), nodes);
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

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges)
{
  hingeNodes = get_hinge_nodes(bulk);

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
//  std::ostringstream os;
  for(const PairwiseSideInfo& info : infoVec) {
    int elem1Idx = find_element_in_groups(groupVec, info.get_element1());
    int elem2Idx = find_element_in_groups(groupVec, info.get_element2());

//    const stk::mesh::BulkData& bulk = info.get_bulk();
//    os << "P" << bulk.parallel_rank() << ": Elements{" << bulk.identifier(info.get_element1()) << "," << bulk.identifier(info.get_element2()) << "} : Indices{" << elem1Idx << "," << elem2Idx << "}" << std::endl;

    populate_group(info, elem1Idx, elem2Idx, groupVec);
  }
//  std::cout << os.str();
}

void insert_into_group(const PairwiseSideInfoVector& node1InfoVec, const PairwiseSideInfoVector& node2InfoVec,
                       stk::mesh::EntityVector commonElem, HingeGroupVector& groupVec)
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
    ThrowAssert(groupVec.size() != 0);
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
    ThrowRequire(part != nullptr);
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

void snip_all_hinges_between_blocks(stk::mesh::BulkData& bulk, bool debug)
{
  HingeNodeVector hingeNodes = get_hinge_nodes(bulk);
  HingeEdgeVector hingeEdges = get_hinge_edges(bulk,hingeNodes);
  HingeGroupVector hingeGroups;
  std::vector<stk::mesh::EntityId> newNodeIdVec;
  std::pair<stk::mesh::Part*, stk::mesh::PartVector> hingeBlocks;

  LinkInfo info;
  std::ostringstream& os = info.os;

  if(debug) {
    os << "P" << bulk.parallel_rank() << std::endl;
  }

  bulk.modification_begin();

  for(HingeNode& hinge : hingeNodes) {
    hingeGroups = get_convex_groupings(bulk, hinge);

    if(debug) {
      print_hinge_node_info(bulk, hinge, 1u, os);
      print_hinge_group_info(bulk, hingeGroups, 2u, os);
    }

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

        if(debug) {
          os << indent(1u) << "Creating disconnect group for hinge group with id " << group.id() << std::endl;
        }

        add_to_sharing_lookup(bulk, hinge.get_node(), info.sharedInfo);

        if(group.needs_to_communicate()) {
          if(debug) {
            os << indent(1u) << "Group needs to communicate ... adding to node map" << std::endl;
          }
          NodeMapKey key(hinge.get_node(), group);
          info.clonedNodeMap[key] = NodeMapValue(bulk, hinge.get_node());
        } else {
          NodeMapKey key(hinge.get_node(), group);
          info.preservedNodeMap[key] = NodeMapValue(bulk, hinge.get_node());
        }
      }
    }
  }

  create_new_duplicate_node_IDs(bulk, info);

  if(debug) {
    for (auto & nodeMapEntry : info.clonedNodeMap) {
      os << indent(1u) << "Disconnect group with id " << nodeMapEntry.first.disconnectedGroup.id()
                       << " has duplicate node id: " << nodeMapEntry.second.newNodeId << std::endl;
    }
  }

  communicate_shared_node_information(bulk, info);

  for(auto mapEntry = info.clonedNodeMap.begin(); mapEntry != info.clonedNodeMap.end(); ++mapEntry) {
    const NodeMapKey& key = mapEntry->first;

    if(debug) {
      print_disconnect_group_info(key.disconnectedGroup, 1u, os);
      os << indent(2u) << "New node id: " << mapEntry->second.newNodeId << std::endl;
    }
    disconnect_elements(bulk, key, mapEntry->second, info);
  }
  info.flush(std::cout);
  bulk.modification_end();
}

} } }
