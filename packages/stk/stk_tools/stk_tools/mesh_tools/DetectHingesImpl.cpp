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
#include "stk_tools/mesh_tools/DetectHinges.hpp"
#include "stk_tools/mesh_tools/EdgeMidNodeDetector.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_util/util/GraphCycleDetector.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
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

template <class ForwardIt, class T>
ForwardIt lower_bound_search(ForwardIt first, ForwardIt last, const T &value)
{
  first = std::lower_bound(first, last, value);
  if (!(first == last) && !(value < *first))
    return first;

  return last;
}

template <class ForwardIt, class T, class Compare>
ForwardIt lower_bound_search(ForwardIt first, ForwardIt last, const T &value, Compare comp)
{
  first = std::lower_bound(first, last, value, comp);
  if (!(first == last) && !(comp(value, *first)))
    return first;

  return last;
}

bool is_vertex_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem, const stk::mesh::ConnectivityOrdinal ordinal)
{
  auto topology = bulk.bucket(elem).topology();
  auto numVertices = topology.num_vertices();

  if(ordinal >= numVertices) {
    return false;
  }

  return true;
}

void populate_pairwise_side_info(const stk::mesh::BulkData& bulk, const PairwiseNodeElementRelation& rel,
                                 PairwiseSideInfoVector& infoVec)
{

  infoVec.emplace_back(bulk, rel.elem1, rel.elem2);

  if(infoVec.back().get_common_nodes().empty()) {
    infoVec.pop_back();
  }
}

bool is_auto_hinge(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  unsigned numConnectedElem = bulk.num_elements(node);
  const stk::mesh::Entity* elems = bulk.begin_elements(node);
  const stk::mesh::ConnectivityOrdinal * ordinals = bulk.begin_ordinals(node, stk::topology::ELEM_RANK);

  for(unsigned i = 0; i < numConnectedElem-1; i++) {
    for(unsigned j = i+1; j < numConnectedElem; j++) {
      if(is_vertex_node(bulk, elems[i], ordinals[i]) != is_vertex_node(bulk, elems[j], ordinals[j])) return true;
    }
  }

  return false;
}

void fill_pairwise_info_for_node(const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node,
                                 PairwiseSideInfoVector& infoVec, bool& autoHinge)
{
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();
  autoHinge = false;

  unsigned numConnectedElem = bulk.num_elements(node);
  if(numConnectedElem == 0) { return; }
  const stk::mesh::Entity* elems = bulk.begin_elements(node);
  const stk::mesh::ConnectivityOrdinal * ordinals = bulk.begin_ordinals(node, stk::topology::ELEM_RANK);

  for(unsigned i = 0; i < numConnectedElem-1; i++) {
    for(unsigned j = i+1; j < numConnectedElem; j++) {
      PairwiseNodeElementRelation rel(node, elems[i], elems[j]);

      populate_pairwise_side_info(bulk, rel, infoVec);

      autoHinge |= (is_vertex_node(bulk, elems[i], ordinals[i]) != is_vertex_node(bulk, elems[j], ordinals[j]));
    }
  }
}

void insert_connected_elements_uniquely(const stk::mesh::EntityVector& inputElems, stk::mesh::EntityVector& connectedElems)
{
  for(stk::mesh::Entity elem : inputElems) {
    stk::util::insert_keep_sorted_and_unique(elem, connectedElems);
  }
}

void insert_connected_elements_uniquely(const stk::mesh::BulkData& bulk, stk::mesh::Entity node, stk::mesh::EntityVector& connectedElems)
{
  unsigned numConnectedElem = bulk.num_elements(node);

  const stk::mesh::Entity* elems = bulk.begin_elements(node);
  for(unsigned i = 0; i < numConnectedElem; i++) {
    stk::util::insert_keep_sorted_and_unique(elems[i], connectedElems);
  }
}

void fill_pairwise_info_for_edge(const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node,
                                 PairwiseSideInfoVector& infoVec, bool& autoHinge)
{
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  stk::mesh::EntityVector connectedElems;
  stk::mesh::EntityVector commonElements;

  const std::vector<stk::tools::EdgeVertices>& edges = midNodeDetector.get_edge_vertices(node);
  for(const stk::tools::EdgeVertices& edge : edges) {
    stk::mesh::Entity nodes[] = {edge.first, edge.second};
    stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEM_RANK, 2, nodes, commonElements);

    insert_connected_elements_uniquely(commonElements, connectedElems);

    const stk::mesh::EntityVector& midNodes = midNodeDetector.get_mid_nodes(edge);
    for(stk::mesh::Entity midNode : midNodes) {
      insert_connected_elements_uniquely(bulk, midNode, connectedElems);
    }
  }

  autoHinge = is_auto_hinge(bulk, node);
  autoHinge |= (midNodeDetector.get_edge_vertices(node).size() > 1);

  unsigned numConnectedElem = connectedElems.size();
  if(numConnectedElem == 0) { return; }

  for(unsigned i = 0; i < numConnectedElem-1; i++) {
    for(unsigned j = i+1; j < numConnectedElem; j++) {
      PairwiseNodeElementRelation rel(node, connectedElems[i], connectedElems[j]);

      populate_pairwise_side_info(bulk, rel, infoVec);
    }
  }
}

void fill_common_nodes_for_connected_elems(const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node,
                                           PairwiseSideInfoVector& infoVec, bool& autoHinge)
{
  if(midNodeDetector.is_mid_edge_node(node)) {
    fill_pairwise_info_for_edge(midNodeDetector, node, infoVec, autoHinge);
  } else {
    fill_pairwise_info_for_node(midNodeDetector, node, infoVec, autoHinge);
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

bool common_nodes_are_part_of_an_edge(const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node1, stk::mesh::Entity node2)
{
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();
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

std::pair<PairwiseSideInfoVector, bool>
get_hinge_info_vec(const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node)
{
  PairwiseSideInfoVector infoVec;
  bool autoHinge = false;

  fill_common_nodes_for_connected_elems(midNodeDetector, node, infoVec, autoHinge);

  return std::make_pair(infoVec, autoHinge);
}

HingeNode convert_to_hinge_node (const EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node)
{
  PairwiseSideInfoVector infoVec;
  bool autoHinge{false};

  std::tie(infoVec, autoHinge) = get_hinge_info_vec(midNodeDetector, node);

  if(autoHinge) {
    // Vertex to non-vertex connection: automatic hinge
    return HingeNode(node, infoVec);
  }

  if(!autoHinge && midNodeDetector.is_mid_edge_node(node)) {
    // Mid node is part of exactly one unique edge ... pruning for hinge edge now depends on edge vertex status
    return HingeNode();
  }

  HingeGroupVector groupVec;
  insert_into_group(infoVec, groupVec);

  if(groupVec.size() >= 2) {
    return HingeNode(node, infoVec);
  }

  return HingeNode();
}

HingeNodeVector get_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector, const stk::mesh::EntityVector& nodes)
{
  HingeNodeVector hingeNodes;
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  for(stk::mesh::Entity node : nodes) {
    HingeNode hingeNode = convert_to_hinge_node(midNodeDetector, node);
    if(hingeNode.is_a_hinge()) {
      hingeNode.set_is_owned( hinge_node_is_locally_owned(bulk, hingeNode) );
      hingeNodes.push_back(hingeNode);
    }
  }

  std::sort(hingeNodes.begin(), hingeNodes.end());
  return hingeNodes;
}

bool is_connected_to_solid_element(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  stk::mesh::ConnectedEntities elems = bulk.get_connected_entities(node, stk::topology::ELEM_RANK);
  for(unsigned i=0; i<elems.size(); ++i) {
    if (stk::is_solid_element(bulk.bucket(elems[i]).topology())) {
      return true;
    }
  }

  return false;
}

stk::mesh::EntityVector get_mesh_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocks, bool onlyIfConnectedToSolidElements = false)
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

  if (onlyIfConnectedToSolidElements) {
    stk::mesh::EntityVector::iterator newEnd =
      std::remove_if(nodes.begin(), nodes.end(),
        [&](stk::mesh::Entity node) { return !is_connected_to_solid_element(bulk, node); });
    nodes.erase(newEnd, nodes.end());
  }

  return nodes;
}

HingeNodeVector get_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector)
{
  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(midNodeDetector.get_bulk_data(), stk::topology::NODE_RANK, nodes);

  HingeNodeVector hingeNodes = get_hinge_nodes(midNodeDetector, nodes);

  return hingeNodes;
}

HingeNodeVector get_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector,
                                const std::vector<std::string>& blocksToDetect,
                                bool onlyIfConnectedToSolidElements)
{
  stk::mesh::EntityVector nodes = get_mesh_nodes(midNodeDetector.get_bulk_data(), blocksToDetect, onlyIfConnectedToSolidElements);

  HingeNodeVector hingeNodes = get_hinge_nodes(midNodeDetector, nodes);

  return hingeNodes;
}

void fill_hinge_edges_for_hinge_node(const EdgeMidNodeDetector& midNodeDetector,
                                     const HingeNodeVector& hingeNodes,
                                     const HingeNode& hingeNode,
                                     HingeEdgeVector& hingeEdges)
{
  const PairwiseSideInfoVector& infoVec = hingeNode.get_info();
  stk::mesh::Entity currentNode = hingeNode.get_node();

  for(const PairwiseSideInfo& info : infoVec) {
    const stk::mesh::EntityVector& commonNodes = info.get_common_nodes();

    if(commonNodes.size() == 2) {
      stk::mesh::Entity otherEdgeNode = (commonNodes[0] == currentNode) ? commonNodes[1] : commonNodes[0];

      if(currentNode < otherEdgeNode) {

        if(common_nodes_are_part_of_an_edge(midNodeDetector, currentNode, otherEdgeNode)) {
          auto iter = lower_bound_search(hingeNodes.begin(), hingeNodes.end(), otherEdgeNode);
          if(iter != hingeNodes.end()) {
            stk::util::insert_keep_sorted_and_unique( {hingeNode, *iter}, hingeEdges);
          }
        }
      }
    }
  }
}


void fill_hinge_edges_from_node_list(const EdgeMidNodeDetector& midNodeDetector,
                                     const HingeNodeVector& hingeNodes,
                                     HingeEdgeVector& hingeEdges)
{
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();
  for(const HingeNode& hingeNode : hingeNodes) {
    fill_hinge_edges_for_hinge_node(bulk, hingeNodes, hingeNode, hingeEdges);
  }
}

void fill_hinge_edges_from_edge_list(const EdgeMidNodeDetector& midNodeDetector,
                                     const HingeNodeVector& hingeNodes,
                                     HingeEdgeVector& hingeEdges)
{
  std::vector<EdgeVertices> allEdges;
  midNodeDetector.fill_all_edges(allEdges);

  bool lastVertexIsHinge = false;
  stk::mesh::Entity lastVertex;

  for(const EdgeVertices& edge : allEdges) {
    if((edge.first == lastVertex) && !lastVertexIsHinge) {
      continue;
    }
    auto iter1 = lower_bound_search(hingeNodes.begin(), hingeNodes.end(), edge.first);

    if(iter1 != hingeNodes.end()) {
      const stk::mesh::EntityVector& midNodes = midNodeDetector.get_mid_nodes(edge);
      size_t numEdgeNodes = 2u + midNodes.size();

      auto iter2 = lower_bound_search(hingeNodes.begin(), hingeNodes.end(), edge.second);
      if(iter2 != hingeNodes.end()) {
        const PairwiseSideInfoVector& infoVec = iter1->get_info();
        for(const PairwiseSideInfo& info : infoVec) {
          const stk::mesh::EntityVector& commonNodes = info.get_common_nodes();

          if(commonNodes.size() == numEdgeNodes) {
            stk::mesh::Entity vertex2 = (commonNodes[0] == edge.first) ? commonNodes[1] : commonNodes[0];

            if(vertex2 == edge.second) {
              stk::util::insert_keep_sorted_and_unique( {*iter1, *iter2}, hingeEdges);
            }
          }
        }
      }

      lastVertexIsHinge = true;
    } else {
      lastVertexIsHinge = false;
    }

    lastVertex = edge.first;
  }
}

HingeEdgeVector get_hinge_edges(const EdgeMidNodeDetector& midNodeDetector, const HingeNodeVector& hingeNodes)
{
  HingeEdgeVector hingeEdges;

  fill_hinge_edges_from_edge_list(midNodeDetector, hingeNodes, hingeEdges);

  return hingeEdges;
}

void remove_entity_from_sorted_list(stk::mesh::EntityVector& entityVec, stk::mesh::Entity entity)
{
  auto it = lower_bound_search(entityVec.begin(), entityVec.end(), entity);
  if(it != entityVec.end()) {
    entityVec.erase(it);
  }
}

void prune_hinge_edge_node(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElems,
                           const HingeNode& hingeNode, HingeNodeVector& hingeNodes)
{
  stk::mesh::Entity node = hingeNode.get_node();
  stk::mesh::EntityVector nodeElems(bulk.begin_elements(node), bulk.begin_elements(node)+bulk.num_elements(node));
  std::sort(nodeElems.begin(), nodeElems.end());

  for(const PairwiseSideInfo& info : hingeNode.get_info()) {
    if(info.is_adjacent()) {
      remove_entity_from_sorted_list(nodeElems, info.get_element1());
      remove_entity_from_sorted_list(nodeElems, info.get_element2());
    }
  }

  for(stk::mesh::Entity elem : commonElems) {
    remove_entity_from_sorted_list(nodeElems, elem);
  }

  if(nodeElems.size() == 0) {
    auto it = lower_bound_search(hingeNodes.begin(), hingeNodes.end(), hingeNode);
    if(it != hingeNodes.end()) {
      hingeNodes.erase(it);
    }
  }
}

void prune_hinge_edge(const EdgeMidNodeDetector& midNodeDetector, const HingeEdge& hingeEdge, HingeNodeVector& hingeNodes)
{
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  stk::mesh::Entity node1 = hingeEdge.first.get_node();
  stk::mesh::Entity node2 = hingeEdge.second.get_node();
  stk::mesh::EntityVector commonElems = get_common_elements(bulk, node1, node2);

  prune_hinge_edge_node(bulk, commonElems, hingeEdge.first, hingeNodes);
  prune_hinge_edge_node(bulk, commonElems, hingeEdge.second, hingeNodes);

  stk::tools::EdgeVertices edge;
  if(bulk.identifier(node1) < bulk.identifier(node2)) {
    edge.first = node1;
    edge.second = node2;
  } else {
    edge.first = node2;
    edge.second = node1;
  }

  stk::mesh::EntityVector midNodes;
  midNodeDetector.fill_mid_nodes(edge, midNodes);

  for(stk::mesh::Entity node : midNodes) {
    auto it = lower_bound_search(hingeNodes.begin(), hingeNodes.end(), node);
    if(it != hingeNodes.end()) {
      hingeNodes.erase(it);
    }
  }
}

void prune_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector, HingeNodeVector& hingeNodes, const HingeEdgeVector& hingeEdges)
{
  for(const HingeEdge& edge : hingeEdges) {
    prune_hinge_edge(midNodeDetector, edge, hingeNodes);
  }
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

std::pair<unsigned, unsigned> get_hinge_count(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges)
{
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

std::pair<unsigned, unsigned> get_hinge_count(const EdgeMidNodeDetector& midNodeDetector)
{
  HingeNodeVector hingeNodes;
  HingeEdgeVector hingeEdges;

  impl::fill_mesh_hinges(midNodeDetector, std::vector<std::string>{}, hingeNodes, hingeEdges);

  return get_hinge_count(midNodeDetector.get_bulk_data(), hingeNodes, hingeEdges);
}

std::pair<unsigned, unsigned> get_hinge_count(const stk::mesh::BulkData& bulk)
{
  HingeNodeVector hingeNodes;
  HingeEdgeVector hingeEdges;

  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  return get_hinge_count(bulk, hingeNodes, hingeEdges);
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

void insert_into_group(const PairwiseSideInfoVector& node1InfoVec, const PairwiseSideInfoVector& /*node2InfoVec*/,
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


HingeGroupVector get_convex_groupings(const stk::tools::EdgeMidNodeDetector& midNodeDetector, stk::mesh::Entity node)
{
  HingeGroupVector groupVec;
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  if(!bulk.is_valid(node)) {
    return groupVec;
  }
  HingeNode hingeNode = convert_to_hinge_node(midNodeDetector, node);
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

HingeGroupVector get_convex_groupings(const stk::tools::EdgeMidNodeDetector& midNodeDetector, const HingeNode& node1, const HingeNode& node2)
{
  HingeGroupVector groupVec;
  const stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  if(!bulk.is_valid(node1.get_node()) || !bulk.is_valid(node2.get_node())) {
    return groupVec;
  }
  const PairwiseSideInfoVector& infoVec1 = node1.get_info();
  const PairwiseSideInfoVector& infoVec2 = node2.get_info();
  stk::mesh::EntityVector commonElem = get_common_elements(bulk, node1.get_node(), node2.get_node());
  insert_into_group(infoVec1, infoVec2, commonElem, groupVec);

  return groupVec;
}

HingeGroupVector get_convex_groupings(const stk::tools::EdgeMidNodeDetector& midNodeDetector, const HingeNode& node)
{
  return get_convex_groupings(midNodeDetector, node.get_node());
}

HingeGroupVector get_convex_groupings(const stk::tools::EdgeMidNodeDetector& midNodeDetector, const HingeEdge& edge)
{
  return get_convex_groupings(midNodeDetector, edge.first, edge.second);
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

void snip_all_hinges(stk::tools::EdgeMidNodeDetector& midNodeDetector, HingeNodeVector& hingeNodes)
{
  HingeGroupVector hingeGroups;
  std::vector<stk::mesh::EntityId> newNodeIdVec;
  std::pair<stk::mesh::Part*, stk::mesh::PartVector> hingeBlocks;

  LinkInfo info;

  stk::mesh::BulkData& bulk = midNodeDetector.get_bulk_data();

  bulk.modification_begin();

  for(HingeNode& hinge : hingeNodes) {
    hingeGroups = get_convex_groupings(midNodeDetector, hinge);

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

HingeNodeVector get_cyclic_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector, HingeNodeVector& hingeNodes)
{
  HingeEdgeVector hingeEdges = get_hinge_edges(midNodeDetector, hingeNodes);
  stk::mesh::impl::LocalIdMapper localIdMapper;
  localIdMapper.set_size(midNodeDetector.get_bulk_data());

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

void prune_hinge_nodes(HingeNodeVector& hingeNodes, const HingeNodeVector& excludedHingeNodes)
{
  for(auto node : excludedHingeNodes) {
    auto it = std::find(hingeNodes.begin(), hingeNodes.end(), node);

    if(it != hingeNodes.end()) {
      hingeNodes.erase(it);
    }
  }
}

void prune_cyclic_hinge_nodes(const EdgeMidNodeDetector& midNodeDetector, HingeNodeVector& hingeNodes)
{
  HingeNodeVector cyclicHingeNodes = get_cyclic_hinge_nodes(midNodeDetector, hingeNodes);
  prune_hinge_nodes(hingeNodes, cyclicHingeNodes);
}

void snip_all_hinges_for_input_nodes(EdgeMidNodeDetector& midNodeDetector,
                                     const stk::mesh::EntityVector nodes,
                                     const HingeNodeVector& preservedHingeNodes)
{
  HingeNodeVector hingeNodes = get_hinge_nodes(midNodeDetector, nodes);
  // prune_cyclic_hinge_nodes(midNodeDetector, hingeNodes);
  prune_hinge_nodes(hingeNodes, preservedHingeNodes);
  snip_all_hinges(midNodeDetector, hingeNodes);
}

void snip_all_hinges_for_input_nodes(EdgeMidNodeDetector& midNodeDetector, const stk::mesh::EntityVector nodes)
{
  snip_all_hinges_for_input_nodes(midNodeDetector, nodes, HingeNodeVector{});
}

void snip_all_hinges_between_blocks(stk::mesh::BulkData& bulk)
{
  EdgeMidNodeDetector midNodeDetector(bulk);
  HingeNodeVector hingeNodes = get_hinge_nodes(midNodeDetector);
  // prune_cyclic_hinge_nodes(midNodeDetector, hingeNodes);
  snip_all_hinges(midNodeDetector, hingeNodes);
}

void fill_mesh_hinges(const EdgeMidNodeDetector& midNodeDetector, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes)
{
  hingeNodes = impl::get_hinge_nodes(midNodeDetector, blocksToDetect);
}

void fill_mesh_hinges(const EdgeMidNodeDetector& midNodeDetector, const std::vector<std::string>& blocksToDetect,
                      HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges, bool onlyIfConnectedToSolidElements)
{
  hingeNodes = impl::get_hinge_nodes(midNodeDetector, blocksToDetect, onlyIfConnectedToSolidElements);

  if(hingeNodes.size() != 0) {
    hingeEdges = impl::get_hinge_edges(midNodeDetector, hingeNodes);
  }

  impl::prune_hinge_nodes(midNodeDetector, hingeNodes, hingeEdges);
}

} } }
