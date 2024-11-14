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

#ifndef _DetectHingesImpl_hpp_
#define _DetectHingesImpl_hpp_

#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <vector>

namespace stk {
namespace mesh { class BulkData; }
namespace tools {
namespace impl {

bool common_nodes_are_part_of_a_side (const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElements, const stk::mesh::EntityVector& commonNodes);
std::pair<stk::mesh::EntityVector,bool> get_pairwise_common_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem1, stk::mesh::Entity elem2);

class PairwiseSideInfo {

public:
  PairwiseSideInfo(const stk::mesh::BulkData& bulk_, stk::mesh::Entity elem1_, stk::mesh::Entity elem2_)
    : bulk(&bulk_), elem1(elem1_), elem2(elem2_)
  {
    std::tie(commonNodes, hasAdjacentFace) = get_pairwise_common_nodes(*bulk, elem1, elem2);
  }

  stk::mesh::Entity get_element1() const { return elem1; }
  stk::mesh::Entity get_element2() const { return elem2; }
  const stk::mesh::EntityVector& get_common_nodes() { return commonNodes; }
  bool is_adjacent() const { return hasAdjacentFace; }
  void set_adjacency(bool adjacent) { hasAdjacentFace = adjacent; }

  const stk::mesh::BulkData& get_bulk() const {
    STK_ThrowRequire(nullptr != bulk);
    return *bulk;
  }

private:
  const stk::mesh::BulkData* bulk = nullptr;
  stk::mesh::Entity elem1;
  stk::mesh::Entity elem2;
  stk::mesh::EntityVector commonNodes;
  bool hasAdjacentFace = false;
};

typedef std::vector<PairwiseSideInfo> PairwiseSideInfoVector;

class HingeNode {
public:
  HingeNode() { }
  HingeNode(const HingeNode& h) : node(h.node), info(h.info), isAHinge(h.isAHinge), isOwned(h.isOwned) {}
  HingeNode(stk::mesh::Entity node_, const PairwiseSideInfoVector& info_) : node(node_), info(info_), isAHinge(true), isOwned(false) { }

  stk::mesh::Entity get_node() const { return node; }
  bool is_a_hinge() const { return isAHinge; }
  const PairwiseSideInfoVector& get_info() const { return info; }
  bool is_owned() const { return isOwned; }

  void set_is_owned(bool owned) { isOwned = owned; }

  bool operator() (const HingeNode& h1, const HingeNode& h2) const { return h1.node < h2.node; }
  bool operator== (const stk::mesh::Entity entity) const { return node == entity; }
  bool operator== (const HingeNode& h) const { return node == h.node; }
  bool operator<  (const HingeNode& h) const { return node < h.node; }
  HingeNode& operator=  (const HingeNode& h) { node = h.node; info = h.info; isAHinge = h.isAHinge; isOwned = h.isOwned; return *this; }

private:
  stk::mesh::Entity node;
  PairwiseSideInfoVector info;
  bool isAHinge = false;
  bool isOwned = false;
};

typedef std::vector<HingeNode> HingeNodeVector;
typedef std::pair<HingeNode, HingeNode> HingeEdge;
typedef std::vector<HingeEdge> HingeEdgeVector;
typedef stk::mesh::EntityVector HingeGroup;
typedef std::vector<HingeGroup> HingeGroupVector;

// Detecting Hinges
bool fill_pairwise_common_nodes_by_face(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem, stk::mesh::EntityVector& commonNodes);

stk::mesh::EntityVector get_common_elements(const stk::mesh::BulkData& bulk, stk::mesh::Entity node1, stk::mesh::Entity node2);

void populate_pairwise_side_info(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem1,
                                 stk::mesh::Entity elem2, PairwiseSideInfoVector& infoVec);

void fill_common_nodes_for_connected_elems(const stk::mesh::BulkData& bulk, stk::mesh::Entity node,
                                           PairwiseSideInfoVector& infoVec);

bool common_nodes_are_part_of_an_edge(const stk::mesh::BulkData& bulk, stk::mesh::Entity node1, stk::mesh::Entity node2);

unsigned get_side_count(const PairwiseSideInfoVector& infoVec) ;

PairwiseSideInfoVector get_hinge_info_vec(const stk::mesh::BulkData& bulk, stk::mesh::Entity node);

HingeNode convert_to_hinge_node (const stk::mesh::BulkData& bulk, stk::mesh::Entity node);

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& nodes);

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk);

HingeNodeVector get_hinge_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect);

void fill_hinge_edges_for_hinge_node(const stk::mesh::BulkData& bulk, const HingeNodeVector& hingeNodes, const HingeNode& hingeNode, HingeEdgeVector& hingeEdges);

HingeEdgeVector get_hinge_edges(const stk::mesh::BulkData& bulk, const HingeNodeVector& hingeNodes);

void remove_entity_from_list(stk::mesh::EntityVector& entityVec, stk::mesh::Entity entity);

void prune_hinge_edge_node(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& commonElems, const HingeNode& hingeNode, HingeNodeVector& hingeNodes);

void prune_hinge_edge(const stk::mesh::BulkData& bulk, const HingeEdge& hingeEdge, HingeNodeVector& hingeNodes);

void prune_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, const HingeEdgeVector& hingeEdges);

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges);

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes);

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges);

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes);

bool hinge_node_is_locally_owned(const stk::mesh::BulkData& bulk, const HingeNode& node);

bool hinge_edge_is_locally_owned(const stk::mesh::BulkData& bulk, const HingeEdge& edge);

std::pair<unsigned, unsigned> get_hinge_count(const stk::mesh::BulkData& bulk);

HingeNodeVector get_cyclic_hinge_nodes(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes);

// Convex Groupings
void merge_groups(HingeGroupVector& groupVec, int idx1, int idx2, const stk::mesh::EntityLess& compare);

int find_element_in_groups(const HingeGroupVector& groupVec, stk::mesh::Entity elem);

void populate_group(const PairwiseSideInfo& info, const int elem1Idx, const int elem2Idx, HingeGroupVector& groupVec);

void insert_into_group(const PairwiseSideInfoVector& infoVec, HingeGroupVector& groupVec);

void insert_into_group(const PairwiseSideInfoVector& node1InfoVec, const PairwiseSideInfoVector& node2InfoVec,
                       const stk::mesh::EntityVector& commonElem, HingeGroupVector& groupVec);

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, stk::mesh::Entity node);

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeNode& node1, const HingeNode& node2);

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeNode& node);

HingeGroupVector get_convex_groupings(const stk::mesh::BulkData& bulk, const HingeEdge& edge);


// Hinge Snipping
void snip_all_hinges_for_input_nodes(stk::mesh::BulkData& bulk, const stk::mesh::EntityVector nodes);

void snip_all_hinges_for_input_nodes(stk::mesh::BulkData& bulk, const stk::mesh::EntityVector nodes,
                                     const HingeNodeVector& preservedHingeNodes);

void snip_all_hinges_between_blocks(stk::mesh::BulkData& bulk);

}}}

#endif
