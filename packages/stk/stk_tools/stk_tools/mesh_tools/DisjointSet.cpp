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

#include "stk_tools/mesh_tools/DisjointSet.hpp"

#include <sstream>

namespace stk::experimental
{

NodeElemKey::operator std::string() const
{
  std::stringstream ss;
  ss << "(n" << get_node() << ", e" << get_elem() << ")";
  return ss.str();
}

void DisjointSet::fill_set(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
{
  auto elemBuckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, selector);
  for (const auto* bucket : elemBuckets) {
    for (const auto& elem : *bucket) {
      auto elemNodes = bulk.get_connected_entities(elem, stk::topology::NODE_RANK);
      for (auto n = 0U; n < elemNodes.size(); ++n) {
        insert(NodeElemKey(elemNodes[n], elem));
      }
    }
  }
  set_first_node_ids(bulk);
}

void DisjointSet::set_first_node_ids(const stk::mesh::BulkData& bulk)
{
  if (m_disjointNodes.empty()) return;

  auto* first = &(m_disjointNodes.begin()->second);
  first->nodeId = bulk.identifier(first->get_node());
  for (auto& [key, node] : m_disjointNodes) {
    if (first->get_node() != node.get_node()) {
      node.nodeId = bulk.identifier(node.get_node());
      first = &node;
    }
  }
}

void DisjointSet::merge_nodes(const NodeElemKey& a, const NodeElemKey& b)
{
  STK_ThrowRequireMsg(
      a.get_node() == b.get_node(), "Attempting to merge two disjoint sub-trees with different originating nodes!");

  auto& aRoot = find_root(a);
  auto& bRoot = find_root(b);
  if (aRoot < bRoot) {
    bRoot.parent = &aRoot;
  } else if (bRoot < aRoot) {
    aRoot.parent = &bRoot;
  }
}

std::size_t DisjointSet::count_trees() const
{
  std::size_t numTrees = 0U;
  for (const auto& [key, node] : m_disjointNodes) {
    if (node.parent == nullptr) {
      ++numTrees;
    }
  }
  return numTrees;
}

}  // namespace stk::experimental