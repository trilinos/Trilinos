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

#ifndef DISJOINT_SET_HPP
#define DISJOINT_SET_HPP

#include <cstddef>

#include "stk_mesh/base/BulkData.hpp"

namespace stk::experimental
{

struct NodeElemKey {
  constexpr NodeElemKey() = default;
  constexpr NodeElemKey(const stk::mesh::Entity& node_, const stk::mesh::Entity& elem_) : node(node_), elem(elem_) {}

  const auto& get_node() const { return node; }
  const auto& get_elem() const { return elem; }

  operator std::string() const;

  stk::mesh::Entity node{stk::mesh::Entity::InvalidEntity};
  stk::mesh::Entity elem{stk::mesh::Entity::InvalidEntity};
};

inline auto operator<(const NodeElemKey& a, const NodeElemKey& b)
{
  if (a.node == b.node) return a.elem < b.elem;
  return a.node < b.node;
}

inline auto operator==(const NodeElemKey& a, const NodeElemKey& b)
{
  return (a.node == b.node && a.elem == b.elem);
}

class DisjointSet
{
  struct DisjointNode {
    DisjointNode(const stk::mesh::Entity& node, const stk::mesh::Entity& elem) : key(node, elem), origEntity(node) {}
    DisjointNode(const NodeElemKey& key_) : key(key_), origEntity(key_.node) {}

    const auto& get_node() const { return key.get_node(); }
    const auto& get_elem() const { return key.get_elem(); }

    stk::mesh::EntityId get_node_id() const { return nodeId; }

    DisjointNode* parent = nullptr;
    NodeElemKey key{};
    stk::mesh::Entity origEntity{stk::mesh::Entity::InvalidEntity};
    stk::mesh::Entity entity{stk::mesh::Entity::InvalidEntity};
    stk::mesh::EntityId nodeId{stk::mesh::InvalidEntityId};
    bool isNew{false};
  };
  friend bool operator<(const DisjointNode& a, const DisjointNode& b);
  friend bool operator==(const DisjointNode& a, const DisjointNode& b);

  template <typename NodeT>
  static auto find_root_impl(NodeT&& node) -> decltype(auto)
  {
    auto current = &node;
    auto next = node.parent;
    while (next != nullptr) {
      current = next;
      next = next->parent;
    }
    return *current;
  }

 public:
  using container = std::map<NodeElemKey, DisjointNode>;
  using value_type = container::value_type;

  auto begin() { return m_disjointNodes.begin(); }
  auto end() { return m_disjointNodes.end(); }
  auto begin() const { return m_disjointNodes.begin(); }
  auto end() const { return m_disjointNodes.end(); }
  auto cbegin() const { return std::cbegin(m_disjointNodes); }
  auto cend() const { return std::cend(m_disjointNodes); }

  void insert(const NodeElemKey& key)
  {
    auto success = m_disjointNodes.insert(container::value_type(key, DisjointNode(key)));
    STK_ThrowRequireMsg(
        success.second, "Unable to insert key " + std::string(key) + "!  Disjoint node already exists.");
  }

  template <typename ContainerT>
  void fill_set(const stk::mesh::BulkData& bulk, const ContainerT& elementList)
  {
    for (const auto& elem : elementList) {
      auto elemNodes = bulk.get_connected_entities(elem, stk::topology::NODE_RANK);
      for (auto n = 0U; n < elemNodes.size(); ++n) {
        insert(NodeElemKey(elemNodes[n], elem));
      }
    }
    set_first_node_ids(bulk);
  }
  void fill_set(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector);
  void set_first_node_ids(const stk::mesh::BulkData& bulk);

  const DisjointNode& operator[](const NodeElemKey& key) const
  {
    auto iter = m_disjointNodes.find(key);
    STK_ThrowRequireMsg(iter != m_disjointNodes.end(), key_not_found_string(key));
    return iter->second;
  }

  DisjointNode& operator[](const NodeElemKey& key)
  {
    auto iter = m_disjointNodes.find(key);
    STK_ThrowRequireMsg(iter != m_disjointNodes.end(), key_not_found_string(key));
    return iter->second;
  }

  const DisjointNode& find_root(const NodeElemKey& key) const
  {
    const auto& node = operator[](key);
    return find_root_impl(node);
  }

  DisjointNode& find_root(const NodeElemKey& key)
  {
    auto& node = operator[](key);
    return find_root_impl(node);
  }

  void merge_nodes(const NodeElemKey& a, const NodeElemKey& b);

  std::size_t count_trees() const;

 private:
  static std::string key_not_found_string(const NodeElemKey& key) { return "Key {} not found!" + std::string(key); }

  container m_disjointNodes;
};

inline bool operator<(const DisjointSet::DisjointNode& a, const DisjointSet::DisjointNode& b)
{
  return a.key < b.key;
}

inline bool operator==(const DisjointSet::DisjointNode& a, const DisjointSet::DisjointNode& b)
{
  return a.key == b.key;
}

}  // namespace stk::experimental

#endif
