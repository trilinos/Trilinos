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

#ifndef STK_STK_UTIL_STK_UTIL_UTIL_TREETRAVERSER_HPP_
#define STK_STK_UTIL_STK_UTIL_UTIL_TREETRAVERSER_HPP_

#include <algorithm>
#include <functional>
#include <string>  // for string
#include <unordered_map>
#include <vector>  // for vector

namespace stk
{
namespace util
{
template <typename Key, typename Iterator>
class TreeTraverser
{
 public:
  virtual Iterator begin() const = 0;

  virtual Iterator end() const = 0;

  virtual const Key& get_node(Iterator iter) const = 0;

  virtual size_t size() const = 0;

  virtual bool node_exists(const Key& node) const = 0;

  virtual const std::vector<Key>& children(const Key& node) const = 0;

  TreeTraverser() {}

  virtual ~TreeTraverser() {}

  bool is_cyclic(const Key& node) const
  {
    initialize();
    return check_for_cycle(node);
  }

  bool is_cyclic() const
  {
    if (begin() != end()) {
      return is_cyclic(get_node(begin()));
    }
    return false;
  }

  std::vector<Key>&& get_forward_traversal_list(const Key& node) const
  {
    initialize();
    fill_traversal(node);

    return std::move(m_traversalList);
  }

  std::vector<Key>&& get_reverse_traversal_list(const Key& node) const
  {
    initialize();
    fill_traversal(node);
    std::reverse(m_traversalList.begin(), m_traversalList.end());

    return std::move(m_traversalList);
  }

 private:
  TreeTraverser(const TreeTraverser&) = delete;

  void fill_traversal(const Key& node) const
  {
    if (node_exists(node)) {
      if (m_visitedNodes[node] == false) {
        m_visitedNodes[node] = true;
        m_traversalList.push_back(node);

        for (const Key& child : children(node)) {
          fill_traversal(child);
        }
      }
    }
  }

  bool check_for_cycle(const Key& node) const
  {
    bool isCyclic = false;
    if (node_exists(node)) {
      if (m_visitedNodes[node] == true) {
        isCyclic = true;
      } else {
        m_visitedNodes[node] = true;

        for (const Key& child : children(node)) {
          isCyclic |= check_for_cycle(child);
        }
      }
    }
    return isCyclic;
  }

  void initialize() const
  {
    m_traversalList.clear();
    m_traversalList.reserve(size());

    for (Iterator iter = begin(); iter != end(); ++iter) {
      const Key& node = get_node(iter);
      m_visitedNodes[node] = false;
    }
  }

  mutable std::unordered_map<Key, bool> m_visitedNodes;
  mutable std::vector<Key> m_traversalList;
};

}  // namespace util
}  // namespace stk

#endif /* STK_STK_UTIL_STK_UTIL_UTIL_TREETRAVERSER_HPP_ */
