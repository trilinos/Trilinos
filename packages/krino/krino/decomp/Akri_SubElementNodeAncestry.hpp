// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_SubElementNodeAncestry_h
#define Akri_SubElementNodeAncestry_h

#include <algorithm>
#include <array>
#include <vector>

#include <Akri_SubElement.hpp>

namespace krino {

class SubElementNodeAncestry {
public:
  SubElementNodeAncestry() = default;
  SubElementNodeAncestry(const SubElementNode * node) : my_node(node) {}

  // Mesh nodes come before child nodes, and are ordered by compare (of SubElementNodes); child nodes
  // by their parents' ancestries, each sorted and compared lexicographically with the same order, at
  // every level.  The orders to use are in Akri_SubElementNodeOrder.hpp.
  template<class LESS>
  static bool compare(const SubElementNodeAncestry & x, const SubElementNodeAncestry & y, const LESS & compare);

  void print(std::ostream & os) const;

private:
  const SubElementNode * my_node = nullptr;
};

// A child node's parents as ancestries, sorted by less: in place for up to 4 parents (an edge, face or
// interior node), so that comparing ancestries does not allocate, and in a vector for more.
class SortedAncestriesOfParents
{
public:
  template<class ANCESTRY_LESS>
  SortedAncestriesOfParents(const SubElementNode & node, const ANCESTRY_LESS & less) : mySize(node.get_num_parents())
  {
    if (mySize > myInPlace.size()) myOnHeap.resize(mySize);
    SubElementNodeAncestry * parents = data();
    for (size_t i = 0; i < mySize; ++i)
      parents[i] = SubElementNodeAncestry(node.get_parent(i));
    std::sort(parents, parents + mySize, less);
  }
  const SubElementNodeAncestry * begin() const { return data(); }
  const SubElementNodeAncestry * end() const { return data() + mySize; }

private:
  SubElementNodeAncestry * data() { return (mySize > myInPlace.size()) ? myOnHeap.data() : myInPlace.data(); }
  const SubElementNodeAncestry * data() const { return (mySize > myInPlace.size()) ? myOnHeap.data() : myInPlace.data(); }
  size_t mySize;
  std::array<SubElementNodeAncestry, 4> myInPlace;
  std::vector<SubElementNodeAncestry> myOnHeap;
};

template<class LESS>
bool SubElementNodeAncestry::compare(const SubElementNodeAncestry & x, const SubElementNodeAncestry & y, const LESS & compare)
{
  // Decided without walking the ancestry where it can be, since walking it again at each level of
  // nesting costs twice as much per level: comparing the parents lexicographically compares a
  // shared parent with itself both ways, and sorting a node's parents compares a mesh node with a
  // child node.
  if (x.my_node == y.my_node) return false;
  const bool xIsMeshNode = x.my_node->is_mesh_node();
  const bool yIsMeshNode = y.my_node->is_mesh_node();
  if (xIsMeshNode || yIsMeshNode)
    return (xIsMeshNode && yIsMeshNode) ? compare(*(x.my_node), *(y.my_node)) : xIsMeshNode;
  const auto less = [&compare](const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
  { return SubElementNodeAncestry::compare(a, b, compare); };
  const SortedAncestriesOfParents xParents(*x.my_node, less);
  const SortedAncestriesOfParents yParents(*y.my_node, less);
  return std::lexicographical_compare(xParents.begin(), xParents.end(), yParents.begin(), yParents.end(), less);
}

inline std::ostream & operator << (std::ostream & os, const SubElementNodeAncestry & ancestry) { ancestry.print(os); return os; }

}

#endif // Akri_SubElementNodeAncestry_h
