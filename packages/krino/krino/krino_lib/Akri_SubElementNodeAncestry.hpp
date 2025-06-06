// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_SubElementNodeAncestry_h
#define Akri_SubElementNodeAncestry_h

#include <Akri_SubElement.hpp>

namespace krino {

class SubElementNodeAncestry {
public:
  SubElementNodeAncestry() = default;
  SubElementNodeAncestry(const SubElementNodeAncestry & rhs) = default;
  SubElementNodeAncestry(const SubElementNode * node) : my_node(node) {}

  // Friend declaration is required to build in c++20 mode or else the vector<SubElementNodeAncestry>::operator<
  // call in compare below fails to compile.
  friend inline bool operator<(const SubElementNodeAncestry & x, const SubElementNodeAncestry & y);

  template<class LESS>
  static bool compare(const SubElementNodeAncestry & x, const SubElementNodeAncestry & y, const LESS & compare)
  {
    const std::vector<SubElementNodeAncestry> xParents = x.get_parents();
    const std::vector<SubElementNodeAncestry> yParents = y.get_parents();
    if (xParents.empty())
    {
      if (yParents.empty()) return compare(*(x.my_node), *(y.my_node));
      else return true;
    }
    else if (yParents.empty()) return false;
    else return xParents < yParents;
  }

  void print(std::ostream & os) const;

private:
  std::vector<SubElementNodeAncestry> get_parents() const
  {
    std::vector<SubElementNodeAncestry> parents;
    if (!my_node->is_mesh_node())
    {
      const NodeVec nodeParents = my_node->get_parents();
      parents.reserve(nodeParents.size());
      for (auto && nodeParent : nodeParents)
        parents.emplace_back(nodeParent);
      std::sort(parents.begin(), parents.end());
    }
    return parents;
  }

  const SubElementNode * my_node = nullptr;
};

inline bool operator<(const SubElementNodeAncestry & x, const SubElementNodeAncestry & y)
{
  return SubElementNodeAncestry::compare(x,y,SubElementNode::less_by_entity_id);
}

inline std::ostream & operator << (std::ostream & os, const SubElementNodeAncestry & ancestry) { ancestry.print(os); return os; }

}

#endif // Akri_SubElementNodeAncestry_h
