// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_SubElementChildNodeAncestry_h
#define Akri_SubElementChildNodeAncestry_h

#include <stk_mesh/base/EntityKey.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <unordered_map>
#include <vector>

namespace krino
{
class CDMesh;
class SubElementNode;

class SubElementChildNodeAncestry {
public:
  SubElementChildNodeAncestry(const SubElementNode * node) { build_ancestry(node); }
  SubElementChildNodeAncestry( stk::CommBuffer & b );

  static bool is_shared(const stk::mesh::BulkData & mesh, const SubElementNode * node);
  void pack_into_buffer(stk::CommBuffer & b) const;
  const SubElementNode * find_subelement_node(CDMesh & mesh) const;
  const SubElementNode * find_subelement_node(CDMesh & mesh, unsigned & ancestry_index) const;
  void get_parent_node_keys(std::vector<stk::mesh::EntityKey> & parent_node_keys) const;
  void build_missing_child_nodes(CDMesh & mesh) const;
  std::vector<SubElementChildNodeAncestry> get_constrained_node_ancestries(const std::unordered_map<stk::mesh::EntityId, std::vector<stk::mesh::EntityId> > & constrained_node_map) const;

private:
  struct Cut {
    Cut(const std::vector<stk::mesh::EntityId> & parentIDs, const std::vector<double> & parentWeights) : myParentIDs(parentIDs), myParentWeights(parentWeights) {}
    std::vector<stk::mesh::EntityId> myParentIDs;
    std::vector<double> myParentWeights;
  };
  void build_ancestry(const SubElementNode * in_node);
  const SubElementNode * build_missing_child_nodes(CDMesh & mesh, unsigned & ancestry_index) const;

  std::vector<Cut> myAncestry;
};

}

#endif // Akri_SubElementChildNodeAncestry_h
