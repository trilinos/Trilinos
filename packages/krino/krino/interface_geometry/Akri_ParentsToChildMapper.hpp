// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_PARENTSTOCHILDMAPPER_H_
#define AKRI_PARENTSTOCHILDMAPPER_H_

#include <stk_mesh/base/BulkData.hpp>

#include <Akri_OrderedIdPair.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino { class FieldRef; }

namespace krino {

class CDFEM_Support;

class ParentsToChildMapper
{
public:
  ParentsToChildMapper() = default;

  void build_map(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const CDFEM_Support & cdfemSupport, const bool addHigherOrderMidSideNodes);
  void build_map(const stk::mesh::BulkData & mesh, const FieldRef & parent_ids_field, const stk::mesh::Selector & elementSelector, const bool add_higher_order_midside_nodes);

  bool get_child_id(const stk::mesh::EntityId parent0,
      const stk::mesh::EntityId parent1,
      stk::mesh::EntityId & child) const
  {
    const auto it = my_parents_to_child_map.find({parent0, parent1});
    if(it == my_parents_to_child_map.end()) return false;
    child = it->second;
    return true;
  }

  stk::mesh::Entity get_child(const stk::mesh::BulkData & mesh, const stk::mesh::Entity parent0,
      const stk::mesh::Entity parent1) const
  {
    const auto it =
        my_parents_to_child_map.find({mesh.identifier(parent0), mesh.identifier(parent1)});
    if(it == my_parents_to_child_map.end()) return stk::mesh::Entity();
    return mesh.get_entity(stk::topology::NODE_RANK, it->second);
  }

private:
  bool need_to_update_map(const stk::mesh::BulkData & mesh, const bool addHigherOrderMidSideNodes) const;
  void mark_map_as_up_to_date(const stk::mesh::BulkData & mesh, const bool addHigherOrderMidSideNodes);
  typedef OrderedIdPair ParentsToChildKey;
  typedef std::map<ParentsToChildKey,stk::mesh::EntityId> ParentsToChildMap;
  ParentsToChildMap my_parents_to_child_map;
  size_t myMeshSyncCount{0};
  bool myHaveHigherOrderMidSideNodes{false};
};

void fill_edge_nodes(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node0,
    const stk::mesh::Entity node1,
    const ParentsToChildMapper & parentToChildMapper,
    std::vector<stk::mesh::Entity> & edgeNodes);

void fill_edge_nodes_and_positions(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node0, stk::mesh::Entity node1,
    const ParentsToChildMapper & parent_child_mapper,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions);

void fill_linear_edge_nodes_and_positions(stk::mesh::Entity node0,
    stk::mesh::Entity node1,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions);

} // namespace krino

#endif /* AKRI_PARENTSTOCHILDMAPPER_H_ */
