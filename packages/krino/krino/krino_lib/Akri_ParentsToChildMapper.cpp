// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ParentsToChildMapper.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <Akri_MasterElementDeterminer.hpp>

namespace krino {

static bool determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  // Assume we can check master element on one bucket and get the right answer

  for ( auto && bucket_ptr : mesh.get_buckets(stk::topology::ELEMENT_RANK, phaseSupport.get_all_decomposed_blocks_selector()) )
  {
    if (bucket_ptr->topology() != bucket_ptr->topology().base())
    {
      for (int ls_index = 0; ls_index < cdfemSupport.num_ls_fields(); ++ls_index)
      {
        stk::topology fieldTopology = MasterElementDeterminer::get_field_topology(*bucket_ptr, cdfemSupport.ls_field(ls_index).isovar);
        if (fieldTopology != fieldTopology.base())
          return true;
      }
    }
  }
  return false;
}

void fill_edge_nodes(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node0,
    const stk::mesh::Entity node1,
    const ParentsToChildMapper & parentToChildMapper,
    std::vector<stk::mesh::Entity> & edgeNodes)
{
  stk::mesh::Entity child = parentToChildMapper.get_child(mesh, node0, node1);

  if (mesh.is_valid(child))
  {
    ThrowRequire(child != node0);
    ThrowRequire(child != node1);
    fill_edge_nodes(mesh, node0,child, parentToChildMapper, edgeNodes);
    fill_edge_nodes(mesh, child,node1, parentToChildMapper, edgeNodes);
  }
  else
  {
    if (edgeNodes.empty())
    {
      edgeNodes.push_back(node0);
    }
    edgeNodes.push_back(node1);
  }
}

static void fill_edge_nodes_and_positions(
    const double pos0, stk::mesh::Entity node0,
    const double pos1, stk::mesh::Entity node1,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions)
{
  if (edgeNodePositions.empty())
  {
    edgeNodePositions.push_back(pos0);
    edgeNodes.push_back(node0);
  }
  edgeNodePositions.push_back(pos1);
  edgeNodes.push_back(node1);
}

static void fill_edge_nodes_and_positions(
    const stk::mesh::BulkData & mesh,
    const double pos0, stk::mesh::Entity node0,
    const double pos1, stk::mesh::Entity node1,
    const ParentsToChildMapper & parent_child_mapper,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions)
{
  stk::mesh::Entity child = parent_child_mapper.get_child(mesh, node0, node1);

  if (mesh.is_valid(child))
  {
    ThrowRequire(child != node0);
    ThrowRequire(child != node1);
    const double child_rel_pos = compute_child_position(mesh, child, node0, node1);
    const double child_abs_pos = pos0 * (1.-child_rel_pos) + pos1 * child_rel_pos;
    fill_edge_nodes_and_positions(mesh,pos0,node0,child_abs_pos,child, parent_child_mapper, edgeNodes, edgeNodePositions);
    fill_edge_nodes_and_positions(mesh,child_abs_pos,child,pos1,node1, parent_child_mapper, edgeNodes, edgeNodePositions);
  }
  else
  {
    fill_edge_nodes_and_positions(pos0, node0, pos1, node1, edgeNodes, edgeNodePositions);
  }
}

void fill_edge_nodes_and_positions(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node0, stk::mesh::Entity node1,
    const ParentsToChildMapper & parent_child_mapper,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions)
{
  edgeNodes.clear();
  edgeNodePositions.clear();
  fill_edge_nodes_and_positions(mesh, 0.0, node0, 1.0, node1, parent_child_mapper, edgeNodes, edgeNodePositions);
}

void fill_linear_edge_nodes_and_positions(stk::mesh::Entity node0, stk::mesh::Entity node1,
    std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<double> & edgeNodePositions)
{
  edgeNodes.clear();
  edgeNodePositions.clear();
  fill_edge_nodes_and_positions(0.0, node0, 1.0, node1, edgeNodes, edgeNodePositions);
}


void ParentsToChildMapper::build_map(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  const bool addHigherOrderMidsideNodes = stk::is_true_on_any_proc(mesh.parallel(), determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(mesh, cdfemSupport, phaseSupport));

  stk::mesh::Selector elementSelectorForParentChildMapper = activePart & mesh.mesh_meta_data().locally_owned_part();

  build_map(mesh, cdfemSupport.get_parent_node_ids_field(), elementSelectorForParentChildMapper, addHigherOrderMidsideNodes);
}

void ParentsToChildMapper::build_map(const stk::mesh::BulkData & mesh, const FieldRef & parent_ids_field, const stk::mesh::Selector & elementSelector, bool add_higher_order_midside_nodes)
{
  my_parents_to_child_map.clear();

  if (parent_ids_field.valid())
  {
    const auto & field_selector = stk::mesh::selectField(parent_ids_field);

    ThrowRequire(!mesh.in_modifiable_state());

    const auto & buckets = mesh.get_buckets(stk::topology::NODE_RANK, field_selector);
    for(auto && b_ptr : buckets)
    {
      for(unsigned i=0; i < b_ptr->size(); ++i)
      {
        const stk::mesh::Entity child_node = (*b_ptr)[i];
        const auto parent_ids = get_edge_node_parent_ids(mesh, parent_ids_field, child_node);
        my_parents_to_child_map[{parent_ids[0], parent_ids[1]}] = mesh.identifier(child_node);
      }
    }
  }

  if (add_higher_order_midside_nodes)
  {
    // Add higher order edge "children"
    for(auto && b_ptr : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
    {
      const stk::topology elem_topology = b_ptr->topology();
      if (elem_topology == stk::topology::TRIANGLE_6_2D || elem_topology == stk::topology::TETRAHEDRON_10)
      {
        for (auto elem : *b_ptr)
        {
          const stk::mesh::Entity * elem_nodes = mesh.begin_nodes(elem);
          for (unsigned edge_i=0; edge_i<elem_topology.num_edges(); ++edge_i) // higher order nodes
          {
            const unsigned * edge_lnn = get_edge_node_ordinals(elem_topology, edge_i);
            const stk::mesh::EntityId parent1 = mesh.identifier(elem_nodes[edge_lnn[0]]);
            const stk::mesh::EntityId parent2 = mesh.identifier(elem_nodes[edge_lnn[1]]);
            const stk::mesh::EntityId midside = mesh.identifier(elem_nodes[edge_lnn[2]]);
            if (my_parents_to_child_map.find({parent1,parent2}) == my_parents_to_child_map.end())
              my_parents_to_child_map[{parent1, parent2}] = midside;
          }
        }
      }
    }
  }
}

} // namespace krino
