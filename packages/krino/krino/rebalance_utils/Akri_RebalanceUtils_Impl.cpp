// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RebalanceUtils_Impl.hpp>

#include <Akri_CDMesh.hpp>
#include <Akri_Element.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_RefinementManager.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <utility>
#include <vector>


namespace krino {
namespace rebalance_utils {
namespace impl {

void set_family_tree_destinations(stk::balance::DecompositionChangeList & decomp_changes, const RefinementManager& refinement, const stk::mesh::BulkData & bulk_data)
{
  // At this point decomp_changes has all adaptivity parent + child elements that are moving
  // and their dest procs match.
  // First iterate the element family trees and set their destinations.
  auto changes = decomp_changes.get_all_partition_changes();
  for (auto && change : changes)
  {
    const auto elem = change.first;
    const auto dest_proc = change.second;

    const auto * fts = bulk_data.begin(elem, stk::topology::CONSTRAINT_RANK);
    const int num_fts = bulk_data.num_connectivity(elem, stk::topology::CONSTRAINT_RANK);
    for (int i = 0; i < num_fts; ++i)
    {
      if (bulk_data.bucket(elem).owned())
      {
        STK_ThrowRequire(!decomp_changes.has_entity(fts[i]) ||
            decomp_changes.get_entity_destination(fts[i]) == dest_proc);
        decomp_changes.set_entity_destination(fts[i], dest_proc);
      }
    }
  }

  std::vector<stk::mesh::Entity> child_sides;
  const auto side_rank = bulk_data.mesh_meta_data().side_rank();
  changes = decomp_changes.get_all_partition_changes();
  for (auto && change : changes)
  {
    const auto elem = change.first;
    const auto dest_proc = change.second;

    if (bulk_data.entity_rank(elem) != stk::topology::ELEMENT_RANK ||
        refinement.is_child(elem) || !refinement.is_parent(elem))
      continue;

    const auto * root_sides = bulk_data.begin(elem, side_rank);
    const int num_root_sides = bulk_data.num_connectivity(elem, side_rank);
    for (int i = 0; i < num_root_sides; ++i)
    {
      if (!bulk_data.bucket(root_sides[i]).owned()) continue;

      fill_all_children(refinement, root_sides[i], child_sides);
      child_sides.push_back(root_sides[i]);

      for (auto && child_side : child_sides)
      {
        STK_ThrowRequire(bulk_data.bucket(child_side).owned());
        const auto * fts = bulk_data.begin(child_side, stk::topology::CONSTRAINT_RANK);
        const int num_fts = bulk_data.num_connectivity(child_side, stk::topology::CONSTRAINT_RANK);
        for (int j = 0; j < num_fts; ++j)
        {
          const auto ft = fts[j];
          if (!decomp_changes.has_entity(ft) ||
              decomp_changes.get_entity_destination(ft) > dest_proc)
          {
            decomp_changes.set_entity_destination(ft, dest_proc);
          }
        }
      }
    }
  }
}

void
update_rebalance_for_adaptivity(stk::balance::DecompositionChangeList & decomp_changes,
    const RefinementManager& refinement,
    const stk::mesh::BulkData & bulk_data)
{
  auto all_changes = decomp_changes.get_all_partition_changes();

  // First pass remove all constrained entities from the list of changes.
  // Second pass will set their destinations according to their dependence
  for(auto && change : all_changes)
  {
    stk::mesh::Entity entity = change.first;
    if(refinement.has_rebalance_constraint(entity))
    {
      decomp_changes.delete_entity(entity);
    }
  }

  all_changes = decomp_changes.get_all_partition_changes();
  stk::mesh::EntityVector adapt_children;
  for(auto && change : all_changes)
  {
    stk::mesh::Entity entity = change.first;
    const auto dest = change.second;

    refinement.fill_dependents(entity, adapt_children);
    for(auto && child : adapt_children)
    {
      decomp_changes.set_entity_destination(child, dest);
    }
  }

  set_family_tree_destinations(decomp_changes, refinement, bulk_data);

  all_changes = decomp_changes.get_all_partition_changes();
  for(auto && change : all_changes)
  {
    stk::mesh::Entity entity = change.first;
    const auto dest = change.second;

    if (bulk_data.entity_rank(entity) == stk::topology::CONSTRAINT_RANK)
    {
      auto side_rank = bulk_data.mesh_meta_data().side_rank();
      const auto * conn_sides = bulk_data.begin(entity, side_rank);
      const int num_sides = bulk_data.num_connectivity(entity, side_rank);
      for (int i = 0; i < num_sides; ++i)
      {
        STK_ThrowRequire(!decomp_changes.has_entity(conn_sides[i]));
        const auto * conn_fts = bulk_data.begin(conn_sides[i], stk::topology::CONSTRAINT_RANK);
        const int num_fts =
            bulk_data.num_connectivity(conn_sides[i], stk::topology::CONSTRAINT_RANK);
        for (int j = 0; j < num_fts; ++j)
        {
          STK_ThrowErrorMsgIf(decomp_changes.get_entity_destination(conn_fts[j]) != dest,
              "Family Tree:\n" << debug_entity(bulk_data, conn_fts[j])
              << "\nHas destination proc = "
              << decomp_changes.get_entity_destination(conn_fts[j])
              << ", which is different than " << dest
              << " the destination for family tree:\n"
              << debug_entity(bulk_data, entity) << "\n\n");
        }
      }
    }
  }
}

void
update_rebalance_for_cdfem(stk::balance::DecompositionChangeList & decomp_changes,
    const stk::mesh::BulkData & bulk_data,
    const CDMesh & cdmesh)
{
  auto all_changes = decomp_changes.get_all_partition_changes();
  std::vector<const SubElement *> subelements;
  const auto & cdfem_parent_part = cdmesh.get_parent_part();
  const auto & cdfem_child_part = cdmesh.get_child_part();
  for(auto && change : all_changes)
  {
    stk::mesh::Entity entity = change.first;
    const auto dest = change.second;

    if(bulk_data.bucket(entity).member(cdfem_parent_part))
    {
      const auto * mesh_elem = cdmesh.find_mesh_element(bulk_data.identifier(entity));
      STK_ThrowRequire(mesh_elem);

      subelements.clear();
      mesh_elem->get_subelements(subelements);

      for(auto && subelem : subelements)
      {
        auto subelem_entity = subelem->entity();
        STK_ThrowRequire(bulk_data.is_valid(subelem_entity));

        decomp_changes.set_entity_destination(subelem_entity, dest);
      }
    }
    else if(bulk_data.bucket(entity).member(cdfem_child_part))
    {
      const auto parent_elem = cdmesh.get_parent_element(entity);
      if(!decomp_changes.has_entity(parent_elem)) decomp_changes.delete_entity(entity);
    }
  }
}

void
accumulate_cdfem_child_weights_to_parents(const stk::mesh::BulkData & bulk_data,
    stk::mesh::Field<double> & element_weights_field,
    const CDMesh & cdmesh)
{
  auto selector =
      stk::mesh::selectField(element_weights_field) & cdmesh.get_child_part() &
      bulk_data.mesh_meta_data().locally_owned_part();
  const auto & buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, selector);
  for(auto && b_ptr : buckets)
  {
    double * weight_data = stk::mesh::field_data(element_weights_field, *b_ptr);
    const int b_size = b_ptr->size();
    for(int i=0; i < b_size; ++i)
    {
      const auto child_elem = (*b_ptr)[i];
      const auto parent = cdmesh.get_parent_element(child_elem);
      double * parent_weight_data = stk::mesh::field_data(element_weights_field, parent);
      STK_ThrowRequire(parent_weight_data);
      *parent_weight_data += weight_data[i];
      weight_data[i] = 0.;
    }
  }
}

void accumulate_adaptivity_child_weights_to_parents(
    const stk::mesh::BulkData & /*bulk_data*/, const RefinementManager& refinement, stk::mesh::Field<double> & element_weights_field)
{
  refinement.update_element_rebalance_weights_incorporating_parallel_owner_constraints(element_weights_field);
}

namespace
{
bool is_owned(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity)
{
  const bool owned_part = mesh.bucket(entity).member(mesh.mesh_meta_data().locally_owned_part());
  const bool owner_rank = mesh.parallel_rank() == mesh.parallel_owner_rank(entity);
  STK_ThrowErrorMsgIf(owned_part != owner_rank,
      "Mismatch between parallel_owner_rank and locally_owned_part membership for:\n"
      << debug_entity(mesh, entity) << "\n");
  return owner_rank;
}
}

bool check_family_tree_element_and_side_ownership(const stk::mesh::BulkData & bulk_data)
{
  const auto & meta = bulk_data.mesh_meta_data();
  bool passed = true;

  const auto & locally_owned_part = meta.locally_owned_part();
  const auto & elem_buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_owned_part);
  for(auto && b_ptr : elem_buckets)
  {
    for(auto && elem : *b_ptr)
    {
      const auto * fts = bulk_data.begin(elem, stk::topology::CONSTRAINT_RANK);
      const int num_fts = bulk_data.num_connectivity(elem, stk::topology::CONSTRAINT_RANK);
      for(int i=0; i < num_fts; ++i)
      {
        if(!is_owned(bulk_data, fts[i]))
        {
          krinolog << "Found non-owned family tree connected to owned element.\n";
          krinolog << "Family tree:\n" << debug_entity(bulk_data, fts[i]) << "\n";
          passed = false;
        }
      }
    }
  }
  const auto side_rank = meta.side_rank();
  const auto & side_buckets = bulk_data.get_buckets(side_rank, locally_owned_part);
  for(auto && b_ptr : side_buckets)
  {
    for(auto && side : *b_ptr)
    {
      const auto * fts = bulk_data.begin(side, stk::topology::CONSTRAINT_RANK);
      const int num_fts = bulk_data.num_connectivity(side, stk::topology::CONSTRAINT_RANK);
      for(int i=0; i < num_fts; ++i)
      {
        if(!is_owned(bulk_data, fts[i]))
        {
          krinolog << "Found non-owned family tree connected to owned side.\n";
          krinolog << "Family tree:\n" << debug_entity(bulk_data, fts[i]) << "\n";
          passed = false;
        }
      }
    }
  }

  const auto & ft_buckets = bulk_data.get_buckets(stk::topology::CONSTRAINT_RANK, locally_owned_part);
  for(auto && b_ptr : ft_buckets)
  {
    for(auto && ft : *b_ptr)
    {
      const auto * elems = bulk_data.begin(ft, stk::topology::ELEMENT_RANK);
      const int num_elems = bulk_data.num_connectivity(ft, stk::topology::ELEMENT_RANK);
      for(int i=0; i < num_elems; ++i)
      {
        if(!is_owned(bulk_data, elems[i]))
        {
          krinolog << "Found non-owned element connected to owned family tree.\n";
          krinolog << "Element ID:\n" << bulk_data.identifier(elems[i]) << "\n";
          krinolog << "Family tree:\n" << debug_entity(bulk_data, ft) << "\n";
          passed = false;
        }
      }
      const auto * sides = bulk_data.begin(ft, side_rank);
      const int num_sides = bulk_data.num_connectivity(ft, side_rank);
      for(int i=0; i < num_sides; ++i)
      {
        if(!is_owned(bulk_data, sides[i]))
        {
          krinolog << "Found non-owned side connected to owned family tree.\n";
          krinolog << "Element ID:\n" << bulk_data.identifier(sides[i]) << "\n";
          krinolog << "Family tree:\n" << debug_entity(bulk_data, ft) << "\n";
          passed = false;
        }
      }
    }
  }

  return passed;
}

}
}
}
