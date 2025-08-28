// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_DistanceSweeper.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino{

namespace {

bool node_has_other_sign(stk::mesh::Entity node, const FieldRef distance_field, const double signed_narrow_band)
{
  double * distance = field_data<double>(distance_field, node);
  STK_ThrowAssert(nullptr != distance);
  return (signed_narrow_band < 0.) ? (*distance > 0.) : (*distance < 0.);
}

bool node_is_inside_narrow_band(stk::mesh::Entity node, const FieldRef distance_field, const double signed_narrow_band)
{
  double * distance = field_data<double>(distance_field, node);
  STK_ThrowAssert(nullptr != distance);
  return (signed_narrow_band < 0.) ? (*distance > signed_narrow_band) : (*distance < signed_narrow_band);
}

void get_nodes_ready_to_update(const AuxMetaData & aux_meta, const stk::mesh::BulkData& mesh, const FieldRef distance_field, const double signed_narrow_band, std::set<stk::mesh::Entity> & nodes_ready_to_update)
{
  stk::mesh::Selector field_not_ghost = aux_meta.active_not_ghost_selector() & stk::mesh::selectField(distance_field);
  const stk::mesh::BucketVector & buckets = mesh.get_buckets( stk::topology::ELEMENT_RANK, field_not_ghost );
  std::vector<stk::mesh::Entity> elem_nodes_to_update;

  for ( auto && bucket : buckets )
  {
    stk::topology elem_topology = bucket->topology();
    const unsigned num_nodes = elem_topology.num_nodes();

    for ( auto && elem : *bucket )
    {
      elem_nodes_to_update.clear();
      const stk::mesh::Entity* nodes = mesh.begin_nodes(elem);
      bool have_node_with_other_sign = false;
      for ( unsigned i = 0; i < num_nodes; ++i )
      {
        stk::mesh::Entity node = nodes[i];
        if (has_field_data(distance_field, node))
        {
          if (node_has_other_sign(node, distance_field, signed_narrow_band))
          {
            have_node_with_other_sign = true;
          }
          else if (!node_is_inside_narrow_band(node, distance_field, signed_narrow_band))
          {
            elem_nodes_to_update.push_back(node);
          }
        }
      }
      if (have_node_with_other_sign)
      {
        nodes_ready_to_update.insert(elem_nodes_to_update.begin(), elem_nodes_to_update.end());
      }
    }
  }
}

void get_neighbor_nodes_ready_to_update(const AuxMetaData & /*aux_meta*/,
    const stk::mesh::BulkData& mesh,
    const FieldRef distance_field,
    const double signed_narrow_band,
    const bool check_if_nbr_is_outside_band,
    stk::mesh::Entity node,
    std::set<stk::mesh::Entity> & nodes_ready_to_update)
{
  const stk::mesh::Entity* node_elems = mesh.begin_elements(node);
  const unsigned num_node_elements = mesh.num_elements(node);

  for (unsigned ielem = 0; ielem < num_node_elements; ++ ielem)
  {
    stk::mesh::Entity elem = node_elems[ielem];
    const unsigned num_nodes = mesh.num_nodes(elem);
    const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);
    for ( unsigned inode = 0; inode < num_nodes; ++inode )
    {
      stk::mesh::Entity nbr_node = elem_nodes[inode];
      if (nbr_node != node &&
          has_field_data(distance_field, nbr_node) &&
          !node_has_other_sign(nbr_node, distance_field, signed_narrow_band) &&
          (!check_if_nbr_is_outside_band || !node_is_inside_narrow_band(nbr_node, distance_field, signed_narrow_band)))
      {
        nodes_ready_to_update.insert(nbr_node);
      }
    }
  }
}

void get_shared_nodes_ready_to_update(const AuxMetaData & aux_meta, const stk::mesh::BulkData& mesh, const FieldRef distance_field, const double signed_narrow_band, std::set<stk::mesh::Entity> & nodes_ready_to_update)
{
  stk::mesh::Selector field_globally_shared = mesh.mesh_meta_data().globally_shared_part() & stk::mesh::selectField(distance_field);
  std::vector< stk::mesh::Entity> shared_nodes;
  stk::mesh::get_selected_entities( field_globally_shared, mesh.buckets( stk::topology::NODE_RANK ), shared_nodes );
  const bool check_if_nbr_is_outside_band = true;

  for (auto && node : shared_nodes)
  {
    if (node_has_other_sign(node, distance_field, signed_narrow_band))
    {
      get_neighbor_nodes_ready_to_update(aux_meta, mesh, distance_field, signed_narrow_band, check_if_nbr_is_outside_band, node, nodes_ready_to_update);
    }
  }
}

} // anonymous

void DistanceSweeper::fix_sign_by_sweeping(const stk::mesh::BulkData& mesh, const FieldRef distance_field, const double signed_narrow_band)
{
  const AuxMetaData & aux_meta = AuxMetaData::get(mesh.mesh_meta_data());
  if (signed_narrow_band == 0.0) return;
  const int sign_to_fix = (signed_narrow_band < 0.) ? -1 : 1;

  std::set<stk::mesh::Entity> nodes_ready_to_update;

  get_nodes_ready_to_update(aux_meta, mesh, distance_field, signed_narrow_band, nodes_ready_to_update);

  bool done = false;
  while (!done)
  {
    size_t num_nodes_updated = 0;

    const bool check_if_nbr_is_outside_band = false;
    while (!nodes_ready_to_update.empty())
    {
      ++num_nodes_updated;
      stk::mesh::Entity node = *nodes_ready_to_update.begin();
      nodes_ready_to_update.erase(nodes_ready_to_update.begin());

      double * distance = field_data<double>(distance_field, node);
      STK_ThrowAssert(nullptr != distance);
      *distance = -signed_narrow_band;

      get_neighbor_nodes_ready_to_update(aux_meta, mesh, distance_field, signed_narrow_band, check_if_nbr_is_outside_band, node, nodes_ready_to_update);
    }

    done = true;
    if (mesh.parallel_size() > 1)
    {
      const size_t local_num_nodes_updated = num_nodes_updated;
      stk::all_reduce_sum(mesh.parallel(), &local_num_nodes_updated, &num_nodes_updated, 1);

      if (num_nodes_updated > 0)
      {
        done = false;

        std::vector<const stk::mesh::FieldBase *> parallel_fields(1, &distance_field.field());
        if (sign_to_fix > 0) stk::mesh::parallel_min(mesh, parallel_fields);
        else  stk::mesh::parallel_max(mesh, parallel_fields);

        get_shared_nodes_ready_to_update(aux_meta, mesh, distance_field, signed_narrow_band, nodes_ready_to_update);
      }
    }
  }
}

//----------------------------------------------------------------
} // namespace krino
