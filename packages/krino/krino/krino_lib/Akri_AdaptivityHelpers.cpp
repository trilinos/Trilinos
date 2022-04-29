// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldState.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <functional>
#include <string>
#include <vector>

namespace krino
{

stk::mesh::Part & get_refinement_active_part(const stk::mesh::MetaData & meta, stk::mesh::EntityRank rank)
{
  const std::string active_part_name = "refine_active_elements_part_"+std::to_string((int)rank);
  stk::mesh::Part* active_part = meta.get_part(active_part_name);
  ThrowRequireMsg(nullptr != active_part, "Active part not found: " << active_part_name);
  return *active_part;
}
stk::mesh::Part & get_refinement_inactive_part(const stk::mesh::MetaData & meta, stk::mesh::EntityRank rank)
{
  const std::string inactive_part_name = "refine_inactive_elements_part_"+std::to_string((int)rank);
  stk::mesh::Part* inactive_part = meta.get_part(inactive_part_name);
  ThrowRequireMsg(nullptr != inactive_part, "Inactive part not found: " << inactive_part_name);
  return *inactive_part;
}

void filter_refinement_marker(const stk::mesh::BulkData & mesh, FieldRef elem_marker, const stk::mesh::Selector & do_not_refine_or_unrefine_selector)
{
  const auto & perceptParentPart = get_refinement_inactive_part(mesh.mesh_meta_data(), stk::topology::ELEMENT_RANK);

  for (auto && bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part()))
  {
    int * markers = field_data<int>(elem_marker, *bucketPtr);
    const int size = bucketPtr->size();
    if (do_not_refine_or_unrefine_selector(*bucketPtr))
    {
      for (int i = 0; i < size; ++i)
        if (markers[i] == Refinement_Marker::REFINE || markers[i] == Refinement_Marker::COARSEN)
          markers[i] = Refinement_Marker::NOTHING;
    }
    else if (bucketPtr->member(perceptParentPart))
    {
      for (int i = 0; i < size; ++i)
        if (markers[i] == Refinement_Marker::REFINE)
          markers[i] = Refinement_Marker::NOTHING;
    }
  }
}

stk::mesh::Selector cdfem_do_not_refine_or_unrefine_selector(const CDFEM_Support & cdfem_support)
{
  const stk::mesh::Selector parent_or_child_selector =
      cdfem_support.get_child_part() | cdfem_support.get_parent_part();
  const stk::mesh::Selector decomposed_blocks_selector =
      krino::Phase_Support::get(cdfem_support.get_mesh_meta()).get_all_decomposed_blocks_selector();
  const stk::mesh::Selector do_not_refine_selector = (!decomposed_blocks_selector) | parent_or_child_selector;
  return do_not_refine_selector;
}


void perform_multilevel_adaptivity(stk::mesh::BulkData & mesh,
    const std::string & marker_field_name,
    const std::function<void(const std::string &, int)> & marker_function,
    const std::function<void(const std::string &, int)> & adapt_function,
    const stk::mesh::Selector & do_not_refine_selector)
{
  Tracespec trace__("perform_multilevel_adaptivity()");

  const auto & aux_meta = AuxMetaData::get(mesh.mesh_meta_data());

  const FieldRef elem_marker = aux_meta.get_field(
      stk::topology::ELEMENT_RANK, marker_field_name, stk::mesh::StateNew);

  const stk::mesh::Selector active_selector = aux_meta.active_part();
  const stk::mesh::Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> children;

  bool done = false;
  int num_refinements = 0;
  while (!done)
  {
    marker_function(marker_field_name, num_refinements);
    filter_refinement_marker(mesh, elem_marker, do_not_refine_selector);

    const auto & local_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, locally_owned_selector);
    unsigned num_marked_refine = 0;
    for (auto && b_ptr : local_buckets)
    {
      int * markers = field_data<int>(elem_marker, *b_ptr);
      for (size_t i = 0; i < b_ptr->size(); ++i)
      {
        if (markers[i] == Refinement_Marker::REFINE) ++num_marked_refine;
      }
    }

    // Only worth cost of adaptive refinement if any elements are marked for refinement
    unsigned global_num_marked_refine = 0;
    stk::all_reduce_sum(mesh.parallel(), &num_marked_refine, &global_num_marked_refine, 1);
    if (global_num_marked_refine > 0)
    {
      const int debug_level = 0;
      adapt_function(marker_field_name, debug_level);

      ++num_refinements;
    }
    else
    {
      done = true;
    }
  }

  // This probably should not be needed.
  CDMesh::fixup_adapted_element_parts(mesh);

  // This probably should not be needed.
  attach_sides_to_elements(mesh);
}
}
