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
#include <Akri_Phase_Support.hpp>
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

#include "Akri_RefinementInterface.hpp"
namespace krino
{

void filter_refinement_marker(const RefinementInterface & refinement, const stk::mesh::BulkData & mesh, FieldRef elem_marker, const stk::mesh::Selector & do_not_refine_or_unrefine_selector)
{
  const auto & parentPart = refinement.parent_part();

  for (auto && bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part()))
  {
    int * markers = field_data<int>(elem_marker, *bucketPtr);
    const int size = bucketPtr->size();
    if (do_not_refine_or_unrefine_selector(*bucketPtr))
    {
      for (int i = 0; i < size; ++i)
        if (markers[i] == static_cast<int>(Refinement_Marker::REFINE) || markers[i] == static_cast<int>(Refinement_Marker::COARSEN))
          markers[i] = static_cast<int>(Refinement_Marker::NOTHING);
    }
    else if (bucketPtr->member(parentPart))
    {
      for (int i = 0; i < size; ++i)
        if (markers[i] == static_cast<int>(Refinement_Marker::REFINE))
          markers[i] = static_cast<int>(Refinement_Marker::NOTHING);
    }
  }
}


void perform_multilevel_adaptivity(RefinementInterface & refinement,
    stk::mesh::BulkData & mesh,
    const std::function<void(int)> & marker_function,
    const stk::mesh::Selector & do_not_refine_selector)
{
  Tracespec trace__("perform_multilevel_adaptivity()");

  const auto & aux_meta = AuxMetaData::get(mesh.mesh_meta_data());

  const FieldRef elem_marker = refinement.get_marker_field_and_sync_to_host();

  const stk::mesh::Selector active_selector = aux_meta.active_part();
  const stk::mesh::Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> children;

  bool done = false;
  int num_refinements = 0;
  while (!done)
  {
    marker_function(num_refinements);
    filter_refinement_marker(refinement, mesh, elem_marker, do_not_refine_selector);

    const auto & local_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, locally_owned_selector);
    unsigned num_marked_refine = 0;
    for (auto && b_ptr : local_buckets)
    {
      int * markers = field_data<int>(elem_marker, *b_ptr);
      for (size_t i = 0; i < b_ptr->size(); ++i)
      {
        if (markers[i] == static_cast<int>(Refinement_Marker::REFINE)) ++num_marked_refine;
      }
    }

    // Only worth cost of adaptive refinement if any elements are marked for refinement
    unsigned global_num_marked_refine = 0;
    stk::all_reduce_sum(mesh.parallel(), &num_marked_refine, &global_num_marked_refine, 1);
    if (global_num_marked_refine > 0)
    {
      const int debug_level = 0;
      refinement.do_refinement(debug_level);

      ++num_refinements;
    }
    else
    {
      krinolog << "Skipping/Terminating refinement because no elements are marked for refinement.\n";
      done = true;
    }
  }

  if (refinement.require_post_refinement_fixups())
  {
    // This probably should not be needed.
    CDMesh::fixup_adapted_element_parts(mesh);

    // This probably should not be needed.
    attach_sides_to_elements(mesh);
  }
}
}
