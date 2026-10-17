// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MetricEdgeMarking.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Edge.hpp>
#include <Akri_EdgeMarker.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSizeMetric.hpp>
#include <Akri_Refinement.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <array>
#include <set>
#include <vector>

namespace krino {

size_t mark_edges_too_long(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const FieldRef coordsField,
    const MeshSizeMetric & sizeMetric,
    const double maxRelativeEdgeLength,
    EdgeBasedEdgeMarker & edgeMarker)
{
  const unsigned spatialDim = mesh.mesh_meta_data().spatial_dimension();
  const stk::mesh::Selector ownedNonParent = mesh.mesh_meta_data().locally_owned_part() & !refinement.parent_part();
  std::set<Edge> checkedEdges;
  std::vector<Edge> elemEdges;
  size_t numMarked = 0;
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedNonParent))
  {
    for (auto elem : *bucketPtr)
    {
      fill_entity_edges(mesh, elem, elemEdges);
      for (auto && edge : elemEdges)
      {
        if (!checkedEdges.insert(edge).second)
          continue;
        const std::array<stk::mesh::Entity,2> edgeNodes = get_edge_nodes(edge);
        const stk::math::Vector3d edgeVec = get_vector_field(mesh, coordsField, edgeNodes[1], spatialDim) - get_vector_field(mesh, coordsField, edgeNodes[0], spatialDim);
        if (sizeMetric.compute_relative_edge_length(edgeNodes[0], edgeNodes[1], edgeVec) > maxRelativeEdgeLength)
        {
          edgeMarker.mark_edge_for_refinement(edge);
          ++numMarked;
        }
      }
    }
  }
  return numMarked;
}

void update_size_for_new_nodes(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const MeshSizeMetric & sizeMetric)
{
  refinement.for_each_newly_refined_node(
      [&sizeMetric](const stk::mesh::Entity newNode, const stk::mesh::Entity * parentBegin, const stk::mesh::Entity * parentEnd)
      { sizeMetric.update_size_at_new_node(newNode, parentBegin, parentEnd); });
  sizeMetric.communicate_updated_sizes(mesh);
}

static size_t mark_edges_too_long_on_all_procs(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const FieldRef coordsField,
    const MeshSizeMetric & sizeMetric,
    const double maxRelativeEdgeLength,
    EdgeBasedEdgeMarker & edgeMarker)
{
  edgeMarker.clear_marked_edges();
  return stk::get_global_sum(mesh.parallel(), mark_edges_too_long(mesh, refinement, coordsField, sizeMetric, maxRelativeEdgeLength, edgeMarker));
}

size_t refine_edges_too_long(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const FieldRef coordsField,
    const MeshSizeMetric & sizeMetric,
    const double maxRelativeEdgeLength,
    const unsigned maxNumPasses)
{
  EdgeBasedEdgeMarker edgeMarker(mesh, refinement);
  for (unsigned pass=0; pass<maxNumPasses; ++pass)
  {
    const size_t numTooLong = mark_edges_too_long_on_all_procs(mesh, refinement, coordsField, sizeMetric, maxRelativeEdgeLength, edgeMarker);
    if (numTooLong == 0)
      return 0;
    krinolog << "Refinement pass " << pass << ": refining " << numTooLong << " edges too long in the size metric." << stk::diag::dendl;
    refinement.do_refinement(edgeMarker);
    update_size_for_new_nodes(mesh, refinement, sizeMetric);
  }
  const size_t numTooLong = mark_edges_too_long_on_all_procs(mesh, refinement, coordsField, sizeMetric, maxRelativeEdgeLength, edgeMarker);
  edgeMarker.clear_marked_edges();
  if (numTooLong > 0)
    krinolog << "After " << maxNumPasses << " refinement passes, " << numTooLong << " edges are still too long in the size metric." << stk::diag::dendl;
  return numTooLong;
}

} // namespace krino
