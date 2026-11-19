// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SmoothMeshRespectingMeshBasedGeometry.hpp>

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_MeshTopologySharpFeatures.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMesh.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

stk::mesh::Selector select_nodes_on_any_boundary(const stk::mesh::MetaData & meta)
{
  stk::mesh::PartVector boundaryParts;
  for (auto * part : meta.get_parts())
  {
    const stk::mesh::EntityRank rank = part->primary_entity_rank();
    if ((rank == stk::topology::FACE_RANK || rank == stk::topology::EDGE_RANK || rank == stk::topology::NODE_RANK) &&
        !stk::mesh::is_topology_root_part(*part))
      boundaryParts.push_back(part);
  }
  return stk::mesh::selectUnion(boundaryParts);
}

using SavedNodeLocations = std::vector<std::pair<stk::mesh::Entity, stk::math::Vector3d>>;

static SavedNodeLocations get_node_locations(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector)
{
  SavedNodeLocations nodeLocations;
  for (auto * bucket : mesh.get_buckets(stk::topology::NODE_RANK, elemSelector))
    for (auto node : *bucket)
      nodeLocations.emplace_back(node, get_vector_field(mesh, coordsField, node, coordsField.dim()));
  return nodeLocations;
}

static void set_node_locations(const CoordinatesFieldRef coordsField, const SavedNodeLocations & nodeLocations)
{
  auto coordsData = coordsField.data<double, stk::mesh::ReadWrite>();
  for (auto && [node, location] : nodeLocations)
  {
    auto nodeCoords = coordsData.entity_values(node);
    for (unsigned j = 0; j < coordsField.dim(); ++j)
      nodeCoords(stk::mesh::ComponentIdx(j)) = location[j];
  }
}

// The average length of the edges of the selected elements (each edge counted once per element), over
// all processors.
static double average_edge_length(const stk::mesh::BulkData & mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Selector & elemSelector)
{
  double sumLength = 0.;
  double numEdges = 0.;
  std::array<unsigned, 3> edgeOrdinals;
  for (auto * bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, elemSelector & mesh.mesh_meta_data().locally_owned_part()))
  {
    const stk::topology topology = bucket->topology();
    for (auto elem : *bucket)
    {
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      for (unsigned iEdge = 0; iEdge < topology.num_edges(); ++iEdge)
      {
        topology.edge_node_ordinals(iEdge, edgeOrdinals.data());
        sumLength += (get_vector_field(mesh, coordsField, elemNodes[edgeOrdinals[0]], coordsField.dim()) -
            get_vector_field(mesh, coordsField, elemNodes[edgeOrdinals[1]], coordsField.dim())).length();
        numEdges += 1.;
      }
    }
  }
  const std::array<double, 2> local{{sumLength, numEdges}};
  std::array<double, 2> global;
  stk::all_reduce_sum(mesh.parallel(), local.data(), global.data(), 2);
  return (global[1] > 0.) ? global[0] / global[1] : 1.;
}

// Smooths the selected elements by optimizing the mean-ratio norm of their node locations with the
// search directions filtered, restoring the node locations if the worst mean ratio drops.  The
// optimizer's default step and gradient tolerances are taken relative to the average edge length,
// so that the smoothing does not depend on the units of the coordinates.
static bool smooth_mesh_with_gradient_filter(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const GradientFilter & gradientFilter)
{
  const MeanRatioQualityMetric meanRatio;
  const double worstMeanRatioBefore = compute_mesh_quality(mesh, elemSelector, coordsField, meanRatio);
  const SavedNodeLocations nodeLocationsBefore = get_node_locations(mesh, coordsField, elemSelector);

  const double edgeLength = average_edge_length(mesh, coordsField, elemSelector);
  const OptimizerFunction optimizer = [edgeLength](const ObjectiveInterface & objFn, DistributedVector & nodeLocations)
  {
    const double initialObjective = objFn.compute_value(nodeLocations);
    lbfgs(objFn, nodeLocations, OptimizationDefaults::stepTol * edgeLength, OptimizationDefaults::gradTol / edgeLength, initialObjective);
  };
  const MeanRatioNormObjective smoothingObj;
  smooth_mesh_by_optimization(mesh, coordsField, elemSelector, optimizer, smoothingObj, gradientFilter, SolutionProjector());

  const double worstMeanRatioAfter = compute_mesh_quality(mesh, elemSelector, coordsField, meanRatio);
  if (worstMeanRatioAfter < worstMeanRatioBefore)
  {
    set_node_locations(coordsField, nodeLocationsBefore);
    krinolog << "Smoothing would lower the worst mean ratio from " << worstMeanRatioBefore << " to "
             << worstMeanRatioAfter << ", so the mesh is left as it was." << stk::diag::dendl;
    return false;
  }
  krinolog << "Smoothing changed the worst mean ratio from " << worstMeanRatioBefore << " to "
           << worstMeanRatioAfter << "." << stk::diag::dendl;
  return true;
}

bool smooth_mesh_keeping_boundary_nodes_fixed(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector)
{
  const stk::mesh::Selector boundaryNodeSelector = select_nodes_on_any_boundary(mesh.mesh_meta_data());
  const GradientFilter keepBoundaryNodesFixed =
      [&mesh, &boundaryNodeSelector](const stk::mesh::Entity node, const stk::math::Vector3d &, stk::math::Vector3d & dir)
  {
    if (boundaryNodeSelector(mesh.bucket(node))) dir = stk::math::Vector3d::ZERO;
  };
  return smooth_mesh_with_gradient_filter(mesh, coordsField, elemSelector, keepBoundaryNodesFixed);
}

// The direction restricted as the node's smoothing direction says: within the plane (less its normal
// component), along the tangent, or zero.
static void restrict_to_smoothing_direction(const SmoothingDirection & smoothingDirection, stk::math::Vector3d & dir)
{
  switch (smoothingDirection.type)
  {
    case SmoothingDirection::SMOOTH_IN_PLANE:
      dir -= Dot(dir, smoothingDirection.direction) * smoothingDirection.direction;
      return;
    case SmoothingDirection::SMOOTH_ALONG_TANGENT:
      dir = Dot(dir, smoothingDirection.direction) * smoothingDirection.direction;
      return;
    case SmoothingDirection::FIXED:
      dir = stk::math::Vector3d::ZERO;
      return;
  }
}

bool smooth_mesh_respecting_mesh_based_geometry(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const std::optional<double> cosFeatureAngle)
{
  STK_ThrowRequireMsg(mesh.parallel_size() == 1 || mesh.is_automatic_aura_on(),
      "Smoothing respecting the mesh-based geometry on more than one processor needs the aura, for the sides of the shared nodes.");
  const std::vector<SmoothingDirection> smoothingDirections =
      compute_smoothing_directions_from_mesh_topology(mesh, coordsField, elemSelector, cosFeatureAngle);
  const GradientFilter respectMeshBasedGeometry =
      [&mesh, &smoothingDirections](const stk::mesh::Entity node, const stk::math::Vector3d &, stk::math::Vector3d & dir)
  {
    const stk::mesh::EntityId nodeId = mesh.identifier(node);
    const auto iter = std::lower_bound(smoothingDirections.begin(), smoothingDirections.end(), nodeId,
        [&mesh](const SmoothingDirection & smoothingDirection, const stk::mesh::EntityId id)
        { return mesh.identifier(smoothingDirection.node) < id; });
    if (iter != smoothingDirections.end() && iter->node == node)
      restrict_to_smoothing_direction(*iter, dir);
  };
  return smooth_mesh_with_gradient_filter(mesh, coordsField, elemSelector, respectMeshBasedGeometry);
}

}
