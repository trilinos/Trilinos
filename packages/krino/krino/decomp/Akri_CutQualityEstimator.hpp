// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_AKRI_CUTQUALITYESTIMATOR_HPP_
#define KRINO_KRINO_DECOMP_AKRI_CUTQUALITYESTIMATOR_HPP_
#include <map>
#include <utility>
#include <vector>

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino
{
class IntersectionPoint;
class QualityMetric;

// For each node, the indices of the intersection points that use it, and whether the node may be
// snapped to each.
typedef std::map<stk::mesh::Entity, std::vector<std::pair<size_t,bool>>> mapFromEntityToIntPtIndexAndSnapAllowed;

stk::math::Vector3d compute_intersection_point_location(const int dim, const FieldRef coordsField, const IntersectionPoint & intersectionPoint);

// Raises the weights below tol to tol, scaling the others down so that the sum is unchanged, so
// that an estimated cut is not closer to a node than tol.
void clip_intersection_point_weights(std::vector<double> & intPtWeights, const double tol);

// Estimates the quality of the elements made by cutting the mesh at intersection points, without
// cutting it: what snapping a node to an intersection point is compared against.
class CutQualityEstimator
{
public:
  CutQualityEstimator(const stk::mesh::BulkData & mesh,
      const stk::mesh::Selector & elementSelector,
      const FieldRef coordsField,
      const QualityMetric & qualityMetric,
      const double minIntPtWeightForEstimatingCutQuality,
      const bool globalIDsAreParallelConsistent);

  // For each set of domains of the intersection points of each node, the worst estimated quality of
  // cutting the elements those points intersect.
  std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> estimate_quality_per_node_per_domain(
      const std::vector<IntersectionPoint> & intersectionPoints,
      const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed) const;

  // The worst estimated quality of cutting the elements intersected by the points, in order, with
  // the node the original node of each.
  double estimate_quality_of_cutting_node_elements(const stk::mesh::Entity node,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<size_t> & sortedIntersectionPointIndices) const;

  // The worst estimated quality of cutting one element at the points, in order: each point whose
  // nodes are all in the element splits it at each of those nodes, and the pieces are cut at the
  // rest.  A piece with no points left counts, at its own quality, only if two or more of its nodes
  // are on the interface or are the original node: a piece that snapping the original node would
  // collapse, as the sliver between it and a nearby point.  Other pieces count as the best value.
  double estimate_quality_of_cutting_element(const std::vector<stk::mesh::Entity> & elemNodes,
      const std::vector<stk::math::Vector3d> & elemNodeCoords,
      const std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<size_t> & intersectionPointIndices) const;

  // The indices of the node's points whose domains are all in domains, in the order of cutting:
  // more domains first, then the smaller weight at the node.
  std::vector<size_t> sorted_intersection_point_indices_for_node_for_domains(
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<std::pair<size_t,bool>> & nodeIntersectionPointIndicesAndWhichSnapsAllowed,
      const stk::mesh::Entity node,
      const std::vector<int> & domains) const;

private:
  const stk::mesh::BulkData & myMesh;
  stk::mesh::Selector myElementSelector;
  FieldRef myCoordsField;
  const QualityMetric & myQualityMetric;
  double myMinIntPtWeightForEstimatingCutQuality;
  bool myGlobalIDsAreParallelConsistent;
};

}

#endif /* KRINO_KRINO_DECOMP_AKRI_CUTQUALITYESTIMATOR_HPP_ */
