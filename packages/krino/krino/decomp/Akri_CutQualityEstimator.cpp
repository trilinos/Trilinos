// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CutQualityEstimator.hpp>
#include <Akri_EntityOrder.hpp>

#include <algorithm>
#include <set>

#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace krino
{

static stk::math::Vector3d compute_intersection_point_location(
    const int dim,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & intPtNodes,
    const std::vector<double> & intPtWeights)
{
  stk::math::Vector3d snapLocation = stk::math::Vector3d::ZERO;
  auto coordsData = coordsField.data<double>();
  for (size_t i=0; i<intPtNodes.size(); ++i)
  {
    const stk::math::Vector3d nodeLocation(coordsData.entity_values(intPtNodes[i]), dim);
    snapLocation += intPtWeights[i] * nodeLocation;
  }
  return snapLocation;
}

stk::math::Vector3d compute_intersection_point_location(const int dim, const FieldRef coordsField, const IntersectionPoint & intersectionPoint)
{
  return compute_intersection_point_location(dim, coordsField, intersectionPoint.get_nodes(), intersectionPoint.get_weights());
}

void clip_intersection_point_weights(std::vector<double> & intPtWeights, const double tol)
{
  double adjustment = 0.;
  double sumGood = 0.;
  for (auto && intPtWeight : intPtWeights)
  {
    if (intPtWeight < tol)
      adjustment += tol - intPtWeight;
    else
      sumGood += intPtWeight;
  }

  if (adjustment > 0.)
  {
    adjustment = adjustment/sumGood;
    for (auto && intPtWeight : intPtWeights)
    {
      if (intPtWeight < tol)
        intPtWeight = tol;
      else
        intPtWeight -= intPtWeight*adjustment;
    }
  }
}

static bool element_has_all_nodes(const std::vector<stk::mesh::Entity> & elemNodes, const std::vector<stk::mesh::Entity> & nodesToFind)
{
  for (auto && nodeToFind : nodesToFind)
    if (std::find(elemNodes.begin(), elemNodes.end(), nodeToFind) == elemNodes.end())
      return false;
  return true;
}

static bool element_will_go_away_if_node_is_snapped(const std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal)
{
  int numElemNodeOnInterfaceOrIsOriginal = 0;
  for (bool isNodeOnInterfaceOrIsOriginal : isElemNodeOnInterfaceOrIsOriginal)
  {
    if (isNodeOnInterfaceOrIsOriginal)
      if (++numElemNodeOnInterfaceOrIsOriginal > 1)
        return true;
  }
  return false;
}

static double get_node_intersection_point_weight(const IntersectionPoint & intersectionPoint, stk::mesh::Entity node)
{
  const std::vector<stk::mesh::Entity> & nodes = intersectionPoint.get_nodes();
  const auto iter = std::find(nodes.begin(), nodes.end(), node);
  STK_ThrowRequire(iter != nodes.end());
  const auto index = std::distance(nodes.begin(), iter);
  return intersectionPoint.get_weights()[index];
}

static std::vector<stk::mesh::EntityId> get_sorted_node_ids(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & nodes)
{
  std::vector<stk::mesh::EntityId> nodeIds;
  nodeIds.reserve(nodes.size());
  for (auto && node : nodes)  nodeIds.push_back(mesh.identifier(node));
  std::sort(nodeIds.begin(), nodeIds.end());
  return nodeIds;
}

static void sort_intersection_points_for_cutting(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const stk::mesh::Entity node,
    const bool globalIDsAreParallelConsistent,
    std::vector<size_t> & sortedIntersectionPointIndices)
{
  // This sorter is designed to match the priority used by the cutting algorithm with CUT_QUADS_BY_MINIMIZING_ANGLES, especially in 3D
  auto sorter = [&intersectionPoints, &mesh, &coordsField, node, globalIDsAreParallelConsistent](const size_t intPtIndex0, const size_t intPtIndex1)
    {
      const IntersectionPoint & intPt0 = intersectionPoints[intPtIndex0];
      const IntersectionPoint & intPt1 = intersectionPoints[intPtIndex1];
      const size_t numDomains0 = intPt0.get_sorted_domains().size();
      const size_t numDomains1 = intPt1.get_sorted_domains().size();
      if (numDomains0 != numDomains1)
        return numDomains0 > numDomains1;
      // reduce precision to float to handle "ties"
      const float wt0 = get_node_intersection_point_weight(intPt0, node);
      const float wt1 = get_node_intersection_point_weight(intPt1, node);
      if (wt0 != wt1)
        return wt0 < wt1;

      // Then by location, or, where the ids are parallel consistent, by the stencils' sorted ids,
      // larger first.
      if (globalIDsAreParallelConsistent)
      {
        const std::vector<stk::mesh::EntityId> sortedNodes0 = get_sorted_node_ids(mesh, intPt0.get_nodes());
        const std::vector<stk::mesh::EntityId> sortedNodes1 = get_sorted_node_ids(mesh, intPt1.get_nodes());
        return less_by_global_id_because_ids_are_parallel_consistent(sortedNodes1, sortedNodes0);
      }

      const stk::math::Vector3d x0 = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPt0);
      const stk::math::Vector3d x1 = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPt1);
      return less_by_coordinates_because_it_impacts_structure(x0, x1);
    };
  std::sort(sortedIntersectionPointIndices.begin(), sortedIntersectionPointIndices.end(), sorter);
}

static std::set<stk::mesh::Entity> get_intersected_elements(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<size_t> & intersectionPointIndices)
{
  std::vector<stk::mesh::Entity> cutElems;
  std::set<stk::mesh::Entity> intersectedElements;
  for (size_t intPtIndex : intersectionPointIndices)
  {
    stk::mesh::get_entities_through_relations(mesh, intersectionPoints[intPtIndex].get_nodes(), stk::topology::ELEMENT_RANK, cutElems);
    for (auto && cutElem : cutElems)
      if (elementSelector(mesh.bucket(cutElem)))
        intersectedElements.insert(cutElem);
  }
  return intersectedElements;
}

static void assign_is_elem_node_original_node(const std::vector<stk::mesh::Entity> & elemNodes, const stk::mesh::Entity originalNode, std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal)
{
  isElemNodeOnInterfaceOrIsOriginal.clear();
  std::for_each(elemNodes.begin(), elemNodes.end(), [originalNode, &isElemNodeOnInterfaceOrIsOriginal](stk::mesh::Entity n) { isElemNodeOnInterfaceOrIsOriginal.push_back(n == originalNode);});
}

CutQualityEstimator::CutQualityEstimator(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const QualityMetric & qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const bool globalIDsAreParallelConsistent)
  : myMesh(mesh),
    myElementSelector(elementSelector),
    myCoordsField(coordsField),
    myQualityMetric(qualityMetric),
    myMinIntPtWeightForEstimatingCutQuality(minIntPtWeightForEstimatingCutQuality),
    myGlobalIDsAreParallelConsistent(globalIDsAreParallelConsistent)
{
}

double CutQualityEstimator::estimate_quality_of_cutting_element(const std::vector<stk::mesh::Entity> & elemNodes,
    const std::vector<stk::math::Vector3d> & elemNodeCoords,
    const std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<size_t> & intersectionPointIndices) const
{
  if (intersectionPointIndices.empty())
  {
    if (element_will_go_away_if_node_is_snapped(isElemNodeOnInterfaceOrIsOriginal))
      return std::max(0., myQualityMetric.get_element_quality_metric(myMesh.mesh_meta_data().spatial_dimension(), elemNodeCoords));
    return myQualityMetric.get_best_value_for_metric();
  }

  // apply the front intersection point to element and recursively call with remaining intersection points
  const IntersectionPoint & intPtToApply = intersectionPoints[*intersectionPointIndices.begin()];
  const std::vector<size_t> remainingIntersectionPointIndices(intersectionPointIndices.begin()+1, intersectionPointIndices.end());

  const auto & intPtNodes = intPtToApply.get_nodes();
  std::vector<double> intPtWts = intPtToApply.get_weights();
  clip_intersection_point_weights(intPtWts, myMinIntPtWeightForEstimatingCutQuality);

  const stk::math::Vector3d intPtLocation = compute_intersection_point_location(myMesh.mesh_meta_data().spatial_dimension(), myCoordsField, intPtNodes, intPtWts);

  double qualityAfterCut = myQualityMetric.get_best_value_for_metric();
  if (element_has_all_nodes(elemNodes, intPtNodes))
  {
    std::vector<stk::mesh::Entity> cutElemNodes;
    std::vector<stk::math::Vector3d> cutElemNodeCoords;
    std::vector<bool> isCutElemNodeOnInterfaceOrIsOriginal;

    for (auto intPtNode : intPtNodes)
    {
      cutElemNodes.clear();
      cutElemNodeCoords.clear();
      isCutElemNodeOnInterfaceOrIsOriginal.clear();
      for (size_t nodeIndex=0; nodeIndex<elemNodes.size(); ++nodeIndex)
      {
        if (elemNodes[nodeIndex] == intPtNode)
        {
          cutElemNodes.push_back(stk::mesh::Entity());
          cutElemNodeCoords.push_back(intPtLocation);
          isCutElemNodeOnInterfaceOrIsOriginal.push_back(true);
        }
        else
        {
          cutElemNodes.push_back(elemNodes[nodeIndex]);
          cutElemNodeCoords.push_back(elemNodeCoords[nodeIndex]);
          isCutElemNodeOnInterfaceOrIsOriginal.push_back(isElemNodeOnInterfaceOrIsOriginal[nodeIndex]);
        }
      }

      const double elemQualityAfterCuts = estimate_quality_of_cutting_element(cutElemNodes, cutElemNodeCoords, isCutElemNodeOnInterfaceOrIsOriginal, intersectionPoints, remainingIntersectionPointIndices);
      if (myQualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
        qualityAfterCut = elemQualityAfterCuts;
    }
  }
  else
  {
    qualityAfterCut = estimate_quality_of_cutting_element(elemNodes, elemNodeCoords, isElemNodeOnInterfaceOrIsOriginal, intersectionPoints, remainingIntersectionPointIndices);
  }

  return qualityAfterCut;
}

std::vector<size_t> CutQualityEstimator::sorted_intersection_point_indices_for_node_for_domains(
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<std::pair<size_t,bool>> & nodeIntersectionPointIndicesAndWhichSnapsAllowed,
    const stk::mesh::Entity node,
    const std::vector<int> & domains) const
{
  std::vector<size_t> sortedIntersectionPointIndices;
  for (auto && intPtIndexAndIsSnapAllowed : nodeIntersectionPointIndicesAndWhichSnapsAllowed)
  {
    const size_t intPtIndex = intPtIndexAndIsSnapAllowed.first;
    if (first_sorted_vector_of_domains_contains_all_domains_in_second_vector(domains, intersectionPoints[intPtIndex].get_sorted_domains()))
      sortedIntersectionPointIndices.push_back(intPtIndex);
  }

  sort_intersection_points_for_cutting(myMesh, myCoordsField, intersectionPoints, node, myGlobalIDsAreParallelConsistent, sortedIntersectionPointIndices);
  return sortedIntersectionPointIndices;
}

double CutQualityEstimator::estimate_quality_of_cutting_node_elements(const stk::mesh::Entity node,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<size_t> & sortedIntersectionPointIndices) const
{
  const int dim = myMesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> elemNodes;
  std::vector<stk::math::Vector3d> elemNodeCoords;
  std::vector<bool> isElemNodeOnInterfaceOrIsOriginal;

  const std::set<stk::mesh::Entity> intersectedElements = get_intersected_elements(myMesh, myElementSelector, intersectionPoints, sortedIntersectionPointIndices);

  double qualityAfterCut = myQualityMetric.get_best_value_for_metric();
  for (auto && elem : intersectedElements)
  {
    elemNodes.assign(myMesh.begin_nodes(elem), myMesh.end_nodes(elem));
    fill_node_locations(dim, myCoordsField, elemNodes, elemNodeCoords);
    assign_is_elem_node_original_node(elemNodes, node, isElemNodeOnInterfaceOrIsOriginal);
    const double elemQualityAfterCuts = estimate_quality_of_cutting_element(elemNodes, elemNodeCoords, isElemNodeOnInterfaceOrIsOriginal, intersectionPoints, sortedIntersectionPointIndices);

    if (myQualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
      qualityAfterCut = elemQualityAfterCuts;
  }
  return qualityAfterCut;
}

std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> CutQualityEstimator::estimate_quality_per_node_per_domain(
    const std::vector<IntersectionPoint> & intersectionPoints,
    const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed) const
{
  std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> domainsToNodesToQuality;
  for (auto entry : nodeToIntPtIndicesAndWhichSnapsAllowed)
  {
    stk::mesh::Entity node = entry.first;
    const auto nodeIntersectionPointIndicesAndWhichSnapsAllowed = entry.second;

    std::set<std::vector<int>> nodeIntPtDomains;
    for (auto && intPtIndexAndIsSnapAllowed : nodeIntersectionPointIndicesAndWhichSnapsAllowed)
      nodeIntPtDomains.insert(intersectionPoints[intPtIndexAndIsSnapAllowed.first].get_sorted_domains());

    for (auto && intPtDomains : nodeIntPtDomains)
    {
      const std::vector<size_t> sortedIntersectionPointIndices = sorted_intersection_point_indices_for_node_for_domains(intersectionPoints, nodeIntersectionPointIndicesAndWhichSnapsAllowed, node, intPtDomains);
      domainsToNodesToQuality[intPtDomains][myMesh.identifier(node)] = estimate_quality_of_cutting_node_elements(node, intersectionPoints, sortedIntersectionPointIndices);
    }
  }

  return domainsToNodesToQuality;
}

}
