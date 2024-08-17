// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AuxMetaData.hpp>
#include <Akri_Snap.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SharpFeature.hpp>
#include <Akri_SnapIndependentSetFinder.hpp>
#include <Akri_SnapInfo.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <memory>

namespace krino
{

static stk::math::Vector3d compute_intersection_point_location(
    const int dim,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & intPtNodes,
    const std::vector<double> & intPtWeights)
{
  stk::math::Vector3d snapLocation = stk::math::Vector3d::ZERO;
  for (size_t i=0; i<intPtNodes.size(); ++i)
  {
    const stk::math::Vector3d nodeLocation(field_data<double>(coordsField, intPtNodes[i]), dim);
    snapLocation += intPtWeights[i] * nodeLocation;
  }
  return snapLocation;
}

stk::math::Vector3d compute_intersection_point_location(const int dim, const FieldRef coordsField, const IntersectionPoint & intersectionPoint)
{
  return compute_intersection_point_location(dim, coordsField, intersectionPoint.get_nodes(), intersectionPoint.get_weights());
}

static void fill_global_ids_of_elements_using_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    stk::mesh::Entity node,
    std::vector<size_t> & globalIdsOfSnapNodeElems)
{
  globalIdsOfSnapNodeElems.clear();
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    if (elementSelector(mesh.bucket(elem)))
      globalIdsOfSnapNodeElems.push_back(mesh.identifier(elem));
}

double compute_quality_if_node_is_snapped_terminating_early_if_below_threshold(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    stk::mesh::Entity node,
    const stk::math::Vector3d & snapLocation,
    const QualityMetric &qualityMetric,
    const double qualityThreshold)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  double qualityAfterSnap = qualityMetric.get_best_value_for_metric();
  std::vector<stk::math::Vector3d> nodeLocations;

  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (elementSelector(mesh.bucket(elem)))
    {
      nodeLocations.clear();
      for (auto elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
      {
        if (elemNode == node)
          nodeLocations.push_back(snapLocation);
        else
          nodeLocations.emplace_back(field_data<double>(coordsField, elemNode), dim);
      }

      const double elemQualityAfterSnap = qualityMetric.get_element_quality_metric(nodeLocations);

      if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterSnap, elemQualityAfterSnap))
      {
        qualityAfterSnap = elemQualityAfterSnap;
        if (qualityMetric.is_first_quality_metric_better_than_second(qualityThreshold, qualityAfterSnap))
          return qualityAfterSnap;
      }
    }
  }
  return qualityAfterSnap;
}

static bool element_has_all_nodes(const std::vector<stk::mesh::Entity> & elemNodes, const std::vector<stk::mesh::Entity> & nodesToFind)
{
  for (auto && nodeToFind : nodesToFind)
    if (std::find(elemNodes.begin(), elemNodes.end(), nodeToFind) == elemNodes.end())
      return false;
  return true;
}

static void clip_intersection_point_weights(std::vector<double> & intPtWeights, const double tol)
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

static double estimate_quality_of_cutting_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elemNodes,
    const std::vector<stk::math::Vector3d> & elemNodeCoords,
    const std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal,
    const std::vector<size_t> intersectionPointIndices,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric &qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality)
{
  if (intersectionPointIndices.empty())
  {
    if (element_will_go_away_if_node_is_snapped(isElemNodeOnInterfaceOrIsOriginal))
      return std::max(0., qualityMetric.get_element_quality_metric(elemNodeCoords));
    return qualityMetric.get_best_value_for_metric();
  }

  // apply the front intersection point to element and recursively call with remaining intersection points
  const IntersectionPoint & intPtToApply = intersectionPoints[*intersectionPointIndices.begin()];
  const std::vector<size_t> remainingIntersectionPointIndices(intersectionPointIndices.begin()+1, intersectionPointIndices.end());

  const auto & intPtNodes = intPtToApply.get_nodes();
  std::vector<double> intPtWts = intPtToApply.get_weights();
  clip_intersection_point_weights(intPtWts, minIntPtWeightForEstimatingCutQuality);

  const stk::math::Vector3d intPtLocation = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPtNodes, intPtWts);

  double qualityAfterCut = qualityMetric.get_best_value_for_metric();
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

      const double elemQualityAfterCuts = estimate_quality_of_cutting_intersection_points(mesh, coordsField, cutElemNodes, cutElemNodeCoords, isCutElemNodeOnInterfaceOrIsOriginal, remainingIntersectionPointIndices, intersectionPoints, qualityMetric, minIntPtWeightForEstimatingCutQuality);
      if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
        qualityAfterCut = elemQualityAfterCuts;
    }
  }
  else
  {
    qualityAfterCut = estimate_quality_of_cutting_intersection_points(mesh, coordsField, elemNodes, elemNodeCoords, isElemNodeOnInterfaceOrIsOriginal, remainingIntersectionPointIndices, intersectionPoints, qualityMetric, minIntPtWeightForEstimatingCutQuality);
  }

  return qualityAfterCut;
}

static double get_node_intersection_point_weight(const IntersectionPoint & intersectionPoint, stk::mesh::Entity node)
{
  const std::vector<stk::mesh::Entity> & nodes = intersectionPoint.get_nodes();
  const auto iter = std::find(nodes.begin(), nodes.end(), node);
  STK_ThrowRequire(iter != nodes.end());
  const auto index = std::distance(nodes.begin(), iter);
  return intersectionPoint.get_weights()[index];
}

std::vector<stk::mesh::EntityId> get_sorted_node_ids(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & nodes)
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

      if (globalIDsAreParallelConsistent)
      {
        const std::vector<stk::mesh::EntityId> sortedNodes0 = get_sorted_node_ids(mesh, intPt0.get_nodes());
        const std::vector<stk::mesh::EntityId> sortedNodes1 = get_sorted_node_ids(mesh, intPt1.get_nodes());
        return sortedNodes1 < sortedNodes0;
      }

      const stk::math::Vector3d x0 = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPt0);
      const stk::math::Vector3d x1 = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPt1);
      return is_less_than_in_x_then_y_then_z(x0, x1);
    };
  std::sort(sortedIntersectionPointIndices.begin(), sortedIntersectionPointIndices.end(), sorter);
}

static void fill_sorted_intersection_point_indices_for_node_for_domains(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<std::pair<size_t,bool>> & nodeIntersectionPointIndicesAndWhichSnapsAllowed,
    const stk::mesh::Entity node,
    const std::vector<int> & domains,
    const bool globalIDsAreParallelConsistent,
    std::vector<size_t> & sortedIntersectionPointIndices)
{
  sortedIntersectionPointIndices.clear();
  for (auto && intPtIndexAndIsSnapAllowed : nodeIntersectionPointIndicesAndWhichSnapsAllowed)
  {
    const size_t intPtIndex = intPtIndexAndIsSnapAllowed.first;
    if (first_sorted_vector_of_domains_contains_all_domains_in_second_vector(domains, intersectionPoints[intPtIndex].get_sorted_domains()))
      sortedIntersectionPointIndices.push_back(intPtIndex);
  }

  sort_intersection_points_for_cutting(mesh, coordsField, intersectionPoints, node, globalIDsAreParallelConsistent, sortedIntersectionPointIndices);
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

static void filter_which_intersection_point_nodes_are_compatible_for_snapping_based_on_sharp_features(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const std::vector<stk::mesh::Entity> & intPtNodes,
    std::vector<bool> & whichSnapsAreAllowed)
{
  if (nullptr != sharpFeatureInfo)
  {
    STK_ThrowAssert(intPtNodes.size() == whichSnapsAreAllowed.size());
    for (size_t iNode=0; iNode<intPtNodes.size(); ++iNode)
    {
      if (whichSnapsAreAllowed[iNode])
      {
        whichSnapsAreAllowed[iNode] = is_intersection_point_node_compatible_for_snapping_based_on_sharp_features(*sharpFeatureInfo, intPtNodes[iNode], intPtNodes);
        if (false == whichSnapsAreAllowed[iNode])
        {
          krinolog << "Blocked snap of node " << mesh.identifier(intPtNodes[iNode]) << " to int pt with nodes ";
          for (auto && intPtNode : intPtNodes)
            krinolog << mesh.identifier(intPtNode) << " ";
          krinolog << stk::diag::dendl;
        }
      }
    }
  }
}

static std::vector<bool> which_intersection_point_nodes_are_allowed_for_snapping(const stk::mesh::BulkData & mesh,
    const AuxMetaData & auxMeta,
    const Phase_Support & phaseSupport,
    const SharpFeatureInfo * sharpFeatureInfo,
    const double maxSnapForEdges,
    const std::vector<stk::mesh::Entity> & intPtNodes,
    const std::vector<double> & intPtWts)
{
  bool anySnapsAreAllowedBasedOnWeight = true;
  std::vector<bool> whichSnapsAreAllowed(intPtNodes.size(), true);
  if (2 == intPtNodes.size() && maxSnapForEdges < 1.0)
  {
    whichSnapsAreAllowed[0] = intPtWts[0] > 1.-maxSnapForEdges;
    whichSnapsAreAllowed[1] = intPtWts[1] > 1.-maxSnapForEdges;
    anySnapsAreAllowedBasedOnWeight = whichSnapsAreAllowed[0] || whichSnapsAreAllowed[1];
  }

  if (anySnapsAreAllowedBasedOnWeight)
  {
    filter_which_intersection_point_nodes_are_compatible_for_snapping(mesh, auxMeta, phaseSupport, intPtNodes, whichSnapsAreAllowed);
    filter_which_intersection_point_nodes_are_compatible_for_snapping_based_on_sharp_features(mesh, sharpFeatureInfo, intPtNodes, whichSnapsAreAllowed);
  }

  return whichSnapsAreAllowed;
}

mapFromEntityToIntPtIndexAndSnapAllowed get_node_to_intersection_point_indices_and_which_snaps_allowed(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const double maxSnapForEdges,
    const std::vector<IntersectionPoint> & intersectionPoints)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());

  mapFromEntityToIntPtIndexAndSnapAllowed nodeToIntPtIndicesAndWhichSnapsAllowed;
  for (size_t intersectionPointIndex=0; intersectionPointIndex<intersectionPoints.size(); ++intersectionPointIndex)
  {
    const auto & intPtNodes = intersectionPoints[intersectionPointIndex].get_nodes();
    const auto & intPtWts = intersectionPoints[intersectionPointIndex].get_weights();
    const std::vector<bool> whichSnapsAreAllowed = which_intersection_point_nodes_are_allowed_for_snapping(mesh, auxMeta, phaseSupport, sharpFeatureInfo, maxSnapForEdges, intPtNodes, intPtWts);
    for (size_t iNode=0; iNode<intPtNodes.size(); ++iNode)
    {
      stk::mesh::Entity node = intPtNodes[iNode];
      if (mesh.bucket(node).owned())
        nodeToIntPtIndicesAndWhichSnapsAllowed[node].emplace_back(intersectionPointIndex,whichSnapsAreAllowed[iNode]);
    }
  }
  return nodeToIntPtIndicesAndWhichSnapsAllowed;
}

static void assign_is_elem_node_original_node(const std::vector<stk::mesh::Entity> & elemNodes, const stk::mesh::Entity originalNode, std::vector<bool> & isElemNodeOnInterfaceOrIsOriginal)
{
  isElemNodeOnInterfaceOrIsOriginal.clear();
  std::for_each(elemNodes.begin(), elemNodes.end(), [originalNode, &isElemNodeOnInterfaceOrIsOriginal](stk::mesh::Entity n) { isElemNodeOnInterfaceOrIsOriginal.push_back(n == originalNode);});
}

std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> determine_quality_per_node_per_domain(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed,
    const QualityMetric &qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const bool globalIDsAreParallelConsistent)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<size_t> sortedIntersectionPointIndices;
  std::vector<stk::mesh::Entity> elemNodes;
  std::vector<stk::math::Vector3d> elemNodeCoords;
  std::vector<bool> isElemNodeOnInterfaceOrIsOriginal;

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
      fill_sorted_intersection_point_indices_for_node_for_domains(mesh, coordsField, intersectionPoints, nodeIntersectionPointIndicesAndWhichSnapsAllowed, node, intPtDomains, globalIDsAreParallelConsistent, sortedIntersectionPointIndices);
      const std::set<stk::mesh::Entity> intersectedElements = get_intersected_elements(mesh, elementSelector, intersectionPoints, sortedIntersectionPointIndices);

      double qualityAfterCut = qualityMetric.get_best_value_for_metric();
      for (auto && elem : intersectedElements)
      {
        elemNodes.assign(mesh.begin_nodes(elem), mesh.end_nodes(elem));
        fill_node_locations(dim, coordsField, elemNodes, elemNodeCoords);
        assign_is_elem_node_original_node(elemNodes, node, isElemNodeOnInterfaceOrIsOriginal);
        const double elemQualityAfterCuts = estimate_quality_of_cutting_intersection_points(mesh, coordsField, elemNodes, elemNodeCoords, isElemNodeOnInterfaceOrIsOriginal, sortedIntersectionPointIndices, intersectionPoints, qualityMetric, minIntPtWeightForEstimatingCutQuality);

        if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
          qualityAfterCut = elemQualityAfterCuts;
      }

      domainsToNodesToQuality[intPtDomains][mesh.identifier(node)] = qualityAfterCut;
    }
  }

  return domainsToNodesToQuality;
}

static void
append_snap_infos_from_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed,
    const QualityMetric &qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const bool globalIDsAreParallelConsistent,
    std::vector<SnapInfo> & snapInfos)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<int> procsThatNeedToKnowAboutThisInfo;
  std::vector<size_t> globalIdsOfSnapNodeElems;

  int owner = mesh.parallel_rank();

  const auto domainsToNodesToQuality = determine_quality_per_node_per_domain(mesh, elementSelector, coordsField, intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed, qualityMetric, minIntPtWeightForEstimatingCutQuality, globalIDsAreParallelConsistent);
  const double minQualityThatIsForSureNotInverted = 0.;

  for (auto entry : nodeToIntPtIndicesAndWhichSnapsAllowed)
  {
    stk::mesh::Entity node = entry.first;
    const auto nodeIntersectionPointIndicesAndWhichSnapsAllowed = entry.second;

    if (mesh.bucket(node).owned())
    {
      const stk::math::Vector3d nodeLocation(field_data<double>(coordsField, node), dim);

      for (auto && intPtIndexAndIsSnapAllowed : nodeIntersectionPointIndicesAndWhichSnapsAllowed)
      {
        const size_t intPtIndex = intPtIndexAndIsSnapAllowed.first;
        const bool isSnapAllowed = intPtIndexAndIsSnapAllowed.second;
        const IntersectionPoint & intersectionPoint = intersectionPoints[intPtIndex];

        const auto & intPtNodes = intersectionPoint.get_nodes();

        if (isSnapAllowed && domains_already_snapped_to_node_are_also_at_intersection_point(nodesToCapturedDomains, node, intersectionPoint.get_sorted_domains()))
        {
          const stk::math::Vector3d snapLocation = compute_intersection_point_location(dim, coordsField, intersectionPoint);
          const double cutQualityEstimate = domainsToNodesToQuality.at(intersectionPoint.get_sorted_domains()).at(mesh.identifier(node));

          // For face and volume cuts, allow quality to go down to acceptable_value_for_metric because estimate is not that good
          //const double minAcceptableQuality = (nodes.size() == 2) ? cutQualityEstimate : std::min(qualityMetric.get_acceptable_value_for_metric(), cutQualityEstimate);
          const double minAcceptableQuality = std::max(minQualityThatIsForSureNotInverted, cutQualityEstimate);

          const double postSnapQuality = compute_quality_if_node_is_snapped_terminating_early_if_below_threshold(mesh, elementSelector, coordsField, node, snapLocation, qualityMetric, minAcceptableQuality);
          if (qualityMetric.is_first_quality_metric_better_than_second(postSnapQuality, minAcceptableQuality))
          {
            const size_t nodeGlobalId = mesh.identifier(node);

            fill_global_ids_of_elements_using_node(mesh, elementSelector, node, globalIdsOfSnapNodeElems);
            fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutThisInfo);

            snapInfos.emplace_back(nodeGlobalId, intPtIndex, nodeLocation, owner, procsThatNeedToKnowAboutThisInfo, globalIdsOfSnapNodeElems, postSnapQuality, snapLocation, intPtNodes.size());
          }
          else if (krinolog.shouldPrint(LOG_DEBUG))
          {
            krinolog << "Skipping snap of " << mesh.identifier(node) << " to " << snapLocation << " at " << debug_output(mesh, intersectionPoint) << " with snap quality at or below " << postSnapQuality << " and estimated cut quality " << cutQualityEstimate << stk::diag::dendl;
          }
        }
      }
    }
  }
}

std::vector<SnapInfo>
build_snap_infos_from_intersection_points(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric &qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxSnapForEdges,
    const bool globalIDsAreParallelConsistent)
{
  std::vector<SnapInfo> snapInfos;

  const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed(mesh, sharpFeatureInfo, maxSnapForEdges, intersectionPoints);
  append_snap_infos_from_intersection_points(mesh, elementSelector, nodesToCapturedDomains, intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed, qualityMetric, minIntPtWeightForEstimatingCutQuality, globalIDsAreParallelConsistent, snapInfos);

  return snapInfos;
}

void interpolate_nodal_field(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node, const FieldRef field,
    const std::vector<stk::mesh::Entity> & interpNodes,
    const std::vector<double> & interpWeights,
    std::vector<double> & scratch)
{
  const unsigned fieldLength = field.length();

  double * val = field_data<double>(field, node);
  if (nullptr == val) return;

  scratch.assign(fieldLength, 0.0);

  for (size_t iNode=0; iNode<interpNodes.size(); ++iNode)
  {
    const double * nodeVal = field_data<double>(field, interpNodes[iNode]);
    if (nullptr == nodeVal)
    {
      krinolog << "When snapping node " << mesh.identifier(node) << ", the field " << field.name() << " is missing on interpolating node " << mesh.identifier(interpNodes[iNode]) << stk::diag::dendl;
      krinolog << "Should the field " << field.name() << " be an interpolation field?" << stk::diag::dendl;
      STK_ThrowRequireMsg(false, "Interpolation field missing on interpolation node " << mesh.identifier(interpNodes[iNode]));
    }

    for (unsigned i=0; i<fieldLength; ++i)
      scratch[i] += interpWeights[iNode] * nodeVal[i];
  }

  for (unsigned i=0; i<fieldLength; ++i)
    val[i] = scratch[i];
}

void snap_nodes(const stk::mesh::BulkData & mesh,
    const FieldSet & interpolationFieldSet,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<SnapInfo> & snapInfos,
    NodeToCapturedDomainsMap & nodesToCapturedDomains )
{
  STK_ThrowRequire(mesh.parallel_size() == 1 || mesh.is_automatic_aura_on());

  std::vector< const stk::mesh::FieldBase *> interpFieldVec;
  for(auto && field : interpolationFieldSet)
    interpFieldVec.push_back(&field.field());
  stk::mesh::communicate_field_data(mesh, interpFieldVec);

  std::vector<double> scratch;
  std::vector<stk::mesh::Entity> snapNodes;
  snapNodes.reserve(snapInfos.size());

  for (auto && snapInfo : snapInfos)
  {
    if (snapInfo.get_owner() == mesh.parallel_rank())
    {
      const size_t intersectionPointIndex = snapInfo.get_intersection_point_index();
      stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id());
      snapNodes.push_back(snapNode);
      const IntersectionPoint & intersectionPoint = intersectionPoints[intersectionPointIndex];

      nodesToCapturedDomains[snapNode] = intersectionPoint.get_sorted_domains();

      const auto & nodes = intersectionPoint.get_nodes();
      const auto & weights = intersectionPoint.get_weights();

      if (krinolog.shouldPrint(LOG_DEBUG))
      {
        krinolog << "Snapping node " << snapInfo.get_node_global_id() << " to " << debug_output(mesh, intersectionPoint) << stk::diag::dendl;
      }

      for(auto && field : interpolationFieldSet)
        interpolate_nodal_field(mesh, snapNode, field, nodes, weights, scratch);
    }
  }

  stk::mesh::communicate_field_data(mesh, interpFieldVec);

  communicate_node_captured_domains_for_given_nodes(mesh, snapNodes, nodesToCapturedDomains);
}

static double interpolate_nodal_field_component(const stk::mesh::BulkData & mesh, const FieldRef field, const unsigned component, const stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & interpNodes, const std::vector<double> & interpWeights)
{
  constexpr double wtTol = 1.e-11;
  double interpVal = 0.;
  double * val = field_data<double>(field, node);
  if (nullptr != val)
  {
    for (size_t iInterpNode=0; iInterpNode<interpNodes.size(); ++iInterpNode)
    {
      const double * nodeVal = field_data<double>(field, interpNodes[iInterpNode]);
      if (nullptr == nodeVal)
      {
        if (interpWeights[iInterpNode] > wtTol)
        {
          krinolog << "When snapping/unsnapping node " << mesh.identifier(node) << " via interpolation with weight " << interpWeights[iInterpNode]
            << ", the field " << field.name() << " is missing on interpolating node " << mesh.identifier(interpNodes[iInterpNode]) << stk::diag::dendl;
          krinolog << "Should the field " << field.name() << " be an interpolation field?" << stk::diag::dendl;
          STK_ThrowRequireMsg(false, "Interpolation field missing on interpolation node " << mesh.identifier(interpNodes[iInterpNode]));
        }
      }
      else
      {
        interpVal += interpWeights[iInterpNode] * nodeVal[component];
      }
    }
  }

  return interpVal;
}

static void interpolate_field_component_on_potentially_conflicting_nodes(const stk::mesh::BulkData & mesh,
    const FieldRef field,
    const unsigned component,
    const std::vector<stk::mesh::Entity> & snapNodes,
    const std::vector<InterpolationPoint> & interpolationPoints,
    std::vector<double> & scratch)
{
  scratch.resize(snapNodes.size());
  for (size_t iSnapNode=0; iSnapNode<snapNodes.size(); ++iSnapNode)
    scratch[iSnapNode] = interpolate_nodal_field_component(mesh, field, component, snapNodes[iSnapNode], interpolationPoints[iSnapNode].get_nodes(), interpolationPoints[iSnapNode].get_weights());

  for (size_t iSnapNode=0; iSnapNode<snapNodes.size(); ++iSnapNode)
  {
    double * val = field_data<double>(field, snapNodes[iSnapNode]);
    if (nullptr != val)
    {
      val[component] = scratch[iSnapNode];
    }
  }
}

static stk::math::Vector3d compute_element_parametric_coords_at_location(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity element, const stk::math::Vector3d & location)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodeCoords;

  for (auto node : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
    nodeCoords.emplace_back(get_vector_field(mesh, coordsField, node, dim));

  return get_parametric_coordinates_of_point(nodeCoords, location);
}

static void fill_interpolation_nodes_in_element_at_parametric_coords(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity containingElem,
    const stk::math::Vector3d & containingElementParametricCoords,
    std::vector<stk::mesh::Entity> & interpNodes,
    std::vector<double> & interpWeights)
{
  const MasterElement & masterElem = MasterElementDeterminer::getMasterElement(mesh.bucket(containingElem).topology());

  interpNodes.assign(mesh.begin_nodes(containingElem), mesh.end_nodes(containingElem));
  interpWeights.assign(interpNodes.size(), 0.);
  masterElem.shape_fcn(1, containingElementParametricCoords.data(), interpWeights.data());
}

static void fill_interplation_nodes_and_weights_at_location(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & location,
    std::vector<stk::mesh::Entity> & interpNodes,
    std::vector<double> & interpWeights)
{
  stk::mesh::Entity containingElem;
  stk::math::Vector3d containingElementParametricCoords;

  double minSqrDist = std::numeric_limits<double>::max();
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (mesh.bucket(elem).member(activePart))
    {
      const stk::math::Vector3d elemParamCoords = compute_element_parametric_coords_at_location(mesh, coordsField, elem, location);
      const double elemParamSqrDist = compute_parametric_square_distance(elemParamCoords);
      if (elemParamSqrDist < minSqrDist)
      {
        minSqrDist = elemParamSqrDist;
        containingElem = elem;
        containingElementParametricCoords = elemParamCoords;
      }
    }
  }

  fill_interpolation_nodes_in_element_at_parametric_coords(mesh, containingElem, containingElementParametricCoords, interpNodes, interpWeights);
}

static std::vector<InterpolationPoint> build_interpolation_points_for_unsnapping(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField,  const FieldRef cdfemSnapField, const std::vector<stk::mesh::Entity> & snapNodes)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> interpNodes;
  std::vector<double> interpWeights;

  std::vector<InterpolationPoint> interpolationPoints;
  interpolationPoints.reserve(snapNodes.size());

  for (auto && node : snapNodes)
  {
    const stk::math::Vector3d oldSnap = get_vector_field(mesh, oldSnapDisplacements, node, dim);
    const stk::math::Vector3d currentLoc = get_vector_field(mesh, coordsField, node, dim);
    const stk::math::Vector3d unsnappedLoc = currentLoc - oldSnap;
    fill_interplation_nodes_and_weights_at_location(mesh, activePart, coordsField, node, unsnappedLoc, interpNodes, interpWeights);
    interpolationPoints.emplace_back(interpNodes, interpWeights);
  }

  return interpolationPoints;
}

static void interpolate_fields(const stk::mesh::BulkData & mesh, const FieldSet & interpFields, const std::vector<stk::mesh::Entity> & ownedNodesToInterpolate, const std::vector<InterpolationPoint> & interpolationPoints)
{
  std::vector<double> scratch;
  for (auto && field : interpFields)
    for (unsigned i=0; i<field.length(); ++i)
      interpolate_field_component_on_potentially_conflicting_nodes(mesh, field, i, ownedNodesToInterpolate, interpolationPoints, scratch);

  std::vector<const stk::mesh::FieldBase *> const_fields;
  for (auto && f : interpFields)
    const_fields.push_back(&f.field());
  stk::mesh::communicate_field_data(mesh, const_fields);
}

std::vector<stk::mesh::Entity> get_owned_nodes_to_unsnap(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> ownedSnapNodes;

  stk::mesh::Selector ownedWithField = stk::mesh::selectField(cdfemSnapField) & mesh.mesh_meta_data().locally_owned_part();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedWithField))
  {
    for(const auto & node : *bucketPtr)
    {
      const stk::math::Vector3d oldSnap = get_vector_field(mesh, oldSnapDisplacements, node, dim);
      if (oldSnap.length_squared() > 0)
      {
        ownedSnapNodes.push_back(node);
      }
    }
  }

  return ownedSnapNodes;
}

void undo_previous_snaps_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & snapFields)
{
  const std::vector<stk::mesh::Entity> ownedNodesToUnsnap = get_owned_nodes_to_unsnap(mesh, activePart, coordsField, cdfemSnapField);

  const std::vector<InterpolationPoint> interpolationPoints = build_interpolation_points_for_unsnapping(mesh, activePart, coordsField, cdfemSnapField, ownedNodesToUnsnap);

  interpolate_fields(mesh, snapFields, ownedNodesToUnsnap, interpolationPoints);
}

static stk::math::Vector3d compute_element_parametric_coords_in_previous_configuration(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef cdfemSnapField, const stk::mesh::Entity element, const stk::math::Vector3d & currentLocation)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> previousNodeCoords;

  for (auto node : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
  {
    const stk::math::Vector3d currentNodeLocation = get_vector_field(mesh, coordsField, node, dim);
    const stk::math::Vector3d oldSnap = get_vector_field(mesh, oldSnapDisplacements, node, dim);
    const stk::math::Vector3d newSnap = get_vector_field(mesh, cdfemSnapField, node, dim);
    previousNodeCoords.push_back(currentNodeLocation - newSnap + oldSnap);
  }

  return get_parametric_coordinates_of_point(previousNodeCoords, currentLocation);
}

static void fill_interpolation_nodes_and_weights_at_node_location_in_previous_configuration(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef cdfemSnapField,
    const stk::mesh::Entity node,
    std::vector<stk::mesh::Entity> & interpNodes,
    std::vector<double> & interpWeights)
{
  stk::mesh::Entity containingElem;
  stk::math::Vector3d containingElementParametricCoords;

  const int dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::math::Vector3d currentLocation = get_vector_field(mesh, coordsField, node, dim);

  double minSqrDist = std::numeric_limits<double>::max();
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (mesh.bucket(elem).member(activePart))
    {
      const stk::math::Vector3d elemParamCoords = compute_element_parametric_coords_in_previous_configuration(mesh, coordsField, cdfemSnapField, elem, currentLocation);
      const double elemParamSqrDist = compute_parametric_square_distance(elemParamCoords);
      if (elemParamSqrDist < minSqrDist)
      {
        minSqrDist = elemParamSqrDist;
        containingElem = elem;
        containingElementParametricCoords = elemParamCoords;
      }
    }
  }

  fill_interpolation_nodes_in_element_at_parametric_coords(mesh, containingElem, containingElementParametricCoords, interpNodes, interpWeights);
}

static std::vector<InterpolationPoint> build_interpolation_points_for_snapping(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, const FieldRef cdfemSnapField, const std::vector<stk::mesh::Entity> & snapNodes)
{
  stk::mesh::Entity containingElem;
  std::vector<stk::mesh::Entity> interpNodes;
  std::vector<double> interpWeights;

  std::vector<InterpolationPoint> interpolationPoints;
  interpolationPoints.reserve(snapNodes.size());

  for (auto && node : snapNodes)
  {
    fill_interpolation_nodes_and_weights_at_node_location_in_previous_configuration(mesh, activePart, coordsField, cdfemSnapField, node, interpNodes, interpWeights);
    interpolationPoints.emplace_back(interpNodes, interpWeights);
  }

  return interpolationPoints;
}

std::vector<stk::mesh::Entity> get_owned_nodes_to_snap(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> ownedSnapNodes;

  stk::mesh::Selector ownedWithField = stk::mesh::selectField(cdfemSnapField) & mesh.mesh_meta_data().locally_owned_part();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedWithField))
  {
    for(const auto & node : *bucketPtr)
    {
      const stk::math::Vector3d oldSnap(field_data<double>(oldSnapDisplacements, node), dim);
      const stk::math::Vector3d newSnap(field_data<double>(cdfemSnapField, node), dim);
      if (oldSnap.length_squared() > 0 || newSnap.length_squared() > 0)
      {
        ownedSnapNodes.push_back(node);
      }
    }
  }

  return ownedSnapNodes;
}

void snap_fields_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & interpFields)
{
  const std::vector<stk::mesh::Entity> ownedNodesToSnap = get_owned_nodes_to_snap(mesh, activePart, coordsField, cdfemSnapField);
  const std::vector<InterpolationPoint> interpolationPoints = build_interpolation_points_for_snapping(mesh, activePart, coordsField, cdfemSnapField, ownedNodesToSnap);
  interpolate_fields(mesh, interpFields, ownedNodesToSnap, interpolationPoints);
}

template<class INFO>
size_t get_global_num_infos(const std::vector<INFO> & infos, stk::ParallelMachine comm)
{
  size_t numInfos = 0;
  int rank{stk::parallel_machine_rank(comm)};
  for (const auto &info : infos)
    if (info.get_owner() == rank)
      ++numInfos;
  const size_t localNumInfos = numInfos;
  stk::all_reduce_sum(comm, &localNumInfos, &numInfos, 1);
  return numInfos;
}

void pack_owned_snap_infos_that_other_procs_need_to_know_about(std::vector<SnapInfo> &snapInfos, stk::CommSparse &commSparse)
{
    stk::pack_and_communicate(commSparse,[&]()
        {
            for(const auto &snapInfo : snapInfos)
            {
                if ( snapInfo.get_owner() == commSparse.parallel_rank() )
                {
                    for ( const int procId : snapInfo.get_procs_that_need_to_know_about_this_info())
                    {
                        if ( procId != commSparse.parallel_rank())
                        {
                            commSparse.send_buffer(procId).pack<size_t>(snapInfo.get_node_global_id());
                            commSparse.send_buffer(procId).pack<size_t>(snapInfo.get_intersection_point_index());
                            stk::pack_vector_to_proc(commSparse, snapInfo.get_procs_that_need_to_know_about_this_info(), procId);
                            stk::pack_vector_to_proc(commSparse, snapInfo.get_conflicting_ids(), procId);
                            commSparse.send_buffer(procId).pack<double>(snapInfo.get_post_worst_quality());
                            commSparse.send_buffer(procId).pack<stk::math::Vector3d>(snapInfo.get_node_location());
                            commSparse.send_buffer(procId).pack<stk::math::Vector3d>(snapInfo.get_snap_location());
                            commSparse.send_buffer(procId).pack<int>(snapInfo.get_snap_rank());
                        }
                    }
                }
            }
        });
}

void receive_snap_infos_that_this_proc_need_to_know_about_and_ghost(std::vector<SnapInfo> &snapInfos, stk::CommSparse &commSparse)
{
    stk::unpack_communications(commSparse, [&commSparse, &snapInfos](int procId)
    {
        size_t globalNodeId{0};
        commSparse.recv_buffer(procId).unpack<size_t>(globalNodeId);

        size_t intersectionPointIndex{0};
        commSparse.recv_buffer(procId).unpack<size_t>(intersectionPointIndex);

        std::vector<int> procsThatNeedToKnowAboutThisInfo;
        std::vector<size_t> globalIdsOfSnapNodeElems;
        stk::unpack_vector_from_proc(commSparse, procsThatNeedToKnowAboutThisInfo, procId);
        stk::unpack_vector_from_proc(commSparse, globalIdsOfSnapNodeElems, procId);

        double postSnapQuality{0};
        commSparse.recv_buffer(procId).unpack<double>(postSnapQuality);

        stk::math::Vector3d nodeLocation;
        commSparse.recv_buffer(procId).unpack<stk::math::Vector3d>(nodeLocation);

        stk::math::Vector3d snapLocation;
        commSparse.recv_buffer(procId).unpack<stk::math::Vector3d>(snapLocation);

        int snapRank;
        commSparse.recv_buffer(procId).unpack<int>(snapRank);

        snapInfos.emplace_back(globalNodeId,
            intersectionPointIndex,
            nodeLocation,
            procId,
            procsThatNeedToKnowAboutThisInfo,
            globalIdsOfSnapNodeElems,
            postSnapQuality,
            snapLocation,
            snapRank);
    });
}

void communicate_snap_infos_that_other_procs_need_to_know_about(std::vector<SnapInfo> &snapInfos, stk::ParallelMachine comm)
{
    stk::CommSparse commSparse(comm);

    pack_owned_snap_infos_that_other_procs_need_to_know_about(snapInfos, commSparse);
    receive_snap_infos_that_this_proc_need_to_know_about_and_ghost(snapInfos, commSparse);
}

std::vector<stk::mesh::Entity> get_sorted_nodes_modified_in_current_snapping_iteration(const stk::mesh::BulkData & mesh, const std::vector<SnapInfo> & iterationSnapInfos)
{
  std::vector<stk::mesh::Entity> sortedSnappedNodes;
  for (auto && snapInfo : iterationSnapInfos)
    sortedSnappedNodes.push_back(mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id()));
  std::sort(sortedSnappedNodes.begin(), sortedSnappedNodes.end());
  return sortedSnappedNodes;
}

std::vector<stk::mesh::EntityId> get_sorted_ids_of_owned_nodes_of_elements_of_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Entity> & nodes)
{
  std::vector<stk::mesh::Entity> nbrNodes;
  for (auto node : nodes)
    for (auto element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
      if (elementSelector(mesh.bucket(element)))
        for (auto nbrNode : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
          nbrNodes.push_back(nbrNode);

  std::vector<stk::mesh::EntityId> nbrNodeIds;
  for (auto && nbrNode : nbrNodes)
    if (mesh.bucket(nbrNode).owned())
      nbrNodeIds.push_back(mesh.identifier(nbrNode));
  stk::util::sort_and_unique(nbrNodeIds);

  return nbrNodeIds;
}

void fill_entity_ids(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & entities,
    std::vector<stk::mesh::EntityId> & entityIds)
{
  entityIds.clear();
  for (auto && entity : entities)
    entityIds.push_back(mesh.identifier(entity));
}

static void prune_snap_infos_modified_by_snap_iteration(const stk::mesh::BulkData & mesh,
    const std::vector<size_t> & oldToNewIntPts,
    const std::vector<stk::mesh::EntityId> & sortedIdsOfNodesThatNeedNewSnapInfos,
    std::vector<SnapInfo> & snapInfos)
{
  const size_t badIndex = std::numeric_limits<size_t>::max();
  int procId = mesh.parallel_rank();
  size_t newNumSnapInfos=0;
  for (auto && snapInfo : snapInfos)
  {
    if (snapInfo.get_owner() == procId)
    {
      const size_t newIntPtIndex = oldToNewIntPts[snapInfo.get_intersection_point_index()];
      if (newIntPtIndex != badIndex &&
          !std::binary_search(sortedIdsOfNodesThatNeedNewSnapInfos.begin(), sortedIdsOfNodesThatNeedNewSnapInfos.end(), snapInfo.get_node_global_id()))
      {
        snapInfo.set_intersection_point_index(newIntPtIndex);
        std::swap(snapInfo, snapInfos[newNumSnapInfos++]);
      }
    }
  }
  snapInfos.erase(snapInfos.begin()+newNumSnapInfos, snapInfos.end());
}

static mapFromEntityToIntPtIndexAndSnapAllowed get_node_to_intersection_point_indices_and_which_snaps_allowed_for_nodes_that_need_new_snap_infos(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const double maxSnapForEdges,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const std::vector<stk::mesh::EntityId> & sortedIdsOfNodesThatNeedNewSnapInfos)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());

  mapFromEntityToIntPtIndexAndSnapAllowed nodeToIntPtIndicesAndWhichSnapsAllowed;
  for (size_t intPtIndex=0; intPtIndex<intersectionPoints.size(); ++intPtIndex)
  {
    const IntersectionPoint & intPt = intersectionPoints[intPtIndex];
    const auto & intPtNodes = intPt.get_nodes();
    const auto & intPtWts = intPt.get_weights();
    const std::vector<bool> whichSnapsAreAllowed = which_intersection_point_nodes_are_allowed_for_snapping(mesh, auxMeta, phaseSupport, sharpFeatureInfo, maxSnapForEdges, intPtNodes, intPtWts);
    for (size_t iNode=0; iNode<intPtNodes.size(); ++iNode)
    {
      stk::mesh::Entity intPtNode = intPtNodes[iNode];
      if (mesh.bucket(intPtNode).owned() && std::binary_search(sortedIdsOfNodesThatNeedNewSnapInfos.begin(), sortedIdsOfNodesThatNeedNewSnapInfos.end(), mesh.identifier(intPtNode)))
        nodeToIntPtIndicesAndWhichSnapsAllowed[intPtNode].emplace_back(intPtIndex,whichSnapsAreAllowed[iNode]);
    }
  }
  return nodeToIntPtIndicesAndWhichSnapsAllowed;
}

void update_intersection_points_and_snap_infos_after_snap_iteration(const stk::mesh::BulkData & mesh,
    const InterfaceGeometry & geometry,
    const SharpFeatureInfo * sharpFeatureInfo,
    const std::vector<stk::mesh::Entity> & iterationSortedSnapNodes,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const stk::mesh::Selector & elementSelector,
    const ScaledJacobianQualityMetric & qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxSnapForEdges,
    const bool globalIDsAreParallelConsistent,
    std::vector<IntersectionPoint> & intersectionPoints,
    std::vector<SnapInfo> & snapInfos)
{
  const std::vector<size_t> oldToNewIntPts = update_intersection_points_after_snap_iteration(mesh, elementSelector, geometry, iterationSortedSnapNodes, nodesToCapturedDomains, intersectionPoints);

  const std::vector<stk::mesh::EntityId> sortedIdsOfNodesThatNeedNewSnapInfos = get_sorted_ids_of_owned_nodes_of_elements_of_nodes(mesh, elementSelector, iterationSortedSnapNodes);

  prune_snap_infos_modified_by_snap_iteration(mesh, oldToNewIntPts, sortedIdsOfNodesThatNeedNewSnapInfos, snapInfos);

  const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed_for_nodes_that_need_new_snap_infos(mesh, sharpFeatureInfo, maxSnapForEdges, intersectionPoints, sortedIdsOfNodesThatNeedNewSnapInfos);

  append_snap_infos_from_intersection_points(mesh, elementSelector, nodesToCapturedDomains, intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed, qualityMetric, minIntPtWeightForEstimatingCutQuality, globalIDsAreParallelConsistent, snapInfos);
}

NodeToCapturedDomainsMap snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent,
    const double snappingSharpFeatureAngleInDegrees,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxSnapForEdges)
{/* %TRACE[ON]% */ Trace trace__("krino::snap_as_much_as_possible_while_maintaining_quality()"); /* %TRACE% */

    const ScaledJacobianQualityMetric qualityMetric;
    size_t iteration{0};
    NodeToCapturedDomainsMap nodesToCapturedDomains;
    stk::ParallelMachine comm = mesh.parallel();
    std::unique_ptr<SharpFeatureInfo> sharpFeatureInfo;
    if (snappingSharpFeatureAngleInDegrees > 0.)
    {
      sharpFeatureInfo = std::make_unique<SharpFeatureInfo>();
      const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
      sharpFeatureInfo->find_sharp_features(mesh, coordsField, elementSelector, std::cos(snappingSharpFeatureAngleInDegrees*M_PI/180.));
    }

    std::vector<IntersectionPoint> intersectionPoints;
    geometry.store_phase_for_uncut_elements(mesh);
    intersectionPoints = build_all_intersection_points(mesh, elementSelector, geometry, nodesToCapturedDomains);
    std::vector<SnapInfo> snapInfos = build_snap_infos_from_intersection_points(mesh, sharpFeatureInfo.get(), elementSelector, nodesToCapturedDomains, intersectionPoints, qualityMetric, minIntPtWeightForEstimatingCutQuality, maxSnapForEdges, globalIDsAreParallelConsistent);

    while (true)
    {
      krinolog << "Snapping To Geometry Iteration " << std::to_string(++iteration) << stk::diag::dendl;

      bool done = stk::is_true_on_all_procs(comm, snapInfos.empty());
      if ( done )
          break;

      communicate_snap_infos_that_other_procs_need_to_know_about(snapInfos, comm);

      const std::vector<SnapInfo> independentSnapInfos = find_snap_info_independent_sets(snapInfos, qualityMetric, comm);

      krinolog << "  Snapping " << get_global_num_infos(independentSnapInfos, comm) << " of " << get_global_num_infos(snapInfos, comm) << " snap candidates." << stk::diag::dendl;

      geometry.store_phase_for_elements_that_will_be_uncut_after_snapping(mesh, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);
      snap_nodes(mesh, interpolationFields, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);

      const std::vector<stk::mesh::Entity> iterationSortedSnapNodes = get_sorted_nodes_modified_in_current_snapping_iteration(mesh, independentSnapInfos);

      update_intersection_points_and_snap_infos_after_snap_iteration(mesh, geometry, sharpFeatureInfo.get(), iterationSortedSnapNodes, nodesToCapturedDomains, elementSelector, qualityMetric, minIntPtWeightForEstimatingCutQuality, maxSnapForEdges, globalIDsAreParallelConsistent, intersectionPoints, snapInfos);
    }

    krinolog << "After snapping quality is " << compute_mesh_quality(mesh, elementSelector, qualityMetric) << stk::diag::dendl;

    return nodesToCapturedDomains;
}
}



