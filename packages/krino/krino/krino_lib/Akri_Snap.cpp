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
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SnapIndependentSetFinder.hpp>
#include <Akri_SnapInfo.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <memory>
#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"

namespace krino
{

static stk::math::Vector3d compute_snap_location(const std::vector<stk::math::Vector3d> & nodeLocations, const std::vector<double> & weights)
{
  stk::math::Vector3d snapLocation{stk::math::Vector3d::ZERO};
  for (size_t i=0; i<nodeLocations.size(); ++i)
    snapLocation += weights[i] * nodeLocations[i];
  return snapLocation;
}

static void fill_node_locations(const int dim, const FieldRef coordsField, const std::vector<stk::mesh::Entity> & nodes, std::vector<stk::math::Vector3d> & nodeLocations)
{
  nodeLocations.clear();
  for (auto node : nodes)
    nodeLocations.emplace_back(field_data<double>(coordsField, node), dim);
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

static double compute_quality_if_node_is_snapped_terminating_early_if_below_threshold(const stk::mesh::BulkData & mesh,
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

static stk::math::Vector3d compute_intersection_point_location(const int dim, const FieldRef coordsField, const IntersectionPoint & intersectionPoint)
{
  std::vector<stk::math::Vector3d> intPtNodeLocations;
  fill_node_locations(dim, coordsField, intersectionPoint.get_nodes(), intPtNodeLocations);
  return compute_snap_location(intPtNodeLocations, intersectionPoint.get_weights());
}

static bool element_has_all_nodes(const std::vector<stk::mesh::Entity> & elemNodes, const std::vector<stk::mesh::Entity> & nodesToFind)
{
  for (auto && nodeToFind : nodesToFind)
    if (std::find(elemNodes.begin(), elemNodes.end(), nodeToFind) == elemNodes.end())
      return false;
  return true;
}

static double estimate_quality_of_cutting_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elemNodes,
    const std::vector<Vector3d> & elemNodeCoords,
    const std::vector<size_t> intersectionPointIndices,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric &qualityMetric)
{
  if (intersectionPointIndices.empty())
    return std::max(0., qualityMetric.get_element_quality_metric(elemNodeCoords));

  // apply the front intersection point to element and recursively call with remaining intersection points
  const IntersectionPoint & intPtToApply = intersectionPoints[*intersectionPointIndices.begin()];
  const std::vector<size_t> remainingIntersectionPointIndices(intersectionPointIndices.begin()+1, intersectionPointIndices.end());

  const stk::math::Vector3d intPtLocation = compute_intersection_point_location(mesh.mesh_meta_data().spatial_dimension(), coordsField, intPtToApply);

  const auto & intPtNodes = intPtToApply.get_nodes();
  double qualityAfterCut = qualityMetric.get_best_value_for_metric();
  if (element_has_all_nodes(elemNodes, intPtNodes))
  {
    std::vector<stk::mesh::Entity> cutElemNodes;
    std::vector<stk::math::Vector3d> cutElemNodeCoords;

    for (auto intPtNode : intPtNodes)
    {
      cutElemNodes.clear();
      cutElemNodeCoords.clear();
      for (size_t nodeIndex=0; nodeIndex<elemNodes.size(); ++nodeIndex)
      {
        if (elemNodes[nodeIndex] == intPtNode)
        {
          cutElemNodes.push_back(stk::mesh::Entity());
          cutElemNodeCoords.push_back(intPtLocation);
        }
        else
        {
          cutElemNodes.push_back(elemNodes[nodeIndex]);
          cutElemNodeCoords.push_back(elemNodeCoords[nodeIndex]);
        }
      }

      const double elemQualityAfterCuts = estimate_quality_of_cutting_intersection_points(mesh, coordsField, cutElemNodes, cutElemNodeCoords, remainingIntersectionPointIndices, intersectionPoints, qualityMetric);
      if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
        qualityAfterCut = elemQualityAfterCuts;
    }
  }
  else
  {
    qualityAfterCut = estimate_quality_of_cutting_intersection_points(mesh, coordsField, elemNodes, elemNodeCoords, remainingIntersectionPointIndices, intersectionPoints, qualityMetric);
  }

  return qualityAfterCut;
}

static bool parts_are_compatible_for_snapping(const stk::mesh::BulkData & mesh, const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & interpNodes)
{
  for (auto && interpNode : interpNodes)
    if (interpNode != node && !parts_are_compatible_for_snapping_when_ignoring_phase(mesh, auxMeta, phaseSupport, node, interpNode))
      return false;
  return true;
}

static double get_node_intersection_point_weight(const IntersectionPoint & intersectionPoint, stk::mesh::Entity node)
{
  const std::vector<stk::mesh::Entity> & nodes = intersectionPoint.get_nodes();
  const auto iter = std::find(nodes.begin(), nodes.end(), node);
  ThrowRequire(iter != nodes.end());
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
    const std::vector<size_t> & candidatesIntersectionPointIndices,
    const stk::mesh::Entity node,
    const std::vector<int> & domains,
    const bool globalIDsAreParallelConsistent,
    std::vector<size_t> & sortedIntersectionPointIndices)
{
  sortedIntersectionPointIndices.clear();
  for (auto && intPtIndex : candidatesIntersectionPointIndices)
  {
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

static std::map<stk::mesh::Entity, std::vector<size_t>> get_node_to_intersection_point_indices(const stk::mesh::BulkData & mesh,
    const std::vector<IntersectionPoint> & intersectionPoints)
{
  std::map<stk::mesh::Entity, std::vector<size_t>> nodeToInsersectionPointIndices;
  for (size_t intersectionPointIndex=0; intersectionPointIndex<intersectionPoints.size(); ++intersectionPointIndex)
    for (auto && node : intersectionPoints[intersectionPointIndex].get_nodes())
      if (mesh.bucket(node).owned())
        nodeToInsersectionPointIndices[node].push_back(intersectionPointIndex);
  return nodeToInsersectionPointIndices;
}

std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> determine_quality_per_node_per_domain(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric &qualityMetric,
    const bool globalIDsAreParallelConsistent)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const auto nodeToInsersectionPointIndices = get_node_to_intersection_point_indices(mesh, intersectionPoints);

  std::vector<size_t> sortedIntersectionPointIndices;
  std::vector<stk::mesh::Entity> elemNodes;
  std::vector<stk::math::Vector3d> elemNodeCoords;

  std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> domainsToNodesToQuality;
  for (auto entry : nodeToInsersectionPointIndices)
  {
    stk::mesh::Entity node = entry.first;
    const auto nodeIntersectionPointIndices = entry.second;

    std::set<std::vector<int>> nodeIntPtDomains;
    for (auto && intPtIndex : nodeIntersectionPointIndices)
      nodeIntPtDomains.insert(intersectionPoints[intPtIndex].get_sorted_domains());

    for (auto && intPtDomains : nodeIntPtDomains)
    {
      fill_sorted_intersection_point_indices_for_node_for_domains(mesh, coordsField, intersectionPoints, nodeIntersectionPointIndices, node, intPtDomains, globalIDsAreParallelConsistent, sortedIntersectionPointIndices);
      const std::set<stk::mesh::Entity> intersectedElements = get_intersected_elements(mesh, elementSelector, intersectionPoints, sortedIntersectionPointIndices);

      double qualityAfterCut = qualityMetric.get_best_value_for_metric();
      for (auto && elem : intersectedElements)
      {
        elemNodes.assign(mesh.begin_nodes(elem), mesh.end_nodes(elem));
        fill_node_locations(dim, coordsField, elemNodes, elemNodeCoords);
        const double elemQualityAfterCuts = estimate_quality_of_cutting_intersection_points(mesh, coordsField, elemNodes, elemNodeCoords, sortedIntersectionPointIndices, intersectionPoints, qualityMetric);

        if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterCut, elemQualityAfterCuts))
          qualityAfterCut = elemQualityAfterCuts;
      }

      domainsToNodesToQuality[intPtDomains][mesh.identifier(node)] = qualityAfterCut;
    }
  }

  return domainsToNodesToQuality;
}

std::vector<SnapInfo>
build_snap_infos_from_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric &qualityMetric,
    const bool globalIDsAreParallelConsistent)
{
  std::vector<SnapInfo> snapInfos;

  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodeLocations;
  std::vector<int> procsThatNeedToKnowAboutThisInfo;
  std::vector<size_t> globalIdsOfSnapNodeElems;

  int owner = mesh.parallel_rank();

  const auto domainsToNodesToQuality = determine_quality_per_node_per_domain(mesh, elementSelector, coordsField, intersectionPoints, qualityMetric, globalIDsAreParallelConsistent);

  for (size_t intersectionPointIndex=0; intersectionPointIndex<intersectionPoints.size(); ++intersectionPointIndex)
  {
    const IntersectionPoint & intersectionPoint = intersectionPoints[intersectionPointIndex];
    const auto & nodes = intersectionPoint.get_nodes();
    fill_node_locations(dim, coordsField, nodes, nodeLocations);
    const stk::math::Vector3d snapLocation = compute_snap_location(nodeLocations, intersectionPoint.get_weights());
    const auto & nodesToQualityIter = domainsToNodesToQuality.find(intersectionPoint.get_sorted_domains());
    for (size_t nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
    {
      stk::mesh::Entity node = nodes[nodeIndex];

      if (mesh.bucket(node).owned() &&
          domains_already_snapped_to_node_are_also_at_intersection_point(nodesToCapturedDomains, node, intersectionPoint.get_sorted_domains()) &&
          parts_are_compatible_for_snapping(mesh, auxMeta, phaseSupport, node, nodes))
      {
        ThrowAssert(nodesToQualityIter != domainsToNodesToQuality.end());
        const auto & nodesToQuality = nodesToQualityIter->second;
        const double cutQualityEstimate = nodesToQuality.at(mesh.identifier(node));

        // For face and volume cuts, allow quality to go down to acceptable_value_for_metric because estimate is not that good
        //const double minAcceptableQuality = (nodes.size() == 2) ? cutQualityEstimate : std::min(qualityMetric.get_acceptable_value_for_metric(), cutQualityEstimate);
        const double minAcceptableQuality = cutQualityEstimate;

        const double postSnapQuality = compute_quality_if_node_is_snapped_terminating_early_if_below_threshold(mesh, elementSelector, coordsField, node, snapLocation, qualityMetric, minAcceptableQuality);
        if (qualityMetric.is_first_quality_metric_better_than_second(postSnapQuality, minAcceptableQuality))
        {
          const size_t nodeGlobalId = mesh.identifier(node);

          fill_global_ids_of_elements_using_node(mesh, elementSelector, node, globalIdsOfSnapNodeElems);
          fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutThisInfo);

          snapInfos.emplace_back(nodeGlobalId, intersectionPointIndex, nodeLocations[nodeIndex], owner, procsThatNeedToKnowAboutThisInfo, globalIdsOfSnapNodeElems, postSnapQuality, snapLocation, nodes.size());
        }
        else if (krinolog.shouldPrint(LOG_DEBUG))
        {
          krinolog << "Skipping snap of " << mesh.identifier(node) << " to " << snapLocation << " at " << debug_output(mesh, intersectionPoint) << " with snap quality at or below " << postSnapQuality << " and estimated cut quality " << cutQualityEstimate << stk::diag::dendl;
        }
      }
    }
  }

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
      ThrowRequireMsg(false, "Interpolation field missing on interpolation node " << mesh.identifier(interpNodes[iNode]));
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
  ThrowRequire(mesh.parallel_size() == 1 || mesh.is_automatic_aura_on());

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

double determine_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const QualityMetric &qualityMetric)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::vector<stk::math::Vector3d> nodeLocations;

  double quality = qualityMetric.get_best_value_for_metric();
  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
  {
    if (bucket->topology().base() == stk::topology::TETRAHEDRON_4 || bucket->topology().base() == stk::topology::TRIANGLE_3_2D)
    {
      for (auto && element : *bucket)
      {
        fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);
        const double elementQuality = qualityMetric.get_element_quality_metric(nodeLocations);
        quality = std::min(quality, elementQuality);
      }
    }
  }

  const double localQuality = quality;
  stk::all_reduce_min(mesh.parallel(), &localQuality, &quality, 1);

  return quality;
}

NodeToCapturedDomainsMap snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent)
{/* %TRACE[ON]% */ Trace trace__("krino::snap_as_much_as_possible_while_maintaining_quality()"); /* %TRACE% */

    const ScaledJacobianQualityMetric qualityMetric;
    size_t iteration{0};
    NodeToCapturedDomainsMap nodesToCapturedDomains;
    stk::ParallelMachine comm = mesh.parallel();

    std::vector<IntersectionPoint> intersectionPoints;
    geometry.store_phase_for_uncut_elements(mesh);
    intersectionPoints = build_all_intersection_points(mesh, geometry, nodesToCapturedDomains);

    while (true)
    {
      krinolog << "Snapping To Geometry Iteration " << std::to_string(++iteration) << stk::diag::dendl;

      std::vector<SnapInfo> snapInfos = build_snap_infos_from_intersection_points(mesh, elementSelector, nodesToCapturedDomains, intersectionPoints, qualityMetric, globalIDsAreParallelConsistent);

      bool done = stk::is_true_on_all_procs(comm, snapInfos.empty());
      if ( done )
          break;

      communicate_snap_infos_that_other_procs_need_to_know_about(snapInfos, comm);

      const std::vector<SnapInfo> independentSnapInfos = find_snap_info_independent_sets(snapInfos, qualityMetric, comm);

      krinolog << "  Snapping " << get_global_num_infos(independentSnapInfos, comm) << " of " << get_global_num_infos(snapInfos, comm) << " snap candidates." << stk::diag::dendl;

      geometry.store_phase_for_elements_that_will_be_uncut_after_snapping(mesh, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);
      snap_nodes(mesh, interpolationFields, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);

      const std::vector<stk::mesh::Entity> iterationSortedSnapNodes = get_sorted_nodes_modified_in_current_snapping_iteration(mesh, independentSnapInfos);
      update_intersection_points_after_snap_iteration(mesh, geometry, iterationSortedSnapNodes, nodesToCapturedDomains, intersectionPoints);
    }

    krinolog << "After snapping quality is " << determine_quality(mesh, elementSelector, qualityMetric) << stk::diag::dendl;

    return nodesToCapturedDomains;
}
}



