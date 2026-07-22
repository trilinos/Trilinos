/*
 * Akri_Smooth.cpp
 *
 *  Created on: Nov 8, 2024
 *      Author: drnoble
 */
#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_QualityMetricWithSensitivities.hpp>
#include <Akri_Smooth.hpp>
#include <Akri_SmoothIndependentSetFinder.hpp>
#include <Akri_SmoothInfo.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

double calculate_quality_before_smoothing(const stk::mesh::BulkData &mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Selector & elemSelector, const stk::mesh::Entity node, const QualityMetric &qualityMetric)
{
  StkMeshEntities nodeElements{mesh.begin_elements(node), mesh.end_elements(node)};
  return compute_minimum_element_quality(mesh, coordsField, elemSelector, nodeElements, qualityMetric);
}

static double compute_node_delta_squared(const stk::mesh::BulkData & mesh,
                                         const CoordinatesFieldRef coordsField,
    const stk::mesh::Entity node,
    const stk::math::Vector3d nodeLoc,
    const stk::mesh::Entity elem)
{
  double sqrDeltaNode = 0.;
  for (auto elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
  {
    if (elemNode != node)
    {
      const stk::math::Vector3d elemNodeLoc = get_vector_field(mesh, coordsField, elemNode, coordsField.dim());
      sqrDeltaNode += (elemNodeLoc - nodeLoc).length_squared();
    }
  }
  return sqrDeltaNode;
}

std::tuple<double, std::array<double,3>> tet_volume_and_sensitivity_to_nth_node(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity elem, const int sensitivityNodeOrdinal)
{
  const auto elementNodeCoords = gather_tet_coordinates(mesh, elem, coordsField);
  return tet_volume_and_sensitivity_to_nth_vertex(elementNodeCoords, sensitivityNodeOrdinal);
}

stk::math::Vector3d smoothed_node_location_using_ODT(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLoc)
{
  double patchVol = 0.;
  stk::math::Vector3d delta = stk::math::Vector3d::ZERO;
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    const stk::mesh::Bucket & elemBucket = mesh.bucket(elem);
    if (elemBucket.topology() == stk::topology::TETRAHEDRON_4 && elemSelector(elemBucket))
    {
      const auto & [vol, sensitivity] = tet_volume_and_sensitivity_to_nth_node(mesh, coordsField, elem, get_entity_node_ordinal(mesh, elem, node));
      patchVol += vol;
      const double sqrDeltaNode = compute_node_delta_squared(mesh, coordsField, node, nodeLoc, elem);
      delta += sqrDeltaNode * stk::math::Vector3d(sensitivity[0], sensitivity[1], sensitivity[2]);
    }
  }
  if (patchVol != 0.)
    delta *= (-0.5/patchVol);

  return nodeLoc+delta;
}

std::tuple<double, std::array<double,3>> tet_mean_ratio_quality_and_sensitivity_to_nth_node(const stk::mesh::BulkData & mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Entity elem, const int sensitivityNodeOrdinal, const stk::math::Vector3d & nodeLocation)
{
  auto elementNodeCoords = gather_tet_coordinates(mesh, elem, coordsField);
  elementNodeCoords[sensitivityNodeOrdinal] = nodeLocation;
  return MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivity_to_nth_vertex(elementNodeCoords, sensitivityNodeOrdinal);
}

double mean_ratio_element_quality(const stk::mesh::BulkData & mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Entity elem, const stk::mesh::Entity node, const stk::math::Vector3d & nodeLocation)
{
  MeanRatioQualityMetric qualityMetric;
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  std::vector<stk::math::Vector3d> nodeLocations;
  nodeLocations.reserve(elemNodes.size());
  for (auto elemNode : elemNodes)
  {
    if (elemNode == node)
      nodeLocations.push_back(nodeLocation);
    else
      nodeLocations.emplace_back(field_data<double>(coordsField, elemNode), coordsField.dim());
  }
  return qualityMetric.get_element_quality_metric(coordsField.dim(), nodeLocations);
}

double compute_worst_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLocation)
{
  double obj = 0.;
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  for (auto nodeElem : nodeElems)
  {
    if (elemSelector(mesh.bucket(nodeElem)))
    {
      const double quality = mean_ratio_element_quality(mesh, coordsField, nodeElem, node, nodeLocation);
      obj = std::max(obj, 1.-quality);
    }
  }
  return obj;
}

#if 0
std::tuple<double, std::array<double,3>> compute_worst_mean_ratio_quality_objective_and_sensitivity(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLocation)
{
  stk::mesh::Entity worstElem;
  double obj = 0.;
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  for (auto nodeElem : nodeElems)
  {
    if (elemSelector(mesh.bucket(nodeElem)))
    {
      const double quality = mean_ratio_element_quality(mesh, coordsField, nodeElem, node, nodeLocation);
      if (1.-quality > obj || (1.-quality == obj && (!mesh.is_valid(worstElem) || mesh.identifier(nodeElem) < mesh.identifier(worstElem))))
      {
        obj = 1.-quality;
        worstElem = nodeElem;
      }
    }
  }

  std::array<double,3> sens{0., 0., 0.};

  if (mesh.is_valid(worstElem))
  {
    const auto & [quality, sensitivity] = tet_mean_ratio_quality_and_sensitivity_to_nth_node(mesh, coordsField, worstElem, get_entity_node_ordinal(mesh, worstElem, node), nodeLocation);
    for (int d=0; d<3; ++d)
      sens[d] = -sensitivity[d];
  }
  return std::make_tuple(obj, sens);
}
#endif
#if 1
std::tuple<double, stk::math::Vector3d> compute_worst_mean_ratio_quality_objective_and_sensitivity(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLocation)
{
  double worstQuality = 1.;
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  std::vector<std::pair<stk::mesh::Entity,double>> elemsAndQualities; // good candidate for small vector
  elemsAndQualities.reserve(nodeElems.size());
  for (auto nodeElem : nodeElems)
  {
    if (elemSelector(mesh.bucket(nodeElem)))
    {
      const double quality = mean_ratio_element_quality(mesh, coordsField, nodeElem, node, nodeLocation);
      elemsAndQualities.push_back(std::make_pair(nodeElem, quality));
      worstQuality = std::min(worstQuality, quality);
    }
  }

  stk::math::Vector3d sens = stk::math::Vector3d::ZERO;

  // If multiple elements have the same worst quality, we need to include them in the sensitivity or else the gradient might not point in the direction of improvement
  unsigned elemCount = 0;
  const double qualityCloseEnoughToIncludeInSens = 1.01 * worstQuality;
  for (const auto & [nodeElem, elemQuality] : elemsAndQualities)
  {
    if (elemQuality < qualityCloseEnoughToIncludeInSens)
    {
      elemCount++;
      const auto & [qual, sensitivity] = tet_mean_ratio_quality_and_sensitivity_to_nth_node(mesh, coordsField, nodeElem, get_entity_node_ordinal(mesh, nodeElem, node), nodeLocation);
      for (int d=0; d<3; ++d)
        sens[d] -= sensitivity[d];
    }
  }

  if (elemCount > 1)
    sens /= elemCount;

  return std::make_tuple(1.-worstQuality, sens);
}
#endif

NodeObjFn build_function_for_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node)
{
  auto fn = [&mesh, coordsField, elemSelector, node](const stk::math::Vector3d & nodeLoc)
  {
    return compute_worst_mean_ratio_quality_objective(mesh, coordsField, elemSelector, node, nodeLoc);
  };
  return fn;
}

NodeObjFnSens build_function_for_mean_ratio_quality_objective_sensitivity(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const NodeSearchDirectionFilter & search_direction_filter)
{
  auto fn = [&mesh, coordsField, elemSelector, node, &search_direction_filter](const stk::math::Vector3d & nodeLoc, stk::math::Vector3d & grad)
  {
    const auto & [origObj, sensitivity] = compute_worst_mean_ratio_quality_objective_and_sensitivity(mesh, coordsField, elemSelector, node, nodeLoc);
    grad = sensitivity;
    search_direction_filter(node, grad);
  };
  return fn;
}

void smoothed_node_location_using_optimized_mean_ratio_smoothing(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const OptimizeNodeLocation & optimize_node_location,
    const NodeSearchDirectionFilter & search_direction_filter,
    stk::math::Vector3d & nodeLocation)
{
  const auto objective_at_loc = build_function_for_mean_ratio_quality_objective(mesh, coordsField, elemSelector, node);
  const auto sensitivity_at_loc = build_function_for_mean_ratio_quality_objective_sensitivity(mesh, coordsField, elemSelector, node, search_direction_filter);

  optimize_node_location(objective_at_loc, sensitivity_at_loc, nodeLocation);
}

void set_node_coordinates(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & newCoords)
{
  double * nodeCoords = get_field_data(mesh, coordsField, node);
  for (unsigned d=0; d<coordsField.dim(); ++d)
    nodeCoords[d] = newCoords[d];
}

void smooth_nodes(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<SmoothInfo> & smoothInfos)
{
  for (auto && smoothInfo : smoothInfos)
  {
    const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, smoothInfo.get_node_global_id());
    if (mesh.bucket(node).owned())
      set_node_coordinates(mesh, coordsField, mesh.get_entity(stk::topology::NODE_RANK, smoothInfo.get_node_global_id()), smoothInfo.get_post_smooth_location());
  }
  parallel_sync_fields(mesh, {&coordsField.field()});
}

static std::vector<stk::mesh::Entity> get_selected_node_neighbors(const stk::mesh::BulkData &mesh,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Selector & nodeSelector,
    const std::vector<stk::mesh::Entity> & nodes)
{
  std::vector<stk::mesh::Entity> selectedNodeNbrs;
  for (const auto &node : nodes)
    for (auto nodeElem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
      if (elemSelector(mesh.bucket(nodeElem)))
        for (auto elemNode : StkMeshEntities{mesh.begin_nodes(nodeElem), mesh.end_nodes(nodeElem)})
          if (elemNode != node && nodeSelector(mesh.bucket(elemNode)))
            selectedNodeNbrs.push_back(elemNode);

  stk::util::sort_and_unique(selectedNodeNbrs);
  return selectedNodeNbrs;
}

std::vector<stk::mesh::Entity> get_nodes_to_consider_smoothing_for_next_iteration(const stk::mesh::BulkData &mesh,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Selector & owndNodeSelector,
    const std::vector<SmoothInfo> &smoothInfosJustExecuted)
{
  std::vector<stk::mesh::Entity> smoothedNodes;
  smoothedNodes.reserve(smoothInfosJustExecuted.size());
  for (const auto &info : smoothInfosJustExecuted)
  {
    const stk::mesh::Entity smoothedNode = mesh.get_entity(stk::topology::NODE_RANK, info.get_node_global_id());
    STK_ThrowAssert(mesh.is_valid(smoothedNode));
    smoothedNodes.push_back(smoothedNode);
  }

  const std::vector<stk::mesh::Entity> nodesToConsider = get_selected_node_neighbors(mesh, elemSelector, owndNodeSelector, smoothedNodes);

  return nodesToConsider;
}

void pack_owned_smoothing_infos_that_other_procs_need_to_know_about(const std::vector<SmoothInfo> &smoothInfos, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for(const auto &smoothInfo : smoothInfos)
    {
      for ( const int procId : smoothInfo.get_procs_that_need_to_know_about_this_info())
      {
        if ( procId != commSparse.parallel_rank())
        {
          commSparse.send_buffer(procId).pack(smoothInfo.get_node_global_id());
          commSparse.send_buffer(procId).pack(smoothInfo.get_owner());
          commSparse.send_buffer(procId).pack<stk::math::Vector3d>(smoothInfo.get_pre_smooth_location());
          commSparse.send_buffer(procId).pack<stk::math::Vector3d>(smoothInfo.get_post_smooth_location());
          commSparse.send_buffer(procId).pack<double>(smoothInfo.get_improvement());
          const auto & elementIds = smoothInfo.get_ids_of_elements_impacted_by_smoothing();
          commSparse.send_buffer(procId).pack<size_t>(elementIds.size());
          for (auto elementId : elementIds)
            commSparse.send_buffer(procId).pack(elementId);
          const auto & procsThatNeedToKnowAboutInfo = smoothInfo.get_procs_that_need_to_know_about_this_info();
          commSparse.send_buffer(procId).pack<size_t>(procsThatNeedToKnowAboutInfo.size());
          for (auto procThatNeedsToKnowAboutThisInfo : procsThatNeedToKnowAboutInfo)
            commSparse.send_buffer(procId).pack(procThatNeedsToKnowAboutThisInfo);
        }
      }
    }
  });
}

void receive_smoothing_infos_that_this_proc_need_to_know_about(std::vector<SmoothInfo> &smoothInfos, stk::CommSparse& commSparse)
{
  stk::unpack_communications(commSparse, [&commSparse, &smoothInfos](int procId)
  {
    stk::mesh::EntityId nodeId{0};
    commSparse.recv_buffer(procId).unpack(nodeId);

    int owner = -1;
    commSparse.recv_buffer(procId).unpack(owner);

    stk::math::Vector3d preSmoothedLocation;
    commSparse.recv_buffer(procId).unpack(preSmoothedLocation);

    stk::math::Vector3d postSmoothedLocation;
    commSparse.recv_buffer(procId).unpack(postSmoothedLocation);

    double improvement = 0.;
    commSparse.recv_buffer(procId).unpack(improvement);

    size_t numElements = 0;
    commSparse.recv_buffer(procId).unpack(numElements);

    std::vector<stk::mesh::EntityId> elementIds(numElements);
    for ( size_t i=0; i<numElements; ++i )
        commSparse.recv_buffer(procId).unpack(elementIds[i]);

    size_t numProcsThatNeedToKnowAboutThisInfo = 0;
    commSparse.recv_buffer(procId).unpack(numProcsThatNeedToKnowAboutThisInfo);

    std::vector<int> procsThatNeedToKnowAboutInfo(numProcsThatNeedToKnowAboutThisInfo);
    for ( size_t i=0; i<numProcsThatNeedToKnowAboutThisInfo; ++i )
        commSparse.recv_buffer(procId).unpack(procsThatNeedToKnowAboutInfo[i]);

    smoothInfos.emplace_back(nodeId, owner, preSmoothedLocation, postSmoothedLocation, improvement, elementIds, procsThatNeedToKnowAboutInfo);
  });
}

void communicate_smoothing_infos_that_other_procs_need_to_know_about(std::vector<SmoothInfo> &smoothInfos, stk::ParallelMachine comm)
{
    stk::CommSparse commSparse(comm);

    pack_owned_smoothing_infos_that_other_procs_need_to_know_about(smoothInfos, commSparse);
    receive_smoothing_infos_that_this_proc_need_to_know_about(smoothInfos, commSparse);
}

double compute_quality_improvement(const stk::mesh::BulkData &mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Selector & elemSelector, const stk::mesh::Entity node, const stk::math::Vector3d & postSmoothedLocation, const QualityMetric &qualityMetric)
{
  const double qualityBeforeSmoothing = calculate_quality_before_smoothing(mesh, coordsField, elemSelector, node, qualityMetric);
  const double qualityAfterSmoothing = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(mesh, elemSelector, coordsField, node, postSmoothedLocation, qualityMetric, qualityBeforeSmoothing);
  if (qualityMetric.is_first_quality_metric_better_than_second(qualityAfterSmoothing, qualityBeforeSmoothing))
    return (1.-qualityBeforeSmoothing)*(1.-qualityBeforeSmoothing) - (1.-qualityAfterSmoothing)*(1.-qualityAfterSmoothing);
  return 0.;
}

using FindSmoothedNodeLocation =
    std::function<stk::math::Vector3d(const stk::mesh::BulkData &, const CoordinatesFieldRef,const stk::mesh::Selector &, const stk::mesh::Entity, const stk::math::Vector3d &)>;

std::vector<SmoothInfo> find_nodes_to_smooth(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const std::vector<stk::mesh::Entity> & nodesToConsider,
    const QualityMetric &qualityMetric,
    const FindSmoothedNodeLocation & find_smoothed_node_location)
{
  std::vector<SmoothInfo> smoothInfos;

  std::vector<stk::mesh::EntityId> elementIds;
  std::vector<int> procsThatNeedToKnowAboutThisInfo;

  for(const auto & node : nodesToConsider)
  {
    const stk::math::Vector3d preSmoothedLocation = get_vector_field(mesh, coordsField, node, coordsField.dim());
    const stk::math::Vector3d postSmoothedLocation = find_smoothed_node_location(mesh, coordsField, elemSelector, node, preSmoothedLocation);
    const double qualityImprovement = compute_quality_improvement(mesh, coordsField, elemSelector, node, postSmoothedLocation, qualityMetric);

    if (qualityImprovement > 0.)
    {
      elementIds.clear();
      for (auto && element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
        if (elemSelector(mesh.bucket(element)))
          elementIds.push_back(mesh.identifier(element));

      fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutThisInfo);

      smoothInfos.emplace_back(mesh.identifier(node), mesh.parallel_rank(), preSmoothedLocation, postSmoothedLocation, qualityImprovement, elementIds, procsThatNeedToKnowAboutThisInfo);
    }
  }

  return smoothInfos;
}

void improve_quality_by_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & nodeSelector,
    const stk::mesh::Selector & elemSelector,
    const FindSmoothedNodeLocation & find_smoothed_node_location,
    const QualityMetric & qualityMetric,
    const size_t maxNumSmoothIterations)
{
  if (2 == coordsField.dim())
  {
    krinolog << "Smoothing is only supported in 3d." << stk::diag::dendl;
    return;
  }
  STK_ThrowRequireMsg(mesh.parallel_size() == 1 || mesh.is_automatic_aura_on(), "Smoothing requires aura because method to find quality when node is moved requires it.");

  parallel_sync_fields(mesh, {&coordsField.field()});

  krinolog << "Before smoothing quality is " << compute_mesh_quality(mesh, elemSelector, coordsField, qualityMetric) << stk::diag::dendl;

  const stk::mesh::Selector ownedNodeSelector = mesh.mesh_meta_data().locally_owned_part() & nodeSelector;

  std::vector< stk::mesh::Entity> nodesToConsider;
  stk::mesh::get_selected_entities( ownedNodeSelector, mesh.buckets( stk::topology::NODE_RANK ), nodesToConsider );

  for ( size_t i=0; i<maxNumSmoothIterations; i++)
  {
    std::vector<SmoothInfo> smoothInfos = find_nodes_to_smooth(mesh, coordsField, elemSelector, nodesToConsider, qualityMetric, find_smoothed_node_location);

    if (stk::is_true_on_all_procs(mesh.parallel(), smoothInfos.empty()))
      break;

    communicate_smoothing_infos_that_other_procs_need_to_know_about(smoothInfos, mesh.parallel());

    smoothInfos = find_smooth_info_independent_sets(smoothInfos, mesh.parallel());

    if (stk::is_true_on_all_procs(mesh.parallel(), smoothInfos.empty()))
      break;

    smooth_nodes(mesh, coordsField, smoothInfos);

    nodesToConsider = get_nodes_to_consider_smoothing_for_next_iteration(mesh, elemSelector, ownedNodeSelector, smoothInfos);

    if (i%10 == 0)
      krinolog << "After " << i << " iterations, smoothing quality is " << compute_mesh_quality(mesh, elemSelector, coordsField, qualityMetric) << stk::diag::dendl;
  }

  krinolog << "After smoothing quality is " << compute_mesh_quality(mesh, elemSelector, coordsField, qualityMetric) << stk::diag::dendl;
}

void improve_quality_by_ODT_smoothing_on_interior(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interiorNodeSelector = elemSelector & !auxMeta.exposed_boundary_part() & !auxMeta.block_boundary_part();
  const MeanRatioQualityMetric qualityMetric;
  improve_quality_by_smoothing(mesh, coordsField, interiorNodeSelector, elemSelector, smoothed_node_location_using_ODT, qualityMetric, maxNumSmoothIterations);
}

FindSmoothedNodeLocation build_function_for_smoothed_node_location_using_optimized_mean_ratio_smoothing(const OptimizeNodeLocation & optimize_node_location,
    const NodeSearchDirectionFilter & node_search_direction_filter)
{
  auto fn = [&](const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const stk::mesh::Entity node,
      const stk::math::Vector3d & preSmoothedLocation)
  {
    stk::math::Vector3d postSmoothedLocation = preSmoothedLocation;
    smoothed_node_location_using_optimized_mean_ratio_smoothing(mesh, coordsField, elemSelector, node, optimize_node_location, node_search_direction_filter, postSmoothedLocation);
    return postSmoothedLocation;
  };
  return fn;
}

void improve_quality_by_optimized_mean_ratio_smoothing_on_interior(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeNodeLocation & optimize_node_location,
    const size_t maxNumSmoothIterations)
{
  static const NodeSearchDirectionFilter pass_through_filter = [](stk::mesh::Entity, stk::math::Vector3d &) {};
  const auto find_smoothed_node_location = build_function_for_smoothed_node_location_using_optimized_mean_ratio_smoothing(optimize_node_location, pass_through_filter);

  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interiorNodeSelector = elemSelector & !auxMeta.exposed_boundary_part() & !auxMeta.block_boundary_part();
  const MeanRatioQualityMetric qualityMetric;

  improve_quality_by_smoothing(mesh, coordsField, interiorNodeSelector, elemSelector, find_smoothed_node_location, qualityMetric, maxNumSmoothIterations);
}

void improve_quality_by_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeNodeLocation & optimize_node_location,
    const NodeSearchDirectionFilter & node_search_direction_filter,
    const size_t maxNumSmoothIterations)
{
  const auto find_smoothed_node_location = build_function_for_smoothed_node_location_using_optimized_mean_ratio_smoothing(optimize_node_location, node_search_direction_filter);

  const stk::mesh::Selector nodeSelector = elemSelector;
  const MeanRatioQualityMetric qualityMetric;

  improve_quality_by_smoothing(mesh, coordsField, nodeSelector, elemSelector, find_smoothed_node_location, qualityMetric, maxNumSmoothIterations);
}

void fill_tet_node_locations(const DistributedVector& nodeLocs, const std::array<unsigned,4> & elemNodeIndices, std::array<stk::math::Vector3d,4> & elemNodeLocs)
{
  for (unsigned i=0; i<4; ++i)
    for (unsigned j=0; j<3; ++j)
      elemNodeLocs[i][j] = nodeLocs[3*elemNodeIndices[i]+j];
}

const double qualityEpsilon = 0.01;

double global_objective_element_contribution(const double elemQuality)
{
  if (elemQuality >= qualityEpsilon)
    return 1./elemQuality;
  return 2./qualityEpsilon - elemQuality/(qualityEpsilon*qualityEpsilon);
}

double global_objective_sensitivity_wrt_quality(const double elemQuality)
{
  if (elemQuality >= qualityEpsilon)
    return -1./(elemQuality*elemQuality);
  return -1./(qualityEpsilon*qualityEpsilon);
}

double compute_global_tet_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const std::vector<std::array<unsigned,4>> & elementsNodeIndices,
    const DistributedVector& nodeLocs)
{
  MeanRatioQualityMetric qualityMetric;
  std::array<stk::math::Vector3d,4> elemNodeLocs;
  double obj = 0.;
  for (auto & elemNodeIndices : elementsNodeIndices)
  {
    fill_tet_node_locations(nodeLocs, elemNodeIndices, elemNodeLocs);
    const double quality = qualityMetric.tet_mean_ratio(elemNodeLocs);
    obj += global_objective_element_contribution(quality);
  }

  all_reduce_sum(mesh.parallel(), obj);
  return obj;
}

void fill_global_tet_mean_ratio_quality_objective_sensitivity(const stk::mesh::BulkData & mesh,
    const std::vector<std::array<unsigned,4>> & elementsNodeIndices,
    const DistributedVector& nodeLocs,
    DistributedVector & sens)
{
  MeanRatioQualityMetricWithSensitivities qualityMetric;
  std::array<stk::math::Vector3d,4> elemNodeLocs;
  sens.assign(nodeLocs.comm(), nodeLocs.size(), nodeLocs.local_size(), 0.);
  for (auto & elemNodeIndices : elementsNodeIndices)
  {
    fill_tet_node_locations(nodeLocs, elemNodeIndices, elemNodeLocs);
    const auto & [quality, elemSens] = MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities(elemNodeLocs);
    const double dObj_dQual = global_objective_sensitivity_wrt_quality(quality);
    for (unsigned i=0; i<4; ++i)
      for (unsigned j=0; j<3; ++j)
        sens[3*elemNodeIndices[i]+j] += dObj_dQual * elemSens[3*i+j];
  }
}

unsigned get_index_in_sorted_entities(const stk::mesh::BulkData &mesh,
  const std::vector<stk::mesh::Entity> & sortedEntities,
  const stk::mesh::Entity entity)
{
  auto iter = std::lower_bound(sortedEntities.begin(), sortedEntities.end(), entity, stk::mesh::EntityLess(mesh));
  STK_ThrowAssert(iter != sortedEntities.end() && *iter == entity);
  return std::distance(sortedEntities.begin(), iter);
}

static void pack_owned_node_sensitivities_for_sharers(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = 3;
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (unsigned i=0; i<ownedNodes.size(); ++i)
    {
      if (mesh.bucket(ownedNodes[i]).shared())
      {
        const double * nodeSens = sens.data() + dim*i;
        const stk::mesh::EntityId nodeId = mesh.identifier(ownedNodes[i]);
        mesh.comm_shared_procs(ownedNodes[i], nodeSharedProcs);
        for (int procId : nodeSharedProcs)
        {
          commSparse.send_buffer(procId).pack(nodeId);
          commSparse.send_buffer(procId).pack(nodeSens, dim);
        }
      }
    }
  });
}

static void unpack_shared_node_sensitivities(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & sortedOwnedNodes,
    const std::vector<stk::mesh::Entity> & sortedSharedUnownedNodes,
    DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = 3;
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(mesh.is_valid(node));
      const unsigned index = get_index_in_sorted_entities(mesh, sortedSharedUnownedNodes, node);
      double * nodeSens = sens.data() + dim*(sortedOwnedNodes.size()+index);
      commSparse.recv_buffer(procId).unpack(nodeSens, dim);
    }
  });
}

static void pack_shared_node_sensitivity_contributions_for_owners(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes,
    const DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = 3;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (unsigned i=0; i<sharedUnownedNodes.size(); ++i)
    {
      const double * nodeSens = sens.data() + dim*(ownedNodes.size()+i);
      const int owner = mesh.parallel_owner_rank(sharedUnownedNodes[i]);

      commSparse.send_buffer(owner).pack(mesh.identifier(sharedUnownedNodes[i]));
      commSparse.send_buffer(owner).pack(nodeSens, dim);
    }
  });
}

static void unpack_shared_node_sensitivity_contributions(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & sortedOwnedNodes,
    const std::vector<stk::mesh::Entity> & sortedSharedUnownedNodes,
    DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = 3;
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(mesh.is_valid(node));
      const unsigned index = get_index_in_sorted_entities(mesh, sortedOwnedNodes, node);
      std::array<double, 3> contrib;
      commSparse.recv_buffer(procId).unpack(contrib.data(), dim);
      for (unsigned i=0; i<dim; ++i)
        sens[dim*index + i] += contrib[i];
    }
  });
}

static void communicate_shared_node_sensitivity_contributions_to_owners(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes,
    DistributedVector & sens)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_shared_node_sensitivity_contributions_for_owners(mesh, ownedNodes, sharedUnownedNodes, sens, commSparse);
  unpack_shared_node_sensitivity_contributions(mesh, ownedNodes, sharedUnownedNodes, sens, commSparse);
}

static void communicate_owned_node_sensitivities_to_sharers(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes,
    DistributedVector & sens)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_owned_node_sensitivities_for_sharers(mesh, ownedNodes, sens, commSparse);
  unpack_shared_node_sensitivities(mesh, ownedNodes, sharedUnownedNodes, sens, commSparse);
}

void communicate_sensitivities(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes,
    DistributedVector & sens)
{
  if (mesh.parallel_size() > 1)
  {
    communicate_shared_node_sensitivity_contributions_to_owners(mesh, ownedNodes, sharedUnownedNodes, sens);
    communicate_owned_node_sensitivities_to_sharers(mesh, ownedNodes, sharedUnownedNodes, sens);
  }
}

MeshNodesObjFn build_function_for_global_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const std::vector<std::array<unsigned,4>> & elemNodeIndices)
{
  auto fn = [&mesh, &elemNodeIndices](const DistributedVector& nodeLocs)
  {
    return compute_global_tet_mean_ratio_quality_objective(mesh, elemNodeIndices, nodeLocs);
  };
  return fn;
}

void filter_node_sensitivity(const NodeSearchDirectionFilter & node_search_direction_filter, const stk::mesh::Entity node, double * sens)
{
  stk::math::Vector3d nodeDir(sens);
  node_search_direction_filter(node, nodeDir);
  for (unsigned i=0; i<3; ++i)
    sens[i] = nodeDir[i];
}

MeshNodesObjFnSens build_function_for_global_mean_ratio_quality_objective_sensitivity(const stk::mesh::BulkData & mesh,
    const std::vector<std::array<unsigned,4>> & elemNodeIndices,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes,
    const NodeSearchDirectionFilter & node_search_direction_filter)
{
  auto fn = [&](const DistributedVector& nodeLocs, DistributedVector & sens)
  {
    fill_global_tet_mean_ratio_quality_objective_sensitivity(mesh, elemNodeIndices, nodeLocs, sens);
    if (node_search_direction_filter)
    {
      const unsigned numOwnedNodes = ownedNodes.size();
      for (unsigned i=0; i<numOwnedNodes; ++i)
        filter_node_sensitivity(node_search_direction_filter, ownedNodes[i], &sens[3*i]);
      for (unsigned i=0; i<sharedUnownedNodes.size(); ++i)
        filter_node_sensitivity(node_search_direction_filter, sharedUnownedNodes[i], &sens[3*(i+numOwnedNodes)]);
    }
    communicate_sensitivities(mesh, ownedNodes, sharedUnownedNodes, sens);
  };
  return fn;
}


DistributedVector gather_node_coordinates(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & ownedNodes,
    const std::vector<stk::mesh::Entity> & sharedUnownedNodes)
{
  const unsigned numOwnedNodes = ownedNodes.size();
  const unsigned numNodes = numOwnedNodes + sharedUnownedNodes.size();
  DistributedVector nodesCoords(mesh.parallel(), 3*numNodes, 3*numOwnedNodes);
  for (unsigned i=0; i<numNodes; ++i)
  {
    const stk::mesh::Entity node = i<numOwnedNodes ? ownedNodes[i] : sharedUnownedNodes[i-numOwnedNodes];
    const double * nodeCoordData = get_field_data(mesh, coordsField, node);
    for (unsigned j=0; j<3; ++j)
      nodesCoords[3*i+j] = nodeCoordData[j];
  }
  return nodesCoords;
}

void set_node_coordinates(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & nodes,
    const DistributedVector & nodesCoords)
{
  for (unsigned i=0; i<nodes.size(); ++i)
  {
    double * nodeCoordData = get_field_data(mesh, coordsField, nodes[i]);
    for (unsigned j=0; j<3; ++j)
      nodeCoordData[j] = nodesCoords[3*i+j];
  }
  parallel_sync_fields(mesh, {&coordsField.field()});
}

unsigned get_index_of_node_in_sorted_owned_or_shared_nodes(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> & sortedOwnedNodes,
    const std::vector<stk::mesh::Entity> & sortedSharedUnownedNodes,
    const stk::mesh::Entity node)
{
  return (mesh.bucket(node).owned()) ?
      get_index_in_sorted_entities(mesh, sortedOwnedNodes, node) :
      (sortedOwnedNodes.size() + get_index_in_sorted_entities(mesh, sortedSharedUnownedNodes, node));
}

std::array<unsigned,4> get_node_indices_for_tet(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> & sortedOwnedNodes,
    const std::vector<stk::mesh::Entity> & sortedSharedUnownedNodes,
    const stk::mesh::Entity elem)
{
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  STK_ThrowAssert(elemNodes.size() == 4);
  std::array<unsigned,4> tetNodeIndices;
  for (unsigned i=0; i<4; ++i)
    tetNodeIndices[i] = get_index_of_node_in_sorted_owned_or_shared_nodes(mesh, sortedOwnedNodes, sortedSharedUnownedNodes, elemNodes[i]);
  return tetNodeIndices;
}

std::vector<std::array<unsigned,4>> get_element_node_indices(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> & elems,
    const std::vector<stk::mesh::Entity> & sortedOwnedNodes,
    const std::vector<stk::mesh::Entity> & sortedSharedUnownedNodes)
{
  std::vector<std::array<unsigned,4>> elemNodeIndices;
  elemNodeIndices.reserve(elems.size());
  for (auto & elem : elems)
    elemNodeIndices.push_back(get_node_indices_for_tet(mesh, sortedOwnedNodes, sortedSharedUnownedNodes, elem));
  return elemNodeIndices;
}

std::tuple<std::vector<stk::mesh::Entity>,std::vector<stk::mesh::Entity>> get_sorted_owned_nodes_and_unowned_shared_nodes(const stk::mesh::BulkData &mesh, const stk::mesh::Selector & elemSelector)
{
  const bool doSortById = true;
  const stk::mesh::Selector ownedSelector = elemSelector & mesh.mesh_meta_data().locally_owned_part();
  std::vector<stk::mesh::Entity> ownedNodes;
  stk::mesh::get_entities(mesh, stk::topology::NODE_RANK, ownedSelector, ownedNodes, doSortById);
  const stk::mesh::Selector sharedUnownedSelector = elemSelector & mesh.mesh_meta_data().globally_shared_part() & !mesh.mesh_meta_data().locally_owned_part();
  std::vector<stk::mesh::Entity> sharedUnownedNodes;
  stk::mesh::get_entities(mesh, stk::topology::NODE_RANK, sharedUnownedSelector, sharedUnownedNodes, doSortById);
  return std::make_tuple(ownedNodes, sharedUnownedNodes);
}

void improve_quality_by_simultaneous_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeMeshNodeLocations & optimize_node_locations,
    const NodeSearchDirectionFilter & node_search_direction_filter)
{
  MeanRatioQualityMetric qualityMetric;

  const auto & [ownedNodes, sharedUnownedNodes] = get_sorted_owned_nodes_and_unowned_shared_nodes(mesh, elemSelector);

  std::vector<stk::mesh::Entity> ownedElems;
  const stk::mesh::Selector ownedElementSelector = elemSelector & mesh.mesh_meta_data().locally_owned_part();
  stk::mesh::get_entities(mesh, stk::topology::ELEMENT_RANK, ownedElementSelector, ownedElems);
  const std::vector<std::array<unsigned,4>> elemNodeIndices = get_element_node_indices(mesh, ownedElems, ownedNodes, sharedUnownedNodes);

  DistributedVector nodeCoords = gather_node_coordinates(mesh, coordsField, ownedNodes, sharedUnownedNodes);

  const auto objective_at_loc = build_function_for_global_mean_ratio_quality_objective(mesh, elemNodeIndices);
  const auto sensitivity_at_loc = build_function_for_global_mean_ratio_quality_objective_sensitivity(mesh, elemNodeIndices, ownedNodes, sharedUnownedNodes, node_search_direction_filter);

  krinolog << "Before smoothing quality is " << compute_mesh_quality(mesh, ownedElementSelector, coordsField, qualityMetric) << stk::diag::dendl;

  optimize_node_locations(objective_at_loc, sensitivity_at_loc, nodeCoords);

  set_node_coordinates(mesh, coordsField, ownedNodes, nodeCoords);

  krinolog << "After smoothing quality is " << compute_mesh_quality(mesh, ownedElementSelector, coordsField, qualityMetric) << stk::diag::dendl;
}
}


