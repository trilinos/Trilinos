/*
 * Akri_Smooth.cpp
 *
 *  Created on: Nov 8, 2024
 *      Author: drnoble
 */
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_QualityMetricWithSensitivities.hpp>
#include <Akri_SmoothIndependentSetFinder.hpp>
#include <Akri_SmoothInfo.hpp>
#include <stk_mesh/base/BulkData.hpp>
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

std::pair<bool,stk::math::Vector3d> smoothed_node_location_using_conditional_ODT(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & preSmoothedLocation,
    const QualityMetric &qualityMetric)
{
  const double qualityBeforeSmoothing = calculate_quality_before_smoothing(mesh, coordsField, elemSelector, node, qualityMetric);
  const stk::math::Vector3d postSmoothedLocation = smoothed_node_location_using_ODT(mesh, coordsField, elemSelector, node, preSmoothedLocation);
  const double minSqrDisplacement = 0.; // Should be larger?

  if ((postSmoothedLocation-preSmoothedLocation).length_squared() > minSqrDisplacement)
  {
    const double qualityAfterSmoothing = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(mesh, elemSelector, coordsField, node, postSmoothedLocation, qualityMetric, qualityBeforeSmoothing);
    if (!qualityMetric.is_first_quality_metric_better_than_second(qualityBeforeSmoothing, qualityAfterSmoothing))
      return {true, postSmoothedLocation};
  }
  return {false, preSmoothedLocation};
}

stk::math::Vector3d compute_optimal_location(const std::function<double(const stk::math::Vector3d&)> & objective_at_loc,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d&)> & sensitivity_at_loc,
    const stk::math::Vector3d & startingLocation)
{
  const stk::math::Vector3d loc = bfgs(objective_at_loc, sensitivity_at_loc, startingLocation);
  return loc;
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

double compute_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLocation)
{
  double obj = 0.;
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  const double normalize = 1./nodeElems.size();
  for (auto nodeElem : nodeElems)
  {
    if (elemSelector(mesh.bucket(nodeElem)))
    {
      const double quality = mean_ratio_element_quality(mesh, coordsField, nodeElem, node, nodeLocation);
      obj += normalize*(1.-quality)*(1.-quality);
    }
  }
  return obj;
}

std::tuple<double, std::array<double,3>> compute_mean_ratio_quality_objective_and_sensitivity(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & nodeLocation)
{
  double obj = 0.;
  std::array<double,3> sens{0., 0., 0.};
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  const double normalize = 1./nodeElems.size();
  for (auto nodeElem : nodeElems)
  {
    if (elemSelector(mesh.bucket(nodeElem)))
    {
      const auto & [quality, sensitivity] = tet_mean_ratio_quality_and_sensitivity_to_nth_node(mesh, coordsField, nodeElem, get_entity_node_ordinal(mesh, nodeElem, node), nodeLocation);
      obj += normalize*(1.-quality)*(1.-quality);
      for (int d=0; d<3; ++d)
        sens[d] -= 2. * normalize * (1.-quality) * sensitivity[d];
    }
  }
  return std::make_tuple(obj, sens);
}

std::function<double(const stk::math::Vector3d&)> build_function_for_mean_ratio_quality_objective(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node)
{
  auto fn = [&mesh, coordsField, elemSelector, node](const stk::math::Vector3d & nodeLoc)
  {
    return compute_mean_ratio_quality_objective(mesh, coordsField, elemSelector, node, nodeLoc);
  };
  return fn;
}

std::function<stk::math::Vector3d(const stk::math::Vector3d&)> build_function_for_mean_ratio_quality_objective_sensitivity(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node)
{
  auto fn = [&mesh, coordsField, elemSelector, node](const stk::math::Vector3d & nodeLoc)
  {
    const auto & [origObj, sensitivity] = compute_mean_ratio_quality_objective_and_sensitivity(mesh, coordsField, elemSelector, node, nodeLoc);
    return stk::math::Vector3d(sensitivity[0], sensitivity[1], sensitivity[2]);
  };
  return fn;
}

std::pair<bool,stk::math::Vector3d> smoothed_node_location_using_optimized_mean_ratio_smoothing(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & preSmoothedLocation,
    const QualityMetric &/*qualityMetric*/)
{
  const auto objective_at_loc = build_function_for_mean_ratio_quality_objective(mesh, coordsField, elemSelector, node);
  const auto sensitivity_at_loc = build_function_for_mean_ratio_quality_objective_sensitivity(mesh, coordsField, elemSelector, node);

  const stk::math::Vector3d postSmoothedLocation = compute_optimal_location(objective_at_loc, sensitivity_at_loc, preSmoothedLocation);

  const double minSqrDisplacement = 0.; // Should be larger?

  if ((postSmoothedLocation-preSmoothedLocation).length_squared() > minSqrDisplacement)
    return {true, postSmoothedLocation};
  return {false, preSmoothedLocation};
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
    if (owndNodeSelector(mesh.bucket(smoothedNode)))
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

    smoothInfos.emplace_back(nodeId, owner, preSmoothedLocation, postSmoothedLocation, elementIds, procsThatNeedToKnowAboutInfo);
  });
}

void communicate_smoothing_infos_that_other_procs_need_to_know_about(std::vector<SmoothInfo> &smoothInfos, stk::ParallelMachine comm)
{
    stk::CommSparse commSparse(comm);

    pack_owned_smoothing_infos_that_other_procs_need_to_know_about(smoothInfos, commSparse);
    receive_smoothing_infos_that_this_proc_need_to_know_about(smoothInfos, commSparse);
}

using FindSmoothedNodeLocation =
    std::function<std::pair<bool,stk::math::Vector3d>(const stk::mesh::BulkData &, const CoordinatesFieldRef,const stk::mesh::Selector &, const stk::mesh::Entity, const stk::math::Vector3d &,const QualityMetric &)>;

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
    const auto & [doSmoothNode, postSmoothedLocation] = find_smoothed_node_location(mesh, coordsField, elemSelector, node, preSmoothedLocation, qualityMetric);

    if (doSmoothNode)
    {
      elementIds.clear();
      for (auto && element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
        if (elemSelector(mesh.bucket(element)))
          elementIds.push_back(mesh.identifier(element));

      fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutThisInfo);

      smoothInfos.emplace_back(mesh.identifier(node), mesh.parallel_rank(), preSmoothedLocation, postSmoothedLocation, elementIds, procsThatNeedToKnowAboutThisInfo);
    }
  }

  return smoothInfos;
}

void improve_quality_by_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
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

  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const stk::mesh::Selector ownedNodeSelector = mesh.mesh_meta_data().locally_owned_part() & auxMeta.active_part() & !auxMeta.exposed_boundary_part() & !auxMeta.block_boundary_part();

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
  }

  krinolog << "After smoothing quality is " << compute_mesh_quality(mesh, elemSelector, coordsField, qualityMetric) << stk::diag::dendl;
}

void improve_quality_by_ODT_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations)
{
  const MeanRatioQualityMetric qualityMetric;
  improve_quality_by_smoothing(mesh, coordsField, elemSelector, smoothed_node_location_using_conditional_ODT, qualityMetric, maxNumSmoothIterations);
}

void improve_quality_by_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations)
{
  const MeanRatioQualityMetric qualityMetric;
  improve_quality_by_smoothing(mesh, coordsField, elemSelector, smoothed_node_location_using_optimized_mean_ratio_smoothing, qualityMetric, maxNumSmoothIterations);
}

}


