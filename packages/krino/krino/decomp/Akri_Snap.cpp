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
#include <Akri_CutQualityEstimator.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_MeshQuality.hpp>
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
        if (false == whichSnapsAreAllowed[iNode] && krinolog.shouldPrint(LOG_DEBUG))
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
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const std::vector<stk::mesh::Entity> & intPtNodes,
    const std::vector<double> & intPtWts)
{
  // A node may snap only to the points in its sector: where its weight exceeds rho times the weight
  // of every other node of the point, rho = (1 - maxEdgeSnap)/maxEdgeSnap.  On an edge, that is a
  // snap of at most maxEdgeSnap of the edge.  On a face or in an element it limits the snap the same
  // way near the edges, and changes continuously as a point moves off an edge or face; for
  // maxEdgeSnap = 0.5, the sectors of the nodes fill the edge, face or element.  The node's weight
  // must also exceed minWeightAtNode: with less, the point is within the edge tolerance of the entity
  // of the other nodes, where cutting would put it.
  bool anySnapsAreAllowedBasedOnWeight = true;
  std::vector<bool> whichSnapsAreAllowed(intPtNodes.size(), true);
  if (maxEdgeSnap < 1.0 || minWeightAtNode > 0.)
  {
    anySnapsAreAllowedBasedOnWeight = false;
    for (size_t iNode = 0; iNode < intPtNodes.size(); ++iNode)
    {
      double maxWeightAtOtherNodes = 0.;
      for (size_t jNode = 0; jNode < intPtNodes.size(); ++jNode)
        if (jNode != iNode)
          maxWeightAtOtherNodes = std::max(maxWeightAtOtherNodes, intPtWts[jNode]);
      const bool isInSector = maxEdgeSnap > 0. &&
          maxEdgeSnap * intPtWts[iNode] > (1. - maxEdgeSnap) * maxWeightAtOtherNodes;
      whichSnapsAreAllowed[iNode] = isInSector && intPtWts[iNode] > minWeightAtNode;
      anySnapsAreAllowedBasedOnWeight =
          anySnapsAreAllowedBasedOnWeight || whichSnapsAreAllowed[iNode];
    }
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
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const std::vector<IntersectionPoint> & intersectionPoints)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());

  mapFromEntityToIntPtIndexAndSnapAllowed nodeToIntPtIndicesAndWhichSnapsAllowed;
  for (size_t intersectionPointIndex=0; intersectionPointIndex<intersectionPoints.size(); ++intersectionPointIndex)
  {
    const auto & intPtNodes = intersectionPoints[intersectionPointIndex].get_nodes();
    const auto & intPtWts = intersectionPoints[intersectionPointIndex].get_weights();
    const std::vector<bool> whichSnapsAreAllowed = which_intersection_point_nodes_are_allowed_for_snapping(mesh, auxMeta, phaseSupport, sharpFeatureInfo, maxEdgeSnap, minWeightAtNode, intPtNodes, intPtWts);
    for (size_t iNode=0; iNode<intPtNodes.size(); ++iNode)
    {
      stk::mesh::Entity node = intPtNodes[iNode];
      if (mesh.bucket(node).owned())
        nodeToIntPtIndicesAndWhichSnapsAllowed[node].emplace_back(intersectionPointIndex,whichSnapsAreAllowed[iNode]);
    }
  }
  return nodeToIntPtIndicesAndWhichSnapsAllowed;
}

// For an edge intersection point, the node moves along the edge as the level set increases, at the
// intersection point's speed.  Other intersection points get no rate or direction, so they are
// ordered by distance, then location.
static void
set_snap_distance_rate_and_direction_for_predictable_selection(const stk::mesh::BulkData & mesh,
    const InterfaceGeometry & geometry,
    const stk::mesh::Entity node,
    const IntersectionPoint & intersectionPoint,
    SnapInfo & snapInfo)
{
  const auto & intPtNodes = intersectionPoint.get_nodes();
  if (intPtNodes.size() != 2) return;
  const stk::mesh::Entity otherNode = (intPtNodes[0] == node) ? intPtNodes[1] : intPtNodes[0];
  const stk::math::Vector3d velocity =
      geometry.edge_intersection_point_velocity(mesh, node, otherNode);
  const double speed = velocity.length();
  if (speed > 0.) snapInfo.set_snap_distance_rate_and_direction(speed, velocity / speed);
}

static void append_snap_infos_from_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed,
    const QualityMetric & qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const bool globalIDsAreParallelConsistent,
    const InterfaceGeometry * geometryForPredictableSelection,
    std::vector<SnapInfo> & snapInfos)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<int> procsThatNeedToKnowAboutThisInfo;
  std::vector<size_t> globalIdsOfSnapNodeElems;

  int owner = mesh.parallel_rank();

  const CutQualityEstimator cutQualityEstimator(mesh,
      elementSelector,
      coordsField,
      qualityMetric,
      minIntPtWeightForEstimatingCutQuality,
      globalIDsAreParallelConsistent);
  const auto domainsToNodesToQuality = cutQualityEstimator.estimate_quality_per_node_per_domain(
      intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed);
  const double minQualityThatIsForSureNotInverted = 0.;

  auto coordsData = coordsField.data<double>();
  for (auto entry : nodeToIntPtIndicesAndWhichSnapsAllowed)
  {
    stk::mesh::Entity node = entry.first;
    const auto nodeIntersectionPointIndicesAndWhichSnapsAllowed = entry.second;

    if (mesh.bucket(node).owned())
    {
      const stk::math::Vector3d nodeLocation(coordsData.entity_values(node), dim);

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

          const double postSnapQuality = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(mesh, elementSelector, coordsField, node, snapLocation, qualityMetric, minAcceptableQuality);
          if (qualityMetric.is_first_quality_metric_better_than_second(postSnapQuality, minAcceptableQuality))
          {
            const size_t nodeGlobalId = mesh.identifier(node);

            fill_global_ids_of_elements_using_node(mesh, elementSelector, node, globalIdsOfSnapNodeElems);
            fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutThisInfo);

            snapInfos.emplace_back(nodeGlobalId, intPtIndex, nodeLocation, owner, procsThatNeedToKnowAboutThisInfo, globalIdsOfSnapNodeElems, postSnapQuality, snapLocation, intPtNodes.size());
            if (nullptr != geometryForPredictableSelection)
              set_snap_distance_rate_and_direction_for_predictable_selection(mesh,
                  *geometryForPredictableSelection,
                  node,
                  intersectionPoint,
                  snapInfos.back());
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

std::vector<SnapInfo> build_snap_infos_from_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const SharpFeatureInfo * sharpFeatureInfo,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const QualityMetric & qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const bool globalIDsAreParallelConsistent,
    const InterfaceGeometry * geometryForPredictableSelection)
{
  std::vector<SnapInfo> snapInfos;

  const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed(mesh, sharpFeatureInfo, maxEdgeSnap, minWeightAtNode, intersectionPoints);
  append_snap_infos_from_intersection_points(mesh,
      coordsField,
      elementSelector,
      nodesToCapturedDomains,
      intersectionPoints,
      nodeToIntPtIndicesAndWhichSnapsAllowed,
      qualityMetric,
      minIntPtWeightForEstimatingCutQuality,
      globalIDsAreParallelConsistent,
      geometryForPredictableSelection,
      snapInfos);

  return snapInfos;
}

void interpolate_nodal_field(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node, const FieldRef field,
    const std::vector<stk::mesh::Entity> & interpNodes,
    const std::vector<double> & interpWeights,
    std::vector<double> & scratch)
{
  const unsigned fieldLength = field.length();

  if (!has_field_data(field, node)) return;

  scratch.assign(fieldLength, 0.0);

  auto fieldReadData = field.data<double>();
  for (size_t iNode=0; iNode<interpNodes.size(); ++iNode)
  {
    if (!has_field_data(field, interpNodes[iNode]))
    {
      krinolog << "When snapping node " << mesh.identifier(node) << ", the field " << field.name() << " is missing on interpolating node " << mesh.identifier(interpNodes[iNode]) << stk::diag::dendl;
      krinolog << "Should the field " << field.name() << " be an interpolation field?" << stk::diag::dendl;
      STK_ThrowRequireMsg(false, "Interpolation field missing on interpolation node " << mesh.identifier(interpNodes[iNode]));
    }

    auto nodeVal = fieldReadData.entity_values(interpNodes[iNode]);
    for (unsigned i=0; i<fieldLength; ++i)
      scratch[i] += interpWeights[iNode] * nodeVal(stk::mesh::ScalarIdx(i));
  }

  auto val = field.data<double, stk::mesh::ReadWrite>().entity_values(node);
  for (unsigned i=0; i<fieldLength; ++i)
    val(stk::mesh::ScalarIdx(i)) = scratch[i];
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

      for(auto && field : interpolationFieldSet)
        interpolate_nodal_field(mesh, snapNode, field, nodes, weights, scratch);

      if (krinolog.shouldPrint(LOG_DEBUG))
      {
        krinolog << "Snapping node " << snapInfo.get_node_global_id() << " to " << snapInfo.get_snap_location() << " at " << debug_output(mesh, intersectionPoint) << ", quality " << snapInfo.get_post_worst_quality() << stk::diag::dendl;
      }
    }
  }

  stk::mesh::communicate_field_data(mesh, interpFieldVec);

  communicate_node_captured_domains_for_given_nodes(mesh, snapNodes, nodesToCapturedDomains);
}

static double interpolate_nodal_field_component(const stk::mesh::BulkData & mesh, const FieldRef field, const unsigned component, const stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & interpNodes, const std::vector<double> & interpWeights)
{
  constexpr double wtTol = 1.e-11;
  double interpVal = 0.;
  if (has_field_data(field, node))
  {
    auto fieldReadData = field.data<double>();
    for (size_t iInterpNode=0; iInterpNode<interpNodes.size(); ++iInterpNode)
    {
      if (!has_field_data(field, interpNodes[iInterpNode]))
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
        auto nodeVal = fieldReadData.entity_values(interpNodes[iInterpNode]);
        interpVal += interpWeights[iInterpNode] * nodeVal(stk::mesh::ScalarIdx(component));
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

  auto fieldData = field.data<double, stk::mesh::ReadWrite>();
  for (size_t iSnapNode=0; iSnapNode<snapNodes.size(); ++iSnapNode)
  {
    if (has_field_data(field, snapNodes[iSnapNode]))
    {
      fieldData.entity_values(snapNodes[iSnapNode])(stk::mesh::ScalarIdx(component)) = scratch[iSnapNode];
    }
  }
}

static void fill_interpolation_nodes_in_element_at_parametric_coords(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity containingElem,
    const stk::math::Vector3d & containingElementParametricCoords,
    std::vector<stk::mesh::Entity> & interpNodes,
    std::vector<double> & interpWeights)
{
  const ParametricElement & parametricElem = AuxMetaData::get(mesh.mesh_meta_data()).parametric_element_library().get_parametric_element(mesh.bucket(containingElem).topology());

  interpNodes.assign(mesh.begin_nodes(containingElem), mesh.end_nodes(containingElem));
  interpWeights.assign(interpNodes.size(), 0.);
  parametricElem.shape_fcn(1, containingElementParametricCoords.data(), interpWeights.data());
}

void fill_interpolation_nodes_and_weights_for_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const stk::mesh::Entity node,
    const std::function<stk::math::Vector3d(const stk::mesh::BulkData &, const stk::mesh::Entity node)> & get_node_location,
    const std::function<void(const stk::mesh::BulkData &, const stk::mesh::Entity, std::vector<stk::math::Vector3d> &)> & fill_element_node_locations,
    std::vector<stk::mesh::Entity> & interpNodes,
    std::vector<double> & interpWeights)
{
  stk::mesh::Entity containingElem;
  stk::math::Vector3d containingElementParametricCoords;
  std::vector<stk::math::Vector3d> elementNodeCoords;

  const stk::math::Vector3d nodeLocation = get_node_location(mesh, node);

  double minSqrDist = std::numeric_limits<double>::max();
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (mesh.bucket(elem).member(activePart))
    {
      fill_element_node_locations(mesh, elem, elementNodeCoords);
      const stk::math::Vector3d elemParamCoords = get_parametric_coordinates_of_point(elementNodeCoords, nodeLocation);
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

std::function<stk::math::Vector3d(const stk::mesh::BulkData &, const stk::mesh::Entity node)>
build_get_unsnapped_node_location(const unsigned dim, const FieldRef coordsField, const FieldRef snapDisplacements)
{
  auto fn = [=](const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
  {
    return get_vector_field(mesh, coordsField, node, dim) - get_vector_field(mesh, snapDisplacements, node, dim);
  };
  return fn;
}

std::function<void(const stk::mesh::BulkData &, const stk::mesh::Entity, std::vector<stk::math::Vector3d> &)>
build_fill_current_element_node_locations(const unsigned dim, const FieldRef coordsField)
{
  auto fn = [=](const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem, std::vector<stk::math::Vector3d> & elemNodeCoords)
  {
    const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
    elemNodeCoords.resize(elemNodes.size());

    for (size_t i=0; i<elemNodes.size(); ++i)
      fill_vector_from_field(mesh, coordsField, elemNodes[i], dim, elemNodeCoords[i]);
  };
  return fn;
}

static std::vector<InterpolationPoint> build_interpolation_points_for_unsnapping(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField,  const FieldRef cdfemSnapField, const std::vector<stk::mesh::Entity> & snapNodes)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> interpNodes;
  std::vector<double> interpWeights;

  std::vector<InterpolationPoint> interpolationPoints;
  interpolationPoints.reserve(snapNodes.size());

  const auto & get_node_location = build_get_unsnapped_node_location(dim, coordsField, oldSnapDisplacements);
  const auto & fill_element_node_locations = build_fill_current_element_node_locations(dim, coordsField);

  for (auto && node : snapNodes)
  {
    fill_interpolation_nodes_and_weights_for_node(mesh, activePart, node, get_node_location, fill_element_node_locations, interpNodes, interpWeights);
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

std::vector<stk::mesh::Entity> get_owned_nodes_to_unsnap(const stk::mesh::BulkData & mesh, const stk::mesh::Part & /*activePart*/, const FieldRef /*coordsField*/, FieldRef cdfemSnapField)
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

std::vector<stk::mesh::Entity> get_owned_nodes_to_snap(const stk::mesh::BulkData & mesh, const stk::mesh::Part & /*activePart*/, const FieldRef /*coordsField*/, FieldRef cdfemSnapField)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> ownedSnapNodes;

  auto oldSnapData = oldSnapDisplacements.data<double>();
  auto newSnapData = cdfemSnapField.data<double>();
  stk::mesh::Selector ownedWithField = stk::mesh::selectField(cdfemSnapField) & mesh.mesh_meta_data().locally_owned_part();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedWithField))
  {
    for(const auto & node : *bucketPtr)
    {
      const stk::math::Vector3d oldSnap(oldSnapData.entity_values(node), dim);
      const stk::math::Vector3d newSnap(newSnapData.entity_values(node), dim);
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
                            commSparse.send_buffer(procId).pack<double>(
                                snapInfo.get_snap_distance_rate());
                            commSparse.send_buffer(procId).pack<stk::math::Vector3d>(
                                snapInfo.get_snap_direction());
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

        double snapDistanceRate{0.};
        commSparse.recv_buffer(procId).unpack<double>(snapDistanceRate);

        stk::math::Vector3d snapDirection;
        commSparse.recv_buffer(procId).unpack<stk::math::Vector3d>(snapDirection);

        snapInfos.emplace_back(globalNodeId,
            intersectionPointIndex,
            nodeLocation,
            procId,
            procsThatNeedToKnowAboutThisInfo,
            globalIdsOfSnapNodeElems,
            postSnapQuality,
            snapLocation,
            snapRank);
        snapInfos.back().set_snap_distance_rate_and_direction(snapDistanceRate, snapDirection);
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
    const double maxEdgeSnap,
    const double minWeightAtNode,
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
    const std::vector<bool> whichSnapsAreAllowed = which_intersection_point_nodes_are_allowed_for_snapping(mesh, auxMeta, phaseSupport, sharpFeatureInfo, maxEdgeSnap, minWeightAtNode, intPtNodes, intPtWts);
    for (size_t iNode=0; iNode<intPtNodes.size(); ++iNode)
    {
      stk::mesh::Entity intPtNode = intPtNodes[iNode];
      if (mesh.bucket(intPtNode).owned() && std::binary_search(sortedIdsOfNodesThatNeedNewSnapInfos.begin(), sortedIdsOfNodesThatNeedNewSnapInfos.end(), mesh.identifier(intPtNode)))
        nodeToIntPtIndicesAndWhichSnapsAllowed[intPtNode].emplace_back(intPtIndex,whichSnapsAreAllowed[iNode]);
    }
  }
  return nodeToIntPtIndicesAndWhichSnapsAllowed;
}

void update_intersection_points_and_snap_infos_after_snap_iteration(
    const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const InterfaceGeometry & geometry,
    const SharpFeatureInfo * sharpFeatureInfo,
    const std::vector<stk::mesh::Entity> & iterationSortedSnapNodes,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const stk::mesh::Selector & potentialParentElementSelector,
    const stk::mesh::Selector & decomposedParentElementSelector,
    const ScaledJacobianQualityMetric & qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const bool globalIDsAreParallelConsistent,
    const InterfaceGeometry * geometryForPredictableSelection,
    std::vector<IntersectionPoint> & intersectionPoints,
    std::vector<SnapInfo> & snapInfos)
{
  const std::vector<size_t> oldToNewIntPts = update_intersection_points_after_snap_iteration(mesh, decomposedParentElementSelector, geometry, iterationSortedSnapNodes, nodesToCapturedDomains, intersectionPoints);

  const std::vector<stk::mesh::EntityId> sortedIdsOfNodesThatNeedNewSnapInfos = get_sorted_ids_of_owned_nodes_of_elements_of_nodes(mesh, potentialParentElementSelector, iterationSortedSnapNodes);

  prune_snap_infos_modified_by_snap_iteration(mesh, oldToNewIntPts, sortedIdsOfNodesThatNeedNewSnapInfos, snapInfos);

  const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed_for_nodes_that_need_new_snap_infos(mesh, sharpFeatureInfo, maxEdgeSnap, minWeightAtNode, intersectionPoints, sortedIdsOfNodesThatNeedNewSnapInfos);

  append_snap_infos_from_intersection_points(mesh,
      coordsField,
      potentialParentElementSelector,
      nodesToCapturedDomains,
      intersectionPoints,
      nodeToIntPtIndicesAndWhichSnapsAllowed,
      qualityMetric,
      minIntPtWeightForEstimatingCutQuality,
      globalIDsAreParallelConsistent,
      geometryForPredictableSelection,
      snapInfos);
}

namespace {
class IntersectQuadraticElementsAsLinearWhileInScope
{
public:
  explicit IntersectQuadraticElementsAsLinearWhileInScope(const InterfaceGeometry & geometry) : myGeometry(geometry) { myGeometry.set_intersect_quadratic_elements_as_linear(true); }
  ~IntersectQuadraticElementsAsLinearWhileInScope() { myGeometry.set_intersect_quadratic_elements_as_linear(false); }
  IntersectQuadraticElementsAsLinearWhileInScope(const IntersectQuadraticElementsAsLinearWhileInScope &) = delete;
  IntersectQuadraticElementsAsLinearWhileInScope & operator=(const IntersectQuadraticElementsAsLinearWhileInScope &) = delete;
private:
  const InterfaceGeometry & myGeometry;
};
}

NodeToCapturedDomainsMap
snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & potentialParentElementSelector,
    const stk::mesh::Selector & decomposedParentElementSelector,
    const FieldRef coordsField,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent,
    const double snappingSharpFeatureAngleInDegrees,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const Snap_Selection snapSelection)
{/* %TRACE[ON]% */ Trace trace__("krino::snap_as_much_as_possible_while_maintaining_quality()"); /* %TRACE% */

    // Snapping moves vertices only, so quadratic elements are intersected as the linear elements of
    // their vertices; their mid-side nodes are made linear after snapping.
    const IntersectQuadraticElementsAsLinearWhileInScope intersectQuadraticElementsAsLinear(geometry);

    const ScaledJacobianQualityMetric qualityMetric;
    krinolog << "Before snapping quality is " << compute_mesh_quality(mesh, potentialParentElementSelector, coordsField, qualityMetric) << stk::diag::dendl;

    size_t iteration{0};
    NodeToCapturedDomainsMap nodesToCapturedDomains;
    stk::ParallelMachine comm = mesh.parallel();
    std::unique_ptr<SharpFeatureInfo> sharpFeatureInfo;
    if (snappingSharpFeatureAngleInDegrees > 0.)
    {
      sharpFeatureInfo = std::make_unique<SharpFeatureInfo>();
      sharpFeatureInfo->find_sharp_features(mesh,
          coordsField,
          potentialParentElementSelector,
          get_snapping_compatibility_side_selector(mesh),
          std::cos(snappingSharpFeatureAngleInDegrees * M_PI / 180.));
      if (krinolog.shouldPrint(LOG_DEBUG))
      {
        for (auto && [node, constraint] : sharpFeatureInfo->get_constraints())
        {
          krinolog << "Node " << mesh.identifier(node) << " is ";
          if (constraint.is_pinned())
          {
            krinolog << "pinned." << stk::diag::dendl;
          }
          else
          {
            const std::array<stk::mesh::Entity, 2> sharpEdgeNbrs =
                constraint.get_sharp_edge_nodes();
            krinolog << "constrained to move along edge between nodes "
                     << mesh.identifier(sharpEdgeNbrs[0]) << " and "
                     << mesh.identifier(sharpEdgeNbrs[1]) << "." << stk::diag::dendl;
          }
        }
      }
    }

    std::vector<IntersectionPoint> intersectionPoints;
    geometry.store_phase_for_uncut_elements(mesh);
    intersectionPoints = build_all_intersection_points(mesh, decomposedParentElementSelector, geometry, nodesToCapturedDomains);
    const bool isPredictable = snapSelection == Snap_Selection::PREDICTABLE;
    const InterfaceGeometry * geometryForPredictableSelection = isPredictable ? &geometry : nullptr;
    std::vector<SnapInfo> snapInfos = build_snap_infos_from_intersection_points(mesh,
        coordsField,
        sharpFeatureInfo.get(),
        potentialParentElementSelector,
        nodesToCapturedDomains,
        intersectionPoints,
        qualityMetric,
        minIntPtWeightForEstimatingCutQuality,
        maxEdgeSnap,
        minWeightAtNode,
        globalIDsAreParallelConsistent,
        geometryForPredictableSelection);

    while (true)
    {
      krinolog << "Snapping To Geometry Iteration " << std::to_string(++iteration) << stk::diag::dendl;

      bool done = stk::is_true_on_all_procs(comm, snapInfos.empty());
      if ( done )
          break;

      communicate_snap_infos_that_other_procs_need_to_know_about(snapInfos, comm);

      const std::vector<SnapInfo> independentSnapInfos =
          find_snap_info_independent_sets(snapInfos, qualityMetric, isPredictable, comm);

      krinolog << "  Snapping " << get_global_num_infos(independentSnapInfos, comm) << " of " << get_global_num_infos(snapInfos, comm) << " snap candidates." << stk::diag::dendl;

      geometry.store_phase_for_elements_that_will_be_uncut_after_snapping(mesh, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);
      snap_nodes(mesh, interpolationFields, intersectionPoints, independentSnapInfos, nodesToCapturedDomains);

      const std::vector<stk::mesh::Entity> iterationSortedSnapNodes = get_sorted_nodes_modified_in_current_snapping_iteration(mesh, independentSnapInfos);

      update_intersection_points_and_snap_infos_after_snap_iteration(mesh,
          coordsField,
          geometry,
          sharpFeatureInfo.get(),
          iterationSortedSnapNodes,
          nodesToCapturedDomains,
          potentialParentElementSelector,
          decomposedParentElementSelector,
          qualityMetric,
          minIntPtWeightForEstimatingCutQuality,
          maxEdgeSnap,
          minWeightAtNode,
          globalIDsAreParallelConsistent,
          geometryForPredictableSelection,
          intersectionPoints,
          snapInfos);
    }

    krinolog << "After snapping quality is " << compute_mesh_quality(mesh, potentialParentElementSelector, coordsField, qualityMetric) << stk::diag::dendl;

    return nodesToCapturedDomains;
}
}



