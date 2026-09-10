// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Objective_MeanRatioAboutNode.hpp>

#include <tuple>
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_QualityMetricWithSensitivities.hpp>
#include <Akri_MeshHelpers.hpp>

namespace krino {

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

double MeanRatioAboutNodeObjective::compute_value(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & x) const
{
  return compute_worst_mean_ratio_quality_objective(mesh, coordsField, elemSelector, node, x);
}

void MeanRatioAboutNodeObjective::fill_gradient(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const GradientFilter & gradientFilter,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & x,
    stk::math::Vector3d & g) const
{
  double obj = 0;
  std::tie(obj, g) = compute_worst_mean_ratio_quality_objective_and_sensitivity(mesh, coordsField, elemSelector, node, x);
  gradientFilter(node, x, g);
}

}

