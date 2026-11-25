// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>
#include <iterator>
#include <tuple>
#include <utility>

#include <Akri_SharpFeatureIntersectionPoints.hpp>
#include <Akri_EntityOrder.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

template <typename NODE>
static void
reduce_stencil(std::vector<NODE> & nodes, std::vector<double> & weights, const double tol)
{
  STK_ThrowAssert(nodes.size() == weights.size() && !nodes.empty());
  const size_t iLargest =
      std::distance(weights.begin(), std::max_element(weights.begin(), weights.end()));
  size_t numKept = 0;
  double sumKept = 0.;
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (weights[i] >= tol || i == iLargest)
    {
      nodes[numKept] = nodes[i];
      weights[numKept] = weights[i];
      sumKept += weights[i];
      ++numKept;
    }
  }
  nodes.resize(numKept);
  weights.resize(numKept);
  for (auto & weight : weights)
    weight /= sumKept;
}

void reduce_intersection_point_stencil(
    std::vector<stk::mesh::Entity> & nodes, std::vector<double> & weights, const double tol)
{
  reduce_stencil(nodes, weights, tol);
}

void fill_reduced_simplex_vertex_weights(const stk::math::Vector3d & paramCoords,
    const unsigned dim,
    const double tol,
    std::vector<int> & vertexOrdinals,
    std::vector<double> & weights)
{
  STK_ThrowAssert(dim == 2 || dim == 3);
  if (dim == 3)
    weights = {1. - paramCoords[0] - paramCoords[1] - paramCoords[2],
        paramCoords[0],
        paramCoords[1],
        paramCoords[2]};
  else
    weights = {1. - paramCoords[0] - paramCoords[1], paramCoords[0], paramCoords[1]};
  vertexOrdinals.resize(weights.size());
  for (size_t i = 0; i < weights.size(); ++i)
    vertexOrdinals[i] = i;
  reduce_stencil(vertexOrdinals, weights, tol);
}

bool node_domains_capture_point_domains(
    const std::vector<int> & sortedNodeDomains, const std::vector<int> & sortedPointDomains)
{
  return std::includes(sortedNodeDomains.begin(),
      sortedNodeDomains.end(),
      sortedPointDomains.begin(),
      sortedPointDomains.end());
}

int sharp_feature_codimension(const std::vector<int> & sortedDomains)
{
  STK_ThrowAssert(!sortedDomains.empty());
  const bool areAllDomainsTheSame = sortedDomains.front() == sortedDomains.back();
  return areAllDomainsTheSame ? int(sortedDomains.size()) : int(sortedDomains.size()) - 1;
}

static std::vector<std::pair<stk::mesh::Entity, double>>
nodes_and_weights_sorted_by_node(const IntersectionPoint & intersectionPoint)
{
  std::vector<std::pair<stk::mesh::Entity, double>> nodesAndWeights;
  for (size_t i = 0; i < intersectionPoint.get_nodes().size(); ++i)
    nodesAndWeights.emplace_back(
        intersectionPoint.get_nodes()[i], intersectionPoint.get_weights()[i]);
  std::sort(nodesAndWeights.begin(), nodesAndWeights.end());
  return nodesAndWeights;
}

bool is_same_intersection_point(
    const IntersectionPoint & a, const IntersectionPoint & b, const double tol)
{
  if (a.get_nodes().size() != b.get_nodes().size() ||
      a.get_sorted_domains() != b.get_sorted_domains())
    return false;
  const auto aNodesAndWeights = nodes_and_weights_sorted_by_node(a);
  const auto bNodesAndWeights = nodes_and_weights_sorted_by_node(b);
  for (size_t i = 0; i < aNodesAndWeights.size(); ++i)
    if (aNodesAndWeights[i].first != bNodesAndWeights[i].first ||
        std::abs(aNodesAndWeights[i].second - bNodesAndWeights[i].second) > tol)
      return false;
  return true;
}

bool is_captured_by_intersection_point(
    const IntersectionPoint & capturing, const IntersectionPoint & point, const double tol)
{
  const std::vector<int> & capturingDomains = capturing.get_sorted_domains();
  const std::vector<int> & pointDomains = point.get_sorted_domains();
  if (capturing.get_nodes().size() != point.get_nodes().size() ||
      !node_domains_capture_point_domains(capturingDomains, pointDomains))
    return false;
  const auto capturingNodesAndWeights = nodes_and_weights_sorted_by_node(capturing);
  const auto pointNodesAndWeights = nodes_and_weights_sorted_by_node(point);
  for (size_t i = 0; i < capturingNodesAndWeights.size(); ++i)
    if (capturingNodesAndWeights[i].first != pointNodesAndWeights[i].first ||
        std::abs(capturingNodesAndWeights[i].second - pointNodesAndWeights[i].second) > tol)
      return false;
  return true;
}

IntersectionPoint reduced_intersection_point_in_canonical_form(
    const stk::mesh::BulkData & mesh, const IntersectionPoint & intersectionPoint, const double tol)
{
  std::vector<stk::mesh::Entity> nodes = intersectionPoint.get_nodes();
  std::vector<double> weights = intersectionPoint.get_weights();
  reduce_intersection_point_stencil(nodes, weights, tol);

  std::vector<std::pair<stk::mesh::EntityId, size_t>> idsAndIndices;
  for (size_t i=0; i<nodes.size(); ++i)
    idsAndIndices.emplace_back(mesh.identifier(nodes[i]), i);
  // The canonical form only represents the point, so its stencil is in id order.
  std::sort(idsAndIndices.begin(),
      idsAndIndices.end(),
      [](const auto & a, const auto & b) { return less_by_global_id_because_it_doesnt_impact_structure(a.first, b.first); });
  std::vector<stk::mesh::Entity> sortedNodes;
  std::vector<double> sortedWeights;
  for (auto & [id, i] : idsAndIndices)
  {
    sortedNodes.push_back(nodes[i]);
    sortedWeights.push_back(weights[i]);
  }
  return IntersectionPoint(intersectionPoint.is_owned(), sortedNodes, sortedWeights, intersectionPoint.get_sorted_domains());
}

namespace {

// The point's stencil nodes, by their model coordinates (compared in single precision, so that the
// order is the same on every build and for any parallel decomposition), with their weights.
std::vector<std::pair<stk::math::Vector3d, double>> stencil_coordinates_and_weights_in_coordinate_order(
    const stk::mesh::BulkData & mesh, const FieldRef coordsField, const IntersectionPoint & intersectionPoint)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<std::pair<stk::math::Vector3d, double>> coordsAndWeights;
  for (size_t i = 0; i < intersectionPoint.get_nodes().size(); ++i)
    coordsAndWeights.emplace_back(get_vector_field(mesh, coordsField, intersectionPoint.get_nodes()[i], dim),
        intersectionPoint.get_weights()[i]);
  std::stable_sort(coordsAndWeights.begin(), coordsAndWeights.end(),
      [](const auto & a, const auto & b) { return less_by_coordinates_because_it_impacts_structure(a.first, b.first); });
  return coordsAndWeights;
}

// Interior points, then face, edge and single-node points; each by its stencil nodes' coordinates,
// then more domains first, then domains, then weights, and only then node ids (never deciding but
// for points at the same nodes to single precision).  Node ids are not used before, because the ids
// that refinement gives the nodes depend on the parallel decomposition and on the build.
bool is_before_in_canonical_order(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const IntersectionPoint & a,
    const IntersectionPoint & b)
{
  if (a.get_nodes().size() != b.get_nodes().size())
    return a.get_nodes().size() > b.get_nodes().size();

  const auto coordsAndWeightsA = stencil_coordinates_and_weights_in_coordinate_order(mesh, coordsField, a);
  const auto coordsAndWeightsB = stencil_coordinates_and_weights_in_coordinate_order(mesh, coordsField, b);
  for (size_t i = 0; i < coordsAndWeightsA.size(); ++i)
  {
    if (less_by_coordinates_because_it_impacts_structure(coordsAndWeightsA[i].first, coordsAndWeightsB[i].first)) return true;
    if (less_by_coordinates_because_it_impacts_structure(coordsAndWeightsB[i].first, coordsAndWeightsA[i].first)) return false;
  }

  if (a.get_sorted_domains().size() != b.get_sorted_domains().size())
    return a.get_sorted_domains().size() > b.get_sorted_domains().size();
  if (a.get_sorted_domains() != b.get_sorted_domains())
    return a.get_sorted_domains() < b.get_sorted_domains();

  for (size_t i = 0; i < coordsAndWeightsA.size(); ++i)
    if (coordsAndWeightsA[i].second != coordsAndWeightsB[i].second)
      return coordsAndWeightsA[i].second < coordsAndWeightsB[i].second;

  std::vector<stk::mesh::EntityId> idsA, idsB;
  for (auto node : a.get_nodes()) idsA.push_back(mesh.identifier(node));
  for (auto node : b.get_nodes()) idsB.push_back(mesh.identifier(node));
  // Only for points the same in all of the above: repeats, so the order changes nothing.
  std::sort(idsA.begin(), idsA.end());
  std::sort(idsB.begin(), idsB.end());
  return less_by_global_id_because_it_doesnt_impact_structure(idsA, idsB);
}

} // namespace

void sort_intersection_points_in_canonical_order(
    const stk::mesh::BulkData & mesh, std::vector<IntersectionPoint> & intersectionPoints)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::stable_sort(intersectionPoints.begin(),
      intersectionPoints.end(),
      [&mesh, coordsField](const IntersectionPoint & a, const IntersectionPoint & b)
      { return is_before_in_canonical_order(mesh, coordsField, a, b); });
}

void remove_repeated_intersection_points(std::vector<IntersectionPoint> & intersectionPoints, const double tol)
{
  std::vector<IntersectionPoint> uniquePoints;
  for (auto & intersectionPoint : intersectionPoints)
  {
    const bool isRepeat = std::any_of(uniquePoints.begin(),
        uniquePoints.end(),
        [&](const IntersectionPoint & unique)
        { return is_captured_by_intersection_point(unique, intersectionPoint, tol); });
    if (!isRepeat)
      uniquePoints.push_back(intersectionPoint);
  }
  intersectionPoints.swap(uniquePoints);
}

std::vector<IntersectionPoint>
reduced_intersection_points_of_codimension(const stk::mesh::BulkData & mesh,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const int codimension,
    const double tol)
{
  std::vector<IntersectionPoint> points;
  for (auto & intersectionPoint : intersectionPoints)
    if (codimension == sharp_feature_codimension(intersectionPoint.get_sorted_domains()))
      points.push_back(reduced_intersection_point_in_canonical_form(mesh, intersectionPoint, tol));
  sort_intersection_points_in_canonical_order(mesh, points);
  // Only the same point found more than once (as on two faces that reduce to the same edge) is a
  // repeat: distinct points near each other are both kept, and whether one is already captured is
  // decided by the domains of the nodes when cutting.
  constexpr double repeatedPointWeightTol = 1.e-12;
  remove_repeated_intersection_points(points, repeatedPointWeightTol);
  return points;
}

} // namespace krino
