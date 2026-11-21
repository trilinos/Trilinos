// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SharpFeatureLeafUtils.hpp>

#include <Akri_SharpFeatureIntersectionPoints.hpp>

namespace krino {

std::vector<double> simplex_weights_at_parametric_coordinates(
    const stk::math::Vector3d & parametricCoords, const unsigned numNodes)
{
  std::vector<double> weights(numNodes, 1.);
  for (unsigned i = 1; i < numNodes; ++i)
  {
    weights[i] = parametricCoords[i - 1];
    weights[0] -= parametricCoords[i - 1];
  }
  return weights;
}

void fill_entity_of_weights_above_tolerance(const std::vector<double> & weights,
    const double tol,
    std::vector<int> & entityOrdinals,
    std::vector<double> & entityWeights)
{
  entityOrdinals.clear();
  entityWeights.clear();
  double sumEntityWeights = 0.;
  for (size_t i = 0; i < weights.size(); ++i)
  {
    if (weights[i] > tol)
    {
      entityOrdinals.push_back(i);
      entityWeights.push_back(weights[i]);
      sumEntityWeights += weights[i];
    }
  }
  for (auto & weight : entityWeights)
    weight /= sumEntityWeights;
}

LocationOnTetFace locate_crossing_on_tet_face(const stk::math::Vector3d & faceCoords, const double edgeTol)
{
  const std::array<double, 3> faceWeights{{1. - faceCoords[0] - faceCoords[1], faceCoords[0], faceCoords[1]}};
  int numWeightsWithinTol = 0;
  int iWeightWithinTol = -1;
  for (int i = 0; i < 3; ++i)
  {
    if (faceWeights[i] <= edgeTol)
    {
      ++numWeightsWithinTol;
      iWeightWithinTol = i;
    }
  }
  if (numWeightsWithinTol == 0) return {TetFaceCrossingLocation::INSIDE, {{-1, -1}}};
  if (numWeightsWithinTol == 1)
    return {TetFaceCrossingLocation::NEAR_EDGE, {{(iWeightWithinTol + 1) % 3, (iWeightWithinTol + 2) % 3}}};
  return {TetFaceCrossingLocation::NEAR_NODE, {{-1, -1}}};
}

bool some_node_captures_point_domains(const std::vector<const std::vector<int> *> & nodesSortedDomains,
    const std::vector<int> & sortedPointDomains)
{
  for (const std::vector<int> * nodeDomains : nodesSortedDomains)
    if (nodeDomains->size() > sortedPointDomains.size() &&
        node_domains_capture_point_domains(*nodeDomains, sortedPointDomains))
      return true;
  return false;
}

bool tet_face_point_is_captured(const std::array<const std::vector<int> *, 3> & faceNodesSortedDomains,
    const stk::math::Vector3d & faceCoords,
    const std::vector<int> & sortedPointDomains,
    const double edgeTol)
{
  if (some_node_captures_point_domains(
          {faceNodesSortedDomains[0], faceNodesSortedDomains[1], faceNodesSortedDomains[2]}, sortedPointDomains))
    return true;
  const std::array<double, 3> faceWeights{{1. - faceCoords[0] - faceCoords[1], faceCoords[0], faceCoords[1]}};
  std::array<bool, 3> nodeHasDomains;
  for (int i = 0; i < 3; ++i)
    nodeHasDomains[i] = node_domains_capture_point_domains(*faceNodesSortedDomains[i], sortedPointDomains);
  for (int i = 0; i < 3; ++i)
  {
    if (nodeHasDomains[i] && faceWeights[i] >= 1. - edgeTol) return true;
    const int j = (i + 1) % 3;
    const int opposite = (i + 2) % 3;
    if (nodeHasDomains[i] && nodeHasDomains[j] && faceWeights[opposite] <= edgeTol) return true;
  }
  return false;
}

} // namespace krino
