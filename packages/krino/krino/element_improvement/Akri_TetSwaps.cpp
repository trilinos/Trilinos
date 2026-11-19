// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>

#include <Akri_QualityMetric.hpp>
#include <Akri_TetSwaps.hpp>

namespace krino
{

double tet_volume(const std::array<stk::math::Vector3d, 4> & nodeCoords)
{
  return Dot(Cross(nodeCoords[1] - nodeCoords[0], nodeCoords[2] - nodeCoords[0]),
             nodeCoords[3] - nodeCoords[0]) /
      6.;
}

int tet_volume_sign(const std::array<stk::math::Vector3d, 4> & nodeCoords)
{
  double longestEdgeSqrLength = 0.;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      longestEdgeSqrLength =
          std::max(longestEdgeSqrLength, (nodeCoords[j] - nodeCoords[i]).length_squared());
  constexpr double relativeTol = 1.e-12;
  const double volumeTol = relativeTol * longestEdgeSqrLength * std::sqrt(longestEdgeSqrLength);
  const double volume = tet_volume(nodeCoords);
  return (volume > volumeTol) ? 1 : ((volume < -volumeTol) ? -1 : 0);
}

double tet_mean_ratio_quality(const std::array<stk::math::Vector3d, 4> & nodeCoords)
{
  static const MeanRatioQualityMetric qualityMetric;
  const std::vector<stk::math::Vector3d> nodeLocations(nodeCoords.begin(), nodeCoords.end());
  return qualityMetric.get_element_quality_metric(3, nodeLocations);
}

size_t num_tet_edge_swap_configurations(const size_t numRingNodes)
{
  if (numRingNodes >
      TetEdgeSwapConfigurations::get_max_num_nodes_around_edge_currently_supported_for_swapping())
    return 0;
  return TetEdgeSwapConfigurations::get_possible_tet_configurations_for_num_tets_around_edge(
      numRingNodes)
      .size();
}

} // namespace krino
