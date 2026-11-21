// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MultipleEdgeCrossings.hpp>
#include <algorithm>

#include <Akri_SharpFeatureIntersectionPoints.hpp>

namespace krino {

// Remove the close crossings of an edge (positions along it): each two neighboring crossings closer
// than tol to each other, where the edge crosses the surface and back within the tolerance, as
// beside a sharp feature moved within the tolerance onto the edge.  Removing them in pairs keeps what
// the signs at the ends show.
static void remove_close_crossings(std::vector<double> & positions, const double tol)
{
  std::sort(positions.begin(), positions.end());
  std::vector<double> kept;
  for (size_t i = 0; i < positions.size(); ++i)
  {
    if (i + 1 < positions.size() && positions[i + 1] - positions[i] <= tol)
      ++i;
    else
      kept.push_back(positions[i]);
  }
  positions.swap(kept);
}

std::vector<double> multiple_crossing_positions_to_cut(const std::vector<double> & positionsFromNode0,
    const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const double closeCrossingTol,
    const double atEndTol)
{
  // The crossings are paired from the end with the lexicographically smaller coordinates, so that
  // every element with the edge, on any processor, removes the same close crossings.
  const bool reversed =
      std::lexicographical_compare(x1.data(), x1.data() + 3, x0.data(), x0.data() + 3);
  const double length = (x1 - x0).length();

  std::vector<double> positions;
  positions.reserve(positionsFromNode0.size());
  for (const double position : positionsFromNode0)
    positions.push_back(reversed ? 1. - position : position);
  remove_close_crossings(positions, closeCrossingTol / length);

  for (auto & position : positions)
  {
    if (position * length <= atEndTol)
      position = 0.;
    else if ((1. - position) * length <= atEndTol)
      position = 1.;
    if (reversed) position = 1. - position;
  }
  std::sort(positions.begin(), positions.end());
  return positions;
}

double multiple_crossing_end_tolerance(const bool surfaceHasSharpFeatures, const double edgeTol)
{
  constexpr double endTolForSurfacesWithoutSharpFeatures = 0.33;
  return surfaceHasSharpFeatures ? edgeTol : endTolForSurfacesWithoutSharpFeatures;
}

bool edge_is_chord_of_interface(const std::vector<int> & end0SortedDomains,
    const std::vector<int> & end1SortedDomains,
    const std::vector<int> & interfaceSortedDomains)
{
  return node_domains_capture_point_domains(end0SortedDomains, interfaceSortedDomains) &&
      node_domains_capture_point_domains(end1SortedDomains, interfaceSortedDomains);
}

std::vector<double> multiple_crossing_positions_to_insert(const std::vector<double> & positionsFromNode0,
    const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const std::vector<int> & end0SortedDomains,
    const std::vector<int> & end1SortedDomains,
    const std::vector<int> & interfaceSortedDomains,
    const double closeCrossingTol,
    const double atEndTol)
{
  std::vector<double> positions =
      multiple_crossing_positions_to_cut(positionsFromNode0, x0, x1, closeCrossingTol, atEndTol);
  std::erase_if(positions,
      [&](const double position)
      {
        return (position == 0. && node_domains_capture_point_domains(end0SortedDomains, interfaceSortedDomains)) ||
            (position == 1. && node_domains_capture_point_domains(end1SortedDomains, interfaceSortedDomains));
      });
  return positions;
}

}
