/*
 * Akri_AdaptiveContourTet.cpp
 *
 *  Created on: Sep 30, 2025
 *      Author: drnoble
 */
#include <Akri_AdaptiveContourTet.hpp>

#include <Akri_AdaptiveContourUtils.hpp>
#include <Akri_ContourTet.hpp>
#include <Akri_ContourUtils.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_MOAB_TetRefiner.hpp>
#include <Akri_RefinerUtils.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <stk_topology/topology.hpp>

using krino::ContourTet;

namespace krino {

static constexpr double snapTol = 1.e-6;
static constexpr double nonlinearTol = 1.e-2;

std::array<double,10> apply_tet_snapping_and_clipping(const std::array<double,10> & unfilteredTet10Dist, const double snapDistTol)
{
  std::array<double,10> tet10Dist;
  for (int n=0; n<4; ++n)
    tet10Dist[n] = (std::abs(unfilteredTet10Dist[n]) < snapDistTol) ? 0.0 : unfilteredTet10Dist[n];
  tet10Dist[4] = clip_midedge_distance(tet10Dist[0], tet10Dist[1], unfilteredTet10Dist[4]);
  tet10Dist[5] = clip_midedge_distance(tet10Dist[1], tet10Dist[2], unfilteredTet10Dist[5]);
  tet10Dist[6] = clip_midedge_distance(tet10Dist[2], tet10Dist[0], unfilteredTet10Dist[6]);
  tet10Dist[7] = clip_midedge_distance(tet10Dist[0], tet10Dist[3], unfilteredTet10Dist[7]);
  tet10Dist[8] = clip_midedge_distance(tet10Dist[1], tet10Dist[3], unfilteredTet10Dist[8]);
  tet10Dist[9] = clip_midedge_distance(tet10Dist[2], tet10Dist[3], unfilteredTet10Dist[9]);
  return tet10Dist;
}

std::array<stk::math::Vector3d,10> get_child_tet10_coords(const std::array<stk::math::Vector3d,10>  & parentCoords,
    const std::array<int, 4> & childElemVertexIndices)
{
  std::array<stk::math::Vector3d,10> childTet10Coords =
  {{
      parentCoords[childElemVertexIndices[0]], parentCoords[childElemVertexIndices[1]], parentCoords[childElemVertexIndices[2]], parentCoords[childElemVertexIndices[3]],
      0.5*(parentCoords[childElemVertexIndices[0]]+parentCoords[childElemVertexIndices[1]]),
      0.5*(parentCoords[childElemVertexIndices[1]]+parentCoords[childElemVertexIndices[2]]),
      0.5*(parentCoords[childElemVertexIndices[2]]+parentCoords[childElemVertexIndices[0]]),
      0.5*(parentCoords[childElemVertexIndices[0]]+parentCoords[childElemVertexIndices[3]]),
      0.5*(parentCoords[childElemVertexIndices[1]]+parentCoords[childElemVertexIndices[3]]),
      0.5*(parentCoords[childElemVertexIndices[2]]+parentCoords[childElemVertexIndices[3]])
  }};
  return childTet10Coords;
}

std::array<stk::math::Vector3d,10> get_tet10_coords_on_tet4(const std::array<stk::math::Vector3d,4>  & tet4Coords)
{
  std::array<stk::math::Vector3d,10> tet10Coords =
  {{
      tet4Coords[0], tet4Coords[1], tet4Coords[2], tet4Coords[3],
      0.5*(tet4Coords[0]+tet4Coords[1]),
      0.5*(tet4Coords[1]+tet4Coords[2]),
      0.5*(tet4Coords[2]+tet4Coords[0]),
      0.5*(tet4Coords[0]+tet4Coords[3]),
      0.5*(tet4Coords[1]+tet4Coords[3]),
      0.5*(tet4Coords[2]+tet4Coords[3])
  }};
  return tet10Coords;
}

std::array<double,10> get_tet10_distance_on_tet4(const std::array<stk::math::Vector3d,10>  & tet10Coords,
    const std::array<double,4> & tet4Dist,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point)
{
  std::array<double,10> tet10Dist =
  {{
      tet4Dist[0], tet4Dist[1], tet4Dist[2], tet4Dist[3],
      distance_at_point(tet10Coords[4]),
      distance_at_point(tet10Coords[5]),
      distance_at_point(tet10Coords[6]),
      distance_at_point(tet10Coords[7]),
      distance_at_point(tet10Coords[8]),
      distance_at_point(tet10Coords[9])
  }};
  return tet10Dist;
}

std::array<double,25> interpolate_subtet(const std::array<double,10> & tet10Dist)
{
  std::array<double,25> subTetDist =
  {{
      0.375*tet10Dist[0] + 0.75*tet10Dist[4] - 0.125*tet10Dist[1], //0-4
      0.375*tet10Dist[1] + 0.75*tet10Dist[4] - 0.125*tet10Dist[0], //1-4
      0.375*tet10Dist[1] + 0.75*tet10Dist[5] - 0.125*tet10Dist[2], //1-5
      0.375*tet10Dist[2] + 0.75*tet10Dist[5] - 0.125*tet10Dist[1], //2-5
      0.375*tet10Dist[2] + 0.75*tet10Dist[6] - 0.125*tet10Dist[0], //2-6
      0.375*tet10Dist[0] + 0.75*tet10Dist[6] - 0.125*tet10Dist[2], //3-6
      0.5*tet10Dist[5] + 0.5*tet10Dist[6] - 0.125*tet10Dist[0] + 0.25*tet10Dist[4] - 0.125*tet10Dist[1], //5-6
      0.5*tet10Dist[4] + 0.5*tet10Dist[6] - 0.125*tet10Dist[1] + 0.25*tet10Dist[5] - 0.125*tet10Dist[2], //4-6
      0.5*tet10Dist[4] + 0.5*tet10Dist[5] - 0.125*tet10Dist[0] + 0.25*tet10Dist[6] - 0.125*tet10Dist[2], //4-5
      0.375*tet10Dist[0] + 0.75*tet10Dist[7] - 0.125*tet10Dist[3], //0-7
      0.5*tet10Dist[4] + 0.5*tet10Dist[7] - 0.125*tet10Dist[1] + 0.25*tet10Dist[8] - 0.125*tet10Dist[3], //4-7
      0.5*tet10Dist[4] + 0.5*tet10Dist[8] - 0.125*tet10Dist[0] + 0.25*tet10Dist[7] - 0.125*tet10Dist[3], //4-8
      0.375*tet10Dist[1] + 0.75*tet10Dist[8] - 0.125*tet10Dist[3], //1-8
      0.5*tet10Dist[5] + 0.5*tet10Dist[8] - 0.125*tet10Dist[2] + 0.25*tet10Dist[9] - 0.125*tet10Dist[3], //5-8
      0.5*tet10Dist[5] + 0.5*tet10Dist[9] - 0.125*tet10Dist[1] + 0.25*tet10Dist[8] - 0.125*tet10Dist[3], //5-9
      0.375*tet10Dist[2] + 0.75*tet10Dist[9] - 0.125*tet10Dist[3], //2-9
      0.5*tet10Dist[6] + 0.5*tet10Dist[9] - 0.125*tet10Dist[0] + 0.25*tet10Dist[7] - 0.125*tet10Dist[3], //6-9
      0.5*tet10Dist[6] + 0.5*tet10Dist[7] - 0.125*tet10Dist[2] + 0.25*tet10Dist[9] - 0.125*tet10Dist[3], //6-7
      0.5*tet10Dist[7] + 0.5*tet10Dist[8] - 0.125*tet10Dist[0] + 0.25*tet10Dist[4] - 0.125*tet10Dist[1], //7-8
      0.5*tet10Dist[8] + 0.5*tet10Dist[9] - 0.125*tet10Dist[1] + 0.25*tet10Dist[5] - 0.125*tet10Dist[2], //8-9
      0.5*tet10Dist[7] + 0.5*tet10Dist[9] - 0.125*tet10Dist[0] + 0.25*tet10Dist[6] - 0.125*tet10Dist[2], //7-9
      0.375*tet10Dist[3] + 0.75*tet10Dist[7] - 0.125*tet10Dist[0], //3-7
      0.375*tet10Dist[3] + 0.75*tet10Dist[8] - 0.125*tet10Dist[1], //3-8
      0.375*tet10Dist[3] + 0.75*tet10Dist[9] - 0.125*tet10Dist[2], //3-9
      -0.125*tet10Dist[0] - 0.125*tet10Dist[1] - 0.125*tet10Dist[2] - 0.125*tet10Dist[3] +
      0.25*tet10Dist[4] + 0.25*tet10Dist[5] + 0.25*tet10Dist[6] + 0.25*tet10Dist[7] + 0.25*tet10Dist[8] + 0.25*tet10Dist[9] // centroid
  }};
  return subTetDist;
}

std::array<double,25> interpolate_and_snap_subtet(const std::array<double,10> & tet10Dist, const double snapDistTol)
{
  std::array<double,25> subTetDist = interpolate_subtet(tet10Dist);
  for (int n=0; n<25; ++n)
    subTetDist[n] = (std::abs(subTetDist[n]) < snapDistTol) ? 0.0 : subTetDist[n];
  return subTetDist;
}

std::array<double,10> get_child_tet10_distance(const std::array<double,10> & tet10Dist,
  const std::array<double,25> & subTetDist,
  const std::array<int, 4> & subElemVertexIndices,
  const std::array<int, 6> & subElemMidsideIndices)
{
  const std::array<double,10> subTet10Dist = {{ tet10Dist[subElemVertexIndices[0]], tet10Dist[subElemVertexIndices[1]], tet10Dist[subElemVertexIndices[2]], tet10Dist[subElemVertexIndices[3]],
      subTetDist[subElemMidsideIndices[0]], subTetDist[subElemMidsideIndices[1]], subTetDist[subElemMidsideIndices[2]],
      subTetDist[subElemMidsideIndices[3]], subTetDist[subElemMidsideIndices[4]], subTetDist[subElemMidsideIndices[5]] }};
  return subTet10Dist;
}

void append_facets_for_subtet_of_converged_tet(const std::array<stk::math::Vector3d,10> & tet10Coords,
  const std::array<double,10> & tet10Dist,
  const std::array<double,25> & subTetDist,
  FacetedSurfaceBase & facets,
  const std::array<int, 4> & subElemVertexIndices,
  const std::array<int, 6> & subElemMidsideIndices)
{
  const std::array<double,10> subTet10Dist = get_child_tet10_distance(tet10Dist, subTetDist, subElemVertexIndices, subElemMidsideIndices);
  ContourTet::append_facets_for_tet4_with_tet10_distance(subarray(tet10Coords, subElemVertexIndices), subTet10Dist, facets);
}

int determine_tet_edge_refinement_case_id(const std::array<double,10> & tet10Dist,
  const double lengthScale,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  if (currentDepth < minDepth)
    return 63;
  if (currentDepth == maxDepth)
    return 0;
  const double nonlinearDistTol = nonlinearTol*lengthScale;
  int caseId = 0;
  if (!is_edge_converged(tet10Dist[0], tet10Dist[1], tet10Dist[4], nonlinearDistTol)) caseId += 1;
  if (!is_edge_converged(tet10Dist[1], tet10Dist[2], tet10Dist[5], nonlinearDistTol)) caseId += 2;
  if (!is_edge_converged(tet10Dist[2], tet10Dist[0], tet10Dist[6], nonlinearDistTol)) caseId += 4;
  if (!is_edge_converged(tet10Dist[0], tet10Dist[3], tet10Dist[7], nonlinearDistTol)) caseId += 8;
  if (!is_edge_converged(tet10Dist[1], tet10Dist[3], tet10Dist[8], nonlinearDistTol)) caseId += 16;
  if (!is_edge_converged(tet10Dist[2], tet10Dist[3], tet10Dist[9], nonlinearDistTol)) caseId += 32;
  return caseId;
}

void append_facets_for_refined_subtet_using_interpolated_distance(const std::array<stk::math::Vector3d,10> & tet10Coords,
  const std::array<stk::math::Vector3d,10> & tet10DepartureCoords,
  const std::array<double,10> & tet10Dist,
  const std::array<double,25> & subTetDist,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const std::array<int, 4> & subElemVertexIndices,
  const std::array<int, 6> & subElemMidsideIndices)
{
  if (!have_possibly_cut_edge(tet10DepartureCoords, tet10Dist, subElemVertexIndices))
    return;

  const std::array<stk::math::Vector3d,10> subTet10Coords = get_child_tet10_coords(tet10Coords, subElemVertexIndices);
  const std::array<double,10> subTet10Dist = get_child_tet10_distance(tet10Dist, subTetDist, subElemVertexIndices, subElemMidsideIndices);

  std::array<double,25> subSubTetDist = interpolate_subtet(subTet10Dist);
  snap_distance(subSubTetDist, snapTol*lengthScale);

  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{0,4,6,7}}, {{0,7,5,9,10,17}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{1,5,4,8}}, {{2,8,1,12,13,11}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{2,6,5,9}}, {{4,6,3,15,16,14}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{7,8,9,3}}, {{18,19,20,21,22,23}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{4,5,6,8}}, {{8,6,7,11,13,24}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{6,5,9,8}}, {{6,14,16,24,13,19}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{6,9,7,8}}, {{16,20,17,24,19,18}});
  append_facets_for_subtet_of_converged_tet(subTet10Coords, subTet10Dist, subSubTetDist, facets, {{6,7,4,8}}, {{17,10,7,24,18,11}});
}

void adaptively_append_facets_for_subtet_using_semilagrangian_distance(const std::array<stk::math::Vector3d,10> & parentCoords,
  const std::array<stk::math::Vector3d,10> & parentDepartureCoords,
  const std::array<double,10> & parentDistance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth,
  const std::array<int, 4> & childElemNodeIndices);

void append_facets_for_converged_tet10(const std::array<stk::math::Vector3d,10> & coords,
  const std::array<stk::math::Vector3d,10> & departureCoords,
  const std::array<double,10> & distance,
  const double lengthScale,
  FacetedSurfaceBase & facets)
{
  const std::array<double,10> filteredDistance = apply_tet_snapping_and_clipping(distance, snapTol*lengthScale);
  std::array<double,25> subTetDist = interpolate_subtet(filteredDistance);

  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{0,4,6,7}}, {{0,7,5,9,10,17}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{1,5,4,8}}, {{2,8,1,12,13,11}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{2,6,5,9}}, {{4,6,3,15,16,14}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{7,8,9,3}}, {{18,19,20,21,22,23}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{4,5,6,8}}, {{8,6,7,11,13,24}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{6,5,9,8}}, {{6,14,16,24,13,19}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{6,9,7,8}}, {{16,20,17,24,19,18}});
  append_facets_for_refined_subtet_using_interpolated_distance(coords, departureCoords, filteredDistance, subTetDist, lengthScale, facets, {{6,7,4,8}}, {{17,10,7,24,18,11}});
}

void adaptively_append_facets_for_tet10_using_semilagrangian_distance(const std::array<stk::math::Vector3d,10> & coords,
  const std::array<stk::math::Vector3d,10> & departureCoords,
  const std::array<double,10> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  const int refinementCaseId = determine_tet_edge_refinement_case_id(distance, lengthScale, currentDepth, minDepth, maxDepth);

  if (0 == refinementCaseId)
  {
    append_facets_for_converged_tet10(coords, departureCoords, distance, lengthScale, facets);
  }
  else
  {
    const LengthRatioMetricForTetRefinement qualityMetric;
    const std::array<int,4> parentNodeRank= get_rank_of_nodes_based_on_coordinates<4>(coords);
    const std::vector<std::array<int,4>> newTetNodes = moab::SimplexTemplateRefiner::refinement_child_nodes_tet4(qualityMetric, refinementCaseId, coords, parentNodeRank);

    for (auto & newTet : newTetNodes)
      adaptively_append_facets_for_subtet_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, newTet);
  }
}

void adaptively_append_facets_for_subtet_using_semilagrangian_distance(const std::array<stk::math::Vector3d,10> & parentCoords,
  const std::array<stk::math::Vector3d,10> & parentDepartureCoords,
  const std::array<double,10> & parentDistance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth,
  const std::array<int, 4> & childElemNodeIndices)
{
  if (!have_possibly_cut_edge(parentDepartureCoords, parentDistance, childElemNodeIndices))
    return;

  const std::array<stk::math::Vector3d,10> tet10Coords = get_child_tet10_coords(parentCoords, childElemNodeIndices);
  const std::array<stk::math::Vector3d,10> tet10DepartureCoords = get_child_tet10_coords(parentDepartureCoords, childElemNodeIndices);
  const std::array<double,10> tet10Dist = get_tet10_distance_on_tet4(tet10DepartureCoords, subarray(parentDistance, childElemNodeIndices), distance_at_point);

  adaptively_append_facets_for_tet10_using_semilagrangian_distance(tet10Coords, tet10DepartureCoords, tet10Dist, distance_at_point, lengthScale, facets, currentDepth+1, minDepth, maxDepth);
}

void adaptively_append_facets_for_tet_using_semilagrangian_distance(const std::array<stk::math::Vector3d,4> & coords,
  const std::array<stk::math::Vector3d,4> & departureCoords,
  const std::array<double,4> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  if (!have_possibly_cut_edge(departureCoords, distance))
    return;

  const std::array<stk::math::Vector3d,10> tet10Coords = get_tet10_coords_on_tet4(coords);
  const std::array<stk::math::Vector3d,10> tet10DepartureCoords = get_tet10_coords_on_tet4(departureCoords);
  const std::array<double,10> tet10Dist = get_tet10_distance_on_tet4(tet10DepartureCoords, distance, distance_at_point);

  adaptively_append_facets_for_tet10_using_semilagrangian_distance(tet10Coords, tet10DepartureCoords, tet10Dist, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth);
}

}

