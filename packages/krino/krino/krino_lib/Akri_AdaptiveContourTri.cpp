/*
 * Akri_AdaptiveElementContour.cpp
 *
 *  Created on: Mar 20, 2024
 *      Author: drnoble
 */

#include <Akri_AdaptiveContourTri.hpp>

#include <Akri_AdaptiveContourUtils.hpp>
#include <Akri_ContourUtils.hpp>
#include <Akri_DiagWriter.hpp>

#include <stk_math/StkVector.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

static constexpr double snapTol = 1.e-6;
static constexpr double nonlinearTol = 1.e-2;

std::array<double,6> apply_tri_snapping_and_clipping(const std::array<double,6> & unfilteredTri6Dist, const double snapDistTol)
{
  std::array<double,6> tri6Dist;
  for (int n=0; n<3; ++n)
    tri6Dist[n] = (std::abs(unfilteredTri6Dist[n]) < snapDistTol) ? 0.0 : unfilteredTri6Dist[n];
  tri6Dist[3] = clip_midedge_distance(tri6Dist[0], tri6Dist[1], unfilteredTri6Dist[3]);
  tri6Dist[4] = clip_midedge_distance(tri6Dist[1], tri6Dist[2], unfilteredTri6Dist[4]);
  tri6Dist[5] = clip_midedge_distance(tri6Dist[2], tri6Dist[0], unfilteredTri6Dist[5]);
  return tri6Dist;
}

std::array<stk::math::Vector3d,6> get_child_tri6_coords(const std::array<stk::math::Vector3d,6>  & parentCoords,
    const std::array<int, 3> & childElemVertexIndices)
{
  std::array<stk::math::Vector3d,6> childTri6Coords =
  {{
      parentCoords[childElemVertexIndices[0]], parentCoords[childElemVertexIndices[1]], parentCoords[childElemVertexIndices[2]],
      0.5*(parentCoords[childElemVertexIndices[0]]+parentCoords[childElemVertexIndices[1]]),
      0.5*(parentCoords[childElemVertexIndices[1]]+parentCoords[childElemVertexIndices[2]]),
      0.5*(parentCoords[childElemVertexIndices[2]]+parentCoords[childElemVertexIndices[0]])
  }};
  return childTri6Coords;
}

std::array<stk::math::Vector3d,6> get_tri6_coords_on_tri3(const std::array<stk::math::Vector3d,3>  & tri3Coords)
{
  std::array<stk::math::Vector3d,6> tri6Coords =
  {{
      tri3Coords[0], tri3Coords[1], tri3Coords[2],
      0.5*(tri3Coords[0]+tri3Coords[1]),
      0.5*(tri3Coords[1]+tri3Coords[2]),
      0.5*(tri3Coords[2]+tri3Coords[0])
  }};
  return tri6Coords;
}

std::array<double,6> get_tri6_distance_on_tri3(const std::array<stk::math::Vector3d,6>  & tri6Coords,
    const std::array<double,3> & tri3Dist,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point)
{
  std::array<double,6> tri6Dist =
  {{
      tri3Dist[0], tri3Dist[1], tri3Dist[2],
      distance_at_point(tri6Coords[3]),
      distance_at_point(tri6Coords[4]),
      distance_at_point(tri6Coords[5]),
  }};
  return tri6Dist;
}

std::array<double,9> interpolate_subtri(const std::array<double,6> & tri6Dist)
{
  std::array<double,9> subTriDist =
  {{
      0.375*tri6Dist[0] + 0.75*tri6Dist[3] - 0.125*tri6Dist[1],
      0.375*tri6Dist[1] + 0.75*tri6Dist[3] - 0.125*tri6Dist[0],
      0.375*tri6Dist[1] + 0.75*tri6Dist[4] - 0.125*tri6Dist[2],
      0.375*tri6Dist[2] + 0.75*tri6Dist[4] - 0.125*tri6Dist[1],
      0.375*tri6Dist[2] + 0.75*tri6Dist[5] - 0.125*tri6Dist[0],
      0.375*tri6Dist[0] + 0.75*tri6Dist[5] - 0.125*tri6Dist[2],
      0.5*tri6Dist[4] + 0.5*tri6Dist[5] - 0.125*tri6Dist[0] + 0.25*tri6Dist[3] - 0.125*tri6Dist[1],
      0.5*tri6Dist[3] + 0.5*tri6Dist[5] - 0.125*tri6Dist[1] + 0.25*tri6Dist[4] - 0.125*tri6Dist[2],
      0.5*tri6Dist[3] + 0.5*tri6Dist[4] - 0.125*tri6Dist[0] + 0.25*tri6Dist[5] - 0.125*tri6Dist[2]
  }};
  return subTriDist;
}

std::array<double,9> interpolate_and_snap_subtri(const std::array<double,6> & tri6Dist, const double snapDistTol)
{
  std::array<double,9> subTriDist = interpolate_subtri(tri6Dist);
  for (int n=0; n<9; ++n)
    subTriDist[n] = (std::abs(subTriDist[n]) < snapDistTol) ? 0.0 : subTriDist[n];
  return subTriDist;
}

std::array<double,6> get_child_tri6_distance(const std::array<double,6> & tri6Dist,
  const std::array<double,9> & subTriDist,
  const std::array<int, 3> & subElemVertexIndices,
  const std::array<int, 3> & subElemMidsideIndices)
{
  const std::array<double,6> subTri6Dist = {{ tri6Dist[subElemVertexIndices[0]], tri6Dist[subElemVertexIndices[1]], tri6Dist[subElemVertexIndices[2]],
      subTriDist[subElemMidsideIndices[0]], subTriDist[subElemMidsideIndices[1]], subTriDist[subElemMidsideIndices[2]] }};
  return subTri6Dist;
}

void append_facets_for_converged_tri(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<double,6> & tri6Dist,
  const double /*lengthScale*/,
  FacetedSurfaceBase & facets)
{
  const int caseId = ContourTri::compute_case_id({{compute_node_sign(tri6Dist[0]), compute_node_sign(tri6Dist[1]), compute_node_sign(tri6Dist[2])}});

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0
      caseId == 26)  // ls[0]>0 && ls[1]>0 && ls[2]>0
    return;

  const std::array<unsigned,6> & i = ContourTri::get_permuted_node_ordinals(caseId);
  const int permutedCaseId = ContourTri::get_permuted_case_id(caseId);

  switch (permutedCaseId)
  {
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0
    case 22: // ls[0]=0 && ls[1]=0 && ls[2]>0
    case 25: // ls[0]=0 && ls[1]>0 && ls[2]>0
      // empty
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0
    case 24: // ls[0]<0 && ls[1]>0 && ls[2]>0
    {
      const stk::math::Vector3d x3 = compute_quadratic_edge_crossing(coords, tri6Dist, i[0], i[1], i[3]);
      const stk::math::Vector3d x5 = compute_quadratic_edge_crossing(coords, tri6Dist, i[2], i[0], i[5]);
      if (permutedCaseId == 2)
        facets.emplace_back_2d(x5, x3);
      else
        facets.emplace_back_2d(x3, x5);
    }
    break;

    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0
    {
      facets.emplace_back_2d(coords[i[0]], coords[i[1]]);
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0
    case 21: // ls[0]<0 && ls[1]=0 && ls[2]>0
    {
      const stk::math::Vector3d x5 = compute_quadratic_edge_crossing(coords, tri6Dist, i[2], i[0], i[5]);
      if (permutedCaseId == 5)
        facets.emplace_back_2d(x5, coords[i[1]]);
      else
        facets.emplace_back_2d(coords[i[1]], x5);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId,permutedCaseId=" << caseId << "," << permutedCaseId);
  }
}

void append_facets_for_subtri_of_converged_tri(const std::array<stk::math::Vector3d,6> & tri6Coords,
  const std::array<double,6> & tri6Dist,
  const std::array<double,9> & subTriDist,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const std::array<int, 3> & subElemVertexIndices,
  const std::array<int, 3> & subElemMidsideIndices)
{
  const std::array<double,6> subTri6Dist = get_child_tri6_distance(tri6Dist, subTriDist, subElemVertexIndices, subElemMidsideIndices);
  append_facets_for_converged_tri(subarray(tri6Coords, subElemVertexIndices), subTri6Dist, lengthScale, facets);
}

int determine_tri_edge_refinement_case_id(const std::array<double,6> & tri6Dist,
  const double lengthScale,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  if (currentDepth < minDepth)
    return 7;
  if (currentDepth == maxDepth)
    return 0;
  const double nonlinearDistTol = nonlinearTol*lengthScale;
  int caseId = 0;
  if (!is_edge_converged(tri6Dist[0], tri6Dist[1], tri6Dist[3], nonlinearDistTol)) caseId += 1;
  if (!is_edge_converged(tri6Dist[1], tri6Dist[2], tri6Dist[4], nonlinearDistTol)) caseId += 2;
  if (!is_edge_converged(tri6Dist[2], tri6Dist[0], tri6Dist[5], nonlinearDistTol)) caseId += 4;
  return caseId;
}

void append_facets_for_refined_subtri_using_interpolated_distance(const std::array<stk::math::Vector3d,6> & tri6Coords,
  const std::array<stk::math::Vector3d,6> & tri6DepartureCoords,
  const std::array<double,6> & tri6Dist,
  const std::array<double,9> & subTriDist,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const std::array<int, 3> & subElemVertexIndices,
  const std::array<int, 3> & subElemMidsideIndices)
{
  if (!have_possibly_cut_edge(tri6DepartureCoords, tri6Dist, subElemVertexIndices))
    return;

  const std::array<stk::math::Vector3d,6> subTri6Coords = get_child_tri6_coords(tri6Coords, subElemVertexIndices);
  const std::array<double,6> subTri6Dist = get_child_tri6_distance(tri6Dist, subTriDist, subElemVertexIndices, subElemMidsideIndices);

  std::array<double,9> subSubTriDist = interpolate_subtri(subTri6Dist);
  snap_distance(subSubTriDist, snapTol*lengthScale);

  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{0,3,5}}, {{0,7,5}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{1,4,3}}, {{2,8,1}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{2,5,4}}, {{4,6,3}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{3,4,5}}, {{8,6,7}});
}

void append_facets_for_converged_tri6(const std::array<stk::math::Vector3d,6> & coords,
  const std::array<stk::math::Vector3d,6> & departureCoords,
  const std::array<double,6> & distance,
  const double lengthScale,
  FacetedSurfaceBase & facets)
{
  const std::array<double,6> filteredDistance = apply_tri_snapping_and_clipping(distance, snapTol*lengthScale);
  std::array<double,9> subTriDist = interpolate_subtri(filteredDistance);

  append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{0,3,5}}, {{0,7,5}});
  append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{1,4,3}}, {{2,8,1}});
  append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{2,5,4}}, {{4,6,3}});
  append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{3,4,5}}, {{8,6,7}});
}

void adaptively_append_facets_for_subtri_using_semilagrangian_distance(const std::array<stk::math::Vector3d,6> & parentCoords,
  const std::array<stk::math::Vector3d,6> & parentDepartureCoords,
  const std::array<double,6> & parentDistance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth,
  const std::array<int, 3> & childElemNodeIndices);

void adaptively_append_facets_for_tri6_using_semilagrangian_distance(const std::array<stk::math::Vector3d,6> & coords,
  const std::array<stk::math::Vector3d,6> & departureCoords,
  const std::array<double,6> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  const int refinementCaseId = determine_tri_edge_refinement_case_id(distance, lengthScale, currentDepth, minDepth, maxDepth);

  switch (refinementCaseId)
  {
    case 0:
    {
      append_facets_for_converged_tri6(coords, departureCoords, distance, lengthScale, facets);
    }
    break;

    case 1:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,3,2}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,2,3}});
    }
    break;

    case 2:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,1,4}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,4,2}});
    }
    break;

    case 3:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,4,3}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,3,2}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{2,3,4}});
    }
    break;

    case 4:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,2,5}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,5,0}});
    }
    break;

    case 5:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,3,5}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{2,5,1}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,5,3}});
    }
    break;

    case 6:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{2,5,4}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,4,0}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,4,5}});
    }
    break;

    case 7:
    {
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{0,3,5}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{1,4,3}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{2,5,4}});
      adaptively_append_facets_for_subtri_using_semilagrangian_distance(coords, departureCoords, distance, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth, {{3,4,5}});
    }
    break;

    default: ThrowRuntimeError("Missing refinement case id =" << refinementCaseId);
  }
}

void adaptively_append_facets_for_subtri_using_semilagrangian_distance(const std::array<stk::math::Vector3d,6> & parentCoords,
  const std::array<stk::math::Vector3d,6> & parentDepartureCoords,
  const std::array<double,6> & parentDistance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth,
  const std::array<int, 3> & childElemNodeIndices)
{
  if (!have_possibly_cut_edge(parentDepartureCoords, parentDistance, childElemNodeIndices))
    return;

  const std::array<stk::math::Vector3d,6> tri6Coords = get_child_tri6_coords(parentCoords, childElemNodeIndices);
  const std::array<stk::math::Vector3d,6> tri6DepartureCoords = get_child_tri6_coords(parentDepartureCoords, childElemNodeIndices);
  const std::array<double,6> tri6Dist = get_tri6_distance_on_tri3(tri6DepartureCoords, subarray(parentDistance, childElemNodeIndices), distance_at_point);

  adaptively_append_facets_for_tri6_using_semilagrangian_distance(tri6Coords, tri6DepartureCoords, tri6Dist, distance_at_point, lengthScale, facets, currentDepth+1, minDepth, maxDepth);
}

void adaptively_append_facets_for_tri_using_semilagrangian_distance(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<stk::math::Vector3d,3> & departureCoords,
  const std::array<double,3> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth)
{
  if (!have_possibly_cut_edge(departureCoords, distance))
    return;

  const std::array<stk::math::Vector3d,6> tri6Coords = get_tri6_coords_on_tri3(coords);
  const std::array<stk::math::Vector3d,6> tri6DepartureCoords = get_tri6_coords_on_tri3(departureCoords);
  const std::array<double,6> tri6Dist = get_tri6_distance_on_tri3(tri6DepartureCoords, distance, distance_at_point);

  adaptively_append_facets_for_tri6_using_semilagrangian_distance(tri6Coords, tri6DepartureCoords, tri6Dist, distance_at_point, lengthScale, facets, currentDepth, minDepth, maxDepth);
}

}

