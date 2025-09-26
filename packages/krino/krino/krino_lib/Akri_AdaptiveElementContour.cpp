/*
 * Akri_AdaptiveElementContour.cpp
 *
 *  Created on: Mar 20, 2024
 *      Author: drnoble
 */

#include <Akri_AdaptiveElementContour.hpp>
#include <Akri_DiagWriter.hpp>

#include <stk_math/StkVector.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_Sign.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

constexpr double snapTol = 1.e-6;
constexpr double nonlinearTol = 1.e-2;

template <size_t NVERT, size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NNODES> & coords,
    const std::array<double,NNODES> & dist)
{
  for (size_t n1=0; n1<NVERT; ++n1)
  {
    for (size_t n2=n1; n2<NVERT; ++n2)
    {
      const double sqrLen = (coords[n1] - coords[n2]).length_squared();
      if (dist[n1]*dist[n1] <= 0.25*sqrLen || dist[n2]*dist[n2] <= 0.25*sqrLen)
        return true;
    }
  }
  return false;
}

template <size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NNODES> & coords,
    const std::array<double,NNODES> & dist)
{
  for (size_t n1=0; n1<NNODES; ++n1)
  {
    for (size_t n2=n1; n2<NNODES; ++n2)
    {
      const double sqrLen = (coords[n1] - coords[n2]).length_squared();
      if (dist[n1]*dist[n1] <= 0.25*sqrLen || dist[n2]*dist[n2] <= 0.25*sqrLen)
        return true;
    }
  }
  return false;
}

template <size_t NPARENTNODES, size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NPARENTNODES> & coords,
    const std::array<double,NPARENTNODES> & dist,
    const std::array<int, NNODES> & subElemVertexIndices)
{
  for (size_t n1=0; n1<NNODES; ++n1)
  {
    const int pn1 = subElemVertexIndices[n1];
    for (size_t n2=n1; n2<NNODES; ++n2)
    {
      const int pn2 = subElemVertexIndices[n2];
      const double sqrLen = (coords[pn1] - coords[pn2]).length_squared();
      if (dist[pn1]*dist[pn1] <= 0.25*sqrLen || dist[pn2]*dist[pn2] <= 0.25*sqrLen)
        return true;
    }
  }
  return false;
}

double clip_midedge_distance(const double d0, const double d1, const double d2)
{
  const double d25 = 0.75*d0+0.25*d1;
  const double d75 = 0.25*d0+0.75*d1;
  if (d0 < d1)
    return std::max(std::min(d2, d75), d25);
  return std::max(std::min(d2, d25), d75);
}

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

template <size_t NNODES>
void snap_distance(std::array<double,NNODES> & dist, const double snapDistTol)
{
  for (int n=0; n<9; ++n)
    dist[n] = (std::abs(dist[n]) < snapDistTol) ? 0.0 : dist[n];
}

std::array<double,9> interpolate_and_snap_subtri(const std::array<double,6> & tri6Dist, const double snapDistTol)
{
  std::array<double,9> subTriDist = interpolate_subtri(tri6Dist);
  for (int n=0; n<9; ++n)
    subTriDist[n] = (std::abs(subTriDist[n]) < snapDistTol) ? 0.0 : subTriDist[n];
  return subTriDist;
}

int compute_node_caseId(const double dist)
{
  return (dist == 0) ? 1 : ((dist < 0) ? 0 : 2);
}

template <size_t NCOORDS, size_t NDIST>
stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NDIST> & distance,
  const unsigned i0, const unsigned i1, const unsigned i2)
{
  const double loc = find_quadratic_crossing(distance[i0], distance[i1], distance[i2]);
  return (1.-loc) * coords[i0] + loc * coords[i1];
}

void append_facets_for_converged_tri(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<double,6> & tri6Dist,
  const double /*lengthScale*/,
  FacetedSurfaceBase & facets)
{
  const int caseId = compute_node_caseId(tri6Dist[0]) + 3*compute_node_caseId(tri6Dist[1]) + 9*compute_node_caseId(tri6Dist[2]);

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0
      caseId == 26)  // ls[0]>0 && ls[1]>0 && ls[2]>0
    return;

  static const unsigned case_permutations[] =
    { 0, 0, 0, 2, 0, 0, 2, 1, 1, 1, //  0-9
      1, 2, 2, 0, 2, 2, 1, 1, 1, 1, //  10-19
      2, 0, 0, 2, 0, 0, 0 };        //  20-26
  static const unsigned permute_case_ids[] =
    { 0, 1, 2, 1, 4, 5, 2,21,24, 1, //  0-9
      4,21, 4,13,22, 5,22,25, 2, 5, //  10-19
     24,21,22,25,24,25,26 };        //  20-26

  stk::topology topo = stk::topology::TRIANGLE_6_2D;
  std::vector<unsigned> permute(6);
  topo.permutation_node_ordinals(case_permutations[caseId], permute.begin());

  const int permutedCaseId = permute_case_ids[caseId];

  const unsigned i0 = permute[0];
  const unsigned i1 = permute[1];
  const unsigned i2 = permute[2];
  const unsigned i3 = permute[3];
  const unsigned i5 = permute[5];

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
      const stk::math::Vector3d x3 = compute_quadratic_edge_crossing(coords, tri6Dist, i0, i1, i3);
      const stk::math::Vector3d x5 = compute_quadratic_edge_crossing(coords, tri6Dist, i2, i0, i5);
      if (permutedCaseId == 2)
        facets.emplace_back_2d(x5, x3);
      else
        facets.emplace_back_2d(x3, x5);
    }
    break;

    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0
    {
      facets.emplace_back_2d(coords[i0], coords[i1]);
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0
    case 21: // ls[0]<0 && ls[1]=0 && ls[2]>0
    {
      const stk::math::Vector3d x5 = compute_quadratic_edge_crossing(coords, tri6Dist, i2, i0, i5);
      if (permutedCaseId == 5)
        facets.emplace_back_2d(x5, coords[i1]);
      else
        facets.emplace_back_2d(coords[i1], x5);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId,permutedCaseId=" << caseId << "," << permutedCaseId);
  }
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

template <size_t NNODES, size_t NSUBNODES, typename T>
std::array<T, NSUBNODES> subarray(const std::array<T, NNODES> & a, const std::array<int,NSUBNODES> & indices)
{
  std::array<T, NSUBNODES> subarray;
  for (size_t n=0; n<NSUBNODES; ++n)
    subarray[n] = a[indices[n]];
  return subarray;
}

void append_facets_for_subtri_of_converged_tri(const std::array<stk::math::Vector3d,6> & tri6Coords,
  const std::array<double,6> & tri6Dist,
  const std::array<double,9> & subTriDist,
  const double lengthScale,
  FacetedSurfaceBase & facets,
  const std::array<int, 3> & subElemVertexIndices,
  const std::array<int, 3> & subElemMidsideIndices)
{
  const std::array<double,6> subTri6Dist = {{ tri6Dist[subElemVertexIndices[0]], tri6Dist[subElemVertexIndices[1]], tri6Dist[subElemVertexIndices[2]],
      subTriDist[subElemMidsideIndices[0]], subTriDist[subElemMidsideIndices[1]], subTriDist[subElemMidsideIndices[2]] }};
  append_facets_for_converged_tri(subarray(tri6Coords, subElemVertexIndices), subTri6Dist, lengthScale, facets);
}

bool is_edge_converged(const double d0, const double d1, const double d2, const double nonlinearDistTol)
{
  return (std::abs(d2 - 0.5*(d0+d1)) < nonlinearDistTol);
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
  const std::array<double,6> subTri6Dist = {{ tri6Dist[subElemVertexIndices[0]], tri6Dist[subElemVertexIndices[1]], tri6Dist[subElemVertexIndices[2]],
            subTriDist[subElemMidsideIndices[0]], subTriDist[subElemMidsideIndices[1]], subTriDist[subElemMidsideIndices[2]] }};

  std::array<double,9> subSubTriDist = interpolate_subtri(subTri6Dist);
  snap_distance(subSubTriDist, snapTol*lengthScale);

  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{0,3,5}}, {{0,7,5}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{1,4,3}}, {{2,8,1}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{2,5,4}}, {{4,6,3}});
  append_facets_for_subtri_of_converged_tri(subTri6Coords, subTri6Dist, subSubTriDist, lengthScale, facets, {{3,4,5}}, {{8,6,7}});
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
      const std::array<double,6> filteredDistance = apply_tri_snapping_and_clipping(distance, snapTol*lengthScale);
      std::array<double,9> subTriDist = interpolate_subtri(filteredDistance);

      append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{0,3,5}}, {{0,7,5}});
      append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{1,4,3}}, {{2,8,1}});
      append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{2,5,4}}, {{4,6,3}});
      append_facets_for_refined_subtri_using_interpolated_distance(coords, departureCoords, filteredDistance, subTriDist, lengthScale, facets, {{3,4,5}}, {{8,6,7}});
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

