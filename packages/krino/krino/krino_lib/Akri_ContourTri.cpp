// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourTri.hpp>

#include <stk_math/StkVector.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_AdaptiveContourUtils.hpp>
#include <Akri_ContourUtils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

unsigned ContourTri::compute_case_id(const std::array<int,3> & nodeSigns)
{
  return (nodeSigns[0]+1) +
      (nodeSigns[1]+1)*3 +
      (nodeSigns[2]+1)*9;
}

unsigned ContourTri::get_permutation_for_case(const unsigned caseId)
{
  static constexpr std::array<unsigned,27> casePermutations =
      { 0, 0, 0, 2, 0, 0, 2, 1, 1, 1, //  0-9
        1, 2, 2, 0, 2, 2, 1, 1, 1, 1, //  10-19
        2, 0, 0, 2, 0, 0, 0 };        //  20-26
  return casePermutations[caseId];
}

std::array<unsigned, 3> ContourTri::get_permuted_tri3_node_ordinals(const unsigned caseId)
{
  stk::topology topo = stk::topology::TRIANGLE_3_2D;
  std::array<unsigned,3> permutedNodeOrdinals;
  topo.permutation_node_ordinals(get_permutation_for_case(caseId), permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

std::array<unsigned, 6> ContourTri::get_permuted_tri6_node_ordinals(const unsigned caseId)
{
  stk::topology topo = stk::topology::TRIANGLE_6_2D;
  std::array<unsigned,6> permutedNodeOrdinals;
  topo.permutation_node_ordinals(get_permutation_for_case(caseId), permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

unsigned ContourTri::get_permuted_case_id(const unsigned caseId)
{
  static constexpr std::array<unsigned,27> permutedCaseIds =
      { 0, 1, 2, 1, 4, 5, 2,21,24, 1, //  0-9
        4,21, 4,13,22, 5,22,25, 2, 5, //  10-19
       24,21,22,25,24,25,26 };        //  20-26

  return permutedCaseIds[caseId];
}

void ContourTri::append_facets_for_tri3(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<double,3> & dist,
  FacetedSurfaceBase & facets)
{
  const int caseId = compute_case_id({{compute_node_sign(dist[0]), compute_node_sign(dist[1]), compute_node_sign(dist[2])}});

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0
      caseId == 26)  // ls[0]>0 && ls[1]>0 && ls[2]>0
    return;

  const std::array<unsigned,3> & i = get_permuted_tri3_node_ordinals(caseId);
  const int permutedCaseId = get_permuted_case_id(caseId);

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
      const stk::math::Vector3d x3 = compute_linear_edge_crossing(coords, dist, i[0], i[1]);
      const stk::math::Vector3d x5 = compute_linear_edge_crossing(coords, dist, i[2], i[0]);
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
      const stk::math::Vector3d x5 = compute_linear_edge_crossing(coords, dist, i[2], i[0]);
      if (permutedCaseId == 5)
        facets.emplace_back_2d(x5, coords[i[1]]);
      else
        facets.emplace_back_2d(coords[i[1]], x5);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId,permutedCaseId=" << caseId << "," << permutedCaseId);
  }
}

void ContourTri::append_facets_for_tri3_with_tri6_distance(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<double,6> & tri6Dist,
  FacetedSurfaceBase & facets)
{
  const int caseId = compute_case_id({{compute_node_sign(tri6Dist[0]), compute_node_sign(tri6Dist[1]), compute_node_sign(tri6Dist[2])}});

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0
      caseId == 26)  // ls[0]>0 && ls[1]>0 && ls[2]>0
    return;

  const std::array<unsigned,6> & i = get_permuted_tri6_node_ordinals(caseId);
  const int permutedCaseId = get_permuted_case_id(caseId);

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

}
