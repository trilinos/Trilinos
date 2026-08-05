// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourTet.hpp>

#include <stk_math/StkVector.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_AdaptiveContourUtils.hpp>
#include <Akri_ContourUtils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

unsigned ContourTet::compute_case_id(const std::array<int,4> & nodeSigns)
{
  return (nodeSigns[0]+1) +
      (nodeSigns[1]+1)*3 +
      (nodeSigns[2]+1)*9 +
      (nodeSigns[3]+1)*27;
}

unsigned ContourTet::get_permutation_for_case(const unsigned caseId)
{
  static constexpr std::array<unsigned,81> casePermutations =
      { 0, 0, 0, 1, 0, 0, 1, 5, 0, 2, //  0-9
        2, 6, 1, 0, 0, 1, 1,10, 2, 2, //  10-19
        2,11, 2, 4, 1, 8, 4, 4, 3, 3, //  20-29
        4, 3, 3, 9, 5, 7, 7, 6, 6, 9, //  30-39
        0, 9, 9, 6, 7, 7, 7, 9,11, 3, //  40-49
        4, 3, 3, 4, 4, 8, 3, 4, 4,11, //  50-59
        4, 2, 2,10, 8, 1,10, 0, 1, 6, //  60-69
        2, 2, 7, 5, 1, 0, 0, 1, 0, 0, //  70-79
        0 };                          //  80
  return casePermutations[caseId];
}

std::array<unsigned, 4> ContourTet::get_permuted_tet4_node_ordinals_for_permutation(const unsigned permutation)
{
  stk::topology topo = stk::topology::TETRAHEDRON_4;
  std::array<unsigned,4> permutedNodeOrdinals;
  topo.permutation_node_ordinals(permutation, permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

std::array<unsigned, 4> ContourTet::get_permuted_tet4_node_ordinals(const unsigned caseId)
{
  return get_permuted_tet4_node_ordinals_for_permutation(get_permutation_for_case(caseId));
}

std::array<unsigned, 10> ContourTet::get_permuted_tet10_node_ordinals_for_permutation(const unsigned permutation)
{
  stk::topology topo = stk::topology::TETRAHEDRON_10;
  std::array<unsigned,10> permutedNodeOrdinals;
  topo.permutation_node_ordinals(permutation, permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

std::array<unsigned, 10> ContourTet::get_permuted_tet10_node_ordinals(const unsigned caseId)
{
  return get_permuted_tet10_node_ordinals_for_permutation(get_permutation_for_case(caseId));
}

unsigned ContourTet::get_permuted_case_id(const unsigned caseId)
{
  static constexpr std::array<unsigned,81> permutedCaseIds =
      { 0, 1, 2, 1, 4, 5, 2, 5, 8, 1, //  0-9
        4, 5, 4,13,14, 5,14,75, 2, 5, //  10-19
        8, 5,14,75, 8,75,78, 1, 4, 5, //  20-29
        4,13,14, 5,14,75, 4,13,14,13, //  30-39
       40,67,14,67,76, 5,14,75,14,67, //  40-49
       76,75,76,79, 2, 5, 8, 5,14,75, //  50-59
        8,75,78, 5,14,75,14,67,76,75, //  60-69
       76,79, 8,75,78,75,76,79,78,79, //  70-79
       80 };                          //  80

  return permutedCaseIds[caseId];
}

const std::array<unsigned, 4> & ContourTet::get_permuted_side_ordinals_for_permutation(const unsigned permutation)
{
  static constexpr std::array<std::array<unsigned,4>,12> permutationSideOrdinals =
      {{ {{ 0, 1, 2, 3 }},
         {{ 1, 2, 0, 3 }},
         {{ 2, 0, 1, 3 }},
         {{ 2, 1, 3, 0 }},
         {{ 1, 3, 2, 0 }},
         {{ 3, 2, 1, 0 }},
         {{ 3, 1, 0, 2 }},
         {{ 1, 0, 3, 2 }},
         {{ 0, 3, 1, 2 }},
         {{ 0, 2, 3, 1 }},
         {{ 2, 3, 0, 1 }},
         {{ 3, 0, 2, 1 }}
      }};
  return permutationSideOrdinals[permutation];
}

const std::array<unsigned, 4> & ContourTet::get_permuted_side_ordinals(const unsigned caseId)
{
  return get_permuted_side_ordinals_for_permutation(get_permutation_for_case(caseId));
}

void ContourTet::append_facets_for_tet4(const std::array<stk::math::Vector3d,4> & coords,
  const std::array<double,4> & dist,
  FacetedSurfaceBase & facets)
{
  const int caseId = compute_case_id({{compute_node_sign(dist[0]), compute_node_sign(dist[1]), compute_node_sign(dist[2]), compute_node_sign(dist[3])}});

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0 && ls[3]<0
      caseId == 80)  // ls[0]>0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    return;

  const std::array<unsigned,4> & i = get_permuted_tet4_node_ordinals(caseId);
  const int permutedCaseId = get_permuted_case_id(caseId);

  switch (permutedCaseId)
  {
    case 0:  // ls[0]<0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 13: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    case 40: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]=0
    case 67: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]>0
    case 76: // ls[0]=0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    case 79: // ls[0]=0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    case 80: // ls[0]>0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      // empty
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 78: // ls[0]<0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x4 = compute_linear_edge_crossing(coords, dist, i[0], i[1]);
      const stk::math::Vector3d x6 = compute_linear_edge_crossing(coords, dist, i[0], i[2]);
      const stk::math::Vector3d x7 = compute_linear_edge_crossing(coords, dist, i[0], i[3]);

      if (2 == permutedCaseId)
        facets.emplace_back_3d(x6, x4, x7);
      else
        facets.emplace_back_3d(x4, x6, x7);
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 75: // ls[0]<0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x6 = compute_linear_edge_crossing(coords, dist, i[0], i[2]);
      const stk::math::Vector3d x7 = compute_linear_edge_crossing(coords, dist, i[0], i[3]);

      if (5 == permutedCaseId)
        facets.emplace_back_3d(x6, coords[i[1]], x7);
      else
        facets.emplace_back_3d(x7, coords[i[1]], x6);
    }
    break;

    case 8:  // ls[0]>0 && ls[1]>0 && ls[2]<0 && ls[3]<0
    {
      const stk::math::Vector3d x5 = compute_linear_edge_crossing(coords, dist, i[1], i[2]);
      const stk::math::Vector3d x6 = compute_linear_edge_crossing(coords, dist, i[0], i[2]);
      const stk::math::Vector3d x7 = compute_linear_edge_crossing(coords, dist, i[0], i[3]);
      const stk::math::Vector3d x8 = compute_linear_edge_crossing(coords, dist, i[1], i[3]);

      // face 4: true: connect 6 and 8, false: connect 7 and 5
      const bool face4 = will_cutting_quad_from_0to2_cut_largest_angle(x8, x7, x6, x5);

      if (face4)
      {
        facets.emplace_back_3d(x8, x7, x6);
        facets.emplace_back_3d(x8, x6, x5);
      }
      else
      {
        facets.emplace_back_3d(x5, x8, x7);
        facets.emplace_back_3d(x5, x7, x6);
      }
    }
    break;

    case 14: // ls[0]>0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    {
      const stk::math::Vector3d x7 = compute_linear_edge_crossing(coords, dist, i[0], i[3]);

      facets.emplace_back_3d(coords[i[1]], x7, coords[i[2]]);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId,permutedCaseId=" << caseId << "," << permutedCaseId);
  }
}

void ContourTet::append_facets_for_tet4_with_tet10_distance(const std::array<stk::math::Vector3d,4> & coords,
  const std::array<double,10> & tet10Dist,
  FacetedSurfaceBase & facets)
{
  const int caseId = compute_case_id({{compute_node_sign(tet10Dist[0]), compute_node_sign(tet10Dist[1]), compute_node_sign(tet10Dist[2]), compute_node_sign(tet10Dist[3])}});

  if (caseId == 0 || // ls[0]<0 && ls[1]<0 && ls[2]<0 && ls[3]<0
      caseId == 80)  // ls[0]>0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    return;

  const std::array<unsigned,10> & i = get_permuted_tet10_node_ordinals(caseId);
  const int permutedCaseId = get_permuted_case_id(caseId);

  switch (permutedCaseId)
  {
    case 0:  // ls[0]<0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 13: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    case 40: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]=0
    case 67: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]>0
    case 76: // ls[0]=0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    case 79: // ls[0]=0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    case 80: // ls[0]>0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      // empty
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 78: // ls[0]<0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x4 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[1], i[4]);
      const stk::math::Vector3d x6 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[2], i[6]);
      const stk::math::Vector3d x7 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[3], i[7]);

      if (2 == permutedCaseId)
        facets.emplace_back_3d(x6, x4, x7);
      else
        facets.emplace_back_3d(x4, x6, x7);
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 75: // ls[0]<0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x6 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[2], i[6]);
      const stk::math::Vector3d x7 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[3], i[7]);

      if (5 == permutedCaseId)
        facets.emplace_back_3d(x6, coords[i[1]], x7);
      else
        facets.emplace_back_3d(x7, coords[i[1]], x6);
    }
    break;

    case 8:  // ls[0]>0 && ls[1]>0 && ls[2]<0 && ls[3]<0
    {
      const stk::math::Vector3d x5 = compute_quadratic_edge_crossing(coords, tet10Dist, i[1], i[2], i[5]);
      const stk::math::Vector3d x6 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[2], i[6]);
      const stk::math::Vector3d x7 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[3], i[7]);
      const stk::math::Vector3d x8 = compute_quadratic_edge_crossing(coords, tet10Dist, i[1], i[3], i[8]);

      // face 4: true: connect 6 and 8, false: connect 7 and 5
      const bool face4 = will_cutting_quad_from_0to2_cut_largest_angle(x8, x7, x6, x5);

      if (face4)
      {
        facets.emplace_back_3d(x8, x7, x6);
        facets.emplace_back_3d(x8, x6, x5);
      }
      else
      {
        facets.emplace_back_3d(x5, x8, x7);
        facets.emplace_back_3d(x5, x7, x6);
      }
    }
    break;

    case 14: // ls[0]>0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    {
      const stk::math::Vector3d x7 = compute_quadratic_edge_crossing(coords, tet10Dist, i[0], i[3], i[7]);

      facets.emplace_back_3d(coords[i[1]], x7, coords[i[2]]);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId,permutedCaseId=" << caseId << "," << permutedCaseId);
  }
}


}



