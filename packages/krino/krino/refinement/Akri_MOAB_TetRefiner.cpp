/*
 * Akri_MOAB_TetRefiner.cpp
 *
 *  Created on: Oct 18, 2022
 *      Author: drnoble
 */




// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QualityMetric.hpp>
#include "Akri_MOAB_TetRefiner.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stack>
#include <string>
#include <stdexcept>

#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

inline double SQR(const double x) { return x*x; }

inline double compute_edge_length_squared(const double *c0, const double *c1)
{
  return SQR(c0[0]-c1[0]) + SQR(c0[1]-c1[1]) + SQR(c0[2]-c1[2]);
}

double MeanRatioMetricForTetRefinement::element_quality(const int *indices, const double *coords[10]) const
{
  const std::array<stk::math::Vector3d,4> nodeLocations{stk::math::Vector3d(coords[indices[0]]), stk::math::Vector3d(coords[indices[1]]), stk::math::Vector3d(coords[indices[2]]), stk::math::Vector3d(coords[indices[3]])};
  return MeanRatioQualityMetric::tet_mean_ratio(nodeLocations);
}

double LengthRatioMetricForTetRefinement::element_quality(const int *indices, const double *coords[10]) const
{
  double edge_min=std::numeric_limits<double>::max();
  double edge_max = 0;
  for (int i=0; i < 3; i++)
    {
      for (int j=i+1; j < 4; j++)
        {
          const double el2 = compute_edge_length_squared(coords[indices[i]], coords[indices[j]]);
          edge_min = std::min(edge_min, el2);
          edge_max = std::max(edge_max, el2);
        }
    }
  return std::sqrt(edge_max/edge_min);
}

namespace moab {

template <size_t N>
int * best_tets_if_lower_quality_is_better( const ElementMetricForTetRefinement & qualityMetric, int * templates, const std::array<int,N> & alternates, const double* coords[10] )
{
  // force better tets to be at least this factor better - helps with tie breaks for structured meshes
  const double factor = 0.95;

  double best_qual=0;
  int * bestTemplate = nullptr;
  for (size_t i=0; i<N; ++i)
    {
      int * current_template = templates + alternates[i];
      // find worst quality element
      double max_qual=0;
      for (int j=0; j < current_template[0]; j++)
        {
          max_qual = std::max(max_qual, qualityMetric.element_quality(current_template + 1 + j*4, coords));
        }
      // find alternates with the best (min) worst quality
      if (i == 0 || max_qual < best_qual*factor)
        {
          best_qual = max_qual;
          bestTemplate = current_template;
        }
    }
  return bestTemplate;
}

template <size_t N>
int * best_tets_if_higher_quality_is_better( const ElementMetricForTetRefinement & qualityMetric, int * templates, const std::array<int,N> & alternates, const double* coords[10] )
{
  // force better tets to be at least this factor better - helps with tie breaks for structured meshes
  const double factor = 1.01;

  double best_qual=0;
  int * bestTemplate = nullptr;
  for (size_t i=0; i<N; ++i)
    {
      int * current_template = templates + alternates[i];
      // find worst quality element
      double min_qual=1.;
      for (int j=0; j < current_template[0]; j++)
        {
          min_qual = std::min(min_qual, qualityMetric.element_quality(current_template + 1 + j*4, coords));
        }
      // find alternates with the best (min) worst quality
      if (i == 0 || min_qual > best_qual*factor)
        {
          best_qual = min_qual;
          bestTemplate = current_template;
        }
    }
  return bestTemplate;
}

template <size_t N>
int * best_tets( const ElementMetricForTetRefinement & qualityMetric, int * templates, const std::array<int,N> & alternates, const double* coords[10] )
{
  if (qualityMetric.is_higher_quality_better())
    return best_tets_if_higher_quality_is_better(qualityMetric, templates, alternates, coords);
  return best_tets_if_lower_quality_is_better(qualityMetric, templates, alternates, coords);
}

unsigned SimplexTemplateRefiner::determine_permuted_case_id_tet4(const unsigned caseId)
{
  static constexpr std::array<unsigned,64> permutedCaseIds{0,1,1,3,1,3,3,7,1,3,10,11,3,13,14,15,1,3,3,13,10,14,11,15,3,7,14,15,11,15,30,31,1,10,3,14,3,11,13,15,3,14,11,30,7,15,15,31,3,11,7,15,14,30,15,31,13,15,15,31,15,31,31,63};
  return permutedCaseIds[caseId];
}

unsigned SimplexTemplateRefiner::determine_permutation_tet4(const unsigned caseId)
{
  static constexpr std::array<unsigned,64> permutations{0,0,1,0,2,2,1,0,3,5,0,0,8,0,0,0,4,4,11,1,1,1,1,1,3,3,5,5,3,3,0,0,7,2,10,2,6,2,2,2,7,4,7,2,6,6,7,2,9,4,9,9,8,1,11,1,4,4,10,4,8,3,7,0};
  return permutations[caseId];
}

unsigned SimplexTemplateRefiner::num_new_child_elements_tet4(const int caseId)
{
  const unsigned permutedCaseId = SimplexTemplateRefiner::determine_permuted_case_id_tet4(caseId);
  switch(permutedCaseId)
  {
    case 0:  // edge marks = 000000, Bits shown left to right
        return 0;
    case 1:  // edge marks = 100000, Bits shown left to right
        return 2;
    case 3:  // edge marks = 110000, Bits shown left to right
        return 3;
    case 7:  // edge marks = 111000, Bits shown left to right
        return 4;
    case 10: // edge marks = 010100, Bits shown left to right
        return 4;
    case 11: // edge marks = 110100, Bits shown left to right
        return 5;
    case 13: // edge marks = 101100, Bits shown left to right
        return 4;
    case 14: // edge marks = 011100, Bits shown left to right
        return 5;
    case 15: // edge marks = 111100, Bits shown left to right
        return 6;
    case 30: // edge marks = 011110, Bits shown left to right
        return 6;
    case 31: // edge marks = 111110, Bits shown left to right
        return 7;
    case 63: // edge marks = 111111, Bits shown left to right
        return 8;
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << caseId << " not supported in num_new_child_elements_tet4.";
      throw std::runtime_error(errorMsg.str());
    }
  }
}

  //p this is from some comments below, and coded in tet_edges.  It's only used for disambiguation of same-length edges.
  //p   * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
  int tet_edges[6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};

  static int determine_tet_side_from_side_nodes(const std::array<int,3> & sideNodes)
  {
    static const std::array<std::vector<int>,10> sortedSidesByNode{{ {0,2,3}, {0,1,3}, {1,2,3}, {0,1,2}, {0,3} ,{1,3}, {2,3}, {0,2}, {0,1}, {1,2} }};

    std::vector<int> nodeSideIntersection;
    for (unsigned i = 0; i < sideNodes.size(); ++i)
    {
      std::vector<int> nodeSides = sortedSidesByNode[sideNodes[i]];
      if (i == 0)
      {
        nodeSideIntersection.swap(nodeSides);
      }
      else
      {
        std::vector<int> workingSet;
        workingSet.swap(nodeSideIntersection);
        std::set_intersection(workingSet.begin(),workingSet.end(),nodeSides.begin(),nodeSides.end(),std::back_inserter(nodeSideIntersection));
      }
    }

    if (!nodeSideIntersection.empty())
    {
      STK_ThrowRequire(1 == nodeSideIntersection.size());
      return nodeSideIntersection[0];
    }

    return -1;
  }

  static std::array<int,4> determine_tet_sides_from_nodes(const std::array<int, 4> & tetNodes)
  {
    std::array<int,4> tetSides;
    stk::topology topo = stk::topology::TETRAHEDRON_4;
    std::array<int,3> tetSideNodes;
    for (int iSide=0; iSide<4; ++iSide)
    {
      topo.face_nodes(tetNodes.data(), iSide, tetSideNodes.data());
      tetSides[iSide] = determine_tet_side_from_side_nodes(tetSideNodes);
    }
    return tetSides;
  }

  std::vector<SimplexTemplateRefiner::TetDescription> SimplexTemplateRefiner::refinement_child_nodes_and_sides_tet4(const ElementMetricForTetRefinement & qualityMetric,
      const unsigned encodedEdgesToRefine,
      const std::array<stk::math::Vector3d,10> & elementNodeCoords,
      const std::array<int,4> & elementNodeScore,
      const bool needSides)
  {
    const std::vector<std::array<int,4>> newTetNodes = moab::SimplexTemplateRefiner::refinement_child_nodes_tet4(qualityMetric, encodedEdgesToRefine, elementNodeCoords, elementNodeScore);

    const unsigned numTets = newTetNodes.size();

    std::vector<SimplexTemplateRefiner::TetDescription> newTets;
    newTets.resize(numTets);
    for (unsigned iTet=0; iTet<numTets; ++iTet)
    {
      for (int iNode=0; iNode<4; ++iNode)
        newTets[iTet].nodeIds[iNode] = newTetNodes[iTet][iNode];

      if (needSides)
      {
        const std::array<int,4> tetSides = determine_tet_sides_from_nodes(newTets[iTet].nodeIds);
        for (int iSide=0; iSide<4; ++iSide)
          newTets[iTet].sideIds[iSide] = tetSides[iSide];
      }
      else
      {
        for (int iSide=0; iSide<4; ++iSide)
          newTets[iTet].sideIds[iSide] = -1; // -1 indicates side is not on any parent side
      }
    }

    return newTets;
  }

  /**\brief Refine a tetrahedron.

   */
  std::vector<TetTupleInt> SimplexTemplateRefiner::refinement_child_nodes_tet4(const ElementMetricForTetRefinement & qualityMetric, const unsigned edge_code, const std::array<stk::math::Vector3d,10> & elementNodeCoords, const std::array<int,4> & elementNodeScore )
  {  // refine_3_simplex

    std::vector<TetTupleInt> new_tets;

    if ( ! edge_code )
      {
        // No edges to subdivide
        return new_tets;
      }

    // Generate tetrahedra that are compatible except when edge
    // lengths are equal on indeterminately subdivided faces.
    const double* permuted_coords[10];
    EntityHandle permuted_hash[4];
    int permuted_local_ids[10];   // this is just a copy of permutations_from_index, for readability
    double permlen[6]; // permuted edge lengths
    int C = SimplexTemplateRefiner::template_index[edge_code][0];
    int P = SimplexTemplateRefiner::template_index[edge_code][1];

    // 1. Permute the tetrahedron into our canonical configuration
    for ( int i = 0; i < 4; ++ i )
      {
        permuted_hash[i] = elementNodeScore[SimplexTemplateRefiner::permutations_from_index[P][i]];
      }
    for ( int i = 0; i < 10; ++ i )
      {
        permuted_local_ids[i] = SimplexTemplateRefiner::permutations_from_index[P][i];
        permuted_coords[i] = elementNodeCoords[SimplexTemplateRefiner::permutations_from_index[P][i]].data();
      }

    for (int i = 0; i < 6; i++)
      {
        const double *c0 = permuted_coords[tet_edges[i][0]];
        const double *c1 = permuted_coords[tet_edges[i][1]];
        if (permuted_hash[tet_edges[i][0]] > permuted_hash[tet_edges[i][1]])
          {
            c0 = permuted_coords[tet_edges[i][1]];
            c1 = permuted_coords[tet_edges[i][0]];
          }
        permlen[i] = compute_edge_length_squared(c0, c1);
      }

    int comparison_bits;
    int comparison_bits_save;
    std::stack<int*> output_tets;
    std::stack<int*> output_perm;
    std::stack<int>  output_sign;

    // cout << "Case " << C << "  Permutation " << P << endl;
    // 2. Generate tetrahedra based on the configuration.
    //    Note that case 0 is handled above (edgeCode == 0).

    switch ( C )
      {
      case 1: // Ruprecht-Muller Case 1
        output_tets.push( SimplexTemplateRefiner::templates + 0 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        break;
      case 2: // Ruprecht-Muller Case 2a
        comparison_bits =
          ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
          0;

        // a tie-breaker for a pair of edges that have the same length - add all the vertex handles
        //   and choose one case or the other depending on if the sum is odd or even

#define V0(iedge) permuted_hash[tet_edges[iedge][0]]
#define V1(iedge) permuted_hash[tet_edges[iedge][1]]

#define VH(iedge) (V0(iedge) + V1(iedge))
#define CMP_VH(ie,je) ( ( std::min(V0(ie),V1(ie)) == std::min(V0(je),V1(je)) ) ? (VH(ie) < VH(je)) : ( std::min(V0(ie),V1(ie)) < std::min(V0(je),V1(je)) ) )

        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(0,1))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }

        output_tets.push( SimplexTemplateRefiner::templates + 9 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        switch ( comparison_bits )
          {
          case 2: // 0>1
            output_tets.push( SimplexTemplateRefiner::templates + 14 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 1: // 1>0
            output_tets.push( SimplexTemplateRefiner::templates + 14 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );
            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 3: // Ruprecht-Muller Case 2b
        output_tets.push( SimplexTemplateRefiner::templates + 23 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        break;
      case 4: // Ruprecht-Muller Case 3a
        //0<3<2<0 ?

        comparison_bits =
          ( permlen[0] <= permlen[3] ? 1 : 0 ) | ( permlen[0] >= permlen[3] ? 2 : 0 ) |
          ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
          ( permlen[0] <= permlen[2] ? 16 : 0 ) | ( permlen[0] >= permlen[2] ? 32 : 0 ) |
          0;
        comparison_bits_save = comparison_bits;
        if ( ( comparison_bits_save & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(0,3))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits_save & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(2,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }
        if ( ( comparison_bits_save & 48 ) == 48 )
          {
            comparison_bits -=  48 ;
            if (CMP_VH(0,2))
              comparison_bits |= 16;
            else
              comparison_bits |= 32;
          }

        output_tets.push( SimplexTemplateRefiner::templates + 40 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        // here
        switch ( comparison_bits )
          {
          case 42: // 0>2>3<0
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 25: // 2>3>0<2
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[11] );
            output_sign.push( 1 );
            break;
          case 37: // 3>0>2<3
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[3] );
            output_sign.push( 1 );
            break;
          case 21: // 3>2>0<3
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[22] );
            output_sign.push( -1 );
            break;
          case 26: // 2>0>3<2
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[12] );
            output_sign.push( -1 );
            break;
          case 38: // 0>3>2<0
            output_tets.push( SimplexTemplateRefiner::templates + 45 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 5: // Ruprecht-Muller Case 3b
        output_tets.push( SimplexTemplateRefiner::templates + 58 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        break;
      case 6: // Ruprecht-Muller Case 3c
        comparison_bits =
          ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
          ( permlen[0] <= permlen[3] ? 4 : 0 ) | ( permlen[0] >= permlen[3] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(0,1))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(0,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }
        switch ( comparison_bits )
          {
          case 10: // 0>1,0>3
            output_tets.push( SimplexTemplateRefiner::templates + 75 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 5: // 1>0,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 96 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 6: // 0>1,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 117 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 9: // 1>0,0>3
            output_tets.push( SimplexTemplateRefiner::templates + 138 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 7: // Ruprecht-Muller Case 3d
        comparison_bits =
          ( permlen[0] <= permlen[2] ? 1 : 0 ) | ( permlen[0] >= permlen[2] ? 2 : 0 ) |
          ( permlen[0] <= permlen[4] ? 4 : 0 ) | ( permlen[0] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(0,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(0,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }

        switch ( comparison_bits )
          {
          case 10: // 0>4,0>2
            output_tets.push( SimplexTemplateRefiner::templates + 159 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 5: // 4>0,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 180 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 0>4,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 201 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 6: // 4>0,0>2
            output_tets.push( SimplexTemplateRefiner::templates + 222 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 8: // Ruprecht-Muller Case 4a
        comparison_bits =
          ( permlen[4] <= permlen[5] ? 1 : 0 ) | ( permlen[4] >= permlen[5] ? 2 : 0 ) |
          ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(4,5))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(3,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }

        output_tets.push( SimplexTemplateRefiner::templates + 243 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        switch ( comparison_bits )
          {
          case 5: // 5>4>3
            output_tets.push( SimplexTemplateRefiner::templates + 252 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 10: // 3>4>5
            output_tets.push( SimplexTemplateRefiner::templates + 252 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );

            break;
          case 6: // 3<4>5
            output_tets.push( SimplexTemplateRefiner::templates + 269 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 3>4<5
            output_tets.push( SimplexTemplateRefiner::templates + 286 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 9: // Ruprecht-Muller Case 4b
        comparison_bits =
          ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
          ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
          ( permlen[3] <= permlen[4] ? 16 : 0 ) | ( permlen[3] >= permlen[4] ? 32 : 0 ) |
          ( permlen[1] <= permlen[4] ? 64 : 0 ) | ( permlen[1] >= permlen[4] ? 128 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(1,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(2,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }
        if ( ( comparison_bits & 48 ) == 48 )
          {
            comparison_bits -=  48 ;
            if (CMP_VH(3,4))
              comparison_bits |= 16;
            else
              comparison_bits |= 32;
          }
        if ( ( comparison_bits & 192 ) == 192 )
          {
            comparison_bits -=  192 ;
            if (CMP_VH(1,4))
              comparison_bits |= 64;
            else
              comparison_bits |= 128;
          }

        switch ( comparison_bits )
          {
          case 85: // 2>1,3>2,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 102: // 1>2,3>2,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            break;
          case 170: // 1>2,2>3,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          case 153: // 2>1,2>3,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            break;
          case 90: // 1>2,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            break;
          case 105: // 2>1,2>3,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            break;
          case 165: // 2>1,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            break;
          case 150: // 1>2,3>2,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 303 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            break;
          case 101: // 2>1,3>2,3>4,4>1
            {
              static constexpr std::array<int,2> alternates{ 328, 353 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 86: // 1>2,3>2,4>3,4>1
            {
              static constexpr std::array<int,2> alternates{ 328, 353 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            break;
          case 154: // 1>2,2>3,4>3,1>4
            {
              static constexpr std::array<int,2> alternates{ 328, 353 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            break;
          case 169: // 2>1,2>3,3>4,1>4
            {
              static constexpr std::array<int,2> alternates{ 328, 353 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          case 89: // 2>1,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 378 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 166: // 1>2,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 378 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 10: // Ruprecht-Muller Case 5
        comparison_bits =
          ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
          ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
            comparison_bits -=  3 ;
            if (CMP_VH(1,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
            comparison_bits -=  12 ;
            if (CMP_VH(3,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
          }

        output_tets.push( SimplexTemplateRefiner::templates + 403 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        switch ( comparison_bits )
          {
          case 10: // 1>2,3>4
            output_tets.push( SimplexTemplateRefiner::templates + 412 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 5: // 2>1,4>3
            output_tets.push( SimplexTemplateRefiner::templates + 412 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );

            break;
          case 6: // 1>2,4>3
            {
              static constexpr std::array<int,2> alternates{ 433, 454 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 2>1,3>4
            {
              static constexpr std::array<int,2> alternates{ 433, 454 };
              output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );

            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 11: // Ruprecht-Muller Case 6
        output_tets.push( SimplexTemplateRefiner::templates + 475 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        {
          static constexpr std::array<int,3> alternates{ 492, 509, 526 };
          output_tets.push( best_tets( qualityMetric, SimplexTemplateRefiner::templates, alternates, permuted_coords ) );
        }
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        break;
      }

    int* tets;
    int  ntets;
    int* perm;
    int  sgn;
    while ( ! output_tets.empty() )
      {
        tets = output_tets.top();
        ntets = *tets;
        tets++;
        perm = output_perm.top();
        sgn = output_sign.top();

        output_tets.pop();
        output_perm.pop();
        output_sign.pop();

        int t;
        if ( sgn > 0 )
          {
            for ( t = 0; t < ntets; ++t )
              {
                TetTupleInt nt = {permuted_local_ids[perm[tets[0]]],
                                  permuted_local_ids[perm[tets[1]]],
                                  permuted_local_ids[perm[tets[2]]],
                                  permuted_local_ids[perm[tets[3]]]};

                new_tets.push_back(nt);

                tets += 4;
              }
          }
        else
          {
            // we have tets that were specified with a negative permutation... reverse the first 2 vertices
            // so the orientation is positive.
            for ( t = 0; t < ntets; ++t )
              {
                TetTupleInt nt = {permuted_local_ids[perm[tets[1]]],
                                  permuted_local_ids[perm[tets[0]]],
                                  permuted_local_ids[perm[tets[2]]],
                                  permuted_local_ids[perm[tets[3]]]};

                new_tets.push_back(nt);


                tets += 4;
              }
          }
      }

    STK_ThrowRequireMsg(num_new_child_elements_tet4(edge_code) == new_tets.size(), "Mismatch of size " << num_new_child_elements_tet4(edge_code) << "  " << new_tets.size() << " for case " << edge_code); // FIXME: remove check

    return new_tets;
  }

  /*
   * The array below is indexed by the edge code for a tetrahedron.
   * Looking up a row with a tet's edge code will return C and P.
   * C is a configuration number and P is a permutation index.
   *
   * C is based on the case number from Ruprecht and
   * Muller's (1998) paper on adaptive tetrahedra. (The case
   * numbers are shown to the left of the row in the column
   * labeled case. The only difference is that we introduce
   * a case 3d which is part of case 3c in the paper.)
   *
   * P is an index into the permutations_from_index array below,
   * and is used to transform the current tetrahedron into
   * the canonical configuration associated with C.
   *
   * The 6-digit binary number to the left (which is shown in
   * the horribly UNconventional LSB->MSB order) is the edge
   * code for the row. The 6 digits correspond to the 6 edges
   * of the tetrahedron; a '0' implies no subdivision while
   * a '1' implies subdivision should occur. The ordering of
   * the bits is
   *
   * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
   *
   * where the numbers are vertices of the tetrahedron 0-1-2-3.
   * Note that Tet 0-1-2-3 must be positive (i.e., the plane
   * specified by Triangle 0-1-2 must have a normal pointing
   * towards vertex 3, and Triangle 0-1-2's normal must be
   * calculated using the cross-product (Edge 0-1) x (Edge 0-2)).
   *
   * ===========
   * References:
   * (Ruprect and Muller, 1998) A Scheme for Edge-based Adaptive
   *   Tetrahedron Subdivision, Mathematical Visualization (eds.
   *   Hege and Polthier), pp. 61--70. Springer-Verlag. 1998.
   */
  int SimplexTemplateRefiner::template_index[64][2] =
    {
      /*      code case      C    P */
      /* 000000  0  0  */ {  0,   0 },
      /* 100000  1  1  */ {  1,   0 },
      /* 010000  2  1  */ {  1,   1 },
      /* 110000  3  2a */ {  2,   0 },
      /* 001000  4  1  */ {  1,   2 },
      /* 101000  5  2a */ {  2,   2 },
      /* 011000  6  2a */ {  2,   1 },
      /* 111000  7  3b */ {  5,  11 },
      /* 000100  8  1  */ {  1,  10 },
      /* 100100  9  2a */ {  2,   5 },
      /* 010100 10  2b */ {  3,   1 },
      /* 110100 11  3c */ {  6,   0 },
      /* 001100 12  2a */ {  2,  10 },
      /* 101100 13  3a */ {  4,   0 },
      /* 011100 14  3d */ {  7,   2 },
      /* 111100 15  4a */ {  8,   6 },
      /* 000010 16  1  */ {  1,   6 },
      /* 100010 17  2a */ {  2,   4 },
      /* 010010 18  2a */ {  2,   8 },
      /* 110010 19  3a */ {  4,   1 },
      /* 001010 20  2b */ {  3,   2 },
      /* 101010 21  3d */ {  7,   0 },
      /* 011010 22  3c */ {  6,   1 },
      /* 111010 23  4a */ {  8,   9 },
      /* 000110 24  2a */ {  2,   3 },
      /* 100110 25  3b */ {  5,   0 },
      /* 010110 26  3d */ {  7,   4 },
      /* 110110 27  4a */ {  8,  11 },
      /* 001110 28  3c */ {  6,  10 },
      /* 101110 29  4a */ {  8,   7 },
      /* 011110 30  4b */ {  9,   0 },
      /* 111110 31  5  */ { 10,   7 },
      /* 000001 32  1  */ {  1,   7 },
      /* 100001 33  2b */ {  3,   0 },
      /* 010001 34  2a */ {  2,   7 },
      /* 110001 35  3d */ {  7,   1 },
      /* 001001 36  2a */ {  2,  11 },
      /* 101001 37  3c */ {  6,   2 },
      /* 011001 38  3a */ {  4,   2 },
      /* 111001 39  4a */ {  8,   3 },
      /* 000101 40  2a */ {  2,   9 },
      /* 100101 41  3d */ {  7,  10 },
      /* 010101 42  3c */ {  6,   7 },
      /* 110101 43  4b */ {  9,   2 },
      /* 001101 44  3b */ {  5,   7 },
      /* 101101 45  4a */ {  8,   8 },
      /* 011101 46  4a */ {  8,   4 },
      /* 111101 47  5  */ { 10,   6 },
      /* 000011 48  2a */ {  2,   6 },
      /* 100011 49  3c */ {  6,   4 },
      /* 010011 50  3b */ {  5,   1 },
      /* 110011 51  4a */ {  8,  10 },
      /* 001011 52  3d */ {  7,   7 },
      /* 101011 53  4b */ {  9,   1 },
      /* 011011 54  4a */ {  8,   5 },
      /* 111011 55  5  */ { 10,  10 },
      /* 000111 56  3a */ {  4,  10 },
      /* 100111 57  4a */ {  8,   1 },
      /* 010111 58  4a */ {  8,   2 },
      /* 110111 59  5  */ { 10,   2 },
      /* 001111 60  4a */ {  8,   0 },
      /* 101111 61  5  */ { 10,   1 },
      /* 011111 62  5  */ { 10,   0 },
      /* 111111 63  6  */ { 11,   0 },
    };


  /* Does this mean anything? If so, then you are either
   * superstitious or much more clever than I (or both?).
   */
  /* permutation index, P:  0  1  2  3  4  5  6  7  8  9 10 11 */
  /* number of references: 12  9  9  3  4  2  5  6  2  3  7  2 */


  /*
   * The array below is a list of all the _positive_
   * permutations of Tetrahedron 0-1-2-3. Given a
   * permutation index, it returns a row of 10 values:
   * these are the vertex numbers of the permuted
   * tetrahedron. The first 4 values are the permuted
   * corner indices, the next 6 values are the
   * permuted edge midpoint indices.
   *
   * There are 24 entries, 6 for each of the 4 faces of
   * the tetrahedron.
   */
  int SimplexTemplateRefiner::permutations_from_index[24][10] =
    {
      /* corners      midpoints          face points   */
      /* POSITIVE ARRANGEMENTS                         */
      { 0, 1, 2, 3,   4, 5, 6, 7, 8, 9 }, /* Face 0-1-2 */
      { 1, 2, 0, 3,   5, 6, 4, 8, 9, 7 },
      { 2, 0, 1, 3,   6, 4, 5, 9, 7, 8 },

      { 0, 3, 1, 2,   7, 8, 4, 6, 9, 5 }, /* Face 0-3-1 */
      { 3, 1, 0, 2,   8, 4, 7, 9, 5, 6 },
      { 1, 0, 3, 2,   4, 7, 8, 5, 6, 9 },

      { 1, 3, 2, 0,   8, 9, 5, 4, 7, 6 }, /* Face 1-3-2 */
      { 3, 2, 1, 0,   9, 5, 8, 7, 6, 4 },
      { 2, 1, 3, 0,   5, 8, 9, 6, 4, 7 },

      { 2, 3, 0, 1,   9, 7, 6, 5, 8, 4 }, /* Face 2-3-0 */
      { 3, 0, 2, 1,   7, 6, 9, 8, 4, 5 },
      { 0, 2, 3, 1,   6, 9, 7, 4, 5, 8 },

      /* NEGATIVE ARRANGEMENTS                         */
      { 0, 2, 1, 3,   6, 5, 4, 7, 9, 8 }, /* Face 0-1-2 */
      { 2, 1, 0, 3,   5, 4, 6, 9, 8, 7 },
      { 1, 0, 2, 3,   4, 6, 5, 8, 7, 9 },

      { 0, 1, 3, 2,   4, 8, 7, 6, 5, 9 }, /* Face 0-3-1 */
      { 1, 3, 0, 2,   8, 7, 4, 5, 9, 6 },
      { 3, 0, 1, 2,   7, 4, 8, 9, 6, 5 },

      { 1, 2, 3, 0,   5, 9, 8, 4, 6, 7 }, /* Face 1-3-2 */
      { 2, 3, 1, 0,   9, 8, 5, 6, 7, 4 },
      { 3, 1, 2, 0,   8, 5, 9, 7, 4, 6 },

      { 2, 0, 3, 1,   6, 7, 9, 5, 4, 8 }, /* Face 2-3-0 */
      { 0, 3, 2, 1,   7, 9, 6, 4, 8, 5 },
      { 3, 2, 0, 1,   9, 6, 7, 8, 5, 4 }
    };

  /*
   * Below is a list of output tetrahedra. The array is
   * generated by TessellatorGenerator.py
   * which also generates the code that references it.
   * Each set of tetrahedra begins with a single integer
   * that is the number of tetrahedra for that particular
   * case. It is followed by 5 integers for each output
   * tetrahedron; the first four numbers on each row are
   * indices of the output tetrahedron. The final number
   * is a bit vector specifying which edges of the
   * tetrahedron are internal to the parent tetrahedron
   * being decomposed.
   *
   * Multiple lists of output tetrahedra may be
   * combined to create the tessellation of a single
   * input tetrahedron.
   */

  int SimplexTemplateRefiner::templates[] =
    {
      // case 1_0, orig start 0, new start 0
      2,
      0,  4,  2,  3,
      4,  1,  2,  3,

      // case 2a_0, orig start 9, new start 9
      1,
      3,  4,  5,  1,

      // case 2a, 0>1, orig start 14, new start 14
      2,
      0,  4,  2,  3,
      4,  5,  2,  3,

      // case 2b_0, orig start 40, new start 23
      4,
      0,  4,  9,  3,
      4,  1,  9,  3,
      0,  4,  2,  9,
      4,  1,  2,  9,

      // case 3a_0, orig start 57, new start 40
      1,
      4,  7,  6,  0,

      // case 3a, 0>2>3<0, orig start 62, new start 45
      3,
      1,  3,  2,  4,
      4,  6,  3,  2,
      4,  6,  7,  3,

      // case 3b_0, orig start 162, new start 58
      4,
      0,  7,  4,  2,
      4,  7,  8,  2,
      4,  8,  1,  2,
      7,  3,  8,  2,

      // case 3c, 0>1,0>3, orig start 179, new start 75
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3c, 1>0,3>0, orig start 200, new start 96
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 0>1,3>0, orig start 221, new start 117
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 1>0,0>3, orig start 242, new start 138
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3d, 0>4,0>2, orig start 362, new start 159
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,2>0, orig start 383, new start 180
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 0>4,2>0, orig start 404, new start 201
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,0>2, orig start 425, new start 222
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 4a_0, orig start 545, new start 243
      2,
      7,  8,  9,  3,
      7,  9,  8,  6,

      // case 4a, 5>4>3, orig start 554, new start 252
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      9,  1,  6,  2,
      9,  1,  8,  6,

      // case 4a, 3<4>5, orig start 571, new start 269
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      8,  2,  6,  9,
      8,  2,  1,  6,

      // case 4a, 3>4<5, orig start 588, new start 286
      4,
      6,  9,  8,  1,
      6,  9,  1,  2,
      6,  7,  0,  1,
      6,  7,  1,  8,

      // case 4b, 2>1,3>2,4>3,4>1, orig start 688, new start 303
      6,
      6,  8,  1,  5,
      6,  8,  0,  1,
      6,  8,  7,  0,
      6,  8,  2,  7,
      7,  8,  2,  3,
      6,  8,  5,  2,

      // case 4b, 2>1,3>2,3>4,4>1, orig start 713, new start 328
      6,
      6,  8,  1,  5,
      6,  8,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      6,  8,  5,  2,
      6,  8,  2,  7,

      // case 4b, 2>1,3>2,3>4,4>1, a, orig start 738, new start 353
      6,
      7,  8,  1,  5,
      6,  5,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      7,  8,  5,  2,
      6,  5,  2,  7,

      // case 4b, 2>1,2>3,4>3,4>1, orig start 763, new start 378
      6,
      6,  8,  5,  2,
      6,  8,  2,  3,
      6,  8,  3,  7,
      6,  8,  7,  0,
      6,  8,  0,  1,
      6,  8,  1,  5,

      // case 5_0, orig start 1107, new start 403
      2,
      7,  8,  9,  3,
      6,  5,  2,  9,

      // case 5, 1>2,3>4, orig start 1116, new start 412
      5,
      5,  7,  1,  8,
      5,  7,  0,  1,
      5,  7,  6,  0,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3, orig start 1137, new start 433
      5,
      0,  5,  6,  7,
      0,  5,  7,  8,
      0,  5,  8,  1,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3, a, orig start 1158, new start 454
      5,
      0,  5,  6,  8,
      0,  6,  7,  8,
      0,  5,  8,  1,
      5,  8,  9,  6,
      6,  7,  8,  9,

      // case 6_0, orig start 1319, new start 475
      4,
      7,  8,  9,  3,
      6,  5,  2,  9,
      4,  1,  5,  8,
      0,  4,  6,  7,

      // case 6_1, orig start 1336, new start 492
      4,
      6,  4,  5,  8,
      6,  5,  9,  8,
      6,  9,  7,  8,
      6,  7,  4,  8,

      // case 6_1, a, orig start 1353, new start 509
      4,
      5,  8,  9,  7,
      5,  9,  6,  7,
      5,  6,  4,  7,
      5,  4,  8,  7,

      // case 6_1, b, orig start 1370, new start 526
      4,
      4,  5,  6,  9,
      4,  6,  7,  9,
      4,  7,  8,  9,
      4,  8,  5,  9,

    };

} // namespace moab
} // namespace krino

