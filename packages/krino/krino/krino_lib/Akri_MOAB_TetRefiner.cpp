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

#include "Akri_MOAB_TetRefiner.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stack>
#include <string>
#include <stdexcept>

#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_DiagWriter.hpp>

namespace krino {
namespace moab {

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
        ThrowRuntimeError("Case " << caseId << " not supported in num_new_child_elements_tet4.");
  }
}

  //p this is from some comments below, and coded in tet_edges.  It's only used for disambiguation of same-length edges.
  //p   * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
  int tet_edges[6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};

  double compute_edge_length_squared(const double *c0, const double *c1)
  {
    return
      (c0[0]-c1[0])*(c0[0]-c1[0])+
      (c0[1]-c1[1])*(c0[1]-c1[1])+
      (c0[2]-c1[2])*(c0[2]-c1[2]);
  }

  inline double SQR(const double x) { return x*x; }

  /// returns tet quality:
  ///   if use_best_quality, then return max edge len/min - 1.0 is ideal, smaller quality is better
  ///   else, return the max edge length (consistent with standard practice heuristics of choosing
  ///      the octagon's shortest diagonal)
  inline double quality(const int *indices, const double *coords[14])
  {
    double edge_min=std::numeric_limits<double>::max();
    double edge_max = 0;
    for (int i=0; i < 3; i++)
      {
        for (int j=i+1; j < 4; j++)
          {
            const double *ci = coords[indices[i]];
            const double *cj = coords[indices[j]];
            const double el2 = SQR(ci[0]-cj[0]) + SQR(ci[1]-cj[1]) + SQR(ci[2]-cj[2]) ;
            edge_min = std::min(edge_min, el2);
            edge_max = std::max(edge_max, el2);
          }
      }
    return std::sqrt(edge_max/edge_min);
  }

  int SimplexTemplateRefiner::best_tets( const int* alternates, const double* coords[14], int, int )
  {
    // force better tets to be at least this factor better - helps with tie breaks for structured meshes
    const double factor = 0.95;
    int nalt = -1;
    for (int i=0; i < 100; i++)
      {
        if (alternates[i] < 0) {
          nalt = i;
          break;
        }
      }
    if (nalt < 0) throw std::runtime_error("SimplexTemplateRefiner::best_tets");
    double best_qual=0;
    int iqual=-1;
    for (int i=0; i < nalt; i++)
      {
        int * current_template = SimplexTemplateRefiner::templates + alternates[i];
        // find worst quality element
        double max_qual=0;
        for (int j=0; j < current_template[0]; j++)
          {
            max_qual = std::max(max_qual, quality(current_template + 1 + j*4, coords));
          }
        // find alternates with the best (min) worst quality
        if (i == 0)
          {
            best_qual = max_qual;
            iqual = 0;
          }
        else if (max_qual < best_qual*factor)
          {
            best_qual = max_qual;
            iqual = i;
          }
      }
    return alternates[iqual];
  }

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
      ThrowRequire(1 == nodeSideIntersection.size());
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

  std::vector<SimplexTemplateRefiner::TetDescription> SimplexTemplateRefiner::refinement_child_nodes_and_sides_tet4(const unsigned encodedEdgesToRefine, const std::array<stk::math::Vector3d,4> & elementNodeCoords, const std::array<int,4> & elementNodeScore, const bool needSides)
  {
    const std::vector<std::array<int,4>> newTetNodes = moab::SimplexTemplateRefiner::refinement_child_nodes_tet4(encodedEdgesToRefine, elementNodeCoords, elementNodeScore );

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
          newTets[iTet].sideIds[iSide] = -1; // -1 indices side is not on any parent side
      }
    }

    return newTets;
  }

  /**\brief Refine a tetrahedron.

   */
  std::vector<TetTupleInt> SimplexTemplateRefiner::refinement_child_nodes_tet4(const unsigned edge_code, const std::array<stk::math::Vector3d,4> & elementNodeCoords, const std::array<int,4> & elementNodeScore )
  {  // refine_3_simplex

    const double* v0 = elementNodeCoords[0].data();
    const double* v1 = elementNodeCoords[1].data();
    const double* v2 = elementNodeCoords[2].data();
    const double* v3 = elementNodeCoords[3].data();
    const EntityHandle h0 = elementNodeScore[0];
    const EntityHandle h1 = elementNodeScore[1];
    const EntityHandle h2 = elementNodeScore[2];
    const EntityHandle h3 = elementNodeScore[3];

    std::vector<TetTupleInt> new_tets;

    double s_midpt0c[3]={0,0,0};
    double s_midpt1c[3]={0,0,0};
    double s_midpt2c[3]={0,0,0};
    double s_midpt3c[3]={0,0,0};
    double s_midpt4c[3]={0,0,0};
    double s_midpt5c[3]={0,0,0};

    double* midpt0c=&s_midpt0c[0];
    double* midpt1c=&s_midpt1c[0];
    double* midpt2c=&s_midpt2c[0];
    double* midpt3c=&s_midpt3c[0];
    double* midpt4c=&s_midpt4c[0];
    double* midpt5c=&s_midpt5c[0];

    EntityHandle midpt0h=0;
    EntityHandle midpt1h=0;
    EntityHandle midpt2h=0;
    EntityHandle midpt3h=0;
    EntityHandle midpt4h=0;
    EntityHandle midpt5h=0;

    for ( int i = 0; i < 3; ++ i )
      {
        midpt0c[i] = ( v0[i] + v1[i] ) * .5;
        midpt1c[i] = ( v1[i] + v2[i] ) * .5;
        midpt2c[i] = ( v2[i] + v0[i] ) * .5;
        midpt3c[i] = ( v0[i] + v3[i] ) * .5;
        midpt4c[i] = ( v1[i] + v3[i] ) * .5;
        midpt5c[i] = ( v2[i] + v3[i] ) * .5;
      }

    if ( ! edge_code )
      {
        // No edges to subdivide
        return new_tets;
      }

    double facept0c[3]={0,0,0};
    double facept1c[3]={0,0,0};
    double facept2c[3]={0,0,0};
    double facept3c[3]={0,0,0};

    const double* vertex_coords[14] = {
      v0, v1, v2, v3,
      midpt0c, midpt1c, midpt2c,
      midpt3c, midpt4c, midpt5c,
      facept0c, facept1c, facept2c, facept3c
    };

    EntityHandle vertex_hash[14] = {
      h0, h1, h2, h3,
      midpt0h, midpt1h, midpt2h,
      midpt3h, midpt4h, midpt5h,
      0, 0, 0, 0
    };

    // Generate tetrahedra that are compatible except when edge
    // lengths are equal on indeterminately subdivided faces.
    const double* permuted_coords[14];
    EntityHandle permuted_hash[14];
    int permuted_local_ids[14];   // this is just a copy of permutations_from_index, for readability
    double permlen[6]; // permuted edge lengths
    int C = SimplexTemplateRefiner::template_index[edge_code][0];
    int P = SimplexTemplateRefiner::template_index[edge_code][1];

    // 1. Permute the tetrahedron into our canonical configuration
    for ( int i = 0; i < 14; ++ i )
      {
        permuted_local_ids[i] = SimplexTemplateRefiner::permutations_from_index[P][i];
        permuted_coords[i] = vertex_coords[SimplexTemplateRefiner::permutations_from_index[P][i]];
        permuted_hash[i] = vertex_hash[SimplexTemplateRefiner::permutations_from_index[P][i]];
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
        output_tets.push( SimplexTemplateRefiner::templates + 40 );
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

        output_tets.push( SimplexTemplateRefiner::templates + 57 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        // here
        switch ( comparison_bits )
          {
          case 42: // 0>2>3<0
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 25: // 2>3>0<2
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[11] );
            output_sign.push( 1 );
            break;
          case 37: // 3>0>2<3
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[3] );
            output_sign.push( 1 );
            break;
          case 21: // 3>2>0<3
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[22] );
            output_sign.push( -1 );
            break;
          case 26: // 2>0>3<2
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[12] );
            output_sign.push( -1 );
            break;
          case 38: // 0>3>2<0
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 5: // Ruprecht-Muller Case 3b
        output_tets.push( SimplexTemplateRefiner::templates + 162 );
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
            output_tets.push( SimplexTemplateRefiner::templates + 179 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 5: // 1>0,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 200 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 6: // 0>1,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 221 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 9: // 1>0,0>3
            output_tets.push( SimplexTemplateRefiner::templates + 242 );
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
            output_tets.push( SimplexTemplateRefiner::templates + 362 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 5: // 4>0,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 383 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 0>4,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 404 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 6: // 4>0,0>2
            output_tets.push( SimplexTemplateRefiner::templates + 425 );
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

        output_tets.push( SimplexTemplateRefiner::templates + 545 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        switch ( comparison_bits )
          {
          case 5: // 5>4>3
            output_tets.push( SimplexTemplateRefiner::templates + 554 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 10: // 3>4>5
            output_tets.push( SimplexTemplateRefiner::templates + 554 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );

            break;
          case 6: // 3<4>5
            output_tets.push( SimplexTemplateRefiner::templates + 571 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 3>4<5
            output_tets.push( SimplexTemplateRefiner::templates + 588 );
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
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 102: // 1>2,3>2,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            break;
          case 170: // 1>2,2>3,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          case 153: // 2>1,2>3,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            break;
          case 90: // 1>2,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            break;
          case 105: // 2>1,2>3,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            break;
          case 165: // 2>1,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            break;
          case 150: // 1>2,3>2,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            break;
          case 101: // 2>1,3>2,3>4,4>1
            {
              int alternates[] = { 713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 0, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 86: // 1>2,3>2,4>3,4>1
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 14, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            break;
          case 154: // 1>2,2>3,4>3,1>4
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 5, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            break;
          case 169: // 2>1,2>3,3>4,1>4
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 15, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            break;
          case 89: // 2>1,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 763 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            break;
          case 166: // 1>2,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 763 );
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

        output_tets.push( SimplexTemplateRefiner::templates + 1107 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        switch ( comparison_bits )
          {
          case 10: // 1>2,3>4
            output_tets.push( SimplexTemplateRefiner::templates + 1116 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 5: // 2>1,4>3
            output_tets.push( SimplexTemplateRefiner::templates + 1116 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );

            break;
          case 6: // 1>2,4>3
            {
              int alternates[] = { 1137, 1158, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 0, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );

            break;
          case 9: // 2>1,3>4
            {
              int alternates[] = {1137, 1158, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 14, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );

            break;
          default:
            throw std::runtime_error("Invalid case");
          }
        break;
      case 11: // Ruprecht-Muller Case 6
        output_tets.push( SimplexTemplateRefiner::templates + 1319 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );

        {
          int alternates[] = { 1336, 1353, 1370, -1 };
          output_tets.push( SimplexTemplateRefiner::templates + best_tets( alternates, permuted_coords, 0, 1 ) );
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
            // we have an inverted tet... reverse the first 2 vertices
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

    ThrowRequireMsg(num_new_child_elements_tet4(edge_code) == new_tets.size(), "Mismatch of size " << num_new_child_elements_tet4(edge_code) << "  " << new_tets.size() << " for case " << edge_code); // FIXME: remove check

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
   * permutation index, it returns a row of 14 values:
   * these are the vertex numbers of the permuted
   * tetrahedron. The first 4 values are the permuted
   * corner indices, the next 6 values are the
   * permuted edge midpoint indices, and the final
   * entries reference mid-face points inserted
   * to maintain a compatible tetrahedralization.
   *
   * There are 24 entries, 6 for each of the 4 faces of
   * the tetrahedron.
   */
  int SimplexTemplateRefiner::permutations_from_index[24][14] =
    {
      /* corners      midpoints          face points   */
      /* POSITIVE ARRANGEMENTS                         */
      { 0, 1, 2, 3,   4, 5, 6, 7, 8, 9,  10, 11, 12, 13 }, /* Face 0-1-2 */
      { 1, 2, 0, 3,   5, 6, 4, 8, 9, 7,  10, 12, 13, 11 },
      { 2, 0, 1, 3,   6, 4, 5, 9, 7, 8,  10, 13, 11, 12 },

      { 0, 3, 1, 2,   7, 8, 4, 6, 9, 5,  11, 13, 12, 10 }, /* Face 0-3-1 */
      { 3, 1, 0, 2,   8, 4, 7, 9, 5, 6,  11, 12, 10, 13 },
      { 1, 0, 3, 2,   4, 7, 8, 5, 6, 9,  11, 10, 13, 12 },

      { 1, 3, 2, 0,   8, 9, 5, 4, 7, 6,  12, 11, 13, 10 }, /* Face 1-3-2 */
      { 3, 2, 1, 0,   9, 5, 8, 7, 6, 4,  12, 13, 10, 11 },
      { 2, 1, 3, 0,   5, 8, 9, 6, 4, 7,  12, 10, 11, 13 },

      { 2, 3, 0, 1,   9, 7, 6, 5, 8, 4,  13, 12, 11, 10 }, /* Face 2-3-0 */
      { 3, 0, 2, 1,   7, 6, 9, 8, 4, 5,  13, 11, 10, 12 },
      { 0, 2, 3, 1,   6, 9, 7, 4, 5, 8,  13, 10, 12, 11 },

      /* NEGATIVE ARRANGEMENTS                         */
      { 0, 2, 1, 3,   6, 5, 4, 7, 9, 8,  10, 13, 12, 11 }, /* Face 0-1-2 */
      { 2, 1, 0, 3,   5, 4, 6, 9, 8, 7,  10, 12, 11, 13 },
      { 1, 0, 2, 3,   4, 6, 5, 8, 7, 9,  10, 11, 13, 12 },

      { 0, 1, 3, 2,   4, 8, 7, 6, 5, 9,  11, 10, 12, 13 }, /* Face 0-3-1 */
      { 1, 3, 0, 2,   8, 7, 4, 5, 9, 6,  11, 12, 13, 10 },
      { 3, 0, 1, 2,   7, 4, 8, 9, 6, 5,  11, 13, 10, 12 },

      { 1, 2, 3, 0,   5, 9, 8, 4, 6, 7,  12, 10, 13, 11 }, /* Face 1-3-2 */
      { 2, 3, 1, 0,   9, 8, 5, 6, 7, 4,  12, 13, 11, 10 },
      { 3, 1, 2, 0,   8, 5, 9, 7, 4, 6,  12, 11, 10, 13 },

      { 2, 0, 3, 1,   6, 7, 9, 5, 4, 8,  13, 10, 11, 12 }, /* Face 2-3-0 */
      { 0, 3, 2, 1,   7, 9, 6, 4, 8, 5,  13, 11, 12, 10 },
      { 3, 2, 0, 1,   9, 6, 7, 8, 5, 4,  13, 12, 10, 11 }
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
      // case 1_0
      2,
      0,  4,  2,  3,
      4,  1,  2,  3,

      // case 2a_0
      1,
      3,  4,  5,  1,

      // case 2a, 0>1
      2,
      0,  4,  2,  3,
      4,  5,  2,  3,

      // case 2a, 0=1
      4,
      10,  3,  0,  4,
      10,  3,  4,  5,
      10,  3,  5,  2,
      10,  3,  2,  0,

      // case 2b_0
      4,
      0,  4,  9,  3,
      4,  1,  9,  3,
      0,  4,  2,  9,
      4,  1,  2,  9,

      // case 3a_0
      1,
      4,  7,  6,  0,

      // case 3a, 0>2>3<0
      3,
      1,  3,  2,  4,
      4,  6,  3,  2,
      4,  6,  7,  3,

      // case 3a, 0=2>3<0
      5,
      4,  6,  7,  3,
      10,  1,  2,  3,
      10,  2,  6,  3,
      10,  6,  4,  3,
      10,  4,  1,  3,

      // case 3a, 3>0=2<3
      5,
      1,  3,  2,  7,
      10,  1,  2,  7,
      10,  2,  6,  7,
      10,  6,  4,  7,
      10,  4,  1,  7,

      // case 3a, 0=2=3=0
      11,
      2,  6, 10, 13,
      3,  7, 13, 11,
      4,  1, 10, 11,
      11,  6, 10,  4,
      11,  6, 13, 10,
      11,  6,  7, 13,
      11,  6,  4,  7,
      2, 10, 11, 13,
      1, 10, 11,  2,
      2, 11,  3, 13,
      3,  2,  1, 11,

      // case 3b_0
      4,
      0,  7,  4,  2,
      4,  7,  8,  2,
      4,  8,  1,  2,
      7,  3,  8,  2,

      // case 3c, 0>1,0>3
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3c, 1>0,3>0
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 0>1,3>0
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 1>0,0>3
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3c, 0=1,0>3
      7,
      4,  1,  5,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  7,  4,  3,
      10,  2,  7,  3,
      10,  5,  2,  3,
      10,  4,  5,  3,

      // case 3c, 3>0,0=1
      7,
      7,  1,  5,  3,
      7,  5,  2,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  5,  2,  7,
      10,  4,  5,  7,
      1,  5,  4,  7,

      // case 3c, 0=1,0=3
      10,
      4,  1,  5, 11,
      11,  1,  5,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  5,  2,  3,
      10,  2,  7,  3,
      10,  7,  4, 11,
      10,  7, 11,  3,
      10,  4,  5, 11,
      10, 11,  5,  3,

      // case 3d, 0>4,0>2
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,2>0
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 0>4,2>0
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,0>2
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 3d, 0=4,0>2
      7,
      4,  1,  2,  8,
      11,  4,  0,  6,
      11,  0,  3,  6,
      11,  2,  4,  6,
      11,  3,  2,  6,
      11,  3,  8,  2,
      11,  8,  4,  2,

      // case 3d, 2>0,0=4
      7,
      6,  2,  8,  1,
      6,  8,  2,  3,
      11,  4,  0,  6,
      11,  0,  3,  6,
      8, 11,  3,  6,
      8,  4, 11,  6,
      1,  6,  4,  8,

      // case 3d, 0=4,0=2
      10,
      4,  1, 10,  8,
      10,  2,  8,  1,
      11,  4,  0,  6,
      11,  0,  3,  6,
      11,  3,  8,  2,
      11,  3,  2,  6,
      11, 10,  4,  6,
      11, 10,  6,  2,
      8,  4, 11, 10,
      11, 10,  2,  8,

      // case 4a_0
      2,
      7,  8,  9,  3,
      7,  9,  8,  6,

      // case 4a, 5>4>3
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      9,  1,  6,  2,
      9,  1,  8,  6,

      // case 4a, 3<4>5
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      8,  2,  6,  9,
      8,  2,  1,  6,

      // case 4a, 3>4<5
      4,
      6,  9,  8,  1,
      6,  9,  1,  2,
      6,  7,  0,  1,
      6,  7,  1,  8,

      // case 4a, 3=4>5
      6,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      1,  2,  6,  8,
      2,  6,  8,  9,

      // case 4a, 5>4,3=4
      6,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      1,  2,  6,  9,
      1,  6,  8,  9,

      // case 4a, 3=4=5
      8,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      6,  1,  2, 12,
      6,  2,  9, 12,
      6,  9,  8, 12,
      6,  8,  1, 12,

      // case 4b, 2>1,3>2,4>3,4>1
      6,
      6,  8,  1,  5,
      6,  8,  0,  1,
      6,  8,  7,  0,
      6,  8,  2,  7,
      7,  8,  2,  3,
      6,  8,  5,  2,

      // case 4b, 2>1,3>2,3>4,4>1
      6,
      6,  8,  1,  5,
      6,  8,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      6,  8,  5,  2,
      6,  8,  2,  7,

      // case 4b, 2>1,3>2,3>4,4>1, a
      6,
      7,  8,  1,  5,
      6,  5,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      7,  8,  5,  2,
      6,  5,  2,  7,

      // case 4b, 2>1,2>3,4>3,4>1
      6,
      6,  8,  5,  2,
      6,  8,  2,  3,
      6,  8,  3,  7,
      6,  8,  7,  0,
      6,  8,  0,  1,
      6,  8,  1,  5,

      // case 4b, 1=2,3>2,3>4,4>1
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  2,  7,  8,
      6,  5,  2,  8,
      7,  8,  2,  3,

      // case 4b, 1=2,2>3,4>3,1>4
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  3,  7,  8,
      6,  5,  3,  8,
      6,  5,  2,  3,

      // case 4b, 1=2,2>3,4>3,4>1
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  3,  7,  8,
      6,  5,  2,  8,
      6,  2,  3,  8,

      // case 4b, 1=2,3>2,3=4,4>1
      11,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  2,  7,  8,
      6,  5,  2,  8,
      7,  8,  2,  3,

      // case 4b, 4>1=2=3,4>3
      12,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3,  8,
      13,  2,  8,  5,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6, 13,  7,  8,
      6,  5, 13,  8,

      // case 4b, 1=2=3>4,1>4
      12,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3,  5,
      13,  3,  8,  5,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6, 13,  7,  8,
      6,  5, 13,  8,

      // case 4b, 1=2=3=4=1
      16,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3, 12,
      13,  2, 12,  5,
      13, 12,  3,  8,
      13, 12,  5,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  5, 13,  8,
      6, 13,  7,  8,

      // case 5_0
      2,
      7,  8,  9,  3,
      6,  5,  2,  9,

      // case 5, 1>2,3>4
      5,
      5,  7,  1,  8,
      5,  7,  0,  1,
      5,  7,  6,  0,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3
      5,
      0,  5,  6,  7,
      0,  5,  7,  8,
      0,  5,  8,  1,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3, a
      5,
      0,  5,  6,  8,
      0,  6,  7,  8,
      0,  5,  8,  1,
      5,  8,  9,  6,
      6,  7,  8,  9,

      // case 5, 1=2,3>4
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      10,  8,  5,  9,
      10,  6,  7,  9,
      10,  7,  8,  9,
      10,  5,  6,  9,

      // case 5, 1=2,3>4, a
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      7,  8,  5,  9,
      10,  6,  7,  5,
      10,  7,  8,  5,
      5,  9,  6,  7,

      // case 5, 1=2,3>4, b
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      6,  8,  5,  9,
      10,  6,  7,  8,
      10,  6,  8,  5,
      8,  9,  6,  7,

      // case 5, 1=2,3=4
      10,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      10,  8,  5,  9,
      10,  6,  7,  9,
      10,  7,  8,  9,
      10,  5,  6,  9,

      // case 6_0
      4,
      7,  8,  9,  3,
      6,  5,  2,  9,
      4,  1,  5,  8,
      0,  4,  6,  7,

      // case 6_1
      4,
      6,  4,  5,  8,
      6,  5,  9,  8,
      6,  9,  7,  8,
      6,  7,  4,  8,

      // case 6_1, a
      4,
      5,  8,  9,  7,
      5,  9,  6,  7,
      5,  6,  4,  7,
      5,  4,  8,  7,

      // case 6_1, b
      4,
      4,  5,  6,  9,
      4,  6,  7,  9,
      4,  7,  8,  9,
      4,  8,  5,  9,

    };

} // namespace moab
} // namespace krino

