// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETEDGESWAPCONFIGURATIONS_HPP_
#define KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETEDGESWAPCONFIGURATIONS_HPP_

#include <array>
#include <cstddef>
#include <vector>

#include <stk_util/util/ReportHandler.hpp>

namespace krino
{

// The ways to re-triangulate the cavity of the tets around an interior edge without that edge, as
// in emend (ngs/emend/emend_lib/ElementEdgeSwapConfigurations.hpp, whose tables are copied here).
// For an edge with n nodes around it in a closed ring, the local node indices are the ring nodes
// 0..n-1, in order around the edge, then the edge's nodes n and n+1.  Each configuration is 2n-4
// tets, with its new edges between ring nodes.
class TetEdgeSwapConfigurations
{
public:
  typedef std::vector<std::array<int, 2>> EdgeGroup;
  typedef std::vector<std::vector<std::vector<int>>> TetsCandidatesArray;
  typedef std::vector<EdgeGroup> InteriorEdgesForTetCandidatesArray;

  static const TetsCandidatesArray &
  get_possible_tet_configurations_for_num_tets_around_edge(const size_t numNodesAroundEdge)
  {
    STK_ThrowRequire(
        numNodesAroundEdge <= get_max_num_nodes_around_edge_currently_supported_for_swapping());
    return allPossibleTetConfigurationsForSwappingEdge[numNodesAroundEdge];
  }

  static const InteriorEdgesForTetCandidatesArray &
  get_interior_edge_configurations_for_num_tets_around_edge(const size_t numNodesAroundEdge)
  {
    STK_ThrowRequire(
        numNodesAroundEdge <= get_max_num_nodes_around_edge_currently_supported_for_swapping());
    return interiorEdgesForConfigurationsForSwappingEdge[numNodesAroundEdge];
  }

  static size_t get_max_num_nodes_around_edge_currently_supported_for_swapping()
  {
    return allPossibleTetConfigurationsForSwappingEdge.size() - 1;
  }

  static constexpr int maxNumHaloNodesInEdgeSwapCavity{8};
  static const std::array<TetsCandidatesArray, maxNumHaloNodesInEdgeSwapCavity>
      allPossibleTetConfigurationsForSwappingEdge;
  static const std::array<std::vector<int>, maxNumHaloNodesInEdgeSwapCavity>
      allPossibleTetConfigurationsIndicesForSwappingEdge;
  static const std::array<InteriorEdgesForTetCandidatesArray, maxNumHaloNodesInEdgeSwapCavity>
      interiorEdgesForConfigurationsForSwappingEdge;
};

} // namespace krino

#endif /* KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETEDGESWAPCONFIGURATIONS_HPP_ */
