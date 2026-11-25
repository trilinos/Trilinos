// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETSWAPS_HPP_
#define KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETSWAPS_HPP_

#include <algorithm>
#include <array>
#include <vector>

#include <Akri_TetEdgeSwapConfigurations.hpp>
#include <stk_math/StkVector.hpp>

namespace krino
{

// Swaps of tets that replace the tets around an interior edge (edge swaps), or the two tets on a
// face (the 2-3 face swap), by others filling the same region.  NODE is any node handle (an index
// or a pointer); a coordinates function maps a node to its location.

double tet_volume(const std::array<stk::math::Vector3d, 4> & nodeCoords);
// The sign of the tet's volume, 0 if it is zero to roundoff: within 1e-12 of the cube of the
// tet's longest edge.  So a tet whose nodes are coplanar, or three of them collinear, is flat
// whatever the roundoff in its volume.
int tet_volume_sign(const std::array<stk::math::Vector3d, 4> & nodeCoords);
double tet_mean_ratio_quality(const std::array<stk::math::Vector3d, 4> & nodeCoords);

// The nodes of the tets around the edge, other than the edge's, in order around it.  Returns
// whether they form a closed ring, as for an interior edge: each tet has both edge nodes, and each
// of its other nodes is shared with one other tet.
template <typename NODE>
bool fill_ring_of_nodes_around_tet_edge(const std::vector<std::array<NODE, 4>> & tetsAroundEdge,
    const NODE edgeNode0,
    const NODE edgeNode1,
    std::vector<NODE> & ring)
{
  ring.clear();
  std::vector<std::array<NODE, 2>> otherNodes;
  for (auto & tet : tetsAroundEdge)
  {
    if (std::count(tet.begin(), tet.end(), edgeNode0) != 1 ||
        std::count(tet.begin(), tet.end(), edgeNode1) != 1)
      return false;
    std::array<NODE, 2> others{};
    int numOthers = 0;
    for (auto & node : tet)
      if (node != edgeNode0 && node != edgeNode1) others[numOthers++] = node;
    otherNodes.push_back(others);
  }
  const size_t numTets = otherNodes.size();
  if (numTets < 3) return false;

  std::vector<bool> used(numTets, false);
  ring = {otherNodes[0][0], otherNodes[0][1]};
  used[0] = true;
  for (size_t numUsed = 1; numUsed < numTets; ++numUsed)
  {
    bool found = false;
    for (size_t i = 0; i < numTets && !found; ++i)
    {
      if (used[i]) continue;
      for (int k = 0; k < 2 && !found; ++k)
      {
        if (otherNodes[i][k] == ring.back())
        {
          used[i] = found = true;
          if (numUsed + 1 < numTets)
            ring.push_back(otherNodes[i][1 - k]);
          else if (otherNodes[i][1 - k] != ring.front())
            return false;
        }
      }
    }
    if (!found) return false;
  }
  if (ring.size() != numTets) return false;
  std::vector<NODE> sortedRing = ring;
  std::sort(sortedRing.begin(), sortedRing.end());
  return std::adjacent_find(sortedRing.begin(), sortedRing.end()) == sortedRing.end();
}

// The number of edge swap configurations for a ring of this many nodes: 0 for fewer than 3, or for
// more than the tables have (TetEdgeSwapConfigurations).
size_t num_tet_edge_swap_configurations(const size_t numRingNodes);

// The tets of an edge swap configuration (TetEdgeSwapConfigurations) of the cavity of the tets
// around the edge with this ring of nodes: 2n-4 tets, none with the edge.  Their orientation
// depends on the direction of the ring (orient_tets_positively).
template <typename NODE>
std::vector<std::array<NODE, 4>> tets_of_tet_edge_swap_configuration(const std::vector<NODE> & ring,
    const NODE edgeNode0,
    const NODE edgeNode1,
    const size_t iConfig)
{
  std::vector<NODE> localNodes = ring;
  localNodes.push_back(edgeNode0);
  localNodes.push_back(edgeNode1);
  std::vector<std::array<NODE, 4>> tets;
  for (auto & tet :
      TetEdgeSwapConfigurations::get_possible_tet_configurations_for_num_tets_around_edge(
          ring.size())[iConfig])
    tets.push_back(
        {{localNodes[tet[0]], localNodes[tet[1]], localNodes[tet[2]], localNodes[tet[3]]}});
  return tets;
}

// The edges between ring nodes that an edge swap configuration makes.
template <typename NODE>
std::vector<std::array<NODE, 2>>
new_edges_of_tet_edge_swap_configuration(const std::vector<NODE> & ring, const size_t iConfig)
{
  std::vector<std::array<NODE, 2>> edges;
  for (auto & edge :
      TetEdgeSwapConfigurations::get_interior_edge_configurations_for_num_tets_around_edge(
          ring.size())[iConfig])
    edges.push_back({{ring[edge[0]], ring[edge[1]]}});
  return edges;
}

// The tets of the 2-3 swap of the tets with the face (faceNodes, p) and (faceNodes, q): the three
// tets around the edge p-q.  Their orientation depends on the order of the face nodes
// (orient_tets_positively).
template <typename NODE>
std::vector<std::array<NODE, 4>>
tets_of_tet_face_swap(const std::array<NODE, 3> & faceNodes, const NODE p, const NODE q)
{
  return {{{faceNodes[0], faceNodes[1], p, q}},
      {{faceNodes[1], faceNodes[2], p, q}},
      {{faceNodes[2], faceNodes[0], p, q}}};
}

// Orient the tets positively: all of them, if their volumes are all positive or all negative (as
// when the ring or face is in the other direction); returns false, leaving them, if the signs are
// mixed or a volume is zero to roundoff (tet_volume_sign), that is, if the tets do not fill the
// region they would replace.
template <typename NODE, typename COORDS>
bool orient_tets_positively(std::vector<std::array<NODE, 4>> & tets, const COORDS & coords)
{
  size_t numPositive = 0;
  size_t numNegative = 0;
  for (auto & tet : tets)
  {
    const int volumeSign =
        tet_volume_sign({{coords(tet[0]), coords(tet[1]), coords(tet[2]), coords(tet[3])}});
    if (volumeSign > 0)
      ++numPositive;
    else if (volumeSign < 0)
      ++numNegative;
  }
  if (numNegative == tets.size())
    for (auto & tet : tets)
      std::swap(tet[0], tet[1]);
  return numPositive == tets.size() || numNegative == tets.size();
}

template <typename NODE, typename COORDS>
double
worst_tet_mean_ratio_quality(const std::vector<std::array<NODE, 4>> & tets, const COORDS & coords)
{
  double worst = 1.;
  for (auto & tet : tets)
    worst = std::min(worst,
        tet_mean_ratio_quality({{coords(tet[0]), coords(tet[1]), coords(tet[2]), coords(tet[3])}}));
  return worst;
}

// Of the edge swap configurations of the cavity around the edge, the valid one (its tets,
// positively oriented, fill the cavity) that accept(tets, iConfig) accepts with the best worst
// mean-ratio quality, or -1 if there is none (as for a ring with more nodes than the tables have).
template <typename NODE, typename COORDS, typename ACCEPT>
int best_tet_edge_swap_configuration(const std::vector<NODE> & ring,
    const NODE edgeNode0,
    const NODE edgeNode1,
    const COORDS & coords,
    const ACCEPT & accept)
{
  int bestConfig = -1;
  double bestQuality = 0.;
  for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(ring.size()); ++iConfig)
  {
    std::vector<std::array<NODE, 4>> tets =
        tets_of_tet_edge_swap_configuration(ring, edgeNode0, edgeNode1, iConfig);
    if (!orient_tets_positively(tets, coords) || !accept(tets, iConfig)) continue;
    const double quality = worst_tet_mean_ratio_quality(tets, coords);
    if (bestConfig < 0 || quality > bestQuality)
    {
      bestConfig = iConfig;
      bestQuality = quality;
    }
  }
  return bestConfig;
}

} // namespace krino

#endif /* KRINO_KRINO_ELEMENT_IMPROVEMENT_AKRI_TETSWAPS_HPP_ */
