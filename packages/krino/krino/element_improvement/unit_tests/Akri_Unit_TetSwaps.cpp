// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_TetEdgeSwapConfigurations.hpp>
#include <Akri_TetSwaps.hpp>
#include <stk_math/StkVector.hpp>

namespace krino
{

namespace
{

// Nodes 0..n-1 on a ring around the edge from node n, at (0,0,-1), to node n+1, at (0,0,1).
class RingAroundEdge
{
public:
  explicit RingAroundEdge(const std::vector<stk::math::Vector3d> & ringCoords)
      : myCoords(ringCoords)
  {
    myCoords.emplace_back(0., 0., -1.);
    myCoords.emplace_back(0., 0., 1.);
  }

  static RingAroundEdge regular(const int n, const double radiusX = 1., const double radiusY = 1.)
  {
    std::vector<stk::math::Vector3d> ringCoords;
    for (int i = 0; i < n; ++i)
      ringCoords.emplace_back(
          radiusX * std::cos(2. * M_PI * i / n), radiusY * std::sin(2. * M_PI * i / n), 0.);
    return RingAroundEdge(ringCoords);
  }

  int num_ring_nodes() const { return myCoords.size() - 2; }
  int edge_node0() const { return num_ring_nodes(); }
  int edge_node1() const { return num_ring_nodes() + 1; }
  std::vector<int> ring() const
  {
    std::vector<int> nodes;
    for (int i = 0; i < num_ring_nodes(); ++i)
      nodes.push_back(i);
    return nodes;
  }
  std::vector<std::array<int, 4>> tets_around_edge() const
  {
    std::vector<std::array<int, 4>> tets;
    for (int i = 0; i < num_ring_nodes(); ++i)
      tets.push_back({{i, (i + 1) % num_ring_nodes(), edge_node0(), edge_node1()}});
    return tets;
  }
  const stk::math::Vector3d & coords(const int node) const { return myCoords[node]; }
  void set_coords(const int node, const stk::math::Vector3d & x) { myCoords[node] = x; }
  auto coords_function() const
  {
    return [this](const int node) { return coords(node); };
  }

private:
  std::vector<stk::math::Vector3d> myCoords;
};

template <typename COORDS>
double total_abs_volume(const std::vector<std::array<int, 4>> & tets, const COORDS & coords)
{
  double volume = 0.;
  for (auto & tet : tets)
    volume +=
        std::abs(tet_volume({{coords(tet[0]), coords(tet[1]), coords(tet[2]), coords(tet[3])}}));
  return volume;
}

bool has_both(const std::array<int, 4> & tet, const int a, const int b)
{
  return std::count(tet.begin(), tet.end(), a) == 1 && std::count(tet.begin(), tet.end(), b) == 1;
}

bool has_edge(const std::vector<std::array<int, 2>> & edges, const int a, const int b)
{
  return std::any_of(edges.begin(),
      edges.end(),
      [&](const std::array<int, 2> & e)
      { return (e[0] == a && e[1] == b) || (e[0] == b && e[1] == a); });
}

const std::array<size_t, 8> expectedNumConfigurations{{0, 0, 0, 1, 2, 5, 14, 42}};

std::pair<int, int> sorted_edge(const int a, const int b)
{
  return {std::min(a, b), std::max(a, b)};
}

bool is_closed_ring(const std::vector<std::array<int, 4>> & tets, const int e0, const int e1)
{
  std::vector<int> ring;
  return fill_ring_of_nodes_around_tet_edge(tets, e0, e1, ring);
}

} // namespace

TEST(TetEdgeSwapConfigurations, allRingSizes_configurationsOfTheCavityWithoutTheEdge)
{
  ASSERT_EQ(7u,
      TetEdgeSwapConfigurations::get_max_num_nodes_around_edge_currently_supported_for_swapping());
  for (size_t n = 3; n <= 7; ++n)
  {
    const auto & configs =
        TetEdgeSwapConfigurations::get_possible_tet_configurations_for_num_tets_around_edge(n);
    const auto & newEdges =
        TetEdgeSwapConfigurations::get_interior_edge_configurations_for_num_tets_around_edge(n);
    ASSERT_EQ(expectedNumConfigurations[n], configs.size()) << "n " << n;
    ASSERT_EQ(configs.size(), newEdges.size()) << "n " << n;
    EXPECT_EQ(expectedNumConfigurations[n], num_tet_edge_swap_configurations(n));
    for (size_t iConfig = 0; iConfig < configs.size(); ++iConfig)
    {
      EXPECT_EQ(2 * n - 4, configs[iConfig].size()) << "n " << n << " config " << iConfig;
      EXPECT_EQ(n - 3, newEdges[iConfig].size()) << "n " << n << " config " << iConfig;
      for (auto & tet : configs[iConfig])
      {
        ASSERT_EQ(4u, tet.size());
        for (int node : tet)
        {
          EXPECT_TRUE(node >= 0 && node < int(n) + 2) << "n " << n << " config " << iConfig;
        }
        EXPECT_FALSE(has_both({{tet[0], tet[1], tet[2], tet[3]}}, n, n + 1))
            << "n " << n << " config " << iConfig;
      }
      for (auto & edge : newEdges[iConfig])
      {
        EXPECT_TRUE(edge[0] < int(n) && edge[1] < int(n)) << "n " << n << " config " << iConfig;
      }
    }
  }
}

// The new edges listed for each configuration are the edges between ring nodes that its tets have,
// other than those around the ring.
TEST(TetEdgeSwapConfigurations, allRingSizes_newEdgesAreThoseOfTheTets)
{
  for (size_t n = 3; n <= 7; ++n)
  {
    const auto & configs =
        TetEdgeSwapConfigurations::get_possible_tet_configurations_for_num_tets_around_edge(n);
    const auto & newEdges =
        TetEdgeSwapConfigurations::get_interior_edge_configurations_for_num_tets_around_edge(n);
    for (size_t iConfig = 0; iConfig < configs.size(); ++iConfig)
    {
      std::set<std::pair<int, int>> edgesOfTets;
      for (auto & tet : configs[iConfig])
        for (int i = 0; i < 4; ++i)
          for (int j = i + 1; j < 4; ++j)
            if (tet[i] < int(n) && tet[j] < int(n)) edgesOfTets.insert(sorted_edge(tet[i], tet[j]));
      for (int i = 0; i < int(n); ++i)
        edgesOfTets.erase(sorted_edge(i, (i + 1) % n));

      std::set<std::pair<int, int>> listedEdges;
      for (auto & edge : newEdges[iConfig])
        listedEdges.insert(sorted_edge(edge[0], edge[1]));
      EXPECT_EQ(n - 3, listedEdges.size()) << "n " << n << " config " << iConfig;
      EXPECT_EQ(edgesOfTets, listedEdges) << "n " << n << " config " << iConfig;
    }
  }
}

TEST(TetEdgeSwapConfigurations, unsupportedRingSizes_noConfigurations)
{
  for (size_t n : {0u, 1u, 2u, 8u, 20u})
  {
    EXPECT_EQ(0u, num_tet_edge_swap_configurations(n)) << "n " << n;
  }
}

TEST(TetSwaps, tetsAroundInteriorEdgeInAnyOrder_closedRingInOrder)
{
  for (int n = 3; n <= 7; ++n)
  {
    const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(n);
    std::vector<std::array<int, 4>> tets = ringAroundEdge.tets_around_edge();
    std::reverse(tets.begin(), tets.end());
    std::rotate(tets.begin(), tets.begin() + 1, tets.end());
    for (const bool flipEdge : {false, true})
    {
      const int e0 = flipEdge ? ringAroundEdge.edge_node1() : ringAroundEdge.edge_node0();
      const int e1 = flipEdge ? ringAroundEdge.edge_node0() : ringAroundEdge.edge_node1();
      std::vector<int> ring;
      ASSERT_TRUE(fill_ring_of_nodes_around_tet_edge(tets, e0, e1, ring)) << "n " << n;
      ASSERT_EQ(size_t(n), ring.size()) << "n " << n;
      for (int i = 0; i < n; ++i)
      {
        const int step = (ring[(i + 1) % n] - ring[i] + n) % n;
        EXPECT_TRUE(step == 1 || step == n - 1) << "n " << n;
      }
    }
  }
}

TEST(TetSwaps, tetsAroundEdgeWithOneMissing_notClosedRing)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(5);
  std::vector<std::array<int, 4>> tets = ringAroundEdge.tets_around_edge();
  tets.erase(tets.begin() + 2);
  std::vector<int> ring;
  EXPECT_FALSE(fill_ring_of_nodes_around_tet_edge(
      tets, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), ring));
}

TEST(TetSwaps, tetWithoutAnEdgeNode_notClosedRing)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(5);
  std::vector<std::array<int, 4>> tets = ringAroundEdge.tets_around_edge();
  tets[2][3] = 4;
  EXPECT_FALSE(is_closed_ring(tets, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1()));
}

// Ring node 1 is in three tets around the edge.
TEST(TetSwaps, ringNodeInThreeTets_notClosedRing)
{
  const int e0 = 5;
  const int e1 = 6;
  const std::vector<std::array<int, 4>> tets{
      {{0, 1, e0, e1}}, {{1, 2, e0, e1}}, {{1, 3, e0, e1}}, {{3, 4, e0, e1}}, {{4, 0, e0, e1}}};
  EXPECT_FALSE(is_closed_ring(tets, e0, e1));
}

// Each tet twice: the ring would wind around the edge twice.
TEST(TetSwaps, repeatedTets_notClosedRing)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(3);
  std::vector<std::array<int, 4>> tets = ringAroundEdge.tets_around_edge();
  const std::vector<std::array<int, 4>> once = tets;
  tets.insert(tets.end(), once.begin(), once.end());
  EXPECT_FALSE(is_closed_ring(tets, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1()));
}

TEST(TetSwaps, twoTetsAroundEdge_notClosedRing)
{
  const int e0 = 2;
  const int e1 = 3;
  EXPECT_FALSE(is_closed_ring({{{0, 1, e0, e1}}, {{1, 0, e0, e1}}}, e0, e1));
}

// The ring closes, but the tables stop at 7 nodes, so there is no configuration to choose.
TEST(TetSwaps, ringOfMoreNodesThanTheTables_noConfiguration)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(8);
  std::vector<int> ring;
  ASSERT_TRUE(fill_ring_of_nodes_around_tet_edge(ringAroundEdge.tets_around_edge(),
      ringAroundEdge.edge_node0(),
      ringAroundEdge.edge_node1(),
      ring));
  EXPECT_EQ(8u, ring.size());
  EXPECT_EQ(-1,
      best_tet_edge_swap_configuration(ring,
          ringAroundEdge.edge_node0(),
          ringAroundEdge.edge_node1(),
          ringAroundEdge.coords_function(),
          [](const std::vector<std::array<int, 4>> &, size_t) { return true; }));
}

// Around a convex ring, every configuration fills the cavity, in either direction of the ring.
TEST(TetSwaps, convexRing_everyConfigurationFillsTheCavity)
{
  for (int n = 3; n <= 7; ++n)
  {
    const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(n);
    const auto coords = ringAroundEdge.coords_function();
    const double cavityVolume = total_abs_volume(ringAroundEdge.tets_around_edge(), coords);
    for (const bool reverseRing : {false, true})
    {
      std::vector<int> ring = ringAroundEdge.ring();
      if (reverseRing) std::reverse(ring.begin(), ring.end());
      for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(n); ++iConfig)
      {
        std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
            ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
        ASSERT_TRUE(orient_tets_positively(tets, coords)) << "n " << n << " config " << iConfig;
        EXPECT_NEAR(cavityVolume, total_abs_volume(tets, coords), 1.e-12)
            << "n " << n << " config " << iConfig;
        for (auto & tet : tets)
        {
          EXPECT_LT(
              0., tet_volume({{coords(tet[0]), coords(tet[1]), coords(tet[2]), coords(tet[3])}}));
          EXPECT_FALSE(has_both(tet, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1()));
        }
      }
    }
  }
}

// Ring nodes perturbed in angle, radius and height, around an edge tilted off the axis: the ring is
// still convex as seen along the edge, so every configuration fills the cavity.
TEST(TetSwaps, irregularConvexRing_everyConfigurationFillsTheCavity)
{
  for (int n = 3; n <= 7; ++n)
  {
    std::vector<stk::math::Vector3d> ringCoords;
    for (int i = 0; i < n; ++i)
    {
      const double angle = 2. * M_PI * i / n + 0.1 * std::sin(i + 1.);
      const double radius = 1. + 0.05 * std::cos(2. * i);
      ringCoords.emplace_back(
          radius * std::cos(angle), radius * std::sin(angle), 0.3 * std::sin(3. * i + 1.));
    }
    RingAroundEdge ringAroundEdge(ringCoords);
    ringAroundEdge.set_coords(ringAroundEdge.edge_node0(), stk::math::Vector3d(0.1, -0.05, -1.5));
    ringAroundEdge.set_coords(ringAroundEdge.edge_node1(), stk::math::Vector3d(-0.05, 0.1, 1.5));
    const auto coords = ringAroundEdge.coords_function();
    std::vector<std::array<int, 4>> tetsAroundEdge = ringAroundEdge.tets_around_edge();
    ASSERT_TRUE(orient_tets_positively(tetsAroundEdge, coords)) << "n " << n;
    const double cavityVolume = total_abs_volume(tetsAroundEdge, coords);
    const std::vector<int> ring = ringAroundEdge.ring();
    for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(n); ++iConfig)
    {
      std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
          ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
      ASSERT_TRUE(orient_tets_positively(tets, coords)) << "n " << n << " config " << iConfig;
      EXPECT_NEAR(cavityVolume, total_abs_volume(tets, coords), 1.e-12)
          << "n " << n << " config " << iConfig;
    }
  }
}

// Ring nodes 0, 1 and 2 are collinear, as seen along the edge, so the new edge (0,2) would make
// tets of zero volume.
TEST(TetSwaps, collinearRingNodes_configurationThroughThemInvalid)
{
  const RingAroundEdge ringAroundEdge({{1., -1., 0.}, {1., 0., 0.}, {1., 1., 0.}, {-1., 0., 0.}});
  const auto coords = ringAroundEdge.coords_function();
  const std::vector<int> ring = ringAroundEdge.ring();
  for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(4); ++iConfig)
  {
    std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
        ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
    const bool hasEdgeThroughCollinearNodes =
        has_edge(new_edges_of_tet_edge_swap_configuration(ring, iConfig), 0, 2);
    EXPECT_EQ(!hasEdgeThroughCollinearNodes, orient_tets_positively(tets, coords))
        << "config " << iConfig;
  }
  const int best = best_tet_edge_swap_configuration(ring,
      ringAroundEdge.edge_node0(),
      ringAroundEdge.edge_node1(),
      coords,
      [](const std::vector<std::array<int, 4>> &, size_t) { return true; });
  ASSERT_LE(0, best);
  EXPECT_TRUE(has_edge(new_edges_of_tet_edge_swap_configuration(ring, best), 1, 3));
}

TEST(TetSwaps, tetVolumeSign_zeroForFlatTetsToRoundoff)
{
  for (const double scale : {1., 1.e-6})
  {
    const stk::math::Vector3d offset(3. * scale, -2. * scale, 1. * scale);
    const auto tet = [&](const std::array<stk::math::Vector3d, 4> & x)
    {
      return std::array<stk::math::Vector3d, 4>{{offset + scale * x[0],
          offset + scale * x[1],
          offset + scale * x[2],
          offset + scale * x[3]}};
    };
    const stk::math::Vector3d o(0., 0., 0.), ex(1., 0., 0.), ey(0., 1., 0.), ez(0., 0., 1.);
    EXPECT_EQ(1, tet_volume_sign(tet({{o, ex, ey, ez}}))) << "scale " << scale;
    EXPECT_EQ(-1, tet_volume_sign(tet({{o, ey, ex, ez}}))) << "scale " << scale;
    EXPECT_EQ(0, tet_volume_sign(tet({{o, ex, ey, ex + ey}}))) << "scale " << scale;
    // Nodes 0, 1 and 2 collinear but for 1e-13: zero to roundoff, though not exactly.
    const stk::math::Vector3d nearlyMidway(0.5, 1.e-13, 0.);
    const auto nearlyFlat = tet({{o, nearlyMidway, ex, ez}});
    EXPECT_NE(0., tet_volume(nearlyFlat)) << "scale " << scale;
    EXPECT_EQ(0, tet_volume_sign(nearlyFlat)) << "scale " << scale;
  }
}

// Ring nodes 0, 1 and 2 collinear but for node 1 pushed out by 1e-13, as seen along the edge: the
// tets through the new edge (0,2) have volumes of one sign, but zero to roundoff, so that
// configuration is not valid; the other is.
TEST(TetSwaps, nearlyCollinearRingNodes_configurationThroughThemInvalid)
{
  const RingAroundEdge ringAroundEdge(
      {{1., -1., 0.}, {1. + 1.e-13, 0., 0.}, {1., 1., 0.}, {-1., 0., 0.}});
  const auto coords = ringAroundEdge.coords_function();
  const std::vector<int> ring = ringAroundEdge.ring();
  for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(4); ++iConfig)
  {
    std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
        ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
    const bool hasEdgeThroughNearlyCollinearNodes =
        has_edge(new_edges_of_tet_edge_swap_configuration(ring, iConfig), 0, 2);
    EXPECT_EQ(!hasEdgeThroughNearlyCollinearNodes, orient_tets_positively(tets, coords))
        << "config " << iConfig;
  }
}

// With ring node 0 pulled in toward the edge, the ring is not convex: configurations with the new
// edge between its neighbours do not fill the cavity, but those from node 0 do.
TEST(TetSwaps, nonConvexRing_someConfigurationsInvalid)
{
  RingAroundEdge ringAroundEdge = RingAroundEdge::regular(5);
  ringAroundEdge.set_coords(0, stk::math::Vector3d(0.2, 0., 0.));
  const auto coords = ringAroundEdge.coords_function();
  const std::vector<int> ring = ringAroundEdge.ring();
  size_t numValid = 0;
  for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(5); ++iConfig)
  {
    std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
        ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
    const bool valid = orient_tets_positively(tets, coords);
    const bool hasEdgeBetweenNeighboursOfNode0 =
        has_edge(new_edges_of_tet_edge_swap_configuration(ring, iConfig), 1, 4);
    EXPECT_EQ(!hasEdgeBetweenNeighboursOfNode0, valid) << "config " << iConfig;
    if (valid) ++numValid;
  }
  EXPECT_LT(0u, numValid);
  EXPECT_GT(num_tet_edge_swap_configurations(5), numValid);

  const int best = best_tet_edge_swap_configuration(ring,
      ringAroundEdge.edge_node0(),
      ringAroundEdge.edge_node1(),
      coords,
      [](const std::vector<std::array<int, 4>> &, size_t) { return true; });
  ASSERT_LE(0, best);
  EXPECT_FALSE(has_edge(new_edges_of_tet_edge_swap_configuration(ring, best), 1, 4));
}

TEST(TetSwaps, acceptRequiringNewEdge_configurationWithIt)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(4);
  const std::vector<int> ring = ringAroundEdge.ring();
  for (const std::array<int, 2> & wanted : {std::array<int, 2>{{0, 2}}, std::array<int, 2>{{1, 3}}})
  {
    const int best = best_tet_edge_swap_configuration(ring,
        ringAroundEdge.edge_node0(),
        ringAroundEdge.edge_node1(),
        ringAroundEdge.coords_function(),
        [&](const std::vector<std::array<int, 4>> &, size_t iConfig)
        {
          return has_edge(
              new_edges_of_tet_edge_swap_configuration(ring, iConfig), wanted[0], wanted[1]);
        });
    ASSERT_LE(0, best);
    EXPECT_TRUE(
        has_edge(new_edges_of_tet_edge_swap_configuration(ring, best), wanted[0], wanted[1]));
  }
}

TEST(TetSwaps, acceptRefusingAll_noConfiguration)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(6);
  EXPECT_EQ(-1,
      best_tet_edge_swap_configuration(ringAroundEdge.ring(),
          ringAroundEdge.edge_node0(),
          ringAroundEdge.edge_node1(),
          ringAroundEdge.coords_function(),
          [](const std::vector<std::array<int, 4>> &, size_t) { return false; }));
}

// The configuration chosen is the one whose worst tet is best.
TEST(TetSwaps, stretchedRing_configurationWithBestWorstQuality)
{
  const RingAroundEdge ringAroundEdge = RingAroundEdge::regular(6, 2., 0.7);
  const auto coords = ringAroundEdge.coords_function();
  const std::vector<int> ring = ringAroundEdge.ring();
  int expectedBest = -1;
  double bestQuality = 0.;
  for (size_t iConfig = 0; iConfig < num_tet_edge_swap_configurations(6); ++iConfig)
  {
    std::vector<std::array<int, 4>> tets = tets_of_tet_edge_swap_configuration(
        ring, ringAroundEdge.edge_node0(), ringAroundEdge.edge_node1(), iConfig);
    ASSERT_TRUE(orient_tets_positively(tets, coords));
    const double quality = worst_tet_mean_ratio_quality(tets, coords);
    if (expectedBest < 0 || quality > bestQuality)
    {
      expectedBest = iConfig;
      bestQuality = quality;
    }
  }
  EXPECT_EQ(expectedBest,
      best_tet_edge_swap_configuration(ring,
          ringAroundEdge.edge_node0(),
          ringAroundEdge.edge_node1(),
          coords,
          [](const std::vector<std::array<int, 4>> &, size_t) { return true; }));
}

// The two tets on the face (0,1,2), with nodes 3 and 4 on either side of it.
TEST(TetSwaps, faceSwapWithEdgeThroughFace_threeTetsAroundIt)
{
  const std::vector<stk::math::Vector3d> nodeCoords{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0.25, 0.25, 1.}, {0.3, 0.2, -1.}};
  const auto coords = [&](const int node) { return nodeCoords[node]; };
  const std::array<int, 3> face{{0, 1, 2}};
  const double volumeBefore =
      total_abs_volume(std::vector<std::array<int, 4>>{{{0, 1, 2, 3}}, {{0, 1, 2, 4}}}, coords);

  std::vector<std::array<int, 4>> tets = tets_of_tet_face_swap(face, 3, 4);
  ASSERT_EQ(3u, tets.size());
  ASSERT_TRUE(orient_tets_positively(tets, coords));
  EXPECT_NEAR(volumeBefore, total_abs_volume(tets, coords), 1.e-12);
  for (auto & tet : tets)
  {
    EXPECT_TRUE(has_both(tet, 3, 4));
  }
}

TEST(TetSwaps, faceSwapWithEdgeMissingFace_invalid)
{
  const std::vector<stk::math::Vector3d> nodeCoords{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0.25, 0.25, 1.}, {2., 2., -1.}};
  const auto coords = [&](const int node) { return nodeCoords[node]; };
  std::vector<std::array<int, 4>> tets = tets_of_tet_face_swap(std::array<int, 3>{{0, 1, 2}}, 3, 4);
  EXPECT_FALSE(orient_tets_positively(tets, coords));
}

} // namespace krino
