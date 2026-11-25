// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <array>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_SharpFeatureLeafUtils.hpp>

namespace krino
{

TEST(SharpFeatureLeafUtils, simplexWeights_firstIsOneMinusTheOthers)
{
  EXPECT_EQ((std::vector<double>{0.25, 0.5, 0.125, 0.125}),
      simplex_weights_at_parametric_coordinates(stk::math::Vector3d(0.5, 0.125, 0.125), 4));
  EXPECT_EQ((std::vector<double>{0.25, 0.5, 0.25}),
      simplex_weights_at_parametric_coordinates(stk::math::Vector3d(0.5, 0.25, 0.), 3));
  EXPECT_EQ((std::vector<double>{0.75, 0.25}),
      simplex_weights_at_parametric_coordinates(stk::math::Vector3d(0.25, 0., 0.), 2));
}

TEST(SharpFeatureLeafUtils, entityOfWeights_dropsWeightsWithinToleranceAndScalesTheRest)
{
  std::vector<int> ordinals;
  std::vector<double> weights;
  fill_entity_of_weights_above_tolerance({0.5, 1.e-12, 0.5, -1.e-12}, 1.e-10, ordinals, weights);
  EXPECT_EQ((std::vector<int>{0, 2}), ordinals);
  EXPECT_EQ((std::vector<double>{0.5, 0.5}), weights);

  fill_entity_of_weights_above_tolerance({1.e-11, 1., 1.e-11}, 1.e-10, ordinals, weights);
  EXPECT_EQ((std::vector<int>{1}), ordinals);
  EXPECT_EQ((std::vector<double>{1.}), weights);
}

TEST(SharpFeatureLeafUtils, locateCrossingOnTetFace_insideNearEdgeOrNearNode)
{
  constexpr double edgeTol = 0.01;
  EXPECT_EQ(TetFaceCrossingLocation::INSIDE, locate_crossing_on_tet_face(stk::math::Vector3d(0.3, 0.3, 0.), edgeTol).location);

  // Weight of node 0 within the tolerance: near the edge of nodes 1 and 2.
  const LocationOnTetFace nearEdge12 = locate_crossing_on_tet_face(stk::math::Vector3d(0.5, 0.495, 0.), edgeTol);
  EXPECT_EQ(TetFaceCrossingLocation::NEAR_EDGE, nearEdge12.location);
  EXPECT_EQ((std::array<int, 2>{{1, 2}}), nearEdge12.edgeNodeOrdinals);

  // Weight of node 2 within the tolerance, by roundoff on the other side: near the edge of nodes 0 and 1.
  const LocationOnTetFace nearEdge01 = locate_crossing_on_tet_face(stk::math::Vector3d(0.4, -1.e-15, 0.), edgeTol);
  EXPECT_EQ(TetFaceCrossingLocation::NEAR_EDGE, nearEdge01.location);
  EXPECT_EQ((std::array<int, 2>{{0, 1}}), nearEdge01.edgeNodeOrdinals);

  EXPECT_EQ(TetFaceCrossingLocation::NEAR_NODE, locate_crossing_on_tet_face(stk::math::Vector3d(0.995, 0.002, 0.), edgeTol).location);
}

namespace
{
const std::vector<int> noDomains;
const std::vector<int> leafUtilsSharpEdgeDomains{0, 0};
const std::vector<int> leafUtilsCornerDomains{0, 0, 0};
} // namespace

TEST(SharpFeatureLeafUtils, tetFacePointIsCaptured_byNodeOfHigherRankAnywhere)
{
  constexpr double edgeTol = 0.01;
  EXPECT_TRUE(tet_face_point_is_captured({{&leafUtilsCornerDomains, &noDomains, &noDomains}},
      stk::math::Vector3d(0.3, 0.3, 0.),
      leafUtilsSharpEdgeDomains,
      edgeTol));
  EXPECT_FALSE(tet_face_point_is_captured({{&leafUtilsSharpEdgeDomains, &noDomains, &noDomains}},
      stk::math::Vector3d(0.3, 0.3, 0.),
      leafUtilsSharpEdgeDomains,
      edgeTol));
}

TEST(SharpFeatureLeafUtils, tetFacePointIsCaptured_byNodeWithItsDomainsWithinTolerance)
{
  constexpr double edgeTol = 0.01;
  const std::array<const std::vector<int> *, 3> node1HasDomains{
      {&noDomains, &leafUtilsSharpEdgeDomains, &noDomains}};
  EXPECT_TRUE(tet_face_point_is_captured(
      node1HasDomains, stk::math::Vector3d(0.995, 0.002, 0.), leafUtilsSharpEdgeDomains, edgeTol));
  EXPECT_FALSE(tet_face_point_is_captured(
      node1HasDomains, stk::math::Vector3d(0.98, 0.01, 0.), leafUtilsSharpEdgeDomains, edgeTol));
  // Within the tolerance of node 0, which does not have the domains.
  EXPECT_FALSE(tet_face_point_is_captured(
      node1HasDomains, stk::math::Vector3d(0.002, 0.002, 0.), leafUtilsSharpEdgeDomains, edgeTol));
}

TEST(SharpFeatureLeafUtils, tetFacePointIsCaptured_byChordWhoseNodesHaveItsDomains)
{
  constexpr double edgeTol = 0.01;
  // Nodes 1 and 2 have the domains: the edge between them, where the weight of node 0 is small.
  const std::array<const std::vector<int> *, 3> chord12{
      {&noDomains, &leafUtilsSharpEdgeDomains, &leafUtilsSharpEdgeDomains}};
  EXPECT_TRUE(tet_face_point_is_captured(
      chord12, stk::math::Vector3d(0.5, 0.495, 0.), leafUtilsSharpEdgeDomains, edgeTol));
  EXPECT_FALSE(tet_face_point_is_captured(
      chord12, stk::math::Vector3d(0.5, 0.4, 0.), leafUtilsSharpEdgeDomains, edgeTol));
  // Near the edge of nodes 0 and 1, which is not a chord.
  EXPECT_FALSE(tet_face_point_is_captured(
      chord12, stk::math::Vector3d(0.5, 0.005, 0.), leafUtilsSharpEdgeDomains, edgeTol));
}

TEST(SharpFeatureLeafUtils, uniqueFacesOfTets_sharedFaceOnceWithBothTetsAndOppositeNodes)
{
  // Two tets sharing the face of nodes 0, 1 and 2.
  const std::vector<std::array<int, 4>> tets{{{0, 1, 2, 3}}, {{0, 2, 1, 4}}};
  const auto faces = unique_faces_of_tets(tets);
  ASSERT_EQ(7u, faces.size());

  size_t numShared = 0;
  for (const auto & face : faces)
  {
    std::array<int, 3> sortedNodes = face.nodes;
    std::sort(sortedNodes.begin(), sortedNodes.end());
    if (sortedNodes == std::array<int, 3>{{0, 1, 2}})
    {
      ++numShared;
      EXPECT_EQ((std::array<int, 3>{{0, 2, 1}}), face.nodes); // as in the first tet's side 3
      EXPECT_EQ((std::vector<std::pair<size_t, int>>{{0, 3}, {1, 4}}), face.tetsAndOppositeNodes);
    }
    else
    {
      EXPECT_EQ(1u, face.tetsAndOppositeNodes.size());
      const auto [iTet, oppositeNode] = face.tetsAndOppositeNodes[0];
      EXPECT_EQ(0, std::count(face.nodes.begin(), face.nodes.end(), oppositeNode));
      EXPECT_EQ(1, std::count(tets[iTet].begin(), tets[iTet].end(), oppositeNode));
    }
  }
  EXPECT_EQ(1u, numShared);
  // The first tet's faces come first, in the order of its sides.
  EXPECT_EQ((std::array<int, 3>{{0, 1, 3}}), faces[0].nodes);
  EXPECT_EQ((std::array<int, 3>{{1, 2, 3}}), faces[1].nodes);
}

} // namespace krino
