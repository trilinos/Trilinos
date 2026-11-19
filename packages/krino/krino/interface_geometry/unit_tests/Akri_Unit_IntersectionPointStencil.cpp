// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <vector>

#include <gtest/gtest.h>

#include <Akri_SharpFeatureIntersectionPoints.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace krino {

namespace {

constexpr double tol = 1.e-3;

const stk::mesh::Entity node0(1), node1(2), node2(3), node3(4);

void expect_reduced_stencil(std::vector<stk::mesh::Entity> nodes,
    std::vector<double> weights,
    const std::vector<stk::mesh::Entity> & goldNodes,
    const std::vector<double> & goldWeights)
{
  reduce_intersection_point_stencil(nodes, weights, tol);
  EXPECT_EQ(goldNodes, nodes);
  ASSERT_EQ(goldWeights.size(), weights.size());
  for (size_t i=0; i<goldWeights.size(); ++i)
  {
    EXPECT_NEAR(goldWeights[i], weights[i], 1.e-14);
  }
}

}

TEST(IntersectionPointStencil, interiorPointAwayFromFaces_unchanged)
{
  expect_reduced_stencil({node0, node1, node2, node3}, {0.1, 0.2, 0.3, 0.4}, {node0, node1, node2, node3}, {0.1, 0.2, 0.3, 0.4});
}

TEST(IntersectionPointStencil, interiorPointNearFace_facePoint)
{
  expect_reduced_stencil({node0, node1, node2, node3}, {0.2, 0.0005, 0.3, 0.4995}, {node0, node2, node3}, {0.2/0.9995, 0.3/0.9995, 0.4995/0.9995});
}

TEST(IntersectionPointStencil, interiorPointNearEdge_edgePoint)
{
  expect_reduced_stencil({node0, node1, node2, node3}, {0.0004, 0.6, 0.0006, 0.399}, {node1, node3}, {0.6/0.999, 0.399/0.999});
}

TEST(IntersectionPointStencil, interiorPointNearNode_capturedByNode)
{
  expect_reduced_stencil({node0, node1, node2, node3}, {0.0002, 0.0003, 0.9991, 0.0004}, {node2}, {1.});
}

TEST(IntersectionPointStencil, facePointNearEdge_edgePoint)
{
  expect_reduced_stencil({node0, node1, node2}, {0.7, 0.2995, 0.0005}, {node0, node1}, {0.7/0.9995, 0.2995/0.9995});
}

TEST(IntersectionPointStencil, facePointNearNode_capturedByNode)
{
  expect_reduced_stencil({node0, node1, node2}, {0.0005, 0.999, 0.0005}, {node1}, {1.});
}

TEST(IntersectionPointStencil, edgePointNearNode_capturedByNode)
{
  expect_reduced_stencil({node0, node1}, {0.9995, 0.0005}, {node0}, {1.});
}

TEST(IntersectionPointStencil, weightAtTolerance_kept)
{
  expect_reduced_stencil({node0, node1}, {1.-tol, tol}, {node0, node1}, {1.-tol, tol});
}

TEST(IntersectionPointStencil, allWeightsBelowTolerance_nodeWithLargestWeightKept)
{
  expect_reduced_stencil({node0, node1, node2}, {0.0001, 0.0004, 0.0002}, {node1}, {1.});
}

TEST(IntersectionPointStencil, sameNodesInDifferentOrder_samePoint)
{
  const IntersectionPoint a(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 0});
  const IntersectionPoint b(false, {node2, node0, node1}, {0.5, 0.2, 0.3}, {0, 0});
  EXPECT_TRUE(is_same_intersection_point(a, b, tol));
}

TEST(IntersectionPointStencil, weightsWithinTolerance_samePoint)
{
  const IntersectionPoint a(true, {node0, node1}, {0.4, 0.6}, {0});
  const IntersectionPoint b(true, {node0, node1}, {0.4005, 0.5995}, {0});
  EXPECT_TRUE(is_same_intersection_point(a, b, tol));
}

TEST(IntersectionPointStencil, weightsDifferByMoreThanTolerance_differentPoints)
{
  const IntersectionPoint a(true, {node0, node1}, {0.4, 0.6}, {0});
  const IntersectionPoint b(true, {node0, node1}, {0.402, 0.598}, {0});
  EXPECT_FALSE(is_same_intersection_point(a, b, tol));
}

TEST(IntersectionPointStencil, differentDomainsOrNodes_differentPoints)
{
  const IntersectionPoint a(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 0});
  EXPECT_FALSE(is_same_intersection_point(a, IntersectionPoint(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 0, 0}), tol));
  EXPECT_FALSE(is_same_intersection_point(a, IntersectionPoint(true, {node0, node1, node3}, {0.2, 0.3, 0.5}, {0, 0}), tol));
  EXPECT_FALSE(is_same_intersection_point(a, IntersectionPoint(true, {node0, node1}, {0.4, 0.6}, {0, 0}), tol));
}

TEST(IntersectionPointStencil, repeatedPoints_firstKept)
{
  std::vector<IntersectionPoint> points{IntersectionPoint(true, {node0, node1}, {0.4, 0.6}, {0}),
      IntersectionPoint(false, {node1, node2}, {0.5, 0.5}, {0}),
      IntersectionPoint(false, {node1, node0}, {0.5995, 0.4005}, {0}),
      IntersectionPoint(false, {node0, node1}, {0.45, 0.55}, {0})};
  remove_repeated_intersection_points(points, tol);
  ASSERT_EQ(3u, points.size());
  EXPECT_TRUE(points[0].is_owned());
  EXPECT_EQ((std::vector<stk::mesh::Entity>{node1, node2}), points[1].get_nodes());
  EXPECT_EQ((std::vector<double>{0.45, 0.55}), points[2].get_weights());
}

TEST(IntersectionPointStencil, moreOrEqualDomains_capture)
{
  const IntersectionPoint corner(true, {node0, node1, node2}, {0.3, 0.3, 0.4}, {0, 0, 0});
  const IntersectionPoint edgeCrossing(true, {node2, node0, node1}, {0.4005, 0.2995, 0.3}, {0, 0});
  EXPECT_TRUE(is_captured_by_intersection_point(corner, edgeCrossing, tol));
  EXPECT_TRUE(is_captured_by_intersection_point(edgeCrossing, edgeCrossing, tol));
  EXPECT_FALSE(is_captured_by_intersection_point(edgeCrossing, corner, tol));
}

TEST(IntersectionPointStencil, differentNodesOrDistantWeights_notCaptured)
{
  const IntersectionPoint corner(true, {node0, node1, node2}, {0.3, 0.3, 0.4}, {0, 0, 0});
  EXPECT_FALSE(is_captured_by_intersection_point(
      corner, IntersectionPoint(true, {node0, node1, node3}, {0.3, 0.3, 0.4}, {0, 0}), tol));
  EXPECT_FALSE(is_captured_by_intersection_point(
      corner, IntersectionPoint(true, {node0, node1, node2}, {0.3, 0.302, 0.398}, {0, 0}), tol));
  EXPECT_FALSE(is_captured_by_intersection_point(
      corner, IntersectionPoint(true, {node0, node1, node2}, {0.3, 0.3, 0.4}, {1, 1}), tol));
}

TEST(IntersectionPointStencil, cornerThenNearbyCrossingsOfItsEdges_onlyCornerKept)
{
  std::vector<IntersectionPoint> points{
      IntersectionPoint(true, {node0, node1, node2}, {0.3, 0.3, 0.4}, {0, 0, 0}),
      IntersectionPoint(true, {node0, node1, node2}, {0.3005, 0.2995, 0.4}, {0, 0}),
      IntersectionPoint(true, {node0, node1, node2}, {0.2995, 0.3, 0.4005}, {0, 0})};
  remove_repeated_intersection_points(points, tol);
  ASSERT_EQ(1u, points.size());
  EXPECT_EQ(3u, points[0].get_sorted_domains().size());
}

TEST(IntersectionPointStencil, featureCodimension_singleSurfaceAndLevelSetsPerPhase)
{
  EXPECT_EQ(1, sharp_feature_codimension({0}));
  EXPECT_EQ(2, sharp_feature_codimension({1, 1}));
  EXPECT_EQ(3, sharp_feature_codimension({2, 2, 2}));
  EXPECT_EQ(1, sharp_feature_codimension({0, 1}));
  EXPECT_EQ(2, sharp_feature_codimension({0, 1, 2}));
  EXPECT_EQ(3, sharp_feature_codimension({0, 1, 2, 3}));

  EXPECT_TRUE(is_corner_intersection_point(
      IntersectionPoint(true, {node0, node1, node2, node3}, {0.25, 0.25, 0.25, 0.25}, {0, 0, 0})));
  EXPECT_TRUE(is_corner_intersection_point(IntersectionPoint(
      true, {node0, node1, node2, node3}, {0.25, 0.25, 0.25, 0.25}, {0, 1, 2, 3})));
  EXPECT_TRUE(is_sharp_edge_intersection_point(
      IntersectionPoint(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 0})));
  EXPECT_TRUE(is_sharp_edge_intersection_point(
      IntersectionPoint(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 1, 2})));
  EXPECT_FALSE(is_corner_intersection_point(
      IntersectionPoint(true, {node0, node1, node2}, {0.2, 0.3, 0.5}, {0, 1, 2})));
}

static void expect_reduced_simplex_weights(const stk::math::Vector3d & paramCoords,
    const unsigned dim,
    const std::vector<int> & goldOrdinals,
    const std::vector<double> & goldWeights)
{
  std::vector<int> ordinals;
  std::vector<double> weights;
  fill_reduced_simplex_vertex_weights(paramCoords, dim, tol, ordinals, weights);
  EXPECT_EQ(goldOrdinals, ordinals);
  ASSERT_EQ(goldWeights.size(), weights.size());
  for (size_t i = 0; i < goldWeights.size(); ++i)
  {
    EXPECT_NEAR(goldWeights[i], weights[i], 1.e-14);
  }
}

TEST(IntersectionPointStencil, reducedTetVertexWeights_entityThePointIsNear)
{
  expect_reduced_simplex_weights({0.2, 0.3, 0.4}, 3, {0, 1, 2, 3}, {0.1, 0.2, 0.3, 0.4});
  expect_reduced_simplex_weights(
      {0.2, 0.3, 0.4995}, 3, {1, 2, 3}, {0.2 / 0.9995, 0.3 / 0.9995, 0.4995 / 0.9995});
  expect_reduced_simplex_weights({0.0004, 0.6, 0.0006}, 3, {0, 2}, {0.399 / 0.999, 0.6 / 0.999});
  expect_reduced_simplex_weights({0.9995, 0.0002, 0.0002}, 3, {1}, {1.});
}

TEST(IntersectionPointStencil, reducedTriVertexWeights_entityThePointIsNear)
{
  expect_reduced_simplex_weights({0.3, 0.5, 0.}, 2, {0, 1, 2}, {0.2, 0.3, 0.5});
  expect_reduced_simplex_weights({0.4995, 0.5, 0.}, 2, {1, 2}, {0.4995 / 0.9995, 0.5 / 0.9995});
}

TEST(IntersectionPointStencil, nodeDomainsCapturePointDomains_countingRepeats)
{
  EXPECT_TRUE(node_domains_capture_point_domains({0, 0, 0}, {0, 0}));
  EXPECT_TRUE(node_domains_capture_point_domains({0, 0}, {0, 0}));
  EXPECT_FALSE(node_domains_capture_point_domains({0}, {0, 0}));
  EXPECT_FALSE(node_domains_capture_point_domains({0, 0}, {0, 0, 0}));
  EXPECT_TRUE(node_domains_capture_point_domains({0, 1, 2, 3}, {1, 2, 3}));
  EXPECT_FALSE(node_domains_capture_point_domains({0, 1, 2}, {1, 2, 3}));
  EXPECT_FALSE(node_domains_capture_point_domains({}, {0}));
}

} // namespace krino
