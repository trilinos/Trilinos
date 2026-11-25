// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <vector>

#include <gtest/gtest.h>

#include <Akri_CutQualityEstimator.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

TEST(ClipIntersectionPointWeights, weightsBelowTolerance_raisedToItWithTheOthersScaledToKeepTheSum)
{
  std::vector<double> weights{0.01, 0.99};
  clip_intersection_point_weights(weights, 0.1);
  EXPECT_NEAR(0.1, weights[0], 1.e-12);
  EXPECT_NEAR(0.9, weights[1], 1.e-12);

  weights = {0.02, 0.08, 0.9};
  clip_intersection_point_weights(weights, 0.1);
  EXPECT_NEAR(0.1, weights[0], 1.e-12);
  EXPECT_NEAR(0.1, weights[1], 1.e-12);
  EXPECT_NEAR(0.8, weights[2], 1.e-12);
}

TEST(ClipIntersectionPointWeights, noWeightBelowTolerance_unchanged)
{
  std::vector<double> weights{0.3, 0.7};
  clip_intersection_point_weights(weights, 0.1);
  EXPECT_EQ(0.3, weights[0]);
  EXPECT_EQ(0.7, weights[1]);
}

// Triangles with nodes 0 (0,0), 1 (1,0), 2 (0,1) and, for TwoRightTris, 3 (-1,0) in the second
// triangle (0,2,3).  The original node, the one that would be snapped, is node 0.
template <typename MESHSPEC> class CutQualityEstimatorTriFixture : public StkMeshTriFixture
{
public:
  CutQualityEstimatorTriFixture()
  {
    set_valid_proc_sizes_for_test({1});
    build_mesh(meshSpec.nodeLocs,
        meshSpec.allElementConn,
        std::vector<unsigned>(meshSpec.allElementConn.size(), 1u));
  }

protected:
  stk::mesh::Entity node(const unsigned index) const { return get_assigned_node_for_index(index); }

  CutQualityEstimator estimator(const double minIntPtWeight = 1.e-6) const
  {
    return CutQualityEstimator(mMesh,
        mMesh.mesh_meta_data().universal_part(),
        mMesh.mesh_meta_data().coordinate_field(),
        qualityMetric,
        minIntPtWeight,
        true);
  }

  stk::math::Vector3d location(const unsigned index) const
  {
    const stk::math::Vec<double, 2> & x = meshSpec.nodeLocs[index];
    return stk::math::Vector3d(x[0], x[1], 0.);
  }

  double quality(const std::vector<stk::math::Vector3d> & triCoords) const
  {
    return qualityMetric.get_element_quality_metric(2, triCoords);
  }

  IntersectionPoint edge_point(const unsigned index0,
      const unsigned index1,
      const double weight0,
      const std::vector<int> & domains = {0}) const
  {
    return IntersectionPoint(true, {node(index0), node(index1)}, {weight0, 1. - weight0}, domains);
  }

  // The first triangle (0,1,2), with only node 0, or also node 1, on the interface or original.
  double estimate_for_first_tri(const CutQualityEstimator & cutQualityEstimator,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<size_t> & intersectionPointIndices,
      const std::vector<bool> & isNodeOnInterfaceOrIsOriginal = {true, false, false}) const
  {
    return cutQualityEstimator.estimate_quality_of_cutting_element({node(0), node(1), node(2)},
        {location(0), location(1), location(2)},
        isNodeOnInterfaceOrIsOriginal,
        intersectionPoints,
        intersectionPointIndices);
  }

  MESHSPEC meshSpec;
  const ScaledJacobianQualityMetric qualityMetric;
};

using CutQualityEstimatorRightTri = CutQualityEstimatorTriFixture<RightTri>;
using CutQualityEstimatorTwoRightTris = CutQualityEstimatorTriFixture<TwoRightTris>;

// With no points, only an element that snapping would collapse (two or more of its nodes on the
// interface or original) counts, at its own quality.
TEST_F(CutQualityEstimatorRightTri, noPoints_bestUnlessTwoNodesOnInterfaceOrOriginal)
{
  if (!is_valid_proc_size_for_test()) return;

  EXPECT_EQ(qualityMetric.get_best_value_for_metric(), estimate_for_first_tri(estimator(), {}, {}));
  EXPECT_DOUBLE_EQ(quality({location(0), location(1), location(2)}),
      estimate_for_first_tri(estimator(), {}, {}, {true, true, false}));
}

// One point on edge 0-1 splits the triangle at it: only the piece with the original node and the
// point counts, the sliver that snapping node 0 to the point would collapse.
TEST_F(CutQualityEstimatorRightTri, onePointOnEdge_qualityOfThePieceWithTheOriginalNode)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{edge_point(0, 1, 0.7)};
  const stk::math::Vector3d pointLocation(0.3, 0., 0.);
  const double expected = quality({location(0), pointLocation, location(2)});
  EXPECT_LT(expected, qualityMetric.get_best_value_for_metric());

  EXPECT_DOUBLE_EQ(expected, estimate_for_first_tri(estimator(), points, {0}));
  EXPECT_DOUBLE_EQ(
      expected, estimator().estimate_quality_of_cutting_node_elements(node(0), points, {0}));
}

// A point closer to a node than the minimum weight is estimated at that weight from the node.
TEST_F(CutQualityEstimatorRightTri, pointNearNode_estimatedAtTheMinimumWeight)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{edge_point(0, 1, 0.99)};
  const double minIntPtWeight = 0.1;
  const stk::math::Vector3d clippedLocation(0.1, 0., 0.);
  EXPECT_DOUBLE_EQ(quality({location(0), clippedLocation, location(2)}),
      estimate_for_first_tri(estimator(minIntPtWeight), points, {0}));
}

// A point whose nodes are not all in the element does not cut it.
TEST_F(CutQualityEstimatorTwoRightTris, pointOnEdgeOfOtherElement_doesNotCutThisOne)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{edge_point(0, 3, 0.7)};
  EXPECT_EQ(
      qualityMetric.get_best_value_for_metric(), estimate_for_first_tri(estimator(), points, {0}));
}

// The estimate for a node is the worst over the elements its points intersect: here only the
// second triangle (0,2,3), which edge 0-3 is in.
TEST_F(CutQualityEstimatorTwoRightTris, nodeElements_worstOverTheIntersectedElements)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{edge_point(0, 3, 0.7)};
  const stk::math::Vector3d pointLocation(-0.3, 0., 0.);
  EXPECT_DOUBLE_EQ(quality({location(0), location(2), pointLocation}),
      estimator().estimate_quality_of_cutting_node_elements(node(0), points, {0}));
}

// The points of a node are cut with more domains first, then the smaller weight at the node, and
// only those whose domains are all in the given ones.
TEST_F(CutQualityEstimatorRightTri, sortedPointsForNode_moreDomainsFirstThenSmallerWeight)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{
      edge_point(0, 1, 0.6, {0}), edge_point(0, 2, 0.2, {0}), edge_point(0, 1, 0.9, {0, 1})};
  const std::vector<std::pair<size_t, bool>> nodePoints{{0, true}, {1, true}, {2, true}};

  EXPECT_EQ((std::vector<size_t>{2, 1, 0}),
      estimator().sorted_intersection_point_indices_for_node_for_domains(
          points, nodePoints, node(0), {0, 1}));
  EXPECT_EQ((std::vector<size_t>{1, 0}),
      estimator().sorted_intersection_point_indices_for_node_for_domains(
          points, nodePoints, node(0), {0}));
}

// Each set of domains of a node's points gets its own estimate, from the points with those domains.
TEST_F(CutQualityEstimatorRightTri, qualityPerNodePerDomain_oneEstimatePerSetOfDomains)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<IntersectionPoint> points{
      edge_point(0, 1, 0.7, {0}), edge_point(0, 2, 0.5, {1})};
  mapFromEntityToIntPtIndexAndSnapAllowed nodeToPoints;
  nodeToPoints[node(0)] = {{0, true}, {1, true}};

  const auto domainsToNodesToQuality =
      estimator().estimate_quality_per_node_per_domain(points, nodeToPoints);

  const stk::mesh::EntityId nodeId = mMesh.identifier(node(0));
  ASSERT_EQ(2u, domainsToNodesToQuality.size());
  EXPECT_DOUBLE_EQ(quality({location(0), stk::math::Vector3d(0.3, 0., 0.), location(2)}),
      domainsToNodesToQuality.at({0}).at(nodeId));
  EXPECT_DOUBLE_EQ(quality({location(0), location(1), stk::math::Vector3d(0., 0.5, 0.)}),
      domainsToNodesToQuality.at({1}).at(nodeId));
}

} // namespace krino
