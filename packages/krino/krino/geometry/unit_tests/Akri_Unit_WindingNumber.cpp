#include <Akri_UnitTestUtils.hpp>
#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <Akri_WindingNumber.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <Akri_WindingNumberCalcs.hpp>
#include <Akri_WindingNumberTopDownTree.hpp>
#include <Akri_WindingNumberBottomUpTree.hpp>

#include <cmath>
#include <iomanip>
#include <limits>

namespace krino {

void expect_winding_number_for_facet(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc, const double goldWindingNumber)
{
  EXPECT_NEAR(goldWindingNumber, compute_exact_facet_winding_number(x0, x1, x2, queryLoc), 1.e-6);
}

void expect_winding_number_in_range_for_point_in_plane_of_facet(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc)
{
  const double windingNumber = compute_exact_facet_winding_number(x0, x1, x2, queryLoc);
  EXPECT_LE(-0.5, windingNumber);
  EXPECT_GE(0.5, windingNumber);
}

TEST(FacetWindingNumber, uniformTetWindingNumberBasedOnAnalyticSolidAngle)
{
  const stk::math::Vector3d x0(1,-1,-1);
  const stk::math::Vector3d x1(-1,-1,1);
  const stk::math::Vector3d x2(-1,1,-1);
  const stk::math::Vector3d queryLoc(1,1,1);

  const double goldWindingNumber = std::acos(23./27.)/(4.*M_PI); // From https://mathworld.wolfram.com/RegularTetrahedron.html
  expect_winding_number_for_facet(x0, x1, x2, queryLoc, goldWindingNumber);
}

TEST(FacetWindingNumber, pointInPlaneOfFacet_getValidWindingNumberInRange)
{
  const stk::math::Vector3d x0(0,0,0);
  const stk::math::Vector3d x1(1,0,0);
  const stk::math::Vector3d x2(0,1,0);
  const stk::math::Vector3d queryLoc(1.,1.,1.e-12);

  expect_winding_number_in_range_for_point_in_plane_of_facet(x0, x1, x2, stk::math::Vector3d(0.2,0.2,0)); // On facet
  expect_winding_number_in_range_for_point_in_plane_of_facet(x0, x1, x2, stk::math::Vector3d(0,0,0)); // On vertex
  expect_winding_number_in_range_for_point_in_plane_of_facet(x0, x1, x2, stk::math::Vector3d(0,0.5,0)); // On edge
  expect_winding_number_in_range_for_point_in_plane_of_facet(x0, x1, x2, stk::math::Vector3d(0.6,0.6,0)); // Outside facet
}

stk::math::Vector3d compute_surface_centroid(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets)
{
  stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
  for (const auto & facetCoords : surfFacets)
    centroid += 1./3. * (facetCoords[0] + facetCoords[1] + facetCoords[2]);
  centroid /= surfFacets.size();
  return centroid;
}

double time_incremental_approximate_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const std::vector<stk::math::Vector3d> & queryLocs, const stk::math::Vector3d & centroid)
{
  const double startTime = stk::wall_time();

  ClusterApproximation approx;
  compute_cluster_approximation(surfFacets, centroid, approx);
  for (const auto & queryLoc: queryLocs)
    compute_approximate_winding_number(approx, queryLoc);

  return stk::wall_time() - startTime;
}

double time_approximate_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const std::vector<stk::math::Vector3d> & queryLocs, const stk::math::Vector3d & centroid)
{
  const double startTime = stk::wall_time();

  FacetClusterApproximation approx;
  compute_cluster_approximation(surfFacets, centroid, approx);
  for (const auto & queryLoc: queryLocs)
    compute_approximate_winding_number(approx, queryLoc);

  return stk::wall_time() - startTime;
}

double time_exact_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const std::vector<stk::math::Vector3d> & queryLocs)
{
  const double startTime = stk::wall_time();

  for (const auto & queryLoc: queryLocs)
    compute_faceted_surface_winding_number(surfFacets, queryLoc);

  return stk::wall_time() - startTime;
}

void test_performance_for_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const std::vector<stk::math::Vector3d> & queryLocs)
{
  const stk::math::Vector3d centroid = compute_surface_centroid(surfFacets);

  std::cout << "Incremental approx time = " << time_incremental_approximate_winding_number(surfFacets, queryLocs, centroid) << std::endl;
  std::cout << "Approx time = " << time_approximate_winding_number(surfFacets, queryLocs, centroid) << std::endl;
  std::cout << "Exact time = " << time_exact_winding_number(surfFacets, queryLocs) << std::endl;
}

double compute_approximate_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc, const stk::math::Vector3d & centroid)
{
  FacetClusterApproximation approx;
  compute_cluster_approximation(surfFacets, centroid, approx);
  return compute_approximate_winding_number(approx, queryLoc);
}

double compute_incremental_approximate_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc, const stk::math::Vector3d & centroid)
{
  ClusterApproximation approx;
  compute_cluster_approximation(surfFacets, centroid, approx);
  return compute_approximate_winding_number(approx, queryLoc);
}

void expect_approximate_winding_number_to_match_exact(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc, const double relativeTol)
{
  const double exactWinding = compute_faceted_surface_winding_number(surfFacets, queryLoc);

  const stk::math::Vector3d centroid = compute_surface_centroid(surfFacets);

  const double approxWinding = compute_approximate_winding_number(surfFacets, queryLoc, centroid);
  EXPECT_NEAR(exactWinding, approxWinding, exactWinding*relativeTol);

  const double incrementalApproxWinding = compute_incremental_approximate_winding_number(surfFacets, queryLoc, centroid);
  EXPECT_NEAR(exactWinding, incrementalApproxWinding, exactWinding*relativeTol);

  EXPECT_NEAR(approxWinding, incrementalApproxWinding, exactWinding*1.e-10);

  std::cout << "For queryLoc " << queryLoc << ", exact = " << exactWinding << ", approx = " << approxWinding << ", incremental approx = " << incrementalApproxWinding << std::endl;
}

void append_refined_facet(const std::array<stk::math::Vector3d,3> & facetCoords, const unsigned numRefine, std::vector<std::array<stk::math::Vector3d,3>> & refinedFacets)
{
  if (numRefine == 0)
  {
    refinedFacets.push_back(facetCoords);
  }
  else
  {
    const stk::math::Vector3d edge0 = 0.5*(facetCoords[0]+facetCoords[1]);
    const stk::math::Vector3d edge1 = 0.5*(facetCoords[1]+facetCoords[2]);
    const stk::math::Vector3d edge2 = 0.5*(facetCoords[2]+facetCoords[0]);
    append_refined_facet({{facetCoords[0], edge0, edge2}}, numRefine-1, refinedFacets);
    append_refined_facet({{facetCoords[1], edge1, edge0}}, numRefine-1, refinedFacets);
    append_refined_facet({{facetCoords[2], edge2, edge1}}, numRefine-1, refinedFacets);
    append_refined_facet({{edge0, edge1, edge2}}, numRefine-1, refinedFacets);
  }
}

std::vector<std::array<stk::math::Vector3d,3>> initialize_nonplanar_facets(const unsigned numRefine)
{
  const std::array<stk::math::Vector3d,3> facet0Coords{{ stk::math::Vector3d(0,1,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0.4,0.5,0.6) }};
  const std::array<stk::math::Vector3d,3> facet1Coords{{ stk::math::Vector3d(0,0,1), stk::math::Vector3d(0,1,0), stk::math::Vector3d(0.4,0.5,0.6) }};
  const std::array<stk::math::Vector3d,3> facet2Coords{{ stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,0,1), stk::math::Vector3d(0.4,0.5,0.6) }};

  std::vector<std::array<stk::math::Vector3d,3>> surfFacets;
  append_refined_facet(facet0Coords, numRefine, surfFacets);
  append_refined_facet(facet1Coords, numRefine, surfFacets);
  append_refined_facet(facet2Coords, numRefine, surfFacets);
  return surfFacets;
}

TEST(approximateWindingNumber, showConvergenceAndAgreementBetweenApproximateMethods)
{
  const std::vector<std::array<stk::math::Vector3d,3>> surfFacets = initialize_nonplanar_facets(3);

  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(1,1,1), 2.e-1);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(2,2,2), 1.e-2);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(4,4,4), 1.e-3);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(8,8,8), 1.e-4);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(1,2,4), 1.e-2);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(4,1,2), 1.e-2);
  expect_approximate_winding_number_to_match_exact(surfFacets, stk::math::Vector3d(2,4,1), 1.e-2);
}

TEST(approximateWindingNumber, compareCPUTimesForApproximateAndExactMethods)
{
  SkipTestMsgIf(is_debug() || stk::EnvData::parallel_size() > 1, "Skipping test debug or in parallel");

  const std::vector<std::array<stk::math::Vector3d,3>> surfFacets = initialize_nonplanar_facets(6);

  const unsigned dim=20;
  std::vector<stk::math::Vector3d> queryLocs;
  for (unsigned i=0; i<dim; ++i)
    for (unsigned j=0; j<dim; ++j)
      for (unsigned k=0; k<dim; ++k)
        queryLocs.emplace_back(-10. + 20.*i/dim, -10. + 20.*j/dim, -10. + 20.*k/dim);

  test_performance_for_winding_number(surfFacets, queryLocs);
}

//
// 2D winding-number tests.
//

// A closed, counter-clockwise polygon (approximating a circle), returned as
// oriented segments. A closed loop has total winding number 1 for interior
// query points and 0 for exterior points.
std::vector<std::array<stk::math::Vector3d,2>> initialize_polygon_segments(const unsigned numSides)
{
  std::vector<stk::math::Vector3d> verts;
  for (unsigned i=0; i<numSides; ++i)
  {
    const double angle = 2.*M_PI*i/numSides;
    verts.emplace_back(std::cos(angle), std::sin(angle), 0.);
  }

  std::vector<std::array<stk::math::Vector3d,2>> segments;
  for (unsigned i=0; i<numSides; ++i)
    segments.push_back({{ verts[i], verts[(i+1)%numSides] }});
  return segments;
}

double compute_exact_polygon_winding_number(const std::vector<std::array<stk::math::Vector3d,2>> & segments, const stk::math::Vector3d & queryLoc)
{
  double windingNumber = 0.;
  for (const auto & seg : segments)
    windingNumber += compute_exact_facet2d_winding_number(seg[0], seg[1], queryLoc);
  return windingNumber;
}

TEST(Facet2dWindingNumber, closedPolygonHasUnitWindingInsideAndZeroOutside)
{
  const auto segments = initialize_polygon_segments(7);

  EXPECT_NEAR(1., compute_exact_polygon_winding_number(segments, stk::math::Vector3d(0,0,0)), 1.e-9);
  EXPECT_NEAR(1., compute_exact_polygon_winding_number(segments, stk::math::Vector3d(0.3,-0.2,0)), 1.e-9);
  EXPECT_NEAR(0., compute_exact_polygon_winding_number(segments, stk::math::Vector3d(5,0,0)), 1.e-9);
  EXPECT_NEAR(0., compute_exact_polygon_winding_number(segments, stk::math::Vector3d(-3,4,0)), 1.e-9);
}

stk::math::Vector3d compute_polygon_centroid(const std::vector<std::array<stk::math::Vector3d,2>> & segments)
{
  stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
  for (const auto & seg : segments)
    centroid += 0.5*(seg[0]+seg[1]);
  centroid /= segments.size();
  return centroid;
}

void expect_approximate_2d_winding_number_to_match_exact(const std::vector<std::array<stk::math::Vector3d,2>> & segments, const stk::math::Vector3d & queryLoc, const double absTol)
{
  const double exactWinding = compute_exact_polygon_winding_number(segments, queryLoc);
  const stk::math::Vector3d centroid = compute_polygon_centroid(segments);

  // Direct segment accumulation about the centroid (top-down build path).
  Cluster2dApproximation directApprox;
  compute_cluster_approximation(segments, centroid, directApprox);
  const double directWinding = compute_approximate_winding_number(directApprox, queryLoc);

  // Same cluster assembled from per-segment approximations merged in (bottom-up
  // build path).
  Cluster2dApproximation incrementalApprox;
  incrementalApprox.clear();
  incrementalApprox.center = centroid;
  for (const auto & seg : segments)
  {
    Cluster2dApproximation segApprox;
    compute_facet_approximation(seg, segApprox);
    accumulate_cluster_approximation(segApprox, incrementalApprox);
  }
  const double incrementalWinding = compute_approximate_winding_number(incrementalApprox, queryLoc);

  EXPECT_NEAR(exactWinding, directWinding, absTol);
  EXPECT_NEAR(exactWinding, incrementalWinding, absTol);
  // The two assembly paths produce the same moments, so they must agree to round-off.
  EXPECT_NEAR(directWinding, incrementalWinding, 1.e-10);
}

TEST(approximate2dWindingNumber, convergesToExactAsQueryPointMovesAway)
{
  const auto segments = initialize_polygon_segments(12);

  // Error of a 2nd-order multipole falls ~ (r_cluster/dist)^3, so tolerance
  // tightens as the query point moves away from the unit-radius polygon.
  expect_approximate_2d_winding_number_to_match_exact(segments, stk::math::Vector3d(3,0,0), 5.e-2);
  expect_approximate_2d_winding_number_to_match_exact(segments, stk::math::Vector3d(6,0,0), 5.e-3);
  expect_approximate_2d_winding_number_to_match_exact(segments, stk::math::Vector3d(12,0,0), 1.e-3);
  expect_approximate_2d_winding_number_to_match_exact(segments, stk::math::Vector3d(5,7,0), 5.e-3);
}

TEST(approximate2dWindingNumber, composableMergeMatchesDirectAccumulationToRoundoff)
{
  const auto segments = initialize_polygon_segments(8);
  const stk::math::Vector3d center(0.1, -0.2, 0.);

  // Direct accumulation of all segments about a single center.
  Cluster2dApproximation direct;
  direct.clear();
  direct.center = center;
  for (const auto & seg : segments)
  {
    Cluster2dApproximation segApprox;
    compute_facet_approximation(seg, segApprox);
    accumulate_cluster_approximation(segApprox, direct);
  }

  // Two-level build: split the segments into two sub-clusters about their own
  // (different) centers, then merge each into the parent about `center`.
  const std::size_t mid = segments.size()/2;
  Cluster2dApproximation left;
  left.clear();
  left.center = stk::math::Vector3d(-0.5, 0.4, 0.);
  for (std::size_t i=0; i<mid; ++i)
  {
    Cluster2dApproximation segApprox;
    compute_facet_approximation(segments[i], segApprox);
    accumulate_cluster_approximation(segApprox, left);
  }

  Cluster2dApproximation right;
  right.clear();
  right.center = stk::math::Vector3d(0.6, -0.3, 0.);
  for (std::size_t i=mid; i<segments.size(); ++i)
  {
    Cluster2dApproximation segApprox;
    compute_facet_approximation(segments[i], segApprox);
    accumulate_cluster_approximation(segApprox, right);
  }

  Cluster2dApproximation incremental;
  incremental.clear();
  incremental.center = center;
  accumulate_cluster_approximation(left, incremental);
  accumulate_cluster_approximation(right, incremental);

  // Incrementally-merged moments must equal the directly-accumulated ones to
  // round-off, and therefore so must the evaluated winding number at any query
  // point.
  for (unsigned i=0; i<2; ++i)
    EXPECT_NEAR(direct.areaN[i], incremental.areaN[i], 1.e-12);
  for (unsigned i=0; i<4; ++i)
    EXPECT_NEAR(direct.areaNiXj[i], incremental.areaNiXj[i], 1.e-12);
  for (unsigned i=0; i<6; ++i)
    EXPECT_NEAR(direct.areaNiDjk[i], incremental.areaNiDjk[i], 1.e-12);

  for (const auto & queryLoc : {stk::math::Vector3d(4,0,0), stk::math::Vector3d(-3,5,0), stk::math::Vector3d(2,-6,0)})
    EXPECT_NEAR(compute_approximate_winding_number(direct, queryLoc),
                compute_approximate_winding_number(incremental, queryLoc), 1.e-12);
}

std::vector<stk::math::Vector3d> build_winding_number_tree_query_points()
{
  // A mix of points: far from the surface (approximation should dominate),
  // moderately close, and close enough to force exact leaf evaluation.
  return {
      stk::math::Vector3d(0.4,0.5,0.6), // Near the shared interior vertex
      stk::math::Vector3d(1,1,1),
      stk::math::Vector3d(2,2,2),
      stk::math::Vector3d(4,4,4),
      stk::math::Vector3d(8,8,8),
      stk::math::Vector3d(-5,3,-2),
      stk::math::Vector3d(1,2,4),
      stk::math::Vector3d(4,1,2),
      stk::math::Vector3d(2,4,1)
  };
}

TEST(WindingNumberTopDownTree, emptyTreeReturnsZero)
{
  WindingNumberFacetArrayTopDownTree::FacetVector noFacets;
  WindingNumberFacetArrayTopDownTree tree(noFacets);

  EXPECT_TRUE(tree.empty());
  EXPECT_EQ(0u, tree.num_nodes());
  EXPECT_EQ(0u, tree.num_facets());
  EXPECT_DOUBLE_EQ(0., tree.compute_winding_number(stk::math::Vector3d(1,1,1)));
  EXPECT_DOUBLE_EQ(0., tree.compute_exact_winding_number(stk::math::Vector3d(1,1,1)));
}

TEST(WindingNumberTopDownTree, treeExactSummationMatchesDirectFacetSummation)
{
  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets = initialize_nonplanar_facets(3);
  WindingNumberFacetArrayTopDownTree tree(surfFacets);

  ASSERT_FALSE(tree.empty());
  EXPECT_EQ(surfFacets.size(), tree.num_facets());

  for (const auto & queryLoc : build_winding_number_tree_query_points())
  {
    const double gold = compute_faceted_surface_winding_number(surfFacets, queryLoc);

    // The tree stores a permutation of the facets, so its exact summation must
    // reproduce the direct summation to round-off regardless of ordering.
    EXPECT_NEAR(gold, tree.compute_exact_winding_number(queryLoc), 1.e-12);

    // theta <= 0 disables approximation everywhere, so the approximate entry
    // point must also reproduce the exact result.
    EXPECT_NEAR(gold, tree.compute_winding_number(queryLoc, 0.), 1.e-12);
  }
}

TEST(WindingNumberTopDownTree, approximationConvergesToExactAsThetaDecreases)
{
  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets = initialize_nonplanar_facets(4);
  WindingNumberFacetArrayTopDownTree tree(surfFacets);

  // Well away from the surface the multipole approximation should be accurate
  // even for a loose theta, and tighten monotonically as theta shrinks.
  const stk::math::Vector3d queryLoc(8,8,8);
  const double gold = tree.compute_exact_winding_number(queryLoc);

  const double errLooseTheta = std::abs(tree.compute_winding_number(queryLoc, 0.8) - gold);
  const double errTightTheta = std::abs(tree.compute_winding_number(queryLoc, 0.1) - gold);
  const double errExact      = std::abs(tree.compute_winding_number(queryLoc, 0.)  - gold);

  EXPECT_LE(errTightTheta, errLooseTheta);
  EXPECT_LE(errExact, errTightTheta);
  EXPECT_NEAR(gold, tree.compute_winding_number(queryLoc, 0.1), std::abs(gold)*1.e-3);
}

TEST(WindingNumberTopDownTree, approximateQueryStaysCloseToExactAcrossQueryPoints)
{
  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets = initialize_nonplanar_facets(4);
  WindingNumberFacetArrayTopDownTree tree(surfFacets);

  const double theta = 0.4;
  for (const auto & queryLoc : build_winding_number_tree_query_points())
  {
    const double gold = tree.compute_exact_winding_number(queryLoc);

    WindingNumberFacetArrayTopDownTree::QueryStats stats;
    const double approx = tree.compute_winding_number(queryLoc, theta, &stats);

    // Absolute tolerance keeps points near the surface (winding number ~ 0)
    // from demanding an unrealistic relative accuracy.
    EXPECT_NEAR(gold, approx, 5.e-3);

    // Every query must at least visit the root, and the traversal accounting
    // should be internally consistent.
    EXPECT_GE(stats.nodesVisited, 1u);
    EXPECT_LE(stats.exactFacetsEvaluated, tree.num_facets());
  }
}

TEST(WindingNumberTopDownTree, singleFacetTreeMatchesFacetWindingNumber)
{
  const WindingNumberFacetArrayTopDownTree::Facet facet{{
      stk::math::Vector3d(0,0,0),
      stk::math::Vector3d(1,0,0),
      stk::math::Vector3d(0,1,0) }};
  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets{facet};
  WindingNumberFacetArrayTopDownTree tree(surfFacets);

  EXPECT_EQ(1u, tree.num_nodes());

  const stk::math::Vector3d queryLoc(0.25,0.25,1.);
  const double gold = compute_exact_facet_winding_number(facet[0], facet[1], facet[2], queryLoc);

  EXPECT_NEAR(gold, tree.compute_exact_winding_number(queryLoc), 1.e-12);
  EXPECT_NEAR(gold, tree.compute_winding_number(queryLoc, 0.5), 1.e-12);
}

TEST(WindingNumberTopDownTree, approximateLeafPolicyStaysAccurateAndUsesLeafApproximations)
{
  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets = initialize_nonplanar_facets(4);

  // Small leaves make it likely that some accepted nodes are leaves, so the
  // policy actually changes traversal behavior.
  const std::size_t leafSize = 4;
  WindingNumberFacetArrayTopDownTree exactLeafTree(surfFacets, leafSize, false);
  WindingNumberFacetArrayTopDownTree approxLeafTree(surfFacets, leafSize, true);

  const double theta = 0.5;
  std::size_t totalApproxLeafAccepted = 0;
  for (const auto & queryLoc : build_winding_number_tree_query_points())
  {
    const double gold = exactLeafTree.compute_exact_winding_number(queryLoc);

    WindingNumberFacetArrayTopDownTree::QueryStats exactStats;
    WindingNumberFacetArrayTopDownTree::QueryStats approxStats;
    const double exactLeafWinding = exactLeafTree.compute_winding_number(queryLoc, theta, &exactStats);
    const double approxLeafWinding = approxLeafTree.compute_winding_number(queryLoc, theta, &approxStats);

    // Both leaf policies must stay accurate; the approximate-leaf policy trades
    // exact facet work for a leaf-level expansion, not accuracy.
    EXPECT_NEAR(gold, exactLeafWinding, 5.e-3);
    EXPECT_NEAR(gold, approxLeafWinding, 5.e-3);

    // The approximate-leaf tree should never do MORE exact facet work than the
    // exact-leaf tree, since any accepted leaf replaces exact facets with one
    // approximation evaluation.
    EXPECT_LE(approxStats.exactFacetsEvaluated, exactStats.exactFacetsEvaluated);

    // The two policies must agree on the exact summation entry point.
    EXPECT_NEAR(exactLeafTree.compute_winding_number(queryLoc, 0.),
                approxLeafTree.compute_winding_number(queryLoc, 0.), 1.e-12);

    if (approxStats.exactFacetsEvaluated < exactStats.exactFacetsEvaluated)
      ++totalApproxLeafAccepted;
  }

  // Sanity check that at least one query actually exercised a leaf-level
  // approximation; otherwise the policy would be untested here.
  EXPECT_GT(totalApproxLeafAccepted, 0u);
}

// Best (least-noisy) wall time over numRepeats sweeps of a query loop. The
// callable returns a checksum so the compiler cannot elide the timed work.
template <typename SweepFn>
double best_sweep_time(const unsigned numRepeats, SweepFn && sweep)
{
  double bestElapsed = std::numeric_limits<double>::max();
  double checksum = 0.;
  for (unsigned r = 0; r < numRepeats; ++r)
  {
    const double startTime = stk::wall_time();
    checksum += sweep();
    bestElapsed = std::min(bestElapsed, stk::wall_time() - startTime);
  }
  EXPECT_TRUE(std::isfinite(checksum));
  return bestElapsed;
}

// Worst absolute error of the tree's approximate query against precomputed gold
// winding numbers (brute-force exact summation).
double worst_abs_error(const WindingNumberFacetArrayTopDownTree & tree, const std::vector<stk::math::Vector3d> & queryLocs, const std::vector<double> & goldWinding, const double theta)
{
  double worst = 0.;
  for (std::size_t i = 0; i < queryLocs.size(); ++i)
    worst = std::max(worst, std::abs(tree.compute_winding_number(queryLocs[i], theta) - goldWinding[i]));
  return worst;
}

TEST(WindingNumberTopDownTree, benchmarkTreePerformance)
{
  // Not gated on is_debug(): every method is timed under identical build
  // settings, so the relative comparisons hold even with asserts enabled
  // (absolute times are inflated in that case). Sized to run in a few seconds.
  SkipTestMsgIf(stk::EnvData::parallel_size() > 1, "Skipping benchmark in parallel");

  const WindingNumberFacetArrayTopDownTree::FacetVector surfFacets = initialize_nonplanar_facets(5); // ~3k facets

  const unsigned dim = 12;
  std::vector<stk::math::Vector3d> queryLocs;
  for (unsigned i=0; i<dim; ++i)
    for (unsigned j=0; j<dim; ++j)
      for (unsigned k=0; k<dim; ++k)
        queryLocs.emplace_back(-10. + 20.*i/dim, -10. + 20.*j/dim, -10. + 20.*k/dim);

  WindingNumberFacetArrayTopDownTree tree(surfFacets);

  // Compute the exact (brute-force) winding number once and reuse it as the
  // accuracy gold for every configuration below. The exact summation is
  // independent of tree layout, so it is valid for all trees here.
  std::vector<double> goldWinding(queryLocs.size());
  for (std::size_t i = 0; i < queryLocs.size(); ++i)
    goldWinding[i] = tree.compute_exact_winding_number(queryLocs[i]);

  // Brute force uses few repeats: each sweep is the expensive one, and its
  // timing is already stable.
  const double bruteForceTime = best_sweep_time(3, [&]()
  {
    double s = 0.;
    for (const auto & queryLoc : queryLocs)
      s += tree.compute_exact_winding_number(queryLoc);
    return s;
  });

  // Tree sweeps are sub-millisecond, so more repeats cost little and reduce noise.
  const unsigned treeRepeats = 50;

  std::cout << "WindingNumberTopDownTree benchmark: " << surfFacets.size()
            << " facets, " << queryLocs.size() << " queries" << std::endl;

  // Part 1: tree (approximate) vs brute force across theta.
  std::cout << "  brute-force (sum all facets) = " << bruteForceTime << " s" << std::endl;
  std::cout << "     theta |    tree (s) | speedup | worst abs err" << std::endl;
  for (const double theta : {0.2, 0.5, 1.0})
  {
    const double treeTime = best_sweep_time(treeRepeats, [&]()
    {
      double s = 0.;
      for (const auto & queryLoc : queryLocs)
        s += tree.compute_winding_number(queryLoc, theta);
      return s;
    });
    const double worstErr = worst_abs_error(tree, queryLocs, goldWinding, theta);

    std::cout << "  " << std::setw(8) << theta
              << " | " << std::setw(11) << treeTime
              << " | " << std::setw(7) << (bruteForceTime / treeTime)
              << " | " << std::setw(13) << worstErr << std::endl;

    // The whole point of the tree: an approximate query must beat summing every
    // facet. Guarded loosely to stay robust across machines and build types.
    EXPECT_LT(treeTime, bruteForceTime);
  }

  // Part 2: exact-leaf vs approximate-leaf policy across leafSize (theta = 0.5).
  const double leafTheta = 0.5;
  std::cout << "  leaf policy (theta = " << leafTheta << "):" << std::endl;
  std::cout << "  leafSize |  exact-leaf (s) | approx-leaf (s) | worst err exact | worst err approx" << std::endl;
  for (const std::size_t leafSize : {std::size_t{4}, std::size_t{8}, std::size_t{16}, std::size_t{32}, std::size_t{64}})
  {
    WindingNumberFacetArrayTopDownTree exactLeafTree(surfFacets, leafSize, false);
    WindingNumberFacetArrayTopDownTree approxLeafTree(surfFacets, leafSize, true);

    const double exactLeafTime = best_sweep_time(treeRepeats, [&]()
    {
      double s = 0.;
      for (const auto & queryLoc : queryLocs)
        s += exactLeafTree.compute_winding_number(queryLoc, leafTheta);
      return s;
    });
    const double approxLeafTime = best_sweep_time(treeRepeats, [&]()
    {
      double s = 0.;
      for (const auto & queryLoc : queryLocs)
        s += approxLeafTree.compute_winding_number(queryLoc, leafTheta);
      return s;
    });

    const double exactErr = worst_abs_error(exactLeafTree, queryLocs, goldWinding, leafTheta);
    const double approxErr = worst_abs_error(approxLeafTree, queryLocs, goldWinding, leafTheta);

    std::cout << "  " << std::setw(8) << leafSize
              << " | " << std::setw(15) << exactLeafTime
              << " | " << std::setw(15) << approxLeafTime
              << " | " << std::setw(15) << exactErr
              << " | " << std::setw(16) << approxErr << std::endl;
  }
}

TEST(WindingNumberBottomUpTree, emptyTreeReturnsZero)
{
  WindingNumberFacetArrayBottomUpTree::FacetVector noFacets;
  WindingNumberFacetArrayBottomUpTree tree(noFacets);

  EXPECT_TRUE(tree.empty());
  EXPECT_EQ(0u, tree.num_nodes());
  EXPECT_EQ(0u, tree.num_facets());
  EXPECT_DOUBLE_EQ(0., tree.compute_winding_number(stk::math::Vector3d(1,1,1)));
  EXPECT_DOUBLE_EQ(0., tree.compute_exact_winding_number(stk::math::Vector3d(1,1,1)));
}

TEST(WindingNumberBottomUpTree, everyNodeCarriesAnApproximation)
{
  const WindingNumberFacetArrayBottomUpTree::FacetVector surfFacets = initialize_nonplanar_facets(4);
  WindingNumberFacetArrayBottomUpTree tree(surfFacets, 4);

  ASSERT_FALSE(tree.empty());
  EXPECT_EQ(surfFacets.size(), tree.num_facets());

  // Unlike the top-down tree there is no exact-leaf policy: the approximation
  // storage is 1:1 with the nodes, so even the smallest leaf can be accepted.
  WindingNumberFacetArrayBottomUpTree::QueryStats stats;
  tree.compute_winding_number(stk::math::Vector3d(8,8,8), 0.5, &stats);
  EXPECT_GT(stats.approximateNodesUsed, 0u);
}

TEST(WindingNumberBottomUpTree, treeExactSummationMatchesDirectFacetSummation)
{
  const WindingNumberFacetArrayBottomUpTree::FacetVector surfFacets = initialize_nonplanar_facets(3);
  WindingNumberFacetArrayBottomUpTree tree(surfFacets);

  ASSERT_FALSE(tree.empty());
  EXPECT_EQ(surfFacets.size(), tree.num_facets());

  for (const auto & queryLoc : build_winding_number_tree_query_points())
  {
    const double gold = compute_faceted_surface_winding_number(surfFacets, queryLoc);

    // The tree stores a permutation of the facets, so its exact summation must
    // reproduce the direct summation to round-off regardless of ordering.
    EXPECT_NEAR(gold, tree.compute_exact_winding_number(queryLoc), 1.e-12);

    // theta <= 0 disables approximation everywhere, so the approximate entry
    // point must also reproduce the exact result.
    EXPECT_NEAR(gold, tree.compute_winding_number(queryLoc, 0.), 1.e-12);
  }
}

TEST(WindingNumberBottomUpTree, approximationConvergesToExactAsThetaDecreases)
{
  const WindingNumberFacetArrayBottomUpTree::FacetVector surfFacets = initialize_nonplanar_facets(4);
  WindingNumberFacetArrayBottomUpTree tree(surfFacets);

  // Well away from the surface the multipole approximation should be accurate
  // even for a loose theta, and tighten monotonically as theta shrinks.
  const stk::math::Vector3d queryLoc(8,8,8);
  const double gold = tree.compute_exact_winding_number(queryLoc);

  const double errLooseTheta = std::abs(tree.compute_winding_number(queryLoc, 0.8) - gold);
  const double errTightTheta = std::abs(tree.compute_winding_number(queryLoc, 0.1) - gold);
  const double errExact      = std::abs(tree.compute_winding_number(queryLoc, 0.)  - gold);

  EXPECT_LE(errTightTheta, errLooseTheta);
  EXPECT_LE(errExact, errTightTheta);
  EXPECT_NEAR(gold, tree.compute_winding_number(queryLoc, 0.1), std::abs(gold)*1.e-3);
}

TEST(WindingNumberBottomUpTree, approximateQueryStaysCloseToExactAcrossQueryPoints)
{
  const WindingNumberFacetArrayBottomUpTree::FacetVector surfFacets = initialize_nonplanar_facets(4);
  WindingNumberFacetArrayBottomUpTree tree(surfFacets);

  const double theta = 0.4;
  for (const auto & queryLoc : build_winding_number_tree_query_points())
  {
    const double gold = tree.compute_exact_winding_number(queryLoc);

    WindingNumberFacetArrayBottomUpTree::QueryStats stats;
    const double approx = tree.compute_winding_number(queryLoc, theta, &stats);

    // Absolute tolerance keeps points near the surface (winding number ~ 0)
    // from demanding an unrealistic relative accuracy.
    EXPECT_NEAR(gold, approx, 5.e-3);

    EXPECT_GE(stats.nodesVisited, 1u);
    EXPECT_LE(stats.exactFacetsEvaluated, tree.num_facets());
  }
}

TEST(WindingNumberBottomUpTree, bottomUpBuildMatchesTopDownTreeToRoundoff)
{
  // The bottom-up (ClusterApproximation, child->parent merge) build must
  // evaluate the SAME expansion as the top-down (FacetClusterApproximation,
  // re-summed) build. Compare both trees node-for-node across query points and
  // acceptance ratios; they should agree to round-off.
  const WindingNumberFacetArrayBottomUpTree::FacetVector surfFacets = initialize_nonplanar_facets(4);

  const std::size_t leafSize = 8;
  WindingNumberFacetArrayBottomUpTree bottomUpTree(surfFacets, leafSize);
  // approximateLeaves = true so both trees can approximate the same leaves.
  WindingNumberFacetArrayTopDownTree topDownTree(surfFacets, leafSize, true);

  ASSERT_EQ(topDownTree.num_nodes(), bottomUpTree.num_nodes());

  for (const auto & queryLoc : build_winding_number_tree_query_points())
    for (const double theta : {0., 0.2, 0.5, 1.0})
      EXPECT_NEAR(topDownTree.compute_winding_number(queryLoc, theta),
                  bottomUpTree.compute_winding_number(queryLoc, theta), 1.e-9);
}

}

