#include <gtest/gtest.h>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_SearchTree.hpp>
#include <Akri_Triangle.hpp>

namespace krino {

class TestFacetedSurface
{
public:
  TestFacetedSurface(const std::vector<stk::math::Vector3d> & nodeLocs, const std::vector<std::array<size_t,3>> & facetNodes)
  : myNodeLocs(nodeLocs), myFacetNodes(facetNodes)
  {
    myFacets.resize(myFacetNodes.size());
    std::iota(myFacets.begin(), myFacets.end(), 0);
  }
  stk::math::Vector3d get_facet_centroid(const size_t facet) const
  {
    const std::array<size_t,3> & facetNodes = myFacetNodes[facet];
    return 1./3.*(myNodeLocs[facetNodes[0]] + myNodeLocs[facetNodes[1]] + myNodeLocs[facetNodes[2]]);
  }
  void insert_facet_into_bounding_box(const size_t facet, BoundingBox & bbox) const
  {
    const std::array<size_t,3> & facetNodes = myFacetNodes[facet];
    bbox.accommodate(myNodeLocs[facetNodes[0]]);
    bbox.accommodate(myNodeLocs[facetNodes[1]]);
    bbox.accommodate(myNodeLocs[facetNodes[2]]);
  }

  std::array<stk::math::Vector3d,3> get_facet_node_locations(const size_t facet) const
  {
    const std::array<size_t,3> & facetNodes = myFacetNodes[facet];
    return {myNodeLocs[facetNodes[0]], myNodeLocs[facetNodes[1]], myNodeLocs[facetNodes[2]]};
  }
  const std::vector<stk::math::Vector3d> & get_node_locations() const { return myNodeLocs; }
  const std::vector<std::array<size_t,3>> & get_facet_nodes() const { return myFacetNodes; }
  const std::vector<size_t> & get_facets() const { return myFacets; }

private:
  std::vector<stk::math::Vector3d> myNodeLocs;
  std::vector<std::array<size_t,3>> myFacetNodes;
  std::vector<size_t> myFacets;
};

// These factories return the lambdas directly (not std::function) so that the
// SearchTree constructor and queries instantiate on the concrete callable type
// and inline it -- this is what the templated-callback interface is designed to
// exploit, and what the performance gauge below measures.
static auto build_insert_facet_into_bounding_box_function(const TestFacetedSurface & surf)
{
    return [&surf](const size_t facet, krino::BoundingBox & bbox) { return surf.insert_facet_into_bounding_box(facet, bbox); };
}

static auto build_get_facet_centroid_function(const TestFacetedSurface & surf)
{
    return [&surf](const size_t facet) { return surf.get_facet_centroid(facet); };
}

static auto build_get_facet_distance_squared(const TestFacetedSurface & surf)
{
    return [&surf](const size_t facet, const stk::math::Vector3d & queryLoc)
    {
      const auto & facetNodeLocs = surf.get_facet_node_locations(facet);
      return CalcTriangle3<double>::distance_squared(facetNodeLocs, queryLoc);
    };
}

TestFacetedSurface build_single_square_surface_with_num_facets_per_side_geometry(const size_t numFacetsPerSide)
{
    std::vector<stk::math::Vector3d> nodeLocs;
    std::vector<std::array<size_t,3>> facetNodes;

    const stk::math::Vector3d minEndPt{-5, 0, -5};
    const stk::math::Vector3d maxEndPt{5, 0, 5};

    for (size_t i=0; i<numFacetsPerSide+1; ++i)
    {
      const double iloc = 1.0*i/numFacetsPerSide;
      const double xLoc = (1.-iloc)*minEndPt[0] + iloc*maxEndPt[0];
      for (size_t j=0; j<numFacetsPerSide+1; ++j)
      {
        const double jloc = 1.0*j/numFacetsPerSide;
        const double zLoc = (1.-jloc)*minEndPt[2] + jloc*maxEndPt[2];
        nodeLocs.emplace_back(xLoc, 0, zLoc);
      }
    }

    for (size_t i=0; i<numFacetsPerSide; ++i)
    {
      for (size_t j=0; j<numFacetsPerSide; ++j)
      {
        const size_t xmzm = i*(numFacetsPerSide+1) + j;
        const size_t xpzm = (i+1)*(numFacetsPerSide+1) + j;
        const size_t xpzp = (i+1)*(numFacetsPerSide+1) + j+1;
        const size_t xmzp = i*(numFacetsPerSide+1) + j+1;
        facetNodes.push_back({xmzm, xpzm, xpzp});
        facetNodes.push_back({xmzm, xpzp, xmzp});
      }
    }

    return TestFacetedSurface(nodeLocs, facetNodes);
}

void test_num_closest_facets_with_num_facets_per_side(const unsigned numFacetsPerSide, const stk::math::Vector3d & queryLoc, const unsigned goldMinNumClosestFacets, const unsigned goldMaxNumClosestFacets)
{
    const TestFacetedSurface surf = build_single_square_surface_with_num_facets_per_side_geometry(numFacetsPerSide);
    SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

    std::vector<size_t> closestFacets;
    facetTree.find_closest_entities(queryLoc, build_get_facet_distance_squared(surf), closestFacets);

    EXPECT_LE(goldMinNumClosestFacets, closestFacets.size());
    EXPECT_GE(goldMaxNumClosestFacets, closestFacets.size());
}

TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, facetedSurfaceWith2FacetsPerSide_getClosestFacets_searchLimitsResultsToActualClosestFacets)
{
    test_num_closest_facets_with_num_facets_per_side(2, stk::math::Vector3d{0,0,0}, 6, 8);
    test_num_closest_facets_with_num_facets_per_side(8, stk::math::Vector3d{0,0,0}, 6, 8);
    test_num_closest_facets_with_num_facets_per_side(2, stk::math::Vector3d{5,0,-5}, 1, 2);
    test_num_closest_facets_with_num_facets_per_side(8, stk::math::Vector3d{5,0,-5}, 1, 2);
    test_num_closest_facets_with_num_facets_per_side(2, stk::math::Vector3d{0,0,-5}, 3, 4);
    test_num_closest_facets_with_num_facets_per_side(8, stk::math::Vector3d{0,0,-5}, 3, 4);
    test_num_closest_facets_with_num_facets_per_side(2, stk::math::Vector3d{-5,0,-5}, 2, 2);
    test_num_closest_facets_with_num_facets_per_side(8, stk::math::Vector3d{-5,0,-5}, 2, 2);
}


void test_num_candidate_intersecting_facets_with_num_facets_per_side(const unsigned numFacetsPerSide, const std::array<stk::math::Vector3d,2> & segmentLocs, const unsigned goldMinNumClosestFacets, const unsigned goldMaxNumClosestFacets)
{
    const TestFacetedSurface surf = build_single_square_surface_with_num_facets_per_side_geometry(numFacetsPerSide);
    SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

    krino::BoundingBox segmentBbox;
    segmentBbox.accommodate(segmentLocs[0]);
    segmentBbox.accommodate(segmentLocs[1]);

    std::vector<size_t> intersectingFacets;
    facetTree.get_intersecting_entities(segmentBbox, intersectingFacets);

    EXPECT_LE(goldMinNumClosestFacets, intersectingFacets.size());
    EXPECT_GE(goldMaxNumClosestFacets, intersectingFacets.size());
}

TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, facetedSurfaceWithNumFacetsPerSide_getIntersectingFacets_searchLimitsResultsToActualIntersectingFacets)
{
    test_num_candidate_intersecting_facets_with_num_facets_per_side(2, {{stk::math::Vector3d{0,-2,0,}, stk::math::Vector3d{0,-1,0}}}, 0, 0);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(4, {{stk::math::Vector3d{0,-2,0,}, stk::math::Vector3d{0,-1,0}}}, 0, 0);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(2, {{stk::math::Vector3d{0,-1,0}, stk::math::Vector3d{0,1,0}}}, 6, 8);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(4, {{stk::math::Vector3d{0,-1,0}, stk::math::Vector3d{0,1,0}}}, 6, 8);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(2, {{stk::math::Vector3d{4,-1,-4}, stk::math::Vector3d{4,1,-4}}}, 1, 2);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(4, {{stk::math::Vector3d{4,-1,-4}, stk::math::Vector3d{4,1,-4}}}, 1, 2);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(2, {{stk::math::Vector3d{0,-1,-4}, stk::math::Vector3d{0,1,-4}}}, 2, 4);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(4, {{stk::math::Vector3d{0,-1,-4}, stk::math::Vector3d{0,1,-4}}}, 2, 4);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(2, {{stk::math::Vector3d{-4,-1,-4}, stk::math::Vector3d{-4,1,-4}}}, 2, 2);
    test_num_candidate_intersecting_facets_with_num_facets_per_side(4, {{stk::math::Vector3d{-4,-1,-4}, stk::math::Vector3d{-4,1,-4}}}, 2, 2);
}

TEST(TestFacetTreeQueriesWithLeafEntityDistance, facetWithDoubleDistanceLessThanBoundingBoxDistance_stillGetClosestFacet)
{
  // This test drives a fix for an interesting corner case:
  // The distance to facet {0,1,2} in double precision is less than the distance to its
  // bounding box, stored with floats.  The fix was to make sure the upperBnd2 is always
  // greater than or equal to the lowerBnd2 in the SearchTree.
  std::vector<stk::math::Vector3d> nodeLocs{
    {-2.49199999999999999289,2.19584979920073175563,-2.64111782888453072715},
    {-2.49199999999999999289,2.19426826798352747616,-2.64711908728139500013},
    {-2.49199999999999999289,2.19672629465853397335,-2.64962435720086819657},
    {-2.49199999999999999289,2.1934,-2.6386126}};
  std::vector<std::array<size_t,3>> facetNodes{{0,1,2},{0,2,3}};
  TestFacetedSurface surf(nodeLocs, facetNodes);
  SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

  stk::math::Vector3d queryLoc{-2.47350389205193943454,2.19619679237152309881,-2.64762050545203919683};
  std::vector<size_t> closestFacets;
  facetTree.find_closest_entities(queryLoc, build_get_facet_distance_squared(surf), closestFacets);

  EXPECT_FALSE(closestFacets.empty());
}

// Brute-force nearest facet, used as ground truth to validate the tree query.
size_t brute_force_closest_facet(const TestFacetedSurface & surf, const stk::math::Vector3d & queryLoc)
{
  size_t closestFacet = 0;
  double closestDist2 = std::numeric_limits<double>::max();
  for (const size_t facet : surf.get_facets())
  {
    const double dist2 = CalcTriangle3<double>::distance_squared(surf.get_facet_node_locations(facet), queryLoc);
    if (dist2 < closestDist2)
    {
      closestDist2 = dist2;
      closestFacet = facet;
    }
  }
  return closestFacet;
}

// Correctness: the candidate set returned by the tree must always contain the
// true nearest facet (as determined by brute force). The existing tests only
// bound the size of the candidate set, not its contents, so this guards the
// actual result against regressions when the bounds logic is changed.
void expect_tree_candidates_contain_true_closest_facet(const TestFacetedSurface & surf,
    SearchTree<size_t> & facetTree,
    const stk::math::Vector3d & queryLoc)
{
  const auto getDist2 = build_get_facet_distance_squared(surf);

  std::vector<size_t> candidateFacets;
  facetTree.find_closest_entities(queryLoc, getDist2, candidateFacets);
  ASSERT_FALSE(candidateFacets.empty());

  const size_t goldClosestFacet = brute_force_closest_facet(surf, queryLoc);
  const double goldClosestDist2 = getDist2(goldClosestFacet, queryLoc);

  EXPECT_NE(candidateFacets.end(), std::find(candidateFacets.begin(), candidateFacets.end(), goldClosestFacet))
      << "Tree candidate set is missing the true closest facet " << goldClosestFacet;

  // The best distance among the candidates must equal the brute-force minimum.
  double treeClosestDist2 = std::numeric_limits<double>::max();
  for (const size_t facet : candidateFacets)
    treeClosestDist2 = std::min(treeClosestDist2, getDist2(facet, queryLoc));
  EXPECT_NEAR(goldClosestDist2, treeClosestDist2, 1.e-10);
}

TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, findClosestEntities_candidateSetAlwaysContainsTrueClosestFacet)
{
  const TestFacetedSurface surf = build_single_square_surface_with_num_facets_per_side_geometry(16);
  SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

  // Sweep query points on, above, and beside the surface, including off the edges.
  std::mt19937 randomEngine(42);
  std::uniform_real_distribution<double> coord(-7.0, 7.0);
  for (int i = 0; i < 200; ++i)
  {
    const stk::math::Vector3d queryLoc(coord(randomEngine), coord(randomEngine), coord(randomEngine));
    expect_tree_candidates_contain_true_closest_facet(surf, facetTree, queryLoc);
  }
}

// Brute-force set of facets whose (per-facet) bounding box intersects the query
// bbox -- ground truth for get_intersecting_entities, which returns a candidate
// superset based on the same per-facet bounding boxes.
std::vector<size_t> brute_force_intersecting_facets(const TestFacetedSurface & surf, const krino::BoundingBox & queryBbox)
{
  std::vector<size_t> intersecting;
  for (const size_t facet : surf.get_facets())
  {
    krino::BoundingBox facetBbox;
    surf.insert_facet_into_bounding_box(facet, facetBbox);
    if (facetBbox.intersects(queryBbox))
      intersecting.push_back(facet);
  }
  return intersecting;
}

// Correctness: the candidate set from get_intersecting_entities must contain
// every facet whose bounding box truly intersects the query box. The existing
// tests only bound the candidate count, not its contents.
void expect_tree_candidates_contain_all_truly_intersecting_facets(const TestFacetedSurface & surf,
    SearchTree<size_t> & facetTree,
    const krino::BoundingBox & queryBbox)
{
  std::vector<size_t> candidateFacets;
  facetTree.get_intersecting_entities(queryBbox, candidateFacets);

  const std::vector<size_t> goldIntersecting = brute_force_intersecting_facets(surf, queryBbox);
  for (const size_t goldFacet : goldIntersecting)
    EXPECT_NE(candidateFacets.end(), std::find(candidateFacets.begin(), candidateFacets.end(), goldFacet))
        << "Tree candidate set is missing truly-intersecting facet " << goldFacet;
}

TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, getIntersectingEntities_candidateSetContainsAllTrulyIntersectingFacets)
{
  const TestFacetedSurface surf = build_single_square_surface_with_num_facets_per_side_geometry(16);
  SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

  std::mt19937 randomEngine(7);
  std::uniform_real_distribution<double> coord(-7.0, 7.0);
  std::uniform_real_distribution<double> halfSpan(0.1, 3.0);
  for (int i = 0; i < 200; ++i)
  {
    const stk::math::Vector3d center(coord(randomEngine), coord(randomEngine), coord(randomEngine));
    const stk::math::Vector3d halfExtent(halfSpan(randomEngine), halfSpan(randomEngine), halfSpan(randomEngine));
    krino::BoundingBox queryBbox;
    queryBbox.accommodate(center - halfExtent);
    queryBbox.accommodate(center + halfExtent);
    expect_tree_candidates_contain_all_truly_intersecting_facets(surf, facetTree, queryBbox);
  }
}

// Non-zero max_search_radius: a radius that excludes the surface entirely must
// yield no candidates, while a generous radius must still find the true closest.
TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, findClosestEntities_maxSearchRadiusLimitsResults)
{
  const TestFacetedSurface surf = build_single_square_surface_with_num_facets_per_side_geometry(8);
  SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));
  const auto getDist2 = build_get_facet_distance_squared(surf);

  // Surface lies in the y=0 plane; query is 10 units away in y.
  const stk::math::Vector3d queryLoc{0, 10, 0};

  // Radius too small to reach the surface: expect nothing.
  std::vector<size_t> closestFacets;
  facetTree.find_closest_entities(queryLoc, getDist2, closestFacets, /*max_search_radius=*/5.0);
  EXPECT_TRUE(closestFacets.empty());

  // Radius large enough to reach the surface: expect the true closest is present.
  facetTree.find_closest_entities(queryLoc, getDist2, closestFacets, /*max_search_radius=*/20.0);
  ASSERT_FALSE(closestFacets.empty());
  const size_t goldClosestFacet = brute_force_closest_facet(surf, queryLoc);
  EXPECT_NE(closestFacets.end(), std::find(closestFacets.begin(), closestFacets.end(), goldClosestFacet));
}

// Single-entity tree: exercises the root-is-leaf path (build with entity_size==1
// and the 0==index special case in recursively_update_tree_bounds).
TEST(TestFacetTreeQueriesUniformSpacedSurfaceFacets, singleFacetTree_findClosestAndIntersecting)
{
  std::vector<stk::math::Vector3d> nodeLocs{{0,0,0}, {1,0,0}, {0,0,1}};
  std::vector<std::array<size_t,3>> facetNodes{{0,1,2}};
  const TestFacetedSurface surf(nodeLocs, facetNodes);
  SearchTree<size_t> facetTree(surf.get_facets(), build_get_facet_centroid_function(surf), build_insert_facet_into_bounding_box_function(surf));

  EXPECT_EQ(1u, facetTree.size());

  std::vector<size_t> closestFacets;
  facetTree.find_closest_entities(stk::math::Vector3d{0.25,3,0.25}, build_get_facet_distance_squared(surf), closestFacets);
  ASSERT_EQ(1u, closestFacets.size());
  EXPECT_EQ(0u, closestFacets[0]);

  krino::BoundingBox queryBbox;
  queryBbox.accommodate(stk::math::Vector3d{-1,-1,-1});
  queryBbox.accommodate(stk::math::Vector3d{2,2,2});
  std::vector<size_t> intersectingFacets;
  facetTree.get_intersecting_entities(queryBbox, intersectingFacets);
  ASSERT_EQ(1u, intersectingFacets.size());
  EXPECT_EQ(0u, intersectingFacets[0]);

  // A query box far from the facet must return nothing.
  krino::BoundingBox farBbox;
  farBbox.accommodate(stk::math::Vector3d{10,10,10});
  farBbox.accommodate(stk::math::Vector3d{11,11,11});
  facetTree.get_intersecting_entities(farBbox, intersectingFacets);
  EXPECT_TRUE(intersectingFacets.empty());
}

}
