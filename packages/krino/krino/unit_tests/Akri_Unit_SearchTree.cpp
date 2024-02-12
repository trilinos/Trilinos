#include <gtest/gtest.h>
#include <numeric>
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

static std::function<void(const size_t, krino::BoundingBox & bbox)> build_insert_facet_into_bounding_box_function(const TestFacetedSurface & surf)
{
    auto fn = [&surf](const size_t facet, krino::BoundingBox & bbox) { return surf.insert_facet_into_bounding_box(facet, bbox); };
    return fn;
}

static std::function<stk::math::Vector3d(const size_t)> build_get_facet_centroid_function(const TestFacetedSurface & surf)
{
    auto fn = [&surf](const size_t facet) { return surf.get_facet_centroid(facet); };
    return fn;
}

static std::function<double(const size_t, const stk::math::Vector3d &)> build_get_facet_distance_squared(const TestFacetedSurface & surf)
{
    auto fn = [&surf](const size_t facet, const stk::math::Vector3d & queryLoc)
    {
      const auto & facetNodeLocs = surf.get_facet_node_locations(facet);
      return CalcTriangle3<double>::distance_squared(facetNodeLocs, queryLoc);
    };
    return fn;
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

}
