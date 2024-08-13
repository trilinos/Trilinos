#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <Akri_WindingNumber.hpp>
#include <stk_util/environment/WallTime.hpp>

namespace krino {

void expect_winding_number_for_facet(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc, const double goldWindingNumber)
{
  EXPECT_NEAR(goldWindingNumber, compute_facet_winding_number(x0, x1, x2, queryLoc), 1.e-6);
}

void expect_winding_number_in_range_for_point_in_plane_of_facet(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc)
{
  const double windingNumber = compute_facet_winding_number(x0, x1, x2, queryLoc);
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
  const std::vector<std::array<stk::math::Vector3d,3>> surfFacets = initialize_nonplanar_facets(6);

  const unsigned dim=20;
  std::vector<stk::math::Vector3d> queryLocs;
  for (unsigned i=0; i<dim; ++i)
    for (unsigned j=0; j<dim; ++j)
      for (unsigned k=0; k<dim; ++k)
        queryLocs.emplace_back(-10. + 20.*i/dim, -10. + 20.*j/dim, -10. + 20.*k/dim);

  test_performance_for_winding_number(surfFacets, queryLocs);
}

}

