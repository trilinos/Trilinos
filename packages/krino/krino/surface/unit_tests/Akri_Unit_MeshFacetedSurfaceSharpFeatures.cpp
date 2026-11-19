// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_MeshFacetedSurface.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <algorithm>
#include <sstream>
#include <string>
#include <array>
#include <vector>

namespace krino
{

namespace
{

// The facets with index i such that i % numProcs == rank, so every facet is owned by one processor.
std::vector<std::array<unsigned, 3>>
facets_owned_by_this_processor(const std::vector<std::array<unsigned, 3>> & allFacets)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<std::array<unsigned, 3>> myFacets;
  for (size_t i = 0; i < allFacets.size(); ++i)
    if (int(i % numProcs) == me) myFacets.push_back(allFacets[i]);
  return myFacets;
}

// The 12 triangles of the box with these 8 vertices, ordered as Cuboid::vertex_location orders
// them (the bottom 4 counterclockwise, then the top 4).
std::vector<std::array<unsigned, 3>> box_facets()
{
  std::vector<std::array<unsigned, 3>> facets{{{0, 2, 1}}, {{0, 3, 2}}, {{4, 5, 6}}, {{4, 6, 7}}};
  for (unsigned i = 0; i < 4; ++i)
  {
    const unsigned j = (i + 1) % 4;
    facets.push_back({{i, j, 4 + j}});
    facets.push_back({{i, 4 + j, 4 + i}});
  }
  return facets;
}

// The tets of a lattice of n^3 cells over the box, each cell cut into 6 tets around its diagonal.
// The lattice is offset by irrational fractions so that no tet edge or face is on the cuboid's
// edges or corners.
std::vector<std::array<stk::math::Vector3d, 4>>
lattice_tets(const BoundingBox & box, const unsigned n)
{
  const stk::math::Vector3d lo(box.get_min()[0], box.get_min()[1], box.get_min()[2]);
  const stk::math::Vector3d hi(box.get_max()[0], box.get_max()[1], box.get_max()[2]);
  auto lattice_point = [&](const unsigned i, const unsigned j, const unsigned k)
  {
    // Within the box: the offsets are less than the 0.1 cells it is shrunk by.
    const double cells = n + 0.1;
    return stk::math::Vector3d(lo[0] + (i + 0.0381966) * (hi[0] - lo[0]) / cells,
        lo[1] + (j + 0.0236068) * (hi[1] - lo[1]) / cells,
        lo[2] + (k + 0.0618034) * (hi[2] - lo[2]) / cells);
  };
  std::vector<std::array<stk::math::Vector3d, 4>> tets;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      for (unsigned k = 0; k < n; ++k)
      {
        auto v = [&](const unsigned a, const unsigned b, const unsigned c)
        { return lattice_point(i + a, j + b, k + c); };
        tets.push_back({{v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(1, 1, 1)}});
        tets.push_back({{v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(1, 1, 1)}});
        tets.push_back({{v(0, 0, 0), v(0, 1, 0), v(1, 1, 0), v(1, 1, 1)}});
        tets.push_back({{v(0, 0, 0), v(0, 1, 0), v(0, 1, 1), v(1, 1, 1)}});
        tets.push_back({{v(0, 0, 0), v(0, 0, 1), v(1, 0, 1), v(1, 1, 1)}});
        tets.push_back({{v(0, 0, 0), v(0, 0, 1), v(0, 1, 1), v(1, 1, 1)}});
      }
  return tets;
}

std::string to_string(const std::vector<stk::math::Vector3d> & points)
{
  std::ostringstream os;
  for (auto & point : points)
    os << " (" << point[0] << ", " << point[1] << ", " << point[2] << ")";
  return os.str();
}

template <size_t N>
void expect_same_parametric_coordinates(const std::vector<stk::math::Vector3d> & expected,
    const std::vector<stk::math::Vector3d> & actual,
    const std::array<stk::math::Vector3d, N> & simplex)
{
  // Each expected point matches a different actual one (in any order: points can differ in the last
  // bits, so sorting them is not reliable).
  bool isSame = expected.size() == actual.size();
  std::vector<bool> isMatched(actual.size(), false);
  for (size_t i = 0; isSame && i < expected.size(); ++i)
  {
    isSame = false;
    for (size_t j = 0; !isSame && j < actual.size(); ++j)
      if (!isMatched[j] && (expected[i] - actual[j]).length() < 1.e-12)
        isSame = isMatched[j] = true;
  }
  EXPECT_TRUE(isSame) << "expected" << to_string(expected) << ", actual" << to_string(actual)
                      << ", for the simplex"
                      << to_string(
                             std::vector<stk::math::Vector3d>(simplex.begin(), simplex.end()));
}

size_t sum_over_processors(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// For each lattice tet in queryBox, the facet box gives the same corners in the tet and feature
// edges crossing its faces as the analytic cuboid.  Returns the numbers of intersections found.
std::array<size_t, 2> expect_same_intersections_as_cuboid(const Cuboid & cuboid,
    const MeshFacetedSurface & surface,
    const std::vector<std::array<stk::math::Vector3d, 4>> & tets,
    const BoundingBox & queryBox)
{
  std::array<size_t, 2> numIntersections{{0, 0}};
  std::vector<stk::math::Vector3d> expected, actual;
  for (auto & tet : tets)
  {
    if (!std::all_of(tet.begin(),
            tet.end(),
            [&](const stk::math::Vector3d & x) { return queryBox.contains(x); }))
      continue;
    cuboid.fill_tetrahedon_intersection_parametric_coordinates(tet, expected);
    surface.fill_tetrahedon_intersection_parametric_coordinates(tet, actual);
    expect_same_parametric_coordinates(expected, actual, tet);
    numIntersections[1] += actual.size();
    for (const std::array<unsigned, 3> & face : {std::array<unsigned, 3>{{0, 1, 3}},
             std::array<unsigned, 3>{{1, 2, 3}},
             std::array<unsigned, 3>{{0, 3, 2}},
             std::array<unsigned, 3>{{0, 2, 1}}})
    {
      const std::array<stk::math::Vector3d, 3> tri{{tet[face[0]], tet[face[1]], tet[face[2]]}};
      cuboid.fill_triangle_intersection_parametric_coordinates(tri, expected);
      surface.fill_triangle_intersection_parametric_coordinates(tri, actual);
      expect_same_parametric_coordinates(expected, actual, tri);
      numIntersections[0] += actual.size();
    }
  }
  return numIntersections;
}

const BoundingBox latticeBox(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
const Cuboid cuboid(stk::math::Vector3d(0.1, -0.05, 0.2), stk::math::Vector3d(1.1, 0.9, 1.3));

MeshFacetedSurface build_facet_box_split_across_processors()
{
  std::vector<stk::math::Vector3d> vertices;
  for (unsigned i = 0; i < 8; ++i)
    vertices.push_back(cuboid.vertex_location(i));
  return MeshFacetedSurface(MPI_COMM_WORLD, 3, vertices, facets_owned_by_this_processor(box_facets()));
}

} // namespace

// A connected box of facets, split across processors, with a feature angle less than its right
// angles: its feature edges cross the faces of a lattice of tets, and its corners are in the tets,
// exactly where those of the same analytic cuboid are.
TEST(MeshFacetedSurfaceSharpFeatures, facetBoxSplitAcrossProcessors_sameIntersectionsAsCuboid)
{
  MeshFacetedSurface surface = build_facet_box_split_across_processors();
  surface.set_sharp_feature_angle(135.);
  EXPECT_TRUE(surface.has_edges());
  EXPECT_TRUE(surface.has_corners());
  surface.prepare_to_compute(0., latticeBox, 0.);

  const auto [numFaceIntersections, numTetIntersections] =
      expect_same_intersections_as_cuboid(cuboid, surface, lattice_tets(latticeBox, 4), latticeBox);
  EXPECT_LT(0u, numFaceIntersections);
  EXPECT_EQ(8u, numTetIntersections); // each corner is in one tet
}

// Each processor queries a slab of the lattice: it is given the features near its slab, including
// those found on other processors, and finds the same intersections as the cuboid there.
TEST(MeshFacetedSurfaceSharpFeatures,
    facetBoxQueriedInSlabs_sameIntersectionsAsCuboidWithFeaturesFromOtherProcessors)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const double slabWidth = 2. / numProcs;
  const BoundingBox slab(stk::math::Vector3d(-1. + me * slabWidth, -1., -1.),
      stk::math::Vector3d(-1. + (me + 1) * slabWidth, 1., 1.));

  MeshFacetedSurface surface = build_facet_box_split_across_processors();
  surface.set_sharp_feature_angle(135.);
  surface.prepare_to_compute(0., slab, 0.);

  const auto [numFaceIntersections, numTetIntersections] =
      expect_same_intersections_as_cuboid(cuboid, surface, lattice_tets(latticeBox, 4), slab);
  EXPECT_LT(0u, sum_over_processors(numFaceIntersections));
  EXPECT_LT(0u, sum_over_processors(numTetIntersections));
}

// Every edge of a lattice of tets, either way round, crosses the facet box where it crosses the
// analytic cuboid, including twice where it clips a feature edge or corner.
TEST(MeshFacetedSurfaceSharpFeatures, facetBoxSplitAcrossProcessors_sameSegmentCrossingsAsCuboid)
{
  MeshFacetedSurface surface = build_facet_box_split_across_processors();
  surface.prepare_to_compute(0., latticeBox, 0.);

  const double edgeCrossingTol = 1.e-10;
  std::vector<double> expected, actual;
  size_t numCrossings = 0, numSegmentsCrossedTwice = 0;
  for (auto & tet : lattice_tets(latticeBox, 4))
    for (unsigned i = 0; i < 4; ++i)
      for (unsigned j = 0; j < 4; ++j)
      {
        if (i == j) continue;
        cuboid.fill_segment_crossings(tet[i], tet[j], edgeCrossingTol, expected);
        surface.fill_segment_crossings(tet[i], tet[j], edgeCrossingTol, actual);
        ASSERT_EQ(expected.size(), actual.size()) << tet[i] << " to " << tet[j];
        for (size_t k = 0; k < expected.size(); ++k)
        {
          EXPECT_NEAR(expected[k], actual[k], 1.e-12) << tet[i] << " to " << tet[j];
        }
        numCrossings += actual.size();
        if (actual.size() == 2) ++numSegmentsCrossedTwice;
      }
  EXPECT_LT(0u, numCrossings);
  EXPECT_LT(0u, numSegmentsCrossedTwice);
}

// A segment through the diagonal of a face, where the two triangles of the face meet, crosses the
// surface once.
TEST(MeshFacetedSurfaceSharpFeatures, facetBoxSegmentThroughFaceDiagonal_oneCrossing)
{
  MeshFacetedSurface surface = build_facet_box_split_across_processors();
  surface.prepare_to_compute(0., latticeBox, 0.);

  const stk::math::Vector3d diagonalMidpoint =
      0.5 * (cuboid.vertex_location(0) + cuboid.vertex_location(2));
  const stk::math::Vector3d pt0 = diagonalMidpoint - stk::math::Vector3d(0., 0., 0.1);
  const stk::math::Vector3d pt1 = diagonalMidpoint + stk::math::Vector3d(0., 0., 0.3);
  std::vector<double> crossings;
  surface.fill_segment_crossings(pt0, pt1, 1.e-10, crossings);
  ASSERT_EQ(1u, crossings.size());
  EXPECT_NEAR(0.25, crossings[0], 1.e-12);
}

TEST(MeshFacetedSurfaceSharpFeatures, withoutFeatureAngle_noSharpFeatures)
{
  MeshFacetedSurface surface = build_facet_box_split_across_processors();
  EXPECT_FALSE(surface.has_edges());
  EXPECT_FALSE(surface.has_corners());
  surface.prepare_to_compute(0., latticeBox, 0.);
  std::vector<stk::math::Vector3d> intParamCoords{stk::math::Vector3d::ZERO};
  surface.fill_tetrahedon_intersection_parametric_coordinates(
      {{cuboid.vertex_location(0) - stk::math::Vector3d(0.1, 0.1, 0.1),
          {1., 0., 0.},
          {0., 1., 0.},
          {0., 0., 1.}}},
      intParamCoords);
  EXPECT_TRUE(intParamCoords.empty());
}

// A square plate's sides cross a triangle only if its open boundary edges are features.
TEST(MeshFacetedSurfaceSharpFeatures, openPlate_boundaryEdgesCrossTrianglesOnlyIfFeatures)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
  const std::vector<std::array<unsigned, 3>> facets{{{0, 1, 2}}, {{0, 2, 3}}};
  const std::array<stk::math::Vector3d, 3> triCrossingSide{
      {{0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {0.5, 0., 0.5}}};
  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(2., 2., 1.));
  for (const bool areOpenBoundaryEdgesFeatures : {false, true})
  {
    MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_owned_by_this_processor(facets));
    surface.set_sharp_feature_angle(135., areOpenBoundaryEdgesFeatures);
    surface.prepare_to_compute(0., box, 0.);
    std::vector<stk::math::Vector3d> intParamCoords;
    surface.fill_triangle_intersection_parametric_coordinates(triCrossingSide, intParamCoords);
    EXPECT_EQ(areOpenBoundaryEdgesFeatures ? 1u : 0u, intParamCoords.size());
  }
}

// In 2D the corners of a square loop of segments are in the triangles that contain them, and it
// has no feature edges.
TEST(MeshFacetedSurfaceSharpFeatures, squareLoop2d_cornersInTriangles)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
  const std::vector<std::array<unsigned, 3>> segments{
      {{0, 1, 0}}, {{1, 2, 0}}, {{2, 3, 0}}, {{3, 0, 0}}};
  MeshFacetedSurface surface(MPI_COMM_WORLD, 2, vertices, facets_owned_by_this_processor(segments));
  surface.set_sharp_feature_angle(135.);
  EXPECT_FALSE(surface.has_edges());
  EXPECT_TRUE(surface.has_corners());
  surface.prepare_to_compute(
      0., BoundingBox(stk::math::Vector3d(-1., -1., 0.), stk::math::Vector3d(2., 2., 0.)), 0.);

  std::vector<stk::math::Vector3d> intParamCoords;
  surface.fill_triangle_intersection_parametric_coordinates(
      {{{0.5, 0.5, 0.}, {1.5, 0.5, 0.}, {0.5, 1.5, 0.}}}, intParamCoords);
  ASSERT_EQ(1u, intParamCoords.size()); // the corner (1,1)
  EXPECT_NEAR(0., (stk::math::Vector3d(0.5, 0.5, 0.) - intParamCoords[0]).length(), 1.e-12);
  surface.fill_triangle_intersection_parametric_coordinates(
      {{{0.2, 0.2, 0.}, {0.8, 0.2, 0.}, {0.2, 0.8, 0.}}}, intParamCoords);
  EXPECT_TRUE(intParamCoords.empty());
}

TEST(MeshFacetedSurfaceSharpFeatures, unconnectedSurface_settingFeatureAngleThrows)
{
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, std::vector<stk::math::Vector3d>{});
  EXPECT_ANY_THROW(surface.set_sharp_feature_angle(135.));
}

} // namespace krino
