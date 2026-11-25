// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_BoundingBox.hpp>
#include <Akri_FacetCoordinates.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

namespace krino
{

namespace
{

constexpr double segmentTol = 1.e-12;

// The facets with index i such that i % numProcs == rank, so every facet is owned by one processor.
std::vector<std::array<unsigned, 3>>
my_facets_for_this_processor(const std::vector<std::array<unsigned, 3>> & allFacets)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<std::array<unsigned, 3>> myFacets;
  for (size_t i = 0; i < allFacets.size(); ++i)
    if (int(i % numProcs) == me) myFacets.push_back(allFacets[i]);
  return myFacets;
}

std::vector<stk::math::Vector3d> segment_facet_vertex_coordinates(const unsigned dim,
    const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facets)
{
  std::vector<stk::math::Vector3d> coords;
  for (auto & facet : facets)
    for (unsigned n = 0; n < dim; ++n)
      coords.push_back(vertices[facet[n]]);
  return coords;
}

// The facets split across processors, as a connected surface and as a soup.
std::vector<std::unique_ptr<MeshFacetedSurface>>
build_split_connected_and_soup_surfaces(const unsigned dim,
    const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & allFacets)
{
  const std::vector<std::array<unsigned, 3>> myFacets = my_facets_for_this_processor(allFacets);
  std::vector<std::unique_ptr<MeshFacetedSurface>> surfaces;
  surfaces.push_back(std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, dim, vertices, myFacets));
  surfaces.push_back(std::make_unique<MeshFacetedSurface>(
      MPI_COMM_WORLD, dim, segment_facet_vertex_coordinates(dim, vertices, myFacets)));
  return surfaces;
}

size_t segment_global_sum(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// Directions spread over the unit sphere (a Fibonacci spiral).
std::vector<stk::math::Vector3d> directions_on_sphere(const unsigned num)
{
  const double goldenAngle = std::numbers::pi * (3. - std::sqrt(5.));
  std::vector<stk::math::Vector3d> directions;
  for (unsigned i = 0; i < num; ++i)
  {
    const double z = 1. - 2. * (i + 0.5) / num;
    const double r = std::sqrt(1. - z * z);
    directions.emplace_back(r * std::cos(goldenAngle * i), r * std::sin(goldenAngle * i), z);
  }
  return directions;
}

// Expects the segment from pt0 to pt1 to cross the surface at its midpoint, with the sign of pt1,
// and the reversed segment to cross there with the opposite sign.
void expect_crossing_at_midpoint(const MeshFacetedSurface & surface,
    const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const int sign)
{
  const auto [forwardSign, forwardLoc] = surface.compute_intersection_with_segment(pt0, pt1, 0.);
  EXPECT_EQ(sign, forwardSign) << pt0 << " to " << pt1;
  EXPECT_NEAR(0.5, forwardLoc, segmentTol) << pt0 << " to " << pt1;
  const auto [reverseSign, reverseLoc] = surface.compute_intersection_with_segment(pt1, pt0, 0.);
  EXPECT_EQ(-sign, reverseSign) << pt1 << " to " << pt0;
  EXPECT_NEAR(0.5, reverseLoc, segmentTol) << pt1 << " to " << pt0;
}

} // namespace

// Each processor queries a slab of x with a narrow band, as CDFEM does, so the facets crossed by
// some of its segments are owned by other processors.  Split across processors, connected or as a
// soup, the surface gives each segment the crossing (or none) it has with all the facets on every
// processor, and the crossings are on the sphere.
TEST(MeshFacetedSurfaceSegmentIntersection,
    sphereSplitAcrossProcessors_crossingsSameAsAllFacetsOnEveryProcessor)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const double slabWidth = 2. * radius / numProcs;
  const double xLo = -radius + me * slabWidth;
  const BoundingBox localQueryBox(
      stk::math::Vector3d(xLo, -1., -1.), stk::math::Vector3d(xLo + slabWidth, 1., 1.));
  const double narrowBand = 0.1;

  // Every facet on every processor: no ghosting is needed.
  MeshFacetedSurface allOnEveryProc(
      MPI_COMM_WORLD, 3, segment_facet_vertex_coordinates(3, vertices, allFacets));
  allOnEveryProc.prepare_to_compute(0., localQueryBox, narrowBand);

  // Radial segments, of length less than the narrow band, that cross the sphere from inside to
  // outside, and that stay inside it.
  std::vector<std::pair<stk::math::Vector3d, stk::math::Vector3d>> segments;
  for (auto & dir : directions_on_sphere(200))
    for (auto [r0, r1] : {std::pair{0.7, 0.8}, std::pair{0.6, 0.65}})
      if (localQueryBox.contains(r0 * dir) && localQueryBox.contains(r1 * dir))
        segments.emplace_back(r0 * dir, r1 * dir);

  for (auto & surface : build_split_connected_and_soup_surfaces(3, vertices, allFacets))
  {
    surface->prepare_to_compute(0., localQueryBox, narrowBand);
    if (numProcs > 1)
    {
      EXPECT_GT(segment_global_sum(surface->num_ghosted_facets()), 0u);
    }

    size_t numCrossings = 0;
    for (auto & [pt0, pt1] : segments)
    {
      const auto [expectedSign, expectedLoc] =
          allOnEveryProc.compute_intersection_with_segment(pt0, pt1, 0.);
      const auto [crossingSign, loc] = surface->compute_intersection_with_segment(pt0, pt1, 0.);
      EXPECT_EQ(expectedSign, crossingSign) << pt0 << " to " << pt1;
      EXPECT_DOUBLE_EQ(expectedLoc, loc) << pt0 << " to " << pt1;
      if (0 != crossingSign)
      {
        ++numCrossings;
        EXPECT_EQ(1, crossingSign);
        EXPECT_NEAR(radius, ((1. - loc) * pt0 + loc * pt1).length(), 0.02);
      }
      else
      {
        EXPECT_LT(pt1.length(), radius);
      }
    }
    EXPECT_LT(0u, numCrossings);
  }
}

// A connected or soup cube split across processors: segments through an edge, a face diagonal
// and a corner cross there, with the sign of their far end.
TEST(MeshFacetedSurfaceSegmentIntersection, cubeEdgesAndCorner_crossingThereWithSignOfTheEnd)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  // Outward facets; each side is split along the diagonal from its bottom-left to top-right corner.
  const std::vector<std::array<unsigned, 3>> allFacets{{{0, 2, 1}},
      {{0, 3, 2}},
      {{4, 5, 6}},
      {{4, 6, 7}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{2, 3, 7}},
      {{2, 7, 6}},
      {{3, 0, 4}},
      {{3, 4, 7}}};
  const BoundingBox box(stk::math::Vector3d(-0.5, -0.5, -0.5), stk::math::Vector3d(1.5, 1.5, 1.5));

  for (auto & surface : build_split_connected_and_soup_surfaces(3, vertices, allFacets))
  {
    surface->prepare_to_compute(0., box, 0.);
    // Through the vertical edge at x = y = 1.
    expect_crossing_at_midpoint(*surface, {0.9, 0.9, 0.5}, {1.1, 1.1, 0.5}, 1);
    // Through the top edge at y = z = 1.
    expect_crossing_at_midpoint(*surface, {0.5, 0.9, 0.9}, {0.5, 1.1, 1.1}, 1);
    // Through the diagonal edge of the x = 1 side, from (1,0,0) to (1,1,1).
    expect_crossing_at_midpoint(*surface, {0.9, 0.5, 0.5}, {1.1, 0.5, 0.5}, 1);
    // Through the corner (1,1,1).
    expect_crossing_at_midpoint(*surface, {0.9, 0.9, 0.9}, {1.1, 1.1, 1.1}, 1);
  }
}

// Two halves of a square pulled apart along their shared diagonal: a segment through the gap
// crosses no facet, but the signed distance changes sign along it.  It is taken to cross where
// the signed distance does when the surface is near there (within 5% of the segment's length), as
// for a small gap, and not otherwise.
TEST(MeshFacetedSurfaceSegmentIntersection, gapBetweenFacets_crossingFromSignedDistanceIfNear)
{
  const stk::math::Vector3d pt0(0.5, 0.5, -0.5);
  const stk::math::Vector3d pt1(0.5, 0.5, 0.5);
  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(2., 2., 1.));
  for (const double gap : {1.e-3, 0.2})
  {
    // Each half moved gap/2 away from the diagonal, so the gap is gap wide.
    const stk::math::Vector3d shift = (0.5 * gap / std::sqrt(2.)) * stk::math::Vector3d(1., 1., 0.);
    const std::vector<stk::math::Vector3d> vertices{stk::math::Vector3d(0, 0, 0) - shift,
        stk::math::Vector3d(1, 0, 0) - shift,
        stk::math::Vector3d(0, 1, 0) - shift,
        stk::math::Vector3d(1, 1, 0) + shift,
        stk::math::Vector3d(0, 1, 0) + shift,
        stk::math::Vector3d(1, 0, 0) + shift};
    const std::vector<std::array<unsigned, 3>> allFacets{{{0, 1, 2}}, {{3, 4, 5}}};

    for (auto & surface : build_split_connected_and_soup_surfaces(3, vertices, allFacets))
    {
      surface->prepare_to_compute(0., box, 0.);
      EXPECT_GT(0., surface->point_signed_distance(pt0));
      EXPECT_LT(0., surface->point_signed_distance(pt1));
      const auto [crossingSign, loc] = surface->compute_intersection_with_segment(pt0, pt1, 0.);
      if (gap < 0.1)
      {
        EXPECT_EQ(1, crossingSign) << "gap " << gap;
        EXPECT_NEAR(0.5, loc, 1.e-3) << "gap " << gap;
      }
      else
      {
        EXPECT_EQ(0, crossingSign) << "gap " << gap;
        EXPECT_GT(0., loc) << "gap " << gap;
      }
    }
  }
}

// A square loop of segments, with a vertex at the middle of each side, split across processors,
// connected or as a soup: segments through a side, a side's middle vertex and a corner cross
// there, with the sign of their far end.
TEST(MeshFacetedSurfaceSegmentIntersection, squareLoop2d_crossingsOfSidesVerticesAndCorner)
{
  // Counterclockwise, so the normals (Facet2d's normal is crossZ(x1-x0)) point out.
  const std::vector<stk::math::Vector3d> vertices{{0, 0, 0},
      {0.5, 0, 0},
      {1, 0, 0},
      {1, 0.5, 0},
      {1, 1, 0},
      {0.5, 1, 0},
      {0, 1, 0},
      {0, 0.5, 0}};
  std::vector<std::array<unsigned, 3>> allSegments;
  for (unsigned i = 0; i < vertices.size(); ++i)
    allSegments.push_back({{i, unsigned((i + 1) % vertices.size()), 0}});
  const BoundingBox box(stk::math::Vector3d(-0.5, -0.5, 0.), stk::math::Vector3d(1.5, 1.5, 0.));

  for (auto & surface : build_split_connected_and_soup_surfaces(2, vertices, allSegments))
  {
    surface->prepare_to_compute(0., box, 0.);
    // Through the bottom side.
    expect_crossing_at_midpoint(*surface, {0.25, 0.1, 0.}, {0.25, -0.1, 0.}, 1);
    // Through the bottom side's middle vertex.
    expect_crossing_at_midpoint(*surface, {0.5, 0.1, 0.}, {0.5, -0.1, 0.}, 1);
    // Through the corner (1,1).
    expect_crossing_at_midpoint(*surface, {0.9, 0.9, 0.}, {1.1, 1.1, 0.}, 1);
  }
}

} // namespace krino
