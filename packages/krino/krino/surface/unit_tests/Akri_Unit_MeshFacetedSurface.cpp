// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_BoundingBox.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_UnitCubitJournalUtils.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetCoordinates.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace krino
{

namespace
{

constexpr double facetedTol = 1.e-12;

const BoundingBox facetedQueryBox(stk::math::Vector3d(-3, -3, -3), stk::math::Vector3d(3, 3, 3));

bool is_serial_run() { return 1 == stk::parallel_machine_size(MPI_COMM_WORLD); }

// Deterministic points in a box, the same on every processor.
std::vector<stk::math::Vector3d>
grid_points_in_box(const BoundingBox & box, const unsigned numPerDirection)
{
  std::vector<stk::math::Vector3d> points;
  const stk::math::Vector3d lo(box.get_min()[0], box.get_min()[1], box.get_min()[2]);
  const stk::math::Vector3d hi(box.get_max()[0], box.get_max()[1], box.get_max()[2]);
  for (unsigned i = 0; i < numPerDirection; ++i)
    for (unsigned j = 0; j < numPerDirection; ++j)
      for (unsigned k = 0; k < numPerDirection; ++k)
      {
        // offset by an irrational fraction so no point lies on a facet edge or vertex
        const double fx = (i + 0.3819660113) / numPerDirection;
        const double fy = (j + 0.2360679775) / numPerDirection;
        const double fz = (k + 0.6180339887) / numPerDirection;
        points.emplace_back(lo[0] + fx * (hi[0] - lo[0]),
            lo[1] + fy * (hi[1] - lo[1]),
            lo[2] + fz * (hi[2] - lo[2]));
      }
  return points;
}

// The facets with index i such that i % numProcs == rank, so every facet is owned by one processor.
std::vector<std::array<unsigned, 3>>
facets_for_this_processor(const std::vector<std::array<unsigned, 3>> & allFacets)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<std::array<unsigned, 3>> myFacets;
  for (size_t i = 0; i < allFacets.size(); ++i)
    if (int(i % numProcs) == me) myFacets.push_back(allFacets[i]);
  return myFacets;
}

// The facets as an unconnected (soup) surface, all of them on every processor: the reference that a
// surface with the facets split across processors, connected or not, is compared to.
std::unique_ptr<MeshFacetedSurface>
build_faceted_surface(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facets)
{
  FacetCoordinates facetCoords(3);
  for (auto & facet : facets)
    facetCoords.emplace_back_3d(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]]);
  return std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, facetCoords.coordinates());
}

size_t faceted_global_sum(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

} // namespace

TEST(MeshFacetedSurface, singleTriangle_signedDistanceClosestPointNormalAndSegmentIntersection)
{
  if (!is_serial_run()) return;

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, {{{0, 1, 2}}});
  surface.prepare_to_compute(0., facetedQueryBox, 0.);

  EXPECT_EQ(FACETED_SURFACE_3D, surface.type());
  EXPECT_EQ(1u, surface.num_owned_facets());
  EXPECT_EQ(0u, surface.num_ghosted_facets());

  EXPECT_NEAR(0.5, surface.point_signed_distance(stk::math::Vector3d(0.25, 0.25, 0.5)), facetedTol);
  EXPECT_NEAR(-0.5, surface.point_signed_distance(stk::math::Vector3d(0.25, 0.25, -0.5)), facetedTol);

  const stk::math::Vector3d closest = surface.closest_point(stk::math::Vector3d(0.25, 0.25, 2.));
  EXPECT_NEAR(0., (closest - stk::math::Vector3d(0.25, 0.25, 0.)).length(), facetedTol);

  const stk::math::Vector3d normal =
      surface.closest_point_normal(stk::math::Vector3d(0.25, 0.25, 2.));
  EXPECT_NEAR(0., (normal - stk::math::Vector3d(0., 0., 1.)).length(), facetedTol);

  const auto [crossingSign, location] = surface.compute_intersection_with_segment(
      stk::math::Vector3d(0.25, 0.25, -1.), stk::math::Vector3d(0.25, 0.25, 1.), 0.);
  EXPECT_EQ(1, crossingSign);
  EXPECT_NEAR(0.5, location, facetedTol);
}

TEST(MeshFacetedSurface, singleSegment_signedDistanceAndClosestPoint)
{
  if (!is_serial_run()) return;

  MeshFacetedSurface surface(MPI_COMM_WORLD, 2, {{0, 0, 0}, {2, 0, 0}}, {{{0, 1, 0}}});
  surface.prepare_to_compute(0., facetedQueryBox, 0.);

  EXPECT_EQ(FACETED_SURFACE_2D, surface.type());
  EXPECT_EQ(1u, surface.num_owned_facets());

  // Facet2d's normal is crossZ(x1-x0), which is -y for this segment.
  EXPECT_NEAR(-0.5, surface.point_signed_distance(stk::math::Vector3d(1., 0.5, 0.)), facetedTol);
  EXPECT_NEAR(0.5, surface.point_signed_distance(stk::math::Vector3d(1., -0.5, 0.)), facetedTol);

  const stk::math::Vector3d closest = surface.closest_point(stk::math::Vector3d(3., 1., 0.));
  EXPECT_NEAR(0., (closest - stk::math::Vector3d(2., 0., 0.)).length(), facetedTol);
}

TEST(MeshFacetedSurface, twoTrianglesApart_journalHasTheFacetInTheRegion)
{
  if (!is_serial_run()) return;

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3,
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {2, 0, 0}, {3, 0, 0}, {2, 1, 0}}, {{{0, 1, 2}}, {{3, 4, 5}}});
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  CubitJournal journal;
  surface.add_to_journal(journal, BoundingBox(stk::math::Vector3d(-0.5, -0.5, -0.5), stk::math::Vector3d(0.5, 0.5, 0.5)));
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  ASSERT_EQ(1u, parsed.tris.size());
  EXPECT_NEAR(0.5, total_area_of_tris(parsed), facetedTol);
}

TEST(MeshFacetedSurface, singleSegment_journalHasTheSegment)
{
  if (!is_serial_run()) return;

  MeshFacetedSurface surface(MPI_COMM_WORLD, 2, {{0, 0, 0}, {2, 0, 0}}, {{{0, 1, 0}}});
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  CubitJournal journal;
  surface.add_to_journal(journal, BoundingBox(stk::math::Vector3d(-0.5, -0.5, 0.), stk::math::Vector3d(0.5, 0.5, 0.)));
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  ASSERT_EQ(1u, parsed.edges.size());
  EXPECT_NEAR(2., total_length_of_edges(parsed), facetedTol);
}

TEST(MeshFacetedSurface, sphereSplitAcrossProcessors_sameDistancesAsSoupWithAllFacets)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  const BoundingBox sphereQueryBox(
      stk::math::Vector3d(-1.5, -1.5, -1.5), stk::math::Vector3d(1.5, 1.5, 1.5));

  // Reference: every facet on every processor, so no communication is involved.
  std::unique_ptr<MeshFacetedSurface> reference = build_faceted_surface(vertices, allFacets);
  reference->prepare_to_compute(0., sphereQueryBox, 0.);

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(0., sphereQueryBox, 0.);

  EXPECT_EQ(allFacets.size(), faceted_global_sum(surface.num_owned_facets()));
  EXPECT_EQ(allFacets.size(),
      surface.num_searched_facets()); // untruncated: every facet is searched, each once

  for (auto & x : grid_points_in_box(sphereQueryBox, 6))
  {
    EXPECT_DOUBLE_EQ(reference->point_signed_distance(x), surface.point_signed_distance(x)) << x;
    EXPECT_NEAR(x.length() - radius, surface.point_signed_distance(x), 0.02) << x;
    EXPECT_DOUBLE_EQ(0., (reference->closest_point(x) - surface.closest_point(x)).length()) << x;
    EXPECT_DOUBLE_EQ(
        0., (reference->closest_point_normal(x) - surface.closest_point_normal(x)).length())
        << x;
  }
}

TEST(MeshFacetedSurface, sphereSplitAcrossProcessors_truncatedQueriesGhostOnlyNearbyFacets)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  // Each processor queries a slab of x near its own rank, with a narrow band, so it needs only the
  // facets near that slab: some of them owned by other processors.
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const double slabWidth = 2. * radius / numProcs;
  const double xLo = -radius + me * slabWidth;
  const BoundingBox localQueryBox(
      stk::math::Vector3d(xLo, -1., -1.), stk::math::Vector3d(xLo + slabWidth, 1., 1.));
  const double narrowBand = 0.1;

  std::unique_ptr<MeshFacetedSurface> reference = build_faceted_surface(vertices, allFacets);
  reference->prepare_to_compute(0., localQueryBox, narrowBand);

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(0., localQueryBox, narrowBand);

  EXPECT_LT(surface.num_searched_facets(), allFacets.size() + 1);
  if (numProcs > 1)
  {
    EXPECT_GT(faceted_global_sum(surface.num_ghosted_facets()), 0u);
  }

  const double farField = 2. * narrowBand;
  for (auto & x : grid_points_in_box(localQueryBox, 5))
    EXPECT_DOUBLE_EQ(reference->possibly_truncated_point_signed_distance(x, narrowBand, farField),
        surface.possibly_truncated_point_signed_distance(x, narrowBand, farField))
        << x;

  // prepare_to_compute again (as when queries move): ghosts are replaced, not accumulated.
  const size_t numSearchedBefore = surface.num_searched_facets();
  surface.prepare_to_compute(0., localQueryBox, narrowBand);
  EXPECT_EQ(numSearchedBefore, surface.num_searched_facets());
}

namespace
{

// A closed prism over a counterclockwise footprint polygon, z in [z0,z1], with outward facets.  Its
// top and bottom are fans from footprint vertex fanApex, from which the whole footprint is visible.
void append_prism(const std::vector<std::array<double, 2>> & footprint,
    const unsigned fanApex,
    const double z0,
    const double z1,
    std::vector<stk::math::Vector3d> & vertices,
    std::vector<std::array<unsigned, 3>> & facets)
{
  const unsigned n = footprint.size();
  const unsigned bottom = vertices.size();
  const unsigned top = bottom + n;
  for (double z : {z0, z1})
    for (auto & xy : footprint)
      vertices.emplace_back(xy[0], xy[1], z);
  for (unsigned i = 1; i + 1 < n; ++i)
  {
    const unsigned a = (fanApex + i) % n;
    const unsigned b = (fanApex + i + 1) % n;
    facets.push_back({{top + fanApex, top + a, top + b}});
    facets.push_back({{bottom + fanApex, bottom + b, bottom + a}});
  }
  for (unsigned i = 0; i < n; ++i)
  {
    const unsigned j = (i + 1) % n;
    facets.push_back({{bottom + i, bottom + j, top + j}});
    facets.push_back({{bottom + i, top + j, top + i}});
  }
}

bool is_inside_polygon(
    const std::vector<std::array<double, 2>> & polygon, const double x, const double y)
{
  bool inside = false;
  for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
    if ((polygon[i][1] > y) != (polygon[j][1] > y) &&
        x < polygon[i][0] +
                (y - polygon[i][1]) * (polygon[j][0] - polygon[i][0]) /
                    (polygon[j][1] - polygon[i][1]))
      inside = !inside;
  return inside;
}

const std::vector<std::array<double, 2>> lFootprint{
    {0., 0.}, {2., 0.}, {2., 1.}, {1., 1.}, {1., 2.}, {0., 2.}}; // counterclockwise
constexpr unsigned lFootprintFanApex = 3;                        // the reentrant corner (1,1)

// The exact signed distance to a closed surface: the distance to its nearest facet, negative
// inside.
template <class FACET, class IS_INSIDE>
double exact_signed_distance(
    const std::vector<FACET> & facets, const IS_INSIDE & is_inside, const stk::math::Vector3d & x)
{
  double minSqrDist = std::numeric_limits<double>::max();
  for (auto & facet : facets)
    minSqrDist = std::min(minSqrDist, facet.point_distance_squared(x));
  return is_inside(x) ? -std::sqrt(minSqrDist) : std::sqrt(minSqrDist);
}

// The signed distance and closest point normal match the exact ones at every point.  Many of the
// points around a polyhedron are closest to one of its edges or vertices, where the sign comes from
// the pseudo normal of the facets connected there.
template <class FACET, class IS_INSIDE, class SURFACE>
void expect_exact_signed_distances_and_normals(const SURFACE & surface,
    const std::vector<FACET> & facets,
    const IS_INSIDE & is_inside,
    const std::vector<stk::math::Vector3d> & points)
{
  for (auto & x : points)
  {
    const double expected = exact_signed_distance(facets, is_inside, x);
    EXPECT_NEAR(expected, surface.point_signed_distance(x), facetedTol) << x;
    const stk::math::Vector3d expectedNormal =
        (expected > 0. ? 1. : -1.) * (x - surface.closest_point(x)).unit_vector();
    EXPECT_NEAR(0., (expectedNormal - surface.closest_point_normal(x)).length(), 1.e-9) << x;
  }
}

std::vector<Facet3d> build_facets(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facets)
{
  std::vector<Facet3d> result;
  for (auto & f : facets)
    result.emplace_back(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
  return result;
}

} // namespace

// Two cubes touching at a corner, with their own vertices there: the sign near the corner comes
// from the facets of each cube separately.
TEST(MeshFacetedSurface, cubesTouchingAtCorner_signFromConnectedFacets)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  append_prism({{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, 0, 0., 1., vertices, allFacets);
  append_prism({{-1., -1.}, {0., -1.}, {0., 0.}, {-1., 0.}}, 0, -1., 0., vertices, allFacets);
  const auto is_inside = [](const stk::math::Vector3d & x)
  {
    const auto in = [&](double lo, double hi)
    { return x[0] > lo && x[0] < hi && x[1] > lo && x[1] < hi && x[2] > lo && x[2] < hi; };
    return in(0., 1.) || in(-1., 0.);
  };

  const BoundingBox box(stk::math::Vector3d(-1.4, -1.4, -1.4), stk::math::Vector3d(1.4, 1.4, 1.4));
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(0., box, 0.);

  std::vector<stk::math::Vector3d> points = grid_points_in_box(box, 9);
  for (auto & x : grid_points_in_box(
           BoundingBox(stk::math::Vector3d(-0.3, -0.3, -0.3), stk::math::Vector3d(0.3, 0.3, 0.3)),
           6))
    points.push_back(x); // near the touching corner
  expect_exact_signed_distances_and_normals(
      surface, build_facets(vertices, allFacets), is_inside, points);
}

// An L-shaped prism has a concave edge, where the points inside nearest to it get their sign from
// the facets on either side of it.
TEST(MeshFacetedSurface, lShapedPrism_signFromConnectedFacetsAtConvexAndConcaveEdges)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  append_prism(lFootprint, lFootprintFanApex, 0., 1., vertices, allFacets);
  const auto is_inside = [](const stk::math::Vector3d & x)
  { return x[2] > 0. && x[2] < 1. && is_inside_polygon(lFootprint, x[0], x[1]); };

  const BoundingBox box(stk::math::Vector3d(-0.4, -0.4, -0.4), stk::math::Vector3d(2.4, 2.4, 1.4));
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(0., box, 0.);

  std::vector<stk::math::Vector3d> points = grid_points_in_box(box, 9);
  for (auto & x : grid_points_in_box(
           BoundingBox(stk::math::Vector3d(0.6, 0.6, 0.1), stk::math::Vector3d(1., 1., 0.9)), 5))
    points.push_back(x); // inside, nearest the concave edge
  expect_exact_signed_distances_and_normals(
      surface, build_facets(vertices, allFacets), is_inside, points);
}

TEST(MeshFacetedSurface, lShapedPolygon2d_signFromConnectedSegmentsAtVertices)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allSegments;
  for (unsigned i = 0; i < lFootprint.size(); ++i)
  {
    vertices.emplace_back(lFootprint[i][0], lFootprint[i][1], 0.);
    allSegments.push_back({{i,
        unsigned((i + 1) % lFootprint.size()),
        0}}); // outward normals (Facet2d's normal is crossZ(x1-x0))
  }
  std::vector<Facet2d> segments;
  for (auto & s : allSegments)
    segments.emplace_back(vertices[s[0]], vertices[s[1]]);
  const auto is_inside = [](const stk::math::Vector3d & x)
  { return is_inside_polygon(lFootprint, x[0], x[1]); };

  const BoundingBox box(stk::math::Vector3d(-0.4, -0.4, 0.), stk::math::Vector3d(2.4, 2.4, 0.));
  MeshFacetedSurface surface(MPI_COMM_WORLD, 2, vertices, facets_for_this_processor(allSegments));
  surface.prepare_to_compute(0., box, 0.);

  std::vector<stk::math::Vector3d> points;
  for (unsigned i = 0; i < 25; ++i)
    for (unsigned j = 0; j < 25; ++j)
      points.emplace_back(
          -0.4 + 2.8 * (i + 0.3819660113) / 25, -0.4 + 2.8 * (j + 0.2360679775) / 25, 0.);
  expect_exact_signed_distances_and_normals(surface, segments, is_inside, points);
}

namespace
{

// A tetrahedron with a sharp corner at the origin, and a cube [-1,-gap]^3 whose corner is gap from
// it along -(1,1,1): two separate bodies with their own vertices.
void fill_sharp_tet_and_cube_separated_at_corners(const double gap,
    std::vector<stk::math::Vector3d> & vertices,
    std::vector<std::array<unsigned, 3>> & facets)
{
  vertices = {{0., 0., 0.}, {1., 0.2, 0.2}, {0.2, 1., 0.2}, {0.2, 0.2, 1.}};
  facets = {{{0, 2, 1}}, {{0, 3, 2}}, {{0, 1, 3}}, {{1, 2, 3}}};
  append_prism(
      {{-1., -1.}, {-gap, -1.}, {-gap, -gap}, {-1., -gap}}, 0, -1., -gap, vertices, facets);
}

} // namespace

// At a point in the gap between the two corners, nearest the tetrahedron's corner, the sign comes
// from the tetrahedron's facets alone.  An unconnected surface, which takes the facets whose
// closest points are within a tolerance of the nearest one, would include the cube's and get the
// sign wrong.
TEST(MeshFacetedSurface, bodiesSeparatedByTinyGap_signFromFacetsConnectedToClosestVertex)
{
  const double gap = 1.e-7;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sharp_tet_and_cube_separated_at_corners(gap, vertices, allFacets);

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(
      0., BoundingBox(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.)), 0.);

  const stk::math::Vector3d x = -0.4 * gap * stk::math::Vector3d(1., 1., 1.);
  EXPECT_NEAR(x.length(), surface.point_signed_distance(x), 1.e-6 * x.length());
  EXPECT_NEAR(0., (x.unit_vector() - surface.closest_point_normal(x)).length(), 1.e-9);
}

// Not taken to be connected, the surface gets its sign from the facets near the closest point, as
// the soup surface does, and so it gets the same (here wrong) sign in the tiny gap.
TEST(MeshFacetedSurface, bodiesSeparatedByTinyGap_notConnected_signFromNearbyFacetsAsSoup)
{
  const double gap = 1.e-7;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sharp_tet_and_cube_separated_at_corners(gap, vertices, allFacets);
  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));

  std::unique_ptr<MeshFacetedSurface> reference = build_faceted_surface(vertices, allFacets);
  reference->prepare_to_compute(0., box, 0.);

  MeshFacetedSurface surface(
      MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets), false);
  surface.prepare_to_compute(0., box, 0.);

  const stk::math::Vector3d x = -0.4 * gap * stk::math::Vector3d(1., 1., 1.);
  EXPECT_LT(reference->point_signed_distance(x), 0.);
  EXPECT_DOUBLE_EQ(reference->point_signed_distance(x), surface.point_signed_distance(x));
  // The facet order, so the last bits of the summed normal, depend on the decomposition.
  EXPECT_NEAR(
      0., (reference->closest_point_normal(x) - surface.closest_point_normal(x)).length(), facetedTol);
}

// Each processor queries its own slab across the concave edge of the L-shaped prism, with a narrow
// band, so the facets sharing an edge or vertex near its points come partly from other processors.
TEST(MeshFacetedSurface, lShapedPrismQueriedInSlabs_signFromConnectedFacetsIncludingGhosts)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  append_prism(lFootprint, lFootprintFanApex, 0., 1., vertices, allFacets);
  const auto is_inside = [](const stk::math::Vector3d & x)
  { return x[2] > 0. && x[2] < 1. && is_inside_polygon(lFootprint, x[0], x[1]); };

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const double zLo = -0.2 + me * 1.4 / numProcs;
  const BoundingBox mySlab(
      stk::math::Vector3d(0.6, 0.6, zLo), stk::math::Vector3d(1.4, 1.4, zLo + 1.4 / numProcs));
  const double narrowBand = 0.5;

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(allFacets));
  surface.prepare_to_compute(0., mySlab, narrowBand);
  if (numProcs > 1)
  {
    EXPECT_GT(faceted_global_sum(surface.num_ghosted_facets()), 0u);
  }

  const std::vector<Facet3d> facets = build_facets(vertices, allFacets);
  for (auto & x : grid_points_in_box(mySlab, 6))
  {
    const double expected = exact_signed_distance(facets, is_inside, x);
    const double expectedTruncated = (std::abs(expected) > narrowBand) ? 2. * narrowBand : expected;
    EXPECT_NEAR(expectedTruncated,
        surface.possibly_truncated_point_signed_distance(x, narrowBand, 2. * narrowBand),
        facetedTol)
        << x;
  }
}

// As sphereSplitAcrossProcessors_sameDistancesAsSoupWithAllFacets, with the facets given
// by their vertices' coordinates (no shared vertices), as contouring gives them, and split across
// processors: the distances, closest points and normals are those with all the facets on every
// processor.
TEST(MeshFacetedSurface,
    unconnectedSphereFromCoordinatesSplitAcrossProcessors_sameAsAllFacetsOnEveryProcessor)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  const BoundingBox sphereQueryBox(
      stk::math::Vector3d(-1.5, -1.5, -1.5), stk::math::Vector3d(1.5, 1.5, 1.5));

  std::unique_ptr<MeshFacetedSurface> reference = build_faceted_surface(vertices, allFacets);
  reference->prepare_to_compute(0., sphereQueryBox, 0.);

  std::vector<stk::math::Vector3d> facetVertexCoords;
  for (auto & facet : facets_for_this_processor(allFacets))
    for (unsigned n = 0; n < 3; ++n)
      facetVertexCoords.push_back(vertices[facet[n]]);
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, facetVertexCoords);
  surface.prepare_to_compute(0., sphereQueryBox, 0.);

  EXPECT_EQ(allFacets.size(), faceted_global_sum(surface.num_owned_facets()));
  EXPECT_FALSE(surface.facet_mesh().has_shared_vertices());

  for (auto & x : grid_points_in_box(sphereQueryBox, 6))
  {
    EXPECT_DOUBLE_EQ(reference->point_signed_distance(x), surface.point_signed_distance(x)) << x;
    EXPECT_DOUBLE_EQ(0., (reference->closest_point(x) - surface.closest_point(x)).length()) << x;
    EXPECT_DOUBLE_EQ(
        0., (reference->closest_point_normal(x) - surface.closest_point_normal(x)).length())
        << x;
  }
}

// A connected sphere split across processors, with a velocity at each vertex that is linear in its
// position: interpolated over the closest facet, the velocity at the closest point is that field
// there.
TEST(MeshFacetedSurface, velocityAtClosestPoint_linearVelocityInterpolatedExactly)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, allFacets);
  auto velocity_at = [](const stk::math::Vector3d & x)
  { return stk::math::Vector3d(x[1] + 2. * x[2], 0.5 - x[0], 3. * x[2] - x[0]); };

  std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> vertexIdsAndCoords;
  std::vector<stk::math::Vector3d> vertexVelocities;
  for (size_t v = 0; v < vertices.size(); ++v)
  {
    vertexIdsAndCoords.emplace_back(v + 1, vertices[v]);
    vertexVelocities.push_back(velocity_at(vertices[v]));
  }
  std::vector<std::array<FacetVertexId, 3>> facetVertexIds;
  for (auto & facet : facets_for_this_processor(allFacets))
    facetVertexIds.push_back({{facet[0] + 1u, facet[1] + 1u, facet[2] + 1u}});

  const BoundingBox facetedQueryBox(
      stk::math::Vector3d(-1.5, -1.5, -1.5), stk::math::Vector3d(1.5, 1.5, 1.5));
  MeshFacetedSurface surface(
      MPI_COMM_WORLD, 3, vertexIdsAndCoords, facetVertexIds, vertexVelocities);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  EXPECT_TRUE(surface.facet_mesh().has_vertex_velocities());

  for (auto & x : grid_points_in_box(facetedQueryBox, 6))
    EXPECT_NEAR(0.,
        (velocity_at(surface.closest_point(x)) - surface.velocity_at_closest_point(x)).length(),
        1.e-12)
        << x;

  MeshFacetedSurface withoutVelocities(MPI_COMM_WORLD, 3, vertexIdsAndCoords, facetVertexIds);
  withoutVelocities.prepare_to_compute(0., facetedQueryBox, 0.);
  EXPECT_ANY_THROW(withoutVelocities.velocity_at_closest_point(stk::math::Vector3d::ZERO));
}

// set_facets replaces the facets at the next prepare_to_compute, switching between a soup (no
// shared vertices) and connected facets: a square of 2 triangles at z=0 (connected), then a
// triangle at z=1 (soup), then the square again.  Processor 0 owns them all.
TEST(MeshFacetedSurface, setFacets_replacesTheFacetsAtTheNextPrepareAsSoupOrConnected)
{
  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  FacetMesh::VertexIdsAndCoords square;
  FacetMesh::FacetVertexIds squareFacets;
  FacetCoordinates triangle(3);
  if (isProc0)
  {
    square = {{1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {0., 1., 0.}}, {4, {1., 1., 0.}}};
    squareFacets = {{{1, 2, 3}}, {{2, 4, 3}}};
    triangle.emplace_back_3d({0., 0., 1.}, {1., 0., 1.}, {0., 1., 1.});
  }
  const stk::math::Vector3d x(0.25, 0.25, 0.2);

  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, square, squareFacets);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  EXPECT_DOUBLE_EQ(0.2, std::abs(surface.point_signed_distance(x)));

  surface.set_facets(triangle);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  EXPECT_FALSE(surface.facet_mesh().has_shared_vertices());
  EXPECT_EQ(1u, faceted_global_sum(surface.num_owned_facets()));
  EXPECT_DOUBLE_EQ(0.8, std::abs(surface.point_signed_distance(x)));
  EXPECT_EQ(triangle.coordinates(), surface.owned_facet_coordinates().coordinates());

  surface.set_facets(square, squareFacets);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);
  EXPECT_TRUE(surface.facet_mesh().has_shared_vertices());
  EXPECT_EQ(4u, surface.global_num_vertices());
  EXPECT_DOUBLE_EQ(0.2, std::abs(surface.point_signed_distance(x)));
}

// With no facets on any processor (a level set with no interface), a soup surface gives the
// far-field value, with and without a narrow band.
TEST(MeshFacetedSurface, emptySoup_farFieldDistances)
{
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, std::vector<stk::math::Vector3d>{});
  const double farField = 0.5;
  for (const double narrowBand : {0., 0.5})
  {
    surface.prepare_to_compute(0., facetedQueryBox, narrowBand);
    for (auto & x : {stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., -2., 0.5)})
    {
      EXPECT_EQ(
          farField, surface.possibly_truncated_point_signed_distance(x, narrowBand, farField));
      EXPECT_EQ(farField, surface.point_unsigned_distance(x, narrowBand, farField));
    }
  }
}

namespace
{

// The facet vertex coordinates of these facets.
std::vector<stk::math::Vector3d>
faceted_vertex_coordinates(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facets)
{
  std::vector<stk::math::Vector3d> coords;
  for (auto & facet : facets)
    for (unsigned n = 0; n < 3; ++n)
      coords.push_back(vertices[facet[n]]);
  return coords;
}

} // namespace

// A soup with a duplicate of one triangle and an opposite of another, each on a different processor
// from its partner when there are several: set_remove_coincident_facets removes them all over the
// processors, leaving the surface of the one triangle.
TEST(MeshFacetedSurface, splitSoupWithRemoveCoincidentFacets_coincidentFacetsOnAnyProcessorsRemoved)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  const std::vector<std::array<unsigned, 3>> allFacets{{{0, 1, 2}},
      {{1, 2, 0}}, // the same orientation as the first: one stays
      {{1, 3, 2}},
      {{1, 2, 3}}}; // opposite orientations: both go

  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  MeshFacetedSurface expected(MPI_COMM_WORLD,
      3,
      faceted_vertex_coordinates(vertices,
          isProc0 ? std::vector<std::array<unsigned, 3>>{allFacets[0]} :
                    std::vector<std::array<unsigned, 3>>{}));
  expected.prepare_to_compute(0., facetedQueryBox, 0.);

  MeshFacetedSurface surface(
      MPI_COMM_WORLD, 3, faceted_vertex_coordinates(vertices, facets_for_this_processor(allFacets)));
  surface.set_remove_coincident_facets(true);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);

  EXPECT_EQ(1u, faceted_global_sum(surface.num_owned_facets()));
  for (auto & x : grid_points_in_box(
           BoundingBox(stk::math::Vector3d(-0.5, -0.5, -0.5), stk::math::Vector3d(1.5, 1.5, 0.5)),
           5))
    EXPECT_DOUBLE_EQ(expected.point_signed_distance(x), surface.point_signed_distance(x)) << x;
}

// A tetrahedron's surface with its bottom face missing, split across processors, connected or as a
// soup: set_check_watertight warns (on processor 0) of the open edges around the missing face,
// unless the watertight bounding box leaves them out or puts them on one of its faces.  The closed
// surface gives no warning.
TEST(MeshFacetedSurface, splitSurfacesWithCheckWatertight_warnOnlyOfOpenEdgesInTheBox)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  const std::vector<std::array<unsigned, 3>> openFacets{{{0, 1, 3}}, {{1, 2, 3}}, {{2, 0, 3}}};
  std::vector<std::array<unsigned, 3>> closedFacets = openFacets;
  closedFacets.push_back({{0, 2, 1}});

  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  const BoundingBox containingBox(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(2., 2., 2.));
  const BoundingBox boxAbove(stk::math::Vector3d(-1., -1., 0.5), stk::math::Vector3d(2., 2., 2.));
  const BoundingBox boxWithOpenEdgesOnItsFace(
      stk::math::Vector3d(-1., -1., 0.), stk::math::Vector3d(2., 2., 2.));
  const std::vector<
      std::tuple<const std::vector<std::array<unsigned, 3>> *, std::optional<BoundingBox>, bool>>
      cases{{&openFacets, std::nullopt, true},
          {&openFacets, containingBox, true},
          {&openFacets, boxAbove, false},
          {&openFacets, boxWithOpenEdgesOnItsFace, false},
          {&closedFacets, std::nullopt, false}};
  for (auto & [facets, watertightBox, expectWarning] : cases)
  {
    const std::vector<std::array<unsigned, 3>> myFacets = facets_for_this_processor(*facets);
    for (const bool connected : {true, false})
    {
      std::unique_ptr<MeshFacetedSurface> surface = connected ?
          std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, vertices, myFacets) :
          std::make_unique<MeshFacetedSurface>(
              MPI_COMM_WORLD, 3, faceted_vertex_coordinates(vertices, myFacets));
      surface->set_check_watertight(true);
      if (watertightBox) surface->set_watertight_bounding_box(*watertightBox);
      const unsigned numWarningsBefore = stk::get_warning_count();
      surface->prepare_to_compute(0., facetedQueryBox, 0.);
      if (isProc0)
      {
        EXPECT_EQ(expectWarning, stk::get_warning_count() > numWarningsBefore)
            << (connected ? "connected" : "soup") << ", " << facets->size() << " facets, "
            << (watertightBox ? "with" : "without") << " box";
      }
    }
  }
}

namespace
{

std::string expected_print_sizes(const size_t numFacets,
    const size_t localNumFacets,
    const double minArea,
    const double maxArea,
    const double totalArea)
{
  std::ostringstream out;
  out << "\t Global sizes: { " << numFacets << " facets on }\n"
      << "\t Local sizes: { " << localNumFacets << " facets on }\n"
      << "\t Global areas: { Min = " << minArea << ", Max = " << maxArea
      << ", Total = " << totalArea << " }\n";
  return out.str();
}

} // namespace

// Two triangles of areas 0.5 and 2 (segments of lengths 1 and 3 in 2D) split across processors,
// connected or as a soup: print_sizes gives the global number of facets and their areas, and this
// processor's number of facets.
TEST(MeshFacetedSurface, printSizes_globalAndLocalFacetCountsAndGlobalAreas)
{
  const std::vector<std::array<unsigned, 3>> twoFacets{{{0, 1, 2}}, {{3, 4, 5}}};
  const std::vector<std::array<unsigned, 3>> myFacets = facets_for_this_processor(twoFacets);
  const size_t localNumFacets = myFacets.size();

  const std::vector<stk::math::Vector3d> triangleVertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}, {2., 0., 1.}, {0., 2., 1.}};
  const std::vector<stk::math::Vector3d> segmentVertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 0., 0.}, {0., 1., 0.}, {3., 1., 0.}, {0., 0., 0.}};
  for (const unsigned dim : {3u, 2u})
  {
    const auto & vertices = (dim == 3) ? triangleVertices : segmentVertices;
    std::vector<stk::math::Vector3d> soupCoords;
    for (auto & facet : myFacets)
      for (unsigned n = 0; n < dim; ++n)
        soupCoords.push_back(vertices[facet[n]]);
    MeshFacetedSurface connected(MPI_COMM_WORLD, dim, vertices, myFacets);
    MeshFacetedSurface soup(MPI_COMM_WORLD, dim, soupCoords);
    const std::string expected = (dim == 3) ?
        expected_print_sizes(2, localNumFacets, 0.5, 2., 2.5) :
        expected_print_sizes(2, localNumFacets, 1., 3., 4.);
    for (MeshFacetedSurface * surface : {&connected, &soup})
    {
      surface->prepare_to_compute(0., facetedQueryBox, 0.);
      EXPECT_EQ(expected, surface->print_sizes()) << "dimension " << dim;
    }
  }
}

TEST(MeshFacetedSurface, setFacetsOfOtherDimension_throws)
{
  MeshFacetedSurface surface(MPI_COMM_WORLD, 3, std::vector<stk::math::Vector3d>{});
  EXPECT_ANY_THROW(surface.set_facets(FacetCoordinates(2)));
}

// set_facets replaces the facets with ones without vertex velocities, connected or as a soup.
TEST(MeshFacetedSurface, setFacets_dropsTheVertexVelocities)
{
  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  FacetMesh::VertexIdsAndCoords triangle;
  FacetMesh::FacetVertexIds triangleFacets;
  std::vector<stk::math::Vector3d> velocities;
  FacetCoordinates soupTriangle(3);
  if (isProc0)
  {
    triangle = {{1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {0., 1., 0.}}};
    triangleFacets = {{{1, 2, 3}}};
    velocities.assign(3, stk::math::Vector3d(1., 2., 3.));
    soupTriangle.emplace_back_3d({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.});
  }
  const stk::math::Vector3d x(0.25, 0.25, 0.5);

  for (const bool setSoup : {false, true})
  {
    MeshFacetedSurface surface(MPI_COMM_WORLD, 3, triangle, triangleFacets, velocities);
    surface.prepare_to_compute(0., facetedQueryBox, 0.);
    EXPECT_DOUBLE_EQ(
        0., (stk::math::Vector3d(1., 2., 3.) - surface.velocity_at_closest_point(x)).length());

    if (setSoup)
      surface.set_facets(soupTriangle);
    else
      surface.set_facets(triangle, triangleFacets);
    surface.prepare_to_compute(0., facetedQueryBox, 0.);
    EXPECT_DOUBLE_EQ(0.5, surface.point_signed_distance(x));
    EXPECT_ANY_THROW(surface.velocity_at_closest_point(x)) << (setSoup ? "soup" : "connected");
  }
}

// A connected square of two triangles, with a degenerate (collinear) triangle along one of its
// edges through a vertex at the edge's middle: near that edge and its end vertices, the signed
// distances and normals are those without the degenerate triangle.  Only degenerate facets near a
// point throw.
TEST(MeshFacetedSurface, connectedWithDegenerateFacet_signAndNormalIgnoreIt)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {0.5, 0., 0.}};
  const std::vector<std::array<unsigned, 3>> square{{{0, 1, 3}}, {{0, 3, 2}}};
  std::vector<std::array<unsigned, 3>> squareAndDegenerate = square;
  squareAndDegenerate.push_back({{0, 4, 1}});

  MeshFacetedSurface expected(MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(square));
  MeshFacetedSurface surface(
      MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(squareAndDegenerate));
  expected.prepare_to_compute(0., facetedQueryBox, 0.);
  surface.prepare_to_compute(0., facetedQueryBox, 0.);

  for (const double z : {0.1, -0.1})
    for (auto & x : {stk::math::Vector3d(0.3, -0.2, z),
             stk::math::Vector3d(0.5, -0.2, z),
             stk::math::Vector3d(-0.2, -0.2, z),
             stk::math::Vector3d(1.2, -0.3, z)})
    {
      EXPECT_DOUBLE_EQ(expected.point_signed_distance(x), surface.point_signed_distance(x)) << x;
      EXPECT_DOUBLE_EQ(
          0., (expected.closest_point_normal(x) - surface.closest_point_normal(x)).length())
          << x;
    }

  const std::vector<std::array<unsigned, 3>> onlyDegenerate{{{0, 4, 1}}};
  MeshFacetedSurface degenerate(
      MPI_COMM_WORLD, 3, vertices, facets_for_this_processor(onlyDegenerate));
  degenerate.prepare_to_compute(0., facetedQueryBox, 0.);
  const stk::math::Vector3d x(0.3, -0.2, 0.1);
  EXPECT_ANY_THROW(degenerate.point_signed_distance(x));
  EXPECT_ANY_THROW(degenerate.closest_point_normal(x));
}

} // namespace krino
