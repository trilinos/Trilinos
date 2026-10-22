// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Faceted_Surface.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <krino/mesh_utils/Akri_AllReduce.hpp>

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_math/StkVector.hpp>

#include <algorithm>
#include <memory>
#include <optional>

namespace krino {

// Serial coverage for FacetedSurfaceBase / Faceted_Surface built directly from
// facets (no mesh, no STL).  We exercise the query paths that
// prepare_to_compute wires up: signed / unsigned point_distance (including
// narrow-band truncation), closest_point, closest_point_normal, the empty
// surface, the facet-container accessors, swap, and the segment-intersection
// result.  All gold values are analytic for a single reference facet:
//   3d: triangle (0,0,0),(1,0,0),(0,1,0) in z=0, facet normal +z.
//   2d: segment (0,0,0)->(2,0,0), facet normal crossZ(x1-x0)=(0,-1,0).

namespace {
const double tol = 1.e-9;

// A generous query box that comfortably contains both the reference facets and
// every query point; the padding computed in prepare_to_compute then pulls the
// facets into the search tree.
const BoundingBox queryBox(stk::math::Vector3d(-3,-3,-3), stk::math::Vector3d(3,3,3));

bool is_serial() { return 1 == stk::parallel_machine_size(MPI_COMM_WORLD); }

std::unique_ptr<FacetedSurfaceBase> build_single_triangle_surface()
{
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  surface->prepare_to_compute(queryBox, 0.);
  return surface;
}

std::unique_ptr<FacetedSurfaceBase> build_single_segment_surface()
{
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(2);
  surface->emplace_back_2d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(2,0,0));
  surface->prepare_to_compute(queryBox, 0.);
  return surface;
}
}

//--------------------------------------------------------------------------
// facet-container accessors and type
//--------------------------------------------------------------------------

TEST(FacetedSurface, build3d_reportsTypeAndFacetCount)
{
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  EXPECT_EQ(FACETED_SURFACE_3D, surface->type());

  surface->emplace_back_3d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  EXPECT_EQ(1u, surface->size());
  EXPECT_EQ(1u, surface->get_facets_3d().size());
  EXPECT_GT(surface->storage_size(), 0u);
}

TEST(FacetedSurface, build2d_reportsTypeAndFacetCount)
{
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(2);
  EXPECT_EQ(FACETED_SURFACE_2D, surface->type());

  surface->emplace_back_2d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(2,0,0));
  EXPECT_EQ(1u, surface->size());
  EXPECT_EQ(1u, surface->get_facets_2d().size());
}

TEST(FacetedSurface, clear_emptiesLocalFacets)
{
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  ASSERT_EQ(1u, surface->size());

  surface->clear();
  EXPECT_EQ(0u, surface->size());
}

TEST(FacetedSurface, swap_exchangesFacets)
{
  std::unique_ptr<FacetedSurfaceBase> a = FacetedSurfaceBase::build(3);
  a->emplace_back_3d(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  a->emplace_back_3d(stk::math::Vector3d(0,0,1), stk::math::Vector3d(1,0,1), stk::math::Vector3d(0,1,1));

  std::unique_ptr<FacetedSurfaceBase> b = FacetedSurfaceBase::build(3);

  a->swap(*b);
  EXPECT_EQ(0u, a->size());
  EXPECT_EQ(2u, b->size());
}

//--------------------------------------------------------------------------
// 3d signed / unsigned distance
//--------------------------------------------------------------------------

TEST(FacetedSurface, signedDistance3d_positiveOnNormalSideNegativeOnOther)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  // Facet normal is +z; the projection of these query points lands on the
  // triangle interior, so the distance is exactly the z gap.
  EXPECT_NEAR(0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,0.5)), tol);
  EXPECT_NEAR(-0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,-0.5)), tol);
  EXPECT_NEAR(0., surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,0.)), tol);
}

TEST(FacetedSurface, unsignedDistance3d_isMagnitudeAndTruncatesToFarField)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  const double farField = 99.;
  // Inside the narrow band: true (unsigned) distance is returned.
  EXPECT_NEAR(0.1, surface->point_unsigned_distance(stk::math::Vector3d(0.25,0.25,0.1), 0.25, farField), tol);
  // Outside the narrow band (0.5 > 0.25): far-field value is returned instead.
  EXPECT_NEAR(farField, surface->point_unsigned_distance(stk::math::Vector3d(0.25,0.25,0.5), 0.25, farField), tol);
}

TEST(FacetedSurface, closestPoint3d_projectsOntoFacetInterior)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  expect_eq(stk::math::Vector3d(0.25,0.25,0.), surface->closest_point(stk::math::Vector3d(0.25,0.25,0.5)));
  expect_eq(stk::math::Vector3d(0.25,0.25,0.), surface->closest_point(stk::math::Vector3d(0.25,0.25,-0.7)));
}

TEST(FacetedSurface, closestPointNormal3d_isFacetNormalForInteriorProjection)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  EXPECT_TRUE(surface->can_approximate_closest_point_normal());
  expect_eq(stk::math::Vector3d(0,0,1), surface->closest_point_normal(stk::math::Vector3d(0.25,0.25,0.5)));
  expect_eq(stk::math::Vector3d(0,0,1), surface->closest_point_normal(stk::math::Vector3d(0.25,0.25,-0.5)));
}

//--------------------------------------------------------------------------
// 2d signed distance / closest point (crossZ normal convention)
//--------------------------------------------------------------------------

TEST(FacetedSurface, signedDistance2d_followsCrossZNormalConvention)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_segment_surface();

  // Facet normal is crossZ(x1-x0) = (0,-1,0), so the -y side is positive.
  EXPECT_NEAR(0.5, surface->point_signed_distance(stk::math::Vector3d(1,-0.5,0)), tol);
  EXPECT_NEAR(-0.5, surface->point_signed_distance(stk::math::Vector3d(1,0.5,0)), tol);
}

TEST(FacetedSurface, closestPoint2d_projectsOntoSegment)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_segment_surface();

  expect_eq(stk::math::Vector3d(1,0,0), surface->closest_point(stk::math::Vector3d(1,0.5,0)));
}

//--------------------------------------------------------------------------
// empty surface
//--------------------------------------------------------------------------

TEST(FacetedSurface, emptySurface_reportsEmptyAndReturnsFarFieldAndZeroClosestPoint)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->prepare_to_compute(queryBox, 0.);

  EXPECT_TRUE(surface->empty());
  EXPECT_EQ(0u, surface->size());

  const double farField = 42.;
  EXPECT_NEAR(farField, surface->point_unsigned_distance(stk::math::Vector3d(1,2,3), 0., farField), tol);
  expect_eq(stk::math::Vector3d::ZERO, surface->closest_point(stk::math::Vector3d(1,2,3)));
}

//--------------------------------------------------------------------------
// segment intersection
//--------------------------------------------------------------------------

TEST(FacetedSurface, computeIntersectionWithSegment_crossingReturnsSignAndLocation)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  // Segment pierces the z=0 triangle at its midpoint (loc 0.5); the far end is
  // on the +z (positive-distance) side, so the crossing sign is +1.
  const auto [sign, loc] = surface->compute_intersection_with_segment(
      stk::math::Vector3d(0.25,0.25,-1.), stk::math::Vector3d(0.25,0.25,1.), 1.e-6);
  EXPECT_EQ(1, sign);
  EXPECT_NEAR(0.5, loc, tol);
}

TEST(FacetedSurface, computeIntersectionWithSegment_noCrossingReturnsNoIntersection)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = build_single_triangle_surface();

  // Segment stays far to the side of the triangle: no crossing.
  const auto [sign, loc] = surface->compute_intersection_with_segment(
      stk::math::Vector3d(5,5,-1.), stk::math::Vector3d(5,5,1.), 1.e-6);
  EXPECT_EQ(0, sign);
  EXPECT_LT(loc, 0.);
}

//--------------------------------------------------------------------------
// coincident opposite-normal facet removal
//--------------------------------------------------------------------------

namespace {
// A triangle and its reverse-wound twin (same three vertices, opposite winding),
// so their normals are +z and -z and cancel in the average-normal calculation.
const stk::math::Vector3d triV0(0,0,0), triV1(1,0,0), triV2(0,1,0);
}

TEST(FacetedSurface, removeCoincident_deletesOpposingPair)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(triV0, triV1, triV2); // normal +z
  surface->emplace_back_3d(triV0, triV2, triV1); // reverse winding, normal -z
  surface->set_remove_coincident_facets(true);
  surface->prepare_to_compute(queryBox, 0.);

  // Balanced orientations cancel entirely.
  EXPECT_EQ(0u, surface->size());
  EXPECT_TRUE(surface->empty());
}

TEST(FacetedSurface, removeCoincident_offByDefaultKeepsPair)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(triV0, triV1, triV2);
  surface->emplace_back_3d(triV0, triV2, triV1);
  // Flag left at its default (off).
  surface->prepare_to_compute(queryBox, 0.);

  EXPECT_EQ(2u, surface->size());
}

TEST(FacetedSurface, removeCoincident_collapsesSameNormalDuplicateToOne)
{
  if (!is_serial()) return;
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(triV0, triV1, triV2); // normal +z
  surface->emplace_back_3d(triV0, triV1, triV2); // identical winding, normal +z
  surface->set_remove_coincident_facets(true);
  surface->prepare_to_compute(queryBox, 0.);

  // Same-normal duplicates now collapse to a single representative.
  EXPECT_EQ(1u, surface->size());

  // The surviving facet still yields a well-defined signed distance.
  EXPECT_NEAR(0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,0.5)), tol);
  EXPECT_NEAR(-0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,-0.5)), tol);
}

TEST(FacetedSurface, removeCoincident_signedDistanceIsWellDefinedAfterRemoval)
{
  if (!is_serial()) return;
  // A single triangle plus its opposing twin: with the twin removed, the query
  // resolves to the surviving +z facet instead of the indeterminate zero normal.
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(triV0, triV1, triV2); // normal +z
  surface->emplace_back_3d(triV0, triV2, triV1); // opposing twin
  surface->emplace_back_3d(triV0, triV1, triV2); // majority +z side -> one representative survives
  surface->set_remove_coincident_facets(true);
  surface->prepare_to_compute(queryBox, 0.);

  ASSERT_EQ(1u, surface->size());
  EXPECT_NEAR(0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,0.5)), tol);
  EXPECT_NEAR(-0.5, surface->point_signed_distance(stk::math::Vector3d(0.25,0.25,-0.5)), tol);
}

TEST(FacetedSurface, removeCoincidentFreeFunction_collapsesSameNormalDuplicates)
{
  // Three identically-wound copies collapse to one; count removed is 2.
  std::vector<Facet3d> facets;
  facets.emplace_back(triV0, triV1, triV2);
  facets.emplace_back(triV0, triV1, triV2);
  facets.emplace_back(triV0, triV1, triV2);
  const size_t numRemoved = remove_coincident_facets(facets);
  EXPECT_EQ(2u, numRemoved);
  ASSERT_EQ(1u, facets.size());
  // Surviving copy keeps the original +z orientation.
  EXPECT_GT(Dot(facets[0].facet_normal(), stk::math::Vector3d(0,0,1)), 0.);
}

TEST(FacetedSurface, removeCoincidentFreeFunction_unbalancedMixKeepsMajoritySide)
{
  // Two +z copies and one -z copy: one opposing pair cancels, leaving one +z copy.
  std::vector<Facet3d> facets;
  facets.emplace_back(triV0, triV1, triV2); // +z
  facets.emplace_back(triV0, triV1, triV2); // +z
  facets.emplace_back(triV0, triV2, triV1); // -z
  const size_t numRemoved = remove_coincident_facets(facets);
  EXPECT_EQ(2u, numRemoved);
  ASSERT_EQ(1u, facets.size());
  EXPECT_GT(Dot(facets[0].facet_normal(), stk::math::Vector3d(0,0,1)), 0.);
}

TEST(FacetedSurface, removeCoincidentFreeFunction_toleranceCatchesNearCoincident)
{
  // Twin perturbed by far less than relTol*facetSize (~1e-6) is still treated as
  // coincident and removed; a twin displaced well beyond tolerance is kept.
  {
    std::vector<Facet3d> facets;
    facets.emplace_back(triV0, triV1, triV2);
    facets.emplace_back(triV0 + stk::math::Vector3d(1.e-10,1.e-10,1.e-10), triV2, triV1);
    const size_t numRemoved = remove_coincident_facets(facets);
    EXPECT_EQ(2u, numRemoved);
    EXPECT_TRUE(facets.empty());
  }
  {
    std::vector<Facet3d> facets;
    facets.emplace_back(triV0, triV1, triV2);
    facets.emplace_back(triV0 + stk::math::Vector3d(0.1,0.1,0.1), triV2, triV1);
    const size_t numRemoved = remove_coincident_facets(facets);
    EXPECT_EQ(0u, numRemoved);
    EXPECT_EQ(2u, facets.size());
  }
}

//--------------------------------------------------------------------------
// watertight check
//--------------------------------------------------------------------------

namespace {
// A closed tetrahedron surface: 4 triangular faces over vertices a,b,c,d.  Every one of the 6
// edges is shared by exactly two faces, so the surface is watertight.  face(i) omits vertex i,
// so dropping any single face exposes exactly 3 boundary edges.
const stk::math::Vector3d tetA(0,0,0), tetB(1,0,0), tetC(0,1,0), tetD(0,0,1);

std::vector<Facet3d> build_tetra_faces()
{
  std::vector<Facet3d> facets;
  facets.emplace_back(tetB, tetC, tetD); // face opposite A
  facets.emplace_back(tetA, tetC, tetD); // face opposite B
  facets.emplace_back(tetA, tetB, tetD); // face opposite C
  facets.emplace_back(tetA, tetB, tetC); // face opposite D
  return facets;
}

double serial_edge_cell_size(const std::vector<Facet3d> & facets)
{
  return compute_global_edge_quantization_cell_size(MPI_COMM_WORLD, facets);
}
}

TEST(FacetedSurface, watertight_closedTetraShell_reportsNoBadEdges)
{
  if (!is_serial()) return;
  const std::vector<Facet3d> facets = build_tetra_faces();
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), std::nullopt);
  EXPECT_TRUE(badEdges.empty());
}

TEST(FacetedSurface, watertight_openShell_reportsBoundaryEdges)
{
  if (!is_serial()) return;
  std::vector<Facet3d> facets = build_tetra_faces();
  facets.pop_back(); // drop face opposite D -> its 3 edges (A-B, B-C, C-A) become boundary edges
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), std::nullopt);
  EXPECT_EQ(3u, badEdges.size());
}

TEST(FacetedSurface, watertight_nonManifoldEdge_flagged)
{
  if (!is_serial()) return;
  // Three triangles sharing the single edge A-B: that edge is incident to 3 facets (non-manifold).
  std::vector<Facet3d> facets;
  facets.emplace_back(tetA, tetB, tetC);
  facets.emplace_back(tetA, tetB, tetD);
  facets.emplace_back(tetA, tetB, stk::math::Vector3d(0,-1,0));
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), std::nullopt);

  // The shared edge A-B (count 3) plus the six free outer edges (count 1) are all non-watertight.
  // Verify the shared edge is among them.
  const bool sharedEdgeFlagged = std::any_of(badEdges.begin(), badEdges.end(),
      [](const std::array<stk::math::Vector3d,2> & e)
      {
        auto matches = [](const stk::math::Vector3d & p, const stk::math::Vector3d & q)
          { return (p-q).length() < 1.e-9; };
        return (matches(e[0],tetA) && matches(e[1],tetB)) || (matches(e[0],tetB) && matches(e[1],tetA));
      });
  EXPECT_TRUE(sharedEdgeFlagged);
}

TEST(FacetedSurface, watertight_boundingBoxRestrictsFlaggedEdges)
{
  if (!is_serial()) return;
  std::vector<Facet3d> facets = build_tetra_faces();
  facets.pop_back(); // boundary edges A-B, B-C, C-A

  // Box tightly around edge A-B = (0,0,0)->(1,0,0); C=(0,1,0) is excluded, so edges B-C and C-A
  // (which touch C) fail the both-endpoints test and only A-B remains.
  const std::optional<BoundingBox> box(BoundingBox(stk::math::Vector3d(-0.1,-0.1,-0.1), stk::math::Vector3d(1.1,0.1,0.1)));
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), box);
  EXPECT_EQ(1u, badEdges.size());
}

TEST(FacetedSurface, watertight_edgesOnBoundingBoxFaceAreNotFlagged)
{
  if (!is_serial()) return;
  std::vector<Facet3d> facets = build_tetra_faces();
  facets.pop_back(); // boundary edges A-B, B-C, C-A all lie in the z=0 plane

  // Box whose min-z face is exactly the z=0 plane and that contains all three boundary edges.
  // Those edges are the open boundary of the clip against this box, not cracks -> none flagged.
  const std::optional<BoundingBox> box(BoundingBox(stk::math::Vector3d(-0.1,-0.1,0.), stk::math::Vector3d(1.1,1.1,1.0)));
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), box);
  EXPECT_TRUE(badEdges.empty());
}

TEST(FacetedSurface, watertight_edgeInsideBoxButOffItsFacesIsStillFlagged)
{
  if (!is_serial()) return;
  std::vector<Facet3d> facets = build_tetra_faces();
  facets.pop_back(); // boundary edges A-B, B-C, C-A in the z=0 plane

  // Same z=0 boundary edges, but lift the box's min-z face below the plane so the edges sit in the
  // box interior rather than on a face: they are genuine cracks and must still be flagged.
  const std::optional<BoundingBox> box(BoundingBox(stk::math::Vector3d(-0.1,-0.1,-1.0), stk::math::Vector3d(1.1,1.1,1.0)));
  const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, serial_edge_cell_size(facets), box);
  EXPECT_EQ(3u, badEdges.size());
}

TEST(FacetedSurface, watertight_2d_closedLoopVsOpenChain)
{
  if (!is_serial()) return;
  // Closed triangle loop: each endpoint shared by exactly two segments -> watertight.
  {
    std::vector<Facet2d> facets;
    facets.emplace_back(tetA, tetB);
    facets.emplace_back(tetB, tetC);
    facets.emplace_back(tetC, tetA);
    const double cell = compute_global_edge_quantization_cell_size(MPI_COMM_WORLD, facets);
    const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, cell, std::nullopt);
    EXPECT_TRUE(badEdges.empty());
  }
  // Open chain: the two free ends (A and C) are each shared by only one segment.
  {
    std::vector<Facet2d> facets;
    facets.emplace_back(tetA, tetB);
    facets.emplace_back(tetB, tetC);
    const double cell = compute_global_edge_quantization_cell_size(MPI_COMM_WORLD, facets);
    const auto badEdges = find_non_watertight_edges(facets, MPI_COMM_WORLD, cell, std::nullopt);
    EXPECT_EQ(2u, badEdges.size());
  }
}

TEST(FacetedSurface, watertight_flagWiresThroughPrepareToCompute)
{
  if (!is_serial()) return;
  // An open shell with the flag on must complete prepare_to_compute (emitting a warning) and
  // remain queryable; with the flag off (default) nothing changes.  We assert it runs cleanly.
  std::unique_ptr<FacetedSurfaceBase> surface = FacetedSurfaceBase::build(3);
  surface->emplace_back_3d(tetB, tetC, tetD);
  surface->emplace_back_3d(tetA, tetC, tetD);
  surface->emplace_back_3d(tetA, tetB, tetD); // open shell (missing face opposite D)
  surface->set_check_watertight(true);
  surface->prepare_to_compute(queryBox, 0.);

  EXPECT_EQ(3u, surface->size());
}

// Exercises the parallel edge-key exchange: the tetra faces are partitioned disjointly across
// ranks (mimicking myLocalFacets after parallel_distribute_facets), so a shared edge is generally
// contributed by facets on different ranks and must be reconciled through the CommSparse owner
// exchange.  The global bad-edge count must match the serial result regardless of proc count.
TEST(FacetedSurface, watertight_parallelEdgeExchange_countsAreProcCountInvariant)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);

  auto local_partition = [numProcs, me](const std::vector<Facet3d> & all)
  {
    std::vector<Facet3d> mine;
    for (size_t i=0; i<all.size(); ++i)
      if (static_cast<int>(i % numProcs) == me)
        mine.push_back(all[i]);
    return mine;
  };

  auto global_bad_count = [](const std::vector<Facet3d> & localFacets)
  {
    const double cell = compute_global_edge_quantization_cell_size(MPI_COMM_WORLD, localFacets);
    const auto badEdges = find_non_watertight_edges(localFacets, MPI_COMM_WORLD, cell, std::nullopt);
    size_t count = badEdges.size();
    krino::all_reduce_sum(MPI_COMM_WORLD, count);
    return count;
  };

  // Closed tetra: watertight -> 0 bad edges on any proc count.
  {
    const std::vector<Facet3d> mine = local_partition(build_tetra_faces());
    EXPECT_EQ(0u, global_bad_count(mine));
  }
  // Open shell (drop one face) -> exactly 3 boundary edges, independent of proc count.
  {
    std::vector<Facet3d> all = build_tetra_faces();
    all.pop_back();
    const std::vector<Facet3d> mine = local_partition(all);
    EXPECT_EQ(3u, global_bad_count(mine));
  }
}

// default_watertight_bounding_box_from_domain fills in the restriction box only for faceted
// surfaces that have the watertight check on and no box of their own; it leaves user-set boxes,
// check-off surfaces, and non-faceted surfaces untouched.
TEST(FacetedSurface, watertight_defaultBoundingBoxFromDomainOnlyFillsUnsetCheckedFacets)
{
  if (!is_serial()) return;

  const BoundingBox domainBox(stk::math::Vector3d(-1,-1,-1), stk::math::Vector3d(1,1,1));
  const BoundingBox userBox(stk::math::Vector3d(-2,-2,-2), stk::math::Vector3d(2,2,2));

  auto make_facet_surface = [](){
    std::shared_ptr<FacetedSurfaceBase> s = FacetedSurfaceBase::build(3);
    s->emplace_back_3d(tetA, tetB, tetC);
    return s;
  };

  // (a) check on, no box -> gets the domain box.
  std::shared_ptr<FacetedSurfaceBase> checkedNoBox = make_facet_surface();
  checkedNoBox->set_check_watertight(true);
  // (b) check on, own box -> keeps its box.
  std::shared_ptr<FacetedSurfaceBase> checkedWithBox = make_facet_surface();
  checkedWithBox->set_check_watertight(true);
  checkedWithBox->set_watertight_bounding_box(userBox);
  // (c) check off -> stays unset.
  std::shared_ptr<FacetedSurfaceBase> uncheckedNoBox = make_facet_surface();
  // (d) non-faceted surface -> skipped without error.
  std::shared_ptr<Surface> sphere = std::make_shared<Sphere>(stk::math::Vector3d(0,0,0), 1.0);

  std::vector<std::shared_ptr<Surface>> surfaces = {checkedNoBox, checkedWithBox, uncheckedNoBox, sphere};
  default_watertight_bounding_box_from_domain(surfaces, domainBox);

  auto expect_box_eq = [](const BoundingBox & expected, const BoundingBox & actual){
    for (int i=0; i<3; ++i)
    {
      EXPECT_FLOAT_EQ(expected.get_min()[i], actual.get_min()[i]);
      EXPECT_FLOAT_EQ(expected.get_max()[i], actual.get_max()[i]);
    }
  };

  ASSERT_TRUE(checkedNoBox->get_watertight_bounding_box().has_value());
  expect_box_eq(domainBox, *checkedNoBox->get_watertight_bounding_box());

  ASSERT_TRUE(checkedWithBox->get_watertight_bounding_box().has_value());
  expect_box_eq(userBox, *checkedWithBox->get_watertight_bounding_box());

  EXPECT_FALSE(uncheckedNoBox->get_watertight_bounding_box().has_value());
}

}
