// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_AdaptiveContourTri.hpp>
#include <Akri_AdaptiveContourTet.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Faceted_Surface.hpp>

#include <array>
#include <functional>

namespace krino {

// The adaptive contour drivers recursively subdivide a tri/tet and contour each
// piece against a semi-Lagrangian distance sampled from a functor. The exact
// number of facets produced is an implementation detail of the refinement, but
// two properties must always hold and are what these tests pin:
//   1. For a LINEAR distance field, every emitted facet vertex lies exactly on
//      the zero isosurface, regardless of how the element is subdivided.
//   2. A field that is strictly one sign over the element produces no facets.
//
// The drivers treat the field as a true (unit-gradient) signed distance -- the
// narrow-band "possibly cut" test compares a node's distance against half the
// edge length -- so the fields below are axis-aligned planes x - c whose
// gradient magnitude is exactly 1. That same narrow-band test also means the
// subdivision depth cannot be pushed arbitrarily high: once the subelements
// shrink below the band, they are skipped. minDepth=1 forces genuine
// subdivision while staying above that threshold. Here departureCoords ==
// coords (no advection), so the functor is sampled at the true point locations.

namespace {

// Signed distance to the plane x = 0.25: positive for x > 0.25.
double planar_distance(const stk::math::Vector3d & p) { return p[0] - 0.25; }

// A planar field that is NOT a unit-gradient signed distance: the gradient
// magnitude here is 2, so node values overstate the true distance to the
// interface (the plane x = 0.3). On the reference elements the endpoints of a
// crossing edge can both fall outside the narrow band yet straddle zero, so the
// narrow-band magnitude test alone would skip the crossing; the sign-change
// guard in have_possibly_cut_edge keeps it. The plane x = 0.3 is chosen to
// avoid coinciding with any midside sample point.
double scaled_nonunit_distance(const stk::math::Vector3d & p) { return 2.0 * (p[0] - 0.3); }

std::array<stk::math::Vector3d,3> tri_coords()
{
  RightTri spec;
  return {{ stk::math::Vector3d(spec.nodeLocs[0][0], spec.nodeLocs[0][1], 0.0),
            stk::math::Vector3d(spec.nodeLocs[1][0], spec.nodeLocs[1][1], 0.0),
            stk::math::Vector3d(spec.nodeLocs[2][0], spec.nodeLocs[2][1], 0.0) }};
}

std::array<stk::math::Vector3d,4> tet_coords()
{
  RightTet spec;
  return {{ spec.nodeLocs[0], spec.nodeLocs[1], spec.nodeLocs[2], spec.nodeLocs[3] }};
}

}

//--------------------------------------------------------------------------
// 2D driver
//--------------------------------------------------------------------------

// A planar field cutting the triangle produces facets, and every facet vertex
// sits on the isosurface to within tolerance -- even though minDepth forces
// several levels of subdivision.
TEST(AdaptiveContourTri, planarField_allFacetVerticesLieOnInterface)
{
  const std::array<stk::math::Vector3d,3> coords = tri_coords();
  const std::array<double,3> dist = {{ planar_distance(coords[0]), planar_distance(coords[1]), planar_distance(coords[2]) }};

  Faceted_Surface<Facet2d> facets;
  const double lengthScale = 1.0;
  adaptively_append_facets_for_tri_using_semilagrangian_distance(
      coords, coords, dist, planar_distance, lengthScale, facets,
      /*currentDepth*/ 0, /*minDepth*/ 1, /*maxDepth*/ 4);

  ASSERT_FALSE(facets.get_facets().empty());
  for (const Facet2d & facet : facets.get_facets())
    for (int v = 0; v < 2; ++v)
      EXPECT_NEAR(0.0, planar_distance(facet.facet_vertex(v)), 1.e-9)
          << "facet vertex must lie on the linear isosurface";
}

// A field that is strictly positive at every point of the triangle never cuts
// it, so no facets are produced.
TEST(AdaptiveContourTri, strictlyPositiveField_producesNoFacets)
{
  const std::array<stk::math::Vector3d,3> coords = tri_coords();
  const auto positive = [](const stk::math::Vector3d & p) { return p[0] + 5.0; };
  const std::array<double,3> dist = {{ positive(coords[0]), positive(coords[1]), positive(coords[2]) }};

  Faceted_Surface<Facet2d> facets;
  adaptively_append_facets_for_tri_using_semilagrangian_distance(
      coords, coords, dist, positive, 1.0, facets, 0, 0, 4);

  EXPECT_TRUE(facets.get_facets().empty());
}

// Robustness: a non-unit-gradient linear field that genuinely cuts the triangle
// must still produce facets (the sign-change guard prevents the narrow-band
// test from dropping the crossing). The field is spatially linear, so linear
// edge interpolation still lands each vertex exactly on the zero plane.
TEST(AdaptiveContourTri, nonUnitGradientField_stillProducesFacetsOnInterface)
{
  const std::array<stk::math::Vector3d,3> coords = tri_coords();
  const std::array<double,3> dist = {{ scaled_nonunit_distance(coords[0]), scaled_nonunit_distance(coords[1]),
                                       scaled_nonunit_distance(coords[2]) }};

  Faceted_Surface<Facet2d> facets;
  adaptively_append_facets_for_tri_using_semilagrangian_distance(
      coords, coords, dist, scaled_nonunit_distance, 1.0, facets, 0, 0, 2);

  ASSERT_FALSE(facets.get_facets().empty());
  for (const Facet2d & facet : facets.get_facets())
    for (int v = 0; v < 2; ++v)
      EXPECT_NEAR(0.0, scaled_nonunit_distance(facet.facet_vertex(v)), 1.e-9);
}

//--------------------------------------------------------------------------
// 3D driver
//--------------------------------------------------------------------------

TEST(AdaptiveContourTet, planarField_allFacetVerticesLieOnInterface)
{
  const std::array<stk::math::Vector3d,4> coords = tet_coords();
  const std::array<double,4> dist = {{ planar_distance(coords[0]), planar_distance(coords[1]),
                                       planar_distance(coords[2]), planar_distance(coords[3]) }};

  Faceted_Surface<Facet3d> facets;
  const double lengthScale = 1.0;
  adaptively_append_facets_for_tet_using_semilagrangian_distance(
      coords, coords, dist, planar_distance, lengthScale, facets,
      /*currentDepth*/ 0, /*minDepth*/ 1, /*maxDepth*/ 4);

  ASSERT_FALSE(facets.get_facets().empty());
  for (const Facet3d & facet : facets.get_facets())
    for (int v = 0; v < 3; ++v)
      EXPECT_NEAR(0.0, planar_distance(facet.facet_vertex(v)), 1.e-9)
          << "facet vertex must lie on the linear isosurface";
}

// 3D counterpart of the tri correctness test: a non-unit-gradient linear field
// must still yield facets that lie exactly on the zero plane. (Unlike the tri
// case this does not by itself require the sign-change guard -- the tet's
// longer diagonal edges give a wider narrow band that already catches this
// crossing -- but it pins the correctness of the tet path for such a field.)
TEST(AdaptiveContourTet, nonUnitGradientField_stillProducesFacetsOnInterface)
{
  const std::array<stk::math::Vector3d,4> coords = tet_coords();
  const std::array<double,4> dist = {{ scaled_nonunit_distance(coords[0]), scaled_nonunit_distance(coords[1]),
                                       scaled_nonunit_distance(coords[2]), scaled_nonunit_distance(coords[3]) }};

  Faceted_Surface<Facet3d> facets;
  adaptively_append_facets_for_tet_using_semilagrangian_distance(
      coords, coords, dist, scaled_nonunit_distance, 1.0, facets, 0, 0, 0);

  ASSERT_FALSE(facets.get_facets().empty());
  for (const Facet3d & facet : facets.get_facets())
    for (int v = 0; v < 3; ++v)
      EXPECT_NEAR(0.0, scaled_nonunit_distance(facet.facet_vertex(v)), 1.e-9);
}

// Regression: a corner-cutting plane x+y+z-0.5 places the interface exactly
// through the three midside nodes of the edges off node 0. When the adaptive
// driver decomposes the converged tet10, the corner subtet has three vertices
// exactly on the interface (caseId 13). That case must emit the zero face --
// before the fix it was silently dropped, leaving a hole. The base linear
// ContourTet already contours this field with a single facet, so the adaptive
// path must produce a non-empty surface whose vertices lie on the plane.
TEST(AdaptiveContourTet, interfaceThroughMidsideNodes_producesFacetsOnInterface)
{
  const auto corner_distance = [](const stk::math::Vector3d & p) { return p[0] + p[1] + p[2] - 0.5; };
  const std::array<stk::math::Vector3d,4> coords = tet_coords();
  const std::array<double,4> dist = {{ corner_distance(coords[0]), corner_distance(coords[1]),
                                       corner_distance(coords[2]), corner_distance(coords[3]) }};

  Faceted_Surface<Facet3d> facets;
  adaptively_append_facets_for_tet_using_semilagrangian_distance(
      coords, coords, dist, corner_distance, 1.0, facets, 0, 0, 0);

  ASSERT_FALSE(facets.get_facets().empty());
  for (const Facet3d & facet : facets.get_facets())
    for (int v = 0; v < 3; ++v)
      EXPECT_NEAR(0.0, corner_distance(facet.facet_vertex(v)), 1.e-9);
}

TEST(AdaptiveContourTet, strictlyNegativeField_producesNoFacets)
{
  const std::array<stk::math::Vector3d,4> coords = tet_coords();
  const auto negative = [](const stk::math::Vector3d & p) { return p[0] - 5.0; };
  const std::array<double,4> dist = {{ negative(coords[0]), negative(coords[1]),
                                       negative(coords[2]), negative(coords[3]) }};

  Faceted_Surface<Facet3d> facets;
  adaptively_append_facets_for_tet_using_semilagrangian_distance(
      coords, coords, dist, negative, 1.0, facets, 0, 0, 4);

  EXPECT_TRUE(facets.get_facets().empty());
}

}
