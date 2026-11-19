// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_Unit_SurfaceDecomposition.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

// A rotated cuboid well inside the domain: with its edges and corners captured, including where
// mesh edges cross it twice beside its sharp edges, the decomposition's negative phase is exactly
// the cuboid, so its volume and surface area are exact.
TEST(RotatedCuboidDecomposition, snapping_volumeAndAreaOfTheCuboid)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain, 0.2);

  const stk::math::Vector3d dimensions(1.0, 0.7, 0.8);
  const Cuboid cuboid(
      stk::math::Vector3d(0.03, -0.02, 0.05), dimensions, stk::math::Vector3d(0.3, 0.5, 0.2));
  const DecompositionMeasures measures = decomposition.decompose_and_measure(cuboid);

  const double exactVolume = dimensions[0] * dimensions[1] * dimensions[2];
  const double exactArea = 2. *
      (dimensions[0] * dimensions[1] + dimensions[0] * dimensions[2] +
          dimensions[1] * dimensions[2]);
  EXPECT_NEAR(exactVolume, measures.negativeVolume, 1.e-10 * exactVolume);
  EXPECT_NEAR(exactArea, measures.interfaceArea, 1.e-10 * exactArea);
  EXPECT_LT(0., measures.worstQuality);
}

namespace
{
const stk::math::Vector3d rotatedCuboidDimensions(1.0, 0.7, 0.8);
const double rotatedCuboidVolume =
    rotatedCuboidDimensions[0] * rotatedCuboidDimensions[1] * rotatedCuboidDimensions[2];
const double rotatedCuboidArea = 2. *
    (rotatedCuboidDimensions[0] * rotatedCuboidDimensions[1] +
        rotatedCuboidDimensions[0] * rotatedCuboidDimensions[2] +
        rotatedCuboidDimensions[1] * rotatedCuboidDimensions[2]);

DecompositionMeasures
decompose_rotated_cuboid(const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::optional<bool> cutSharpFeaturesBeforeSurfaces,
    const std::optional<double> edgeTol = std::nullopt,
    const stk::topology elemTopology = stk::topology::TETRAHEDRON_4)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain,
      0.2,
      edgeDegeneracyHandling,
      cutSharpFeaturesBeforeSurfaces,
      elemTopology,
      edgeTol);
  const Cuboid cuboid(stk::math::Vector3d(0.03, -0.02, 0.05),
      rotatedCuboidDimensions,
      stk::math::Vector3d(0.3, 0.5, 0.2));
  return decomposition.decompose_and_measure(cuboid);
}
} // namespace

// Only cutting, with its sharp features cut before it: the cuboid's corners and edges are captured
// to within the edge tolerance, so its volume and area are close to exact, and no leaf is a sliver.
TEST(RotatedCuboidDecomposition, onlyCuttingWithSharpFeaturesFirst_volumeAndAreaNearlyExact)
{
  const DecompositionMeasures measures = decompose_rotated_cuboid(SNAP_TO_NODE, true);
  EXPECT_NEAR(rotatedCuboidVolume, measures.negativeVolume, 1.e-5 * rotatedCuboidVolume);
  EXPECT_NEAR(rotatedCuboidArea, measures.interfaceArea, 1.e-5 * rotatedCuboidArea);
  EXPECT_LT(1.e-5, measures.worstQuality);
}

// Only cutting, with the sharp features not cut first: the cuboid's edges are cut off.
TEST(RotatedCuboidDecomposition, onlyCuttingWithoutSharpFeaturesFirst_sharpFeaturesNotCaptured)
{
  const DecompositionMeasures measures = decompose_rotated_cuboid(SNAP_TO_NODE, false);
  EXPECT_LT(1.e-3 * rotatedCuboidVolume, std::abs(rotatedCuboidVolume - measures.negativeVolume));
  EXPECT_LT(0., measures.worstQuality);
}

// Snapping, with the sharp features not cut first: it decomposes, but the edges are cut off.
TEST(RotatedCuboidDecomposition, snappingWithoutSharpFeaturesFirst_sharpFeaturesNotCaptured)
{
  const DecompositionMeasures measures =
      decompose_rotated_cuboid(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, false);
  EXPECT_LT(1.e-3 * rotatedCuboidVolume, std::abs(rotatedCuboidVolume - measures.negativeVolume));
  EXPECT_LT(0., measures.worstQuality);
}

// Snapping, with the sharp features captured only by snapping to them (none is cut): it decomposes,
// near the cuboid, and right after, it does not need updating.
TEST(RotatedCuboidDecomposition, snappingCapturingSharpFeaturesOnlyDuringSnapping_nearTheCuboidAndNoUpdateNeeded)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain, 0.2);
  decomposition.capture_sharp_features_only_during_snapping();
  const Cuboid cuboid(stk::math::Vector3d(0.03, -0.02, 0.05),
      rotatedCuboidDimensions,
      stk::math::Vector3d(0.3, 0.5, 0.2));
  const DecompositionMeasures measures = decomposition.decompose_and_measure(cuboid);
  EXPECT_NEAR(rotatedCuboidVolume, measures.negativeVolume, 3.e-3 * rotatedCuboidVolume);
  EXPECT_NEAR(rotatedCuboidArea, measures.interfaceArea, 3.e-3 * rotatedCuboidArea);
  EXPECT_LT(0., measures.worstQuality);
  EXPECT_FALSE(decomposition.decomposition_needs_update(cuboid));
}

// Nothing has moved since the decomposition, so it does not need updating, although nodes at the
// cuboid's edges and corners, and where mesh edges cross it twice beside its edges, are on it
// without a crossing of their edges there.
TEST(RotatedCuboidDecomposition, unmovedCuboid_decompositionDoesNotNeedUpdate)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain, 0.2);
  const Cuboid cuboid(stk::math::Vector3d(0.03, -0.02, 0.05),
      stk::math::Vector3d(1.0, 0.7, 0.8),
      stk::math::Vector3d(0.3, 0.5, 0.2));
  decomposition.decompose_and_measure(cuboid);
  EXPECT_FALSE(decomposition.decomposition_needs_update(cuboid));
}

// A box between the planes x = xLo and x = xHi, reaching beyond the domain in y and z, as a
// connected faceted surface without sharp features (no sharp feature angle).  The facets are all
// owned by processor 0.
static std::unique_ptr<MeshFacetedSurface> thin_faceted_slab_normal_to_x(const double xLo, const double xHi)
{
  constexpr double w = 5.;
  const std::vector<stk::math::Vector3d> vertices{{xLo, -w, -w},
      {xHi, -w, -w},
      {xHi, w, -w},
      {xLo, w, -w},
      {xLo, -w, w},
      {xHi, -w, w},
      {xHi, w, w},
      {xLo, w, w}};
  std::vector<std::array<unsigned, 3>> facets;
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
    facets = {{{0, 4, 7}},
        {{0, 7, 3}},
        {{1, 2, 6}},
        {{1, 6, 5}},
        {{0, 1, 5}},
        {{0, 5, 4}},
        {{2, 3, 7}},
        {{2, 7, 6}},
        {{0, 3, 2}},
        {{0, 2, 1}},
        {{4, 5, 6}},
        {{4, 6, 7}}};
  return std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, vertices, facets);
}

// A slab without sharp features thinner than the tolerance crosses each mesh edge through it twice
// within the tolerance: close crossings, which are not cut, so the decomposition is the mesh.
// Nothing has moved since the decomposition, so it does not need updating.
TEST(RotatedCuboidDecomposition, unmovedSlabThinnerThanTolerance_decompositionDoesNotNeedUpdate)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain, 0.2);
  const auto slab = thin_faceted_slab_normal_to_x(0.013, 0.013 + 2.e-5);
  ASSERT_FALSE(slab->has_sharp_features());
  decomposition.decompose_and_measure(*slab);
  EXPECT_FALSE(decomposition.decomposition_needs_update(*slab));
}

TEST(RotatedCuboidDecomposition, movedCuboid_decompositionNeedsUpdate)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  SurfaceDecomposition decomposition(domain, 0.2);
  const stk::math::Vector3d dimensions(1.0, 0.7, 0.8);
  const stk::math::Vector3d rotation(0.3, 0.5, 0.2);
  decomposition.decompose_and_measure(
      Cuboid(stk::math::Vector3d(0.03, -0.02, 0.05), dimensions, rotation));
  EXPECT_TRUE(decomposition.decomposition_needs_update(
      Cuboid(stk::math::Vector3d(0.08, -0.02, 0.05), dimensions, rotation)));
}

namespace
{
const double rotatedRectangleArea = 1.0 * 0.7;
const double rotatedRectanglePerimeter = 2. * (1.0 + 0.7);

DecompositionMeasures
decompose_rotated_rectangle(const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::optional<bool> cutSharpFeaturesBeforeSurfaces,
    const std::optional<double> edgeTol = std::nullopt,
    const stk::topology elemTopology = stk::topology::TRIANGLE_3_2D)
{
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-1., -1., 0.), stk::math::Vector3d(1., 1., 0.));
  SurfaceDecomposition decomposition(domain,
      0.2,
      edgeDegeneracyHandling,
      cutSharpFeaturesBeforeSurfaces,
      elemTopology,
      edgeTol);
  const Cuboid rectangle(stk::math::Vector3d(0.03, -0.02, 0.),
      stk::math::Vector3d(1.0, 0.7, 10.),
      stk::math::Vector3d(0., 0., 0.3));
  return decomposition.decompose_and_measure(rectangle);
}
} // namespace

// In 2D, a rotated rectangle well inside the domain: with its corners captured, the decomposition's
// negative phase is exactly the rectangle, so its area and perimeter are exact.
TEST(RotatedRectangleDecomposition, snapping_areaAndPerimeterOfTheRectangle)
{
  const DecompositionMeasures measures = decompose_rotated_rectangle(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, std::nullopt);
  EXPECT_NEAR(rotatedRectangleArea, measures.negativeVolume, 1.e-10 * rotatedRectangleArea);
  EXPECT_NEAR(
      rotatedRectanglePerimeter, measures.interfaceArea, 1.e-10 * rotatedRectanglePerimeter);
  EXPECT_LT(0., measures.worstQuality);
}

// In 2D, only cutting, with its sharp features cut before it: the rectangle's corners are captured
// to within the edge tolerance, so its area and perimeter are close to exact.
TEST(RotatedRectangleDecomposition, onlyCuttingWithSharpFeaturesFirst_areaAndPerimeterNearlyExact)
{
  const DecompositionMeasures measures = decompose_rotated_rectangle(SNAP_TO_NODE, true);
  EXPECT_NEAR(rotatedRectangleArea, measures.negativeVolume, 1.e-5 * rotatedRectangleArea);
  EXPECT_NEAR(rotatedRectanglePerimeter, measures.interfaceArea, 1.e-5 * rotatedRectanglePerimeter);
  EXPECT_LT(1.e-3, measures.worstQuality);
}

// In 2D, only cutting, with the sharp features not cut first: the corners are cut off.
TEST(RotatedRectangleDecomposition, onlyCuttingWithoutSharpFeaturesFirst_cornersNotCaptured)
{
  const DecompositionMeasures measures = decompose_rotated_rectangle(SNAP_TO_NODE, false);
  EXPECT_LT(1.e-3 * rotatedRectangleArea, std::abs(rotatedRectangleArea - measures.negativeVolume));
  EXPECT_LT(0., measures.worstQuality);
}

// With a large edge tolerance, the sharp features are moved by up to 15% of an edge onto the
// entities they are near: the decomposition is no longer exact, but stays close to the shape, with
// no slivers, whether snapping or only cutting (with the sharp features cut first).
constexpr double largeEdgeTol = 0.15;

TEST(RotatedCuboidDecomposition, snappingWithLargeEdgeTolerance_volumeAndAreaNearTheCuboid)
{
  const DecompositionMeasures measures = decompose_rotated_cuboid(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, std::nullopt, largeEdgeTol);
  EXPECT_NEAR(rotatedCuboidVolume, measures.negativeVolume, 0.02 * rotatedCuboidVolume);
  EXPECT_NEAR(rotatedCuboidArea, measures.interfaceArea, 0.03 * rotatedCuboidArea);
  EXPECT_LT(1.e-3, measures.worstQuality);
}

// Only cutting, a sharp edge moved within the tolerance leaves the leaf edges that stand for it
// beside it, crossing the cuboid's faces near it twice within the tolerance: those close crossings
// are not cut, which would make edges much shorter than the tolerance (worst quality 7e-3).
TEST(RotatedCuboidDecomposition,
    onlyCuttingWithSharpFeaturesFirstAndLargeEdgeTolerance_volumeAndAreaNearTheCuboid)
{
  const DecompositionMeasures measures = decompose_rotated_cuboid(SNAP_TO_NODE, true, largeEdgeTol);
  EXPECT_NEAR(rotatedCuboidVolume, measures.negativeVolume, 0.02 * rotatedCuboidVolume);
  EXPECT_NEAR(rotatedCuboidArea, measures.interfaceArea, 0.03 * rotatedCuboidArea);
  EXPECT_LT(1.e-2, measures.worstQuality);
}

TEST(RotatedRectangleDecomposition, snappingWithLargeEdgeTolerance_areaAndPerimeterNearTheRectangle)
{
  const DecompositionMeasures measures = decompose_rotated_rectangle(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, std::nullopt, largeEdgeTol);
  EXPECT_NEAR(rotatedRectangleArea, measures.negativeVolume, 0.02 * rotatedRectangleArea);
  EXPECT_NEAR(rotatedRectanglePerimeter, measures.interfaceArea, 0.03 * rotatedRectanglePerimeter);
  EXPECT_LT(1.e-3, measures.worstQuality);
}

TEST(RotatedRectangleDecomposition,
    onlyCuttingWithSharpFeaturesFirstAndLargeEdgeTolerance_areaAndPerimeterNearTheRectangle)
{
  const DecompositionMeasures measures =
      decompose_rotated_rectangle(SNAP_TO_NODE, true, largeEdgeTol);
  EXPECT_NEAR(rotatedRectangleArea, measures.negativeVolume, 0.02 * rotatedRectangleArea);
  EXPECT_NEAR(rotatedRectanglePerimeter, measures.interfaceArea, 0.03 * rotatedRectanglePerimeter);
  EXPECT_LT(1.e-3, measures.worstQuality);
}

namespace
{
// Quadratic elements are decomposed as the linear elements of their vertices: the same leaves, the
// same phase volume and interface area.
void expect_same_decomposition(const DecompositionMeasures & linear, const DecompositionMeasures & quadratic)
{
  EXPECT_EQ(linear.numElements, quadratic.numElements);
  EXPECT_NEAR(linear.negativeVolume, quadratic.negativeVolume, 1.e-12 * linear.negativeVolume);
  EXPECT_NEAR(linear.interfaceArea, quadratic.interfaceArea, 1.e-12 * linear.interfaceArea);
}
}

TEST(RotatedCuboidDecomposition, tet10OnlyCuttingWithoutSharpFeaturesFirst_sameAsTet4)
{
  expect_same_decomposition(decompose_rotated_cuboid(SNAP_TO_NODE, false),
      decompose_rotated_cuboid(SNAP_TO_NODE, false, std::nullopt, stk::topology::TETRAHEDRON_10));
}

TEST(RotatedCuboidDecomposition, tet10OnlyCuttingWithSharpFeaturesFirst_sameAsTet4)
{
  expect_same_decomposition(decompose_rotated_cuboid(SNAP_TO_NODE, true),
      decompose_rotated_cuboid(SNAP_TO_NODE, true, std::nullopt, stk::topology::TETRAHEDRON_10));
}

TEST(RotatedCuboidDecomposition, tet10Snapping_sameAsTet4)
{
  const Edge_Degeneracy_Handling snapping = SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE;
  expect_same_decomposition(decompose_rotated_cuboid(snapping, std::nullopt),
      decompose_rotated_cuboid(snapping, std::nullopt, std::nullopt, stk::topology::TETRAHEDRON_10));
}

TEST(RotatedRectangleDecomposition, tri6OnlyCuttingWithoutSharpFeaturesFirst_sameAsTri3)
{
  expect_same_decomposition(decompose_rotated_rectangle(SNAP_TO_NODE, false),
      decompose_rotated_rectangle(SNAP_TO_NODE, false, std::nullopt, stk::topology::TRIANGLE_6_2D));
}

TEST(RotatedRectangleDecomposition, tri6OnlyCuttingWithSharpFeaturesFirst_sameAsTri3)
{
  expect_same_decomposition(decompose_rotated_rectangle(SNAP_TO_NODE, true),
      decompose_rotated_rectangle(SNAP_TO_NODE, true, std::nullopt, stk::topology::TRIANGLE_6_2D));
}

TEST(RotatedRectangleDecomposition, tri6Snapping_sameAsTri3)
{
  const Edge_Degeneracy_Handling snapping = SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE;
  expect_same_decomposition(decompose_rotated_rectangle(snapping, std::nullopt),
      decompose_rotated_rectangle(snapping, std::nullopt, std::nullopt, stk::topology::TRIANGLE_6_2D));
}

} // namespace krino
