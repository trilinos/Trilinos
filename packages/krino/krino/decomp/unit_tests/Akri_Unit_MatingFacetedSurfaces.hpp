// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_MATINGFACETEDSURFACES_HPP_
#define KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_MATINGFACETEDSURFACES_HPP_

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino
{

// Faceted surfaces of two halves of a box (or a rectangle in 2D) that mate along the middle face (or
// segment) between them, as the bounding surfaces of two phases that meet, with the elsewhere phase
// around them.  The facets are all owned by processor 0.

// The 12 triangles of a box with these 8 vertices, ordered as Cuboid::vertex_location orders them
// (the bottom 4 counterclockwise, then the top 4), facing out.
inline std::vector<std::array<unsigned, 3>> mating_box_facets()
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

// The bottom and top halves of the cuboid, split halfway between its bottom and top faces, as connected
// faceted surfaces that mate: they share the 4 vertices of the middle face (the same doubles).  With
// matchingMiddleFaces they also share its 2 triangles, in opposite orientations; otherwise the top
// half splits the middle face along its other diagonal, so the halves coincide there without matching
// facets.
inline std::array<std::unique_ptr<MeshFacetedSurface>, 2> mating_halves_of_cuboid(
    const Cuboid & cuboid, const bool matchingMiddleFaces = true)
{
  std::array<stk::math::Vector3d, 8> corners;
  for (unsigned i = 0; i < 8; ++i)
    corners[i] = cuboid.vertex_location(i);
  std::array<stk::math::Vector3d, 4> middle;
  for (unsigned i = 0; i < 4; ++i)
    middle[i] = 0.5 * (corners[i] + corners[4 + i]);

  const std::vector<stk::math::Vector3d> bottomVertices{
      corners[0], corners[1], corners[2], corners[3], middle[0], middle[1], middle[2], middle[3]};
  const std::vector<stk::math::Vector3d> topVertices{
      middle[0], middle[1], middle[2], middle[3], corners[4], corners[5], corners[6], corners[7]};
  const bool ownsFacets = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  const std::vector<std::array<unsigned, 3>> bottomFacets =
      ownsFacets ? mating_box_facets() : std::vector<std::array<unsigned, 3>>{};
  std::vector<std::array<unsigned, 3>> topFacets = bottomFacets;
  if (ownsFacets && !matchingMiddleFaces)
  {
    // The bottom face of the top half, facing down, split along the diagonal from vertex 1 to 3.
    topFacets[0] = {{0, 3, 1}};
    topFacets[1] = {{1, 3, 2}};
  }
  return {std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, bottomVertices, bottomFacets),
      std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, topVertices, topFacets)};
}

// The vertex of a rectangle of these dimensions, centered at center and rotated counterclockwise by
// angle, at these coordinates in units of its half dimensions.
inline stk::math::Vector3d mating_rectangle_point(const stk::math::Vector3d & center,
    const std::array<double, 2> & dimensions,
    const double angle,
    const double x,
    const double y)
{
  const double u = 0.5 * x * dimensions[0];
  const double v = 0.5 * y * dimensions[1];
  return center + stk::math::Vector3d(std::cos(angle) * u - std::sin(angle) * v, std::sin(angle) * u + std::cos(angle) * v, 0.);
}

// In 2D, the bottom and top halves of the rotated rectangle, split halfway between its bottom and top
// sides, as connected loops of segments, counterclockwise (facing out), that mate: they share the 2
// vertices of the middle segment (the same doubles).  With matchingMiddleSegments both have that
// segment, in opposite directions; otherwise the top half splits it at its midpoint, so the halves
// coincide there without matching segments.
inline std::array<std::unique_ptr<MeshFacetedSurface>, 2> mating_halves_of_rectangle(const stk::math::Vector3d & center,
    const std::array<double, 2> & dimensions,
    const double angle,
    const bool matchingMiddleSegments = true)
{
  auto point = [&](const double x, const double y)
  { return mating_rectangle_point(center, dimensions, angle, x, y); };
  const stk::math::Vector3d middleLeft = point(-1., 0.);
  const stk::math::Vector3d middleRight = point(1., 0.);

  const std::vector<stk::math::Vector3d> bottomVertices{point(-1., -1.), point(1., -1.), middleRight, middleLeft};
  std::vector<stk::math::Vector3d> topVertices{middleLeft, middleRight, point(1., 1.), point(-1., 1.)};
  const bool ownsFacets = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<std::array<unsigned, 3>> bottomSegments, topSegments;
  if (ownsFacets)
  {
    bottomSegments = {{{0, 1, 0}}, {{1, 2, 0}}, {{2, 3, 0}}, {{3, 0, 0}}};
    topSegments = {{{0, 1, 0}}, {{1, 2, 0}}, {{2, 3, 0}}, {{3, 0, 0}}};
  }
  if (!matchingMiddleSegments)
  {
    topVertices.push_back(0.5 * (middleLeft + middleRight));
    if (ownsFacets)
    {
      topSegments[0] = {{0, 4, 0}};
      topSegments.push_back({{4, 1, 0}});
    }
  }
  return {std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 2, bottomVertices, bottomSegments),
      std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 2, topVertices, topSegments)};
}

} // namespace krino

#endif /* KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_MATINGFACETEDSURFACES_HPP_ */
