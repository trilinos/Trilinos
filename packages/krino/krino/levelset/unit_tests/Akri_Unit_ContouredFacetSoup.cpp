// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_ClosestPointRedistance.hpp>
#include <Akri_FacetCoordinates.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_Unit_BoundingBoxMesh.hpp>
#include <Akri_UnitMeshUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace krino
{

namespace
{

// Points in the box, the same on every processor, offset by irrational fractions so that none is on
// a facet edge or vertex.
std::vector<stk::math::Vector3d>
grid_points_in_box(const BoundingBox & box, const unsigned dim, const unsigned n)
{
  std::vector<stk::math::Vector3d> points;
  const auto & lo = box.get_min();
  const auto & hi = box.get_max();
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      for (unsigned k = 0; k < (dim == 3 ? n : 1u); ++k)
      {
        const double fx = (i + 0.3819660113) / n;
        const double fy = (j + 0.2360679775) / n;
        const double fz = (k + 0.6180339887) / n;
        points.emplace_back(lo[0] + fx * (hi[0] - lo[0]),
            lo[1] + fy * (hi[1] - lo[1]),
            (dim == 3) ? lo[2] + fz * (hi[2] - lo[2]) : 0.);
      }
  return points;
}

// The distance from x to the nearest of all processors' facets, found by checking each one.
template <class FACET>
double brute_force_distance_to_facets(
    const stk::ParallelMachine comm, const FacetCoordinates & facets, const stk::math::Vector3d & x)
{
  double localMinSqrDist = std::numeric_limits<double>::max();
  for (const FACET & facet : facets.make_facets<FACET>())
    localMinSqrDist = std::min(localMinSqrDist, facet.point_distance_squared(x));
  double minSqrDist = 0.;
  stk::all_reduce_min(comm, &localMinSqrDist, &minSqrDist, 1);
  return std::sqrt(minSqrDist);
}

// Contours the level set of distance from a sphere (circle in 2D) on the bounding box mesh into a
// facet soup, and expects the soup as an unconnected MeshFacetedSurface to give the distance to the
// nearest facet (truncated to the narrow band), a closest point at that distance, and away from the
// interface the sign of the level set.
void expect_contoured_soup_surface_distances_to_be_those_of_its_facets(BoundingBoxMesh & bboxMesh)
{
  stk::mesh::MetaData & meta = bboxMesh.meta_data();
  const unsigned dim = meta.spatial_dimension();
  const double meshSize = 0.2;
  populate_bounding_box_mesh_and_activate(
      bboxMesh, {0., 0., 0.}, {1., 1., (dim == 3) ? 1. : 0.}, meshSize);
  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();

  meta.enable_late_fields();
  auto & distField = meta.declare_field<double>(stk::topology::NODE_RANK, "distance");
  stk::mesh::put_field_on_mesh(distField, meta.universal_part(), 1, 1, nullptr);
  const FieldRef coords = AuxMetaData::get(meta).get_current_coordinates();
  const stk::math::Vector3d center(0.45, 0.5, (dim == 3) ? 0.55 : 0.);
  const double radius = 0.3;
  auto levelset = [&](const stk::math::Vector3d & x) { return (x - center).length() - radius; };
  {
    auto coordsData = coords.data<double>();
    auto distData = distField.data<stk::mesh::ReadWrite>();
    for (auto * bucket : mesh.buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
        distData.entity_values(node)() =
            levelset(stk::math::Vector3d(coordsData.entity_values(node), dim));
  }

  std::vector<stk::mesh::Entity> elems;
  stk::mesh::get_selected_entities(AuxMetaData::get(meta).active_part() & meta.locally_owned_part(),
      mesh.buckets(stk::topology::ELEMENT_RANK),
      elems);
  FacetCoordinates soup(dim);
  ClosestPointRedistance::build_facets_for_elements(
      mesh, coords, FieldRef(distField), elems, meshSize, soup);
  size_t numFacets = soup.num_facets(), globalNumFacets = 0;
  stk::all_reduce_sum(mesh.parallel(), &numFacets, &globalNumFacets, 1);
  EXPECT_LT(0u, globalNumFacets);

  MeshFacetedSurface surface(mesh.parallel(), dim, soup.coordinates());
  // The same query points on every processor, for the brute-force distance over all of them.
  const BoundingBox queryBox(
      stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., (dim == 3) ? 1. : 0.));
  for (const double narrowBand : {0., 0.15})
  {
    const double farField = narrowBand;
    surface.prepare_to_compute(0., queryBox, narrowBand);
    for (auto & x : grid_points_in_box(queryBox, dim, 7))
    {
      const double distance = (dim == 3) ?
          brute_force_distance_to_facets<Facet3d>(mesh.parallel(), soup, x) :
          brute_force_distance_to_facets<Facet2d>(mesh.parallel(), soup, x);
      const bool isInBand = 0. == narrowBand || distance <= narrowBand;
      const double expected = isInBand ? distance : farField;
      EXPECT_NEAR(expected, surface.point_unsigned_distance(x, narrowBand, farField), 1.e-12)
          << "narrow band " << narrowBand << " at " << x;
      const double signedDistance =
          surface.possibly_truncated_point_signed_distance(x, narrowBand, farField);
      if (isInBand)
      {
        EXPECT_NEAR(expected, std::abs(signedDistance), 1.e-12)
            << "narrow band " << narrowBand << " at " << x;
        if (std::abs(levelset(x)) > meshSize)
        {
          EXPECT_EQ(levelset(x) > 0., signedDistance > 0.)
              << "narrow band " << narrowBand << " at " << x;
        }
      }
      if (0. == narrowBand)
      {
        EXPECT_NEAR(distance, (surface.closest_point(x) - x).length(), 1.e-12) << x;
      }
    }
  }
}

} // namespace

TEST(ContouredFacetSoup, sphereLevelSetOnTets_distancesAreToTheNearestFacet)
{
  BoundingBoxMeshTet4 tets;
  expect_contoured_soup_surface_distances_to_be_those_of_its_facets(tets);
}

TEST(ContouredFacetSoup, circleLevelSetOnTris_distancesAreToTheNearestFacet)
{
  BoundingBoxMeshTri3 tris;
  expect_contoured_soup_surface_distances_to_be_those_of_its_facets(tris);
}

} // namespace krino
