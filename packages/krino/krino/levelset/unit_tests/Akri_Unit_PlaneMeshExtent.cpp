// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_CSG_Surface.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_PlaneMeshExtent.hpp>
#include <Akri_StkMeshFixture.hpp>

namespace krino {

class CubeOf12TetsPlaneMeshExtent : public StkMeshTetFixture
{
protected:
  void build_mesh_spread_over_procs()
  {
    // Spread the elements over the processors, so the extent must be reduced over them.
    std::vector<int> elemProcOwners;
    for (size_t i=0; i<meshSpec.allElementConn.size(); ++i)
      elemProcOwners.push_back(i % parallel_size());
    build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, std::vector<unsigned>(meshSpec.allElementConn.size(), 1u), elemProcOwners);
  }

  FieldRef coords_field() const { return FieldRef(get_coordinates_field()); }

  CubeOf12Tets meshSpec; // unit cube centered on the origin
};

TEST_F(CubeOf12TetsPlaneMeshExtent, minMaxAlongDirection_globalOverAllProcs)
{
  build_mesh_spread_over_procs();

  const double tol = 1.e-12;
  const auto [minY, maxY] = compute_min_max_of_nodes_along_direction(get_mesh(), coords_field(), stk::math::Vector3d(0,1,0));
  EXPECT_NEAR(-0.5, minY, tol);
  EXPECT_NEAR(0.5, maxY, tol);

  // Along the cube diagonal the extremes are the corners (-0.5,-0.5,-0.5) and (0.5,0.5,0.5).
  const stk::math::Vector3d diagonal = stk::math::Vector3d(1,1,1).unit_vector();
  const auto [minDiag, maxDiag] = compute_min_max_of_nodes_along_direction(get_mesh(), coords_field(), diagonal);
  EXPECT_NEAR(-0.5*std::sqrt(3.), minDiag, tol);
  EXPECT_NEAR(0.5*std::sqrt(3.), maxDiag, tol);
}

TEST_F(CubeOf12TetsPlaneMeshExtent, setPlaneOffsets_planesAtFractionOfExtentPlacedSameOnAllProcs)
{
  build_mesh_spread_over_procs();

  auto plane = Plane::at_fraction_of_extent(stk::math::Vector3d(0,2,0), 0.75, 1.);
  set_plane_offsets_from_mesh_extent(get_mesh(), coords_field(), {plane.get()});

  EXPECT_FALSE(plane->needs_offset_from_extent());
  const double tol = 1.e-12;
  EXPECT_NEAR(0., plane->point_signed_distance(stk::math::Vector3d(0.3,0.25,-0.1)), tol); // y = -0.5 + 0.75*1
  EXPECT_NEAR(-0.75, plane->point_signed_distance(stk::math::Vector3d(0,-0.5,0)), tol);
}

TEST_F(CubeOf12TetsPlaneMeshExtent, setPlaneOffsets_planesSharingANormalAndPlanesWithOtherNormals)
{
  build_mesh_spread_over_procs();

  // Two planes along y (normals scaled differently, the same once unitized) and two along the diagonal,
  // interleaved, with one inside a CSG surface.
  auto bottomY = Plane::at_fraction_of_extent(stk::math::Vector3d(0,1,0), 0., 1.);
  auto cornerDiag = Plane::at_fraction_of_extent(stk::math::Vector3d(1,1,1), 1., 1.);
  auto topY = Plane::at_fraction_of_extent(stk::math::Vector3d(0,3,0), 1., 1.);
  auto centerDiag = Plane::at_fraction_of_extent(stk::math::Vector3d(2,2,2), 0.5, 1.);
  CSG_Surface csgWithCenterDiag;
  csgWithCenterDiag.set_tree(CSGTree(centerDiag));

  set_plane_offsets_from_mesh_extent(get_mesh(), coords_field(), {bottomY.get(), cornerDiag.get(), topY.get(), &csgWithCenterDiag});

  const double tol = 1.e-12;
  EXPECT_NEAR(0., bottomY->point_signed_distance(stk::math::Vector3d(0.2,-0.5,0.1)), tol);
  EXPECT_NEAR(0., topY->point_signed_distance(stk::math::Vector3d(0.2,0.5,0.1)), tol);
  EXPECT_NEAR(0., cornerDiag->point_signed_distance(stk::math::Vector3d(0.5,0.5,0.5)), tol);
  EXPECT_NEAR(0., centerDiag->point_signed_distance(stk::math::Vector3d(0.,0.,0.)), tol);
}

TEST_F(CubeOf12TetsPlaneMeshExtent, setPlaneOffsets_leavesOtherPlanesUnchanged)
{
  build_mesh_spread_over_procs();

  Plane plane(stk::math::Vector3d(0,0,1), -7.);
  set_plane_offsets_from_mesh_extent(get_mesh(), coords_field(), {&plane});
  EXPECT_DOUBLE_EQ(-7., plane.offset());
}
}
