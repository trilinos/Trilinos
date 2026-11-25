// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <functional>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalMotionConstraint.hpp>
#include <Akri_SharpFeature.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

using Vec = stk::math::Vector3d;

// Planar level sets on a bounding box mesh; the walls are found from the domain side parts, or,
// with those ignored, from the sharp edges between the exposed boundary sides.
class NodalMotionConstraintOnBox : public ::testing::Test
{
protected:
  void build(const stk::topology topology, const double meshSize)
  {
    myMesh = std::make_unique<BoundingBoxMesh>(topology, MPI_COMM_WORLD);
    myDim = topology.dimension();
    myLevelSetField = AuxMetaData::get(myMesh->meta_data())
                          .register_field("LS",
                              FieldType::REAL,
                              stk::topology::NODE_RANK,
                              1,
                              1,
                              myMesh->meta_data().universal_part());
    const double zMin = (myDim == 3) ? -1. : 0.;
    const double zMax = (myDim == 3) ? 1. : 0.;
    populate_bounding_box_mesh_and_activate(*myMesh, {-1., -1., zMin}, {1., 1., zMax}, meshSize);
  }

  // Bends the top wall (z = 1) to z = 1 + a (x^2 + y^2), keeping the bottom fixed.
  void bend_top_wall(const double a)
  {
    stk::mesh::BulkData & mesh = myMesh->bulk_data();
    auto coordsData = coords_field().data<double, stk::mesh::ReadWrite>();
    for (auto * bucketPtr :
        mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().universal_part()))
      for (auto node : *bucketPtr)
      {
        auto x = coordsData.entity_values(node);
        const double bump = x(stk::mesh::ComponentIdx(0)) * x(stk::mesh::ComponentIdx(0)) +
            x(stk::mesh::ComponentIdx(1)) * x(stk::mesh::ComponentIdx(1));
        x(stk::mesh::ComponentIdx(2)) += a * 0.5 * (1. + x(stk::mesh::ComponentIdx(2))) * bump;
      }
  }

  void set_level_set(const std::function<double(const Vec &)> & levelSet)
  {
    stk::mesh::BulkData & mesh = myMesh->bulk_data();
    auto lsData = myLevelSetField.data<double, stk::mesh::ReadWrite>();
    for (auto * bucketPtr :
        mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().universal_part()))
      for (auto node : *bucketPtr)
      {
        const double phi = levelSet(get_vector_field(mesh, coords_field(), node, myDim));
        lsData.entity_values(node)() = (std::abs(phi) < 1.e-12) ? 0. : phi;
      }
  }

  void set_plane(const Vec & normal, const Vec & pointOnPlane)
  {
    const Vec n = normal.unit_vector();
    set_level_set([n, pointOnPlane](const Vec & x) { return Dot(n, x - pointOnPlane); });
  }

  FieldRef coords_field() const { return FieldRef(myMesh->meta_data().coordinate_field()); }
  stk::mesh::Part & exposed_part() const
  {
    return AuxMetaData::get(myMesh->meta_data()).exposed_boundary_part();
  }

  std::vector<const stk::mesh::Part *> wall_parts(const stk::mesh::Entity side) const
  {
    std::vector<const stk::mesh::Part *> parts;
    if (myUseWallParts)
      for (const char * name : {"Xminus", "Xplus", "Yminus", "Yplus", "Zminus", "Zplus"})
        if (stk::mesh::Part * part = myMesh->meta_data().get_part(name))
          if (myMesh->bulk_data().bucket(side).member(*part)) parts.push_back(part);
    return parts;
  }

  bool are_sides_in_same_wall(const stk::mesh::Entity node,
      const stk::mesh::Entity side0,
      const stk::mesh::Entity side1) const
  {
    if (wall_parts(side0) != wall_parts(side1)) return false;
    if (myUseWallParts) return true;
    const stk::mesh::BulkData & mesh = myMesh->bulk_data();
    const double cosFeatureAngle = std::cos(135. * M_PI / 180.);
    if (myDim == 2)
      return !is_sharp_between_line_sides(
          mesh, coords_field(), node, side0, side1, cosFeatureAngle);
    for (auto n0 : StkMeshEntities{mesh.begin_nodes(side0), mesh.end_nodes(side0)})
      for (auto n1 : StkMeshEntities{mesh.begin_nodes(side1), mesh.end_nodes(side1)})
        if (n0 != node && n0 == n1)
          return !is_sharp_between_tri_sides(
              mesh, coords_field(), {{node, n0}}, side0, side1, cosFeatureAngle);
    return true;
  }

  NodalMotionConstraint constraint(const stk::mesh::Entity node) const
  {
    return compute_zero_levelset_boundary_node_constraint(myMesh->bulk_data(),
        coords_field(),
        myLevelSetField,
        node,
        exposed_part(),
        myMesh->meta_data().universal_part(),
        [&](const stk::mesh::Entity side0, const stk::mesh::Entity side1)
        { return are_sides_in_same_wall(node, side0, side1); });
  }

  // The locally owned node at x, if any.
  stk::mesh::Entity owned_node_at(const Vec & x) const
  {
    const stk::mesh::BulkData & mesh = myMesh->bulk_data();
    for (auto * bucketPtr :
        mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part()))
      for (auto node : *bucketPtr)
        if ((get_vector_field(mesh, coords_field(), node, myDim) - x).length() < 1.e-9) return node;
    return stk::mesh::Entity();
  }

  // Checks the constraint of the node at x on the processor that owns it.
  void expect_constraint_at(
      const Vec & x, const NodalMotionConstraint::Type goldType, const Vec & goldDirection)
  {
    size_t numChecked = 0;
    const stk::mesh::Entity node = owned_node_at(x);
    if (myMesh->bulk_data().is_valid(node))
    {
      ++numChecked;
      const NodalMotionConstraint c = constraint(node);
      EXPECT_EQ(goldType, c.type) << x;
      if (goldType == NodalMotionConstraint::SURFACE)
      {
        EXPECT_NEAR(1., Dot(c.direction, goldDirection.unit_vector()), 1.e-12) << x;
      }
      else if (goldType == NodalMotionConstraint::LINE)
      {
        EXPECT_NEAR(1., std::abs(Dot(c.direction, goldDirection.unit_vector())), 1.e-12) << x;
      }
    }
    all_reduce_sum(myMesh->bulk_data().parallel(), numChecked);
    EXPECT_EQ(1u, numChecked) << x;
  }

  void expect_plane_cases()
  {
    // a face node, the interface crossing its wall: within the wall
    set_plane({1., 2., 4.}, {0., 0., 0.});
    expect_constraint_at({1., 0.5, -0.5}, NodalMotionConstraint::SURFACE, {1., 0., 0.});
    // an interior node: free
    expect_constraint_at({0., 0., 0.}, NodalMotionConstraint::FREE, Vec::ZERO);
    // a box edge node, the interface crossing the edge: along the edge
    set_plane({1., 1., 1.}, {1., 1., 0.5});
    expect_constraint_at({1., 1., 0.5}, NodalMotionConstraint::LINE, {0., 0., 1.});
    // an isolated corner: along the box edge whose neighbor's sign differs from the other two
    set_plane({1., 1., -1.}, {1., 1., 1.});
    expect_constraint_at({1., 1., 1.}, NodalMotionConstraint::LINE, {0., 0., 1.});
    set_plane({3., -1., -1.}, {1., 1., 1.});
    expect_constraint_at({1., 1., 1.}, NodalMotionConstraint::LINE, {1., 0., 0.});
    // a box edge in the interface: within the wall on the negative side (y = 1), and at the corners
    // at its ends, along the box edge where that wall meets the third
    set_plane({1., -2., 0.}, {1., 1., 0.});
    expect_constraint_at({1., 1., 0.}, NodalMotionConstraint::SURFACE, {0., 1., 0.});
    expect_constraint_at({1., 1., 1.}, NodalMotionConstraint::LINE, {1., 0., 0.});
  }

  std::unique_ptr<BoundingBoxMesh> myMesh;
  FieldRef myLevelSetField;
  unsigned myDim = 3;
  bool myUseWallParts = true;
};

} // namespace

TEST_F(NodalMotionConstraintOnBox, tetsWithWallParts_followsTheInterfaceContact)
{
  build(stk::topology::TETRAHEDRON_4, 0.5);
  expect_plane_cases();
}

TEST_F(
    NodalMotionConstraintOnBox, tetsWithSharpFeaturesInsteadOfWallParts_followsTheInterfaceContact)
{
  build(stk::topology::TETRAHEDRON_4, 0.5);
  myUseWallParts = false;
  expect_plane_cases();
}

// When every boundary side is in one wall, a box edge node moves within it, normal to the mean of
// the two faces.
TEST_F(NodalMotionConstraintOnBox, tetsWithOneWall_boxEdgeNodeOnTheMeanWall)
{
  build(stk::topology::TETRAHEDRON_4, 0.5);
  const stk::mesh::BulkData & mesh = myMesh->bulk_data();
  set_plane({1., 1., 1.}, {1., 1., 0.5});
  const stk::mesh::Entity node = owned_node_at({1., 1., 0.5});
  if (!mesh.is_valid(node)) return;
  const NodalMotionConstraint c = compute_zero_levelset_boundary_node_constraint(mesh,
      coords_field(),
      myLevelSetField,
      node,
      exposed_part(),
      mesh.mesh_meta_data().universal_part(),
      [](const stk::mesh::Entity, const stk::mesh::Entity) { return true; });
  EXPECT_EQ(NodalMotionConstraint::SURFACE, c.type);
  EXPECT_GT(c.direction[0], 0.);
  EXPECT_GT(c.direction[1], 0.);
}

// On a curved top wall, a node away from the box edges moves within the wall, normal to it to the
// accuracy of the sides' mean normal.
TEST_F(NodalMotionConstraintOnBox, tetsWithCurvedTopWall_surfaceNearTheAnalyticNormal)
{
  const double a = 0.3;
  build(stk::topology::TETRAHEDRON_4, 1. / 3.);
  bend_top_wall(a);
  myUseWallParts = false;
  set_level_set([](const Vec & x) { return x[0] - x[1]; });
  const Vec x(1. / 3., 1. / 3., 1. + a * 2. / 9.);
  const stk::mesh::Entity node = owned_node_at(x);
  size_t numChecked = 0;
  if (myMesh->bulk_data().is_valid(node))
  {
    ++numChecked;
    const NodalMotionConstraint c = constraint(node);
    ASSERT_EQ(NodalMotionConstraint::SURFACE, c.type);
    const Vec analyticNormal = Vec(-2. * a * x[0], -2. * a * x[1], 1.).unit_vector();
    EXPECT_GT(Dot(c.direction, analyticNormal), 0.999);
  }
  all_reduce_sum(myMesh->bulk_data().parallel(), numChecked);
  EXPECT_EQ(1u, numChecked);
}

TEST_F(NodalMotionConstraintOnBox, trisWithWallParts_sideAndCornerNodes)
{
  build(stk::topology::TRIANGLE_3_2D, 0.25);
  // a side node: along its wall
  const Vec sideNode(1., -1. / 9., 0.); // the lattice rows are at y = -1 + 2k/9
  set_plane({1., 2., 0.}, sideNode);
  expect_constraint_at(sideNode, NodalMotionConstraint::SURFACE, {1., 0., 0.});
  // a corner: along the wall with the negative neighbor (y = 1)
  set_plane({1., -2., 0.}, {1., 1., 0.});
  expect_constraint_at({1., 1., 0.}, NodalMotionConstraint::SURFACE, {0., 1., 0.});
}

TEST(NodalMotionConstraint, constrainKeepingNormalSpeed)
{
  const Vec v(1., 2., 3.);
  const NodalMotionConstraint surface{NodalMotionConstraint::SURFACE, Vec(0., 0., 1.)};
  const Vec onSurface = constrain_keeping_normal_speed(surface, v);
  EXPECT_NEAR(0., onSurface[2], 1.e-14);
  EXPECT_NEAR(Dot(v, v), Dot(v, onSurface), 1.e-12);

  const NodalMotionConstraint line{NodalMotionConstraint::LINE, Vec(0., 1., 0.)};
  const Vec onLine = constrain_keeping_normal_speed(line, v);
  EXPECT_NEAR(0., onLine[0], 1.e-14);
  EXPECT_NEAR(0., onLine[2], 1.e-14);
  EXPECT_NEAR(Dot(v, v), Dot(v, onLine), 1.e-12);

  EXPECT_EQ(v, constrain_keeping_normal_speed({NodalMotionConstraint::FREE, Vec::ZERO}, v));

  // The constrained part would be orthogonal to v.
  const Vec normalToSurface(0., 0., 2.);
  EXPECT_EQ(normalToSurface, constrain_keeping_normal_speed(surface, normalToSurface));
}

} // namespace krino
