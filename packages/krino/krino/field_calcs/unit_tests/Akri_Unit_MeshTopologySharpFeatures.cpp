// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshTopologySharpFeatures.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

namespace {

// A flat-walled generated mesh of [-1,1]^dim, each wall in its own sideset (Xminus, Xplus, ...) unless
// requested otherwise, and the axes along which a node is on a wall.
class MeshTopologySharpFeaturesOfBox
{
public:
  // Without wall sidesets, the walls' sides are only in the exposed-boundary part: every wall is in one surface.
  MeshTopologySharpFeaturesOfBox(const stk::topology elemTopology, const BoundingBoxMeshStructureType structure, const bool withWallSidesets = true)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD)
  {
    myBBoxMesh.set_mesh_structure_type(structure);
    myBBoxMesh.set_create_domain_sidesets(withWallSidesets);
    const unsigned dim = elemTopology.dimension();
    const stk::math::Vector3d corner(1., 1., (dim == 3) ? 1. : 0.);
    populate_bounding_box_mesh_and_activate(myBBoxMesh, -corner, corner, 0.5);
  }

  const stk::mesh::BulkData & mesh() const { return myBBoxMesh.bulk_data(); }
  FieldRef coords() const { return FieldRef(myBBoxMesh.meta_data().coordinate_field()); }
  std::vector<SmoothingDirection> smoothing_directions(const std::optional<double> cosFeatureAngle = std::nullopt) const
  {
    return compute_smoothing_directions_from_mesh_topology(
        mesh(), coords(), AuxMetaData::get(myBBoxMesh.meta_data()).active_part(), cosFeatureAngle);
  }

  std::vector<unsigned> wall_axes(const stk::mesh::Entity node) const
  {
    const unsigned dim = mesh().mesh_meta_data().spatial_dimension();
    const stk::math::Vector3d x = get_vector_field(mesh(), coords(), node, dim);
    std::vector<unsigned> axes;
    for (unsigned axis = 0; axis < dim; ++axis)
      if (std::abs(std::abs(x[axis]) - 1.) < 1.e-12) axes.push_back(axis);
    return axes;
  }

  size_t num_owned_or_shared_nodes_on_a_wall() const
  {
    const stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    size_t count = 0;
    for (auto * bucket : mesh().get_buckets(stk::topology::NODE_RANK, meta.locally_owned_part() | meta.globally_shared_part()))
      for (auto node : *bucket)
        if (!wall_axes(node).empty()) ++count;
    return count;
  }

private:
  BoundingBoxMesh myBBoxMesh;
};

stk::math::Vector3d unit_axis(const unsigned axis)
{
  stk::math::Vector3d v = stk::math::Vector3d::ZERO;
  v[axis] = 1.;
  return v;
}

void expect_parallel_to(const stk::math::Vector3d & expected, const stk::math::Vector3d & direction)
{
  EXPECT_NEAR(1., std::abs(Dot(expected, direction)), 1.e-12) << "direction " << direction << ", expected along " << expected;
}

void expect_sorted_by_id(const stk::mesh::BulkData & mesh, const std::vector<SmoothingDirection> & directions)
{
  EXPECT_TRUE(std::is_sorted(directions.begin(), directions.end(),
      [&mesh](const SmoothingDirection & a, const SmoothingDirection & b) { return mesh.identifier(a.node) < mesh.identifier(b.node); }));
}

// Every wall node of the 3D box: within the plane of its wall, along the edge of two walls, or fixed at a
// corner; no interior node.
void expect_box_walls_edges_and_corners_3D(const MeshTopologySharpFeaturesOfBox & box, const std::vector<SmoothingDirection> & directions)
{
  EXPECT_EQ(box.num_owned_or_shared_nodes_on_a_wall(), directions.size());
  expect_sorted_by_id(box.mesh(), directions);
  for (auto && direction : directions)
  {
    const std::vector<unsigned> axes = box.wall_axes(direction.node);
    if (axes.size() == 1)
    {
      EXPECT_EQ(SmoothingDirection::SMOOTH_IN_PLANE, direction.type);
      expect_parallel_to(unit_axis(axes[0]), direction.direction);
    }
    else if (axes.size() == 2)
    {
      EXPECT_EQ(SmoothingDirection::SMOOTH_ALONG_TANGENT, direction.type);
      expect_parallel_to(unit_axis(3 - axes[0] - axes[1]), direction.direction);
    }
    else
    {
      EXPECT_EQ(3u, axes.size());
      EXPECT_EQ(SmoothingDirection::FIXED, direction.type);
      EXPECT_EQ(stk::math::Vector3d::ZERO, direction.direction);
    }
  }
}

// Every wall node of the 2D box: along its wall, or fixed at a corner; no interior node.
void expect_box_walls_and_corners_2D(const MeshTopologySharpFeaturesOfBox & box, const std::vector<SmoothingDirection> & directions)
{
  EXPECT_EQ(box.num_owned_or_shared_nodes_on_a_wall(), directions.size());
  expect_sorted_by_id(box.mesh(), directions);
  for (auto && direction : directions)
  {
    const std::vector<unsigned> axes = box.wall_axes(direction.node);
    if (axes.size() == 1)
    {
      EXPECT_EQ(SmoothingDirection::SMOOTH_ALONG_TANGENT, direction.type);
      expect_parallel_to(unit_axis(1 - axes[0]), direction.direction);
    }
    else
    {
      EXPECT_EQ(2u, axes.size());
      EXPECT_EQ(SmoothingDirection::FIXED, direction.type);
      EXPECT_EQ(stk::math::Vector3d::ZERO, direction.direction);
    }
  }
}

// The box's edges and corners are 90 degrees: sharp for a feature angle of 135; its walls are flat.
const double cosOf135Degrees = std::cos(135. * M_PI / 180.);

} // namespace

TEST(MeshTopologySharpFeatures, box3D_wallNodesInPlaneEdgeNodesAlongTangentCornersFixedInteriorNodesAbsent)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TETRAHEDRON_4, FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  expect_box_walls_edges_and_corners_3D(box, box.smoothing_directions());
}

TEST(MeshTopologySharpFeatures, box2D_sideNodesAlongTangentCornersFixedInteriorNodesAbsent)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TRIANGLE_3_2D, FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH);
  expect_box_walls_and_corners_2D(box, box.smoothing_directions());
}

TEST(MeshTopologySharpFeatures, box3DWallsInOneSurface_withoutFeatureAngle_everyWallNodeInPlane)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TETRAHEDRON_4, FLAT_WALLED_BCC_BOUNDING_BOX_MESH, false);
  const std::vector<SmoothingDirection> directions = box.smoothing_directions();

  EXPECT_EQ(box.num_owned_or_shared_nodes_on_a_wall(), directions.size());
  for (auto && direction : directions)
    EXPECT_EQ(SmoothingDirection::SMOOTH_IN_PLANE, direction.type);
}

TEST(MeshTopologySharpFeatures, box3DWallsInOneSurface_withFeatureAngle_edgesAndCornersAsWithASidesetPerWall)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TETRAHEDRON_4, FLAT_WALLED_BCC_BOUNDING_BOX_MESH, false);
  expect_box_walls_edges_and_corners_3D(box, box.smoothing_directions(cosOf135Degrees));
}

TEST(MeshTopologySharpFeatures, box2DWallsInOneSurface_withoutFeatureAngle_everyWallNodeAlongTangent)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TRIANGLE_3_2D, FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH, false);
  const std::vector<SmoothingDirection> directions = box.smoothing_directions();

  EXPECT_EQ(box.num_owned_or_shared_nodes_on_a_wall(), directions.size());
  for (auto && direction : directions)
    EXPECT_EQ(SmoothingDirection::SMOOTH_ALONG_TANGENT, direction.type);
}

TEST(MeshTopologySharpFeatures, box2DWallsInOneSurface_withFeatureAngle_cornersAsWithASidesetPerWall)
{
  const MeshTopologySharpFeaturesOfBox box(stk::topology::TRIANGLE_3_2D, FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH, false);
  expect_box_walls_and_corners_2D(box, box.smoothing_directions(cosOf135Degrees));
}

} // namespace krino
