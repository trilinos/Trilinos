// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <map>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMeshRespectingMeshBasedGeometry.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

namespace {

// A flat-walled BCC tet4 mesh of [-scale,scale]^3, each wall in its own sideset (Xminus, Xplus, ...), with a
// node inside the wall x = scale moved along the wall and a node on the edge x = y = scale moved along the
// edge, which makes poor tets beside them.
class PerturbedBoxForSmoothing
{
public:
  explicit PerturbedBoxForSmoothing(const double scale = 1.) : myScale(scale), myBBoxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD)
  {
    myBBoxMesh.set_mesh_structure_type(FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
    populate_bounding_box_mesh_and_activate(myBBoxMesh, -scale * stk::math::Vector3d(1., 1., 1.), scale * stk::math::Vector3d(1., 1., 1.), 0.5 * scale);
    move_node_at(scale * stk::math::Vector3d(1., 0., 0.), scale * stk::math::Vector3d(1., 0.2, 0.15));
    move_node_at(scale * stk::math::Vector3d(1., 1., 0.), scale * stk::math::Vector3d(1., 1., 0.2));
  }

  const stk::mesh::BulkData & mesh() const { return myBBoxMesh.bulk_data(); }
  CoordinatesFieldRef coords() const { return CoordinatesFieldRef(myBBoxMesh.meta_data().coordinate_field(), 3); }
  const stk::mesh::Selector activeSelector() const { return AuxMetaData::get(myBBoxMesh.meta_data()).active_part(); }
  double worst_mean_ratio() const { return compute_mesh_quality(mesh(), activeSelector(), coords(), MeanRatioQualityMetric()); }

  double scale() const { return myScale; }

  // The local nodes on a wall, with their locations.
  std::map<stk::mesh::Entity, stk::math::Vector3d> wall_node_locations() const
  {
    std::map<stk::mesh::Entity, stk::math::Vector3d> locations;
    for (auto * bucket : mesh().buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
      {
        const stk::math::Vector3d x = get_vector_field(mesh(), coords(), node, 3);
        if (!walls_of_location(x / myScale).empty()) locations[node] = x;
      }
    return locations;
  }

  static std::vector<unsigned> walls_of_location(const stk::math::Vector3d & x)
  {
    std::vector<unsigned> walls;
    for (unsigned axis = 0; axis < 3; ++axis)
      if (std::abs(std::abs(x[axis]) - 1.) < 1.e-12) walls.push_back(axis);
    return walls;
  }

private:
  // Moves the node at oldLocation on every processor that has it.
  void move_node_at(const stk::math::Vector3d & oldLocation, const stk::math::Vector3d & newLocation)
  {
    auto coordsData = coords().data<double, stk::mesh::ReadWrite>();
    for (auto * bucket : myBBoxMesh.bulk_data().buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
      {
        if ((get_vector_field(mesh(), coords(), node, 3) - oldLocation).length() > 1.e-12 * myScale) continue;
        auto nodeCoords = coordsData.entity_values(node);
        for (unsigned j = 0; j < 3; ++j)
          nodeCoords(stk::mesh::ComponentIdx(j)) = newLocation[j];
      }
  }

  double myScale;
  BoundingBoxMesh myBBoxMesh;
};

} // namespace

static void expect_nodes_stay_on_their_walls_corners_stay_and_quality_rises(PerturbedBoxForSmoothing & box)
{
  const double worstBefore = box.worst_mean_ratio();
  const auto locationsBefore = box.wall_node_locations();

  EXPECT_TRUE(smooth_mesh_respecting_mesh_based_geometry(box.mesh(), box.coords(), box.activeSelector()));

  EXPECT_LT(worstBefore, box.worst_mean_ratio());
  for (auto && [node, before] : locationsBefore)
  {
    const stk::math::Vector3d after = get_vector_field(box.mesh(), box.coords(), node, 3);
    const std::vector<unsigned> walls = PerturbedBoxForSmoothing::walls_of_location(before / box.scale());
    for (const unsigned axis : walls)
      EXPECT_NEAR(before[axis], after[axis], 1.e-12 * box.scale()) << "node " << box.mesh().identifier(node) << " left its wall";
    if (walls.size() == 3)
    {
      EXPECT_EQ(before, after) << "corner node " << box.mesh().identifier(node) << " moved";
    }
  }
}

TEST(SmoothMeshRespectingMeshBasedGeometry, perturbedNodesOnWallAndEdge_nodesStayOnTheirWallsCornersStayAndQualityRises)
{
  PerturbedBoxForSmoothing box;
  expect_nodes_stay_on_their_walls_corners_stay_and_quality_rises(box);
}

// The optimizer's tolerances are relative to the edge length: the same smoothing in a box a millionth the
// size (as a mesh in meters of a micron-scale microstructure).
TEST(SmoothMeshRespectingMeshBasedGeometry, perturbedNodesOnWallAndEdgeOfSmallBox_nodesStayOnTheirWallsCornersStayAndQualityRises)
{
  PerturbedBoxForSmoothing box(1.e-6);
  expect_nodes_stay_on_their_walls_corners_stay_and_quality_rises(box);
}

TEST(SmoothMeshRespectingMeshBasedGeometry, perturbedNodesOnWallAndEdge_betterThanKeepingTheBoundaryFixed)
{
  PerturbedBoxForSmoothing boundaryFixed;
  smooth_mesh_keeping_boundary_nodes_fixed(boundaryFixed.mesh(), boundaryFixed.coords(), boundaryFixed.activeSelector());

  PerturbedBoxForSmoothing respectingGeometry;
  smooth_mesh_respecting_mesh_based_geometry(respectingGeometry.mesh(), respectingGeometry.coords(), respectingGeometry.activeSelector());

  EXPECT_LT(boundaryFixed.worst_mean_ratio(), respectingGeometry.worst_mean_ratio());
}

} // namespace krino
