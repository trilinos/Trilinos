// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <type_traits>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetsFromSides.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

// The interface velocity, linear in the coordinates, at states N and NP1.
stk::math::Vector3d velocity_n(const stk::math::Vector3d & x)
{
  return stk::math::Vector3d(1., x[1], -x[0]);
}

stk::math::Vector3d velocity_np1(const stk::math::Vector3d & x)
{
  return stk::math::Vector3d(x[1] + 1., 2. * x[0], x[2] - x[0]);
}

// Two elements in blocks 1 and 2, whose block boundary side is the interface.
template <stk::topology::topology_t TOPO, typename MESH_SPEC>
class InterfaceSideFacets : public StkMeshFixture<TOPO>
{
protected:
  static constexpr unsigned DIM = StkMeshAndBuilder<TOPO>::NUM_DIM;

  // The interface of the two elements is the side x=0 between block_1 (x>0) and block_2 (x<0).
  void test_facet_mesh_from_interface_sides()
  {
    this->set_valid_proc_sizes_for_test({1, 2});
    if (!this->is_valid_proc_size_for_test()) return;

    stk::mesh::BulkData & mesh = this->mMesh;
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    this->build_mesh_from_spec_in_blocks(MESH_SPEC(), {1, 2}, {{2, {0, 1}}});

    const stk::mesh::Selector interfaceSelector = this->get_aux_meta().block_boundary_part();
    const stk::mesh::Selector negativeSideBlockSelector = *meta.get_part("block_1");
    const FacetMesh facets = build_facet_mesh_from_interface_sides(mesh,
        {{interfaceSelector, negativeSideBlockSelector}},
        meta.universal_part(),
        this->get_coordinates_field());

    const unsigned numOwnedSides = stk::mesh::count_selected_entities(
        interfaceSelector & meta.locally_owned_part(), mesh.buckets(meta.side_rank()));
    EXPECT_EQ(numOwnedSides, facets.num_owned_facets());
    EXPECT_EQ(DIM * numOwnedSides, facets.num_owned_facet_vertices());
    for (unsigned f = 0; f < facets.num_owned_facets(); ++f)
    {
      for (unsigned v : facets.facet_vertices(f))
      {
        const stk::mesh::Entity node =
            mesh.get_entity(stk::topology::NODE_RANK, facets.vertex_id(v));
        ASSERT_TRUE(mesh.is_valid(node));
        EXPECT_EQ(this->get_node_coordinates(node), facets.vertex_coordinates(v));
      }
      using FACET = std::conditional_t<DIM == 3, Facet3d, Facet2d>;
      const stk::math::Vector3d normal = facets.make_facet<FACET>(f).facet_normal();
      EXPECT_NEAR(0., (normal - stk::math::Vector3d(-1., 0., 0.)).length(), 1.e-12);
    }

    // A side of two surfaces is two facets, with the same vertices.
    const FacetMesh twice = build_facet_mesh_from_interface_sides(mesh,
        {{interfaceSelector, negativeSideBlockSelector},
            {interfaceSelector, negativeSideBlockSelector}},
        meta.universal_part(),
        this->get_coordinates_field());
    EXPECT_EQ(2 * numOwnedSides, twice.num_owned_facets());
    EXPECT_EQ(facets.num_owned_facet_vertices(), twice.num_owned_facet_vertices());
  }

  // With a two-state nodal interface velocity of velocity_n and velocity_np1, the surface built
  // from the interface sides gives, at the closest point to any node, the average of N and NP1.
  void test_velocity_at_closest_point_and_max_velocity()
  {
    this->set_valid_proc_sizes_for_test({1, 2});
    if (!this->is_valid_proc_size_for_test()) return;

    stk::mesh::BulkData & mesh = this->mMesh;
    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    auto & velocityField =
        meta.declare_field<double>(stk::topology::NODE_RANK, "interface_velocity", 2);
    stk::mesh::put_field_on_mesh(velocityField, meta.universal_part(), DIM, nullptr);
    this->build_mesh_from_spec_in_blocks(MESH_SPEC(), {1, 2}, {{2, {0, 1}}});

    const FieldRef velocity(velocityField);
    for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    {
      for (auto node : *bucketPtr)
      {
        const stk::math::Vector3d x = this->get_node_coordinates(node);
        double * np1 = field_data<double>(velocity.field(), node);
        double * n = field_data<double>(velocity.field_state(stk::mesh::StateOld).field(), node);
        for (unsigned d = 0; d < DIM; ++d)
        {
          np1[d] = velocity_np1(x)[d];
          n[d] = velocity_n(x)[d];
        }
      }
    }
    auto expected_velocity = [](const stk::math::Vector3d & x)
    {
      stk::math::Vector3d v = 0.5 * (velocity_n(x) + velocity_np1(x));
      if (DIM == 2) v[2] = 0.;
      return v;
    };

    // The fixture's sides are not in the active part, so every side counts as active here.
    const stk::mesh::Part & activePart = meta.universal_part();
    const stk::mesh::Selector interfaceSelector = this->get_aux_meta().block_boundary_part();
    const stk::mesh::Selector negativeSideBlockSelector = *meta.get_part("block_1");
    const std::unique_ptr<MeshFacetedSurface> surface =
        build_mesh_faceted_surface_with_velocity_from_interface_sides(mesh,
            interfaceSelector,
            negativeSideBlockSelector,
            activePart,
            this->get_coordinates_field(),
            velocity,
            2);
    BoundingBox nodeBbox;
    for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
      for (auto node : *bucketPtr)
        nodeBbox.accommodate(this->get_node_coordinates(node));
    surface->prepare_to_compute(0., nodeBbox, 0.);

    for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    {
      for (auto node : *bucketPtr)
      {
        const stk::math::Vector3d x = this->get_node_coordinates(node);
        const stk::math::Vector3d closestPt = surface->closest_point(x);
        EXPECT_NEAR(0.,
            (expected_velocity(closestPt) - surface->velocity_at_closest_point(x)).length(),
            1.e-12)
            << x;
      }
    }

    // The largest velocity at the nodes of this processor's owned interface sides.
    double expectedMax = 0.;
    for (auto * bucketPtr : mesh.get_buckets(
             meta.side_rank(), interfaceSelector & activePart & meta.locally_owned_part()))
      for (auto side : *bucketPtr)
        for (unsigned n = 0; n < mesh.num_nodes(side); ++n)
          expectedMax = std::max(expectedMax,
              expected_velocity(this->get_node_coordinates(mesh.begin_nodes(side)[n])).length());
    EXPECT_NEAR(expectedMax,
        compute_max_interface_velocity_magnitude(
            mesh, interfaceSelector, negativeSideBlockSelector, activePart, velocity, 2),
        1.e-12);
  }
};

using InterfaceSideFacetsTwoTets = InterfaceSideFacets<stk::topology::TETRAHEDRON_4, TwoRightTets>;
using InterfaceSideFacetsTwoTris = InterfaceSideFacets<stk::topology::TRIANGLE_3_2D, TwoRightTris>;

} // namespace

TEST_F(InterfaceSideFacetsTwoTets, buildFacetMesh_nodeIdsCoordinatesAndOrientation)
{
  test_facet_mesh_from_interface_sides();
}

TEST_F(InterfaceSideFacetsTwoTris, buildFacetMesh_nodeIdsCoordinatesAndOrientation)
{
  test_facet_mesh_from_interface_sides();
}

TEST_F(
    InterfaceSideFacetsTwoTets, velocityAtClosestPoint_averageOfNAndNp1AndMaxAtOwnedInterfaceNodes)
{
  test_velocity_at_closest_point_and_max_velocity();
}

TEST_F(
    InterfaceSideFacetsTwoTris, velocityAtClosestPoint_averageOfNAndNp1AndMaxAtOwnedInterfaceNodes)
{
  test_velocity_at_closest_point_and_max_velocity();
}

} // namespace krino
