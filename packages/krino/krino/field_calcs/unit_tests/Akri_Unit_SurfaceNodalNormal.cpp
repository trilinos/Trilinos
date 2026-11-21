// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_SurfaceNodalNormal.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <string>

namespace krino
{

// FourRightTris in blocks {1,2,2,1}: the block boundary is the edges 1-3 (between elements 0 and 1)
// and 3-2 (between elements 1 and 2), so only edge 1-3 separates the blocks; edge 3-2 is interior
// to block 2.  Built without aura to show the sum doesn't need it.
class SurfaceNodalNormal : public StkMeshTriFixture
{
public:
  SurfaceNodalNormal()
  {
    set_valid_proc_sizes_for_test({1, 2, 3, 4});
    mMesh.set_automatic_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    mBuilder.create_sideset_part(surfaceId);
    mNormal = register_surface_nodal_normal_field(mMesh.mesh_meta_data(), "normal", surface_part());
    mElement1Part = &mMesh.mesh_meta_data().declare_part("element_1", stk::topology::ELEMENT_RANK);
    build_mesh_from_spec_in_blocks(meshSpec, {1, 2, 2, 1});
    if (is_valid_proc_size_for_test())
    {
      create_surface_sides();
      std::vector<stk::mesh::Entity> elems;
      if (is_owned_element(1)) elems.push_back(get_assigned_element_for_index(1));
      mMesh.batch_change_entity_parts(elems,
          std::vector<stk::mesh::PartVector>(elems.size(), {mElement1Part}),
          std::vector<stk::mesh::PartVector>(elems.size()));
    }
  }

protected:
  static constexpr unsigned surfaceId = 1;
  stk::mesh::Part & surface_part() const
  {
    return *mMesh.mesh_meta_data().get_part("surface_" + std::to_string(surfaceId));
  }
  stk::mesh::Part & block_part(const unsigned blockId) const
  {
    return *mMesh.mesh_meta_data().get_part("block_" + std::to_string(blockId));
  }

  void update(const stk::mesh::Selector & elementSelector)
  {
    update_surface_nodal_normal(
        mMesh, mMesh.mesh_meta_data().coordinate_field(), mNormal, surface_part(), elementSelector);
  }

  // The interior sides 1-3 (side 1 of element 0) and 3-2 (side 1 of element 1).  They are created,
  // if they don't exist yet, by the owner of that element (the element graph attaches the element
  // on the other side), then the owner of each side adds it to the surface part.
  void create_surface_sides()
  {
    std::vector<SideDescription> sideRequests;
    for (unsigned elemIndex : {0u, 1u})
      if (is_owned_element(elemIndex))
        sideRequests.emplace_back(
            get_assigned_element_for_index(elemIndex), 1, stk::mesh::PartVector{});
    batch_create_sides(mMesh, sideRequests);

    std::vector<stk::mesh::Entity> surfaceSides;
    for (const std::array<unsigned, 2> sideNodeIndices :
        {std::array<unsigned, 2>{{1, 3}}, std::array<unsigned, 2>{{3, 2}}})
    {
      if (!has_node(sideNodeIndices[0]) || !has_node(sideNodeIndices[1])) continue;
      const stk::mesh::Entity side =
          mBuilder.get_side_with_nodes({get_assigned_node_for_index(sideNodeIndices[0]),
              get_assigned_node_for_index(sideNodeIndices[1])});
      if (mMesh.is_valid(side) && mMesh.bucket(side).owned()) surfaceSides.push_back(side);
    }
    mMesh.batch_change_entity_parts(surfaceSides,
        std::vector<stk::mesh::PartVector>(surfaceSides.size(), {&surface_part()}),
        std::vector<stk::mesh::PartVector>(surfaceSides.size()));
  }

  bool has_node(const unsigned nodeIndex) const { return is_node_present_on_proc(nodeIndex); }
  stk::math::Vector3d normal_at(const unsigned nodeIndex) const
  {
    return get_vector_field(mMesh, mNormal, get_assigned_node_for_index(nodeIndex), 2);
  }
  void expect_normals(const std::array<stk::math::Vector3d, 3> & goldAtNodes1To3) const
  {
    for (unsigned n : {1u, 2u, 3u})
      if (has_node(n)) expect_near_absolute(goldAtNodes1To3[n - 1], normal_at(n), 1.e-14);
  }

  FourRightTris meshSpec;
  FieldRef mNormal;
  stk::mesh::Part * mElement1Part{nullptr};
};

TEST_F(SurfaceNodalNormal, block1_outOfElement0AcrossEdge13)
{
  if (!is_valid_proc_size_for_test()) return;
  update(block_part(1));

  // Edge 1-3 has length sqrt 2 and outward normal (1,1)/sqrt 2 from element 0.
  const stk::math::Vector3d areaVec13(1., 1., 0.);
  expect_normals({areaVec13, stk::math::Vector3d::ZERO, areaVec13});

  for (unsigned n : {0u, 4u, 5u})
  {
    if (has_node(n))
    {
      EXPECT_FALSE(has_field_data(mNormal, get_assigned_node_for_index(n)));
    }
  }
}

TEST_F(SurfaceNodalNormal, block2_oppositeOrientationAndInteriorSideCancels)
{
  if (!is_valid_proc_size_for_test()) return;
  update(block_part(2));

  // Edge 1-3 out of element 1 is -(1,1); edge 3-2 is attached to elements 1 and 2, both in block 2,
  // whose outward normals (0,1) and (0,-1) cancel.
  const stk::math::Vector3d areaVec13(-1., -1., 0.);
  expect_normals({areaVec13, stk::math::Vector3d::ZERO, areaVec13});
}

TEST_F(SurfaceNodalNormal, oneElementOnTwoSides_sumsTheSidesAtTheSharedNode)
{
  if (!is_valid_proc_size_for_test()) return;
  update(*mElement1Part);

  // Out of element 1 (nodes 1, 2, 3): edge 1-3 is -(1,1), edge 3-2 (length 1) is (0,1).
  const stk::math::Vector3d areaVec13(-1., -1., 0.), areaVec32(0., 1., 0.);
  expect_normals({areaVec13, areaVec32, areaVec13 + areaVec32});
}

TEST_F(SurfaceNodalNormal, allElements_bothOrientationsOfEachSideCancel)
{
  if (!is_valid_proc_size_for_test()) return;
  update(mMesh.mesh_meta_data().universal_part());
  expect_normals({stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO});
}

} // namespace krino
