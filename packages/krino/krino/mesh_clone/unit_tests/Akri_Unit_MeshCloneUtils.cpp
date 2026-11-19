// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <memory>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_MeshCloneUtils.hpp>
#include <Akri_Unit_MeshCloneTestUtils.hpp>

namespace krino
{

class MeshCloneUtils : public StkMeshFixture<stk::topology::TETRAHEDRON_4>
{
};

TEST_F(MeshCloneUtils, cloneMetaDataAndBulkData_equalMeshAndSharingByIdsRestored)
{
  if (!build_two_block_tet_mesh(*this)) return;

  const std::unique_ptr<stk::mesh::BulkData> clone = clone_mesh(mMesh);

  EXPECT_FALSE(clone->use_entity_ids_for_resolving_sharing());
  expect_equal_meta_data_parts_and_fields(mMesh.mesh_meta_data(), clone->mesh_meta_data());
  expect_equal_meshes(mMesh, *clone);
}

TEST_F(MeshCloneUtils, cloneWithAddedPartAndField_equalMeshWithAddedPartAndField)
{
  if (!build_two_block_tet_mesh(*this)) return;
  const stk::mesh::MetaData & meta = mMesh.mesh_meta_data();

  std::shared_ptr<stk::mesh::MetaData> cloneMeta = stk::mesh::MeshBuilder(mComm).create_meta_data();
  clone_meta_data_parts_and_fields(meta, *cloneMeta);
  stk::mesh::Part & addedPart = cloneMeta->declare_part("added_part", stk::topology::NODE_RANK);
  auto & addedField = cloneMeta->declare_field<double>(stk::topology::NODE_RANK, "added_field");
  stk::mesh::put_field_on_mesh(addedField, cloneMeta->universal_part(), nullptr);
  std::unique_ptr<stk::mesh::BulkData> clone = stk::mesh::MeshBuilder(mComm).create(cloneMeta);
  cloneMeta->commit();
  clone_bulk_data(mMesh, *clone);

  expect_equal_meta_data_parts_and_fields(meta, *cloneMeta);
  expect_equal_meshes(mMesh, *clone);
  EXPECT_LE(meta.get_parts().size(), addedPart.mesh_meta_data_ordinal());
  EXPECT_LE(meta.get_fields().size(), addedField.mesh_meta_data_ordinal());
}

TEST(MeshCloneUtilsIo,
    cloneIoAttributesAndSurfaceToBlockMapping_rolesOutputTypesDistributionFactorsAndMapCopied)
{
  std::shared_ptr<stk::mesh::MetaData> meta =
      stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create_meta_data();
  stk::mesh::Part & block =
      meta->declare_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  stk::mesh::Part & surface =
      meta->declare_part_with_topology("surface_1", stk::topology::TRIANGLE_3);
  stk::io::put_io_part_attribute(block);
  stk::io::put_io_part_attribute(surface);
  meta->set_surface_to_block_mapping(&surface, {&block});

  auto & coords = meta->declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh(coords, meta->universal_part(), 3, nullptr);
  stk::io::set_field_role(coords, Ioss::Field::MESH);
  stk::io::set_field_output_type(coords, stk::io::FieldOutputType::VECTOR_3D);
  auto & distFact = meta->declare_field<double>(stk::topology::NODE_RANK, "distribution_factors");
  stk::mesh::put_field_on_mesh(distFact, surface, nullptr);
  stk::io::set_distribution_factor_field(surface, distFact);

  stk::mesh::MetaData clone;
  clone_meta_data_parts_and_fields(*meta, clone);
  clone_io_attributes_and_surface_to_block_mapping(*meta, clone);

  expect_equal_meta_data_parts_and_fields(*meta, clone);
  expect_equal_io_attributes_and_surface_to_block_mapping(*meta, clone);
  ASSERT_EQ(1u, clone.get_surfaces_in_surface_to_block_map().size());
}

} // namespace krino
