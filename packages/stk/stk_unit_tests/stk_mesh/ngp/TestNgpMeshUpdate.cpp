#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpAtomics.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "NgpUnitTestUtils.hpp"

namespace {

class UpdateNgpMesh : public stk::unit_test_util::MeshFixture
{
public:
  void setup_test_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(UpdateNgpMesh, lazyAutoUpdate)
{
  setup_test_mesh();

  // Don't store persistent pointers/references if you want automatic updates
  // when acquiring an NgpMesh from BulkData
  stk::mesh::NgpMesh * ngpMesh = &get_bulk().get_updated_ngp_mesh();

  get_bulk().modification_begin();
  get_bulk().modification_end();

#ifdef STK_USE_DEVICE_MESH
  EXPECT_FALSE(ngpMesh->is_up_to_date());
  ngpMesh = &get_bulk().get_updated_ngp_mesh();  // Trigger update
  EXPECT_TRUE(ngpMesh->is_up_to_date());
#else
  EXPECT_TRUE(ngpMesh->is_up_to_date());
  ngpMesh = &get_bulk().get_updated_ngp_mesh();
  EXPECT_TRUE(ngpMesh->is_up_to_date());
#endif
}

TEST_F(UpdateNgpMesh, manualUpdate)
{
  setup_test_mesh();

  // If storing a persistent reference, call the update_ngp_mesh() method
  // to ensure that it is synchronized with BulkData
  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();

  get_bulk().modification_begin();
  get_bulk().modification_end();

#ifdef STK_USE_DEVICE_MESH
  EXPECT_FALSE(ngpMesh.is_up_to_date());
  get_bulk().update_ngp_mesh();
  EXPECT_TRUE(ngpMesh.is_up_to_date());
#else
  EXPECT_TRUE(ngpMesh.is_up_to_date());
  get_bulk().update_ngp_mesh();
  EXPECT_TRUE(ngpMesh.is_up_to_date());
#endif
}

TEST_F(UpdateNgpMesh, OnlyOneDeviceMesh_InternalAndExternal)
{
  setup_test_mesh();

  // Create first NgpMesh inside BulkData
  get_bulk().get_updated_ngp_mesh();

#ifdef STK_USE_DEVICE_MESH
  EXPECT_THROW(stk::mesh::NgpMesh secondNgpMesh(get_bulk()), std::logic_error);
#else
  EXPECT_NO_THROW(stk::mesh::NgpMesh secondNgpMesh(get_bulk()));
#endif
}

TEST_F(UpdateNgpMesh, OnlyOneDeviceMesh_TwoExternal)
{
  setup_test_mesh();

  stk::mesh::NgpMesh firstNgpMesh(get_bulk());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_THROW(stk::mesh::NgpMesh secondNgpMesh(get_bulk()), std::logic_error);
#else
  EXPECT_NO_THROW(stk::mesh::NgpMesh secondNgpMesh(get_bulk()));
#endif
}

class BucketLayoutModification : public stk::unit_test_util::MeshFixture
{
public:
  void setup_mesh_3hex_3block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_3hex_3block(get_bulk(), bucketCapacity);
  }

  void setup_mesh_3hex_2block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_3hex_2block(get_bulk(), bucketCapacity);
  }
};


//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, DeleteBucketInMiddle)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 2;
  setup_mesh_3hex_3block(bucketCapacity);

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_2|block_1|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, AddBucketInMiddle)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 1;
  setup_mesh_3hex_3block(bucketCapacity);

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_1|block_3|        |block_1|block_3|block_3|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, ChangeBucketContents)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 2;
  setup_mesh_3hex_2block(bucketCapacity);

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_3", {2,3}}});
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, DeleteBucketInMiddle_WithCopy)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 2;
  setup_mesh_3hex_3block(bucketCapacity);

  stk::mesh::NgpMesh ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_2|block_1|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, AddBucketInMiddle_WithCopy)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 1;
  setup_mesh_3hex_3block(bucketCapacity);

  stk::mesh::NgpMesh ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_1|block_3|        |block_1|block_3|block_3|
//   -------------------------        -------------------------
//
TEST_F(BucketLayoutModification, ChangeBucketContents_WithCopy)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 2;
  setup_mesh_3hex_2block(bucketCapacity);

  stk::mesh::NgpMesh ngpMesh = get_bulk().get_updated_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_3", {2,3}}});
}

}
