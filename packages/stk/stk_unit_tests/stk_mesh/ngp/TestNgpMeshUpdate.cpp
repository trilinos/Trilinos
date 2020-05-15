#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpAtomics.hpp>
#include <stk_mesh/base/NgpMultistateField.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

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

struct BucketContents
{
  std::string partName;
  std::vector<stk::mesh::EntityId> elements;
};

class BucketLayoutModification : public stk::unit_test_util::MeshFixture
{
public:
  void setup_mesh_3hex_3block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
                           "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_3hex_2block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                           "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void check_bucket_layout(const std::vector<BucketContents> & expectedBucketLayout)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
    size_t numBuckets = buckets.size();
    ASSERT_EQ(numBuckets, expectedBucketLayout.size());

    size_t numElemsAcrossBuckets = 0;
    for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
      const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
      const BucketContents & bucketContents = expectedBucketLayout[bucketIdx];

      const stk::mesh::Part & expectedPart = *get_meta().get_part(bucketContents.partName);
      EXPECT_TRUE(bucket.member(expectedPart));

      numElemsAcrossBuckets += bucket.size();
      ASSERT_EQ(bucket.size(), bucketContents.elements.size());
      for (unsigned i = 0; i < bucket.size(); ++i) {
        EXPECT_EQ(get_bulk().identifier(bucket[i]), bucketContents.elements[i]);
      }
    }

    using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::mesh::MemSpace>;
    BucketPartOrdinalType bucketPartOrdinals("bucketPartOrdinals", numBuckets);
    BucketPartOrdinalType::HostMirror hostBucketPartOrdinals = Kokkos::create_mirror_view(bucketPartOrdinals);
    for (size_t i = 0; i < buckets.size(); ++i) {
      hostBucketPartOrdinals[i] = get_meta().get_part(expectedBucketLayout[i].partName)->mesh_meta_data_ordinal();
    }
    Kokkos::deep_copy(bucketPartOrdinals, hostBucketPartOrdinals);

    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(size_t /*index*/) {
                           NGP_ASSERT_EQ(ngpMesh.num_buckets(stk::topology::ELEM_RANK), numBuckets);
                           for (unsigned i = 0; i < numBuckets; ++i) {
                             NGP_EXPECT_TRUE(ngpMesh.get_bucket(stk::topology::ELEM_RANK, i).member(bucketPartOrdinals[i]));
                           }
                         });

    using BucketEntitiesType = Kokkos::View<stk::mesh::EntityId*, stk::mesh::MemSpace>;
    BucketEntitiesType bucketEntities("bucketEntities", numElemsAcrossBuckets+numBuckets);
    BucketEntitiesType::HostMirror hostBucketEntities = Kokkos::create_mirror_view(bucketEntities);
    size_t index = 0;
    for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
      const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
      hostBucketEntities[index++] = bucket.size();
      for (unsigned i = 0; i < bucket.size(); ++i) {
        hostBucketEntities[index++] = get_bulk().identifier(bucket[i]);
      }
    }
    Kokkos::deep_copy(bucketEntities, hostBucketEntities);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(size_t /*index*/) {
                           NGP_ASSERT_EQ(ngpMesh.num_buckets(stk::topology::ELEM_RANK), numBuckets);
                           size_t idx = 0;
                           for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                                                                                                bucketIdx);
                             NGP_ASSERT_EQ(bucket.size(), bucketEntities[idx++]);
                             for (size_t i = 0; i < bucket.size(); ++i) {
                               NGP_EXPECT_EQ(ngpMesh.identifier(bucket[i]), bucketEntities[idx++]);
                             }
                           }
                         });
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

  check_bucket_layout({{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  check_bucket_layout({{"block_1", {1,2}}, {"block_3", {3}}});
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

  check_bucket_layout({{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  check_bucket_layout({{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
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

  check_bucket_layout({{"block_1", {1,2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  ngpMesh.update_mesh();

  check_bucket_layout({{"block_1", {1}}, {"block_3", {2,3}}});

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

  check_bucket_layout({{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  check_bucket_layout({{"block_1", {1,2}}, {"block_3", {3}}});
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

  check_bucket_layout({{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  check_bucket_layout({{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
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

  check_bucket_layout({{"block_1", {1,2}}, {"block_3", {3}}});

  get_bulk().modification_begin();
  stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
  stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
  get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
  get_bulk().modification_end();

  get_bulk().update_ngp_mesh();

  check_bucket_layout({{"block_1", {1}}, {"block_3", {2,3}}});
}

}
