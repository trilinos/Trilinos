#ifndef NgpUnitTestUtils_hpp
#define NgpUnitTestUtils_hpp

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_ngp_test/ngp_test.hpp>

namespace ngp_unit_test_utils {

struct BucketContents
{
  std::string partName;
  std::vector<stk::mesh::EntityId> entities;
};

template<typename DualViewType>
DualViewType create_dualview(const std::string& name, unsigned size)
{
  DualViewType result(name, size);

  Kokkos::deep_copy(result.h_view, 0);
  result.template modify<typename DualViewType::host_mirror_space>();
  result.template sync<typename DualViewType::execution_space>();

  return result;
}

inline void setup_mesh_4hex_4block(stk::mesh::BulkData& bulk, unsigned bucketCapacity)
{
  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(4);
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_3hex_3block(stk::mesh::BulkData& bulk, unsigned bucketCapacity)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
                         "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_3hex_2block(stk::mesh::BulkData& bulk, unsigned bucketCapacity)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                         "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_2hex_2block(stk::mesh::BulkData& bulk, unsigned bucketCapacity)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

struct CheckPartMembership {
  using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;
  CheckPartMembership(
      const stk::mesh::NgpMesh& _ngpMesh, BucketPartOrdinalType _bucketPartOrdinals, size_t _numBuckets,
      const stk::topology::rank_t _bucketRank)
    : ngpMesh(_ngpMesh),
      bucketPartOrdinals(_bucketPartOrdinals),
      numBuckets(_numBuckets),
      bucketRank(_bucketRank)
  {
  }

  KOKKOS_FUNCTION
  void operator()(size_t) const
  {
    for (unsigned i = 0; i < numBuckets; ++i) {
      NGP_EXPECT_TRUE(ngpMesh.get_bucket(bucketRank, i).member(bucketPartOrdinals[i]));
    }
  }

private:
  stk::mesh::NgpMesh ngpMesh;
  BucketPartOrdinalType bucketPartOrdinals;
  size_t numBuckets;
  stk::topology::rank_t bucketRank;
};

inline void check_bucket_layout(const stk::mesh::BulkData& bulk,  const std::vector<BucketContents> & expectedBucketLayout,
                                const stk::topology::rank_t bucketRank = stk::topology::ELEM_RANK)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector & buckets = bulk.buckets(bucketRank);
  size_t numBuckets = buckets.size();
  ASSERT_EQ(numBuckets, expectedBucketLayout.size());

  size_t numEntitiesAcrossBuckets = 0;
  for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    const BucketContents & bucketContents = expectedBucketLayout[bucketIdx];

    const stk::mesh::Part & expectedPart = *meta.get_part(bucketContents.partName);
    EXPECT_TRUE(bucket.member(expectedPart));

    numEntitiesAcrossBuckets += bucket.size();
    ASSERT_EQ(bucket.size(), bucketContents.entities.size());
    for (unsigned i = 0; i < bucket.size(); ++i) {
      EXPECT_EQ(bulk.identifier(bucket[i]), bucketContents.entities[i]);
    }
  }

  using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;
  BucketPartOrdinalType bucketPartOrdinals("bucketPartOrdinals", numBuckets);
  BucketPartOrdinalType::HostMirror hostBucketPartOrdinals = Kokkos::create_mirror_view(bucketPartOrdinals);
  for (size_t i = 0; i < buckets.size(); ++i) {
    hostBucketPartOrdinals[i] = meta.get_part(expectedBucketLayout[i].partName)->mesh_meta_data_ordinal();
  }
  Kokkos::deep_copy(bucketPartOrdinals, hostBucketPartOrdinals);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  CheckPartMembership checkElementMembership(ngpMesh, bucketPartOrdinals, numBuckets, bucketRank);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), checkElementMembership);

  ASSERT_EQ(ngpMesh.num_buckets(bucketRank), numBuckets);

  for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::NgpMesh::BucketType & ngpBucket = ngpMesh.get_bucket(bucketRank, bucketIdx);
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    ASSERT_EQ(bucket.size(), ngpBucket.size());
  }

  using BucketEntitiesType = Kokkos::View<stk::mesh::EntityId*, stk::ngp::MemSpace>;
  BucketEntitiesType bucketEntities("bucketEntities", numEntitiesAcrossBuckets);
  BucketEntitiesType::HostMirror hostBucketEntities = Kokkos::create_mirror_view(bucketEntities);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(size_t /*index*/) {
      size_t idx = 0;
      for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
        const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(bucketRank, bucketIdx);
        for (size_t i = 0; i < bucket.size(); ++i) {
          bucketEntities[idx++] = ngpMesh.identifier(bucket[i]);
        }
      }
    });

  Kokkos::deep_copy(hostBucketEntities, bucketEntities);

  size_t index = 0;
  for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    for (unsigned i = 0; i < bucket.size(); ++i) {
      EXPECT_EQ(bulk.identifier(bucket[i]), hostBucketEntities[index++]);
    }
  }
}

} // ngp_unit_test_utils

#endif

