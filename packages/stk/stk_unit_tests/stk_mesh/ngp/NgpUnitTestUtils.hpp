#ifndef NgpUnitTestUtils_hpp
#define NgpUnitTestUtils_hpp

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

namespace ngp_unit_test_utils {

struct BucketContents
{
  std::string partName;
  std::vector<stk::mesh::EntityId> elements;
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

inline void setup_mesh_2hex_3block(stk::mesh::BulkData& bulk, unsigned bucketCapacity)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void check_bucket_layout(const stk::mesh::BulkData& bulk, const std::vector<BucketContents> & expectedBucketLayout)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::ELEM_RANK);
  size_t numBuckets = buckets.size();
  ASSERT_EQ(numBuckets, expectedBucketLayout.size());

  size_t numElemsAcrossBuckets = 0;
  for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    const BucketContents & bucketContents = expectedBucketLayout[bucketIdx];

    const stk::mesh::Part & expectedPart = *meta.get_part(bucketContents.partName);
    EXPECT_TRUE(bucket.member(expectedPart));

    numElemsAcrossBuckets += bucket.size();
    ASSERT_EQ(bucket.size(), bucketContents.elements.size());
    for (unsigned i = 0; i < bucket.size(); ++i) {
      EXPECT_EQ(bulk.identifier(bucket[i]), bucketContents.elements[i]);
    }
  }

  using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::mesh::MemSpace>;
  BucketPartOrdinalType bucketPartOrdinals("bucketPartOrdinals", numBuckets);
  BucketPartOrdinalType::HostMirror hostBucketPartOrdinals = Kokkos::create_mirror_view(bucketPartOrdinals);
  for (size_t i = 0; i < buckets.size(); ++i) {
    hostBucketPartOrdinals[i] = meta.get_part(expectedBucketLayout[i].partName)->mesh_meta_data_ordinal();
  }
  Kokkos::deep_copy(bucketPartOrdinals, hostBucketPartOrdinals);

  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
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
      hostBucketEntities[index++] = bulk.identifier(bucket[i]);
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


} // ngp_unit_test_utils

#endif

