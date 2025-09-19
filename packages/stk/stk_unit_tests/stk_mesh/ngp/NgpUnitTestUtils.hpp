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
#include "Kokkos_Macros.hpp"

namespace ngp_unit_test_utils {

constexpr unsigned MaxNumParts = 4;

struct BucketContents
{
  std::vector<std::string> partNames;
  std::vector<stk::mesh::EntityId> entities;
};

template<typename DualViewType>
DualViewType create_dualview(const std::string& name, unsigned size)
{
  DualViewType result(name, size);

  Kokkos::deep_copy(result.view_host(), 0);
  result.template modify<typename DualViewType::host_mirror_space>();
  result.template sync<typename DualViewType::execution_space>();

  return result;
}

inline void setup_mesh_4hex_4block(stk::mesh::BulkData& bulk, unsigned /*bucketCapacity*/)
{
  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(4);
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_3hex_3block(stk::mesh::BulkData& bulk, unsigned /*bucketCapacity*/)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
                         "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_3hex_2block(stk::mesh::BulkData& bulk, unsigned /*bucketCapacity*/)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                         "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

inline void setup_mesh_2hex_2block(stk::mesh::BulkData& bulk, unsigned /*bucketCapacity*/)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
}

struct CheckBucketParts {
  using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;
  CheckBucketParts(
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
    for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
      for (unsigned partIdx = 0; partIdx < MaxNumParts; ++partIdx) {
        unsigned partOffset = bucketIdx*MaxNumParts + partIdx;
        if (bucketPartOrdinals[partOffset] != stk::mesh::InvalidOrdinal) {
          NGP_EXPECT_TRUE(ngpMesh.get_bucket(bucketRank, bucketIdx).member(bucketPartOrdinals[partOffset]));
        }
      }
    }
  }

private:
  stk::mesh::NgpMesh ngpMesh;
  BucketPartOrdinalType bucketPartOrdinals;
  size_t numBuckets;
  stk::topology::rank_t bucketRank;
};

inline void check_bucket_layout(const stk::mesh::BulkData& bulk,
                                const std::vector<BucketContents> & expectedBucketLayout,
                                const stk::topology::rank_t bucketRank = stk::topology::ELEM_RANK)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector & buckets = bulk.buckets(bucketRank);
  size_t numBuckets = buckets.size();
  ASSERT_EQ(numBuckets, expectedBucketLayout.size()) << "Found " << numBuckets << " Host Buckets when expecting "
                                                     << expectedBucketLayout.size();

  size_t numEntitiesAcrossBuckets = 0;
  for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    const BucketContents & bucketContents = expectedBucketLayout[bucketIdx];

    for (const std::string& partName : bucketContents.partNames) {
      const stk::mesh::Part& expectedPart = *meta.get_part(partName);
      EXPECT_TRUE(bucket.member(expectedPart)) << "Host Bucket " << bucket.bucket_id() << " not a member of Part "
                                               << expectedPart.name();
    }

    numEntitiesAcrossBuckets += bucket.size();
    ASSERT_EQ(bucket.size(), bucketContents.entities.size()) << "Found " << bucket.size()
                                                             << " Entities in Host Bucket when expecting "
                                                             << bucketContents.entities.size();
    for (unsigned i = 0; i < bucket.size(); ++i) {
      EXPECT_EQ(bulk.identifier(bucket[i]), bucketContents.entities[i]) << "Found " << bucket[i]
                                                                        << " in Host Bucket when expecting "
                                                                        << bucketContents.entities[i];
    }
  }

  using BucketPartOrdinalType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;
  BucketPartOrdinalType bucketPartOrdinals("bucketPartOrdinals", numBuckets*MaxNumParts);
  BucketPartOrdinalType::host_mirror_type hostBucketPartOrdinals = Kokkos::create_mirror_view(bucketPartOrdinals);
  Kokkos::deep_copy(hostBucketPartOrdinals, stk::mesh::InvalidOrdinal);
  for (size_t bucketIdx = 0; bucketIdx < buckets.size(); ++bucketIdx) {
    const unsigned numExpectedParts = expectedBucketLayout[bucketIdx].partNames.size();
    STK_ThrowRequireMsg(numExpectedParts <= MaxNumParts, "Checking more Parts than test fixture supports");
    for (size_t partIdx = 0; partIdx < numExpectedParts; ++partIdx) {
      const std::string& partName = expectedBucketLayout[bucketIdx].partNames[partIdx];
      unsigned partOffset = bucketIdx*MaxNumParts + partIdx;
      hostBucketPartOrdinals[partOffset] = meta.get_part(partName)->mesh_meta_data_ordinal();
    }
  }
  Kokkos::deep_copy(bucketPartOrdinals, hostBucketPartOrdinals);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  CheckBucketParts checkBucketParts(ngpMesh, bucketPartOrdinals, numBuckets, bucketRank);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), checkBucketParts);

  ASSERT_EQ(ngpMesh.num_buckets(bucketRank), numBuckets) << "Found " << ngpMesh.num_buckets(bucketRank)
                                                         << " Device Buckets when expecting " << numBuckets;

  for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::NgpMesh::BucketType & ngpBucket = ngpMesh.get_bucket(bucketRank, bucketIdx);
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    ASSERT_EQ(bucket.size(), ngpBucket.size());
  }

  using BucketEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
  BucketEntitiesType bucketEntities("bucketEntities", numEntitiesAcrossBuckets);
  BucketEntitiesType::host_mirror_type hostBucketEntities = Kokkos::create_mirror_view(bucketEntities);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(size_t /*index*/) {
      size_t idx = 0;
      for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
        const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(bucketRank, bucketIdx);
        for (size_t i = 0; i < bucket.size(); ++i) {
          bucketEntities[idx++] = bucket[i];
        }
      }
    });

  Kokkos::deep_copy(hostBucketEntities, bucketEntities);

  size_t index = 0;
  for (size_t bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    const stk::mesh::Bucket & bucket = *buckets[bucketIdx];
    for (unsigned i = 0; i < bucket.size(); ++i) {
      const stk::mesh::Entity deviceEntity = hostBucketEntities[index++];
      EXPECT_EQ(bucket[i], deviceEntity) << "Found " << deviceEntity << " in Device Bucket when expecting "
                                         << bucket[i];
    }
  }
}

template <typename DeviceMesh, typename DeviceEntityViewType, typename DevicePartOrdinalsViewType>
inline void check_entity_parts_on_device(DeviceMesh const& deviceMesh,
                                         DeviceEntityViewType const& entities,
                                         DevicePartOrdinalsViewType const& addPartOrdinals,
                                         DevicePartOrdinalsViewType const& removePartOrdinals,
                                         const stk::topology::rank_t bucketRank = stk::topology::ELEM_RANK)
{
  using TeamType = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<> teamPolicy(entities.extent(0), Kokkos::AUTO);

  Kokkos::View<int*, stk::ngp::UVMMemSpace> result = Kokkos::View<int*, stk::ngp::UVMMemSpace>("", 1);

  Kokkos::parallel_for(teamPolicy,
    KOKKOS_LAMBDA(TeamType const& team) {
      auto idx = team.league_rank();
      auto fastMeshIndex = deviceMesh.fast_mesh_index(entities(idx));
      auto bucketId = fastMeshIndex.bucket_id;
      auto& bucket = deviceMesh.get_bucket(bucketRank, bucketId);
      auto& bucketPartOrdinals = bucket.get_part_ordinals();
      int failCount = 0;

      for (unsigned i = 0; i < addPartOrdinals.extent(0); ++i) {
        auto ordinal = addPartOrdinals(i);
        auto it = Kokkos::Experimental::find(team, bucketPartOrdinals, ordinal);
        if (it == Kokkos::Experimental::end(bucketPartOrdinals))
          failCount++;
      }

      for (unsigned i = 0; i < removePartOrdinals.extent(0); ++i) {
        auto ordinal = removePartOrdinals(i);
        auto it = Kokkos::Experimental::find(team, bucketPartOrdinals, ordinal);
        if (it != Kokkos::Experimental::end(bucketPartOrdinals))
          failCount++;
      }

      Kokkos::single(Kokkos::PerTeam(team), [=]() { result(0) = failCount; });
    }
  );
  Kokkos::fence();

  EXPECT_EQ(result(0), 0);
}

} // ngp_unit_test_utils

#endif

