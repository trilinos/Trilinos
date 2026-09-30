// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/FieldBLAS.hpp"
#include "stk_mesh/base/DumpMeshInfo.hpp"
#include "stk_ngp_test/ngp_test.hpp"
#include "stk_topology/topology.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_unit_test_utils/DeviceBucketTestUtils.hpp"
#include "stk_unit_test_utils/GetMeshSpec.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

#ifdef STK_USE_DEVICE_MESH

enum class Status 
{
  isOwned,
  notOwned,
  isShared,
  notShared,
  isInAura,
  notInAura
};

struct OwnershipStatus
{
  Status owned;
  Status shared;
  Status inAura;
};

struct CheckIsOwned {
  template <typename Bucket>
  KOKKOS_INLINE_FUNCTION
  bool operator()(const Bucket& bucket) const { return bucket->is_owned(); }
};

struct CheckIsShared {
  template <typename Bucket>
  KOKKOS_INLINE_FUNCTION
  bool operator()(const Bucket& bucket) const { return bucket->is_shared(); }
};

struct CheckIsInAura {
  template <typename Bucket>
  KOKKOS_INLINE_FUNCTION
  bool operator()(const Bucket& bucket) const { return bucket->is_in_aura(); }
};

class DeviceMeshParallelTester : public ::ngp_testing::Test
{
public:
  using PartOrdinalType = unsigned;

  DeviceMeshParallelTester() = default;

  DeviceMeshParallelTester(unsigned spatialDim, stk::mesh::BulkData::AutomaticAuraOption auraOp)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(spatialDim);
    builder.set_initial_bucket_capacity(10u);
    builder.set_maximum_bucket_capacity(10u);
    builder.set_aura_option(auraOp);
    bulk = builder.create();
    meta = &bulk->mesh_meta_data();
  }

  void setup_mesh(std::string& meshDesc)
  {
    stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);
  }

  void check_all_bucket_ownerships(OwnershipStatus const& status = {Status::notOwned, Status::notShared, Status::notInAura})
  {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto& deviceBucketRepo = deviceMesh.get_device_bucket_repository();
    auto rankCount = meta->entity_rank_count();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < rankCount; ++rank) {
        auto numBuckets = deviceBucketRepo.num_buckets(rank);
        
        for (unsigned i = 0; i < numBuckets; ++i) {
          auto bucket = deviceBucketRepo.get_bucket(rank, i);

          NGP_EXPECT_EQ(bucket->is_owned(), (status.owned == Status::isOwned) ? true : false);
          NGP_EXPECT_EQ(bucket->is_shared(), (status.shared == Status::isShared) ? true : false);
          NGP_EXPECT_EQ(bucket->is_in_aura(), (status.inAura == Status::isInAura) ? true : false);
        }
      }
    });
  }

  void declare_new_parts(std::vector<std::string>&& partsToCreate, stk::topology::topology_t topology = stk::topology::topology_t::INVALID_TOPOLOGY)
  {
    for (auto partName : partsToCreate) {
      if (topology == stk::topology::topology_t::INVALID_TOPOLOGY) {
        meta->declare_part(partName);
      } else {
        meta->declare_part_with_topology(partName, topology);
      }
    }
  }

  void check_bucket_ownerships(stk::mesh::EntityRank rank, std::vector<PartOrdinalType>&& bucketIdxToCheck, OwnershipStatus const& status = {})
  {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

    Kokkos::View<unsigned*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> bucketIdx(bucketIdxToCheck.data(), bucketIdxToCheck.size());
    auto bucketIdxOnDevice = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace{}, bucketIdx);
  
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      for (unsigned i = 0; i < bucketIdxOnDevice.extent(0); ++i) {
        auto idx = bucketIdxOnDevice(i);
        auto bucket = deviceBucketRepo.get_bucket(rank, idx);

        NGP_EXPECT_EQ(bucket->is_owned(), (status.owned == Status::isOwned) ? true : false);
        NGP_EXPECT_EQ(bucket->is_shared(), (status.shared == Status::isShared) ? true : false);
        NGP_EXPECT_EQ(bucket->is_in_aura(), (status.inAura == Status::isInAura) ? true : false);
      }
    });
  }

  template <typename EntityView, typename Op>
  bool check_bucket_ownership(EntityView entities, Op const& testOp)
  {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

    bool allOwned = false;
    Kokkos::parallel_reduce(entities.extent(0),
      KOKKOS_LAMBDA(const int i, bool& update) {
        auto entity = entities(i);
        auto fmi = deviceMesh.fast_mesh_index(entity);
        auto rank = deviceMesh.entity_rank(entity);
        auto bucket = deviceBucketRepo.get_bucket(rank, fmi.bucket_id);
        update &= testOp(bucket);
      },
    Kokkos::LAnd<bool>(allOwned));

    return allOwned;
  }

  bool is_owned(stk::mesh::EntityVector&& entitiesToCheck)
  {
    using EntityViewOnHost = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    EntityViewOnHost entitiesOnHost(entitiesToCheck.data(), entitiesToCheck.size());
    auto entities = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, entitiesOnHost);
    return check_bucket_ownership(entities, CheckIsOwned{});
  }

  bool is_shared(stk::mesh::EntityVector&& entitiesToCheck)
  {
    using EntityViewOnHost = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    EntityViewOnHost entitiesOnHost(entitiesToCheck.data(), entitiesToCheck.size());
    auto entities = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, entitiesOnHost);
    return check_bucket_ownership(entities, CheckIsShared{});
  }

  bool is_in_aura(stk::mesh::EntityVector&& entitiesToCheck)
  {
    using EntityViewOnHost = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    EntityViewOnHost entitiesOnHost(entitiesToCheck.data(), entitiesToCheck.size());
    auto entities = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, entitiesOnHost);
    return check_bucket_ownership(entities, CheckIsInAura{});
  }

  void check_entity_ownership(stk::mesh::EntityVector&& entitiesToCheck, OwnershipStatus&& expectedStatus)
  {
    if (expectedStatus.owned == Status::isOwned)
      EXPECT_TRUE(is_owned(std::move(entitiesToCheck)));
    else
      EXPECT_FALSE(is_owned(std::move(entitiesToCheck)));
    if (expectedStatus.shared == Status::isShared)
      EXPECT_TRUE(is_shared(std::move(entitiesToCheck)));
    else
      EXPECT_FALSE(is_shared(std::move(entitiesToCheck)));   
    if (expectedStatus.inAura == Status::isInAura)
      EXPECT_TRUE(is_in_aura(std::move(entitiesToCheck)));
    else
      EXPECT_FALSE(is_in_aura(std::move(entitiesToCheck)));
  }

  void check_device_bucket_part_membership(stk::mesh::EntityVector&& entitiesToCheck, std::vector<std::string>&& includedPartNames, std::vector<std::string>&& excludedPartNames)
  {
    using EntityViewOnHost = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

    Kokkos::View<unsigned*> expectedPartOrdinals("", includedPartNames.size());
    Kokkos::View<unsigned*> excludedPartOrdinals("", excludedPartNames.size());

    fill_view(includedPartNames, expectedPartOrdinals);
    fill_view(excludedPartNames, excludedPartOrdinals);

    EntityViewOnHost entitiesOnHost(entitiesToCheck.data(), entitiesToCheck.size());
    auto entities = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, entitiesOnHost);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      for (unsigned i = 0; i < entities.extent(0); ++i) {
        auto entity = entities(i);
        auto rank = deviceMesh.entity_rank(entity);
        auto fmi = deviceMesh.fast_mesh_index(entity);
        auto bucket = deviceBucketRepo.get_bucket(rank, fmi.bucket_id);

        for (unsigned j = 0; j < expectedPartOrdinals.extent(0); ++j) {
          NGP_EXPECT_TRUE(bucket->member(expectedPartOrdinals(j)));
        }
        for (unsigned j = 0; j < excludedPartOrdinals.extent(0); ++j) {
          NGP_EXPECT_FALSE(bucket->member(excludedPartOrdinals(j)));
        }
      }
    });
  }

  template <typename VectorType, typename ViewType>
  void fill_view(VectorType const& vector, ViewType view)
  {
    auto hostView = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, view);

    for (unsigned i = 0; i < hostView.extent(0); ++i) {
      if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vector)>::value_type, stk::mesh::Part*>) {
        hostView(i) = vector[i]->mesh_meta_data_ordinal();
      } else if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vector)>::value_type, std::string>) {
        hostView(i) = meta->get_part(vector[i])->mesh_meta_data_ordinal();
      } else {
        hostView(i) = vector[i];
      }
    }
    Kokkos::deep_copy(view, hostView);
  }

  template <typename EntityFmiVector, typename ViewType>
  void get_entities(stk::mesh::EntityRank rank, EntityFmiVector&& entityFmi, ViewType view)
  {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);

    Kokkos::View<stk::mesh::FastMeshIndex*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> fmi(entityFmi.data(), entityFmi.size());
    auto deviceFmi = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, fmi);

    Kokkos::parallel_for(entityFmi.size(),
      KOKKOS_LAMBDA(const int i) {
        view(i) = deviceMesh.get_entity(rank, deviceFmi(i));
      }
    );
  }

  void check_bucket_count(stk::mesh::EntityRank rank, unsigned count) {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

    auto numBuckets = deviceBucketRepo.num_buckets(rank);
    EXPECT_EQ(numBuckets, count);
  }

  template <typename DeviceMeshType, typename FmiVectorType, typename ViewType>
  void fill_entities(stk::mesh::EntityRank rank, DeviceMeshType const& deviceMesh, FmiVectorType&& entityFmiVec, ViewType entityView)
  {
    Kokkos::View<stk::mesh::FastMeshIndex*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> hostFmi(entityFmiVec.data(), entityFmiVec.size());
    auto deviceFmi = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace{}, hostFmi);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int) {
      for (unsigned i = 0; i < deviceFmi.extent(0); ++i) {
        auto fmi = deviceFmi(i);
        auto entity = deviceMesh.get_entity(rank, fmi);
        entityView(i) = entity;
      }
    });
  }

  void change_entity_parts_in_parallel(stk::mesh::EntityVector const& entities, std::vector<std::string>&& addPartStrs, std::vector<std::string>&& removePartStrs)
  {
    Kokkos::View<stk::mesh::Entity*> entityView("", entities.size());
    Kokkos::View<unsigned*> addPartOrdinals("", addPartStrs.size());
    Kokkos::View<unsigned*> removePartOrdinals("", removePartStrs.size());

    fill_view(entities, entityView);
    fill_view(addPartStrs, addPartOrdinals);
    fill_view(removePartStrs, removePartOrdinals);

    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    deviceMesh.batch_change_entity_parts(entityView, addPartOrdinals, removePartOrdinals);
  }

  stk::mesh::MetaData* meta;
  std::shared_ptr<stk::mesh::BulkData> bulk;
};

class DeviceMeshParallelTester_NoAura : public DeviceMeshParallelTester
{
public:
  DeviceMeshParallelTester_NoAura()
    : DeviceMeshParallelTester(3, stk::mesh::BulkData::NO_AUTO_AURA)
  {}
};

NGP_TEST_F(DeviceMeshParallelTester_NoAura, check_one_proc_bucket_ownership)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  unsigned numBlocks = 2;
  unsigned numProcs = 1;
  auto desc = stk::unit_test_util::get_many_block_mesh_desc(numBlocks, numProcs);
  setup_mesh(desc);

  check_all_bucket_ownerships({Status::isOwned, Status::notShared, Status::notInAura});
}

NGP_TEST_F(DeviceMeshParallelTester_NoAura, check_two_procs_bucket_ownership_one_block)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,1,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_FALSE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
  } else {
    auto node15 = bulk->get_entity(stk::topology::NODE_RANK, 15);
    auto node9 = bulk->get_entity(stk::topology::NODE_RANK, 9);
    EXPECT_TRUE(is_owned({node15}));
    EXPECT_FALSE(is_shared({node15}));
    EXPECT_FALSE(is_in_aura({node15}));
    EXPECT_TRUE(is_owned({node9}));
    EXPECT_FALSE(is_shared({node9}));
    EXPECT_FALSE(is_in_aura({node9}));
  }
}

NGP_TEST_F(DeviceMeshParallelTester_NoAura, check_two_procs_bucket_ownership_two_blocks)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  setup_mesh(desc);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_FALSE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
  } else {
    auto node15 = bulk->get_entity(stk::topology::NODE_RANK, 15);
    auto node9 = bulk->get_entity(stk::topology::NODE_RANK, 9);
    EXPECT_TRUE(is_owned({node15}));
    EXPECT_FALSE(is_shared({node15}));
    EXPECT_FALSE(is_in_aura({node15}));
    EXPECT_TRUE(is_owned({node9}));
    EXPECT_FALSE(is_shared({node9}));
    EXPECT_FALSE(is_in_aura({node9}));
  }
}

NGP_TEST_F(DeviceMeshParallelTester_NoAura, check_three_procs_bucket_ownership)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,9,10,11,12,13,14,15,16,block_2\n" \
                     "2,3,HEX_8,13,14,15,16,17,18,19,20,block_3";
  setup_mesh(desc);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_FALSE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
  } else if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 1) {
    auto node9 = bulk->get_entity(stk::topology::NODE_RANK, 9);
    auto node15 = bulk->get_entity(stk::topology::NODE_RANK, 15);
    EXPECT_TRUE(is_owned({node9}));
    EXPECT_FALSE(is_shared({node9}));
    EXPECT_FALSE(is_in_aura({node9}));
    EXPECT_TRUE(is_owned({node15}));
    EXPECT_TRUE(is_shared({node15}));
    EXPECT_FALSE(is_in_aura({node15}));
  } else {
    auto node15 = bulk->get_entity(stk::topology::NODE_RANK, 15);
    auto node20 = bulk->get_entity(stk::topology::NODE_RANK, 20);
    EXPECT_FALSE(is_owned({node15}));
    EXPECT_TRUE(is_shared({node15}));
    EXPECT_FALSE(is_in_aura({node15}));
    EXPECT_TRUE(is_owned({node20}));
    EXPECT_FALSE(is_shared({node20}));
    EXPECT_FALSE(is_in_aura({node20}));
  }
}

class DeviceMeshParallelMeshModTester_NoAura : public DeviceMeshParallelTester_NoAura
{
};

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_unshared_elems_diff_parts)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1", "new_elem_part2"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 1);
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 9);
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  }

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = (rank == 0) ? "new_elem_part1" : "new_elem_part2";

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    change_entity_parts_in_parallel({elem}, {newPart}, {});
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 9);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    change_entity_parts_in_parallel({elem}, {newPart}, {});
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_unshared_elems_same_parts)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  setup_mesh(desc);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 1);
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 9);
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  }

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = (rank == 1) ? "block_1" : "block_2";

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    change_entity_parts_in_parallel({elem}, {newPart}, {});
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 9);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    change_entity_parts_in_parallel({elem}, {newPart}, {});
    EXPECT_TRUE(is_owned({node}));
    EXPECT_FALSE(is_shared({node}));
    EXPECT_FALSE(is_in_aura({node}));
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, require_entity_owner)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1", "new_elem_part2"}, stk::topology::HEX_8);
  declare_new_parts({"new_node_part1", "new_node_part2"}, stk::topology::NODE);

  if (myRank == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem}, {"new_elem_part1"}, {}));
    Kokkos::fence();

    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_NO_THROW(change_entity_parts_in_parallel({node5}, {"new_node_part1"}, {}));
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node9 = bulk->get_entity(stk::topology::NODE_RANK, 9);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem}, {"new_elem_part2"}, {}));
    Kokkos::fence();

    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    EXPECT_NO_THROW(change_entity_parts_in_parallel({node9}, {"new_node_part2"}, {}));
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_only_shared_node_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    check_bucket_count(stk::topology::NODE_RANK, 2);

    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {"new_part1"});
  } else {
    check_bucket_count(stk::topology::NODE_RANK, 2);

    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    EXPECT_TRUE(is_owned({node12}));
    EXPECT_FALSE(is_shared({node12}));
    EXPECT_FALSE(is_in_aura({node12}));
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {"new_part1"});
  }

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node7 = bulk->get_entity(stk::topology::NODE_RANK, 7);
    change_entity_parts_in_parallel({node5}, {"new_part1"}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 3);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    EXPECT_TRUE(is_owned({node7}));
    EXPECT_TRUE(is_shared({node7}));
    EXPECT_FALSE(is_in_aura({node7}));

    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_part1"});
    check_device_bucket_part_membership({node7}, {"block_1", "block_2"}, {"new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_part1"}, {});
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({}, {}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 3);
    EXPECT_TRUE(is_owned({node12}));
    EXPECT_FALSE(is_shared({node12}));
    EXPECT_FALSE(is_in_aura({node12}));
    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    EXPECT_FALSE(is_owned({node8}));
    EXPECT_TRUE(is_shared({node8}));
    EXPECT_FALSE(is_in_aura({node8}));

    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_part1"});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2"}, {"new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_part1"}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_all_nodes_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);

    stk::mesh::EntityVector entitiesToApplyParts;
    bulk->get_entities(stk::topology::NODE_RANK, meta->universal_part(), entitiesToApplyParts);
    change_entity_parts_in_parallel(entitiesToApplyParts, {"new_part1"}, {});

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));

    check_device_bucket_part_membership({node1}, {"block_1", "new_part1"}, {"block_2"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_part1"}, {});
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({}, {}, {});

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node12}));
    EXPECT_FALSE(is_shared({node12}));
    EXPECT_FALSE(is_in_aura({node12}));
    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));

    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_part1"}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_non_inducible_parts_one_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    change_entity_parts_in_parallel({elem}, {"new_part1"}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {});

    check_bucket_count(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(is_owned({elem}));
    EXPECT_FALSE(is_shared({elem}));
    EXPECT_FALSE(is_in_aura({elem}));
    check_device_bucket_part_membership({elem}, {"block_1", "new_part1"}, {"block_2"});
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    change_entity_parts_in_parallel({}, {}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node12}));
    EXPECT_FALSE(is_shared({node12}));
    EXPECT_FALSE(is_in_aura({node12}));
    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {});

    check_bucket_count(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(is_owned({elem}));
    EXPECT_FALSE(is_shared({elem}));
    EXPECT_FALSE(is_in_aura({elem}));
    check_device_bucket_part_membership({elem}, {"block_2"}, {"block_1", "new_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_non_inducible_parts_two_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_part1", "new_part2"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    change_entity_parts_in_parallel({elem}, {"new_part1"}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node1}));
    EXPECT_FALSE(is_shared({node1}));
    EXPECT_FALSE(is_in_aura({node1}));
    EXPECT_TRUE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_part1", "new_part2"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {"new_part1", "new_part2"});

    check_bucket_count(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(is_owned({elem}));
    EXPECT_FALSE(is_shared({elem}));
    EXPECT_FALSE(is_in_aura({elem}));
    check_device_bucket_part_membership({elem}, {"block_1", "new_part1"}, {"block_2", "new_part2"});
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    change_entity_parts_in_parallel({elem}, {"new_part2"}, {});
    Kokkos::fence();

    check_bucket_count(stk::topology::NODE_RANK, 2);
    EXPECT_TRUE(is_owned({node12}));
    EXPECT_FALSE(is_shared({node12}));
    EXPECT_FALSE(is_in_aura({node12}));
    EXPECT_FALSE(is_owned({node5}));
    EXPECT_TRUE(is_shared({node5}));
    EXPECT_FALSE(is_in_aura({node5}));
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_part1", "new_part2"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {"new_part1", "new_part2"});

    check_bucket_count(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(is_owned({elem}));
    EXPECT_FALSE(is_shared({elem}));
    EXPECT_FALSE(is_in_aura({elem}));
    check_device_bucket_part_membership({elem}, {"block_2", "new_part2"}, {"block_1", "new_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_ownership_change_entity_parts_modify_elem_with_inducible_parts_by_proc_owning_shared_nodes)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node}, {Status::isOwned, Status::notShared, Status::notInAura});
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 5);
    check_entity_ownership({node}, {Status::notOwned, Status::isShared, Status::notInAura});
  }

  std::string newPart = "new_elem_part1";
  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem}, {newPart}, {});
    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem}, {}, {});
    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_entity_parts_modify_nodes_with_only_inducible_parts_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_node_part1"}, stk::topology::topology_t::NODE);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 7);
    check_entity_ownership({node}, {Status::isOwned, Status::isShared, Status::notInAura});
  } else {
    auto node = bulk->get_entity(stk::topology::NODE_RANK, 5);
    check_entity_ownership({node}, {Status::notOwned, Status::isShared, Status::notInAura});
  }

  std::string newPart = "new_node_part1";
  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({node5, node8}, {newPart}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_node_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_node_part1"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_node_part1"}, {});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem}, {}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2"}, {"new_node_part1"});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_node_part1"}, {});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_node_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_entity_parts_modify_elem_with_only_inducible_parts_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  std::string newPart = "new_elem_part1";
  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem}, {newPart}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({node1}, {"block_1", "new_elem_part1"}, {"block_2"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_elem_part1"}, {});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem}, {}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", "new_elem_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_only_inducible_parts_by_proc_non_owning_shared_nodes)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  std::string newPart = "new_elem_part1";
  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem}, {}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", "new_elem_part1"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_elem_part1"}, {});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem}, {newPart}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", "new_elem_part1"}, {});
    check_device_bucket_part_membership({node12}, {"block_2", "new_elem_part1"}, {"block_1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_only_inducible_parts_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newPart1 = "new_elem_part1";
  std::string newPart2 = "new_elem_part2";
  declare_new_parts({newPart1, newPart2}, stk::topology::topology_t::HEX_8);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem}, {newPart1}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({node1}, {"block_1", newPart1}, {"block_2", newPart2});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", newPart1, newPart2}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newPart1, newPart2}, {});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem}, {newPart2}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", newPart1, newPart2}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newPart1, newPart2}, {});
    check_device_bucket_part_membership({node12}, {"block_2", newPart2}, {"block_1", newPart1});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_nodes_with_only_inducible_elem_parts_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newPart1 = "new_elem_part1";
  std::string newPart2 = "new_elem_part2";
  declare_new_parts({newPart1, newPart2}, stk::topology::topology_t::HEX_8);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({node1, node5}, {newPart1}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", newPart1, newPart2});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {newPart1, newPart2});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2"}, {newPart1, newPart2});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({node12}, {newPart2}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2"}, {newPart1, newPart2});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2"}, {newPart1, newPart2});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newPart1, newPart2});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_mixed_parts_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newRankedPart = "new_elem_part";
  std::string newUnrankedPart = "new_part";
  declare_new_parts({newRankedPart}, stk::topology::topology_t::HEX_8);
  declare_new_parts({newUnrankedPart});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem}, {newRankedPart, newUnrankedPart}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_1", newRankedPart, newUnrankedPart}, {"block_2"});
    check_device_bucket_part_membership({node1}, {"block_1", newRankedPart}, {"block_2", newUnrankedPart});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({}, {}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_2"}, {"block_1", newRankedPart, newUnrankedPart});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newRankedPart, newUnrankedPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_and_node_with_mixed_parts_by_owning_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newRankedPart = "new_elem_part";
  std::string newUnrankedPart = "new_part";
  declare_new_parts({newRankedPart}, stk::topology::topology_t::HEX_8);
  declare_new_parts({newUnrankedPart});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem, node8}, {newRankedPart, newUnrankedPart}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_1", newRankedPart, newUnrankedPart}, {"block_2"});
    check_device_bucket_part_membership({node1}, {"block_1", newRankedPart}, {"block_2", newUnrankedPart});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({}, {}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_2"}, {"block_1", newRankedPart, newUnrankedPart});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newRankedPart, newUnrankedPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_modify_elem_with_mixed_parts_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }
  
  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newRankedPart1 = "new_elem_part1";
  std::string newRankedPart2 = "new_elem_part2";
  std::string newUnrankedPart1 = "new_part1";
  std::string newUnrankedPart2 = "new_part2";
  declare_new_parts({newRankedPart1, newRankedPart2}, stk::topology::topology_t::HEX_8);
  declare_new_parts({newUnrankedPart1, newUnrankedPart2});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    change_entity_parts_in_parallel({elem, node8}, {newRankedPart1, newUnrankedPart1}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_1", newRankedPart1, newUnrankedPart1}, {"block_2", newRankedPart2, newUnrankedPart2});
    check_device_bucket_part_membership({node1}, {"block_1", newRankedPart1}, {"block_2", newUnrankedPart1, newRankedPart2, newUnrankedPart2});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2", newRankedPart1, newRankedPart2}, {newUnrankedPart1, newUnrankedPart2});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart1, newUnrankedPart1, newRankedPart2}, {newUnrankedPart2});
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    change_entity_parts_in_parallel({elem, node12}, {newRankedPart2, newUnrankedPart2}, {});

    check_entity_ownership({elem}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node8}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {"block_2", newRankedPart2, newUnrankedPart2}, {"block_1", newRankedPart1, newUnrankedPart1});
    check_device_bucket_part_membership({node6}, {"block_1", "block_2", newRankedPart1, newRankedPart2}, {newUnrankedPart1, newUnrankedPart2});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart1, newUnrankedPart1, newRankedPart2}, {newUnrankedPart2});
    check_device_bucket_part_membership({node12}, {"block_2", newUnrankedPart2, newRankedPart2}, {"block_1", newRankedPart1, newUnrankedPart1});
  }
}

class DeviceMeshParallelMeshModTester_NoAura_2D : public DeviceMeshParallelMeshModTester_NoAura
{
public:
  DeviceMeshParallelMeshModTester_NoAura_2D()
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(2);
    builder.set_initial_bucket_capacity(10u);
    builder.set_maximum_bucket_capacity(10u);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    bulk = builder.create();
    meta = &bulk->mesh_meta_data();
  }

  void check_field_data(stk::mesh::Entity entity, stk::mesh::Field<double>& field)
  {
    auto& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
    auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        auto fastMeshIndex = deviceMesh.device_mesh_index(entity);
        auto fieldEntityValues = fieldData.entity_values(fastMeshIndex);

        for (stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
          auto fieldValue = static_cast<double>(fieldEntityValues(i));
          NGP_EXPECT_EQ(fieldValue, 1.0);
        }
      }
    );
  };
};

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura_2D, check_change_entity_parts_keyhole_elem_ranked_and_shared_nodes_unranked)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,QUAD_4_2D,1,2,6,5,block_1\n"
                     "0,2,QUAD_4_2D,2,3,7,6,block_1\n"
                     "0,3,QUAD_4_2D,3,4,8,7,block_1\n"
                     "0,4,QUAD_4_2D,5,6,10,9,block_1\n"
                     "1,5,QUAD_4_2D,6,7,11,10,block_2\n"
                     "0,6,QUAD_4_2D,7,8,12,11,block_1\n"
                     "0,7,QUAD_4_2D,9,10,14,13,block_1\n"
                     "0,8,QUAD_4_2D,10,11,15,14,block_1\n"
                     "0,9,QUAD_4_2D,11,12,16,15,block_1";
  setup_mesh(desc);

  std::string newRankedPart = "new_elem_part";
  std::string newUnrankedPart = "new_part";
  declare_new_parts({newRankedPart}, stk::topology::topology_t::QUAD_4_2D);
  declare_new_parts({newUnrankedPart});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 1) {
    auto elem5 = bulk->get_entity(stk::topology::ELEM_RANK, 5);
    auto node6 = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node11 = bulk->get_entity(stk::topology::NODE_RANK, 11);
    change_entity_parts_in_parallel({elem5}, {newRankedPart}, {});

    check_entity_ownership({elem5}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node6}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem5}, {"block_2", newRankedPart}, {"block_1", newUnrankedPart});
    check_device_bucket_part_membership({node6},  {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
    check_device_bucket_part_membership({node11}, {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
  } else {
    auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
    auto node6  = bulk->get_entity(stk::topology::NODE_RANK, 6);
    auto node7  = bulk->get_entity(stk::topology::NODE_RANK, 7);
    auto node10 = bulk->get_entity(stk::topology::NODE_RANK, 10);
    auto node11 = bulk->get_entity(stk::topology::NODE_RANK, 11);
    change_entity_parts_in_parallel({node6, node7, node10, node11}, {newUnrankedPart}, {});

    check_entity_ownership({node6}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({node1},  {"block_1"}, {"block_2", newRankedPart, newUnrankedPart});
    check_device_bucket_part_membership({node6},  {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
    check_device_bucket_part_membership({node7},  {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
    check_device_bucket_part_membership({node10}, {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
    check_device_bucket_part_membership({node11}, {"block_1", "block_2", newRankedPart, newUnrankedPart}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura_2D, check_change_entity_parts_elem_add_and_remove_ranked_parts)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  auto& block1Part = meta->declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  auto& block2Part = meta->declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);

  stk::mesh::Field<double>& oneField = meta->declare_field<double>(stk::topology::NODE_RANK, "field_of_one");
  double initValZero = 0.0;
  stk::mesh::put_field_on_mesh(oneField, block1Part, &initValZero);

  std::string desc = "0,1,QUAD_4_2D,1,2,5,6,block_1\n"
                     "1,2,QUAD_4_2D,2,3,4,5,block_1\n"
                     "|dimension:2";
  setup_mesh(desc);

  stk::mesh::field_fill(1.0, oneField);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node2 = bulk->get_entity(stk::topology::NODE_RANK, 2);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    change_entity_parts_in_parallel({elem}, {block2Part.name()}, {block1Part.name()});

    check_entity_ownership({node2}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {block2Part.name()}, {block1Part.name()});
    check_device_bucket_part_membership({node2}, {block1Part.name(), block2Part.name()}, {});
    check_device_bucket_part_membership({node5}, {block1Part.name(), block2Part.name()}, {});

    check_field_data(elem, oneField);
    check_field_data(node2, oneField);
    check_field_data(node5, oneField);
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node2 = bulk->get_entity(stk::topology::NODE_RANK, 2);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    change_entity_parts_in_parallel({elem}, {}, {});

    check_entity_ownership({node2}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {block1Part.name()}, {block2Part.name()});
    check_device_bucket_part_membership({node2}, {block1Part.name(), block2Part.name()}, {});
    check_device_bucket_part_membership({node5}, {block1Part.name(), block2Part.name()}, {});

    check_field_data(elem, oneField);
    check_field_data(node2, oneField);
    check_field_data(node5, oneField);
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura_2D, check_change_entity_parts_elem_remove_ranked_part_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  auto& block1Part = meta->declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  auto& block2Part = meta->declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);

  stk::mesh::Field<double>& oneField = meta->declare_field<double>(stk::topology::NODE_RANK, "field_of_one");
  double initValZero = 0.0;
  stk::mesh::put_field_on_mesh(oneField, block1Part, &initValZero);
  stk::mesh::put_field_on_mesh(oneField, block2Part, &initValZero);

  std::string desc = "0,1,QUAD_4_2D,1,2,5,6,block_1\n"
                     "1,2,QUAD_4_2D,2,3,4,5,block_1\n"
                     "|dimension:2";
  setup_mesh(desc);

  stk::mesh::field_fill(1.0, oneField);

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 1);
    auto node2 = bulk->get_entity(stk::topology::NODE_RANK, 2);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    change_entity_parts_in_parallel({elem}, {block2Part.name()}, {block1Part.name()});

    check_entity_ownership({node2}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {block2Part.name()}, {block1Part.name()});
    check_device_bucket_part_membership({node2}, {block2Part.name()}, {block1Part.name()});
    check_device_bucket_part_membership({node5}, {block2Part.name()}, {block1Part.name()});

    check_field_data(node2, oneField);
    check_field_data(node5, oneField);
  } else {
    auto elem = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    auto node2 = bulk->get_entity(stk::topology::NODE_RANK, 2);
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    change_entity_parts_in_parallel({elem}, {block2Part.name()}, {block1Part.name()});

    check_entity_ownership({node2}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_device_bucket_part_membership({elem}, {block2Part.name()}, {block1Part.name()});
    check_device_bucket_part_membership({node2}, {block2Part.name()}, {block1Part.name()});
    check_device_bucket_part_membership({node5}, {block2Part.name()}, {block1Part.name()});

    check_field_data(node2, oneField);
    check_field_data(node5, oneField);
  }

}

NGP_TEST_F(DeviceMeshParallelMeshModTester_NoAura, check_change_entity_parts_ranked_removal_denial_three_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }
  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);
  declare_new_parts({"block_2"}, stk::topology::topology_t::HEX_8);

  auto elem2  = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto elem3  = bulk->get_entity(stk::topology::ELEM_RANK, 3);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8  = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node9  = bulk->get_entity(stk::topology::NODE_RANK, 9);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
  auto node13 = bulk->get_entity(stk::topology::NODE_RANK, 13);

  if (myRank == 1) {
    EXPECT_TRUE(is_owned({node9}));
    EXPECT_TRUE(is_shared({node9}));
    check_device_bucket_part_membership({node9},  {"block_1"}, {"block_2"});
    check_device_bucket_part_membership({node12}, {"block_1"}, {"block_2"});
  } else if (myRank == 2) {
    EXPECT_FALSE(is_owned({node9}));
    EXPECT_TRUE(is_shared({node9}));
    check_device_bucket_part_membership({node9},  {"block_1"}, {"block_2"});
    check_device_bucket_part_membership({node13}, {"block_1"}, {"block_2"});
  }

  if (myRank == 1) {
    change_entity_parts_in_parallel({elem2}, {"block_2"}, {"block_1"});
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  if (myRank == 0) {
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {});
    check_device_bucket_part_membership({node8}, {"block_1", "block_2"}, {});
  } else if (myRank == 1) {
    check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1"});
    check_device_bucket_part_membership({node9},  {"block_1", "block_2"}, {});
    check_device_bucket_part_membership({node12}, {"block_1", "block_2"}, {});
  } else {
    check_device_bucket_part_membership({elem3},  {"block_1"}, {"block_2"});
    check_device_bucket_part_membership({node13}, {"block_1"}, {"block_2"});
  }
}

class DeviceMeshParallelTester_WithAura : public DeviceMeshParallelTester
{
public:
  DeviceMeshParallelTester_WithAura()
    : DeviceMeshParallelTester(3, stk::mesh::BulkData::AUTO_AURA)
  {}
};

class DeviceMeshParallelMeshModTester_WithAura : public DeviceMeshParallelTester_WithAura
{
};

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, check_aura_ownerships)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"});

  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});

    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});

    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    check_entity_ownership({node12}, {Status::notOwned, Status::notShared, Status::isInAura});
  } else {
    auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
    check_entity_ownership({node5}, {Status::notOwned, Status::isShared, Status::notInAura});

    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::isInAura});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, abort_making_part_changes_entities_in_aura)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"});

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";

  if (rank == 0) {
    auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
    check_entity_ownership({node12}, {Status::notOwned, Status::notShared, Status::isInAura});

    EXPECT_ANY_THROW(change_entity_parts_in_parallel({node12}, {newPart}, {}));
  } else {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::isInAura});

    EXPECT_ANY_THROW(change_entity_parts_in_parallel({node1}, {newPart}, {}));
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_by_owning_proc_visible_to_aured_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";

  if (rank == 0) {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});

    EXPECT_NO_THROW(change_entity_parts_in_parallel({node1}, {newPart}, {}));
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});
  } else {
    auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::isInAura});

    change_entity_parts_in_parallel({}, {}, {});
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_remove_unranked_part_by_owning_proc_visible_to_aured_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);

  bulk->modification_begin();
  if (rank == 0) {
    bulk->change_entity_parts(node1, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});

  if (rank == 0) {
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({node1}, {}, {newPart});
  } else {
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({node1}, {"block_1"}, {newPart});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_to_owned_shared_and_aured_entity)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);

  if (myRank == 0) {
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    change_entity_parts_in_parallel({node5}, {newPart}, {});
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({node5}, {"block_1", newPart}, {});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_by_two_owners_visible_to_common_aured_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1", "new_part2"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node16 = bulk->get_entity(stk::topology::NODE_RANK, 16);

  if (myRank == 0) {
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({node1}, {"new_part1"}, {});
  } else if (myRank == 1) {
    check_entity_ownership({node1},  {Status::notOwned, Status::notShared, Status::isInAura});
    check_entity_ownership({node16}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  } else {
    check_entity_ownership({node16}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({node16}, {"new_part2"}, {});
  }
  Kokkos::fence();

  if (myRank == 0) {
    check_device_bucket_part_membership({node1}, {"block_1", "new_part1"}, {"new_part2"});
  } else if (myRank == 1) {
    check_device_bucket_part_membership({node1},  {"block_1", "new_part1"}, {"new_part2"});
    check_device_bucket_part_membership({node16}, {"block_1", "new_part2"}, {"new_part1"});
  } else {
    check_device_bucket_part_membership({node16}, {"block_1", "new_part2"}, {"new_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_not_needed_by_non_neighbor_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);

  if (myRank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({node1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({}, {}, {}));
  }
  Kokkos::fence();

  if (myRank == 0) {
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});
  } else if (myRank == 1) {
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});
  } else {
    check_bucket_count(stk::topology::NODE_RANK, 3u);
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_to_elem_visible_to_both_aured_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1\n" \
                     "3,4,HEX_8,13,14,15,16,17,18,19,20,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1", "new_part2"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto elem4 = bulk->get_entity(stk::topology::ELEM_RANK, 4);

  if (myRank == 1) {
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({elem2}, {"new_part1"}, {});
  } else if (myRank == 3) {
    check_entity_ownership({elem4}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({elem4}, {"new_part2"}, {});
  } else {
    check_entity_ownership({elem2}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  if (myRank == 0 || myRank == 1) {
    check_device_bucket_part_membership({elem2}, {"block_1", "new_part1"}, {"new_part2"});
  } else if (myRank == 2) {
    check_device_bucket_part_membership({elem2}, {"block_1", "new_part1"}, {"new_part2"});
    check_device_bucket_part_membership({elem4}, {"block_1", "new_part2"}, {"new_part1"});
  } else {
    check_device_bucket_part_membership({elem4}, {"block_1", "new_part2"}, {"new_part1"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_parts_to_unshared_elems_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_part1", "new_part2"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);

  if (myRank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem2}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({elem1}, {"new_part1"}, {});
  } else {
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::isInAura});
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({elem2}, {"new_part2"}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1", "new_part1"}, {"block_2", "new_part2"});
  check_device_bucket_part_membership({elem2}, {"block_2", "new_part2"}, {"block_1", "new_part1"});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {"new_part1", "new_part2"});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_unranked_part_to_shared_node_leaves_aura_entities_alone)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_part1"});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  auto& deviceMeshBefore = stk::mesh::get_updated_ngp_mesh(*bulk);
  auto numNodeBucketsBefore = deviceMeshBefore.get_device_bucket_repository().num_buckets(stk::topology::NODE_RANK);

  if (myRank == 0) {
    check_entity_ownership({node1},  {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5},  {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({node5}, {newPart}, {});
  } else {
    check_entity_ownership({node1},  {Status::notOwned, Status::notShared, Status::isInAura});
    check_entity_ownership({node5},  {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({node1},  {"block_1"}, {newPart});
  check_device_bucket_part_membership({node5}, {"block_1", newPart}, {});
  check_device_bucket_part_membership({node12}, {"block_1"}, {newPart});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_ranked_part_to_unshared_elem_induces_onto_aured_nodes)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";
  auto elem1  = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8  = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  if (myRank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({elem1}, {newPart}, {});
  } else {
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::isInAura});
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {"block_2"});
  check_device_bucket_part_membership({node1}, {"block_1", newPart}, {"block_2"});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2", newPart}, {});
  check_device_bucket_part_membership({node8}, {"block_1", "block_2", newPart}, {});
  check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newPart});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_ranked_parts_to_unshared_elems_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newPart1 = "new_elem_part1";
  std::string newPart2 = "new_elem_part2";
  declare_new_parts({newPart1, newPart2}, stk::topology::topology_t::HEX_8);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  auto elem1  = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2  = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8  = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  if (myRank == 0) {
    change_entity_parts_in_parallel({elem1}, {newPart1}, {});
  } else {
    change_entity_parts_in_parallel({elem2}, {newPart2}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1", newPart1}, {"block_2", newPart2});
  check_device_bucket_part_membership({elem2}, {"block_2", newPart2}, {"block_1", newPart1});
  check_device_bucket_part_membership({node1},  {"block_1", newPart1}, {"block_2", newPart2});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2", newPart1, newPart2}, {});
  check_device_bucket_part_membership({node8}, {"block_1", "block_2", newPart1, newPart2}, {});
  check_device_bucket_part_membership({node12}, {"block_2", newPart2}, {"block_1", newPart1});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_mixed_ranked_and_unranked_parts_to_unshared_elem)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  std::string newRankedPart = "new_elem_part";
  std::string newUnrankedPart = "new_part";
  declare_new_parts({newRankedPart}, stk::topology::topology_t::HEX_8);
  declare_new_parts({newUnrankedPart});

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  auto elem1  = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8  = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  if (myRank == 0) {
    change_entity_parts_in_parallel({elem1}, {newRankedPart, newUnrankedPart}, {});
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1", newRankedPart, newUnrankedPart}, {"block_2"});
  check_device_bucket_part_membership({node1}, {"block_1", newRankedPart}, {"block_2", newUnrankedPart});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
  check_device_bucket_part_membership({node8}, {"block_1", "block_2", newRankedPart}, {newUnrankedPart});
  check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newRankedPart, newUnrankedPart});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_remove_ranked_part_from_unshared_elem_uninduces_on_aured_nodes)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";
  auto elem1  = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8  = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  bulk->modification_begin();
  if (myRank == 0) {
    bulk->change_entity_parts(elem1, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {});
  check_device_bucket_part_membership({node1}, {"block_1", newPart}, {});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2", newPart}, {});

  if (myRank == 0) {
    change_entity_parts_in_parallel({elem1}, {}, {newPart});
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1"}, {newPart});
  check_device_bucket_part_membership({node1}, {"block_1"}, {newPart});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {newPart});
  check_device_bucket_part_membership({node8}, {"block_1", "block_2"}, {newPart});
  check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newPart});
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithAura, change_entity_parts_add_ranked_part_to_middle_elem_induces_onto_both_neighbors)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_mesh(desc);
  declare_new_parts({"new_elem_part1"}, stk::topology::topology_t::HEX_8);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";
  auto elem2  = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1  = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5  = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node9  = bulk->get_entity(stk::topology::NODE_RANK, 9);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);
  auto node16 = bulk->get_entity(stk::topology::NODE_RANK, 16);

  if (myRank == 0) {
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node9}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  } else if (myRank == 1) {
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});
    change_entity_parts_in_parallel({elem2}, {newPart}, {});
  } else {
    check_entity_ownership({node9}, {Status::notOwned, Status::isShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::notOwned, Status::notShared, Status::isInAura});
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem2},  {"block_1", newPart}, {});
  check_device_bucket_part_membership({node5},  {"block_1", newPart}, {});
  check_device_bucket_part_membership({node9},  {"block_1", newPart}, {});
  check_device_bucket_part_membership({node12}, {"block_1", newPart}, {});

  if (myRank != 2) {
    check_device_bucket_part_membership({node1}, {"block_1"}, {newPart});
  }
  if (myRank != 0) {
    check_device_bucket_part_membership({node16}, {"block_1"}, {newPart});
  }
}

class DeviceMeshParallelTester_WithGhosting : public DeviceMeshParallelTester
{
public:
  DeviceMeshParallelTester_WithGhosting()
    : DeviceMeshParallelTester(3, stk::mesh::BulkData::NO_AUTO_AURA)
  {}

  stk::mesh::Ghosting& create_custom_ghosting(std::string const& ghostingName,
                                              stk::mesh::EntityProcVec const& entitiesToGhost)
  {
    bulk->modification_begin();
    stk::mesh::Ghosting& ghosting = bulk->create_ghosting(ghostingName);
    bulk->change_ghosting(ghosting, entitiesToGhost);
    bulk->modification_end();

    return ghosting;
  }

  std::string ghosting_part_name(stk::mesh::Ghosting const& ghosting) const
  {
    return bulk->ghosting_part(ghosting).name();
  }

  void check_host_entity_ownership(stk::mesh::Entity entity, OwnershipStatus const& expectedStatus)
  {
    const stk::mesh::Bucket& hostBucket = bulk->bucket(entity);

    EXPECT_EQ(hostBucket.owned(), expectedStatus.owned  == Status::isOwned);
    EXPECT_EQ(hostBucket.shared(), expectedStatus.shared == Status::isShared);
    EXPECT_EQ(hostBucket.in_aura(), expectedStatus.inAura == Status::isInAura);
  }
};

class DeviceMeshParallelMeshModTester_WithGhosting : public DeviceMeshParallelTester_WithGhosting
{
};

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, check_custom_ghost_ownerships)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());

  stk::mesh::EntityProcVec entitiesToGhost;
  if (myRank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);

  if (myRank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});

    check_host_entity_ownership(elem1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node5, {Status::isOwned, Status::isShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1"}, {ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1"}, {ghostPart});
  } else {
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::notOwned, Status::isShared, Status::notInAura});

    check_host_entity_ownership(elem1, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node5, {Status::notOwned, Status::isShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1", ghostPart}, {});
    check_device_bucket_part_membership({node1}, {"block_1", ghostPart}, {});
    check_device_bucket_part_membership({node5}, {"block_1"}, {ghostPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, check_custom_ghost_ownerships_separate_blocks)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());

  stk::mesh::EntityProcVec entitiesToGhost;
  if (myRank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);

  if (myRank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::isOwned, Status::isShared, Status::notInAura});

    check_host_entity_ownership(elem1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node5, {Status::isOwned, Status::isShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1"}, {"block_2", ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", ghostPart});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {ghostPart});
  } else {
    auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node5}, {Status::notOwned, Status::isShared, Status::notInAura});

    check_host_entity_ownership(elem1, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(elem2, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node5, {Status::notOwned, Status::isShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1", ghostPart}, {"block_2"});
    check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1", ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1", ghostPart}, {"block_2"});
    check_device_bucket_part_membership({node5}, {"block_1", "block_2"}, {ghostPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, check_custom_ghost_ownerships_ghosted_by_both_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto myRank = stk::parallel_machine_rank(stk::parallel_machine_world());

  stk::mesh::EntityProcVec entitiesToGhost;
  if (myRank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  if (myRank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem2}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::notOwned, Status::notShared, Status::notInAura});

    check_host_entity_ownership(elem1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(elem2, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node12, {Status::notOwned, Status::notShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1"}, {"block_2", ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", ghostPart}, {"block_1"});
    check_device_bucket_part_membership({node12}, {"block_2", ghostPart}, {"block_1"});
  } else {
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node12}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({node1}, {Status::notOwned, Status::notShared, Status::notInAura});

    check_host_entity_ownership(elem2, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node12, {Status::isOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(elem1, {Status::notOwned, Status::notShared, Status::notInAura});
    check_host_entity_ownership(node1, {Status::notOwned, Status::notShared, Status::notInAura});

    check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1", ghostPart});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", ghostPart});
    check_device_bucket_part_membership({elem1}, {"block_1", ghostPart}, {"block_2"});
    check_device_bucket_part_membership({node1}, {"block_1", ghostPart}, {"block_2"});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_unranked_part_by_owning_proc_visible_to_ghosted_proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({"new_part1"});

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_TRUE(bulk->is_valid(elem1));

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  if (rank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {ghostPart});
  } else {
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_device_bucket_part_membership({elem1}, {"block_1", newPart, ghostPart}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_unranked_part_by_both_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart});

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {newPart}, {}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_entity_ownership({elem1}, {Status::isOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem2}, {Status::notOwned, Status::notShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", newPart, ghostPart}, {});
  } else {
    check_entity_ownership({elem1}, {Status::notOwned, Status::notShared, Status::notInAura});
    check_entity_ownership({elem2}, {Status::isOwned, Status::notShared, Status::notInAura});

    check_device_bucket_part_membership({elem1}, {"block_1", newPart, ghostPart}, {});
    check_device_bucket_part_membership({elem2}, {"block_2", newPart}, {ghostPart});

  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_and_remove_same_unranked_part_by_different_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart});

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);

  bulk->modification_begin();
  if (rank == 1) {
    bulk->change_entity_parts(elem2, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {}, {newPart}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", ghostPart}, {newPart});
  } else {
    check_device_bucket_part_membership({elem2}, {"block_2"}, {newPart, ghostPart});
    check_device_bucket_part_membership({elem1}, {"block_1", newPart, ghostPart}, {});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_remove_unranked_part_by_both_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart});

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {newPart}, {}));
  }
  Kokkos::fence();

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {}, {newPart}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {}, {newPart}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {"block_1"}, {newPart, ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", ghostPart}, {newPart});
  } else {
    check_device_bucket_part_membership({elem2}, {"block_2"}, {newPart, ghostPart});
    check_device_bucket_part_membership({elem1}, {"block_1", ghostPart}, {newPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_ranked_part_by_both_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart}, stk::topology::topology_t::HEX_8);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {newPart}, {}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {"block_2", ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", newPart, ghostPart}, {"block_1"});
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {"block_2", ghostPart});
    check_device_bucket_part_membership({node12}, {"block_2", newPart, ghostPart}, {"block_1"});
  } else {
    check_device_bucket_part_membership({elem1}, {"block_1", newPart, ghostPart}, {"block_2"});
    check_device_bucket_part_membership({elem2}, {"block_2", newPart}, {"block_1", ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1", newPart, ghostPart}, {"block_2"});
    check_device_bucket_part_membership({node12}, {"block_2", newPart}, {"block_1", ghostPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_and_remove_same_ranked_part_by_different_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart}, stk::topology::topology_t::HEX_8);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  bulk->modification_begin();
  if (rank == 1) {
    bulk->change_entity_parts(elem2, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {}, {newPart}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {"block_2", ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", ghostPart}, {"block_1", newPart});
    check_device_bucket_part_membership({node1}, {"block_1", newPart}, {"block_2", ghostPart});
    check_device_bucket_part_membership({node12}, {"block_2", ghostPart}, {"block_1", newPart});
  } else {
    check_device_bucket_part_membership({elem1}, {"block_1", newPart, ghostPart}, {"block_2"});
    check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1", newPart, ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1", newPart, ghostPart}, {"block_2"});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newPart, ghostPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_remove_ranked_part_by_both_owning_procs_ghosted_to_each_other)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  if (rank == 0) {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 1), 1);
  } else {
    entitiesToGhost.emplace_back(bulk->get_entity(stk::topology::ELEM_RANK, 2), 0);
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart}, stk::topology::topology_t::HEX_8);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node12 = bulk->get_entity(stk::topology::NODE_RANK, 12);

  bulk->modification_begin();
  if (rank == 0) {
    bulk->change_entity_parts(elem1, stk::mesh::PartVector{meta->get_part(newPart)});
  } else {
    bulk->change_entity_parts(elem2, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {}, {newPart}));
  } else {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem2}, {}, {newPart}));
  }
  Kokkos::fence();

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {"block_1"}, {"block_2", newPart, ghostPart});
    check_device_bucket_part_membership({elem2}, {"block_2", ghostPart}, {"block_1", newPart});
    check_device_bucket_part_membership({node1}, {"block_1"}, {"block_2", newPart, ghostPart});
    check_device_bucket_part_membership({node12}, {"block_2", ghostPart}, {"block_1", newPart});
  } else {
    check_device_bucket_part_membership({elem1}, {"block_1", ghostPart}, {"block_2", newPart});
    check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1", newPart, ghostPart});
    check_device_bucket_part_membership({node1}, {"block_1", ghostPart}, {"block_2", newPart});
    check_device_bucket_part_membership({node12}, {"block_2"}, {"block_1", newPart, ghostPart});
  }
}

NGP_TEST_F(DeviceMeshParallelMeshModTester_WithGhosting, change_entity_parts_add_and_remove_ranked_part_by_different_owning_procs_ghosted_to_all_procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 3) { GTEST_SKIP(); }

  std::string desc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n" \
                     "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n" \
                     "2,3,HEX_8,9,10,11,12,13,14,15,16,block_3";
  setup_mesh(desc);

  auto rank = stk::parallel_machine_rank(stk::parallel_machine_world());
  std::string newPart = "new_elem_part1";

  stk::mesh::EntityProcVec entitiesToGhost;
  auto ownedElem = bulk->get_entity(stk::topology::ELEM_RANK, rank + 1);
  for (int otherRank = 0; otherRank < 3; ++otherRank) {
    if (otherRank != static_cast<int>(rank)) {
      entitiesToGhost.emplace_back(ownedElem, otherRank);
    }
  }
  const stk::mesh::Ghosting& ghosting = create_custom_ghosting("customGhosting", entitiesToGhost);
  const std::string ghostPart = ghosting_part_name(ghosting);

  declare_new_parts({newPart}, stk::topology::topology_t::HEX_8);

  auto elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  auto elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);
  auto elem3 = bulk->get_entity(stk::topology::ELEM_RANK, 3);
  auto node1 = bulk->get_entity(stk::topology::NODE_RANK, 1);
  auto node5 = bulk->get_entity(stk::topology::NODE_RANK, 5);
  auto node8 = bulk->get_entity(stk::topology::NODE_RANK, 8);
  auto node9 = bulk->get_entity(stk::topology::NODE_RANK, 9);
  auto node13 = bulk->get_entity(stk::topology::NODE_RANK, 13);
  auto node16 = bulk->get_entity(stk::topology::NODE_RANK, 16);

  bulk->modification_begin();
  if (rank == 2) {
    bulk->change_entity_parts(elem3, stk::mesh::PartVector{meta->get_part(newPart)});
  }
  bulk->modification_end();

  if (rank == 0) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem1}, {newPart}, {}));
  } else if (rank == 2) {
    EXPECT_NO_THROW(change_entity_parts_in_parallel({elem3}, {}, {newPart}));
  } else {
    change_entity_parts_in_parallel({}, {}, {});
  }
  Kokkos::fence();

  check_device_bucket_part_membership({elem1}, {"block_1", newPart}, {"block_2", "block_3"});
  check_device_bucket_part_membership({elem2}, {"block_2"}, {"block_1", "block_3", newPart});
  check_device_bucket_part_membership({elem3}, {"block_3"}, {"block_1", "block_2", newPart});

  check_device_bucket_part_membership({node1}, {"block_1", newPart}, {"block_2", "block_3"});
  check_device_bucket_part_membership({node5}, {"block_1", "block_2", newPart}, {"block_3"});
  check_device_bucket_part_membership({node8}, {"block_1", "block_2", newPart}, {"block_3"});
  check_device_bucket_part_membership({node9}, {"block_2", "block_3"}, {"block_1", newPart});
  check_device_bucket_part_membership({node13}, {"block_3"}, {"block_1", "block_2", newPart});
  check_device_bucket_part_membership({node16}, {"block_3"}, {"block_1", "block_2", newPart});

  if (rank == 0) {
    check_device_bucket_part_membership({elem1}, {}, {ghostPart});
    check_device_bucket_part_membership({elem2, elem3}, {ghostPart}, {});
    check_device_bucket_part_membership({node1, node5, node8}, {}, {ghostPart});
    check_device_bucket_part_membership({node9, node13, node16}, {ghostPart}, {});
  } else if (rank == 1) {
    check_device_bucket_part_membership({elem1, elem3}, {ghostPart}, {});
    check_device_bucket_part_membership({elem2}, {}, {ghostPart});
    check_device_bucket_part_membership({node1, node13, node16}, {ghostPart}, {});
    check_device_bucket_part_membership({node5, node8, node9}, {}, {ghostPart});
  } else {
    check_device_bucket_part_membership({elem1, elem2}, {ghostPart}, {});
    check_device_bucket_part_membership({elem3}, {}, {ghostPart});
    check_device_bucket_part_membership({node1, node5, node8}, {ghostPart}, {});
    check_device_bucket_part_membership({node9, node13, node16}, {}, {ghostPart});
  }
}

#endif
