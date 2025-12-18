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

#ifndef STK_MESH_DEVICEMESH_HPP
#define STK_MESH_DEVICEMESH_HPP

#include <stk_util/stk_config.h>
#include "Kokkos_Macros.hpp"
#include "View/Kokkos_ViewCtor.hpp"
#include "stk_mesh/base/DeviceFieldDataManagerBase.hpp"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "Kokkos_Core.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/DeviceFieldDataManager.hpp"
#include <string>

#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "stk_util/util/StkNgpVector.hpp"
#include "stk_util/util/ReportHandler.hpp"

#include "stk_mesh/baseImpl/DeviceMeshViewVector.hpp"
#include "stk_mesh/baseImpl/Partition.hpp"
#include "stk_mesh/baseImpl/NgpMeshHostData.hpp"
#include "stk_mesh/base/DeviceBucket.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"

namespace stk {
namespace mesh {

using DeviceBucket = DeviceBucketT<stk::ngp::MemSpace>;

template<typename NgpMemSpace>
class DeviceMeshT : public NgpMeshBase
{
public:
  typedef NgpMemSpace ngp_mem_space;

  static_assert(Kokkos::is_memory_space_v<NgpMemSpace>);
  using MeshExecSpace     =  typename NgpMemSpace::execution_space;
  using BucketType        = DeviceBucketT<NgpMemSpace>;
  using ConnectedNodes    = typename BucketType::ConnectedNodes;
  using ConnectedEntities = typename BucketType::ConnectedEntities;
  using ConnectedOrdinals = typename BucketType::ConnectedOrdinals;
  using Permutations      = typename BucketType::Permutations;
  using MeshIndex         = FastMeshIndex;

  KOKKOS_FUNCTION
  DeviceMeshT()
    : NgpMeshBase(),
      bulk(nullptr),
      spatial_dimension(0),
      lastBulkDataSynchronizedCount(0),
      synchronizedCount(0),
#ifndef STK_HIDE_DEPRECATED_CODE
      m_needSyncToHost(false),
#endif
      deviceMeshHostData(nullptr)
  {}

  explicit DeviceMeshT(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&const_cast<stk::mesh::BulkData&>(b)),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      lastBulkDataSynchronizedCount(0),
      synchronizedCount(0),
#ifndef STK_HIDE_DEPRECATED_CODE
      m_needSyncToHost(false),
#endif
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      deviceMeshHostData(nullptr),
      m_deviceBucketRepo(this, b.get_initial_bucket_capacity(), b.get_maximum_bucket_capacity()),
      m_deviceBufferOffsets(UnsignedViewType<NgpMemSpace>("deviceBufferOffsets", 1)),
      m_deviceMeshIndicesOffsets(UnsignedViewType<NgpMemSpace>("deviceMeshIndicesOffsets", 1))
  {
    bulk->register_device_mesh();
    deviceMeshHostData = impl::get_ngp_mesh_host_data<NgpMemSpace>(*bulk);
    update_mesh();

    deviceMeshHostData->m_hostBufferOffsets = Kokkos::create_mirror_view(m_deviceBufferOffsets);;
    deviceMeshHostData->m_hostMeshIndicesOffsets = Kokkos::create_mirror_view(m_deviceMeshIndicesOffsets);
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(DeviceMeshT &&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(DeviceMeshT &&) = default;

  KOKKOS_FUNCTION
  virtual ~DeviceMeshT() override {
#ifndef STK_HIDE_DEPRECATED_CODE
    m_needSyncToHost = false;
#endif
  }

  void update_mesh() override;

  void update_bulk_data() override;

  bool need_update_bulk_data() const override {
    return synchronizedCount > bulk->synchronized_count();
  }

  unsigned synchronized_count() const override { return synchronizedCount; }

  KOKKOS_FUNCTION
  unsigned get_spatial_dimension() const
  {
    return spatial_dimension;
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityId identifier(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].id();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].rank();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()];
  }

  KOKKOS_FUNCTION
  unsigned local_id(stk::mesh::Entity entity) const
  {
    return entityLocalIds[entity.local_offset()];
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank,
                               const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return m_deviceBucketRepo.m_buckets[rank][meshIndex.bucket_id][meshIndex.bucket_ord];
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    return m_deviceBucketRepo.m_buckets[rank][entityIndex.bucket_id].get_connected_entities(entityIndex.bucket_ord, connectedRank);
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return m_deviceBucketRepo.m_buckets[rank][entityIndex.bucket_id].get_nodes(entityIndex.bucket_ord);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    return m_deviceBucketRepo.m_buckets[rank][entityIndex.bucket_id].get_connected_ordinals(entityIndex.bucket_ord, connectedRank);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    return m_deviceBucketRepo.m_buckets[rank][entityIndex.bucket_id].get_connected_permutations(entityIndex.bucket_ord, connectedRank);
  }

  KOKKOS_FUNCTION
  Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
  {
    return device_mesh_index(entity);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
  {
    return deviceMeshIndices(entity.local_offset());
  }

  stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
  {
    return stk::mesh::get_bucket_ids(get_bulk_on_host(), rank, selector);
  }

  KOKKOS_FUNCTION
  EntityRank get_end_rank() const
  {
    return endRank;
  }

  KOKKOS_FUNCTION
  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return m_deviceBucketRepo.num_buckets(rank);
  }

  KOKKOS_FUNCTION
  const DeviceBucketT<NgpMemSpace> &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    return m_deviceBucketRepo.m_buckets[rank][index];
  }

  KOKKOS_FUNCTION
  NgpCommMapIndices<NgpMemSpace> volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc,
                                                               bool includeGhosts=false) const
  {
    const size_t dataBegin = volatileFastSharedCommMapOffset[rank][proc];
    const size_t dataEnd   = includeGhosts ? volatileFastSharedCommMapOffset[rank][proc+1]
                                           : dataBegin + volatileFastSharedCommMapNumShared[rank][proc];
    NgpCommMapIndices<NgpMemSpace> buffer = Kokkos::subview(volatileFastSharedCommMap[rank],
                                                            Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
    return buffer;
  }

  void clear()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
      m_deviceBucketRepo.m_buckets[rank] = BucketView();
  }

  stk::mesh::BulkData &get_bulk_on_host()
  {
    STK_ThrowRequireMsg(bulk != nullptr, "DeviceMesh::get_bulk_on_host, bulk==nullptr");
    return *bulk;
  }

  const stk::mesh::BulkData &get_bulk_on_host() const
  {
    STK_ThrowRequireMsg(bulk != nullptr, "DeviceMeshT::get_bulk_on_host, bulk==nullptr");
    return *bulk;
  }

  bool is_up_to_date() const {
    if(bulk == nullptr) { return false; }
    return synchronizedCount == bulk->synchronized_count();
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    using AddPartOrdinalsMemorySpace = typename std::remove_reference<decltype(addPartOrdinals)>::type::memory_space;
    using RemovePartOrdinalsMemorySpace = typename std::remove_reference<decltype(removePartOrdinals)>::type::memory_space;

    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, AddPartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'addPartOrdinals' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, RemovePartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'removePartOrdinals' View is inaccessible from the DeviceMesh execution space");

    bool hasRankedPart = impl::has_ranked_part(get_device_bucket_repository(), addPartOrdinals, removePartOrdinals);

    if (!hasRankedPart) {
      impl_batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
    } else {
      impl_batch_change_entity_parts_with_inducible_parts(entities, addPartOrdinals, removePartOrdinals);
    }

#ifndef STK_HIDE_DEPRECATED_CODE
    m_needSyncToHost = true;
#endif
    increment_synchronized_count();
  }

#ifndef STK_HIDE_DEPRECATED_CODE
  STK_DEPRECATED_MSG("Use update_bulk_data() instead.")
  void sync_to_host() {
    m_needSyncToHost = false;
  }

  STK_DEPRECATED_MSG("Use need_update_bulk_data() instead.")
  bool need_sync_to_host() const override {
    return m_needSyncToHost;
  }
#endif

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    using NewBucketsToAddViewType = Kokkos::View<impl::NumNewBucketsToAddPerPartition*, NgpMemSpace>;

#ifndef NDEBUG
    check_parts_are_not_internal(addPartOrdinals, removePartOrdinals);
#endif

    const unsigned maxCurrentNumPartsPerEntity = impl::get_max_num_parts_per_entity(*this, entities);
    const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*maxNewNumPartsPerEntity);
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(addPartOrdinals);

    // A proxy indices view to part ordinals: (rank, startPtr, length)
    using PartOrdinalsProxyViewType = Kokkos::View<impl::PartOrdinalsProxyIndices*>;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), entities.size());
    impl::set_new_part_list_per_entity(*this, entities, sortedAddPartOrdinals, removePartOrdinals,
                                       maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);

    PartOrdinalsProxyViewType copiedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "copiedPartOrdinalsProxy"), entities.size());
    Kokkos::deep_copy(copiedPartOrdinalsProxy, partOrdinalsProxy);

    impl::sort_and_unique_and_resize(partOrdinalsProxy, MeshExecSpace{});
    Kokkos::fence();

    m_deviceBucketRepo.batch_create_partitions(partOrdinalsProxy);

    using EntitySrcDestView = Kokkos::View<impl::EntitySrcDest*, NgpMemSpace>;
    EntitySrcDestView entitySrcDestView(Kokkos::view_alloc("srcDestPartitionIdPerEntity", Kokkos::WithoutInitializing), entities.size());

    m_deviceBucketRepo.batch_get_partitions(entities, copiedPartOrdinalsProxy, entitySrcDestView);

    NewBucketsToAddViewType numNewBucketsToAddInPartitions("NumNewBucketsToAddInPartitions", entities.size());
    impl::assign_dest_bucket_id_and_ordinal(*this, entitySrcDestView, numNewBucketsToAddInPartitions);

    {
      const MetaData& meta = get_bulk_on_host().mesh_meta_data();
      const FieldVector& fields = meta.get_fields();
      for(FieldBase* field : fields) {
        if (field->type_is<double>()) {
          field->data<double,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<float>()) {
          field->data<float,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<int>()) {
          field->data<int,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<unsigned>()) {
          field->data<unsigned,ReadWrite,stk::ngp::DeviceSpace>();
        }
      }
    }

    m_deviceBucketRepo.batch_create_buckets(numNewBucketsToAddInPartitions);

    set_dest_bucket_ids(*this, entitySrcDestView);

    unsigned maxNumBuckets = 0;
    for (auto rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank) {
      maxNumBuckets = std::max(maxNumBuckets, m_deviceBucketRepo.num_buckets(rank));
    }
    STK_ThrowAssert(maxNumBuckets != std::numeric_limits<unsigned>::max());

    using BucketConnectivitySize2DView = Kokkos::View<impl::BucketConnectivitySizes**, NgpMemSpace>;
    BucketConnectivitySize2DView numNewConnectivityViewSizesInBuckets("numNewConnectivityViewSizesInBuckets", static_cast<unsigned>(stk::topology::NUM_RANKS), maxNumBuckets);

    determine_connectivity_view_sizes_in_bucket(*this, entitySrcDestView, numNewConnectivityViewSizesInBuckets);

    m_deviceBucketRepo.batch_init_bucket_connectivity_views(entitySrcDestView, numNewConnectivityViewSizesInBuckets);

    m_deviceBucketRepo.batch_move_entities(entitySrcDestView);

    m_deviceBucketRepo.sync_from_partitions();

    Kokkos::fence();

    synchronizedCount++;
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts_with_inducible_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                                           const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                                           const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    using EntityWrapperViewType = Kokkos::View<impl::EntityWrapper*, NgpMemSpace>;
    using NewBucketsToAddViewType = Kokkos::View<impl::NumNewBucketsToAddPerPartition*, NgpMemSpace>;

#ifndef NDEBUG
    check_parts_are_not_internal(addPartOrdinals, removePartOrdinals);
#endif

    Kokkos::Profiling::pushRegion("construct_part_ordinal_proxy");
    auto maxNumDownwardConnectedEntities = impl::get_max_num_downward_connected_entities(*this, entities);
    auto entityInterval = maxNumDownwardConnectedEntities + 1;
    auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

    EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
    impl::populate_all_downward_connected_entities_and_wrap_entities(*this, entities, entityInterval, wrappedEntities);
    impl::remove_invalid_entities_sort_unique_and_resize(wrappedEntities, MeshExecSpace{});

    // determine resulting parts per entity including inducible parts
    const unsigned maxCurrentNumPartsPerEntity = impl::get_max_num_parts_per_entity(*this, wrappedEntities);
    const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(addPartOrdinals);

    // create part ordinals proxy, sort and unique it (and realloc)
    // A proxy indices view to part ordinals: (rank, startPtr, length)
    using PartOrdinalsProxyViewType = Kokkos::View<impl::PartOrdinalsProxyIndices*>;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
    impl::set_new_part_list_per_entity_with_induced_parts(*this, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                          maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);

    PartOrdinalsProxyViewType copiedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "copiedPartOrdinalsProxy"), wrappedEntities.size());
    Kokkos::deep_copy(copiedPartOrdinalsProxy, partOrdinalsProxy);

    impl::sort_and_unique_and_resize(partOrdinalsProxy, MeshExecSpace{});
    Kokkos::fence();
    Kokkos::Profiling::popRegion();

    m_deviceBucketRepo.batch_create_partitions(partOrdinalsProxy);

    using EntitySrcDestView = Kokkos::View<impl::EntitySrcDest*, NgpMemSpace>;
    EntitySrcDestView entitySrcDestView(Kokkos::view_alloc("srcDestPartitionIdPerEntity", Kokkos::WithoutInitializing), wrappedEntities.size());

    m_deviceBucketRepo.batch_get_partitions(wrappedEntities, copiedPartOrdinalsProxy, entitySrcDestView);

    NewBucketsToAddViewType numNewBucketsToAddInPartitions("NumNewBucketsToAddInPartitions", wrappedEntities.size());
    impl::assign_dest_bucket_id_and_ordinal(*this, entitySrcDestView, numNewBucketsToAddInPartitions);

    {
      const MetaData& meta = get_bulk_on_host().mesh_meta_data();
      const FieldVector& fields = meta.get_fields();
      for(FieldBase* field : fields) {
        if (field->type_is<double>()) {
          field->data<double,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<float>()) {
          field->data<float,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<int>()) {
          field->data<int,ReadWrite,stk::ngp::DeviceSpace>();
        }
        else if (field->type_is<unsigned>()) {
          field->data<unsigned,ReadWrite,stk::ngp::DeviceSpace>();
        }
      }
    }

    m_deviceBucketRepo.batch_create_buckets(numNewBucketsToAddInPartitions);

    set_dest_bucket_ids(*this, entitySrcDestView);

    unsigned maxNumBuckets = 0;
    for (auto rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank) {
      maxNumBuckets = std::max(maxNumBuckets, m_deviceBucketRepo.num_buckets(rank));
    }
    STK_ThrowAssert(maxNumBuckets != std::numeric_limits<unsigned>::max());

    using BucketConnectivitySize2DView = Kokkos::View<impl::BucketConnectivitySizes**, NgpMemSpace>;
    BucketConnectivitySize2DView numNewConnectivityViewSizesInBuckets("numNewConnectivityViewSizesInBuckets", stk::topology::NUM_RANKS, maxNumBuckets);

    determine_connectivity_view_sizes_in_bucket(*this, entitySrcDestView, numNewConnectivityViewSizesInBuckets);

    m_deviceBucketRepo.batch_init_bucket_connectivity_views(entitySrcDestView, numNewConnectivityViewSizesInBuckets);

    m_deviceBucketRepo.batch_move_entities(entitySrcDestView);

    m_deviceBucketRepo.sync_from_partitions();

    Kokkos::fence();

    synchronizedCount++;
  }

  MeshIndexType<NgpMemSpace>& get_fast_mesh_indices() {
    return deviceMeshIndices;
  }

  auto& get_ngp_parallel_sum_host_buffer_offsets() {
    return deviceMeshHostData->m_hostBufferOffsets;
  }

  auto& get_ngp_parallel_sum_host_mesh_indices_offsets() {
    return deviceMeshHostData->m_hostMeshIndicesOffsets;
  }

  auto& get_ngp_parallel_sum_device_mesh_indices_offsets() {
    return m_deviceMeshIndicesOffsets;
  }

  KOKKOS_INLINE_FUNCTION
  impl::DeviceBucketRepository<NgpMemSpace>& get_device_bucket_repository() {
    return m_deviceBucketRepo;
  }

  KOKKOS_INLINE_FUNCTION
  impl::DeviceBucketRepository<NgpMemSpace> const& get_device_bucket_repository() const {
    return m_deviceBucketRepo;
  }

  bool should_sort_buckets_by_first_entity_identifier() const {
    return bulk->should_sort_buckets_by_first_entity_identifier();
  }

  template <typename AddPartOrdinalsViewType, typename RemovePartOrdinalsViewType>
  void check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals, RemovePartOrdinalsViewType const& removePartOrdinals);

  DeviceFieldDataManagerBase* get_field_data_manager(const stk::mesh::BulkData& bulk_in);

  const DeviceFieldDataManagerBase* get_field_data_manager(const stk::mesh::BulkData& bulk_in) const;


private:
  bool fill_buckets(const stk::mesh::BulkData& bulk_in);

  void copy_entity_keys_to_device();

  void copy_entity_local_ids_to_device();

  void copy_sparse_connectivities_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();

  void update_field_data_manager();

  void update_last_internal_part_ordinal() {
    m_deviceBucketRepo.update_last_internal_part_ordinal();
  }

  void update_field_metadata_host_pointers();

  void copy_all_fields_to_device();

  void check_all_fields_to_synced_device();

  void set_all_synchronized_counts(unsigned count);

  void increment_synchronized_count() { ++synchronizedCount; }

  using BucketView = Kokkos::View<DeviceBucketT<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  stk::mesh::BulkData* bulk;
  unsigned spatial_dimension;
  unsigned lastBulkDataSynchronizedCount;
  unsigned synchronizedCount;

#ifndef STK_HIDE_DEPRECATED_CODE
  bool m_needSyncToHost;
#endif

  stk::mesh::EntityRank endRank;
  impl::NgpMeshHostData<NgpMemSpace>* deviceMeshHostData;

  EntityKeyViewType<NgpMemSpace> entityKeys;
  UnsignedViewType<NgpMemSpace> entityLocalIds;

  impl::DeviceBucketRepository<NgpMemSpace> m_deviceBucketRepo;
  HostMeshIndexType<NgpMemSpace> hostMeshIndices;
  MeshIndexType<NgpMemSpace> deviceMeshIndices;

  UnsignedViewType<NgpMemSpace> volatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];
  UnsignedViewType<NgpMemSpace> volatileFastSharedCommMapNumShared[stk::topology::NUM_RANKS];
  FastSharedCommMapViewType<NgpMemSpace> volatileFastSharedCommMap[stk::topology::NUM_RANKS];

  UnsignedViewType<NgpMemSpace> m_deviceBufferOffsets;
  UnsignedViewType<NgpMemSpace> m_deviceMeshIndicesOffsets;
};

using DeviceMesh = DeviceMeshT<stk::ngp::MemSpace>;

constexpr double RESIZE_FACTOR = 0.05;

template <typename DEVICE_VIEW, typename HOST_VIEW>
inline void reallocate_views(DEVICE_VIEW & deviceView, HOST_VIEW & hostView, size_t requiredSize, double resizeFactor = 0.0)
{
  const size_t currentSize = deviceView.extent(0);
  const size_t shrinkThreshold = currentSize - static_cast<size_t>(2*resizeFactor*currentSize);
  const bool needGrowth = (requiredSize > currentSize);
  const bool needShrink = (requiredSize < shrinkThreshold);

  if (needGrowth || needShrink) {
    const size_t newSize = requiredSize + static_cast<size_t>(resizeFactor*requiredSize);
    deviceView = DEVICE_VIEW(Kokkos::view_alloc(Kokkos::WithoutInitializing, deviceView.label()), newSize);
    hostView = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, deviceView);
  }
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_mesh()
{
  STK_ThrowRequireMsg(!bulk->m_meshModification.in_modifiable_state(),
                      "BulkData cannot be in a mod cycle during device to host synchronization");
  STK_ThrowRequireMsg(synchronizedCount == bulk->m_meshModification.last_device_synchronized_count(),
                      "DeviceMesh has been modified since last update from BulkData.  "
                      "Cannot update DeviceMesh from BulkData without syncing DeviceMesh to BulkData first.");
  if (is_up_to_date()) {
    return;
  }

  require_ngp_mesh_rank_limit(bulk->mesh_meta_data());

  Kokkos::Profiling::pushRegion("DeviceMeshT::update_mesh");
  const bool anyChanges = fill_buckets(*bulk);

  if (anyChanges) {
    Kokkos::Profiling::pushRegion("anyChanges stuff");

    Kokkos::Profiling::pushRegion("entity-keys");
    copy_entity_keys_to_device();
    copy_entity_local_ids_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("volatile-fast-shared-comm-map");
    copy_volatile_fast_shared_comm_map_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("mesh-indices");
    deviceMeshIndices = bulk->get_updated_fast_mesh_indices<NgpMemSpace>();
    Kokkos::Profiling::popRegion();

    update_field_data_manager();

    update_last_internal_part_ordinal();

    Kokkos::Profiling::popRegion();
  }

  set_all_synchronized_counts(bulk->synchronized_count());
  Kokkos::Profiling::popRegion();
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_bulk_data()
{
#ifndef STK_HIDE_DEPRECATED_CODE
  m_needSyncToHost = false;
#endif

  STK_ThrowRequireMsg(lastBulkDataSynchronizedCount == bulk->synchronized_count(),
                      "BulkData has been modified since last DeviceMesh update.  "
                      "Cannot update BulkData from DeviceMesh without syncing BulkData to DeviceMesh first.");

  if (synchronizedCount == bulk->synchronized_count())
  {
    return;
  }
  bulk->modification_begin_for_sync_to_host("sync DeviceMesh to host");
  check_all_fields_to_synced_device();
  m_deviceBucketRepo.sync_to_host(bulk->m_bucket_repository);
  Kokkos::deep_copy(bulk->m_entity_keys.get_view(), entityKeys);
  Kokkos::deep_copy(bulk->m_local_ids.get_view(), entityLocalIds);
  for (int i=0; i < stk::topology::NUM_RANKS; ++i)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMap[i], volatileFastSharedCommMap[i].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMapOffset[i], volatileFastSharedCommMapOffset[i].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMapNumShared[i], volatileFastSharedCommMapNumShared[i].extent(0));

    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMap[i], volatileFastSharedCommMap[i]);
    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMapOffset[i], volatileFastSharedCommMapOffset[i]);
    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMapNumShared[i], volatileFastSharedCommMapNumShared[i]);
  }

  set_all_synchronized_counts(synchronizedCount);
  bulk->modification_end_for_sync_to_host();

  for (EntityRank rank=stk::topology::BEGIN_RANK; rank < bulk->mesh_meta_data().entity_rank_count(); ++rank)
  {
    for (Bucket* bucket : bulk->buckets(rank))
    {
      for (unsigned i=0; i < bucket->size(); ++i)
      {
        bulk->set_mesh_index((*bucket)[i], bucket, i);
      }
    }
  }

  update_field_metadata_host_pointers();
}


template<typename NgpMemSpace>
bool DeviceMeshT<NgpMemSpace>::fill_buckets(const stk::mesh::BulkData& bulk_in)
{
  bool anyBucketChanges = false;

  Kokkos::Profiling::pushRegion("fill_buckets");
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    auto& hostBuckets = bulk_in.buckets(rank);
    auto& hostPartitions = bulk_in.m_bucket_repository.m_partitions[rank];
    m_deviceBucketRepo.copy_buckets_and_partitions_from_host(rank, hostBuckets, hostPartitions, anyBucketChanges);
  }
  Kokkos::Profiling::popRegion();

  return anyBucketChanges;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_keys_to_device()
{
  if (get_bulk_on_host().m_entity_keys.capacity() != entityKeys.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, entityKeys, get_bulk_on_host().m_entity_keys.capacity());
  }

  Kokkos::deep_copy(entityKeys, get_bulk_on_host().m_entity_keys.get_view());
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_local_ids_to_device()
{
  if (get_bulk_on_host().m_local_ids.capacity() != entityLocalIds.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, entityLocalIds, get_bulk_on_host().m_local_ids.capacity());
  }

  Kokkos::deep_copy(entityLocalIds, get_bulk_on_host().m_local_ids.get_view());
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_volatile_fast_shared_comm_map_to_device()
{
  bulk->volatile_fast_shared_comm_map<NgpMemSpace>(stk::topology::NODE_RANK, 0);
  auto& hostVolatileFastSharedCommMapOffset = deviceMeshHostData->hostVolatileFastSharedCommMapOffset;
  auto& hostVolatileFastSharedCommMapNumShared = deviceMeshHostData->hostVolatileFastSharedCommMapNumShared;
  auto& hostVolatileFastSharedCommMap = deviceMeshHostData->hostVolatileFastSharedCommMap;

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapNumShared[rank], hostVolatileFastSharedCommMapNumShared[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank].extent(0));
    Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMapNumShared[rank], hostVolatileFastSharedCommMapNumShared[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
  }
}

template<typename NgpMemSpace>
DeviceFieldDataManagerBase* DeviceMeshT<NgpMemSpace>::get_field_data_manager(const stk::mesh::BulkData& bulk_in)
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = nullptr;

  if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::DeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::DeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::UVMDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::UVMDeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::HostPinnedDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::HostPinnedDeviceSpace>();
  }
  else {
    STK_ThrowErrorMsg("Requested a DeviceFieldDataManager from a DeviceMesh with an unsupported MemorySpace: " <<
                      typeid(NgpMemSpace).name());
  }

  return deviceFieldDataManagerBase;
}

template<typename NgpMemSpace>
const DeviceFieldDataManagerBase* DeviceMeshT<NgpMemSpace>::get_field_data_manager(const stk::mesh::BulkData& bulk_in) const
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = nullptr;

  if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::DeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::DeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::UVMDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::UVMDeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::HostPinnedDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::HostPinnedDeviceSpace>();
  }
  else {
    STK_ThrowErrorMsg("Requested a DeviceFieldDataManager from a DeviceMesh with an unsupported MemorySpace: " <<
                      typeid(NgpMemSpace).name());
  }

  return deviceFieldDataManagerBase;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_field_data_manager()
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = get_field_data_manager(*bulk);
  STK_ThrowRequire(deviceFieldDataManagerBase != nullptr);
  deviceFieldDataManagerBase->update_all_bucket_allocations();
}

template<typename NgpMemSpace>
template <typename AddPartOrdinalsViewType, typename RemovePartOrdinalsViewType>
void DeviceMeshT<NgpMemSpace>::check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals,
                                                            RemovePartOrdinalsViewType const& removePartOrdinals)
{
  Kokkos::parallel_for(addPartOrdinals.extent(0),
    KOKKOS_CLASS_LAMBDA(const int i) {
      if (m_deviceBucketRepo.is_internal_part(addPartOrdinals(i)))
        Kokkos::abort("Cannot add an internal part.\n");
    }
  );

  Kokkos::parallel_for(removePartOrdinals.extent(0),
    KOKKOS_CLASS_LAMBDA(const int i) {
      if (m_deviceBucketRepo.is_internal_part(removePartOrdinals(i)))
        Kokkos::abort("Cannot remove an internal part.\n");
    }
  );
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_field_metadata_host_pointers()
{
  {
    for (FieldBase* field : bulk->mesh_meta_data().get_fields())
    {
      get_field_data_manager(*bulk)->update_host_bucket_pointers(field->mesh_meta_data_ordinal());
      field->modify_on_device();
    }
  }
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_all_fields_to_device()
{
  for (stk::mesh::FieldBase* field : bulk->mesh_meta_data().get_fields())
  {
    if (field->need_sync_to_device())
    {
      // this might be a problem: if built-in fields like coordinates are never
      // accessed on device, they still need a sync to device, but if the DeviceMesh
      // has changed then we can't do the sync
      field->sync_to_device();
    }

    field->modify_on_device();
  }
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::check_all_fields_to_synced_device()
{
  for (stk::mesh::FieldBase* field : bulk->mesh_meta_data().get_fields())
  {
    STK_ThrowRequireMsg(!field->need_sync_to_device(), "fields must have been synced to device prior to mesh modification");
  }
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::set_all_synchronized_counts(unsigned count)
{
  synchronizedCount             = count;
  lastBulkDataSynchronizedCount = count;
  bulk->m_meshModification.set_sync_count(count);
  bulk->m_meshModification.set_last_device_synchronized_count(count);
}

}
}

#endif

