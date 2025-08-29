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
#include "View/Kokkos_ViewCtor.hpp"
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
      synchronizedCount(0),
      m_needSyncToHost(false),
      deviceMeshHostData(nullptr)
  {}

  explicit DeviceMeshT(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&const_cast<stk::mesh::BulkData&>(b)),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      synchronizedCount(0),
      m_needSyncToHost(false),
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      deviceMeshHostData(nullptr),
      m_deviceBucketRepo(this),
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
    m_needSyncToHost = false;
  }

  void update_mesh() override;

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
    return m_deviceBucketRepo.m_buckets[rank](meshIndex.bucket_id)[meshIndex.bucket_ord];
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    return m_deviceBucketRepo.m_buckets[rank](entityIndex.bucket_id).get_connected_entities(entityIndex.bucket_ord, connectedRank);
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return m_deviceBucketRepo.m_buckets[rank](entityIndex.bucket_id).get_nodes(entityIndex.bucket_ord);
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
    return m_deviceBucketRepo.m_buckets[rank](entityIndex.bucket_id).get_connected_ordinals(entityIndex.bucket_ord, connectedRank);
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
    return m_deviceBucketRepo.m_buckets[rank](entityIndex.bucket_id).get_connected_permutations(entityIndex.bucket_ord, connectedRank);
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
  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return m_deviceBucketRepo.num_buckets(rank);
  }

  KOKKOS_FUNCTION
  const DeviceBucketT<NgpMemSpace> &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    m_deviceBucketRepo.m_buckets[rank](index).m_owningMesh = this;
    return m_deviceBucketRepo.m_buckets[rank](index);
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

  // This is an initial crude implementation that brings the device-side Views back to
  // the host and then kicks off a host-side mesh modification.  The modified host mesh
  // is then synchronized back to device.  This will not perform well and the semantics
  // are a little different from the final device-side capability (because the host mesh
  // will not be left in an unsynchronized state), but it can serve as a stand-in for
  // the final device-side mesh modification capability in the meantime.
  //
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

    using HostEntitiesType = typename std::remove_reference<decltype(entities)>::type::HostMirror;
    using HostAddPartOrdinalsType = typename std::remove_reference<decltype(addPartOrdinals)>::type::HostMirror;
    using HostRemovePartOrdinalsType = typename std::remove_reference<decltype(removePartOrdinals)>::type::HostMirror;

    HostEntitiesType copiedEntities = Kokkos::create_mirror_view(entities);
    HostAddPartOrdinalsType copiedAddPartOrdinals = Kokkos::create_mirror_view(addPartOrdinals);
    HostRemovePartOrdinalsType copiedRemovePartOrdinals = Kokkos::create_mirror_view(removePartOrdinals);

    Kokkos::deep_copy(copiedEntities, entities);
    Kokkos::deep_copy(copiedAddPartOrdinals, addPartOrdinals);
    Kokkos::deep_copy(copiedRemovePartOrdinals, removePartOrdinals);

    std::vector<stk::mesh::Entity> hostEntities;
    std::vector<stk::mesh::Part*> hostAddParts;
    std::vector<stk::mesh::Part*> hostRemoveParts;

    hostEntities.reserve(copiedEntities.extent(0));
    for (size_t i = 0; i < copiedEntities.extent(0); ++i) {
      hostEntities.push_back(copiedEntities[i]);
    }

    const stk::mesh::PartVector& parts = bulk->mesh_meta_data().get_parts();

    hostAddParts.reserve(copiedAddPartOrdinals.extent(0));
    for (size_t i = 0; i < copiedAddPartOrdinals.extent(0); ++i) {
      const size_t partOrdinal = copiedAddPartOrdinals[i];
      STK_ThrowRequire(partOrdinal < parts.size());
      hostAddParts.push_back(parts[partOrdinal]);
    }

    hostRemoveParts.reserve(copiedRemovePartOrdinals.extent(0));
    for (size_t i = 0; i < copiedRemovePartOrdinals.extent(0); ++i) {
      const size_t partOrdinal = copiedRemovePartOrdinals[i];
      STK_ThrowRequire(partOrdinal < parts.size());
      hostRemoveParts.push_back(parts[partOrdinal]);
    }

    m_needSyncToHost = false;
    bulk->batch_change_entity_parts(hostEntities, hostAddParts, hostRemoveParts);

    update_mesh();
    m_needSyncToHost = true;
  }

  // This function should be called before doing any host-side mesh operations after a
  // device-side mesh modification, to avoid accessing stale data.  Accessing the host
  // mesh without syncing it first should result in a throw.
  //
  void sync_to_host() {
    m_needSyncToHost = false;
  }

  // This can be used to check if the device-side mesh has been modified without
  // synchronizing it to the host.
  //
  bool need_sync_to_host() const override {
    return m_needSyncToHost;
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
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

    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    const unsigned maxCurrentNumPartsPerEntity = impl::get_max_num_parts_per_entity(*this,entities);
    const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();
    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*(1+maxNewNumPartsPerEntity));
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(addPartOrdinals);
    impl::set_new_part_list_per_entity(*this, entities, sortedAddPartOrdinals, removePartOrdinals,
                                       maxNewNumPartsPerEntity, newPartOrdinalsPerEntity);

    // TODO 
    // Need to check that parts are valid before using on device
    // Either check against host parts or
    // store and maintain a part ordinals list for lookup
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

  impl::DeviceBucketRepository<NgpMemSpace>& get_device_bucket_repository() {
    return m_deviceBucketRepo;
  }
  
  impl::DeviceBucketRepository<NgpMemSpace> const& get_device_bucket_repository() const {
    return m_deviceBucketRepo;
  }

  bool should_sort_buckets_by_first_entity_identifier() const {
    return bulk->should_sort_buckets_by_first_entity_identifier();
  }

private:
  void set_entity_keys(const stk::mesh::BulkData& bulk_in);

  void set_entity_local_ids(const stk::mesh::BulkData& bulk_in);

  bool fill_buckets(const stk::mesh::BulkData& bulk_in);

  void copy_entity_keys_to_device();

  void copy_entity_local_ids_to_device();

  void copy_sparse_connectivities_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();

  void update_field_data_manager(const stk::mesh::BulkData& bulk_in);

  using BucketView = Kokkos::View<DeviceBucketT<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  stk::mesh::BulkData* bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;
  bool m_needSyncToHost;
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
  if (is_up_to_date()) {
    return;
  }

  require_ngp_mesh_rank_limit(bulk->mesh_meta_data());

  Kokkos::Profiling::pushRegion("DeviceMeshT::update_mesh");
  const bool anyChanges = fill_buckets(*bulk);

  if (anyChanges) {
    Kokkos::Profiling::pushRegion("anyChanges stuff");

    Kokkos::Profiling::pushRegion("entity-keys");
    set_entity_keys(*bulk);
    set_entity_local_ids(*bulk);
    copy_entity_keys_to_device();
    copy_entity_local_ids_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("volatile-fast-shared-comm-map");
    copy_volatile_fast_shared_comm_map_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("mesh-indices");
    deviceMeshIndices = bulk->get_updated_fast_mesh_indices<NgpMemSpace>();
    Kokkos::Profiling::popRegion();

    update_field_data_manager(*bulk);

    Kokkos::Profiling::popRegion();
  }

  synchronizedCount = bulk->synchronized_count();
  Kokkos::Profiling::popRegion();
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
  m_deviceBucketRepo.sync_num_buckets_host_to_device();
  m_deviceBucketRepo.sync_num_partitions_host_to_device();
  Kokkos::Profiling::popRegion();

  return anyBucketChanges;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::set_entity_keys(const stk::mesh::BulkData& bulk_in)
{
  unsigned totalNumEntityKeys = bulk_in.get_size_of_entity_index_space();
  auto& hostEntityKeys = deviceMeshHostData->hostEntityKeys;

  reallocate_views(entityKeys, hostEntityKeys, totalNumEntityKeys, RESIZE_FACTOR);

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    for (unsigned i = 0; i < stkBuckets.size(); ++i) {
      const stk::mesh::Bucket & bucket = *stkBuckets[i];
      for (unsigned j = 0; j < bucket.size(); ++j) {
        stk::mesh::Entity entity = bucket[j];
        hostEntityKeys[entity.local_offset()] = bulk_in.entity_key(entity);
      }
    }
  }
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::set_entity_local_ids(const stk::mesh::BulkData& bulk_in)
{
  unsigned totalNumEntityLocalIds = bulk_in.get_size_of_entity_index_space();
  auto& hostEntityLocalIds = deviceMeshHostData->hostEntityLocalIds;

  reallocate_views(entityLocalIds, hostEntityLocalIds, totalNumEntityLocalIds, RESIZE_FACTOR);

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    for (unsigned i = 0; i < stkBuckets.size(); ++i) {
      const stk::mesh::Bucket & bucket = *stkBuckets[i];
      for (unsigned j = 0; j < bucket.size(); ++j) {
        stk::mesh::Entity entity = bucket[j];
        hostEntityLocalIds[entity.local_offset()] = bulk_in.local_id(entity);
      }
    }
  }
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_keys_to_device()
{
  auto& hostEntityKeys = deviceMeshHostData->hostEntityKeys;

  Kokkos::deep_copy(entityKeys, hostEntityKeys);
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_local_ids_to_device()
{
  auto& hostEntityLocalIds = deviceMeshHostData->hostEntityLocalIds;

  Kokkos::deep_copy(entityLocalIds, hostEntityLocalIds);
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
void DeviceMeshT<NgpMemSpace>::update_field_data_manager(const stk::mesh::BulkData& bulk_in)
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<NgpMemSpace>();
  STK_ThrowRequire(deviceFieldDataManagerBase != nullptr);

  deviceFieldDataManagerBase->update_all_bucket_allocations();
}

}
}

#endif

