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
#include <stk_util/util/StridedArray.hpp>
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include <Kokkos_Core.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <string>

#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include "stk_mesh/baseImpl/NgpMeshHostData.hpp"

namespace stk {
namespace mesh {

template<typename NgpMemSpace> class DeviceMeshT;

template<typename BucketNgpMemSpace>
struct DeviceBucketT {
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  KOKKOS_FUNCTION
  DeviceBucketT()
    : m_owningMesh(nullptr),
      m_bucketSize(0),
      m_bucketCapacity(0),
      m_bucketId(0),
      m_bucketTopology(),
      m_entityRank(stk::topology::NODE_RANK)
  {}

  KOKKOS_FUNCTION
  unsigned bucket_id() const { return m_bucketId; }

  KOKKOS_FUNCTION
  size_t size() const { return m_bucketSize; }

  KOKKOS_FUNCTION
  size_t capacity() const { return m_bucketCapacity; }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank() const { return m_entityRank; }

  KOKKOS_FUNCTION
  stk::topology topology() const { return m_bucketTopology; }

  KOKKOS_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_INLINE_FUNCTION
  Permutations get_connected_permutations(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity operator[](unsigned offsetIntoBucket) const {
    return m_entities(offsetIntoBucket);
  }

  KOKKOS_FUNCTION
  bool member(stk::mesh::PartOrdinal partOrdinal) const
  {
    for(unsigned i=0; i<m_partOrdinals.size(); i++) {
      if(m_partOrdinals(i) == partOrdinal) {
        return true;
      }
    }
    return false;
  }

  void initialize_bucket_attributes(const stk::mesh::Bucket &bucket);
  void initialize_fixed_data_from_host(const stk::mesh::Bucket &bucket);
  void update_entity_data_from_host(const stk::mesh::Bucket &bucket);
  void update_sparse_connectivity_from_host(const stk::mesh::Bucket &bucket);

  void resize_device_views(const stk::mesh::Bucket &bucket);
  std::pair<unsigned, unsigned> scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket);

  EntityViewType m_entities;
  BucketConnectivityType m_nodeConnectivity;
  OrdinalViewType m_nodeConnectivityOffsets;

  OrdinalViewType m_nodeOrdinals;

  Unsigned2dViewType m_sparseConnectivityOffsets;
  BucketConnectivityType m_sparseConnectivity;
  OrdinalViewType m_sparseConnectivityOrdinals;
  PermutationViewType m_sparseConnectivityPermutations;

  PartOrdinalViewType m_partOrdinals;

  const stk::mesh::DeviceMeshT<BucketNgpMemSpace>* m_owningMesh;

  unsigned m_bucketSize;
  unsigned m_bucketCapacity;

  unsigned m_bucketId;
  stk::topology m_bucketTopology;
  stk::mesh::EntityRank m_entityRank;
};

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
      deviceMeshHostData(nullptr)
  {
    bulk->register_device_mesh();
    deviceMeshHostData = impl::get_ngp_mesh_host_data<NgpMemSpace>(*bulk);
    update_mesh();
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(DeviceMeshT &&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(DeviceMeshT &&) = default;

  KOKKOS_FUNCTION
  virtual ~DeviceMeshT() override {
    m_needSyncToHost = false;
    clear_buckets_and_views();
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
  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank,
                               const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return buckets[rank](meshIndex.bucket_id)[meshIndex.bucket_ord];
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    return buckets[rank](entityIndex.bucket_id).get_connected_entities(entityIndex.bucket_ord, connectedRank);
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return buckets[rank](entityIndex.bucket_id).get_nodes(entityIndex.bucket_ord);
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
    return buckets[rank](entityIndex.bucket_id).get_connected_ordinals(entityIndex.bucket_ord, connectedRank);
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
    return buckets[rank](entityIndex.bucket_id).get_connected_permutations(entityIndex.bucket_ord, connectedRank);
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
    return buckets[rank].size();
  }

  KOKKOS_FUNCTION
  const DeviceBucketT<NgpMemSpace> &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    buckets[rank](index).m_owningMesh = this;
    return buckets[rank](index);
  }

  KOKKOS_FUNCTION
  NgpCommMapIndicesT<NgpMemSpace> volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc) const
  {
    const size_t dataBegin = volatileFastSharedCommMapOffset[rank][proc];
    const size_t dataEnd   = volatileFastSharedCommMapOffset[rank][proc+1];
    NgpCommMapIndicesT<NgpMemSpace> buffer = Kokkos::subview(volatileFastSharedCommMap[rank], Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
    return buffer;
  }

  void clear()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
      buckets[rank] = BucketView();
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

    using HostEntitiesType = typename std::remove_reference<decltype(entities)>::type::HostMirror;
    using HostAddPartOrdinalsType = typename std::remove_reference<decltype(addPartOrdinals)>::type::HostMirror;
    using HostRemovePartOrdinalsType = typename std::remove_reference<decltype(removePartOrdinals)>::type::HostMirror;
  }

private:
  void set_entity_keys(const stk::mesh::BulkData& bulk_in);

  KOKKOS_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const
  {
    return (buckets[rank].use_count() == 1);
  }

  KOKKOS_FUNCTION
  void clear_buckets_and_views()
  {
    KOKKOS_IF_ON_HOST((
      if (is_last_bucket_reference()) {
        for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
          for (unsigned iBucket = 0; iBucket < buckets[rank].size(); ++iBucket) {
            buckets[rank][iBucket].~DeviceBucket();
          }
        }
      }
    ))
  }

  bool fill_buckets(const stk::mesh::BulkData& bulk_in);

  void fill_mesh_indices(const stk::mesh::BulkData& bulk_in);

  void copy_entity_keys_to_device();

  void copy_mesh_indices_to_device();

  void copy_sparse_connectivities_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();


  using BucketView = Kokkos::View<DeviceBucketT<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  stk::mesh::BulkData* bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;
  bool m_needSyncToHost;
  stk::mesh::EntityRank endRank;
  impl::NgpMeshHostData<NgpMemSpace>* deviceMeshHostData;

  EntityKeyViewType entityKeys;

  BucketView buckets[stk::topology::NUM_RANKS];
  HostMeshIndexType hostMeshIndices;
  MeshIndexType deviceMeshIndices;

  UnsignedViewType volatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];
  FastSharedCommMapViewType volatileFastSharedCommMap[stk::topology::NUM_RANKS];
};

using DeviceMesh = DeviceMeshT<stk::ngp::MemSpace>;

template<typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedEntities
DeviceBucketT<BucketNgpMemSpace>::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    const size_t nodeOffset = m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedEntities(&m_nodeConnectivity(nodeOffset), numNodes, 1);
  }

  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  return ConnectedEntities(&m_sparseConnectivity(offset), length, 1);
}

template<typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedOrdinals
DeviceBucketT<BucketNgpMemSpace>::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedOrdinals(m_nodeOrdinals.data(), numNodes, 1);
  }

  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  return ConnectedOrdinals(&m_sparseConnectivityOrdinals(offset), length, 1);
}

template<typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::Permutations
DeviceBucketT<BucketNgpMemSpace>::get_connected_permutations(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  if (m_sparseConnectivityPermutations.size() <= offset) {
    return Permutations(nullptr, 0);
  }

  return Permutations(&m_sparseConnectivityPermutations(offset), length, 1);
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::initialize_bucket_attributes(const stk::mesh::Bucket &bucket)
{
  m_bucketId = bucket.bucket_id();
  m_bucketCapacity = bucket.capacity();
  m_bucketSize = bucket.size();
  m_entityRank = bucket.entity_rank();
  m_bucketTopology = bucket.topology();
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::initialize_fixed_data_from_host(const stk::mesh::Bucket &bucket)
{
  const stk::mesh::PartVector& parts = bucket.supersets();
  m_partOrdinals = PartOrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"),
                                       parts.size());
  auto hostPartOrdinals = HostPartOrdinalViewType(bucket.superset_part_ordinals().first, parts.size());
  Kokkos::deep_copy(m_partOrdinals, hostPartOrdinals);
}

template<typename BucketNgpMemSpace>
std::pair<unsigned, unsigned>
DeviceBucketT<BucketNgpMemSpace>::scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket)
{
  if (bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
    unsigned maxNodesPerEntity = 0;
    unsigned totalNumConnectedNodes = 0;
    for (unsigned i = 0; i < bucket.size(); ++i) {
      maxNodesPerEntity = std::max(maxNodesPerEntity, bucket.num_nodes(i));
      totalNumConnectedNodes += bucket.num_nodes(i);
    }
    return std::make_pair(maxNodesPerEntity, totalNumConnectedNodes);
  }

  return std::make_pair<unsigned, unsigned>(bucket.topology().num_nodes(),
                                            bucket.topology().num_nodes() * m_bucketCapacity);
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::resize_device_views(const stk::mesh::Bucket & bucket)
{
  Kokkos::Profiling::pushRegion("resize_device_views()");

  Kokkos::Profiling::pushRegion("set node ordinals");

  const auto [maxNodesPerEntity, totalNumConnectedNodes] = scan_entities_for_nodal_connectivity(bucket);

  if (m_nodeOrdinals.size() != maxNodesPerEntity) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeOrdinals, static_cast<size_t>(maxNodesPerEntity));
    OrdinalViewType& nodeOrds = m_nodeOrdinals; //local var to avoid implicit this capture
    Kokkos::parallel_for(Kokkos::RangePolicy<stk::ngp::ExecSpace>(0, maxNodesPerEntity),
      KOKKOS_LAMBDA(const int i) {
        nodeOrds(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
      });
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("bucket entities");
  if (m_entities.size() != m_bucketCapacity) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_entities, m_bucketCapacity);
    STK_ThrowRequireMsg(m_bucketCapacity > 0, "bucket capacity must be greater than 0");
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("nodal connectivity");
  if (m_nodeConnectivity.size() != totalNumConnectedNodes) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeConnectivity, totalNumConnectedNodes);
  }

  if (m_nodeConnectivityOffsets.size() != m_bucketCapacity+1) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeConnectivityOffsets, m_bucketCapacity+1);
  }
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_entity_data_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_entity_data_from_host()");

  m_bucketSize = bucket.size();
  m_bucketCapacity = bucket.capacity();

  resize_device_views(bucket);

  Kokkos::Profiling::pushRegion("filling host-side Views");
  auto hostEntities = HostEntityViewType(bucket.begin(), m_bucketCapacity);
  auto hostNodeConnectivity = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeConnectivity);
  auto hostNodeConnectivityOffsets = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeConnectivityOffsets);
  unsigned nodeOffset = 0;
  for (unsigned iEntity = 0; iEntity < bucket.size(); ++iEntity) {
    const unsigned nodesPerEntity = bucket.num_nodes(iEntity);
    const stk::mesh::Entity * elemNodes = bucket.begin_nodes(iEntity);
    for (unsigned iNode = 0; iNode < nodesPerEntity; ++iNode) {
      hostNodeConnectivity(nodeOffset + iNode) = elemNodes[iNode];
    }
    hostNodeConnectivityOffsets(iEntity) = nodeOffset;
    nodeOffset += nodesPerEntity;
  }
  hostNodeConnectivityOffsets(bucket.size()) = nodeOffset;
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("deep_copy entities/connectivity/offsets");
  Kokkos::deep_copy(m_entities, hostEntities);
  Kokkos::deep_copy(m_nodeConnectivity, hostNodeConnectivity);
  Kokkos::deep_copy(m_nodeConnectivityOffsets, hostNodeConnectivityOffsets);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::popRegion();
}

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

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_sparse_connectivity_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_sparse_connectivity_from_host()");

  Unsigned2dViewType::HostMirror hostConnectivityOffsets("hostConnectivityOffsets", 0,0); 
  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityOffsets, stk::topology::NUM_RANKS, bucket.size()+1);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityOffsets, stk::topology::NUM_RANKS, bucket.size()+1);
  BucketConnectivityType::HostMirror hostConnectivity("hostConnectivity", 0);
  OrdinalViewType::HostMirror hostConnectivityOrdinals("hostConnectivityOrdinals", 0);
  PermutationViewType::HostMirror hostConnectivityPermutations("hostConnectivityPermutations", 0);

  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bucket.mesh().mesh_meta_data().entity_rank_count());

  unsigned offset = 0;
  for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; ++connectedRank) {
    for(unsigned i=0; i<bucket.size(); ++i) {
      hostConnectivityOffsets(connectedRank,i) = offset;
      offset += bucket.num_connectivity(i, connectedRank);
    }
    hostConnectivityOffsets(connectedRank,bucket.size()) = offset;
  }

  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivity, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivity, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityOrdinals, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityOrdinals, offset);
  const bool hasPermutations = bucket.has_permutation(stk::topology::EDGE_RANK)
                            || bucket.has_permutation(stk::topology::FACE_RANK)
                            || bucket.has_permutation(stk::topology::ELEM_RANK);
  if (hasPermutations) {
    Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityPermutations, offset);
    Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityPermutations, offset);
  }

  offset = 0;
  for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; ++connectedRank) {
    for(unsigned i=0; i<bucket.size(); ++i) {
      const unsigned numConn = bucket.num_connectivity(i, connectedRank);
      const Entity* conn = bucket.begin(i, connectedRank);
      const ConnectivityOrdinal* ords = bucket.begin_ordinals(i, connectedRank);
      const Permutation* perms = bucket.begin_permutations(i, connectedRank);

      for(unsigned nc=0; nc<numConn; ++nc) {
        hostConnectivity(offset+nc) = conn[nc];
        hostConnectivityOrdinals(offset+nc) = ords[nc];
        if (hasPermutations) {
          hostConnectivityPermutations(offset+nc) = perms[nc];
        }
      }
      offset += numConn;
    }
  }

  Kokkos::deep_copy(m_sparseConnectivityOffsets, hostConnectivityOffsets);
  Kokkos::deep_copy(m_sparseConnectivity, hostConnectivity);
  Kokkos::deep_copy(m_sparseConnectivityOrdinals, hostConnectivityOrdinals);
  if (hasPermutations) {
    Kokkos::deep_copy(m_sparseConnectivityPermutations, hostConnectivityPermutations);
  }

  Kokkos::Profiling::popRegion();
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
    copy_entity_keys_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("volatile-fast-shared-comm-map");
    copy_volatile_fast_shared_comm_map_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("mesh-indices");
    fill_mesh_indices(*bulk);
    copy_mesh_indices_to_device();
    Kokkos::Profiling::popRegion();

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
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    unsigned numStkBuckets = stkBuckets.size();

    BucketView bucketBuffer(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketBuffer"), numStkBuckets);

    if (numStkBuckets != buckets[rank].size()) {
      anyBucketChanges = true;
    }

    for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket) {
      stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
      const unsigned ngpBucketId = stkBucket.ngp_bucket_id();

      if (ngpBucketId == INVALID_BUCKET_ID) {
        Kokkos::Profiling::pushRegion("new bucket");
        // New bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucketT<NgpMemSpace>();
        bucketBuffer[iBucket].initialize_bucket_attributes(stkBucket);
        bucketBuffer[iBucket].initialize_fixed_data_from_host(stkBucket);
        bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
        bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
        anyBucketChanges = true;
        Kokkos::Profiling::popRegion();
      }
      else {
        Kokkos::Profiling::pushRegion("pre-existing bucket");
        // Pre-existing bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucketT(buckets[rank][ngpBucketId]);
        if (stkBucket.is_modified()) {
          bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
          bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
          anyBucketChanges = true;
        }
        bucketBuffer[iBucket].m_bucketId = stkBucket.bucket_id();
        Kokkos::Profiling::popRegion();
      }

      stkBucket.set_ngp_bucket_id(iBucket);
    }

    if (is_last_bucket_reference(rank)) {
      for (unsigned iBucket = 0; iBucket < buckets[rank].size(); ++iBucket) {
        buckets[rank][iBucket].~DeviceBucketT();
      }
    }

    buckets[rank] = bucketBuffer;
  }
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
void DeviceMeshT<NgpMemSpace>::fill_mesh_indices(const stk::mesh::BulkData& bulk_in)
{
  const size_t indexSpaceSize = bulk->get_size_of_entity_index_space();
  hostMeshIndices = HostMeshIndexType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "host_mesh_indices"), indexSpaceSize);

  for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
    const stk::mesh::BucketVector& bkts = bulk_in.buckets(rank);

    for(const stk::mesh::Bucket* bktptr : bkts)
    {
      const stk::mesh::Bucket& bkt = *bktptr;
      const unsigned bktId = bkt.bucket_id();
      for(unsigned i = 0; i < bkt.size(); ++i)
      {
        hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex{bktId, i};
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
void DeviceMeshT<NgpMemSpace>::copy_mesh_indices_to_device()
{
  unsigned length = hostMeshIndices.size();
  Kokkos::View<stk::mesh::FastMeshIndex*, NgpMemSpace> nonconst_device_mesh_indices(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_dev_mesh_indices"), length);
  Kokkos::deep_copy(nonconst_device_mesh_indices, hostMeshIndices);
  deviceMeshIndices = nonconst_device_mesh_indices;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_volatile_fast_shared_comm_map_to_device()
{
  bulk->volatile_fast_shared_comm_map<NgpMemSpace>(stk::topology::NODE_RANK, 0);
  auto& hostVolatileFastSharedCommMapOffset = deviceMeshHostData->hostVolatileFastSharedCommMapOffset;
  auto& hostVolatileFastSharedCommMap = deviceMeshHostData->hostVolatileFastSharedCommMap;

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank].extent(0));
    Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
  }
}
}
}

#endif

