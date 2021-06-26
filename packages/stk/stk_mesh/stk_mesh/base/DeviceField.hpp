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

#ifndef STK_MESH_DEVICE_FIELD_HPP
#define STK_MESH_DEVICE_FIELD_HPP

#include "stk_util/stk_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"
#include "stk_mesh/base/FieldSyncDebugging.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/EntityFieldData.hpp"
#include "stk_mesh/baseImpl/NgpFieldAux.hpp"
#include "stk_mesh/base/NgpFieldSyncDebugger.hpp"

namespace stk {
namespace mesh {

constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T, template <typename> class NgpDebugger>
class DeviceField : public NgpFieldBase
{
private:
  using FieldDataDeviceUnmanagedViewType = Kokkos::View<T***, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using StkDebugger = typename NgpDebugger<T>::StkFieldSyncDebuggerType;

public:
  using value_type = T;

  KOKKOS_FUNCTION
  DeviceField()
    : NgpFieldBase(),
      rank(stk::topology::NODE_RANK),
      ordinal(INVALID_ORDINAL),
      hostBulk(nullptr),
      hostField(nullptr),
      bucketCapacity(0),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(0),
      userSpecifiedSelector(false),
      syncSelector(nullptr),
      asyncCopyState(impl::AsyncCopyState()),
      fieldSyncDebugger(nullptr)
  {
    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
  }

  DeviceField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      rank(stkField.entity_rank()),
      ordinal(stkField.mesh_meta_data_ordinal()),
      hostBulk(&bulk),
      hostField(&stkField),
      bucketCapacity(0),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(bulk.synchronized_count()),
      copyCounter("copy_counter"),
      userSpecifiedSelector(false),
      syncSelector(new Selector()),
      asyncCopyState(impl::AsyncCopyState()),
      fieldSyncDebugger(nullptr)
  {
    ThrowRequireMsg(isFromGetUpdatedNgpField, "NgpField must be obtained from get_updated_ngp_field()");
    initialize();
    update_field();

    fieldSyncDebugger.initialize_debug_views(this);
  }

  DeviceField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField,
              const ExecSpace& execSpace, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      rank(stkField.entity_rank()),
      ordinal(stkField.mesh_meta_data_ordinal()),
      hostBulk(&bulk),
      hostField(&stkField),
      bucketCapacity(0),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(bulk.synchronized_count()),
      copyCounter("copy_counter"),
      userSpecifiedSelector(false),
      syncSelector(new Selector()),
      asyncCopyState(impl::AsyncCopyState()),
      fieldSyncDebugger(nullptr)
  {
    ThrowRequireMsg(isFromGetUpdatedNgpField, "NgpField must be obtained from get_updated_ngp_field()");
    asyncCopyState.set_state(execSpace, impl::INVALID);
    initialize();
    update_field();

    fieldSyncDebugger.initialize_debug_views(this);
  }

  void initialize()
  {
    hostField->template make_field_sync_debugger<StkDebugger>();
    fieldSyncDebugger = NgpDebugger<T>(&hostField->get_field_sync_debugger<StkDebugger>());

    needSyncToHost = Kokkos::View<bool, HostExecSpace>("needSyncToHost");
    needSyncToDevice = Kokkos::View<bool, HostExecSpace>("needSyncToDevice");

    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
  }

  void set_field_states(DeviceField<T, NgpDebugger>* fields[])
  {
    const unsigned numStates = hostField->number_of_states();
    for (unsigned state = 0; state < numStates; ++state) {
      stateFields[state] = fields[state];
    }
  }

  void set_execution_space(const ExecSpace& execSpace) override
  {
    asyncCopyState.set_execution_space(execSpace);
  }

  void update_field(bool needToSyncAllDataToDevice = false) override
  {
    ProfilingBlock prof("update_field for " + hostField->name());
    stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
    auto hostFieldEntityRank = hostField->entity_rank();
    const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostFieldEntityRank, selector);
    const stk::mesh::BucketVector& allBuckets = hostBulk->buckets(hostFieldEntityRank);
    numBucketsForField = buckets.size();
    maxNumScalarsPerEntity = hostField->max_size(rank);

    if (!buckets.empty()) {
      bucketCapacity = buckets[0]->capacity();
    }

    construct_field_exist_view(allBuckets, selector);
    construct_all_fields_buckets_num_components_per_entity_view(allBuckets);
    construct_field_buckets_num_components_per_entity_view(buckets);
    construct_bucket_sizes_view(buckets);
    construct_new_index_view(allBuckets);
    if (numBucketsForField != deviceData.extent(0)) {
      construct_view(buckets, "deviceData_"+hostField->name(), maxNumScalarsPerEntity);
    } else {
      move_unmodified_buckets(buckets, maxNumScalarsPerEntity);
    }

    if(needToSyncAllDataToDevice) {
      clear_sync_state_flags();
      set_modify_on_host();
      copy_host_to_device();
    } else {
      copy_new_and_modified_buckets_from_host(buckets, maxNumScalarsPerEntity);
    }

    fieldSyncDebugger.update_field(this, needToSyncAllDataToDevice);

    for(auto * bucket : allBuckets)
    {
      bucket->set_ngp_field_bucket_id(get_ordinal(), INVALID_BUCKET_ID);
    }
    for(auto * bucket : buckets)
    {
      bucket->set_ngp_field_bucket_id(get_ordinal(), bucket->bucket_id());
    }

    synchronizedCount = hostBulk->synchronized_count();
    hostSelectedBucketOffset = newHostSelectedBucketOffset;
    deviceSelectedBucketOffset = newDeviceSelectedBucketOffset;
    reset_sync_selector();
  }

  size_t num_syncs_to_host() const override { return hostField->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return hostField->num_syncs_to_device(); }

  void fence() override {
    asyncCopyState.execSpace.fence();

    if(asyncCopyState.syncMode == impl::DEVICE_TO_HOST_ASYNC) {
      Selector selector = selectField(*hostField);
      copy_host_data_to_stk_field_data(selector);
      reset_sync_selector();
    }

    asyncCopyState.reset_state();
  }

  void modify_on_host() override
  {
    set_modify_on_host();
    userSpecifiedSelector = false;
    *syncSelector = Selector(hostBulk->mesh_meta_data().universal_part());
  }

  void modify_on_host(const Selector& selector) override
  {
    set_modify_on_host();
    userSpecifiedSelector = true;
    *syncSelector |= selector;
  }

  void modify_on_device() override
  {
    set_modify_on_device();
    userSpecifiedSelector = false;
    *syncSelector = Selector(hostBulk->mesh_meta_data().universal_part());
  }

  void modify_on_device(const Selector& selector) override
  {
    set_modify_on_device();
    userSpecifiedSelector = true;
    *syncSelector |= selector;
  }

  void clear_sync_state() override
  {
    fieldSyncDebugger.clear_sync_state(this);
    clear_sync_state_flags();
    reset_sync_selector();
  }

  void clear_host_sync_state() override
  {
    fieldSyncDebugger.clear_host_sync_state(this);
    if (need_sync_to_device()) {
      clear_sync_state_flags();
      reset_sync_selector();
    }
  }

  void clear_device_sync_state() override
  {
    fieldSyncDebugger.clear_device_sync_state(this);
    if (need_sync_to_host()) {
      clear_sync_state_flags();
      reset_sync_selector();
    }
  }

  void sync_to_host() override
  {
    asyncCopyState.set_state(Kokkos::DefaultExecutionSpace(), impl::DEVICE_TO_HOST);
    internal_sync_to_host();
    Kokkos::fence();
    asyncCopyState.reset_state();
  }

  void sync_to_host(const ExecSpace& newExecSpace) override
  {
    asyncCopyState.set_state(newExecSpace, impl::DEVICE_TO_HOST_ASYNC);
    internal_sync_to_host();
  }

  void sync_to_device() override
  {
    asyncCopyState.set_state(Kokkos::DefaultExecutionSpace(), impl::HOST_TO_DEVICE);
    internal_sync_to_device();
    Kokkos::fence();
    asyncCopyState.reset_state();
  }

  void sync_to_device(const ExecSpace& newExecSpace) override
  {
    asyncCopyState.set_state(newExecSpace, impl::HOST_TO_DEVICE_ASYNC);
    internal_sync_to_device();
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  KOKKOS_DEFAULTED_FUNCTION DeviceField(const DeviceField &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T, NgpDebugger>& operator=(const DeviceField<T, NgpDebugger>&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T, NgpDebugger>& operator=(DeviceField<T, NgpDebugger>&&) = default;

  KOKKOS_FUNCTION virtual ~DeviceField()
  {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    if (is_last_field_copy()) {
      delete syncSelector;
    }
#endif
  }

  unsigned get_component_stride() const
  {
    unsigned stride = 1;
#ifdef KOKKOS_ENABLE_CUDA
    stride = bucketCapacity;
#endif
    return stride;
  }

  KOKKOS_FUNCTION
  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entityIndex) const {
    const unsigned bucketId = entityIndex.bucket_id;
    return deviceAllFieldsBucketsNumComponentsPerEntity[bucketId];
  }

  unsigned debug_get_bucket_offset(unsigned bucketOrdinal) const override {
    return fieldSyncDebugger.get_bucket_offset(bucketOrdinal);
  }

  template <typename Mesh> KOKKOS_FUNCTION
  T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component, fileName, lineNumber);
  }

  KOKKOS_FUNCTION
  T& get(stk::mesh::FastMeshIndex index, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& get(MeshIndex index, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket->bucket_id()), ORDER_INDICES(index.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  T& operator()(const stk::mesh::FastMeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& operator()(const MeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket->bucket_id()), ORDER_INDICES(index.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index,
                                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, fileName, lineNumber);
    T* dataPtr = &deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, 0));
    const unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_component_stride());
  }

  template <typename Mesh>
  void set_all(const Mesh& ngpMesh, const T& value,
               const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER)
  {
    fieldSyncDebugger.device_stale_access_check(this, fileName, lineNumber);
    Kokkos::deep_copy(deviceData, value);
    modify_on_device();
    fieldSyncDebugger.set_any_potential_device_field_modification(true);
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank get_rank() const { return rank; }

  KOKKOS_FUNCTION
  unsigned get_ordinal() const { return ordinal; }

  const stk::mesh::BulkData& get_bulk() const { return *hostBulk; }

  FieldState state() const { return hostField->state(); }

  void rotate_multistate_data() override
  {
  }

  KOKKOS_FUNCTION
  void swap(DeviceField<T, NgpDebugger> &other)
  {
    swap_views(hostData,   other.hostData);
    swap_views(deviceData, other.deviceData);
    swap_views(needSyncToHost, other.needSyncToHost);
    swap_views(needSyncToDevice, other.needSyncToDevice);
  }

protected:
  bool need_sync_to_host() const
  {
    return needSyncToHost();
  }

  bool need_sync_to_device() const
  {
    return needSyncToDevice();
  }

  unsigned get_contiguous_bucket_offset_end(const BucketVector& buckets, unsigned i)
  {
    unsigned prevOffset = newHostSelectedBucketOffset(buckets[i]->bucket_id());
    unsigned endIndex = i;

    for(unsigned bucketIdx = i+1; bucketIdx < buckets.size(); bucketIdx++) {
       unsigned currOffset = newHostSelectedBucketOffset(buckets[bucketIdx]->bucket_id());
       if(currOffset == prevOffset+1) {
         prevOffset = currOffset;
         endIndex = bucketIdx;
       } else {
         break;
       }
    }

    return endIndex;
  }

  void debug_modification_begin() override
  {
    fieldSyncDebugger.detect_device_field_modification(this);
    hostField->get_field_sync_debugger<StkDebugger>().fill_last_mod_location_field_from_device();
  }

  void debug_modification_end(size_t synchronizationCount) override
  {
    hostField->get_field_sync_debugger<StkDebugger>().clear_last_field_value();
    fieldSyncDebugger.update_debug_storage(this, synchronizationCount);
    hostField->get_field_sync_debugger<StkDebugger>().fill_last_mod_location_view_from_host();
  }

  void debug_detect_device_field_modification() override
  {
    fieldSyncDebugger.detect_device_field_modification(this);
  }

private:

  STK_FUNCTION
  bool is_last_field_copy() const
  {
    return (copyCounter.use_count() == 1);
  }

  void reset_sync_selector()
  {
    userSpecifiedSelector = false;
    *syncSelector = Selector();
  }

  void set_modify_on_host()
  {
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    ThrowRequireMsg(needSyncToHost() == false, "field " << hostField->name());
#endif
    needSyncToDevice() = true;
  }

  void set_modify_on_device()
  {
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    ThrowRequireMsg(needSyncToDevice() == false, "field " << hostField->name());
#endif
    needSyncToHost() = true;
  }

  void clear_sync_state_flags()
  {
    needSyncToHost() = false;
    needSyncToDevice() = false;
  }

  void move_unmodified_buckets_in_range_from_back(const BucketVector& buckets, unsigned numPerEntity, unsigned& currBaseIndex)
  {
    int startIndex = currBaseIndex;
    int endIndex = buckets.size() - 1;
    unsigned oldBucketId, newBucketId;
    unsigned oldBucketOffset, newBucketOffset;

    for(unsigned j = currBaseIndex; j < buckets.size(); j++) {
      oldBucketId = buckets[j]->get_ngp_field_bucket_id(get_ordinal());
      newBucketId = buckets[j]->bucket_id();

      if(oldBucketId == INVALID_BUCKET_ID) {
        endIndex = j - 1;
        break;
      }

      oldBucketOffset = hostSelectedBucketOffset(oldBucketId);
      newBucketOffset = newHostSelectedBucketOffset(newBucketId);

      if(oldBucketOffset >= newBucketOffset || buckets[j]->get_ngp_field_bucket_is_modified(get_ordinal())) {
        endIndex = j - 1;
        break;
      }
    }

    for(int j = endIndex; j >= startIndex; j--) {
      oldBucketId = buckets[j]->get_ngp_field_bucket_id(get_ordinal());
      newBucketId = buckets[j]->bucket_id();

      copy_moved_host_bucket_data<FieldDataHostViewType<T>, UnmanagedHostInnerView<T>>(hostData, hostData, oldBucketId, newBucketId, numPerEntity);
      copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(deviceData, deviceData, oldBucketId, newBucketId, numPerEntity);
    }
    currBaseIndex = endIndex;
  }

  void move_unmodified_buckets(const BucketVector& buckets, unsigned numPerEntity)
  {
    for(unsigned i = 0; i < buckets.size(); i++) {
      if(hostFieldExistsOnBucket(buckets[i]->bucket_id())) {
        unsigned oldBucketId = buckets[i]->get_ngp_field_bucket_id(get_ordinal());
        unsigned newBucketId = buckets[i]->bucket_id();

        if(oldBucketId == INVALID_BUCKET_ID) { continue; }

        unsigned oldBucketOffset = hostSelectedBucketOffset(oldBucketId);
        unsigned newBucketOffset = newHostSelectedBucketOffset(newBucketId);

        if(oldBucketOffset != newBucketOffset && !buckets[i]->get_ngp_field_bucket_is_modified(get_ordinal())) {
          if(oldBucketOffset > newBucketOffset) {
            copy_moved_host_bucket_data<FieldDataHostViewType<T>, UnmanagedHostInnerView<T>>(hostData, hostData, oldBucketId, newBucketId, numPerEntity);
            copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(deviceData, deviceData, oldBucketId, newBucketId, numPerEntity);
          } else {
            move_unmodified_buckets_in_range_from_back(buckets, numPerEntity, i);
          }
        }
      }
    }
  }

  void copy_unmodified_buckets(const BucketVector& buckets, FieldDataDeviceViewType<T> destDevView, FieldDataHostViewType<T> destHostView, unsigned numPerEntity)
  {
    for(unsigned i = 0; i < buckets.size(); i++) {
      if(hostFieldExistsOnBucket(buckets[i]->bucket_id())) {
        unsigned oldBucketId = buckets[i]->get_ngp_field_bucket_id(get_ordinal());
        unsigned newBucketId = buckets[i]->bucket_id();

        if(!buckets[i]->get_ngp_field_bucket_is_modified(get_ordinal())) {
          ThrowRequire(hostData.extent(0) != 0 && hostSelectedBucketOffset.extent(0) != 0);
          copy_moved_host_bucket_data<FieldDataHostViewType<T>, UnmanagedHostInnerView<T>>(destHostView, hostData, oldBucketId, newBucketId, numPerEntity);
          copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(destDevView, deviceData, oldBucketId, newBucketId, numPerEntity);
        }
      }
    }
  }

  void fill_host_view(Bucket* bucket, UnmanagedHostInnerView<T>& unHostInnerView, unsigned numPerEntity)
  {

    T* hostFieldPtr = static_cast<T*>(stk::mesh::field_data<FieldBase, EmptyStkFieldSyncDebugger>(*hostField, *bucket));
    unsigned numEntities = bucket->size();

    for(unsigned i = 0; i < numEntities; i++) {
      for(unsigned j = 0; j < numPerEntity; j++) {
        unHostInnerView(i,j) = hostFieldPtr[i * numPerEntity + j];
      }
    }
  }

  template<typename DataView, typename UnmanagedView>
  void copy_moved_host_bucket_data(DataView destView, DataView srcView, unsigned oldBucketId, unsigned newBucketId, unsigned numPerEntity)
  {
    unsigned oldSelectedOffset = hostSelectedBucketOffset(oldBucketId);
    unsigned newSelectedOffset = newHostSelectedBucketOffset(newBucketId);

    T* srcPtr = srcView.data() + oldSelectedOffset * bucketCapacity * numPerEntity;
    T* destPtr = destView.data() + newSelectedOffset * bucketCapacity * numPerEntity;

    UnmanagedView unInnerSrcView(srcPtr, bucketCapacity, numPerEntity);
    UnmanagedView unInnerDestView(destPtr, bucketCapacity, numPerEntity);
    Kokkos::deep_copy(asyncCopyState.execSpace, unInnerDestView, unInnerSrcView);
  }

  template<typename DataView, typename UnmanagedView>
  void copy_moved_device_bucket_data(DataView destView, DataView srcView, unsigned oldBucketId, unsigned newBucketId, unsigned numPerEntity)
  {
    unsigned oldSelectedOffset = hostSelectedBucketOffset(oldBucketId);
    unsigned newSelectedOffset = newHostSelectedBucketOffset(newBucketId);

    T* srcPtr = srcView.data() + oldSelectedOffset * bucketCapacity * numPerEntity;
    T* destPtr = destView.data() + newSelectedOffset * bucketCapacity * numPerEntity;

    UnmanagedView unInnerSrcView(srcPtr, ORDER_INDICES(bucketCapacity, numPerEntity));
    UnmanagedView unInnerDestView(destPtr, ORDER_INDICES(bucketCapacity, numPerEntity));
    Kokkos::deep_copy(asyncCopyState.execSpace, unInnerDestView, unInnerSrcView);
  }

  void copy_bucket_from_device_to_host(Bucket* bucket, unsigned numPerEntity, unsigned numContiguousBuckets = 1)
  {
    unsigned selectedBucketOffset = newHostSelectedBucketOffset(bucket->bucket_id());

    T* devicePtr = deviceData.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unDeviceInnerView(devicePtr, ORDER_INDICES(bucketCapacity*numContiguousBuckets, numPerEntity));

    T* bufferPtr = bufferData.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unBufferInnerView(bufferPtr, bucketCapacity*numContiguousBuckets, numPerEntity);

    transpose_contiguous_device_data_into_buffer(asyncCopyState.execSpace, bucketCapacity*numContiguousBuckets, numPerEntity, unDeviceInnerView, unBufferInnerView);

    UnmanagedHostInnerView<T> unHostInnerView(&hostData(selectedBucketOffset,0,0),
                                              bucketCapacity*numContiguousBuckets, numPerEntity);
    Kokkos::deep_copy(asyncCopyState.execSpace, unHostInnerView, unBufferInnerView);
  }

  void copy_bucket_from_host_to_device(Bucket* bucket, unsigned maxNumPerEntity, unsigned numContiguousBuckets = 1)
  {
    unsigned selectedBucketOffset = newHostSelectedBucketOffset(bucket->bucket_id());

    UnmanagedHostInnerView<T> unHostInnerView(&hostData(selectedBucketOffset,0,0),
                                              bucketCapacity*numContiguousBuckets, maxNumPerEntity);

    auto bucketNumPerEntity = stk::mesh::field_scalars_per_entity(*hostField, *bucket);
    fill_host_view(bucket, unHostInnerView, bucketNumPerEntity);

    T* bufferPtr = bufferData.data() + selectedBucketOffset * bucketCapacity * maxNumPerEntity;
    UnmanagedDevInnerView<T> unBufferInnerView(bufferPtr, bucketCapacity*numContiguousBuckets, maxNumPerEntity);
    Kokkos::deep_copy(asyncCopyState.execSpace, unBufferInnerView, unHostInnerView);

    T* devicePtr = deviceData.data() + selectedBucketOffset * bucketCapacity * maxNumPerEntity;
    UnmanagedDevInnerView<T> unDeviceInnerView(devicePtr, ORDER_INDICES(bucketCapacity*numContiguousBuckets, maxNumPerEntity));

    transpose_buffer_into_contiguous_device_data(asyncCopyState.execSpace, bucketCapacity*numContiguousBuckets, maxNumPerEntity, unBufferInnerView, unDeviceInnerView);
  }

  void copy_new_and_modified_buckets_from_host(const BucketVector& buckets, unsigned numPerEntity)
  {
    for(unsigned bucketIdx = 0; bucketIdx < buckets.size(); bucketIdx++) {
      Bucket* bucket = buckets[bucketIdx];
      unsigned oldBucketId = bucket->get_ngp_field_bucket_id(get_ordinal());

      if(oldBucketId == INVALID_BUCKET_ID || bucket->get_ngp_field_bucket_is_modified(get_ordinal())) {
        copy_bucket_from_host_to_device(bucket, numPerEntity);
      }
    }
  }

  void construct_view(const BucketVector& buckets, const std::string& name, unsigned numPerEntity)
  {
    unsigned numBuckets = buckets.size();
    FieldDataDeviceViewType<T> tempDataDeviceView = FieldDataDeviceViewType<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets,
                                                                               ORDER_INDICES(bucketCapacity, numPerEntity));
    FieldDataHostViewType<T> tempDataHostView = FieldDataHostViewType<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets,
                                                                         bucketCapacity, numPerEntity);
    fieldSyncDebugger.initialize_view(tempDataDeviceView);
    fieldSyncDebugger.initialize_view(tempDataHostView);

    copy_unmodified_buckets(buckets, tempDataDeviceView, tempDataHostView, numPerEntity);

    deviceData = tempDataDeviceView;
    hostData = tempDataHostView;

    bufferData = FieldDataDeviceViewType<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets, bucketCapacity, numPerEntity);
  }

  void construct_new_index_view(const BucketVector& allBuckets)
  {
    Selector selector(*hostField);
    unsigned bucketIndex = 0;

    newDeviceSelectedBucketOffset = UnsignedViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, hostField->name() + "_bucket_offset"),
                                                     allBuckets.size());
    newHostSelectedBucketOffset = Kokkos::create_mirror_view(Kokkos::HostSpace(), newDeviceSelectedBucketOffset, Kokkos::WithoutInitializing);

    for(unsigned i = 0; i < allBuckets.size(); i++) {
      if(selector(*allBuckets[i])) {
        newHostSelectedBucketOffset(i) = bucketIndex++;
      } else {
        newHostSelectedBucketOffset(i) = INVALID_BUCKET_ID;
      }
    }

    Kokkos::deep_copy(asyncCopyState.execSpace, newDeviceSelectedBucketOffset, newHostSelectedBucketOffset);
  }

  void construct_field_exist_view(const BucketVector& allBuckets, const Selector& selector)
  {
    deviceFieldExistsOnBucket = BoolViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, hostField->name() + "_exists_on_bucket"),
                                             allBuckets.size());
    hostFieldExistsOnBucket = Kokkos::create_mirror_view(Kokkos::HostSpace(), deviceFieldExistsOnBucket, Kokkos::WithoutInitializing);
    Kokkos::deep_copy(asyncCopyState.execSpace, hostFieldExistsOnBucket, false);
    for (size_t i = 0; i < allBuckets.size(); ++i) {
      if(selector(*allBuckets[i])) {
        hostFieldExistsOnBucket(allBuckets[i]->bucket_id()) = true;
      }
    }
    Kokkos::deep_copy(asyncCopyState.execSpace, deviceFieldExistsOnBucket, hostFieldExistsOnBucket);
  }

  void construct_unsigned_bucket_views(const BucketVector & buckets, const std::string& suffix,
                                       typename UnsignedViewType::HostMirror& hostView, UnsignedViewType& deviceView)
  {
    if (deviceView.extent(0) != buckets.size()) {
      deviceView = UnsignedViewType(Kokkos::ViewAllocateWithoutInitializing(hostField->name() + suffix), buckets.size());
      hostView = Kokkos::create_mirror_view(deviceView);
    }
  }

  void construct_all_fields_buckets_num_components_per_entity_view(const BucketVector & allBuckets) {
    construct_unsigned_bucket_views(allBuckets, "_numComponentsPerEntity", hostAllFieldsBucketsNumComponentsPerEntity, deviceAllFieldsBucketsNumComponentsPerEntity);

    for (stk::mesh::Bucket * bucket : allBuckets) {
      hostAllFieldsBucketsNumComponentsPerEntity[bucket->bucket_id()] = stk::mesh::field_scalars_per_entity(*hostField, *bucket);
    }

    Kokkos::deep_copy(asyncCopyState.execSpace, deviceAllFieldsBucketsNumComponentsPerEntity, hostAllFieldsBucketsNumComponentsPerEntity);
  }

  void construct_field_buckets_num_components_per_entity_view(const BucketVector & buckets) {
    construct_unsigned_bucket_views(buckets, "_numFieldBucketComponentsPerEntity", hostFieldBucketsNumComponentsPerEntity, deviceFieldBucketsNumComponentsPerEntity);

    for (size_t i=0; i<buckets.size(); ++i) {
      hostFieldBucketsNumComponentsPerEntity[i] = stk::mesh::field_scalars_per_entity(*hostField, *buckets[i]);
    }

    Kokkos::deep_copy(asyncCopyState.execSpace, deviceFieldBucketsNumComponentsPerEntity, hostFieldBucketsNumComponentsPerEntity);
  }

  void construct_bucket_sizes_view(const BucketVector & buckets) {
    construct_unsigned_bucket_views(buckets, "_bucketSizes", hostBucketSizes, deviceBucketSizes);

    for (size_t i=0; i<buckets.size(); ++i) {
      hostBucketSizes[i] = buckets[i]->size();
    }

    Kokkos::deep_copy(asyncCopyState.execSpace, deviceBucketSizes, hostBucketSizes);
  }

  void copy_contiguous_buckets_from_device_to_host(const BucketVector& buckets, unsigned numPerEntity, unsigned& i)
  {
    copy_bucket_from_device_to_host(buckets[i], numPerEntity, 1);
  }

  void copy_selected_buckets_to_host()
  {
    if(hostField) {
      Selector selector = selectField(*hostField) & *syncSelector;
      const BucketVector& buckets = hostBulk->get_buckets(rank, selector);

      for(unsigned i = 0; i < buckets.size(); i++) {
        copy_contiguous_buckets_from_device_to_host(buckets, maxNumScalarsPerEntity, i);
      }
    }
  }

  void sync_to_host_using_selector()
  {
    if (!userSpecifiedSelector) {
      auto& execSpace = asyncCopyState.execSpace;
      transpose_all_device_data_into_buffer(execSpace, *hostField, deviceData, bufferData, deviceBucketSizes, deviceFieldBucketsNumComponentsPerEntity);

      Kokkos::deep_copy(execSpace, hostData, bufferData);
    }
    else {
      copy_selected_buckets_to_host();
    }
    clear_sync_state_flags();
  }

  void copy_host_data_to_stk_field_data(Selector& selector)
  {
    const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
    hostField->increment_num_syncs_to_host();
    copy_data(buckets, [](T &hostFieldData, T &stkFieldData){
      stkFieldData = hostFieldData;
    }, selector);
  }

  void copy_device_to_host()
  {
    if (hostField) {
      stk::mesh::Selector selector = stk::mesh::selectField(*hostField);

      if(userSpecifiedSelector) {
        selector &= *syncSelector;
      }

      sync_to_host_using_selector();

      if(asyncCopyState.syncMode == impl::DEVICE_TO_HOST) {
        fence();
        copy_host_data_to_stk_field_data(selector);
      }
    }
  }

  void internal_sync_to_host()
  {
    if (need_sync_to_host()) {
      ProfilingBlock prof("copy_to_host for " + hostField->name());
      copy_device_to_host();
      fieldSyncDebugger.sync_to_host(this);
      if (asyncCopyState.syncMode == impl::DEVICE_TO_HOST) {
        reset_sync_selector();
      }
    }
  }

  void copy_contiguous_buckets_from_host_to_device(const BucketVector& buckets, unsigned numPerEntity, unsigned& i)
  {
    copy_bucket_from_host_to_device(buckets[i], numPerEntity, 1);
  }

  void copy_selected_buckets_to_device()
  {
    if(hostField) {
      Selector selector = selectField(*hostField) & *syncSelector;
      const BucketVector& buckets = hostBulk->get_buckets(rank, selector);

      for(unsigned i = 0; i < buckets.size(); i++) {
        copy_contiguous_buckets_from_host_to_device(buckets, maxNumScalarsPerEntity, i);
      }
    }
  }

  void sync_to_device_using_selector()
  {
    if (!userSpecifiedSelector) {
      auto& execSpace = asyncCopyState.execSpace;
      Kokkos::deep_copy(execSpace, bufferData, hostData);

      transpose_buffer_into_all_device_data(execSpace, *hostField, bufferData, deviceData, deviceBucketSizes, deviceFieldBucketsNumComponentsPerEntity);
    }
    else {
      copy_selected_buckets_to_device();
    }
    clear_sync_state_flags();
  }

  void copy_host_to_device()
  {
    if (hostField) {
      stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
      const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);

      if(userSpecifiedSelector) {
        selector &= *syncSelector;
      }

      hostField->increment_num_syncs_to_device();
      copy_data(buckets, [](T &hostFieldData, T &stkFieldData){
        hostFieldData = stkFieldData;
      }, selector);
      sync_to_device_using_selector();
    }
  }

  void internal_sync_to_device()
  {
    bool needToSyncToDevice = need_sync_to_device();

    if (needToSyncToDevice) {
      ProfilingBlock prof("copy_to_device for " + hostField->name());
      if (hostBulk->synchronized_count() == synchronizedCount) {
        copy_host_to_device();
        fieldSyncDebugger.sync_to_device(this);
      }
      else {
        update_field(needToSyncToDevice);
      }
      reset_sync_selector();
    }
  }

  template <typename ViewType>
  KOKKOS_FUNCTION
  void swap_views(ViewType & view1, ViewType & view2)
  {
    ViewType tmpView = view2;
    view2 = view1;
    view1 = tmpView;
  }

  template <typename ViewType>
  void allocate_view(ViewType &view, typename ViewType::HostMirror &host, const std::string &name, size_t size)
  {
    view = ViewType(name, size);
    host = Kokkos::create_mirror_view(view);
  }

  template <typename Assigner>
  void copy_data(const stk::mesh::BucketVector& buckets, const Assigner &assigner, const Selector& selector)
  {
    ProfilingBlock prof("copy_data for " + hostField->name());
    for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
    {
      const stk::mesh::Bucket &bucket = *buckets[iBucket];

      if(selector(bucket) == false) { continue; }

      const unsigned numPerEntity = stk::mesh::field_scalars_per_entity(*hostField, bucket);
      unsigned selectedBucketOffset = newHostSelectedBucketOffset(bucket.bucket_id());

      T* data = static_cast<T*>(stk::mesh::field_data<FieldBase, EmptyStkFieldSyncDebugger>(*hostField, bucket));
      for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
      {
        for(unsigned j=0; j<numPerEntity; j++)
        {
          assigner(hostData(selectedBucketOffset, iEntity, j), data[iEntity*numPerEntity + j]);
        }
      }
    }
  }

  friend NgpDebugger<T>;

  template <typename U, template <typename> class NgpDebuggerType>
  friend void impl::internal_fence_no_sync_to_host(NgpField<U, NgpDebuggerType>& ngpField);

  FieldDataHostViewType<T> hostData;
  FieldDataDeviceViewType<T> deviceData;
  FieldDataDeviceViewType<T> bufferData;
  Kokkos::View<bool, Kokkos::HostSpace> needSyncToHost;
  Kokkos::View<bool, Kokkos::HostSpace> needSyncToDevice;

  typename BoolViewType::HostMirror hostFieldExistsOnBucket;
  BoolViewType deviceFieldExistsOnBucket;

  typename UnsignedViewType::HostMirror hostSelectedBucketOffset;
  UnsignedViewType deviceSelectedBucketOffset;
  typename UnsignedViewType::HostMirror newHostSelectedBucketOffset;
  UnsignedViewType newDeviceSelectedBucketOffset;
  typename UnsignedViewType::HostMirror hostAllFieldsBucketsNumComponentsPerEntity;
  UnsignedViewType deviceAllFieldsBucketsNumComponentsPerEntity;

  stk::mesh::EntityRank rank;
  unsigned ordinal;

  const stk::mesh::BulkData* hostBulk;
  const stk::mesh::FieldBase* hostField;
  unsigned bucketCapacity;
  unsigned numBucketsForField;
  unsigned maxNumScalarsPerEntity;
  size_t synchronizedCount;

  DeviceField<T, NgpDebugger>* stateFields[MaximumFieldStates];

  Kokkos::View<int[1], Kokkos::HostSpace> copyCounter;
  bool userSpecifiedSelector;
  Selector* syncSelector;

  typename UnsignedViewType::HostMirror hostBucketSizes;
  UnsignedViewType deviceBucketSizes;

  typename UnsignedViewType::HostMirror hostFieldBucketsNumComponentsPerEntity;
  UnsignedViewType deviceFieldBucketsNumComponentsPerEntity;

  impl::AsyncCopyState asyncCopyState;

  NgpDebugger<T> fieldSyncDebugger;
};

}
}

#endif

