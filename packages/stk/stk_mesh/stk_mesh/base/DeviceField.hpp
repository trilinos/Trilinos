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

#include <cmath>

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
#include "stk_mesh/base/EntityFieldData.hpp"
#include "stk_mesh/baseImpl/NgpFieldAux.hpp"
#include "stk_mesh/base/NgpFieldSyncDebugger.hpp"

namespace stk {
namespace mesh {

constexpr unsigned NUM_COMPONENTS_INDEX = 0;
constexpr unsigned FIRST_DIMENSION_INDEX = 1;
constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T, template <typename> class NgpDebugger> class DeviceField;

namespace impl {
  constexpr double OVERALLOCATION_FACTOR = 1.1;

  inline long allocation_size(const double size_requested)
  {
    return std::lround(size_requested*OVERALLOCATION_FACTOR);
  }

  template <typename T> const FieldDataDeviceViewType<T> get_device_data(const DeviceField<T>& deviceField);
  template <typename T> FieldDataDeviceViewType<T> get_device_data(DeviceField<T>&);
}

template<typename T, template <typename> class NgpDebugger>
class DeviceField : public NgpFieldBase
{
private:
  using StkDebugger = typename NgpDebugger<T>::StkFieldSyncDebuggerType;

 public:
  using ExecSpace = stk::ngp::ExecSpace;
  using value_type = T;

  KOKKOS_FUNCTION
  DeviceField()
    : NgpFieldBase(),
      rank(stk::topology::INVALID_RANK),
      ordinal(INVALID_ORDINAL),
      hostBulk(nullptr),
      hostField(nullptr),
      bucketCapacity(0),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(0),
      fieldSyncDebugger(nullptr)
  {
  }

  DeviceField(const BulkData& bulk, const FieldBase &stkField, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      rank(stkField.entity_rank()),
      ordinal(stkField.mesh_meta_data_ordinal()),
      hostBulk(&bulk),
      hostField(&stkField),
      bucketCapacity(bulk.get_maximum_bucket_capacity()),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(0),
      fieldSyncDebugger(nullptr)
  {
    STK_ThrowRequireMsg(isFromGetUpdatedNgpField, "NgpField must be obtained from get_updated_ngp_field()");
    initialize();
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceField(const DeviceField<T, NgpDebugger>&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField(DeviceField<T, NgpDebugger>&&) = default;
  KOKKOS_FUNCTION ~DeviceField() {}
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T, NgpDebugger>& operator=(const DeviceField<T, NgpDebugger>&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T, NgpDebugger>& operator=(DeviceField<T, NgpDebugger>&&) = default;

  void initialize()
  {
    hostField->template make_field_sync_debugger<StkDebugger>();
    fieldSyncDebugger = NgpDebugger<T>(&hostField->get_field_sync_debugger<StkDebugger>());
  }

  void set_field_states(DeviceField<T, NgpDebugger>* /*fields*/[])
  {
  }

  void update_field(const ExecSpace& newExecSpace) override
  {
    set_execution_space(newExecSpace);
    update_field();
  }

  void update_field(ExecSpace&& newExecSpace) override
  {
    set_execution_space(std::forward<ExecSpace>(newExecSpace));
    update_field();
  }

  size_t num_syncs_to_host() const override { return hostField->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return hostField->num_syncs_to_device(); }

  void fence() override {
    get_execution_space().fence();
    reset_execution_space();
  }

  void modify_on_host() override
  {
    set_modify_on_host();
  }

  void modify_on_host(const Selector& selector) override
  {
    modify_on_host();
  }

  void modify_on_device() override
  {
    set_modify_on_device();
  }

  void modify_on_device(const Selector& selector) override
  {
    modify_on_device();
  }

  void notify_sync_debugger_clear_sync_state() override
  {
    fieldSyncDebugger.clear_sync_state(this);
  }

  void notify_sync_debugger_clear_host_sync_state() override
  {
    fieldSyncDebugger.clear_host_sync_state(this);
  }

  void notify_sync_debugger_clear_device_sync_state() override
  {
    fieldSyncDebugger.clear_device_sync_state(this);
  }

  void clear_sync_state() override
  {
    hostField->clear_sync_state();
  }

  void clear_host_sync_state() override
  {
    hostField->clear_host_sync_state();
  }

  void clear_device_sync_state() override
  {
    hostField->clear_device_sync_state();
  }

  void sync_to_host() override
  {
    reset_execution_space();
    if (internal_sync_to_host()) {
      Kokkos::fence();
      reset_execution_space();
    }
  }

  void sync_to_host(const ExecSpace& newExecSpace) override
  {
    set_execution_space(newExecSpace);
    internal_sync_to_host();
  }

  void sync_to_host(ExecSpace&& newExecSpace) override
  {
    set_execution_space(std::forward<ExecSpace>(newExecSpace));
    internal_sync_to_host();
  }

  void sync_to_device() override
  {
    reset_execution_space();
    if (internal_sync_to_device()) {
      Kokkos::fence();
      reset_execution_space();
    }
  }

  void sync_to_device(const ExecSpace& newExecSpace) override
  {
    set_execution_space(newExecSpace);
    internal_sync_to_device();
  }

  void sync_to_device(ExecSpace&& newExecSpace) override
  {
    set_execution_space(std::forward<ExecSpace>(newExecSpace));
    internal_sync_to_device();
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  KOKKOS_FUNCTION
  unsigned get_component_stride() const
  {
    unsigned stride = 1;
#ifdef STK_USE_DEVICE_MESH
    stride = bucketCapacity;
#endif
    return stride;
  }

  KOKKOS_FUNCTION
  unsigned get_num_components_per_entity(const FastMeshIndex& entityIndex) const {
    const unsigned bucketId = entityIndex.bucket_id;
    return deviceAllFieldsBucketsLayoutPerEntity(bucketId, NUM_COMPONENTS_INDEX);
  }

  KOKKOS_FUNCTION
  unsigned get_extent0_per_entity(const FastMeshIndex& entityIndex) const {
    const unsigned bucketId = entityIndex.bucket_id;
    return deviceAllFieldsBucketsLayoutPerEntity(bucketId, FIRST_DIMENSION_INDEX);
  }

  KOKKOS_FUNCTION
  unsigned get_extent1_per_entity(const FastMeshIndex& entityIndex) const {
    const unsigned bucketId = entityIndex.bucket_id;
    const unsigned numScalars = deviceAllFieldsBucketsLayoutPerEntity(bucketId, NUM_COMPONENTS_INDEX);
    const unsigned firstDimension = deviceAllFieldsBucketsLayoutPerEntity(bucketId, FIRST_DIMENSION_INDEX);
    return (numScalars != 0) ? numScalars / firstDimension : 0;
  }

  KOKKOS_FUNCTION
  unsigned get_extent_per_entity(const FastMeshIndex& entityIndex, unsigned dimension) const {
    const unsigned bucketId = entityIndex.bucket_id;
    if (dimension == 0) {
      return deviceAllFieldsBucketsLayoutPerEntity(bucketId, FIRST_DIMENSION_INDEX);
    }
    else if (dimension == 1) {
      const unsigned numScalars = deviceAllFieldsBucketsLayoutPerEntity(bucketId, NUM_COMPONENTS_INDEX);
      const unsigned firstDimension = deviceAllFieldsBucketsLayoutPerEntity(bucketId, FIRST_DIMENSION_INDEX);
      return (numScalars != 0) ? numScalars / firstDimension : 0;
    }
    else {
      const unsigned numScalars = deviceAllFieldsBucketsLayoutPerEntity(bucketId, NUM_COMPONENTS_INDEX);
      return (numScalars != 0) ? 1 : 0;
    }
  }

  unsigned debug_get_bucket_offset(unsigned bucketOrdinal) const override {
    return fieldSyncDebugger.get_bucket_offset(bucketOrdinal);
  }

  template <typename Mesh> KOKKOS_FUNCTION
  T& get(const Mesh& ngpMesh, Entity entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component, fileName, lineNumber);
  }

  KOKKOS_FUNCTION
  T& get(FastMeshIndex index, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, component));
  }

  KOKKOS_FUNCTION
  T& operator()(const FastMeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    fieldSyncDebugger.device_stale_access_check(this, index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, component));
  }

  KOKKOS_FUNCTION
  EntityFieldData<T> operator()(const FastMeshIndex& index,
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
    clear_sync_state();
    Kokkos::deep_copy(deviceData, value);
    modify_on_device();
    fieldSyncDebugger.set_any_potential_device_field_modification(true);
  }

  KOKKOS_FUNCTION
  EntityRank get_rank() const { return rank; }

  KOKKOS_FUNCTION
  unsigned get_ordinal() const { return ordinal; }

  const BulkData& get_bulk() const { return *hostBulk; }

  FieldState state() const { return hostField->state(); }

  const FieldBase* get_field_base() const { return hostField; }

  void update_bucket_pointer_view() override
  {
    Selector selector = selectField(*hostField);
    auto hostFieldEntityRank = hostField->entity_rank();
    const BucketVector& buckets = hostBulk->get_buckets(hostFieldEntityRank, selector);
    construct_field_buckets_pointer_view(buckets);
  }

  void swap_field_views(NgpFieldBase *other) override
  {
    DeviceField<T,NgpDebugger>* deviceFieldT = dynamic_cast<DeviceField<T,NgpDebugger>*>(other);
    STK_ThrowRequireMsg(deviceFieldT != nullptr, "DeviceField::swap_field_views called with class that can't dynamic_cast to DeviceField<T,NgpDebugger>");
    swap_views(deviceData, deviceFieldT->deviceData);
  }

  KOKKOS_FUNCTION
  void swap(DeviceField<T, NgpDebugger> &other)
  {
    swap_views(deviceData, other.deviceData);
  }

  bool need_sync_to_host() const override
  {
    return hostField->need_sync_to_host();
  }

  bool need_sync_to_device() const override
  {
    return hostField->need_sync_to_device();
  }

  void debug_initialize_debug_views() override { fieldSyncDebugger.initialize_debug_views(this); }

 protected:

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
 ExecSpace& get_execution_space() const { return hostField->get_execution_space(); }

 void set_execution_space(const ExecSpace& executionSpace) { hostField->set_execution_space(executionSpace); }

 void set_execution_space(ExecSpace&& executionSpace)
 {
   hostField->set_execution_space(std::forward<ExecSpace>(executionSpace));
 }

 void reset_execution_space() { hostField->reset_execution_space(); }

 void set_modify_on_host() { hostField->modify_on_host(); }

 void set_modify_on_device() { hostField->modify_on_device(); }

 void clear_sync_state_flags() { hostField->clear_sync_state(); }

 void update_field()
 {
   STK_ThrowRequireMsg(hostBulk->synchronized_count() >= synchronizedCount,
       "Invalid sync state detected for NgpField: " << hostField->name());
   if (hostBulk->synchronized_count() == synchronizedCount) {
     return;
   }

   ProfilingBlock prof("update_field for " + hostField->name());
   Selector selector = selectField(*hostField);
   auto hostFieldEntityRank = hostField->entity_rank();
   const BucketVector& buckets = hostBulk->get_buckets(hostFieldEntityRank, selector);
   const BucketVector& allBuckets = hostBulk->buckets(hostFieldEntityRank);
   numBucketsForField = buckets.size();
   maxNumScalarsPerEntity = hostField->max_size();

   construct_field_buckets_pointer_view(buckets);
   construct_all_fields_buckets_data_layout_per_entity_view(allBuckets);
   construct_field_buckets_num_components_per_entity_view(buckets);
   construct_bucket_sizes_view(buckets);
   construct_new_index_view(allBuckets);

   if (numBucketsForField != deviceData.extent(0)) {
     construct_view(buckets, "deviceData_" + hostField->name(), maxNumScalarsPerEntity);
   } else {
     move_unmodified_buckets(buckets, maxNumScalarsPerEntity);
   }

   copy_new_and_modified_buckets_from_host(buckets, maxNumScalarsPerEntity);

   fieldSyncDebugger.update_field(this);

   for (auto* bucket : allBuckets) {
     bucket->set_ngp_field_bucket_id(get_ordinal(), INVALID_BUCKET_ID);
   }
   for (auto* bucket : buckets) {
     bucket->set_ngp_field_bucket_id(get_ordinal(), bucket->bucket_id());
   }

   hostField->increment_num_syncs_to_device();
   synchronizedCount = hostBulk->synchronized_count();
   hostSelectedBucketOffset = newHostSelectedBucketOffset;
   deviceSelectedBucketOffset = newDeviceSelectedBucketOffset;
 }

  void construct_view(const BucketVector& buckets, const std::string& name, unsigned numPerEntity)
  {
    unsigned numBuckets = buckets.size();
    FieldDataDeviceViewType<T> tempDataDeviceView = FieldDataDeviceViewType<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets,
                                                    ORDER_INDICES(bucketCapacity, numPerEntity));
    fieldSyncDebugger.initialize_view(tempDataDeviceView);

    copy_unmodified_buckets(buckets, tempDataDeviceView, numPerEntity);

    deviceData = tempDataDeviceView;
  }

  void construct_new_index_view(const BucketVector& allBuckets)
  {
    Selector selector(*hostField);
    unsigned bucketIndex = 0;

    newDeviceSelectedBucketOffset = UnsignedViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, hostField->name() + "_bucket_offset"),
                                                     allBuckets.size());
    newHostSelectedBucketOffset =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, newDeviceSelectedBucketOffset);

    for(unsigned i = 0; i < allBuckets.size(); i++) {
      if(selector(*allBuckets[i])) {
        newHostSelectedBucketOffset(i) = bucketIndex++;
      } else {
        newHostSelectedBucketOffset(i) = INVALID_BUCKET_ID;
      }
    }

    Kokkos::deep_copy(get_execution_space(), newDeviceSelectedBucketOffset, newHostSelectedBucketOffset);
  }

  template<typename ViewType>
  void construct_bucket_views(const BucketVector & buckets, const std::string& suffix,
                                   typename ViewType::HostMirror& hostView, ViewType& deviceView)
  {
    if (buckets.size() > deviceView.extent(0)) {
      deviceView = ViewType(Kokkos::ViewAllocateWithoutInitializing(hostField->name() + suffix), impl::allocation_size(buckets.size()));
      if (hostView.extent(0) != deviceView.extent(0)) {
#if defined STK_USE_DEVICE_MESH && !defined(STK_ENABLE_GPU)
        hostView = Kokkos::create_mirror(deviceView);
#else
        hostView = Kokkos::create_mirror_view(deviceView);
#endif
      }
    }
  }

  template<typename ViewType>
  void construct_2d_bucket_view(const BucketVector & buckets, const std::string& suffix, unsigned secondExtent,
                                ViewType& deviceView)
  {
    if ((buckets.size() > deviceView.extent(0)) || (secondExtent != deviceView.extent(1))) {
      deviceView = ViewType(Kokkos::ViewAllocateWithoutInitializing(hostField->name() + suffix),
                            impl::allocation_size(buckets.size()), secondExtent);
    }
  }

  void construct_all_fields_buckets_data_layout_per_entity_view(const BucketVector & allBuckets)
  {
    construct_2d_bucket_view(allBuckets, "_layoutPerEntity", 2, deviceAllFieldsBucketsLayoutPerEntity);
    typename Unsigned2dViewType::HostMirror hostAllFieldsBucketsLayoutPerEntity =
        Kokkos::create_mirror_view(deviceAllFieldsBucketsLayoutPerEntity);

    for (Bucket * bucket : allBuckets) {
      hostAllFieldsBucketsLayoutPerEntity(bucket->bucket_id(), NUM_COMPONENTS_INDEX)  = field_scalars_per_entity(*hostField, *bucket);
      hostAllFieldsBucketsLayoutPerEntity(bucket->bucket_id(), FIRST_DIMENSION_INDEX) = field_extent0_per_entity(*hostField, *bucket);
    }

    Kokkos::deep_copy(get_execution_space(), deviceAllFieldsBucketsLayoutPerEntity, hostAllFieldsBucketsLayoutPerEntity);
  }

  void construct_field_buckets_num_components_per_entity_view(const BucketVector & buckets)
  {
    construct_bucket_views(buckets, "_numFieldBucketComponentsPerEntity", hostBucketScratchMemory, deviceFieldBucketsNumComponentsPerEntity);

    for (size_t i=0; i<buckets.size(); ++i) {
      hostBucketScratchMemory[i] = field_scalars_per_entity(*hostField, *buckets[i]);
    }

    Kokkos::deep_copy(get_execution_space(), deviceFieldBucketsNumComponentsPerEntity, hostBucketScratchMemory);
  }

  void construct_field_buckets_pointer_view(const BucketVector& buckets)
  {
    if (buckets.size() > hostBucketPtrData.extent(0)) {
      hostBucketPtrData = FieldDataPointerHostViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "HostBucketPtrDataView"), impl::allocation_size(buckets.size()));
      deviceBucketPtrData = FieldDataPointerDeviceViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketPtrDataView"), impl::allocation_size(buckets.size()));
    }

    set_field_buckets_pointer_view(buckets);
  }

  void construct_bucket_sizes_view(const BucketVector & buckets)
  {
    construct_bucket_views(buckets, "_bucketSizes", hostBucketScratchMemory, deviceBucketSizes);

    for (size_t i=0; i<buckets.size(); ++i) {
      hostBucketScratchMemory[i] = buckets[i]->size();
    }

    Kokkos::deep_copy(get_execution_space(), deviceBucketSizes, hostBucketScratchMemory);
  }

  void set_field_buckets_pointer_view(const BucketVector& buckets)
  {
    STK_ThrowRequireMsg(hostBucketPtrData.extent(0) >= buckets.size(), "hostBucketPtrData is not large enough for the selected buckets");
    for(unsigned i = 0; i < buckets.size(); i++) {
      T* hostBucketPtr = reinterpret_cast<T*>(field_data<FieldBase, EmptyStkFieldSyncDebugger>(*hostField, *buckets[i]));
      T* deviceBucketPtr = hostBucketPtr;

#ifdef KOKKOS_ENABLE_CUDA
      cudaError_t status = cudaHostGetDevicePointer((void**)&deviceBucketPtr, (void*)hostBucketPtr, 0);

      STK_ThrowRequireMsg(status == cudaSuccess, "Something went wrong during cudaHostGetDevicePointer: " + std::string(cudaGetErrorString(status)));
#elif defined(KOKKOS_ENABLE_HIP) 
      hipError_t status = hipHostGetDevicePointer((void**)&deviceBucketPtr, (void*)hostBucketPtr, 0);

      STK_ThrowRequireMsg(status == hipSuccess, "Something went wrong during hipHostGetDevicePointer: " + std::string(hipGetErrorString(status)));
#endif

      hostBucketPtrData(i) = reinterpret_cast<uintptr_t>(deviceBucketPtr);
    }

    Kokkos::deep_copy(get_execution_space(), deviceBucketPtrData, hostBucketPtrData);
  }

  void copy_unmodified_buckets(const BucketVector& buckets, FieldDataDeviceViewType<T> destDevView, unsigned numPerEntity)
  {
    for(unsigned i = 0; i < buckets.size(); i++) {
      unsigned oldBucketId = buckets[i]->get_ngp_field_bucket_id(get_ordinal());
      unsigned newBucketId = buckets[i]->bucket_id();

      if(!buckets[i]->get_ngp_field_bucket_is_modified(get_ordinal())) {
        STK_ThrowRequire(deviceData.extent(0) != 0 && deviceSelectedBucketOffset.extent(0) != 0);
        copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(destDevView, deviceData, oldBucketId, newBucketId, numPerEntity);
      }
    }
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
    Kokkos::deep_copy(get_execution_space(), unInnerDestView, unInnerSrcView);
  }

  struct BackwardShiftIndices {
    BackwardShiftIndices(unsigned _oldIndex, unsigned _newIndex)
      : oldIndex(_oldIndex),
        newIndex(_newIndex)
    {}
    unsigned oldIndex;
    unsigned newIndex;
  };

  void move_unmodified_buckets(const BucketVector& buckets, unsigned numPerEntity)
  {
    std::vector<BackwardShiftIndices> backwardShiftList;

    for (unsigned i = 0; i < buckets.size(); ++i) {
      unsigned oldBucketId = buckets[i]->get_ngp_field_bucket_id(get_ordinal());
      unsigned newBucketId = buckets[i]->bucket_id();

      const bool isNewBucket = (oldBucketId == INVALID_BUCKET_ID);
      if (isNewBucket) {
        continue;
      }

      unsigned oldBucketOffset = hostSelectedBucketOffset(oldBucketId);
      unsigned newBucketOffset = newHostSelectedBucketOffset(newBucketId);

      const bool bucketNotForThisField = (newBucketOffset == INVALID_BUCKET_ID);
      const bool bucketHasNotMoved = (oldBucketId == newBucketId);
      const bool bucketIsUnmodified = not buckets[i]->get_ngp_field_bucket_is_modified(get_ordinal());

      if (bucketNotForThisField || (bucketHasNotMoved && bucketIsUnmodified)) {
        continue;
      }

      if (newBucketOffset < oldBucketOffset) {
        shift_bucket_forward(oldBucketId, newBucketId, numPerEntity);
      }
      else {
        backwardShiftList.emplace_back(oldBucketId, newBucketId);
      }
    }

    shift_buckets_backward(backwardShiftList, numPerEntity);
  }

  void shift_bucket_forward(unsigned oldBucketId, unsigned newBucketId, unsigned numPerEntity)
  {
    copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(deviceData, deviceData,
                                                                                        oldBucketId, newBucketId,
                                                                                        numPerEntity);
  }

  void shift_buckets_backward(const std::vector<BackwardShiftIndices> & backwardShiftList, unsigned numPerEntity)
  {
    for (auto it = backwardShiftList.rbegin(); it != backwardShiftList.rend(); ++it) {
      const BackwardShiftIndices& shiftIndices = *it;
      copy_moved_device_bucket_data<FieldDataDeviceViewType<T>, UnmanagedDevInnerView<T>>(deviceData, deviceData,
                                                                                          shiftIndices.oldIndex,
                                                                                          shiftIndices.newIndex,
                                                                                          numPerEntity);
    }
  }

  void copy_new_and_modified_buckets_from_host(const BucketVector& buckets, unsigned numPerEntity)
  {
    construct_bucket_views(buckets, "_bucketSizes", hostBucketScratchMemory, deviceFieldBucketsMarkedModified);

    for(unsigned bucketIdx = 0; bucketIdx < buckets.size(); bucketIdx++) {
      Bucket* bucket = buckets[bucketIdx];
      unsigned oldBucketId = bucket->get_ngp_field_bucket_id(get_ordinal());
      bool isModified = oldBucketId == INVALID_BUCKET_ID || bucket->get_ngp_field_bucket_is_modified(get_ordinal());
      
      hostBucketScratchMemory(bucketIdx) = isModified ? 1 : 0;
    }

    Kokkos::deep_copy(get_execution_space(), deviceFieldBucketsMarkedModified, hostBucketScratchMemory);

    impl::transpose_new_and_modified_buckets_to_device(get_execution_space(), deviceBucketPtrData, deviceData,
                                                       deviceBucketSizes, deviceFieldBucketsNumComponentsPerEntity, deviceFieldBucketsMarkedModified);
  }

  void copy_device_to_host()
  {
    if (hostField) {
      impl::transpose_to_pinned_and_mapped_memory(get_execution_space(), deviceBucketPtrData, deviceData, deviceBucketSizes, deviceFieldBucketsNumComponentsPerEntity);
      clear_device_sync_state();
      hostField->increment_num_syncs_to_host();
    }
  }

  bool internal_sync_to_host()
  {
    if (need_sync_to_host()) {
      ProfilingBlock prof("copy_to_host for " + hostField->name());
      copy_device_to_host();
      fieldSyncDebugger.sync_to_host(this);
      return true;
    }

    return false;
  }

  void copy_host_to_device()
  {
    if (hostField) {
      impl::transpose_from_pinned_and_mapped_memory(get_execution_space(), deviceBucketPtrData, deviceData, deviceBucketSizes, deviceFieldBucketsNumComponentsPerEntity);
      clear_host_sync_state();
      hostField->increment_num_syncs_to_device();
    }
  }

  bool internal_sync_to_device()
  {
    if (need_sync_to_device()) {
      ProfilingBlock prof("copy_to_device for " + hostField->name());
      if(hostBulk->synchronized_count() != synchronizedCount) {
        update_field();
      }
      copy_host_to_device();
      fieldSyncDebugger.sync_to_device(this);

      return true;
    }

    return false;
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

  friend NgpDebugger<T>;
  friend const FieldDataDeviceViewType<T> impl::get_device_data<T>(const DeviceField<T>&);
  friend FieldDataDeviceViewType<T> impl::get_device_data<T>(DeviceField<T>&);

  FieldDataDeviceViewType<T> deviceData;

  FieldDataPointerHostViewType hostBucketPtrData;
  FieldDataPointerDeviceViewType deviceBucketPtrData;

  typename UnsignedViewType::HostMirror hostSelectedBucketOffset;
  UnsignedViewType deviceSelectedBucketOffset;
  typename UnsignedViewType::HostMirror newHostSelectedBucketOffset;
  UnsignedViewType newDeviceSelectedBucketOffset;
  Unsigned2dViewType deviceAllFieldsBucketsLayoutPerEntity;

  EntityRank rank;
  unsigned ordinal;

  const BulkData* hostBulk;
  const FieldBase* hostField;
  unsigned bucketCapacity;
  unsigned numBucketsForField;
  unsigned maxNumScalarsPerEntity;
  size_t synchronizedCount;

  typename UnsignedViewType::HostMirror hostBucketScratchMemory;
  UnsignedViewType deviceBucketSizes;
  UnsignedViewType deviceFieldBucketsNumComponentsPerEntity;
  UnsignedViewType deviceFieldBucketsMarkedModified;

  NgpDebugger<T> fieldSyncDebugger;
};

namespace impl {

//not for public consumption. calling this will void your warranty.
template<typename T>
const FieldDataDeviceViewType<T> get_device_data(const DeviceField<T>& deviceField)
{
  return deviceField.deviceData;
}

template<typename T>
FieldDataDeviceViewType<T> get_device_data(DeviceField<T>& deviceField)
{
  return deviceField.deviceData;
}

} //namespace impl

}
}

#endif

