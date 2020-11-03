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

#ifndef STK_MESH_NGPFIELD_HPP
#define STK_MESH_NGPFIELD_HPP

#include "stk_util/stk_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"
#include "stk_mesh/base/FieldSyncDebugging.hpp"
#include "stk_mesh/base/NgpUtils.hpp"

namespace stk {
namespace mesh {

constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T>
class EntityFieldData {
public:
  KOKKOS_FUNCTION EntityFieldData(T* dataPtr, unsigned length, unsigned componentStride=1)
  : fieldDataPtr(dataPtr), fieldDataLength(length), fieldComponentStride(componentStride)
  {}
  KOKKOS_DEFAULTED_FUNCTION ~EntityFieldData() = default;

  KOKKOS_FUNCTION unsigned size() const { return fieldDataLength; }
  KOKKOS_FUNCTION T& operator[](unsigned idx) { return fieldDataPtr[idx*fieldComponentStride]; }
  KOKKOS_FUNCTION const T& operator[](unsigned idx) const { return fieldDataPtr[idx*fieldComponentStride]; }

private:
  T* fieldDataPtr;
  unsigned fieldDataLength;
  unsigned fieldComponentStride;
};

template<typename T> using UnmanagedHostInnerView = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template<typename T> using UnmanagedDevInnerView = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template<typename T>
class HostField : public NgpFieldBase
{
public:
  typedef T value_type;

  HostField()
    : NgpFieldBase(),
      field(nullptr)
  {
    needSyncToHost = std::make_shared<bool>(false);
    needSyncToDevice = std::make_shared<bool>(false);
  }

  HostField(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      field(&f)
  {
    needSyncToHost = std::make_shared<bool>(false);
    needSyncToDevice = std::make_shared<bool>(false);
  }

  HostField(const HostField<T>&) = default;
  HostField(HostField<T>&&) = default;
  HostField<T>& operator=(const HostField<T>&) = default;
  HostField<T>& operator=(HostField<T>&&) = default;

  virtual ~HostField() = default;

  void update_field(bool needToSyncAllDataToDevice = false) override {}

  void set_field_states(HostField<T>* fields[]) {}

  size_t num_syncs_to_host() const override { return field->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return field->num_syncs_to_device(); }

  unsigned get_component_stride() const { return 1; }

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
  }

#ifdef STK_DEBUG_FIELD_SYNC
  unsigned get_bucket_offset(unsigned bucketOrdinal) const override {
    return bucketOrdinal;
  }

  void detect_device_field_modification() override {
  }

  void update_debug_storage(size_t hostSynchronizedCount) override {
  }

  bool lost_device_field_data() const override {
    return false;
  }

  T& get(const HostMesh& ngpMesh, stk::mesh::Entity entity, int component,
         const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity));
    ThrowAssert(data);
    return data[component];
  }

  T& get(stk::mesh::FastMeshIndex entity, int component,
         const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket_id, entity.bucket_ord));
    ThrowAssert(data);
    return data[component];
  }

  T& get(HostMesh::MeshIndex entity, int component,
         const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T* data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket->bucket_id(), entity.bucketOrd));
    ThrowAssert(data);
    return data[component];
  }

  T& operator()(const stk::mesh::FastMeshIndex& index, int component,
                const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    ThrowAssert(data);
    return data[component];
  }

  T& operator()(const HostMesh::MeshIndex& index, int component,
                const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T* data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket->bucket_id(), index.bucketOrd));
    ThrowAssert(data);
    return data[component];
  }

  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index,
                                const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    unsigned numScalars = stk::mesh::field_scalars_per_entity(*field, index.bucket_id);
    ThrowAssert(data);
    return EntityFieldData<T>(data, numScalars);
  }

#else

  T& get(const HostMesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity));
    ThrowAssert(data);
    return data[component];
  }

  T& get(stk::mesh::FastMeshIndex entity, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket_id, entity.bucket_ord));
    ThrowAssert(data);
    return data[component];
  }

  T& get(HostMesh::MeshIndex entity, int component) const
  {
    T* data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket->bucket_id(), entity.bucketOrd));
    ThrowAssert(data);
    return data[component];
  }

  T& operator()(const stk::mesh::FastMeshIndex& index, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    ThrowAssert(data);
    return data[component];
  }

  T& operator()(const HostMesh::MeshIndex& index, int component) const
  {
    T* data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket->bucket_id(), index.bucketOrd));
    ThrowAssert(data);
    return data[component];
  }

  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    unsigned numScalars = stk::mesh::field_scalars_per_entity(*field, index.bucket_id);
    ThrowAssert(data);
    return EntityFieldData<T>(data, numScalars);
  }
#endif

  void set_all(const HostMesh& ngpMesh, const T& value)
  {
    stk::mesh::for_each_entity_run(ngpMesh, field->entity_rank(), *field, KOKKOS_LAMBDA(const FastMeshIndex& entity) {
                                     T* fieldPtr = static_cast<T*>(stk::mesh::field_data(*field, entity.bucket_id, entity.bucket_ord));
                                     int numScalars = stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
                                     for (int i=0; i<numScalars; i++) {
                                       fieldPtr[i] = value;
                                     }
                             });
  }

  void modify_on_host() override { set_modify_on_host(); }

  void modify_on_device() override { set_modify_on_device(); }

  void clear_sync_state() override
  {
    *needSyncToHost = false;
    *needSyncToDevice = false;
  }

  void modify_on_host(const Selector& selector) override { set_modify_on_host(); }

  void modify_on_device(const Selector& selector) override { set_modify_on_device(); }

  void clear_host_sync_state() override
  {
    *needSyncToDevice = false;
  }

  void clear_device_sync_state() override
  {
    *needSyncToHost = false;
  }

  void sync_to_host() override
  {
    if (need_sync_to_host()) {
      field->increment_num_syncs_to_host();
      clear_sync_state();
    }
  }

  void sync_to_device() override 
  {
    if (need_sync_to_device()) {
      field->increment_num_syncs_to_device();
      clear_sync_state();
    }
  }

  size_t synchronized_count() const override { return field->mesh_meta_data().mesh_bulk_data().synchronized_count(); }

  FieldState state() const { return field->state(); }

  void rotate_multistate_data() override { }

  void swap(HostField<T> &other) {
    needSyncToHost.swap(other.needSyncToHost);
    needSyncToDevice.swap(other.needSyncToDevice);
  }

  stk::mesh::EntityRank get_rank() const { return field->entity_rank(); }

  unsigned get_ordinal() const { return field->mesh_meta_data_ordinal(); }


protected:
  bool need_sync_to_host() const { return *needSyncToHost; }

  bool need_sync_to_device() const { return *needSyncToDevice; }

  unsigned get_contiguous_bucket_offset_end(const BucketVector& buckets, unsigned i) { return i; }

private:
  void set_modify_on_host() {
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    ThrowRequire(*needSyncToHost == false);
#endif
    *needSyncToDevice = true;
  }

  void set_modify_on_device() {
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    ThrowRequire(*needSyncToDevice == false);
#endif
    *needSyncToHost = true;
  }

  void copy_host_to_device() { }

  void copy_device_to_host() { }

  const stk::mesh::FieldBase * field;

  std::shared_ptr<bool> needSyncToHost;
  std::shared_ptr<bool> needSyncToDevice;
};

template<typename T>
class DeviceField : public NgpFieldBase
{
private:
  typedef Kokkos::View<T***, Kokkos::LayoutRight, MemSpace> FieldDataDeviceViewType;
  typedef Kokkos::View<T***, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> FieldDataDeviceUnmanagedViewType;
  typedef Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::HostSpace> FieldDataHostViewType;
#ifdef STK_DEBUG_FIELD_SYNC
  template <typename U> using HostArrayType = Kokkos::View<U*, Kokkos::HostSpace>;
#endif

public:
  typedef T value_type;

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
      syncSelector(nullptr)
  {
    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
  }

  DeviceField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      rank(field.entity_rank()),
      ordinal(field.mesh_meta_data_ordinal()),
      hostBulk(&bulk),
      hostField(&field),
      bucketCapacity(0),
      numBucketsForField(0),
      maxNumScalarsPerEntity(0),
      synchronizedCount(bulk.synchronized_count()),
      copyCounter("copy_counter"),
      userSpecifiedSelector(false),
      syncSelector(new Selector())
  {
    needSyncToHost = Kokkos::View<bool, HostExecSpace>("needSyncToHost");
    needSyncToDevice = Kokkos::View<bool, HostExecSpace>("needSyncToDevice");

    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
    ThrowRequireMsg(isFromGetUpdatedNgpField, "NgpField must be obtained from get_updated_ngp_field()");

    update_field();

#ifdef STK_DEBUG_FIELD_SYNC
    initialize_debug_views();
#endif
  }

  void set_field_states(DeviceField<T>* fields[])
  {
    const unsigned numStates = hostField->number_of_states();
    for (unsigned state = 0; state < numStates; ++state) {
      stateFields[state] = fields[state];
    }
  }

  void update_field(bool needToSyncAllDataToDevice = false) override
  {
    ProfilingBlock prof("update_field for " + hostField->name());
    stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
    const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
    const stk::mesh::BucketVector& allBuckets = hostBulk->buckets(hostField->entity_rank());
    numBucketsForField = buckets.size();
    maxNumScalarsPerEntity = hostField->max_size(rank);

    if (!buckets.empty()) {
      bucketCapacity = buckets[0]->capacity();
    }

    construct_field_exist_view(allBuckets, selector);
    construct_num_components_per_entity_view(allBuckets);
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

#ifdef STK_DEBUG_FIELD_SYNC
    if (lastFieldValue.extent(0) != 0) {
      lostDeviceFieldData() = false;
      anyPotentialDeviceFieldModification() = false;

      Kokkos::deep_copy(lastFieldValue, deviceData);

      if (needToSyncAllDataToDevice) {
        reset_last_modification_state(LastModLocation::HOST_OR_DEVICE);
      }

      hostSynchronizedCount() = hostBulk->synchronized_count();
    }

    isFieldLayoutConsistent = true;
#endif

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

  void modify_on_host() override
  {
    set_modify_on_host();
    userSpecifiedSelector = false;
    *syncSelector = Selector(hostBulk->mesh_meta_data().universal_part());
  }

  void modify_on_host(const Selector& selector) override
  {
#ifdef STK_DEBUG_FIELD_SYNC
    ThrowErrorMsg("Field Sync Debugging is not supported with Selector-based modify_on_host().");
#endif
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
#ifdef STK_DEBUG_FIELD_SYNC
    ThrowErrorMsg("Field Sync Debugging is not supported with Selector-based modify_on_device().");
#endif
    set_modify_on_device();
    userSpecifiedSelector = true;
    *syncSelector |= selector;
  }

  void clear_sync_state() override
  {
#ifdef STK_DEBUG_FIELD_SYNC
    if (hostBulk->synchronized_count() != synchronizedCount) {
      update_field(need_sync_to_device());
    }
    Kokkos::deep_copy(lastFieldValue, deviceData);
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;
    reset_last_modification_state(LastModLocation::HOST_OR_DEVICE);
#endif
    clear_sync_state_flags();
    reset_sync_selector();
  }

  void clear_host_sync_state() override
  {
#ifdef STK_DEBUG_FIELD_SYNC
    set_last_modification_state_bit(LastModLocation::DEVICE);
#endif
    if (need_sync_to_device()) {
      clear_sync_state_flags();
      reset_sync_selector();
    }
  }

  void clear_device_sync_state() override
  {
#ifdef STK_DEBUG_FIELD_SYNC
    if (hostBulk->synchronized_count() != synchronizedCount) {
      update_field(need_sync_to_device());
    }
    Kokkos::deep_copy(lastFieldValue, deviceData);
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;
    set_last_modification_state_bit(LastModLocation::HOST);
#endif
    if (need_sync_to_host()) {
      clear_sync_state_flags();
      reset_sync_selector();
    }
  }

  void sync_to_host() override
  {
    if (need_sync_to_host()) {
      ProfilingBlock prof("copy_to_host for " + hostField->name());
      copy_device_to_host();
#ifdef STK_DEBUG_FIELD_SYNC
      Kokkos::deep_copy(lastFieldValue, deviceData);
      reset_last_modification_state(LastModLocation::HOST_OR_DEVICE);
      lostDeviceFieldData() = false;
      anyPotentialDeviceFieldModification() = false;
#endif
      reset_sync_selector();
    }
  }

  void sync_to_device() override
  {
    bool needToSyncToDevice = need_sync_to_device();
    if (needToSyncToDevice) {
      ProfilingBlock prof("copy_to_device for " + hostField->name());
      if (hostBulk->synchronized_count() == synchronizedCount) {
        copy_host_to_device();
#ifdef STK_DEBUG_FIELD_SYNC
        Kokkos::deep_copy(lastFieldValue, deviceData);
        reset_last_modification_state(LastModLocation::HOST_OR_DEVICE);
        lostDeviceFieldData() = false;
        anyPotentialDeviceFieldModification() = false;
#endif
      }
      else {
        update_field(needToSyncToDevice);
      }
      reset_sync_selector();
    }
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  KOKKOS_DEFAULTED_FUNCTION DeviceField(const DeviceField &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T>& operator=(const DeviceField<T>&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField<T>& operator=(DeviceField<T>&&) = default;

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
    return deviceNumComponentsPerEntity[bucketId];
  }

#ifdef STK_DEBUG_FIELD_SYNC
  unsigned get_bucket_offset(unsigned bucketOrdinal) const override {
    return debugHostSelectedBucketOffset(bucketOrdinal);
  }

  void initialize_debug_views() {
    stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
    const stk::mesh::BucketVector& allBuckets = hostBulk->buckets(hostField->entity_rank());
    const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
    unsigned numPerEntity = hostField->max_size(rank);

    construct_debug_views(allBuckets, buckets, numPerEntity);

    std::strcpy(fieldName.data(), hostField->name().data());
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;

    debugHostSelectedBucketOffset = newHostSelectedBucketOffset;
    debugDeviceSelectedBucketOffset = newDeviceSelectedBucketOffset;

    Kokkos::deep_copy(lastFieldValue, deviceData);

    reset_last_modification_state(LastModLocation::HOST_OR_DEVICE);

    hostSynchronizedCount() = hostBulk->synchronized_count();
    isFieldLayoutConsistent = true;
  }

  void construct_debug_views(const BucketVector & allBuckets, const BucketVector & buckets, unsigned numPerEntity) {
    if (buckets.size() != allBuckets.size()) {
      ThrowErrorMsg("Field Sync Debugging is not supported with partial Field allocations.");
    }
    if (buckets.size() != 0) {
      FieldDataDeviceViewType tempLastFieldValue = FieldDataDeviceViewType(hostField->name()+"_lastValue", buckets.size(),
                                                                           ORDER_INDICES(bucketCapacity, numPerEntity));
      LastFieldModLocationType tempLastFieldModLocation = LastFieldModLocationType(hostField->name()+"_lastModLocation", buckets.size(),
                                                                                   ORDER_INDICES(bucketCapacity, numPerEntity));
      HostArrayType<unsigned> tempLastMeshModBucketId = HostArrayType<unsigned>("lastMeshModBucketId", allBuckets.size());

      lastFieldValue = tempLastFieldValue;
      lastFieldModLocation = tempLastFieldModLocation;
      lastMeshModBucketId = tempLastMeshModBucketId;
    }

    hostField->set_last_modification_view(lastFieldModLocation);

    if (fieldName.extent(0) == 0u) {
      fieldName = FieldNameType(hostField->name()+"_name", hostField->name().size()+1);
      hostSynchronizedCount = ScalarUvmType<size_t>(hostField->name()+"_hostSyncCount");
      lostDeviceFieldData = ScalarUvmType<bool>(hostField->name()+"_anyDeviceFieldMod");
      anyPotentialDeviceFieldModification = ScalarUvmType<bool>(hostField->name()+"_anyPotentialDeviceFieldMod");
    }
  }

  bool lost_device_field_data() const override {
    return lostDeviceFieldData();
  }

  void detect_device_field_modification() {
    if (!anyPotentialDeviceFieldModification()) {
      return;
    }

    if (field_not_updated_after_mesh_mod()) {
      return;
    }

    stk::mesh::NgpMesh & ngpMesh = hostBulk->get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = hostBulk->mesh_meta_data();

    UnsignedViewType & localDeviceNumComponentsPerEntity = deviceNumComponentsPerEntity;
    FieldDataDeviceViewType & localDeviceData = deviceData;
    FieldDataDeviceViewType & localLastFieldValue = lastFieldValue;
    LastFieldModLocationType & localLastFieldModLocation = lastFieldModLocation;
    ScalarUvmType<bool> & localLostDeviceFieldData = lostDeviceFieldData;
    const bool inModCycle = hostBulk->in_modifiable_state();

    stk::mesh::for_each_entity_run(ngpMesh, rank, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& index)
    {
      const unsigned bucketId = index.bucket_id;
      const unsigned numComponents = localDeviceNumComponentsPerEntity(bucketId);
      for (unsigned component = 0; component < numComponents; ++component) {
        if (localDeviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) !=
            localLastFieldValue(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)))
        {
          localLastFieldModLocation(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) = LastModLocation::DEVICE;
          if (inModCycle) {
            localLostDeviceFieldData() = true;
          }
        }
      }
    });

    anyPotentialDeviceFieldModification() = false;
    Kokkos::fence();
  }

  void update_debug_storage(size_t hostSyncCount) {
    hostSynchronizedCount() = hostSyncCount;

    stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
    const stk::mesh::BucketVector & allBuckets = hostBulk->buckets(hostField->entity_rank());
    const stk::mesh::BucketVector & buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
    const unsigned numPerEntity = hostField->max_size(rank);

    construct_new_index_view(allBuckets);
    construct_debug_views(allBuckets, buckets, numPerEntity);

    debugHostSelectedBucketOffset = newHostSelectedBucketOffset;
    debugDeviceSelectedBucketOffset = newDeviceSelectedBucketOffset;
    isFieldLayoutConsistent = false;
  }

  void set_last_modification_state_bit(LastModLocation value) {
    for (size_t i = 0; i < lastFieldModLocation.extent(0); ++i) {
      for (size_t j = 0; j < lastFieldModLocation.extent(1); ++j) {
        for (size_t k = 0; k < lastFieldModLocation.extent(2); ++k) {
          lastFieldModLocation(i, j, k) = static_cast<LastModLocation>(lastFieldModLocation(i, j, k) | value);
        }
      }
    }
  }

  void reset_last_modification_state(LastModLocation value) {
    Kokkos::deep_copy(lastFieldModLocation, value);
  }

  STK_INLINE_FUNCTION
  bool data_is_stale_on_device(const stk::mesh::FastMeshIndex & index, int component) const {
    return !(lastFieldModLocation(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) &
             LastModLocation::DEVICE);
  }

  STK_INLINE_FUNCTION
  bool data_is_stale_on_device(const stk::mesh::DeviceMesh::MeshIndex & index, int component) const {
    return !(lastFieldModLocation(index.bucket->bucket_id(), ORDER_INDICES(index.bucketOrd, component)) &
             LastModLocation::DEVICE);
  }

  STK_INLINE_FUNCTION
  bool data_is_stale_on_device(int bucketId, int bucketOrdinal, int component) const {
    return !(lastFieldModLocation(bucketId, ORDER_INDICES(bucketOrdinal, component)) &
             LastModLocation::DEVICE);
  }

  STK_INLINE_FUNCTION
  void set_last_modification(const stk::mesh::FastMeshIndex & index, int component, LastModLocation location) const {
    lastFieldModLocation(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) = location;
  }

  STK_INLINE_FUNCTION
  bool field_not_updated_after_mesh_mod() const {
    return hostSynchronizedCount() != synchronizedCount;
  }

  STK_INLINE_FUNCTION
  void print_unupdated_field_warning(const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing un-updated Field %s on Device after mesh modification\n", fieldName.data());
    }
    else {
      printf("%s:%i *** WARNING: Accessing un-updated Field %s on Device after mesh modification\n",
             fileName, lineNumber, fieldName.data());
    }
  }

  STK_INLINE_FUNCTION
  void print_stale_data_warning(unsigned bucketId, unsigned bucketOrdinal, int component,
                                const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing stale data on Device for Field %s[%i]=%f\n",
             fieldName.data(), component,
             static_cast<double>(deviceData(bucketId, ORDER_INDICES(bucketOrdinal, component))));
    }
    else {
      printf("%s:%i *** WARNING: Accessing stale data on Device for Field %s[%i]=%f\n",
             fileName, lineNumber, fieldName.data(), component,
             static_cast<double>(deviceData(bucketId, ORDER_INDICES(bucketOrdinal, component))));
    }
  }

  STK_INLINE_FUNCTION
  void print_stale_data_warning_without_field_values(unsigned bucketId, unsigned bucketOrdinal, int component,
                                                     const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing stale data on Device for Field %s[%i]\n", fieldName.data(), component);
    }
    else {
      printf("%s:%i *** WARNING: Accessing stale data on Device for Field %s[%i]\n",
             fileName, lineNumber, fieldName.data(), component);
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const stk::mesh::FastMeshIndex & index, int component, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    if (data_is_stale_on_device(index, component)) {
      print_stale_data_warning(index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const stk::mesh::DeviceMesh::MeshIndex & index, int component, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    if (data_is_stale_on_device(index, component)) {
      print_stale_data_warning(index.bucket->bucket_id(), index.bucketOrd, component, fileName, lineNumber);
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const stk::mesh::FastMeshIndex & index, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    for (unsigned component = 0; component < get_num_components_per_entity(index); ++component) {
      if (data_is_stale_on_device(index, component)) {
        print_stale_data_warning(index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
      }
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    for (size_t i = 0; i < lastFieldModLocation.extent(0); ++i) {
      for (size_t j = 0; j < lastFieldModLocation.extent(1); ++j) {
        for (size_t k = 0; k < lastFieldModLocation.extent(2); ++k) {
          if (data_is_stale_on_device(i, ORDER_INDICES(j, k))) {
            print_stale_data_warning_without_field_values(i, ORDER_INDICES(j, k), fileName, lineNumber);
          }
        }
      }
    }
  }

  template <typename Mesh> KOKKOS_FUNCTION
  T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component, fileName, lineNumber);
  }

  KOKKOS_FUNCTION
  T& get(stk::mesh::FastMeshIndex entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(entity, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(entity.bucket_id), ORDER_INDICES(entity.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& get(MeshIndex entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(entity, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(entity.bucket->bucket_id()), ORDER_INDICES(entity.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  T& operator()(const stk::mesh::FastMeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& operator()(const MeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, component, fileName, lineNumber);
    return deviceData(deviceSelectedBucketOffset(index.bucket->bucket_id()), ORDER_INDICES(index.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index,
                                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, fileName, lineNumber);
    T* dataPtr = &deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, 0));
    const unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_component_stride());
  }

  template <typename Mesh>
  void set_all(const Mesh& ngpMesh, const T& value,
               const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER)
  {
    check_stale_field_access(fileName, lineNumber);
    Kokkos::deep_copy(deviceData, value);
    modify_on_device();
    anyPotentialDeviceFieldModification() = true;
  }

#else

  template <typename Mesh> KOKKOS_FUNCTION
  T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component);
  }

  KOKKOS_FUNCTION
  T& get(stk::mesh::FastMeshIndex entity, int component) const
  {
    return deviceData(deviceSelectedBucketOffset(entity.bucket_id), ORDER_INDICES(entity.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& get(MeshIndex entity, int component) const
  {
    return deviceData(deviceSelectedBucketOffset(entity.bucket->bucket_id()), ORDER_INDICES(entity.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  T& operator()(const stk::mesh::FastMeshIndex& index, int component) const
  {
    unsigned offset = deviceSelectedBucketOffset(index.bucket_id);
    return deviceData(offset, ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> KOKKOS_FUNCTION
  T& operator()(const MeshIndex& index, int component) const
  {
    return deviceData(deviceSelectedBucketOffset(index.bucket->bucket_id()), ORDER_INDICES(index.bucketOrd, component));
  }

  KOKKOS_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index) const
  {
    T* dataPtr = &deviceData(deviceSelectedBucketOffset(index.bucket_id), ORDER_INDICES(index.bucket_ord, 0));
    const unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_component_stride());
  }

  template <typename Mesh>
  void set_all(const Mesh& ngpMesh, const T& value)
  {
    Kokkos::deep_copy(deviceData, value);
    modify_on_device();
  }

#endif

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
  void swap(DeviceField<T> &other)
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
    ThrowRequire(needSyncToHost() == false);
#endif
    needSyncToDevice() = true;
  }

  void set_modify_on_device()
  {
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    ThrowRequire(needSyncToDevice() == false);
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

      copy_moved_host_bucket_data<FieldDataHostViewType, UnmanagedHostInnerView<T>>(hostData, hostData, oldBucketId, newBucketId, numPerEntity);
      copy_moved_device_bucket_data<FieldDataDeviceViewType, UnmanagedDevInnerView<T>>(deviceData, deviceData, oldBucketId, newBucketId, numPerEntity);
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
            copy_moved_host_bucket_data<FieldDataHostViewType, UnmanagedHostInnerView<T>>(hostData, hostData, oldBucketId, newBucketId, numPerEntity);
            copy_moved_device_bucket_data<FieldDataDeviceViewType, UnmanagedDevInnerView<T>>(deviceData, deviceData, oldBucketId, newBucketId, numPerEntity);
          } else {
            move_unmodified_buckets_in_range_from_back(buckets, numPerEntity, i);
          }
        }
      }
    }
  }

  void copy_unmodified_buckets(const BucketVector& buckets, FieldDataDeviceViewType destDevView, FieldDataHostViewType destHostView, unsigned numPerEntity)
  {
    for(unsigned i = 0; i < buckets.size(); i++) {
      if(hostFieldExistsOnBucket(buckets[i]->bucket_id())) {
        unsigned oldBucketId = buckets[i]->get_ngp_field_bucket_id(get_ordinal());
        unsigned newBucketId = buckets[i]->bucket_id();

        if(!buckets[i]->get_ngp_field_bucket_is_modified(get_ordinal())) {
          ThrowRequire(hostData.extent(0) != 0 && hostSelectedBucketOffset.extent(0) != 0);
          copy_moved_host_bucket_data<FieldDataHostViewType, UnmanagedHostInnerView<T>>(destHostView, hostData, oldBucketId, newBucketId, numPerEntity);
          copy_moved_device_bucket_data<FieldDataDeviceViewType, UnmanagedDevInnerView<T>>(destDevView, deviceData, oldBucketId, newBucketId, numPerEntity);
        }
      }
    }
  }

  void fill_host_view(Bucket* bucket, UnmanagedHostInnerView<T>& unHostInnerView, unsigned numPerEntity)
  {
#ifdef STK_DEBUG_FIELD_SYNC
      T* hostFieldPtr = static_cast<T*>(stk::mesh::ngp_debug_field_data(*hostField, *bucket));
#else
      T* hostFieldPtr = static_cast<T*>(stk::mesh::field_data(*hostField, *bucket));
#endif
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
    Kokkos::deep_copy(unInnerDestView, unInnerSrcView);
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
    Kokkos::deep_copy(unInnerDestView, unInnerSrcView);
  }

  void copy_bucket_from_device_to_host(Bucket* bucket, unsigned numPerEntity, unsigned numContiguousBuckets = 1)
  {
    unsigned selectedBucketOffset = newHostSelectedBucketOffset(bucket->bucket_id());

    T* devicePtr = deviceData.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unDeviceInnerView(devicePtr, ORDER_INDICES(bucketCapacity*numContiguousBuckets, numPerEntity));

    T* bufferPtr = deviceBuffer.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unBufferInnerView(bufferPtr, bucketCapacity*numContiguousBuckets, numPerEntity);

    transpose_contiguous_device_data_into_buffer(bucketCapacity*numContiguousBuckets, numPerEntity, unDeviceInnerView, unBufferInnerView);
    Kokkos::fence();

    UnmanagedHostInnerView<T> unHostInnerView(&hostData(selectedBucketOffset,0,0),
                                              bucketCapacity*numContiguousBuckets, numPerEntity);
    Kokkos::deep_copy(unHostInnerView, unBufferInnerView);
  }

  void copy_bucket_from_host_to_device(Bucket* bucket, unsigned numPerEntity, unsigned numContiguousBuckets = 1)
  {
    unsigned selectedBucketOffset = newHostSelectedBucketOffset(bucket->bucket_id());

    UnmanagedHostInnerView<T> unHostInnerView(&hostData(selectedBucketOffset,0,0),
                                              bucketCapacity*numContiguousBuckets, numPerEntity);

    fill_host_view(bucket, unHostInnerView, numPerEntity);

    T* bufferPtr = deviceBuffer.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unBufferInnerView(bufferPtr, bucketCapacity*numContiguousBuckets, numPerEntity);
    Kokkos::deep_copy(unBufferInnerView, unHostInnerView);

    T* devicePtr = deviceData.data() + selectedBucketOffset * bucketCapacity * numPerEntity;
    UnmanagedDevInnerView<T> unDeviceInnerView(devicePtr, ORDER_INDICES(bucketCapacity*numContiguousBuckets, numPerEntity));

    transpose_buffer_into_contiguous_device_data(bucketCapacity*numContiguousBuckets, numPerEntity, unBufferInnerView, unDeviceInnerView);
    Kokkos::fence();
  }

  void reset_device_buffer()
  {
    deviceBuffer = FieldDataDeviceUnmanagedViewType(reinterpret_cast<T*>(hostBulk->get_ngp_field_sync_buffer()),
                                                    hostData.extent(0), hostData.extent(1), hostData.extent(2));
  }

  void copy_new_and_modified_buckets_from_host(const BucketVector& buckets, unsigned numPerEntity)
  {
    reset_device_buffer();

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
#ifdef STK_DEBUG_FIELD_SYNC
    FieldDataDeviceViewType tempDataDeviceView = FieldDataDeviceViewType(name, numBuckets,
                                                                         ORDER_INDICES(bucketCapacity, numPerEntity));
    FieldDataHostViewType tempDataHostView = FieldDataHostViewType(name, numBuckets,
                                                                   bucketCapacity, numPerEntity);
#else
    FieldDataDeviceViewType tempDataDeviceView = FieldDataDeviceViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets,
                                                                         ORDER_INDICES(bucketCapacity, numPerEntity));
    FieldDataHostViewType tempDataHostView = FieldDataHostViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), numBuckets,
                                                                   bucketCapacity, numPerEntity);
#endif

    copy_unmodified_buckets(buckets, tempDataDeviceView, tempDataHostView, numPerEntity);

    deviceData = tempDataDeviceView;
    hostData = tempDataHostView;
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

    Kokkos::deep_copy(newDeviceSelectedBucketOffset, newHostSelectedBucketOffset);
  }

  void construct_field_exist_view(const BucketVector& allBuckets, const Selector& selector)
  {
    deviceFieldExistsOnBucket = BoolViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, hostField->name() + "_exists_on_bucket"),
                                             allBuckets.size());
    hostFieldExistsOnBucket = Kokkos::create_mirror_view(Kokkos::HostSpace(), deviceFieldExistsOnBucket, Kokkos::WithoutInitializing);
    Kokkos::deep_copy(hostFieldExistsOnBucket, false);
    for (size_t i = 0; i < allBuckets.size(); ++i) {
      if(selector(*allBuckets[i])) {
        hostFieldExistsOnBucket(allBuckets[i]->bucket_id()) = true;
      }
    }
    Kokkos::deep_copy(deviceFieldExistsOnBucket, hostFieldExistsOnBucket);
  }

  void construct_num_components_per_entity_view(const BucketVector & allBuckets) {
    if (deviceNumComponentsPerEntity.extent(0) != allBuckets.size()) {
      deviceNumComponentsPerEntity = UnsignedViewType(Kokkos::ViewAllocateWithoutInitializing(hostField->name() + "_numComponentsPerEntity"),
                                                      allBuckets.size());
      hostNumComponentsPerEntity = Kokkos::create_mirror_view(deviceNumComponentsPerEntity);
    }

    for (stk::mesh::Bucket * bucket : allBuckets) {
      hostNumComponentsPerEntity[bucket->bucket_id()] = stk::mesh::field_scalars_per_entity(*hostField, *bucket);
    }

    Kokkos::deep_copy(deviceNumComponentsPerEntity, hostNumComponentsPerEntity);
  }

  void copy_contiguous_buckets_from_device_to_host(const BucketVector& buckets, unsigned numPerEntity, unsigned& i)
  {
    unsigned endIndex = get_contiguous_bucket_offset_end(buckets, i);

    unsigned numCopyBucketCount = endIndex - i + 1;

    copy_bucket_from_device_to_host(buckets[i], numPerEntity, numCopyBucketCount);
    i = endIndex;
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
    reset_device_buffer();

    if (!userSpecifiedSelector) {
      transpose_all_device_data_into_buffer(*hostField, deviceData, deviceBuffer);
      Kokkos::fence();

      Kokkos::deep_copy(hostData, deviceBuffer);
    }
    else {
      copy_selected_buckets_to_host();
    }
    clear_sync_state_flags();
  }

  void copy_device_to_host()
  {
    if (hostField) {
      stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
      const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);

      if(userSpecifiedSelector) {
        selector &= *syncSelector;
      }

      sync_to_host_using_selector();

      hostField->increment_num_syncs_to_host();
      copy_data(buckets, [](T &hostFieldData, T &stkFieldData){
        stkFieldData = hostFieldData;
      }, selector);
    }
  }

  void copy_contiguous_buckets_from_host_to_device(const BucketVector& buckets, unsigned numPerEntity, unsigned& i)
  {
    unsigned endIndex = get_contiguous_bucket_offset_end(buckets, i);

    unsigned numCopyBucketCount = endIndex - i + 1;

    copy_bucket_from_host_to_device(buckets[i], numPerEntity, numCopyBucketCount);
    i = endIndex;
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
    reset_device_buffer();

    if (!userSpecifiedSelector) {
      Kokkos::deep_copy(deviceBuffer, hostData);

      transpose_buffer_into_all_device_data(*hostField, deviceBuffer, deviceData);
      Kokkos::fence();
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

#ifdef STK_DEBUG_FIELD_SYNC
      T* data = static_cast<T*>(stk::mesh::ngp_debug_field_data(*hostField, bucket));
#else
      T* data = static_cast<T*>(stk::mesh::field_data(*hostField, bucket));
#endif
      for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
      {
        for(unsigned j=0; j<numPerEntity; j++)
        {
          assigner(hostData(selectedBucketOffset, iEntity, j), data[iEntity*numPerEntity + j]);
        }
      }
    }
  }

  FieldDataHostViewType hostData;
  FieldDataDeviceViewType deviceData;
  FieldDataDeviceUnmanagedViewType deviceBuffer;
  Kokkos::View<bool, Kokkos::HostSpace> needSyncToHost;
  Kokkos::View<bool, Kokkos::HostSpace> needSyncToDevice;

  typename BoolViewType::HostMirror hostFieldExistsOnBucket;
  BoolViewType deviceFieldExistsOnBucket;

  typename UnsignedViewType::HostMirror hostSelectedBucketOffset;
  UnsignedViewType deviceSelectedBucketOffset;
  typename UnsignedViewType::HostMirror newHostSelectedBucketOffset;
  UnsignedViewType newDeviceSelectedBucketOffset;
  typename UnsignedViewType::HostMirror hostNumComponentsPerEntity;
  UnsignedViewType deviceNumComponentsPerEntity;

  stk::mesh::EntityRank rank;
  unsigned ordinal;

  const stk::mesh::BulkData* hostBulk;
  const stk::mesh::FieldBase* hostField;
  unsigned bucketCapacity;
  unsigned numBucketsForField;
  unsigned maxNumScalarsPerEntity;
  size_t synchronizedCount;

  DeviceField<T>* stateFields[MaximumFieldStates];

  Kokkos::View<int[1], Kokkos::HostSpace> copyCounter;
  bool userSpecifiedSelector;
  Selector* syncSelector;

#ifdef STK_DEBUG_FIELD_SYNC
  using FieldNameType = Kokkos::View<char*, UVMMemSpace>;
  template <typename U> using ScalarUvmType = Kokkos::View<U, UVMMemSpace>;

  FieldNameType fieldName;
  bool isFieldLayoutConsistent;
  ScalarUvmType<size_t> hostSynchronizedCount;
  ScalarUvmType<bool> lostDeviceFieldData;
  ScalarUvmType<bool> anyPotentialDeviceFieldModification;
  HostArrayType<unsigned> lastMeshModBucketId;
  LastFieldModLocationType lastFieldModLocation;
  FieldDataDeviceViewType lastFieldValue;
  typename UnsignedViewType::HostMirror debugHostSelectedBucketOffset;
  UnsignedViewType debugDeviceSelectedBucketOffset;
#endif
};

}
}

#endif
