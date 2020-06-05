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

namespace stk {
namespace mesh {

constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T>
class EntityFieldData {
public:
  STK_FUNCTION EntityFieldData(T* dataPtr, unsigned length, unsigned stride=1)
  : fieldDataPtr(dataPtr), fieldDataLength(length), fieldDataStride(stride)
  {}
  KOKKOS_DEFAULTED_FUNCTION ~EntityFieldData() = default;

  STK_FUNCTION unsigned size() const { return fieldDataLength; }
  STK_FUNCTION T& operator[](unsigned idx) { return fieldDataPtr[idx*fieldDataStride]; }
  STK_FUNCTION const T& operator[](unsigned idx) const { return fieldDataPtr[idx*fieldDataStride]; }

private:
  T* fieldDataPtr;
  unsigned fieldDataLength;
  unsigned fieldDataStride;
};

template<typename T>
class HostField : public NgpFieldBase
{
public:
  typedef T value_type;

  HostField()
    : NgpFieldBase(),
      field(nullptr) { }

  HostField(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f)
    : NgpFieldBase(),
      field(&f) { }

  virtual ~HostField() = default;

  void update_field() override {}

  void set_field_states(HostField<T>* fields[]) {}

  size_t num_syncs_to_host() const override { return field->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return field->num_syncs_to_device(); }

  unsigned get_stride() const { return 1; }

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
  }

#ifdef STK_DEBUG_FIELD_SYNC
  void detect_device_field_modification() override {
  }

  void update_debug_storage(size_t hostSynchronizedCount) override {
  }

  bool any_device_field_modification() const override {
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

  void modify_on_host() override { }

  void modify_on_device() override { }

  void clear_sync_state() override {}

  void sync_to_host() override { }

  void sync_to_device() override { }

  size_t synchronized_count() const override { return field->mesh_meta_data().mesh_bulk_data().synchronized_count(); }

  FieldState state() const { return field->state(); }

  void rotate_multistate_data() override { }

  void swap(HostField<T> &sf) { }

  stk::mesh::EntityRank get_rank() const { return field->entity_rank(); }

  unsigned get_ordinal() const { return field->mesh_meta_data_ordinal(); }

protected:
  bool need_sync_to_host() const { return false; }

  bool need_sync_to_device() const { return false; }

private:
  void copy_host_to_device() { }

  void copy_device_to_host() { }

  const stk::mesh::FieldBase * field;
};


template<typename T>
class DeviceField : public NgpFieldBase
{
public:
  typedef T value_type;

  STK_FUNCTION
  DeviceField()
    : NgpFieldBase(),
      rank(stk::topology::NODE_RANK),
      ordinal(INVALID_ORDINAL),
      hostBulk(nullptr),
      hostField(nullptr),
      bucketCapacity(0),
      synchronizedCount(0)
  {
    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
  }

  DeviceField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : NgpFieldBase(),
      rank(field.entity_rank()),
      ordinal(field.mesh_meta_data_ordinal()),
      hostBulk(&bulk),
      hostField(&field),
      bucketCapacity(0),
      synchronizedCount(bulk.synchronized_count())
  {
    const int maxStates = static_cast<int>(stk::mesh::MaximumFieldStates);
    for (int s=0; s<maxStates; ++s) {
      stateFields[s] = nullptr;
    }
    update_field();
  }

  void construct_view(const std::string& name, unsigned numBuckets, unsigned numPerEntity) {
#ifdef KOKKOS_ENABLE_CUDA
    deviceData = FieldDataType(name, numBuckets, numPerEntity, bucketCapacity);
#else
    deviceData = FieldDataType(name, numBuckets, bucketCapacity, numPerEntity);
#endif
    hostData = Kokkos::create_mirror_view(deviceData);
    fieldData = FieldDataDualViewType(deviceData, hostData);

    deviceFieldExistsOnBucket = BoolViewType(name + "_exists_on_bucket", numBuckets);
    hostFieldExistsOnBucket = Kokkos::create_mirror_view(deviceFieldExistsOnBucket);

#ifdef STK_DEBUG_FIELD_SYNC
    fieldName = FieldNameType(hostField->name()+"_name", hostField->name().size()+1);
    hostSynchronizedCount = ScalarUvmType<size_t>(hostField->name()+"_hostSyncCount");
    anyDeviceFieldModification = ScalarUvmType<bool>(hostField->name()+"_anyDeviceFieldMod");
    anyPotentialDeviceFieldModification = ScalarUvmType<bool>(hostField->name()+"_anyPotentialDeviceFieldMod");
    construct_debug_views(numBuckets, numPerEntity);
#endif
  }

  void set_field_states(DeviceField<T>* fields[])
  {
    const unsigned numStates = hostField->number_of_states();
    for (unsigned state = 0; state < numStates; ++state) {
      stateFields[state] = fields[state];
    }   
  }

  void update_field() override
  {
    ProfilingBlock prof("update_field for " + hostField->name());
    stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
    const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
    if (!buckets.empty()) {
      bucketCapacity = buckets[0]->capacity();
    }

    const stk::mesh::BucketVector& allBuckets = hostBulk->buckets(hostField->entity_rank());
    unsigned numPerEntity = hostField->max_size(rank);

    construct_view("deviceData_"+hostField->name(), allBuckets.size(), numPerEntity);

    hostField->increment_num_syncs_to_device();
    copy_data(buckets, [](T &hostFieldData, T &stkFieldData){hostFieldData = stkFieldData;});
    Kokkos::deep_copy(deviceData, hostData);

#ifdef STK_DEBUG_FIELD_SYNC
    std::strcpy(fieldName.data(), hostField->name().data());
    Kokkos::deep_copy(lastFieldValue, hostData);
    reset_last_modification_state(LastModLocation::HOST | LastModLocation::DEVICE);
    hostSynchronizedCount() = hostBulk->synchronized_count();
    anyDeviceFieldModification() = false;
    anyPotentialDeviceFieldModification() = false;
    isFieldLayoutConsistent = true;
#endif

    Kokkos::deep_copy(hostFieldExistsOnBucket, false);
    for (size_t i = 0; i < buckets.size(); ++i) {
      hostFieldExistsOnBucket(buckets[i]->bucket_id()) = true;
    }
    Kokkos::deep_copy(deviceFieldExistsOnBucket, hostFieldExistsOnBucket);

    synchronizedCount = hostBulk->synchronized_count();
  }

  size_t num_syncs_to_host() const override { return hostField->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return hostField->num_syncs_to_device(); }

  void modify_on_host() override
  {
    fieldData.modify_host();
  }

  void modify_on_device() override
  {
    fieldData.modify_device();
  }

  void clear_sync_state() override
  {
    fieldData.clear_sync_state();
  }

  void sync_to_host() override
  {
    if (need_sync_to_host()) {
      ProfilingBlock prof("copy_to_host for " + hostField->name());
      copy_device_to_host();
#ifdef STK_DEBUG_FIELD_SYNC
      Kokkos::deep_copy(lastFieldValue, deviceData);
      reset_last_modification_state(LastModLocation::HOST | LastModLocation::DEVICE);
      anyDeviceFieldModification() = false;
      anyPotentialDeviceFieldModification() = false;
#endif
    }
  }

  void sync_to_device() override
  {
    if (need_sync_to_device()) {
      ProfilingBlock prof("copy_to_device for " + hostField->name());
      if (hostBulk->synchronized_count() == synchronizedCount) {
        copy_host_to_device();
#ifdef STK_DEBUG_FIELD_SYNC
        Kokkos::deep_copy(lastFieldValue, deviceData);
        reset_last_modification_state(LastModLocation::HOST | LastModLocation::DEVICE);
        anyDeviceFieldModification() = false;
        anyPotentialDeviceFieldModification() = false;
#endif
      }
      else {
        update_field();
      }
    }
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  STK_FUNCTION DeviceField(const DeviceField &) = default;

  STK_FUNCTION virtual ~DeviceField() {}

  unsigned get_stride() const
  {
    unsigned stride = 1;
#ifdef KOKKOS_ENABLE_CUDA
    stride = bucketCapacity;
#endif
    return stride;
  }

  STK_FUNCTION
  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entityIndex) const {
    const unsigned bucketId = entityIndex.bucket_id;
    if (deviceFieldExistsOnBucket[bucketId]) {
#ifdef KOKKOS_ENABLE_CUDA
      return deviceData.extent(1);
#else
      return deviceData.extent(2);
#endif
    }
    else {
      return 0;
    }
  }

#ifdef STK_DEBUG_FIELD_SYNC

  void construct_debug_views(unsigned numBuckets, unsigned numPerEntity) {
#ifdef KOKKOS_ENABLE_CUDA
    lastFieldValue = FieldDataType(hostField->name()+"_lastValue", numBuckets, numPerEntity, bucketCapacity);
    lastFieldModLocation = LastFieldModLocationType(hostField->name()+"_lastModLocation",
                                                    numBuckets, numPerEntity, bucketCapacity);
#else
    lastFieldValue = FieldDataType(hostField->name()+"_lastValue", numBuckets, bucketCapacity, numPerEntity);
    lastFieldModLocation = LastFieldModLocationType(hostField->name()+"_lastModLocation",
                                                    numBuckets, bucketCapacity, numPerEntity);
#endif
    hostField->set_last_modification_view(lastFieldModLocation);
  }

  bool any_device_field_modification() const override {
    return anyDeviceFieldModification();
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

    FieldDataType & localDeviceData = deviceData;
    FieldDataType & localLastFieldValue = lastFieldValue;
    LastFieldModLocationType & localLastFieldModLocation = lastFieldModLocation;

    stk::mesh::for_each_entity_run(ngpMesh, rank, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& index)
    {
#ifdef KOKKOS_ENABLE_CUDA
      const unsigned numComponents = localDeviceData.extent(1);
#else
      const unsigned numComponents = localDeviceData.extent(2);
#endif
      for (unsigned component = 0; component < numComponents; ++component) {
        if (localDeviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) !=
            localLastFieldValue(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)))
        {
          localLastFieldModLocation(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) = LastModLocation::DEVICE;
        }
      }
    });

    anyPotentialDeviceFieldModification() = false;
    Kokkos::fence();
  }

  void detect_any_device_field_modification() {
    FieldDataType & localDeviceData = deviceData;
    FieldDataType & localLastFieldValue = lastFieldValue;
    ScalarUvmType<bool> & localAnyDeviceFieldModification = anyDeviceFieldModification;
    localAnyDeviceFieldModification() = false;

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
      for (size_t i = 0; i < localDeviceData.extent(0); ++i) {
        for (size_t j = 0; j < localDeviceData.extent(1); ++j) {
          for (size_t k = 0; k < localDeviceData.extent(2); ++k) {
            if (localDeviceData(i, j, k) != localLastFieldValue(i, j, k)) {
              localAnyDeviceFieldModification() = true;
              return;
            }
          }
        }
      }
    });

    Kokkos::fence();
  }

  void update_debug_storage(size_t hostSyncCount) {
    hostSynchronizedCount() = hostSyncCount;

    if (isFieldLayoutConsistent) {
      detect_any_device_field_modification();
      isFieldLayoutConsistent = false;
    }

    const unsigned numBuckets = hostBulk->buckets(hostField->entity_rank()).size();
    const unsigned numPerEntity = hostField->max_size(rank);
    construct_debug_views(numBuckets, numPerEntity);

    if (anyDeviceFieldModification()) {
      reset_last_modification_state(LastModLocation::DEVICE);
    }
    else {
      reset_last_modification_state(LastModLocation::HOST | LastModLocation::DEVICE);
    }
  }

  void reset_last_modification_state(uint8_t value) {
    Kokkos::deep_copy(lastFieldModLocation, static_cast<LastModLocation>(value));
  }

  STK_INLINE_FUNCTION
  bool last_modification_not_on_device(const stk::mesh::FastMeshIndex & index, int component) const {
    return !(lastFieldModLocation(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) &
             LastModLocation::DEVICE);
  }

  STK_INLINE_FUNCTION
  bool last_modification_not_on_device(const stk::mesh::DeviceMesh::MeshIndex & index, int component) const {
    return !(lastFieldModLocation(index.bucket->bucket_id(), ORDER_INDICES(index.bucketOrd, component)) &
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
  void print_unsynchronized_field_warning(const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing unsynchronized Field %s on Device after mesh modification\n", fieldName.data());
    }
    else {
      printf("%s:%i *** WARNING: Accessing unsynchronized Field %s on Device after mesh modification\n",
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
  void check_stale_field_access(const stk::mesh::FastMeshIndex & index, int component, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unsynchronized_field_warning(fileName, lineNumber);
      return;
    }

    if (last_modification_not_on_device(index, component)) {
      print_stale_data_warning(index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const stk::mesh::DeviceMesh::MeshIndex & index, int component, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unsynchronized_field_warning(fileName, lineNumber);
      return;
    }

    if (last_modification_not_on_device(index, component)) {
      print_stale_data_warning(index.bucket->bucket_id(), index.bucketOrd, component, fileName, lineNumber);
    }
  }

  STK_INLINE_FUNCTION
  void check_stale_field_access(const stk::mesh::FastMeshIndex & index, const char * fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod()) {
      print_unsynchronized_field_warning(fileName, lineNumber);
      return;
    }

    for (unsigned component = 0; component < get_num_components_per_entity(index); ++component) {
      if (last_modification_not_on_device(index, component)) {
        print_stale_data_warning(index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
      }
    }
  }

  template <typename Mesh> STK_FUNCTION
  T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component, fileName, lineNumber);
  }

  STK_FUNCTION
  T& get(stk::mesh::FastMeshIndex entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(entity, component, fileName, lineNumber);
    return deviceData(entity.bucket_id, ORDER_INDICES(entity.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T& get(MeshIndex entity, int component,
         const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(entity, component, fileName, lineNumber);
    return deviceData(entity.bucket->bucket_id(), ORDER_INDICES(entity.bucketOrd, component));
  }

  STK_FUNCTION
  T& operator()(const stk::mesh::FastMeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, component, fileName, lineNumber);
    return deviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T& operator()(const MeshIndex& index, int component,
                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, component, fileName, lineNumber);
    return deviceData(index.bucket->bucket_id(), ORDER_INDICES(index.bucketOrd, component));
  }

  STK_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index,
                                const char * fileName = DEVICE_DEBUG_FILE_NAME, int lineNumber = DEVICE_DEBUG_LINE_NUMBER) const
  {
    check_stale_field_access(index, fileName, lineNumber);
    T* dataPtr = &deviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, 0));
    const unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_stride());
  }

#else

  template <typename Mesh> STK_FUNCTION
  T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component);
  }

  STK_FUNCTION
  T& get(stk::mesh::FastMeshIndex entity, int component) const
  {
    return deviceData(entity.bucket_id, ORDER_INDICES(entity.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T& get(MeshIndex entity, int component) const
  {
    return deviceData(entity.bucket->bucket_id(), ORDER_INDICES(entity.bucketOrd, component));
  }

  STK_FUNCTION
  T& operator()(const stk::mesh::FastMeshIndex& index, int component) const
  {
    return deviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T& operator()(const MeshIndex& index, int component) const
  {
    return deviceData(index.bucket->bucket_id(), ORDER_INDICES(index.bucketOrd, component));
  }

  STK_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index) const
  {
    T* dataPtr = &deviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, 0));
    const unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_stride());
  }

#endif

  template <typename Mesh>
  void set_all(const Mesh& ngpMesh, const T& value)
  {
    Kokkos::deep_copy(hostData, value);
    Kokkos::deep_copy(deviceData, value);
    clear_sync_state();
#ifdef STK_DEBUG_FIELD_SYNC
    anyPotentialDeviceFieldModification() = true;
#endif
  }

  STK_FUNCTION
  stk::mesh::EntityRank get_rank() const { return rank; }

  STK_FUNCTION
  unsigned get_ordinal() const { return ordinal; }

  const stk::mesh::BulkData& get_bulk() const { return *hostBulk; }

  FieldState state() const { return hostField->state(); }

  void rotate_multistate_data() override
  {
    //This is code brought down from BulkData::update_field_data_states,
    //the 'real' multistate rotation still needs to be implemented for
    //DeviceField. To improve performance, it should be possible to
    //rotate DeviceField views on device, without the modify/sync operation.

    const int numStates = hostField->number_of_states();
    if (numStates > 1) {
      for (int s=0; s<numStates; ++s) {
        FieldState stateEnum = static_cast<FieldState>(s);
        FieldBase* fieldOfState = hostField->field_state(stateEnum);
        fieldOfState->sync_to_host();
        fieldOfState->modify_on_host();
      }
    }
  }

  STK_FUNCTION
  void swap(DeviceField<T> &sf)
  {
    swap_views(hostData,   sf.hostData);
    swap_views(deviceData, sf.deviceData);
    swap_views(fieldData,  sf.fieldData);
  }

protected:
  bool need_sync_to_host() const
  {
    return fieldData.need_sync_host();
  }

  bool need_sync_to_device() const
  {
    return fieldData.need_sync_device();
  }

  void copy_device_to_host()
  {
    fieldData.sync_host();  // New Kokkos API

    if (hostField) {
      stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
      const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
      hostField->increment_num_syncs_to_host();
      copy_data(buckets, [](T &hostFieldData, T &stkFieldData){stkFieldData = hostFieldData;});
    }
  }

  void copy_host_to_device()
  {
    if (hostField) {
      stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
      const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
      hostField->increment_num_syncs_to_device();
      copy_data(buckets, [](T &hostFieldData, T &stkFieldData){hostFieldData = stkFieldData;});
    }

    fieldData.sync_device();  // New Kokkos API
  }

  template <typename ViewType>
  STK_FUNCTION
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

  void create_device_field_data(const T& initialValue)
  {
    Kokkos::deep_copy(hostData, initialValue);
    Kokkos::deep_copy(deviceData, hostData);
  }

  template <typename Assigner>
  void copy_data(const stk::mesh::BucketVector& buckets, const Assigner &assigner)
  {
    unsigned numPerEntity = get_max_num_components_per_entity();
    for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
    {
      const stk::mesh::Bucket &bucket = *buckets[iBucket];

#ifdef STK_DEBUG_FIELD_SYNC
      T* data = static_cast<T*>(stk::mesh::ngp_debug_field_data(*hostField, bucket));
#else
      T* data = static_cast<T*>(stk::mesh::field_data(*hostField, bucket));
#endif
      for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
      {
        for(unsigned j=0; j<numPerEntity; j++)
        {
          assigner(hostData(bucket.bucket_id(), ORDER_INDICES(iEntity,j) ), data[iEntity*numPerEntity + j]);
        }
      }
    }
  }

  unsigned get_max_num_components_per_entity() const {
#ifdef KOKKOS_ENABLE_CUDA
    return hostData.extent(1);
#else
    return hostData.extent(2);
#endif
  }

#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::DualView<T***, Kokkos::LayoutRight,
                           Kokkos::CudaSpace> FieldDataDualViewType;
#else
  typedef Kokkos::DualView<T***, Kokkos::LayoutRight> FieldDataDualViewType;
#endif

  typedef typename FieldDataDualViewType::t_dev FieldDataType;
  typedef typename FieldDataDualViewType::t_host FieldDataHostType;

  FieldDataHostType hostData;
  FieldDataType deviceData;
  FieldDataDualViewType fieldData;

  typename BoolViewType::HostMirror hostFieldExistsOnBucket;
  BoolViewType deviceFieldExistsOnBucket;

  stk::mesh::EntityRank rank;
  unsigned ordinal;

  const stk::mesh::BulkData* hostBulk;
  const stk::mesh::FieldBase* hostField;
  unsigned bucketCapacity;
  size_t synchronizedCount;

  DeviceField<T>* stateFields[MaximumFieldStates];

#ifdef STK_DEBUG_FIELD_SYNC
  using FieldNameType = Kokkos::View<char*, UVMMemSpace>;
  template <typename U> using ScalarUvmType = Kokkos::View<U, UVMMemSpace>;

  FieldNameType fieldName;
  bool isFieldLayoutConsistent;
  ScalarUvmType<size_t> hostSynchronizedCount;
  ScalarUvmType<bool> anyDeviceFieldModification;
  ScalarUvmType<bool> anyPotentialDeviceFieldModification;
  LastFieldModLocationType lastFieldModLocation;
  FieldDataType lastFieldValue;
#endif
};

}
}

#endif
