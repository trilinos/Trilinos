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
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"

namespace stk {
namespace mesh {

constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T> class ConstHostField;
template<typename T> class MultistateField;

class NgpFieldBase
{
public:
  STK_FUNCTION NgpFieldBase() = default;
  STK_FUNCTION virtual ~NgpFieldBase() {}
  virtual void update_field() = 0;
  virtual void modify_on_host() = 0;
  virtual void modify_on_device() = 0;
  virtual void sync_to_host() = 0;
  virtual void sync_to_device() = 0;
  virtual size_t synchronized_count() const = 0;
  virtual size_t num_syncs_to_host() const = 0;
  virtual size_t num_syncs_to_device() const = 0;
};

template<typename T>
class EntityFieldData {
public:
  STK_FUNCTION EntityFieldData(T* dataPtr, unsigned length, unsigned stride=1)
  : fieldDataPtr(dataPtr), fieldDataLength(length), fieldDataStride(stride)
  {}
  STK_FUNCTION ~EntityFieldData() = default;

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

  size_t num_syncs_to_host() const override { return field->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return field->num_syncs_to_device(); }

  unsigned get_stride() const { return 1; }

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
  }

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

  void sync_to_host() override { }

  void sync_to_device() override { }

  size_t synchronized_count() const override { return field->mesh_meta_data().mesh_bulk_data().synchronized_count(); }

  void clear_sync_state() { }

  void swap(ConstHostField<T> &sf) { }

  void swap(HostField<T> &sf) { }

  stk::mesh::EntityRank get_rank() const { return field->entity_rank(); }

  unsigned get_ordinal() const { return field->mesh_meta_data_ordinal(); }

private:
  void copy_host_to_device() { }

  void copy_device_to_host() { }

  bool need_sync_to_host() const { return false; }

  bool need_sync_to_device() const { return false; }

  const stk::mesh::FieldBase * field;

  friend ConstHostField<T>;
  friend MultistateField<T>;
};

template<typename T>
class ConstHostField : public NgpFieldBase
{
public:
  typedef T value_type;

  ConstHostField()
    : NgpFieldBase(),
      stkFieldAdapter() { }

  ConstHostField(const HostField<T> &sfa)
    : NgpFieldBase(),
      stkFieldAdapter(sfa)
  {
  }

  ConstHostField(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f)
    : NgpFieldBase(),
      stkFieldAdapter(b, f)
  {
  }

  virtual ~ConstHostField() = default;

  void update_field() override {}

  size_t num_syncs_to_host() const override { return stkFieldAdapter.num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return stkFieldAdapter.num_syncs_to_device(); }

  unsigned get_stride() const { return stkFieldAdapter.get_stride(); }

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& index) const {
    return stkFieldAdapter.get_num_components_per_entity(index);
  }

  const T& get(const HostMesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    return stkFieldAdapter.get(ngpMesh, entity, component);
  }

  const T& get(stk::mesh::FastMeshIndex entity, int component) const
  {
    return stkFieldAdapter.get(entity, component);
  }

  const T& get(HostMesh::MeshIndex entity, int component) const
  {
    return stkFieldAdapter.get(entity, component);
  }

  const T& operator()(stk::mesh::FastMeshIndex index, int component) const
  {
    return stkFieldAdapter(index, component);
  }

  const T& operator()(HostMesh::MeshIndex index, int component) const
  {
    return stkFieldAdapter(index, component);
  }

  stk::mesh::EntityRank get_rank() const { return stkFieldAdapter.get_rank(); }

  unsigned get_ordinal() const { return stkFieldAdapter.get_ordinal(); }

  void modify_on_host() override { }

  void modify_on_device() override { }

  void sync_to_host() override { }

  void sync_to_device() override { }

  size_t synchronized_count() const override { return stkFieldAdapter.synchronized_count(); }

  void swap(ConstHostField<T> &sf) { }

  void swap(HostField<T> &sf) { }

  void clear_sync_state() { }

private:
  void copy_host_to_device() { }

  void copy_device_to_host() { }

  bool need_sync_to_host() const { return false; }

  bool need_sync_to_device() const { return false; }

  HostField<T> stkFieldAdapter;

  friend MultistateField<T>;
};


#ifdef KOKKOS_ENABLE_CUDA
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

template<typename T> class ConstDeviceField;

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
  { }

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
    update_field();
  }

  void update_field() override
  {
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

  void sync_to_host() override
  {
    if (need_sync_to_host()) {
      ProfilingBlock prof("copy_to_host for " + hostField->name());
      copy_device_to_host();
    }
  }

  void sync_to_device() override
  {
    if (need_sync_to_device()) {
      ProfilingBlock prof("copy_to_device for " + hostField->name());
      if (hostBulk->synchronized_count() == synchronizedCount) {
        copy_host_to_device();
      }
      else {
        update_field();
      }
    }
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  void clear_sync_state()
  {
    fieldData.clear_sync_state();
  }

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

  template <typename Mesh>
  void set_all(const Mesh& ngpMesh, const T& value)
  {
    Kokkos::deep_copy(hostData, value);
    Kokkos::deep_copy(deviceData, value);
    clear_sync_state();
  }

  STK_FUNCTION
  stk::mesh::EntityRank get_rank() const { return rank; }

  STK_FUNCTION
  unsigned get_ordinal() const { return ordinal; }

  const stk::mesh::BulkData& get_bulk() const { return *hostBulk; }

  STK_FUNCTION
  void swap(DeviceField<T> &sf)
  {
    swap_views(hostData,   sf.hostData);
    swap_views(deviceData, sf.deviceData);
    swap_views(fieldData,  sf.fieldData);
  }

private:
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

  bool need_sync_to_host() const
  {
    return fieldData.need_sync_host();  // New Kokkos API
  }

  bool need_sync_to_device() const
  {
    return fieldData.need_sync_device();  // New Kokkos API
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

      T* data = static_cast<T*>(stk::mesh::field_data(*hostField, bucket));
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

  friend ConstDeviceField<T>;
  friend MultistateField<T>;
};

template<typename T>
class ConstDeviceField : public NgpFieldBase
{
public:
  typedef T value_type;

  STK_FUNCTION
  ConstDeviceField()
    : NgpFieldBase()
  {
  }

  ConstDeviceField(const DeviceField<T> &sf)
    : NgpFieldBase(),
      deviceField(sf),
      constDeviceData(deviceField.deviceData)
  {
  }

  ConstDeviceField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : NgpFieldBase(),
      deviceField(bulk, field),
      constDeviceData(deviceField.deviceData)
  {
  }

  STK_FUNCTION ConstDeviceField(const ConstDeviceField &) = default;

  STK_FUNCTION virtual ~ConstDeviceField() {}

  void update_field() override { deviceField.update_field(); }

  size_t num_syncs_to_host() const override { return deviceField.num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return deviceField.num_syncs_to_device(); }

  unsigned get_stride() const
  {
    unsigned stride = 1;
#ifdef KOKKOS_ENABLE_CUDA
    stride = deviceField.bucketCapacity;
#endif
    return stride;
  }

  STK_FUNCTION
  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entityIndex) const {
    return deviceField.get_num_components_per_entity(entityIndex);
  }
  
  template <typename Mesh> STK_FUNCTION
  T get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
    return get(fastIndex, component);
  }

  STK_FUNCTION
  T get(stk::mesh::FastMeshIndex entity, int component) const
  {
    return constDeviceData(entity.bucket_id, ORDER_INDICES(entity.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T get(MeshIndex entity, int component) const
  {
    return constDeviceData(entity.bucket->bucket_id(), ORDER_INDICES(entity.bucketOrd, component));
  }

  STK_FUNCTION
  T operator()(const stk::mesh::FastMeshIndex& index, int component) const
  {
    return constDeviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, component));
  }

  template <typename MeshIndex> STK_FUNCTION
  T operator()(const MeshIndex& index, int component) const
  {
    return constDeviceData(index.bucket->bucket_id(), ORDER_INDICES(index.bucketOrd, component));
  }

  STK_FUNCTION
  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index) const
  {
    T* dataPtr = constDeviceData(index.bucket_id, ORDER_INDICES(index.bucket_ord, 0));
    unsigned numScalars = get_num_components_per_entity(index);
    return EntityFieldData<T>(dataPtr, numScalars, get_stride());
  }

  STK_FUNCTION
  stk::mesh::EntityRank get_rank() const { return deviceField.get_rank(); }

  STK_FUNCTION
  unsigned get_ordinal() const { return deviceField.get_ordinal(); }

  void modify_on_host() override
  {
  }

  void modify_on_device() override
  {
  }

  void sync_to_host() override
  {
    if (need_sync_to_host()) {
      copy_device_to_host();
    }
  }

  void sync_to_device() override
  {
    if (need_sync_to_device()) {
      copy_host_to_device();
    }
  }

  size_t synchronized_count() const override { return deviceField.synchronized_count(); }

  STK_FUNCTION
  void swap(ConstDeviceField<T> &sf)
  {
    swap(sf.deviceField);
  }

  STK_FUNCTION
  void swap(DeviceField<T> &sf)
  {
    deviceField.swap(sf);
    constDeviceData = deviceField.deviceData;
  }

  void clear_sync_state()
  {
    deviceField.clear_sync_state();
  }

private:
  void copy_device_to_host()
  {
    deviceField.modify_on_device();
    deviceField.copy_device_to_host();
  }

  void copy_host_to_device()
  {
    deviceField.modify_on_host();
    deviceField.copy_host_to_device();
  }

  bool need_sync_to_host() const
  {
    return deviceField.need_sync_to_host();
  }

  bool need_sync_to_device() const
  {
    return deviceField.need_sync_to_device();
  }

#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::View<const T * * *, Kokkos::LayoutRight, Kokkos::CudaSpace,
                       Kokkos::MemoryTraits<Kokkos::RandomAccess>> ConstFieldDataType;
#else
  typedef Kokkos::View<const T***, Kokkos::LayoutRight,
                       Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;
#endif

  DeviceField<T> deviceField;
  ConstFieldDataType constDeviceData;

  friend MultistateField<T>;
};

}
}

#endif
