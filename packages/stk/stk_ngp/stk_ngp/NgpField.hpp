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

#ifndef STK_NGP_NGPFIELD_H_
#define STK_NGP_NGPFIELD_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_ngp/NgpForEachEntity.hpp>
#include <stk_ngp/NgpMesh.hpp>
#include <stk_ngp/NgpProfilingBlock.hpp>

namespace ngp {

constexpr unsigned INVALID_ORDINAL = 9999999;

template<typename T> class ConstStkFieldAdapter;
template<typename T> class MultistateField;

class FieldBase
{
public:
  STK_FUNCTION FieldBase() = default;
  STK_FUNCTION virtual ~FieldBase() {}
  virtual void sync_to_host() {}
};


template<typename T>
class StkFieldAdapter : public FieldBase
{
public:
    typedef T value_type;

    StkFieldAdapter()
      : FieldBase(),
        field(nullptr) { }

    StkFieldAdapter(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f)
      : FieldBase(),
        field(&f) { }

    virtual ~StkFieldAdapter() = default;

    unsigned get_num_components_per_entity(stk::mesh::FastMeshIndex entity) const {
        return stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
    }

    T& get(const StkMeshAdapter& ngpMesh, stk::mesh::Entity entity, int component) const
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

    T& get(StkMeshAdapter::MeshIndex entity, int component) const
    {
        T* data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket->bucket_id(), entity.bucketOrd));
        ThrowAssert(data);
        return data[component];
    }

    void set_all(const StkMeshAdapter& ngpMesh, const T& value)
    {
      ngp::for_each_entity_run(ngpMesh, field->entity_rank(), *field, KOKKOS_LAMBDA(const StkMeshAdapter::MeshIndex& entity) {
                                 T* fieldPtr = static_cast<T*>(stk::mesh::field_data(*field, *entity.bucket, entity.bucketOrd));
                                 int numScalars = stk::mesh::field_scalars_per_entity(*field, *entity.bucket);
                                 for (int i=0; i<numScalars; i++) {
                                   fieldPtr[i] = value;
                                 }
                               });
    }

    void sync_to_host() override { }

    void sync_to_device() { }

    void modify_on_host() { }

    void modify_on_device() { }

    void clear_sync_state() { }

    void swap(ConstStkFieldAdapter<T> &sf) { }

    void swap(StkFieldAdapter<T> &sf) { }

    stk::mesh::EntityRank get_rank() const { return field->entity_rank(); }

    unsigned get_ordinal() const { return field->mesh_meta_data_ordinal(); }

#ifdef STK_HIDE_DEPRECATED_CODE
private:
#endif
    void copy_host_to_device() { };

    void copy_device_to_host() { };
#ifndef STK_HIDE_DEPRECATED_CODE
    void copy_host_to_device(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& f) { };

    void copy_device_to_host(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& f) { };
private:
#endif

    bool need_sync_to_host() const { return false; }

    bool need_sync_to_device() const { return false; }

    const stk::mesh::FieldBase * field;

    friend ConstStkFieldAdapter<T>;
    friend MultistateField<T>;
};

template<typename T>
class ConstStkFieldAdapter : public FieldBase
{
public:
    typedef T value_type;

    ConstStkFieldAdapter()
      : FieldBase(),
        stkFieldAdapter() { }

    ConstStkFieldAdapter(const StkFieldAdapter<T> &sfa)
      : FieldBase(),
        stkFieldAdapter(sfa)
    {
    }

    ConstStkFieldAdapter(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f)
      : FieldBase(),
        stkFieldAdapter(b, f)
    {
    }

    virtual ~ConstStkFieldAdapter() = default;

    const T& get(const StkMeshAdapter& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        return stkFieldAdapter.get(ngpMesh, entity, component);
    }

    const T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return stkFieldAdapter.get(entity, component);
    }

    const T& get(StkMeshAdapter::MeshIndex entity, int component) const
    {
        return stkFieldAdapter.get(entity, component);
    }

    stk::mesh::EntityRank get_rank() const { return stkFieldAdapter.get_rank(); }

    unsigned get_ordinal() const { return stkFieldAdapter.get_ordinal(); }

    void sync_to_host() override { }

    void sync_to_device() { }

    void swap(ConstStkFieldAdapter<T> &sf) { }

    void swap(StkFieldAdapter<T> &sf) { }

#ifdef STK_HIDE_DEPRECATED_CODE
private:
#endif
    void copy_host_to_device() { };

    void copy_device_to_host() { };
#ifndef STK_HIDE_DEPRECATED_CODE
    void copy_host_to_device(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { };

    void copy_device_to_host(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { };
private:
#endif

    bool need_sync_to_host() const { return false; }

    bool need_sync_to_device() const { return false; }

    void clear_sync_state() { }

    StkFieldAdapter<T> stkFieldAdapter;

    friend MultistateField<T>;
};


#ifdef KOKKOS_ENABLE_CUDA
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

template<typename T> class ConstStaticField;

template<typename T>
class StaticField : public FieldBase
{
public:
    typedef T value_type;

    STK_FUNCTION
    StaticField()
      : FieldBase(),
        rank(stk::topology::NODE_RANK),
        ordinal(INVALID_ORDINAL),
        hostBulk(nullptr),
        hostField(nullptr) { }

    void construct_view(const std::string& name, unsigned numBuckets, unsigned numPerEntity) {
#ifdef KOKKOS_ENABLE_CUDA
        deviceData = FieldDataType(name, numBuckets, numPerEntity);
#else
        deviceData = FieldDataType(name, numBuckets, bucketSize, numPerEntity);
#endif
        hostData = Kokkos::create_mirror_view(deviceData);
        fieldData = FieldDataDualViewType(deviceData, hostData);

        deviceFieldExistsOnBucket = BoolViewType(name + "_exists_on_bucket", numBuckets);
        hostFieldExistsOnBucket = Kokkos::create_mirror_view(deviceFieldExistsOnBucket);
    }

    StaticField(stk::mesh::EntityRank r, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
      : FieldBase(),
        rank(r),
        ordinal(INVALID_ORDINAL),
        hostBulk(&bulk),
        hostField(nullptr)
    {
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(rank);

        construct_view("no_name", allBuckets.size(), 1);

        create_device_field_data(initialValue);
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
      : FieldBase(),
        rank(field.entity_rank()),
        ordinal(field.mesh_meta_data_ordinal()),
        hostBulk(&bulk),
        hostField(&field)
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(field.entity_rank());
        unsigned numPerEntity=field.max_size(rank);

        construct_view("deviceData_"+field.name(), allBuckets.size(), numPerEntity);

        copy_data(buckets, [](T &hostFieldData, T &stkFieldData){hostFieldData = stkFieldData;});
        Kokkos::deep_copy(deviceData, hostData);

        Kokkos::deep_copy(hostFieldExistsOnBucket, false);
        for (size_t i = 0; i < buckets.size(); ++i) {
          hostFieldExistsOnBucket(buckets[i]->bucket_id()) = true;
        }
        Kokkos::deep_copy(deviceFieldExistsOnBucket, hostFieldExistsOnBucket);
    }

    void sync_to_host() override
    {
        if (need_sync_to_host()) {
            ProfilingBlock prof("copy_to_host for " + hostField->name());
            copy_device_to_host();
        }
    }

    void sync_to_device()
    {
        if (need_sync_to_device()) {
            ProfilingBlock prof("copy_to_device for " + hostField->name());
            copy_host_to_device();
        }
    }

    void modify_on_host()
    {
        fieldData.modify_host();  // New Kokkos API
    }

    void modify_on_device()
    {
        fieldData.modify_device();  // New Kokkos API
    }

    void clear_sync_state()
    {
        fieldData.clear_sync_state();  // New Kokkos API
    }

    STK_FUNCTION StaticField(const StaticField &) = default;

    STK_FUNCTION virtual ~StaticField() {}

    STK_FUNCTION
    unsigned get_num_components_per_entity(stk::mesh::FastMeshIndex entityIndex) const {
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
    void swap(StaticField<T> &sf)
    {
      swap_views(hostData,   sf.hostData);
      swap_views(deviceData, sf.deviceData);
      swap_views(fieldData,  sf.fieldData);
    }

#ifdef STK_HIDE_DEPRECATED_CODE
private:
#endif
    void copy_device_to_host()
    {
        fieldData.sync_host();  // New Kokkos API

        if (hostField) {
          stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
          const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
          copy_data(buckets, [](T &hostFieldData, T &stkFieldData){stkFieldData = hostFieldData;});
        }
    }

    void copy_host_to_device()
    {
        if (hostField) {
          stk::mesh::Selector selector = stk::mesh::selectField(*hostField);
          const stk::mesh::BucketVector& buckets = hostBulk->get_buckets(hostField->entity_rank(), selector);
          copy_data(buckets, [](T &hostFieldData, T &stkFieldData){hostFieldData = stkFieldData;});
        }

        fieldData.sync_device();  // New Kokkos API
    }

#ifndef STK_HIDE_DEPRECATED_CODE
    void copy_host_to_device(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { copy_host_to_device(); };

    void copy_device_to_host(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { copy_device_to_host(); };
private:
#endif
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
    typedef Kokkos::DualView<T**[bucketSize], Kokkos::LayoutRight,
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

    friend ConstStaticField<T>;
    friend MultistateField<T>;
};

template<typename T>
class ConstStaticField : public FieldBase
{
public:
    typedef T value_type;

    STK_FUNCTION
    ConstStaticField()
      : FieldBase()
    {
    }

    ConstStaticField(const StaticField<T> &sf)
      : FieldBase(),
        staticField(sf),
        constDeviceData(staticField.deviceData)
    {
    }

    ConstStaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
      : FieldBase(),
        staticField(bulk, field),
        constDeviceData(staticField.deviceData)
    {
    }

    STK_FUNCTION ConstStaticField(const ConstStaticField &) = default;

    STK_FUNCTION virtual ~ConstStaticField() {}

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
    stk::mesh::EntityRank get_rank() const { return staticField.get_rank(); }

    STK_FUNCTION
    unsigned get_ordinal() const { return staticField.get_ordinal(); }

    void sync_to_host() override
    {
        if (need_sync_to_host()) {
            copy_device_to_host();
        }
    }

    void sync_to_device()
    {
        if (need_sync_to_device()) {
            copy_host_to_device();
        }
    }

    STK_FUNCTION
    void swap(ConstStaticField<T> &sf)
    {
        swap(sf.staticField);
    }

    STK_FUNCTION
    void swap(StaticField<T> &sf)
    {
        staticField.swap(sf);
        constDeviceData = staticField.deviceData;
    }

#ifdef STK_HIDE_DEPRECATED_CODE
private:
#endif
    void copy_device_to_host()
    {
        staticField.modify_on_device();
        staticField.copy_device_to_host();
    }

    void copy_host_to_device()
    {
        staticField.modify_on_host();
        staticField.copy_host_to_device();
    }
#ifndef STK_HIDE_DEPRECATED_CODE
    void copy_host_to_device(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { copy_host_to_device(); };

    void copy_device_to_host(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase& field) { copy_device_to_host(); };
private:
#endif

    bool need_sync_to_host() const
    {
        return staticField.need_sync_to_host();
    }

    bool need_sync_to_device() const
    {
        return staticField.need_sync_to_device();
    }

    void clear_sync_state()
    {
        staticField.clear_sync_state();
    }

#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::View<
      const T * * [bucketSize], Kokkos::LayoutRight, Kokkos::CudaSpace,
      Kokkos::MemoryTraits<Kokkos::RandomAccess>> ConstFieldDataType;
#else
    typedef Kokkos::View<const T***, Kokkos::LayoutRight, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;
#endif

    StaticField<T> staticField;
    ConstFieldDataType constDeviceData;

    friend MultistateField<T>;
};

}


#endif /* STK_NGP_NGPFIELD_H_ */
