#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <ngp/NgpMesh.hpp>
#include <ngp/NgpForEachEntity.hpp>

namespace ngp {

template<typename T> class ConstStkFieldAdapter;

template<typename T>
class StkFieldAdapter
{
public:
    typedef T value_type;

    StkFieldAdapter() : field(nullptr) { }

    StkFieldAdapter(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f) : field(&f) { }

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
            *fieldPtr = value;
        });
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field_in)
    {
    }

    stk::mesh::EntityRank get_rank() const { return field->entity_rank(); }

private:
    const stk::mesh::FieldBase * field;

    friend ConstStkFieldAdapter<T>;
};

template<typename T>
class ConstStkFieldAdapter
{
public:
    typedef T value_type;

    ConstStkFieldAdapter() : stkFieldAdapter() { }

    ConstStkFieldAdapter(const StkFieldAdapter<T> &sfa)
    :   stkFieldAdapter(sfa)
    {
    }

    ConstStkFieldAdapter(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f)
    :   stkFieldAdapter(b, f)
    {
    }

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

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field_in)
    {
    }

    stk::mesh::EntityRank get_rank() const { return stkFieldAdapter.get_rank(); }

    void swap_data(ConstStkFieldAdapter<T> &sf)
    {
    }
    void swap_data(StkFieldAdapter<T> &sf)
    {
    }
private:
    StkFieldAdapter<T> stkFieldAdapter;
};


#ifdef KOKKOS_HAVE_CUDA
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

template<typename T> class ConstStaticField;

template<typename T>
class StaticField {
public:
    typedef T value_type;

    StaticField() : rank(stk::topology::NODE_RANK) { }

    void construct_view(const std::string& name, unsigned numBuckets, unsigned numPerEntity) {
#ifdef KOKKOS_HAVE_CUDA
        deviceData = FieldDataType(name, numBuckets, numPerEntity);
#else
        deviceData = FieldDataType(name, numBuckets, bucketSize, numPerEntity);
#endif
        hostData = Kokkos::create_mirror_view(deviceData);
    }

    StaticField(stk::mesh::EntityRank r, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : deviceData(), rank(r)
    {
      // const stk::mesh::BucketVector& buckets = bulk.get_buckets(rank, selector);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(rank);

        construct_view("no_name", allBuckets.size(), 1);

        create_device_field_data(initialValue);
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : deviceData(), rank(field.entity_rank())
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(field.entity_rank());
        unsigned numPerEntity=field.max_size(rank);

        construct_view("deviceData_"+field.name(), allBuckets.size(), numPerEntity);

        copy_data(buckets, field, [](T &staticData, T &fieldData){staticData = fieldData;});

        Kokkos::deep_copy(deviceData, hostData);
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        Kokkos::deep_copy(hostData, deviceData);

        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);
        copy_data(buckets, field, [](T &staticData, T &fieldData){fieldData = staticData;});
    }

    STK_FUNCTION ~StaticField(){}

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

    template <typename Mesh> STK_FUNCTION
    void set_all(const Mesh& ngpMesh, const T& value)
    {
        Kokkos::deep_copy(hostData, value);
        Kokkos::deep_copy(deviceData, value);
    }

    STK_FUNCTION
    stk::mesh::EntityRank get_rank() const { return rank; }

private:

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
    void copy_data(const stk::mesh::BucketVector& buckets, const stk::mesh::FieldBase &field, const Assigner &assigner)
    {
        unsigned numPerEntity = get_num_components_per_entity();
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
        {
            const stk::mesh::Bucket &bucket = *buckets[iBucket];

            T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
            for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
            {
                for(unsigned j=0; j<numPerEntity; j++)
                {
                    assigner(hostData(bucket.bucket_id(), ORDER_INDICES(iEntity,j) ), data[iEntity*numPerEntity + j]);
                }
            }
        }
    }

    unsigned get_num_components_per_entity() const {
#ifdef KOKKOS_HAVE_CUDA
        return hostData.dimension_1();
#else
        return hostData.dimension_2();
#endif
    }

#ifdef KOKKOS_HAVE_CUDA
    typedef Kokkos::View<T**[bucketSize], Kokkos::LayoutRight> FieldDataType;
#else
    typedef Kokkos::View<T***, Kokkos::LayoutRight> FieldDataType;
#endif

    typename FieldDataType::HostMirror hostData;
    FieldDataType deviceData;

    stk::mesh::EntityRank rank;

    friend ConstStaticField<T>;
};

template<typename T>
class ConstStaticField {
public:
    typedef T value_type;

    ConstStaticField() { }

    ConstStaticField(const StaticField<T> &sf) :
        staticField(sf),
        constDeviceData(staticField.deviceData)
    {
    }

    ConstStaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    :   staticField(bulk, field),
        constDeviceData(staticField.deviceData)
    {
    }

    STK_FUNCTION ~ConstStaticField(){}

    template <typename Mesh> STK_FUNCTION
    const T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        stk::mesh::FastMeshIndex fastIndex = ngpMesh.fast_mesh_index(entity);
        return get(fastIndex, component);
    }

    STK_FUNCTION
    const T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return constDeviceData(entity.bucket_id, ORDER_INDICES(entity.bucket_ord, component));
    }

    template <typename MeshIndex> STK_FUNCTION
    const T& get(MeshIndex entity, int component) const
    {
        return constDeviceData(entity.bucket->bucket_id(), ORDER_INDICES(entity.bucketOrd, component));
    }

    STK_FUNCTION
    stk::mesh::EntityRank get_rank() const { return staticField.rank; }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        staticField.copy_device_to_host(bulk, field);
    }

    STK_FUNCTION
    void swap_data(ConstStaticField<T> &sf)
    {
        swap_data(sf.staticField);
    }
    STK_FUNCTION
    void swap_data(StaticField<T> &sf)
    {
        StaticField<T> tmp = sf;
        sf = staticField;
        staticField = tmp;
        constDeviceData = staticField.deviceData;
    }

private:

#ifdef KOKKOS_HAVE_CUDA
    typedef Kokkos::View<const T**[bucketSize], Kokkos::LayoutRight, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;
#else
    typedef Kokkos::View<const T***, Kokkos::LayoutRight, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;
#endif

    StaticField<T> staticField;
    ConstFieldDataType constDeviceData;
};

}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_ */
