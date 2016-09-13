#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <ngp/NgpMesh.hpp>

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

    stk::mesh::EntityRank get_rank() const { return stkFieldAdapter.get_rank(); }

private:
    StkFieldAdapter<T> stkFieldAdapter;
};



template<typename T> class ConstStaticField;

template<typename T>
class StaticField {
public:
    typedef T value_type;

    StaticField() : rank(stk::topology::NODE_RANK) { }

    StaticField(stk::mesh::EntityRank r, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : deviceData(), rank(r)
    {
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(rank, selector);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(rank);
        fill_num_per_entity(buckets, allBuckets.size(), 1);
        fill_bucket_indices(buckets, allBuckets.size(), "bucketIndices");

        unsigned alloc_size = compute_alloc_size(buckets);
        create_device_field_data(alloc_size, initialValue);
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : deviceData(), rank(field.entity_rank())
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(field.entity_rank());
        fill_num_per_entity(buckets, allBuckets.size(), field);
        fill_bucket_indices(buckets, allBuckets.size(), "bucketIndices_"+field.name());

        allocate_view(deviceData, hostData, "deviceData_"+field.name(), compute_alloc_size(buckets));

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
        unsigned idx = get_index(ngpMesh, entity);
        return deviceData(idx+component);
    }

    STK_FUNCTION
    T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return deviceData(get_index(entity.bucket_id, entity.bucket_ord)+component);
    }

    template <typename MeshIndex> STK_FUNCTION
    T& get(MeshIndex entity, int component) const
    {
        return deviceData(get_index(entity.bucket->bucket_id(), entity.bucketOrd)+component);
    }

    STK_FUNCTION
    stk::mesh::EntityRank get_rank() const { return rank; }

private:
    template <typename Mesh> STK_FUNCTION
    unsigned get_index(const Mesh& ngpMesh, stk::mesh::Entity entity) const
    {
        stk::mesh::FastMeshIndex fast = ngpMesh.fast_mesh_index(entity);
        return get_index(fast.bucket_id, fast.bucket_ord);
    }

    STK_FUNCTION
    unsigned get_index(unsigned bucketId, unsigned bucketOrd) const
    {
        return get_index(bucketIndices(bucketId), bucketOrd, numPerEntity(bucketId));
    }

    unsigned get_index_host(unsigned bucketId, unsigned bucketOrd) const
    {
        return get_index(hostBucketIndices(bucketId), bucketOrd, hostNumPerEntity(bucketId));
    }

    STK_FUNCTION
    unsigned get_index(unsigned bucketIndex, unsigned bucketOrd, unsigned numPer) const
    {
        return bucketIndex + bucketOrd*numPer;
    }

    void fill_bucket_indices(const stk::mesh::BucketVector& buckets, size_t numAllBuckets, const std::string &name)
    {
        allocate_view(bucketIndices, hostBucketIndices, name, numAllBuckets);
        Kokkos::deep_copy(hostBucketIndices, std::numeric_limits<unsigned>::max());

        unsigned runningOffset = 0;
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
        {
            const stk::mesh::Bucket &bucket = *buckets[iBucket];
            unsigned scalarsPerEntity = hostNumPerEntity(bucket.bucket_id());
            hostBucketIndices(bucket.bucket_id()) = runningOffset;
            runningOffset += scalarsPerEntity*bucket.size();
        }

        Kokkos::deep_copy(bucketIndices, hostBucketIndices);
    }

    void fill_num_per_entity(const stk::mesh::BucketVector& buckets, size_t numAllBuckets, const stk::mesh::FieldBase &field)
    {
        allocate_num_per_entity(numAllBuckets);
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            hostNumPerEntity(buckets[iBucket]->bucket_id()) = stk::mesh::field_scalars_per_entity(field, *buckets[iBucket]);
        Kokkos::deep_copy(numPerEntity, hostNumPerEntity);
    }

    void fill_num_per_entity(const stk::mesh::BucketVector& buckets, size_t numAllBuckets, unsigned constNumPer)
    {
        allocate_num_per_entity(numAllBuckets);
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            hostNumPerEntity(buckets[iBucket]->bucket_id()) = constNumPer;
        Kokkos::deep_copy(numPerEntity, hostNumPerEntity);
    }

    void allocate_num_per_entity(size_t numAllBuckets)
    {
        allocate_view(numPerEntity, hostNumPerEntity, "numPerEntity", numAllBuckets);
        Kokkos::deep_copy(hostNumPerEntity, std::numeric_limits<unsigned>::max());
    }

    template <typename ViewType>
    void allocate_view(ViewType &view, typename ViewType::HostMirror &host, const std::string &name, size_t size)
    {
        view = ViewType(name, size);
        host = Kokkos::create_mirror_view(view);
    }

    unsigned compute_alloc_size(const stk::mesh::BucketVector& buckets)
    {
        unsigned alloc_size = 0;
        for(size_t i=0; i<buckets.size(); i++)
        {
            const stk::mesh::Bucket& bkt = *buckets[i];
            alloc_size += bkt.size()*hostNumPerEntity(bkt.bucket_id());
        }
        return alloc_size;
    }

    void create_device_field_data(unsigned allocSize, const T& initialValue)
    {
        deviceData = FieldDataType("deviceData", allocSize);
        hostData = Kokkos::create_mirror_view(deviceData);
        Kokkos::deep_copy(hostData, initialValue);
        Kokkos::deep_copy(deviceData, hostData);
    }

    template <typename Assigner>
    void copy_data(const stk::mesh::BucketVector& buckets, const stk::mesh::FieldBase &field, const Assigner &assigner)
    {
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
        {
            const stk::mesh::Bucket &bucket = *buckets[iBucket];
            unsigned numPer = hostNumPerEntity(bucket.bucket_id());
            T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
            for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
            {
                for(unsigned j=0; j<numPer; j++)
                {
                    unsigned index = get_index_host(bucket.bucket_id(), iEntity) + j;
                    assigner(hostData(index), data[iEntity*numPer + j]);
                }
            }
        }
    }

    typedef Kokkos::View<T*> FieldDataType;
    typedef Kokkos::View<unsigned*> UnsignedType;

    UnsignedType numPerEntity;
    typename UnsignedType::HostMirror hostNumPerEntity;

    UnsignedType bucketIndices;
    typename UnsignedType::HostMirror hostBucketIndices;

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
        unsigned idx = staticField.get_index(ngpMesh, entity);
        return constDeviceData(idx+component);
    }

    STK_FUNCTION
    const T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return constDeviceData(staticField.get_index(entity.bucket_id, entity.bucket_ord)+component);
    }

    template <typename MeshIndex> STK_FUNCTION
    const T& get(MeshIndex entity, int component) const
    {
        return constDeviceData(staticField.get_index(entity.bucket->bucket_id(), entity.bucketOrd)+component);
    }

    STK_FUNCTION
    stk::mesh::EntityRank get_rank() const { return staticField.rank; }

private:
    StaticField<T> staticField;

    typedef Kokkos::View<const T*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;
    ConstFieldDataType constDeviceData;
};

}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPFIELD_H_ */
