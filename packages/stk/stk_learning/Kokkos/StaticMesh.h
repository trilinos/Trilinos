#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkVector.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include <string>

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Field.hpp"

namespace stk {
namespace mesh {
struct ConstMeshIndex
{
    const stk::mesh::Bucket *bucket;
    size_t bucketOrd;
};
}
}

#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif

#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda     ExecSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP   HostExecSpace ;
#elif KOKKOS_HAVE_CUDA
  typedef Kokkos::Serial   HostExecSpace ;
#else
  typedef Kokkos::Serial   HostExecSpace ;
#endif

#ifdef KOKKOS_HAVE_OPENMP
   typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_HAVE_CUDA
   typedef Kokkos::CudaSpace    MemSpace;
#else
   typedef Kokkos::HostSpace    MemSpace;
#endif

#ifdef KOKKOS_HAVE_OPENMP
typedef Kokkos::OpenMP       UVMMemSpace;
#elif KOKKOS_HAVE_CUDA
typedef Kokkos::CudaUVMSpace UVMMemSpace;
#else
typedef Kokkos::HostSpace    UVMMemSpace;
#endif

namespace ngp {

#ifdef KOKKOS_HAVE_CUDA
typedef Kokkos::LayoutLeft   Layout ;
#else
typedef Kokkos::LayoutRight   Layout ;
#endif

typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;

template <typename Mesh, typename AlgorithmPerEntity>
struct ThreadFunctor
{
    STK_FUNCTION
    ThreadFunctor(const typename Mesh::BucketType *b, const AlgorithmPerEntity &f) :
        bucket(b),
        functor(f)
    {}
    STK_FUNCTION
    void operator()(const int& i) const
    {
        functor(typename Mesh::MeshIndex{bucket, static_cast<unsigned>(i)});
    }
    const typename Mesh::BucketType *bucket;
    const AlgorithmPerEntity &functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
struct TeamFunctor
{
    STK_FUNCTION
    TeamFunctor(const Mesh m, const stk::mesh::EntityRank r, stk::Vector<unsigned> b, const AlgorithmPerEntity f) :
        mesh(m),
        rank(r),
        bucketIds(b),
        functor(f)
    {}
    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ScheduleType>::member_type TeamHandleType;
    STK_FUNCTION
    void operator()(const TeamHandleType& team) const
    {
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
        unsigned numElements = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements),
                             ThreadFunctor<Mesh, AlgorithmPerEntity>(&bucket, functor));
    }
    const Mesh mesh;
    const stk::mesh::EntityRank rank;
    stk::Vector<unsigned> bucketIds;
    const AlgorithmPerEntity functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
void for_each_entity_run(Mesh &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const AlgorithmPerEntity &functor)
{
    stk::Vector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
    unsigned numBuckets = bucketIds.size();
    Kokkos::parallel_for(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO),
                         TeamFunctor<Mesh, AlgorithmPerEntity>(mesh, rank, bucketIds, functor));
}

//    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ScheduleType>::member_type TeamHandleType;
//    unsigned numBuckets = mesh.num_buckets(rank);
//    Kokkos::parallel_for(Kokkos::TeamPolicy<MyExecSpace>(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
//    {
//        const int bucketIndex = team.league_rank();
//        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
//        unsigned numElements = bucket.size();
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&](const int& i)
//        {
//            functor(typename Mesh::MeshIndex{&bucket, static_cast<unsigned>(i)});
//        });
//    });

template <typename Mesh, typename AlgorithmPerEntity>
void for_each_entity_run_stkmesh(Mesh &mesh, stk::topology::rank_t rank, const AlgorithmPerEntity &functor)
{
    typedef typename Kokkos::TeamPolicy<HostExecSpace, ScheduleType>::member_type TeamHandleType;
    const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part());
    unsigned numBuckets = elemBuckets.size();
    Kokkos::parallel_for(Kokkos::TeamPolicy<HostExecSpace>(numBuckets, Kokkos::AUTO), [&](const TeamHandleType& team)
    {
        const int bucketIndex = team.league_rank();
        const stk::mesh::Bucket &bucket = *elemBuckets[bucketIndex];
        unsigned numElements = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
        {
            functor(stk::mesh::FastMeshIndex{bucket.bucket_id(), static_cast<unsigned>(i)});
        });
    });
}

template <typename T>
class Entities
{
public:
    STK_FUNCTION
    Entities(T e, unsigned n) : entities(e), num(n)
    {
    }
    STK_FUNCTION
    stk::mesh::Entity operator[](unsigned i) const
    {
        return entities[i];
    }
    STK_FUNCTION
    unsigned size() const
    {
        return num;
    }
private:
    T entities;
    unsigned num;
};


inline stk::Vector<unsigned> get_bucket_ids(const stk::mesh::BulkData &bulk,
                     stk::mesh::EntityRank rank,
                     const stk::mesh::Selector &selector)
{
    const stk::mesh::BucketVector &buckets = bulk.get_buckets(rank, selector);
    stk::Vector<unsigned> bucketIds(buckets.size());
    for(size_t i=0; i<buckets.size(); i++)
        bucketIds[i] = buckets[i]->bucket_id();
    bucketIds.copy_host_to_device();
    return bucketIds;
}


class WrapperMesh
{
public:
    typedef HostExecSpace MeshExecSpace;

    typedef stk::mesh::ConstMeshIndex MeshIndex;
    typedef stk::mesh::Bucket BucketType;
    typedef Entities<const stk::mesh::Entity *> ConnectedNodes;

    WrapperMesh(const stk::mesh::BulkData& b) : bulk(b)
    {
    }

    ConnectedNodes get_nodes(const stk::mesh::ConstMeshIndex &elem) const
    {
        return ConnectedNodes(elem.bucket->begin_nodes(elem.bucketOrd), elem.bucket->num_nodes(elem.bucketOrd));
    }

    stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
    {
        const stk::mesh::MeshIndex &meshIndex = bulk.mesh_index(entity);
        return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), static_cast<unsigned>(meshIndex.bucket_ordinal)};
    }

    stk::Vector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector)
    {
        return ngp::get_bucket_ids(bulk, rank, selector);
    }

    unsigned num_buckets(stk::mesh::EntityRank rank) const
    {
        return bulk.buckets(rank).size();
    }

    const stk::mesh::Bucket & get_bucket(stk::mesh::EntityRank rank, unsigned i) const
    {
        return *bulk.buckets(rank)[i];
    }

private:
    const stk::mesh::BulkData &bulk;
};

template<typename T>
class WrapperField
{
public:
    WrapperField(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f) : field(f)
    {
    }

    T& get(const WrapperMesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        T *data = static_cast<T *>(stk::mesh::field_data(field, entity));
        return data[component];
    }

    const T& const_get(const WrapperMesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        const T *data = static_cast<T *>(stk::mesh::field_data(field, entity));
        return data[component];
    }

    T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        T *data = static_cast<T *>(stk::mesh::field_data(field, entity.bucket_id, entity.bucket_ord));
        return data[component];
    }

    const T& const_get(stk::mesh::FastMeshIndex entity, int component) const
    {
        const T *data = static_cast<T *>(stk::mesh::field_data(field, entity.bucket_id, entity.bucket_ord));
        return data[component];
    }

    T& get(WrapperMesh::MeshIndex entity, int component) const
    {
        T *data = static_cast<T *>(stk::mesh::field_data(field, entity.bucket->bucket_id(), entity.bucketOrd));
        return data[component];
    }

    const T& const_get(WrapperMesh::MeshIndex entity, int component) const
    {
        const T *data = static_cast<T *>(stk::mesh::field_data(field, entity.bucket->bucket_id(), entity.bucketOrd));
        return data[component];
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
    }

private:
    const stk::mesh::FieldBase& field;
};




typedef Kokkos::View<stk::mesh::Entity*, Kokkos::LayoutStride, MemSpace> ConnectedNodesType;
typedef Kokkos::View<stk::mesh::Entity*, MemSpace> EntityViewType;
typedef Kokkos::View<stk::mesh::Entity**, Layout, MemSpace> BucketConnectivityType;

struct StaticBucket {
    STK_FUNCTION
    StaticBucket()
     : bucketId(0), entityRank(stk::topology::NODE_RANK), entities(), connectivity() {}

    void initialize(unsigned bucket_id, stk::mesh::EntityRank rank, unsigned numEntities, unsigned numNodesPerEntity) {
        bucketId = bucket_id;
        entityRank = rank;
        entities = EntityViewType("BucketEntities"+std::to_string(bucket_id), numEntities);
        hostEntities = Kokkos::create_mirror_view(entities);
        connectivity = BucketConnectivityType("BucketConnectivity"+std::to_string(bucket_id), numEntities, numNodesPerEntity);
        hostConnectivity = Kokkos::create_mirror_view(connectivity);
    }

    void copy_to_device()
    {
        Kokkos::deep_copy(entities, hostEntities);
        Kokkos::deep_copy(connectivity, hostConnectivity);
    }

    STK_FUNCTION
    unsigned bucket_id() const { return bucketId; }

    STK_FUNCTION
    size_t size() const { return entities.size(); }

    STK_FUNCTION
    stk::mesh::EntityRank entity_rank() const { return entityRank; }

    STK_FUNCTION
    unsigned get_num_nodes_per_entity() const { return connectivity.extent(1); }

    STK_FUNCTION
    ConnectedNodesType get_nodes(unsigned offsetIntoBucket) const {
        return Kokkos::subview(connectivity, offsetIntoBucket, Kokkos::ALL());
    }

    STK_FUNCTION
    stk::mesh::Entity operator[](unsigned offsetIntoBucket) const {
        return entities(offsetIntoBucket);
    }

    stk::mesh::Entity host_get_entity(unsigned offsetIntoBucket) const {
        return hostEntities(offsetIntoBucket);
    }

    unsigned bucketId;
    stk::mesh::EntityRank entityRank;

    EntityViewType entities;
    EntityViewType::HostMirror hostEntities;

    BucketConnectivityType connectivity;
    BucketConnectivityType::HostMirror hostConnectivity;
};

struct StaticMeshIndex
{
    const StaticBucket *bucket;
    size_t bucketOrd;
};

class StaticMesh
{
public:
    typedef ExecSpace MeshExecSpace;

    typedef StaticMeshIndex MeshIndex;
    typedef StaticBucket BucketType;
    typedef Entities<ConnectedNodesType> ConnectedNodes;

    StaticMesh(const stk::mesh::BulkData& b) : bulk(b)
    {
        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror("host_mesh_indices", bulk.get_size_of_entity_index_space());
        fill_mesh_indices(stk::topology::NODE_RANK, bulk);
        fill_mesh_indices(stk::topology::ELEM_RANK, bulk);
        fill_buckets(bulk, bulk.mesh_meta_data().universal_part());
        copy_mesh_indices_to_device();
    }

    ~StaticMesh()
    {
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(const StaticMeshIndex &elem) const
    {
        ConnectedNodesType nodes = buckets(elem.bucket->bucket_id()).get_nodes(elem.bucketOrd);
        return ConnectedNodes(nodes, nodes.size());
    }

    STK_FUNCTION
    stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
    {
        return device_mesh_index(entity);
    }

    STK_FUNCTION
    const stk::mesh::FastMeshIndex& device_mesh_index(stk::mesh::Entity entity) const
    {
        return deviceMeshIndices(entity.local_offset());
    }

    const stk::mesh::FastMeshIndex& host_mesh_index(stk::mesh::Entity entity) const
    {
        return hostMeshIndices(entity.local_offset());
    }

    stk::Vector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector)
    {
        return ngp::get_bucket_ids(get_bulk_on_host(), rank, selector);
    }

    STK_FUNCTION
    unsigned num_buckets(stk::mesh::EntityRank rank) const
    {
        return buckets.size();
    }

    STK_FUNCTION
    const StaticBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
    {
        return buckets(index);
    }

    void clear()
    {
        buckets = BucketView();
    }


private:

    void fill_buckets(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);
        unsigned numElementBuckets = elemBuckets.size();

        buckets = BucketView("ElementNodeConnectivity", numElementBuckets);

        for (unsigned elemBucketIndex = 0; elemBucketIndex < numElementBuckets; ++elemBucketIndex)
        {
            const stk::mesh::Bucket& bucket = *elemBuckets[elemBucketIndex];
            unsigned nodesPerElem = bucket.topology().num_nodes();
            StaticBucket& bkt = buckets(bucket.bucket_id());
            bkt.initialize(bucket.bucket_id(), stk::topology::ELEM_RANK, bucket.size(), nodesPerElem);
            unsigned bucket_id = bucket.bucket_id();

            for(unsigned elemIndex = 0; elemIndex < bucket.size(); ++elemIndex)
            {
                stk::mesh::Entity element = bucket[elemIndex];
                buckets(bucket_id).hostEntities(elemIndex) = element;
                const stk::mesh::Entity * elemNodes = bulk.begin_nodes(element);
                for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
                {
                    buckets(bucket_id).hostConnectivity(elemIndex, iNode) = elemNodes[iNode];
                }
            }
        }
        for(unsigned b=0; b<buckets.size(); ++b) {
            buckets(b).copy_to_device();
        }
    }

    void fill_mesh_indices(stk::mesh::EntityRank rank, const stk::mesh::BulkData& bulk)
    {
        const stk::mesh::BucketVector& bkts = bulk.buckets(rank);

        for(const stk::mesh::Bucket* bktptr : bkts)
        {
            const stk::mesh::Bucket& bkt = *bktptr;
            for(size_t i = 0; i < bkt.size(); ++i)
            {
                hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex{bkt.bucket_id(), static_cast<unsigned>(i)};
            }
        }
    }

    void copy_mesh_indices_to_device()
    {
        unsigned length = hostMeshIndices.size();
        Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> tmp_device_mesh_indices("tmp_dev_mesh_indices", length);
        Kokkos::deep_copy(tmp_device_mesh_indices, hostMeshIndices);
        deviceMeshIndices = tmp_device_mesh_indices;
    }

    const stk::mesh::BulkData &get_bulk_on_host()
    {
        return bulk;
    }


    typedef Kokkos::View<StaticBucket*, UVMMemSpace> BucketView;
    const stk::mesh::BulkData &bulk;
    BucketView buckets;
    Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
    Kokkos::View<const stk::mesh::FastMeshIndex*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;
    //    Kokkos::View<const stk::mesh::FastMeshIndex*> device_mesh_indices;
};



template<typename T>
class StaticField {
public:
    typedef Kokkos::View<T**> BucketFieldDataType;
    typedef Kokkos::View<BucketFieldDataType*, UVMMemSpace> FieldDataType;
    typedef Kokkos::View<const T**, Kokkos::MemoryTraits<Kokkos::RandomAccess>> ConstBucketFieldDataType;
    typedef Kokkos::View<ConstBucketFieldDataType*, UVMMemSpace> ConstFieldDataType;
//    typedef Kokkos::View<T***> FieldDataType;
//    typedef Kokkos::View<const T***, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstFieldDataType;

    StaticField(stk::mesh::EntityRank rank, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : deviceData()
    {
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(rank);
        if(!allBuckets.empty())
        {
            size_t sizePerBucket = allBuckets[0]->capacity();
            create_device_field_data(allBuckets.size(), sizePerBucket, 1, initialValue);
        }
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : deviceData()
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& allBuckets = bulk.buckets(field.entity_rank());
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
            deviceData = FieldDataType("deviceData_"+field.name(), allBuckets.size());
            constDeviceData = ConstFieldDataType("constDeviceData_"+field.name(), allBuckets.size());

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];

                const std::string name = "bucket_" + std::to_string(iBucket) + "_" + field.name();
                unsigned numPerEntity = stk::mesh::field_scalars_per_entity(field, bucket);
                if(numPerEntity > 0)
                {
                    deviceData(bucket.bucket_id()) = BucketFieldDataType(name, bucket.size(), numPerEntity);
                    typename BucketFieldDataType::HostMirror tmpHost = Kokkos::create_mirror_view(deviceData(bucket.bucket_id()));

                    const T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                    {
                        for(unsigned j=0; j<numPerEntity; j++)
                        {
                            tmpHost(iEntity, j) = data[iEntity*numPerEntity + j];
                        }
                    }
                    Kokkos::deep_copy(deviceData(bucket.bucket_id()), tmpHost);
                    constDeviceData(iBucket) = deviceData(iBucket);
                }
            }
//            constDeviceData = deviceData;
        }
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        Kokkos::fence();

        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
//            Kokkos::deep_copy(hostData, deviceData);

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];
                unsigned numPerEntity = stk::mesh::field_scalars_per_entity(field, bucket);

                if(numPerEntity > 0)
                {
                    typename BucketFieldDataType::HostMirror tmpHost = Kokkos::create_mirror_view(deviceData(bucket.bucket_id()));
                    Kokkos::deep_copy(tmpHost, deviceData(bucket.bucket_id()));

                    T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                    {
                        for(unsigned j=0; j<numPerEntity; j++)
                        {
                            data[iEntity*numPerEntity + j] = tmpHost(iEntity, j);
                        }
                    }
                }
            }
        }
    }

    ~StaticField(){}

    template <typename Mesh> STK_FUNCTION
    T& get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(entity);
        return deviceData(index.bucket_id)(index.bucket_ord, component);
    }

    template <typename Mesh> STK_FUNCTION
    const T& const_get(const Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(entity);
        return constDeviceData(index.bucket_id)(index.bucket_ord, component);
    }

    STK_FUNCTION
    T& get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return deviceData(entity.bucket_id)(entity.bucket_ord, component);
    }

    STK_FUNCTION
    const T& const_get(stk::mesh::FastMeshIndex entity, int component) const
    {
        return constDeviceData(entity.bucket_id)(entity.bucket_ord, component);
    }

    template <typename MeshIndex> STK_FUNCTION
    T& get(MeshIndex entity, int component) const
    {
        return deviceData(entity.bucket->bucket_id())(entity.bucketOrd, component);
    }

    template <typename MeshIndex> STK_FUNCTION
    const T& const_get(MeshIndex entity, int component) const
    {
        return constDeviceData(entity.bucket->bucket_id())(entity.bucketOrd, component);
    }

private:

    STK_FUNCTION
    unsigned get_index(size_t bucketOrd) const
    {
        return bucketOrd;
    }

    void create_device_field_data(size_t numBuckets, size_t sizePerBucket, size_t dim, const T& initialValue)
    {
        deviceData = FieldDataType("deviceData", numBuckets);
        constDeviceData = ConstFieldDataType("constDeviceData", numBuckets);

        for(size_t i=0; i<numBuckets; ++i)
        {
            const std::string name = "bucket_" + std::to_string(i);
            deviceData(i) = BucketFieldDataType(name, sizePerBucket, dim);
            typename BucketFieldDataType::HostMirror tmpHost = Kokkos::create_mirror_view(deviceData(i));
            for(size_t j=0; j<sizePerBucket; ++j)
                for(size_t k=0; k<dim; ++k)
                    tmpHost(j,k) = initialValue;
            Kokkos::deep_copy(deviceData(i), tmpHost);
            constDeviceData(i) = deviceData(i);
        }
    }

    FieldDataType deviceData;
    ConstFieldDataType constDeviceData;
};


#ifdef KOKKOS_HAVE_CUDA
typedef ngp::StaticMesh StkNgpMesh;
typedef ngp::StaticField<double> StkNgpField;
#else
typedef ngp::WrapperMesh StkNgpMesh;
//typedef ngp::StaticField<double> StkNgpField;
typedef ngp::WrapperField<double> StkNgpField;
#endif

}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_ */
