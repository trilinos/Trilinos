#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include <string>

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Field.hpp"

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

typedef Kokkos::View<stk::mesh::Entity*, Kokkos::LayoutStride, MemSpace> ConnectedNodesType;
typedef Kokkos::View<stk::mesh::Entity*, MemSpace> EntityViewType;
typedef Kokkos::View<stk::mesh::Entity**, Layout, MemSpace> BucketConnectivityType;

struct StaticBucket {
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

typedef Kokkos::View<StaticBucket*, UVMMemSpace> BucketsType;

struct StaticMeshIndex
{
    const ngp::StaticBucket &bucket;
    int bucketOrd;
    STK_FUNCTION StaticMeshIndex(const ngp::StaticBucket &bucketIn, int ordinal) : bucket(bucketIn), bucketOrd(ordinal) {}
};

typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;
typedef Kokkos::TeamPolicy<ExecSpace, ScheduleType> TeamPolicyType;
typedef TeamPolicyType::member_type TeamHandleType;
typedef Kokkos::TeamPolicy<HostExecSpace, ScheduleType> HostTeamPolicyType;
typedef HostTeamPolicyType::member_type HostTeamHandleType;

template <typename MESH, typename ALGORITHM_PER_ENTITY>
void for_each_entity_run(MESH &mesh, stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
{
    unsigned numBuckets = mesh.num_buckets(rank);
    Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace>(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const int bucketIndex = team.league_rank();
        const ngp::StaticBucket &bucket = mesh.get_bucket(rank, bucketIndex);
        unsigned numElements = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
        {
            functor(stk::mesh::FastMeshIndex{bucket.bucket_id(), static_cast<unsigned>(i)});
        });
    });
}

template <typename MESH, typename ALGORITHM_PER_ENTITY>
void for_each_entity_run_stkmesh(MESH &mesh, stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
{
    const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part());
    unsigned numBuckets = elemBuckets.size();
    Kokkos::parallel_for(Kokkos::TeamPolicy<HostExecSpace>(numBuckets, Kokkos::AUTO), [&](const HostTeamHandleType& team)
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

class StaticMesh
{
public:
    StaticMesh(const stk::mesh::BulkData& bulk)
    {
        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror("host_mesh_indices", bulk.get_size_of_entity_index_space());
        fill_mesh_indices(stk::topology::NODE_RANK, bulk);
        fill_mesh_indices(stk::topology::ELEM_RANK, bulk);
        fill_buckets(bulk, bulk.mesh_meta_data().locally_owned_part());
        copy_mesh_indices_to_device();
    }

    ~StaticMesh()
    {
    }

    STK_FUNCTION
    ConnectedNodesType get_nodes(const stk::mesh::FastMeshIndex &meshIndex) const
    {
        return buckets(meshIndex.bucket_id).get_nodes(meshIndex.bucket_ord);
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
        buckets = BucketsType();
    }

private:

    void fill_buckets(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);
        unsigned numElementBuckets = elemBuckets.size();

        buckets = BucketsType("ElementNodeConnectivity", numElementBuckets);

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
                hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex(bkt.bucket_id(), i);
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

    BucketsType buckets;
    Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
    Kokkos::View<const stk::mesh::FastMeshIndex*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;
    //    Kokkos::View<const stk::mesh::FastMeshIndex*> device_mesh_indices;
};

template<typename T>
class StaticField {
public:
    typedef Kokkos::View<T*> FieldDataType;

    StaticField(stk::mesh::EntityRank rank, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : deviceData(), numPerEntity(1)
    {
        unsigned alloc_size = compute_alloc_size(rank, bulk, selector);
        create_device_field_data(alloc_size, initialValue);
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : deviceData(), numPerEntity(0)
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
            numPerEntity = field.max_size(field.entity_rank());
            size_t allocationSize = buckets.size()*buckets[0]->capacity()*numPerEntity;
            deviceData = FieldDataType("deviceData"+std::to_string(numPerEntity), allocationSize);

            hostData = Kokkos::create_mirror_view(deviceData);

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];
                const T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                {
                    for(unsigned j=0; j<numPerEntity; j++)
                    {
                        unsigned index = iBucket*bucket.capacity()*numPerEntity + iEntity*numPerEntity + j;
                        hostData(index) = data[iEntity*numPerEntity + j];
                    }
                }
            }

            Kokkos::deep_copy(deviceData, hostData);
        }
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
            Kokkos::deep_copy(hostData, deviceData);

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];
                T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                {
                    for(unsigned j=0; j<numPerEntity; j++)
                    {
                        unsigned index = iBucket*bucket.capacity()*numPerEntity + iEntity*numPerEntity + j;
                        data[iEntity*numPerEntity + j] = hostData(index);
                    }
                }
            }
        }
    }

    ~StaticField(){}

    STK_FUNCTION
    unsigned get_index(const StaticMesh& ngpMesh, stk::mesh::Entity entity) const
    {
        unsigned bkt_id = ngpMesh.device_mesh_index(entity).bucket_id;
        unsigned bkt_ord = ngpMesh.device_mesh_index(entity).bucket_ord;
        unsigned idx = bkt_id*512*numPerEntity + bkt_ord*numPerEntity;
        return idx;
    }

    STK_FUNCTION
    T& get(const StaticMesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        unsigned idx = get_index(ngpMesh, entity);
        return deviceData(idx+component);
    }

private:
    unsigned compute_alloc_size(stk::mesh::EntityRank rank, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    {
        const stk::mesh::BucketVector& bkts = bulk.get_buckets(rank, selector);
        unsigned alloc_size = 0;
        for(const stk::mesh::Bucket* bktptr : bkts) {
            const stk::mesh::Bucket& bkt = *bktptr;
            alloc_size += bkt.capacity();
        }
        return alloc_size*numPerEntity;
    }

    void create_device_field_data(unsigned allocSize, const T& initialValue)
    {
        deviceData = FieldDataType("deviceData", allocSize);
        hostData = Kokkos::create_mirror_view(deviceData);

        for(size_t i=0; i<allocSize; ++i)
            hostData(i) = initialValue;
        Kokkos::deep_copy(deviceData, hostData);
    }

    typename FieldDataType::HostMirror hostData;
    FieldDataType deviceData;
    unsigned numPerEntity;
};

}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_STATICMESH_H_ */
