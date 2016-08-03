#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPMESH_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPMESH_H_

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

#include <ngp/NgpSpaces.hpp>

namespace stk {
namespace mesh {
struct ConstMeshIndex
{
    const stk::mesh::Bucket *bucket;
    size_t bucketOrd;
};
}
}


namespace ngp {

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


class StkMeshAdapter
{
public:
    typedef ngp::HostExecSpace MeshExecSpace;

    typedef stk::mesh::ConstMeshIndex MeshIndex;
    typedef stk::mesh::Bucket BucketType;
    typedef Entities<const stk::mesh::Entity *> ConnectedNodes;

    StkMeshAdapter(const stk::mesh::BulkData& b) : bulk(b)
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




#ifdef KOKKOS_HAVE_CUDA
typedef Kokkos::LayoutLeft   Layout ;
#else
typedef Kokkos::LayoutRight   Layout ;
#endif

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
    typedef ngp::ExecSpace MeshExecSpace;

    typedef StaticMeshIndex MeshIndex;
    typedef StaticBucket BucketType;
    typedef Entities<ConnectedNodesType> ConnectedNodes;

    StaticMesh(const stk::mesh::BulkData& b) : bulk(b)
    {
        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror("host_mesh_indices", bulk.get_size_of_entity_index_space());
        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            fill_buckets(bulk, rank);
            fill_mesh_indices(bulk, rank);
        }
        copy_mesh_indices_to_device();
    }

    ~StaticMesh()
    {
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(const StaticMeshIndex &elem) const
    {
        ConnectedNodesType nodes = buckets[elem.bucket->entity_rank()](elem.bucket->bucket_id()).get_nodes(elem.bucketOrd);
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
        return buckets[rank].size();
    }

    STK_FUNCTION
    const StaticBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
    {
        return buckets[rank](index);
    }

    void clear()
    {
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
            buckets[rank] = BucketView();
    }


private:

    void fill_buckets(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank)
    {
        const stk::mesh::BucketVector& stkBuckets = bulk.buckets(rank);
        unsigned numStkBuckets = stkBuckets.size();

        buckets[rank] = BucketView("Buckets", numStkBuckets);
        BucketView &bucketsOfRank = buckets[rank];

        for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket)
        {
            const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
            unsigned bucketId = stkBucket.bucket_id();
            unsigned nodesPerElem = stkBucket.topology().num_nodes();
            bucketsOfRank(bucketId).initialize(stkBucket.bucket_id(), rank, stkBucket.size(), nodesPerElem);

            for(unsigned iEntity = 0; iEntity < stkBucket.size(); ++iEntity)
            {
                stk::mesh::Entity entity = stkBucket[iEntity];
                bucketsOfRank(bucketId).hostEntities(iEntity) = entity;
                const stk::mesh::Entity * elemNodes = bulk.begin_nodes(entity);
                for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
                {
                    bucketsOfRank(bucketId).hostConnectivity(iEntity, iNode) = elemNodes[iNode];
                }
            }
        }
        for(unsigned b=0; b<bucketsOfRank.size(); ++b) {
            bucketsOfRank(b).copy_to_device();
        }
    }

    void fill_mesh_indices(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank)
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
    BucketView buckets[stk::topology::NUM_RANKS];
    Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
    Kokkos::View<const stk::mesh::FastMeshIndex*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;
    //    Kokkos::View<const stk::mesh::FastMeshIndex*> device_mesh_indices;
};


}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPMESH_H_ */
