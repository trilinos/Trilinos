// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

constexpr int bucketSize = 512;

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
#ifdef KOKKOS_ENABLE_CUDA
        return entities[bucketSize*i];
#else
        return entities[i];
#endif
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

    StkMeshAdapter() : bulk(nullptr)
    {
    }
    StkMeshAdapter(const stk::mesh::BulkData& b) : bulk(&b)
    {
    }

    ConnectedNodes get_nodes(const stk::mesh::ConstMeshIndex &elem) const
    {
        return ConnectedNodes(elem.bucket->begin_nodes(elem.bucketOrd), elem.bucket->num_nodes(elem.bucketOrd));
    }

    ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &elem) const
    {
        const stk::mesh::Bucket &bucket = get_bucket(rank, elem.bucket_id);
        return ConnectedNodes(bucket.begin_nodes(elem.bucket_ord), bucket.num_nodes(elem.bucket_ord));
    }

    stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
    {
        const stk::mesh::MeshIndex &meshIndex = bulk->mesh_index(entity);
        return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), static_cast<unsigned>(meshIndex.bucket_ordinal)};
    }

    stk::Vector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
    {
        return ngp::get_bucket_ids(*bulk, rank, selector);
    }

    unsigned num_buckets(stk::mesh::EntityRank rank) const
    {
        return bulk->buckets(rank).size();
    }

    const stk::mesh::Bucket & get_bucket(stk::mesh::EntityRank rank, unsigned i) const
    {
        return *bulk->buckets(rank)[i];
    }

private:
    const stk::mesh::BulkData *bulk;
};




typedef Kokkos::View<stk::mesh::Entity*, MemSpace> EntityViewType;
typedef Kokkos::View<stk::mesh::Entity**, MemSpace> BucketConnectivityType;

struct StaticBucket {
    typedef Entities<const stk::mesh::Entity *> ConnectedNodes;

    STK_FUNCTION
    StaticBucket()
     : bucketId(0), entityRank(stk::topology::NODE_RANK), entities(), connectivity() {}

    void initialize(unsigned bucket_id_in, stk::mesh::EntityRank rank, unsigned numEntities, unsigned numNodesPerEntity) {
        bucketId = bucket_id_in;
        entityRank = rank;
        entities = EntityViewType("BucketEntities"+std::to_string(bucket_id_in), numEntities);
        hostEntities = Kokkos::create_mirror_view(entities);
        connectivity = BucketConnectivityType("BucketConnectivity"+std::to_string(bucket_id_in), bucketSize, numNodesPerEntity);
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
    ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
        return ConnectedNodes(&connectivity(offsetIntoBucket,0), connectivity.dimension_1());
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
    typedef StaticBucket::ConnectedNodes ConnectedNodes;
    typedef StaticMeshIndex MeshIndex;
    typedef StaticBucket BucketType;

    StaticMesh() : bulk(nullptr) {}

    StaticMesh(const stk::mesh::BulkData& b) : bulk(&b)
    {
        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror("host_mesh_indices", bulk->get_size_of_entity_index_space());
        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            fill_buckets(*bulk, rank);
            fill_mesh_indices(*bulk, rank);
        }
        copy_mesh_indices_to_device();
    }

    ~StaticMesh()
    {
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(const StaticMeshIndex &elem) const
    {
        return buckets[elem.bucket->entity_rank()](elem.bucket->bucket_id()).get_nodes(elem.bucketOrd);
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &elem) const
    {
        return buckets[rank](elem.bucket_id).get_nodes(elem.bucket_ord);
    }

    STK_FUNCTION
    stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
    {
        return device_mesh_index(entity);
    }

    STK_FUNCTION
    stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
    {
        return deviceMeshIndices(entity.local_offset());
    }

    const stk::mesh::FastMeshIndex& host_mesh_index(stk::mesh::Entity entity) const
    {
        return hostMeshIndices(entity.local_offset());
    }

    stk::Vector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
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

    void fill_buckets(const stk::mesh::BulkData& bulk_in, stk::mesh::EntityRank rank)
    {
        const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
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
                const stk::mesh::Entity * elemNodes = bulk_in.begin_nodes(entity);
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

    void fill_mesh_indices(const stk::mesh::BulkData& bulk_in, stk::mesh::EntityRank rank)
    {
        const stk::mesh::BucketVector& bkts = bulk_in.buckets(rank);

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

    const stk::mesh::BulkData &get_bulk_on_host() const
    {
        return *bulk;
    }


    typedef Kokkos::View<StaticBucket*, UVMMemSpace> BucketView;
    const stk::mesh::BulkData *bulk;
    BucketView buckets[stk::topology::NUM_RANKS];
    Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
    Kokkos::View<const stk::mesh::FastMeshIndex*, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;
    //    Kokkos::View<const stk::mesh::FastMeshIndex*> device_mesh_indices;
};


}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPMESH_H_ */
