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

#ifndef STK_NGP_NGPDYNAMICMESH_H_
#define STK_NGP_NGPDYNAMICMESH_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_MemoryPool.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include <string>

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Field.hpp"

#include <stk_ngp/NgpSpaces.hpp>
#include <stk_ngp/NgpMesh.hpp>
#include <stk_util/util/StkNgpVector.hpp>

#ifdef KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda Device;
#else
  typedef Kokkos::Serial Device;
#endif

namespace ngp {


typedef Kokkos::View<unsigned*,MemSpace> UnsignedViewType;
typedef Kokkos::View<stk::mesh::EntityRank*,MemSpace> RankViewType;
typedef Kokkos::View<unsigned*,Kokkos::LayoutRight, MemSpace> PartOrdViewType;
typedef Kokkos::View<unsigned*,Kokkos::LayoutRight, MemSpace, Kokkos::MemoryUnmanaged> UnmanagedPartOrdViewType;
typedef Kokkos::View<stk::mesh::Entity*, MemSpace> EntityViewType;
typedef Kokkos::View<stk::mesh::Entity*, MemSpace, Kokkos::MemoryUnmanaged> UnmanagedEntityViewType;
typedef Kokkos::View<stk::mesh::Entity**, Kokkos::LayoutRight, MemSpace> BktConnectivityType;
typedef Kokkos::View<stk::mesh::Entity**, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryUnmanaged> UnmanagedBktConnectivityType;

STK_FUNCTION
inline
bool all_parts_match(const UnmanagedPartOrdViewType& lhsParts, const UnmanagedPartOrdViewType& rhsParts)
{
    if (lhsParts.size() == 0 && rhsParts.size() == 0) {
        return true;
    }
    if (lhsParts.size() == 0 || rhsParts.size() == 0) {
        return false;
    }
    unsigned numParts = lhsParts(0);
    if (numParts != rhsParts(0)) {
        return false;
    }
    for(unsigned p=1; p<=numParts; ++p) {
        if (lhsParts(p) != rhsParts(p)) {
             return false;
        }
    }
    return true;
}

struct DynamicBucket {
    typedef Entities<const stk::mesh::Entity> ConnectedNodes;

    STK_FUNCTION
    DynamicBucket()
     : bucketId(0), entityRank(stk::topology::NODE_RANK), topo(stk::topology::INVALID_TOPOLOGY),
       partOrds(), entities(), connectivity()
    {
    }

    bool is_empty_on_host() const
    { return hostPartOrds.size() == 0 || hostPartOrds(0) == 0; }

    STK_FUNCTION
    bool is_empty() const { return partOrds.size() == 0 || partOrds(0) == 0; }
    
    void initialize(unsigned bucket_id_in, stk::mesh::EntityRank rank,
                    stk::topology stkTopo, unsigned nParts, unsigned nEntities,
                    unsigned nodesPerEntity)
    {
        bucketId = bucket_id_in;
        entityRank = rank;
        topo = stkTopo;
        numParts = nParts;
        numEntities = nEntities;
        numNodesPerEntity = nodesPerEntity;
    }

    void initialize_host_views_from_device_views()
    {
        hostPartOrds = Kokkos::create_mirror_view(partOrds);
        hostEntities = Kokkos::create_mirror_view(entities);
        hostConnectivity = Kokkos::create_mirror_view(connectivity);
    }

    void copy_to_device()
    {
        Kokkos::deep_copy(partOrds, hostPartOrds);
        Kokkos::deep_copy(entities, hostEntities);
        Kokkos::deep_copy(connectivity, hostConnectivity);
    }

    STK_FUNCTION
    unsigned bucket_id() const { return bucketId; }

    STK_FUNCTION
    size_t size() const { return numEntities; }

    STK_FUNCTION
    stk::mesh::EntityRank entity_rank() const { return entityRank; }

    STK_FUNCTION
    stk::topology topology() const { return topo; }

    unsigned host_get_num_parts() const { return hostPartOrds(0); }

    STK_FUNCTION
    unsigned get_num_parts() const { return partOrds(0); }

    STK_FUNCTION
    UnmanagedPartOrdViewType get_parts() const { return partOrds; }

    STK_FUNCTION
    bool is_member(unsigned partOrd) const
    {
        for(unsigned i=1; i<=partOrds(0); ++i)
        {
            if (partOrds(i) == partOrd)
            {
                return true;
            }
        }
        return false;
    }

    STK_FUNCTION
    unsigned get_num_nodes_per_entity() const { return topo.num_nodes(); }

    STK_FUNCTION
    ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
        return ConnectedNodes(&connectivity(offsetIntoBucket,0), connectivity.extent(1));
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
    stk::topology topo;
    unsigned numParts;
    unsigned numEntities;
    unsigned numNodesPerEntity;

    UnmanagedPartOrdViewType partOrds;
    PartOrdViewType::HostMirror hostPartOrds;

    UnmanagedEntityViewType entities;
    EntityViewType::HostMirror hostEntities;

    UnmanagedBktConnectivityType connectivity;
    BktConnectivityType::HostMirror hostConnectivity;
};

class DynamicMesh
{
public:
    typedef ngp::ExecSpace MeshExecSpace;
    typedef DynamicBucket::ConnectedNodes ConnectedNodes;
    typedef DynamicBucket BucketType;

    DynamicMesh() : bulk(nullptr) {}

    DynamicMesh(const stk::mesh::BulkData& b)
     : bulk(&b),
       maxNumParts(b.mesh_meta_data().get_parts().size()+1),
       bytesPerPartOrdAlloc(maxNumParts*sizeof(unsigned)),
       totalPartOrdBytes(100*bytesPerPartOrdAlloc),
       partOrdPool(Device::memory_space(), totalPartOrdBytes, bytesPerPartOrdAlloc, bytesPerPartOrdAlloc),

       maxNumEntities(512),
       bytesPerEntityViewAlloc(512*sizeof(stk::mesh::Entity)),
       totalEntityViewBytes(100*bytesPerEntityViewAlloc),
       entityViewPool(Device::memory_space(), totalEntityViewBytes, bytesPerEntityViewAlloc, bytesPerEntityViewAlloc),

       minNumConnectivity(maxNumEntities*2),
       maxNumConnectivity(maxNumEntities*8),
       minBytesPerConnectivityAlloc(minNumConnectivity*sizeof(stk::mesh::Entity)),
       maxBytesPerConnectivityAlloc(maxNumConnectivity*sizeof(stk::mesh::Entity)),
       totalConnectivityBytes(100*maxBytesPerConnectivityAlloc),
       connectivityViewPool(Device::memory_space(), totalConnectivityBytes, minBytesPerConnectivityAlloc, maxBytesPerConnectivityAlloc)
    {
        indexOfFirstEmptyBucket = UnsignedViewType("indexOfFirstEmptyBucket",stk::topology::NUM_RANKS);
        hostIndexOfFirstEmptyBucket = Kokkos::create_mirror_view(indexOfFirstEmptyBucket);
        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*,MemSpace>::HostMirror("host_mesh_indices", bulk->get_size_of_entity_index_space());
        hostEntityRanks = RankViewType::HostMirror("host_entity_ranks", bulk->get_size_of_entity_index_space());

        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());

        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            fill_buckets(*bulk, rank);
            fill_mesh_indices_and_ranks(*bulk, rank);
        }
        copy_mesh_indices_and_ranks_to_device();
    }

    ~DynamicMesh()
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
    stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
    {
        return deviceEntityRanks(entity.local_offset());
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

    stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
    {
        return ngp::get_bucket_ids(get_bulk_on_host(), rank, selector);
    }

    STK_FUNCTION
    unsigned num_buckets(stk::mesh::EntityRank rank) const
    {
        return indexOfFirstEmptyBucket(rank);
    }

    STK_FUNCTION
    const DynamicBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
    {
        return buckets[rank](index);
    }

    STK_FUNCTION
    int find_bucket_with_parts(stk::mesh::EntityRank rank, const UnmanagedPartOrdViewType& partOrds) const
    {
        unsigned numBuckets = num_buckets(rank);
        for(unsigned i=0; i<numBuckets; ++i) {
            if (buckets[rank](i).is_empty()) continue;
            if (all_parts_match(partOrds, buckets[rank](i).get_parts())) {
                return i;
            }
        }
        return -1;
    }

    STK_FUNCTION
    int add_bucket(stk::mesh::EntityRank rank, stk::topology topo, const UnmanagedPartOrdViewType& partOrds) const
    {
       int newBucketIndex = num_buckets(rank);
       buckets[rank](newBucketIndex).topo = topo;
       buckets[rank](newBucketIndex).partOrds = partOrds;
       buckets[rank](newBucketIndex).entities = UnmanagedEntityViewType(static_cast<stk::mesh::Entity*>(entityViewPool.allocate(bytesPerEntityViewAlloc)), maxNumEntities);
       unsigned nodesPerEntity = topo.num_nodes();
       unsigned bytesToAllocate = nodesPerEntity*maxNumEntities*sizeof(stk::mesh::Entity);
       buckets[rank](newBucketIndex).connectivity = UnmanagedBktConnectivityType(static_cast<UnmanagedBktConnectivityType::pointer_type>(connectivityViewPool.allocate(bytesToAllocate)), maxNumEntities, nodesPerEntity);
       indexOfFirstEmptyBucket(rank)++;
       return newBucketIndex;
    }

    STK_FUNCTION
    UnmanagedPartOrdViewType get_new_part_ords(unsigned partOrdToAdd, const UnmanagedPartOrdViewType& oldPartOrds) const
    {
         UnmanagedPartOrdViewType newParts = UnmanagedPartOrdViewType(static_cast<unsigned*>(partOrdPool.allocate(bytesPerPartOrdAlloc)), maxNumParts);
         newParts(0) = oldPartOrds.size()==0 ? 1 : oldPartOrds(0) + 1;
         unsigned offset = 1;
         bool alreadyAddedNewPart = false;
         for(unsigned i=1; i<=oldPartOrds(0); ++i) {
             if (oldPartOrds(i) > partOrdToAdd && !alreadyAddedNewPart) {
                 newParts(offset++) = partOrdToAdd;
                 alreadyAddedNewPart = true;
             }
             newParts(offset++) = oldPartOrds(i);
         }
         return newParts;
    }

    STK_FUNCTION
    void change_entity_parts(stk::mesh::Entity entity, unsigned partOrdToAdd) const
    {
        stk::mesh::FastMeshIndex meshIndex = device_mesh_index(entity);
        stk::mesh::EntityRank rank = entity_rank(entity);
        const DynamicBucket& bkt = get_bucket(rank, meshIndex.bucket_id);
        if (!bkt.is_member(partOrdToAdd))
        {
              UnmanagedPartOrdViewType newParts = get_new_part_ords(partOrdToAdd, bkt.get_parts());
              int bucketIndex = find_bucket_with_parts(rank, newParts);
              if (bucketIndex < 0) {
                  bucketIndex = add_bucket(rank, bkt.topology(), newParts);
              }
              stk::mesh::FastMeshIndex fmi = move_entity_to_bucket(rank, entity, bucketIndex);
              deviceMeshIndices(entity.local_offset()) = fmi;
        }
    }

    void clear()
    {
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
            buckets[rank] = BucketView();
    }

    struct alloc_part_ords_tag {};
    struct alloc_entities_tag {};
    struct alloc_connectivity_tag {};

    STK_FUNCTION
    void operator()(alloc_part_ords_tag, const int& bucketIndex) const {
        bucketsOfRank(bucketIndex).partOrds = UnmanagedPartOrdViewType(static_cast<unsigned*>(partOrdPool.allocate(bytesPerPartOrdAlloc)), maxNumParts);
    }

    STK_FUNCTION
    void operator()(alloc_entities_tag, const int& bucketIndex) const {
        bucketsOfRank(bucketIndex).entities = UnmanagedEntityViewType(static_cast<stk::mesh::Entity*>(entityViewPool.allocate(bytesPerEntityViewAlloc)), maxNumEntities);
    }

    STK_FUNCTION
    void operator()(alloc_connectivity_tag, const int& bucketIndex) const {
        unsigned nodesPerEntity = bucketsOfRank(bucketIndex).topology().num_nodes();
        unsigned numEntities = maxNumEntities;
        unsigned connAllocBytes = maxNumEntities*nodesPerEntity*sizeof(stk::mesh::Entity);
        bucketsOfRank(bucketIndex).connectivity = UnmanagedBktConnectivityType(static_cast<UnmanagedBktConnectivityType::pointer_type>(connectivityViewPool.allocate(connAllocBytes)), numEntities, nodesPerEntity);
    }

private:

    STK_FUNCTION
    stk::mesh::FastMeshIndex move_entity_to_bucket(stk::mesh::EntityRank rank,
                                                   stk::mesh::Entity entity,
                                                   unsigned newBucketIndex) const
    {
        stk::mesh::FastMeshIndex oldMeshIndex = device_mesh_index(entity);
        DynamicBucket& bucket = buckets[rank](oldMeshIndex.bucket_id);
        unsigned numEntities = bucket.size();
        unsigned numNodesPerEntity = bucket.get_num_nodes_per_entity();
        for(unsigned i=oldMeshIndex.bucket_ord+1; i<numEntities; ++i) {
            bucket.entities(i-1) = bucket.entities(i);
            for(unsigned j=0; j<numNodesPerEntity; ++j) {
                bucket.connectivity(i-1,j) = bucket.connectivity(i,j);
            }
        }
        bucket.numEntities--;
        DynamicBucket& newBucket = buckets[rank](newBucketIndex);
        unsigned newBucketOrd = newBucket.numEntities;
        newBucket.entities(newBucketOrd) = entity;
        newBucket.numEntities++;
        return stk::mesh::FastMeshIndex{newBucketIndex,newBucketOrd};
    }

    void fill_host_bucket_part_ords(stk::mesh::EntityRank rank, unsigned bucketId, const stk::mesh::Bucket& stkBucket)
    {
        hostBuckets[rank](bucketId).hostPartOrds(0) = hostBuckets[rank](bucketId).numParts;
        unsigned offset = 1;
        for(unsigned iPart = 0; iPart < stkBucket.supersets().size(); ++iPart)
        {
            const stk::mesh::Part* part = stkBucket.supersets()[iPart];
            hostBuckets[rank](bucketId).hostPartOrds(offset++) = part->mesh_meta_data_ordinal();
        }
    }

    void fill_host_bucket_entities(stk::mesh::EntityRank rank, unsigned bucketId, const stk::mesh::Bucket& stkBucket)
    {
        for(unsigned iEntity = 0; iEntity < stkBucket.size(); ++iEntity)
        {
            hostBuckets[rank](bucketId).hostEntities(iEntity) = stkBucket[iEntity];
        }
    }

    void fill_host_bucket_connectivity(stk::mesh::EntityRank rank, unsigned bucketId, const stk::mesh::Bucket& stkBucket)
    {
        for(unsigned iEntity = 0; iEntity < stkBucket.size(); ++iEntity)
        {
            const stk::mesh::Entity * elemNodes = stkBucket.begin_nodes(iEntity);
            unsigned nodesPerElem = stkBucket.topology().num_nodes();
            for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
            {
                hostBuckets[rank](bucketId).hostConnectivity(iEntity, iNode) = elemNodes[iNode];
            }
        }
    }

    void copy_buckets_to_device(stk::mesh::EntityRank rank)
    {
        Kokkos::deep_copy(buckets[rank], hostBuckets[rank]);
        for(unsigned b=0; b<hostBuckets[rank].size(); ++b)
        {
            if (!hostBuckets[rank](b).is_empty_on_host()) {
                hostBuckets[rank](b).copy_to_device();
            }
        }
    }

    void setup_bucket_storage(stk::mesh::EntityRank rank, unsigned numStkBuckets)
    {
        hostIndexOfFirstEmptyBucket(rank) = numStkBuckets;

        Kokkos::deep_copy(buckets[rank], hostBuckets[rank]);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace,alloc_part_ords_tag>(0,numStkBuckets), *this);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace,alloc_entities_tag>(0,numStkBuckets), *this);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace,alloc_connectivity_tag>(0,numStkBuckets), *this);

        Kokkos::deep_copy(indexOfFirstEmptyBucket, hostIndexOfFirstEmptyBucket);
        Kokkos::deep_copy(hostBuckets[rank], buckets[rank]);
    }

    void fill_buckets(const stk::mesh::BulkData& bulk_in, stk::mesh::EntityRank rank)
    {
        const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
        unsigned numStkBuckets = stkBuckets.size();

        unsigned emptyBucketsToAdd = std::max(5u, numStkBuckets/10);
        buckets[rank] = BucketView("Buckets", numStkBuckets+emptyBucketsToAdd);
        hostBuckets[rank] = Kokkos::create_mirror_view(buckets[rank]);

        bucketsOfRank = buckets[rank];

        for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket)
        {
            const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
            unsigned bucketId = stkBucket.bucket_id();
            unsigned nodesPerElem = stkBucket.topology().num_nodes();

            hostBuckets[rank](bucketId).initialize(bucketId, rank, stkBucket.topology(), stkBucket.supersets().size(), stkBucket.size(), nodesPerElem);
        }

        setup_bucket_storage(rank, numStkBuckets);

        for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket)
        {
            const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
            unsigned bucketId = stkBucket.bucket_id();

            hostBuckets[rank](bucketId).initialize_host_views_from_device_views();
            fill_host_bucket_part_ords(rank, bucketId, stkBucket);
            fill_host_bucket_entities(rank, bucketId, stkBucket);
            fill_host_bucket_connectivity(rank, bucketId, stkBucket);
        }

        for(unsigned iBucket=0; iBucket<emptyBucketsToAdd; ++iBucket)
        {
            unsigned bucketId = iBucket+numStkBuckets;
            hostBuckets[rank](iBucket+numStkBuckets).initialize(bucketId, rank, stk::topology::INVALID_TOPOLOGY, 0, 0, 0);
        }

        copy_buckets_to_device(rank);
    }

    void fill_mesh_indices_and_ranks(const stk::mesh::BulkData& bulk_in, stk::mesh::EntityRank rank)
    {
        const stk::mesh::BucketVector& bkts = bulk_in.buckets(rank);

        for(const stk::mesh::Bucket* bktptr : bkts)
        {
            const stk::mesh::Bucket& bkt = *bktptr;
            for(size_t i = 0; i < bkt.size(); ++i)
            {
                hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex{bkt.bucket_id(), static_cast<unsigned>(i)};
                hostEntityRanks[bkt[i].local_offset()] = rank;
            }
        }
    }

    void copy_mesh_indices_and_ranks_to_device()
    {
        unsigned length = hostMeshIndices.size();
        Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> tmp_device_mesh_indices("tmp_dev_mesh_indices", length);
        Kokkos::deep_copy(tmp_device_mesh_indices, hostMeshIndices);
        deviceMeshIndices = tmp_device_mesh_indices;
        RankViewType tmp_device_entity_ranks("tmp_dev_ent_ranks", length);
        Kokkos::deep_copy(tmp_device_entity_ranks, hostEntityRanks);
        deviceEntityRanks = tmp_device_entity_ranks;
    }

    const stk::mesh::BulkData &get_bulk_on_host() const
    {
        return *bulk;
    }


    typedef Kokkos::View<DynamicBucket*, Kokkos::LayoutRight, MemSpace> BucketView;
    const stk::mesh::BulkData *bulk;

    unsigned maxNumParts;
    unsigned bytesPerPartOrdAlloc;
    unsigned totalPartOrdBytes;
    Kokkos::MemoryPool<Device> partOrdPool;

    unsigned maxNumEntities;
    unsigned bytesPerEntityViewAlloc;
    unsigned totalEntityViewBytes;
    Kokkos::MemoryPool<Device> entityViewPool;

    unsigned minNumConnectivity;
    unsigned maxNumConnectivity;
    unsigned minBytesPerConnectivityAlloc;
    unsigned maxBytesPerConnectivityAlloc;
    unsigned totalConnectivityBytes;
    Kokkos::MemoryPool<Device> connectivityViewPool;

    BucketView buckets[stk::topology::NUM_RANKS];
    BucketView bucketsOfRank;
    UnsignedViewType indexOfFirstEmptyBucket;
    UnsignedViewType::HostMirror hostIndexOfFirstEmptyBucket;
    BucketView::HostMirror hostBuckets[stk::topology::NUM_RANKS];
    Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace>::HostMirror hostMeshIndices;
    Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> deviceMeshIndices;
    RankViewType deviceEntityRanks;
    RankViewType::HostMirror hostEntityRanks;
};


}


#endif /* STK_NGP_NGPMESH_H_ */
