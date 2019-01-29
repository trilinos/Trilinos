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

#ifndef _STK_NGP_MESH_HPP_
#define _STK_NGP_MESH_HPP_

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

#include <stk_ngp/NgpSpaces.hpp>
#include <stk_util/util/StkNgpVector.hpp>

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
    Entities() : entities(nullptr), num(0), stride(bucketSize)
    {
    }
    STK_FUNCTION
    Entities(T* e, unsigned n, int stride_in=bucketSize) : entities(e), num(n), stride(stride_in)
    {
    }
    STK_FUNCTION
    T operator[](unsigned i) const
    {
#ifdef __CUDA_ARCH__
        return entities[stride*i];
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
    T* entities;
    unsigned num;
    int stride;
};

inline stk::NgpVector<unsigned> get_bucket_ids(const stk::mesh::BulkData &bulk,
                     stk::mesh::EntityRank rank,
                     const stk::mesh::Selector &selector)
{
    const stk::mesh::BucketVector &buckets = bulk.get_buckets(rank, selector);
    stk::NgpVector<unsigned> bucketIds(buckets.size());
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
    typedef Entities<const stk::mesh::Entity> ConnectedNodes;
    typedef Entities<const stk::mesh::Entity> ConnectedEntities;
    typedef Entities<const stk::mesh::ConnectivityOrdinal> ConnectedOrdinals;
    typedef Entities<const stk::mesh::Permutation> Permutations;

    StkMeshAdapter() : bulk(nullptr)
    {
    }
    StkMeshAdapter(const stk::mesh::BulkData& b) : bulk(&b)
    {
    }

    unsigned get_spatial_dimension() const
    {
        return bulk->mesh_meta_data().spatial_dimension();
    }

    ConnectedNodes get_nodes(const stk::mesh::ConstMeshIndex &elem) const
    {
        return ConnectedNodes(elem.bucket->begin_nodes(elem.bucketOrd), elem.bucket->num_nodes(elem.bucketOrd));
    }

    ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
        return ConnectedEntities(bucket.begin(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
    }

    ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
        return ConnectedOrdinals(bucket.begin_ordinals(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
    }

    ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::NODE_RANK);
    }

    ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::EDGE_RANK);
    }

    ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::FACE_RANK);
    }

    ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::ELEM_RANK);
    }

    ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
    }

    ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
    }

    ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
    }

    ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
    }

    Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
        return Permutations(bucket.begin_permutations(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
    }

    Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::NODE_RANK);
    }

    Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::EDGE_RANK);
    }

    Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::FACE_RANK);
    }

    Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::ELEM_RANK);
    }

    stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
    {
        const stk::mesh::MeshIndex &meshIndex = bulk->mesh_index(entity);
        return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), static_cast<unsigned>(meshIndex.bucket_ordinal)};
    }

    stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
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
typedef Kokkos::View<unsigned*, MemSpace> UnsignedViewType;
typedef Kokkos::View<stk::mesh::ConnectivityOrdinal*, MemSpace> OrdinalViewType;
typedef Kokkos::View<stk::mesh::Permutation*, MemSpace> PermutationViewType;

class StaticMesh;

struct StaticBucket {
    typedef Entities<const stk::mesh::Entity> ConnectedNodes;
    typedef Entities<const stk::mesh::Entity> ConnectedEntities;
    typedef Entities<const stk::mesh::ConnectivityOrdinal> ConnectedOrdinals;
    typedef Entities<const stk::mesh::Permutation> Permutations;

    STK_FUNCTION
    StaticBucket()
     : bucketId(0), entityRank(stk::topology::NODE_RANK), entities(),
       nodeConnectivity(), hostNodeConnectivity(),
       owningMesh(nullptr)
    {}

    void initialize(unsigned bucket_id_in, stk::mesh::EntityRank rank, stk::topology topo,
                    unsigned numEntities, unsigned numNodesPerEntity) {
        bucketId = bucket_id_in;
        entityRank = rank;
        bucketTopology = topo;
        const std::string bktIdStr(std::to_string(bucket_id_in));
        entities = EntityViewType(Kokkos::ViewAllocateWithoutInitializing("BucketEntities"+bktIdStr), numEntities);
        hostEntities = Kokkos::create_mirror_view(entities);
        nodeConnectivity = BucketConnectivityType(Kokkos::ViewAllocateWithoutInitializing("BucketConnectivity"+bktIdStr), bucketSize, numNodesPerEntity);
        hostNodeConnectivity = Kokkos::create_mirror_view(nodeConnectivity);
        nodeOrdinals = OrdinalViewType(Kokkos::ViewAllocateWithoutInitializing("NodeOrdinals"+bktIdStr), numNodesPerEntity);
        hostNodeOrdinals = Kokkos::create_mirror_view(nodeOrdinals);
        for(unsigned i=0; i<numNodesPerEntity; ++i)
        {
            hostNodeOrdinals(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
        }
    }

    void copy_to_device()
    {
        Kokkos::deep_copy(entities, hostEntities);
        Kokkos::deep_copy(nodeConnectivity, hostNodeConnectivity);
        Kokkos::deep_copy(nodeOrdinals, hostNodeOrdinals);
    }

    STK_FUNCTION
    unsigned bucket_id() const { return bucketId; }

    STK_FUNCTION
    size_t size() const { return entities.size(); }

    STK_FUNCTION
    stk::mesh::EntityRank entity_rank() const { return entityRank; }

    STK_FUNCTION
    stk::topology topology() const { return bucketTopology; }

    STK_FUNCTION
    unsigned get_num_nodes_per_entity() const { return nodeConnectivity.extent(1); }

    STK_INLINE_FUNCTION
    ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

    STK_INLINE_FUNCTION
    ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

    STK_FUNCTION
    ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
        return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
    }

    STK_FUNCTION
    ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
        return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
    }

    STK_FUNCTION
    ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
        return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
    }

    STK_FUNCTION
    ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
        return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
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
    stk::topology bucketTopology;

    EntityViewType entities;
    EntityViewType::HostMirror hostEntities;

    BucketConnectivityType nodeConnectivity;
    BucketConnectivityType::HostMirror hostNodeConnectivity;

    OrdinalViewType nodeOrdinals;
    OrdinalViewType::HostMirror hostNodeOrdinals;

    const ngp::StaticMesh* owningMesh;
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
    typedef StaticBucket::ConnectedEntities ConnectedEntities;
    typedef StaticBucket::ConnectedOrdinals ConnectedOrdinals;
    typedef StaticBucket::Permutations Permutations;
    typedef StaticMeshIndex MeshIndex;
    typedef StaticBucket BucketType;

    STK_FUNCTION
    StaticMesh() : bulk(nullptr), spatial_dimension(0) {}

    explicit StaticMesh(const stk::mesh::BulkData& b)
     : bulk(&b), spatial_dimension(b.mesh_meta_data().spatial_dimension())
    {
        set_bucket_entity_offsets(b);
        copy_bucket_entity_offsets_to_device();
        fill_sparse_connectivities(b);
        copy_sparse_connectivities_to_device();

        hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror(Kokkos::ViewAllocateWithoutInitializing("host_mesh_indices"), bulk->get_size_of_entity_index_space());
        const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            fill_buckets(*bulk, rank);
            fill_mesh_indices(*bulk, rank);
        }
        copy_mesh_indices_to_device();
    }

    STK_FUNCTION
    StaticMesh(const StaticMesh &) = default;
    STK_FUNCTION
    ~StaticMesh() {}

    STK_FUNCTION
    unsigned get_spatial_dimension() const
    {
        return spatial_dimension;
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(const StaticMeshIndex &entity) const
    {
        return buckets[entity.bucket->entity_rank()](entity.bucket->bucket_id()).get_nodes(entity.bucketOrd);
    }

    STK_FUNCTION
    ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        if (connectedRank == stk::topology::NODE_RANK)
        {
            return buckets[rank](entity.bucket_id).get_connected_entities(entity.bucket_ord, connectedRank);
        }

        int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
        int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
        size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
                            - connectivityOffset;
        ConnectedEntities connectedEntities(nullptr, 0);
        if (numConnected > 0)
        {
            int stride = 1;
            connectedEntities = ConnectedEntities(&sparseConnectivity[rank][connectedRank](connectivityOffset), numConnected, stride);
        }
        return connectedEntities;
    }

    STK_FUNCTION
    ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return buckets[rank](entity.bucket_id).get_nodes(entity.bucket_ord);
    }

    STK_FUNCTION
    ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::EDGE_RANK);
    }

    STK_FUNCTION
    ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::FACE_RANK);
    }

    STK_FUNCTION
    ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_entities(rank, entity, stk::topology::ELEM_RANK);
    }

    STK_FUNCTION
    ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        if (connectedRank == stk::topology::NODE_RANK)
        {
            return buckets[rank](entity.bucket_id).get_connected_ordinals(entity.bucket_ord, connectedRank);
        }

        int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
        int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
        size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
                            - connectivityOffset;
        ConnectedOrdinals connectedOrdinals(nullptr, 0);
        if (numConnected > 0)
        {
            int stride = 1;
            connectedOrdinals = ConnectedOrdinals(&sparseConnectivityOrdinals[rank][connectedRank](connectivityOffset), numConnected, stride);
        }
        return connectedOrdinals;
    }

    STK_FUNCTION
    ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
    }

    STK_FUNCTION
    ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
    }

    STK_FUNCTION
    ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
    }

    STK_FUNCTION
    ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
    }

    STK_FUNCTION
    Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
    {
        Permutations permutations(nullptr, 0);
        if (connectedRank == stk::topology::NODE_RANK)
        {
            return permutations;
        }

        int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
        int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
        size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
                            - connectivityOffset;
        if (numConnected > 0)
        {
            int stride = 1;
            permutations = Permutations(&sparsePermutations[rank][connectedRank](connectivityOffset), numConnected, stride);
        }
        return permutations;
    }

    STK_FUNCTION
    Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::NODE_RANK);
    }

    STK_FUNCTION
    Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::EDGE_RANK);
    }

    STK_FUNCTION
    Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::FACE_RANK);
    }

    STK_FUNCTION
    Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
    {
        return get_permutations(rank, entity, stk::topology::ELEM_RANK);
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
        return buckets[rank].size();
    }

    STK_FUNCTION
    const StaticBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
    {
        buckets[rank](index).owningMesh = this;
        return buckets[rank](index);
    }

    void clear()
    {
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
            buckets[rank] = BucketView();
    }


private:

    void set_bucket_entity_offsets(const stk::mesh::BulkData& bulk_in)
    {
        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
            bucketEntityOffsets[rank] = Kokkos::View<int*,MemSpace>(Kokkos::ViewAllocateWithoutInitializing("bucketEntityOffsets"), stkBuckets.size()+1);
            hostBucketEntityOffsets[rank] = Kokkos::create_mirror_view(bucketEntityOffsets[rank]);

            int bucketOffsetIntoEntities = 0;
            for(unsigned i=0; i<stkBuckets.size(); ++i)
            {
                hostBucketEntityOffsets[rank](i) = bucketOffsetIntoEntities;
                bucketOffsetIntoEntities += stkBuckets[i]->size();
            }
            hostBucketEntityOffsets[rank](stkBuckets.size()) = bucketOffsetIntoEntities;
        }
    }

    void fill_sparse_connectivities(const stk::mesh::BulkData& bulk_in)
    {
        unsigned totalNumConnectedEntities[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}};
        unsigned totalNumPermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}};
        const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
            for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
            {
                const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
                const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
                for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
                {
                    const bool hasPermutation = stkBucket.has_permutation(connectedRank);
                    for(unsigned iEntity=0; iEntity<stkBucket.size(); ++iEntity)
                    {
                        const unsigned numConnected = stkBucket.num_connectivity(iEntity, connectedRank);
                        totalNumConnectedEntities[rank][connectedRank] += numConnected;
                        totalNumPermutations[rank][connectedRank] += hasPermutation ? numConnected : 0;
                    }
                }
            }

            const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
            for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
            {
                size_t numEntities = hostBucketEntityOffsets[rank](hostBucketEntityOffsets[rank].size()-1);
                entityConnectivityOffset[rank][connectedRank] = UnsignedViewType(Kokkos::ViewAllocateWithoutInitializing("FacesRowMap"), numEntities+1);
                hostEntityConnectivityOffset[rank][connectedRank] = Kokkos::create_mirror_view(entityConnectivityOffset[rank][connectedRank]);

                sparseConnectivity[rank][connectedRank] = EntityViewType(Kokkos::ViewAllocateWithoutInitializing("SparseConnectivity"), totalNumConnectedEntities[rank][connectedRank]);
                hostSparseConnectivity[rank][connectedRank] = Kokkos::create_mirror_view(sparseConnectivity[rank][connectedRank]);

                sparseConnectivityOrdinals[rank][connectedRank] = OrdinalViewType(Kokkos::ViewAllocateWithoutInitializing("SparseConnectivityOrdinals"), totalNumConnectedEntities[rank][connectedRank]);
                hostSparseConnectivityOrdinals[rank][connectedRank] = Kokkos::create_mirror_view(sparseConnectivityOrdinals[rank][connectedRank]);

                sparsePermutations[rank][connectedRank] = PermutationViewType(Kokkos::ViewAllocateWithoutInitializing("SparsePermutations"), totalNumPermutations[rank][connectedRank]);
                hostSparsePermutations[rank][connectedRank] = Kokkos::create_mirror_view(sparsePermutations[rank][connectedRank]);
            }
            int entriesOffsets[stk::topology::NUM_RANKS] = {0};
            int myOffset = 0;
            for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
            {
                const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
                int bucketEntityOffset = hostBucketEntityOffsets[rank](stkBucket.bucket_id());
                for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
                {
                    const bool hasPermutation = stkBucket.has_permutation(connectedRank);
                    const stk::mesh::Entity* connectedEntities = stkBucket.begin(0, connectedRank);
                    const stk::mesh::ConnectivityOrdinal* connectedOrdinals = stkBucket.begin_ordinals(0, connectedRank);
                    const stk::mesh::Permutation* permutations = hasPermutation ? stkBucket.begin_permutations(0, connectedRank) : nullptr;
                    for(unsigned iEntity=0; iEntity<stkBucket.size(); ++iEntity)
                    {
                        myOffset = bucketEntityOffset + iEntity;
                        unsigned numConnected = stkBucket.num_connectivity(iEntity, connectedRank);

                        int entriesOffset = entriesOffsets[connectedRank];
                        hostEntityConnectivityOffset[rank][connectedRank](myOffset) = entriesOffset;
                        for(unsigned i=0; i<numConnected; ++i)
                        {
                            hostSparseConnectivity[rank][connectedRank](entriesOffset+i) = connectedEntities[i];
                            hostSparseConnectivityOrdinals[rank][connectedRank](entriesOffset+i) = connectedOrdinals[i];
                            if (hasPermutation) {
                                hostSparsePermutations[rank][connectedRank](entriesOffset+i) = permutations[i];
                            }
                        }
                        entriesOffsets[connectedRank] = entriesOffset + numConnected;
                        connectedEntities += numConnected;
                        connectedOrdinals += numConnected;
                        if (hasPermutation) permutations += numConnected;
                    }
                }
            }
            for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
            {
                hostEntityConnectivityOffset[rank][connectedRank](myOffset+1) = entriesOffsets[connectedRank];
            }
        }
    }

    void fill_buckets(const stk::mesh::BulkData& bulk_in, stk::mesh::EntityRank bucketRank)
    {
        const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(bucketRank);
        unsigned numStkBuckets = stkBuckets.size();

        buckets[bucketRank] = BucketView("Buckets", numStkBuckets);
        Kokkos::fence();
        BucketView &bucketsOfRank = buckets[bucketRank];

        for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket)
        {
            const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
            unsigned bucketId = stkBucket.bucket_id();
            unsigned nodesPerElem = stkBucket.topology().num_nodes();
            BucketType& bucket = bucketsOfRank(bucketId);
            bucket.initialize(stkBucket.bucket_id(), bucketRank, stkBucket.topology(), stkBucket.size(), nodesPerElem);

            for(unsigned iEntity = 0; iEntity < stkBucket.size(); ++iEntity)
            {
                stk::mesh::Entity entity = stkBucket[iEntity];
                bucket.hostEntities(iEntity) = entity;
                const stk::mesh::Entity * elemNodes = stkBucket.begin_nodes(iEntity);
                for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
                {
                    bucket.hostNodeConnectivity(iEntity, iNode) = elemNodes[iNode];
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
            const unsigned bktId = bkt.bucket_id();
            for(unsigned i = 0; i < bkt.size(); ++i)
            {
                hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex{bktId, i};
            }
        }
    }

    void copy_mesh_indices_to_device()
    {
        unsigned length = hostMeshIndices.size();
        Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> nonconst_device_mesh_indices(Kokkos::ViewAllocateWithoutInitializing("tmp_dev_mesh_indices"), length);
        Kokkos::deep_copy(nonconst_device_mesh_indices, hostMeshIndices);
        deviceMeshIndices = nonconst_device_mesh_indices;
    }

    void copy_bucket_entity_offsets_to_device()
    {
        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::EDGE_RANK; rank<endRank; rank++)
        {
            Kokkos::deep_copy(bucketEntityOffsets[rank], hostBucketEntityOffsets[rank]);
        }
    }

    void copy_sparse_connectivities_to_device()
    {
        const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
        {
            const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
            for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
            {
                Kokkos::deep_copy(entityConnectivityOffset[rank][connectedRank], hostEntityConnectivityOffset[rank][connectedRank]);
                Kokkos::deep_copy(sparseConnectivity[rank][connectedRank], hostSparseConnectivity[rank][connectedRank]);
                Kokkos::deep_copy(sparseConnectivityOrdinals[rank][connectedRank], hostSparseConnectivityOrdinals[rank][connectedRank]);
                Kokkos::deep_copy(sparsePermutations[rank][connectedRank], hostSparsePermutations[rank][connectedRank]);
            }
        }
    }

    const stk::mesh::BulkData &get_bulk_on_host() const
    {
        return *bulk;
    }


    typedef Kokkos::View<StaticBucket*, UVMMemSpace> BucketView;
    const stk::mesh::BulkData *bulk;
    unsigned spatial_dimension;
    BucketView buckets[stk::topology::NUM_RANKS];
    Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
    Kokkos::View<const stk::mesh::FastMeshIndex*, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;

    Kokkos::View<int*,MemSpace> bucketEntityOffsets[stk::topology::NUM_RANKS];
    Kokkos::View<int*,MemSpace>::HostMirror hostBucketEntityOffsets[stk::topology::NUM_RANKS];

    UnsignedViewType entityConnectivityOffset[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
    UnsignedViewType::HostMirror hostEntityConnectivityOffset[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

    EntityViewType sparseConnectivity[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
    EntityViewType::HostMirror hostSparseConnectivity[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

    OrdinalViewType sparseConnectivityOrdinals[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
    OrdinalViewType::HostMirror hostSparseConnectivityOrdinals[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

    PermutationViewType sparsePermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
    PermutationViewType::HostMirror hostSparsePermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
};

STK_INLINE_FUNCTION
StaticBucket::ConnectedEntities
StaticBucket::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
    NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
    if (connectedRank == stk::topology::NODE_RANK) {
        return ConnectedEntities(&nodeConnectivity(offsetIntoBucket,0), nodeConnectivity.extent(1));
    }
    NGP_ThrowAssert(owningMesh != nullptr);
    stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
    return owningMesh->get_connected_entities(entity_rank(), meshIndex, connectedRank);
}

STK_INLINE_FUNCTION
StaticBucket::ConnectedOrdinals
StaticBucket::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
    NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
    if (connectedRank == stk::topology::NODE_RANK) {
        return ConnectedOrdinals(&nodeOrdinals(0), nodeOrdinals.extent(0));
    }
    NGP_ThrowAssert(owningMesh != nullptr);
    stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
    return owningMesh->get_connected_ordinals(entity_rank(), meshIndex, connectedRank);
}

}


#endif /* STK_NGP_MESH_HPP_ */
