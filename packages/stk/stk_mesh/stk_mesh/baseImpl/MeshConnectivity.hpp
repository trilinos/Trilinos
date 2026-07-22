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

#ifndef STK_MESH_IMPL_MESHCONNECTIVITY_HPP
#define STK_MESH_IMPL_MESHCONNECTIVITY_HPP

#include <Kokkos_Core.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/StridedArray.hpp>
#include <sstream>
#include <iostream>

namespace stk {
namespace mesh {
namespace impl {

using ConnectedEntities = stk::util::StridedArray<const stk::mesh::Entity>;
using ConnectedOrdinals = stk::util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
using ConnectedPermutations = stk::util::StridedArray<const stk::mesh::Permutation>;

class EntityConnectivity
{
public:

  KOKKOS_DEFAULTED_FUNCTION EntityConnectivity() = default;
  KOKKOS_DEFAULTED_FUNCTION EntityConnectivity(const EntityConnectivity&) = default;
  KOKKOS_DEFAULTED_FUNCTION EntityConnectivity(EntityConnectivity&&) = default;
  KOKKOS_DEFAULTED_FUNCTION EntityConnectivity& operator=(const EntityConnectivity&) = default;
  KOKKOS_DEFAULTED_FUNCTION EntityConnectivity& operator=(EntityConnectivity&&) = default;

  KOKKOS_FUNCTION
  ~EntityConnectivity(){}

  KOKKOS_INLINE_FUNCTION
  uint16_t get_capacity() const { return m_capacity; }

  KOKKOS_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank connectedRank) const
  {
    const uint16_t len = m_offsets[connectedRank+1]-m_offsets[connectedRank];
    return len==0 ? ConnectedEntities(nullptr, 0) : ConnectedEntities(&m_conn[m_offsets[connectedRank]], len);
  }

  KOKKOS_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank connectedRank) const
  {
    const uint16_t len = m_offsets[connectedRank+1]-m_offsets[connectedRank];
    return len==0 ? ConnectedOrdinals(nullptr, 0) : ConnectedOrdinals(&m_ord[m_offsets[connectedRank]], len);
  }

  KOKKOS_INLINE_FUNCTION
  ConnectedPermutations get_connected_permutations(stk::mesh::EntityRank connectedRank) const
  {
    const uint16_t len = m_offsets[connectedRank+1]-m_offsets[connectedRank];
    return (m_perm==nullptr || len==0) ? ConnectedPermutations(nullptr, 0) : ConnectedPermutations(&m_perm[m_offsets[connectedRank]], len); 
  }

  template<typename PoolAllocator>
  void add_connectivity(PoolAllocator& poolAlloc,
                        stk::mesh::EntityRank connRank,
                        unsigned numConnectivity, 
                        const stk::mesh::Entity* entities,
                        const stk::mesh::ConnectivityOrdinal* ords,
                        const stk::mesh::Permutation* perms,
                        bool addPermutations)
  {
    if (numConnectivity == 0) {
      return;
    }
    uint16_t curSize = m_offsets[stk::topology::NUM_RANKS];
    uint16_t avail = m_capacity - curSize;
    if (avail < numConnectivity) {
      grow_capacity(poolAlloc, numConnectivity, addPermutations);
    }     

    uint16_t destOffset = m_offsets[connRank+1];

    STK_ThrowRequireMsg(m_conn != nullptr,"m_conn is nullptr");
    STK_ThrowRequireMsg(m_ord != nullptr,"m_ord is nullptr");
    if (addPermutations && m_perm == nullptr) {
      add_permutations(poolAlloc);
    }

    stk::mesh::EntityRank lastRank = static_cast<stk::mesh::EntityRank>(stk::topology::NUM_RANKS-1);
    for(stk::mesh::EntityRank rank = lastRank; rank > connRank; --rank) {
      uint16_t offset = m_offsets[rank];
      uint16_t num = m_offsets[rank+1] - offset;
      for(uint16_t i=0; i<num; ++i) {
        m_conn[offset+numConnectivity+i] = m_conn[offset+i];
        m_ord[offset+numConnectivity+i] = m_ord[offset+i];
        if (addPermutations) {
          m_perm[offset+numConnectivity+i] = m_perm[offset+i];
        }
      }

      m_offsets[rank+1] += numConnectivity;
    }

    m_offsets[connRank+1] += numConnectivity;

    for(uint16_t i=0; i<numConnectivity; ++i) {
      m_conn[destOffset+i] = entities[i];
      m_ord[destOffset+i] = ords[i];
      if (addPermutations) {
        m_perm[destOffset+i] = perms != nullptr ? perms[i] : INVALID_PERMUTATION;
      }
    }
  }

  template<typename PoolAllocator>
  void reset(PoolAllocator& poolAlloc, unsigned newCapacity, bool addPermutations)
  {
    if (newCapacity == 0) {
      const size_t unusedAllocSize = 0;
      poolAlloc.deallocate(m_data, unusedAllocSize);
      m_conn = nullptr;
      m_ord = nullptr;
      m_perm = nullptr;
      m_data = nullptr;
      m_capacity = 0;
    }

    if (newCapacity > m_capacity) {
      uint16_t sizeofPerm = addPermutations ? sizeof(stk::mesh::Permutation) : 0;
      size_t newNumBytes = newCapacity*(sizeof(stk::mesh::Entity) + sizeof(stk::mesh::ConnectivityOrdinal) + sizeofPerm);
      void* newData = poolAlloc.allocate(newNumBytes);
      throw_err_if_nullptr(newData, newNumBytes, poolAlloc);

      stk::mesh::Entity* newConn = (stk::mesh::Entity*)newData;
      stk::mesh::Entity* endConn = newConn+newCapacity;
      stk::mesh::ConnectivityOrdinal* newOrds = (stk::mesh::ConnectivityOrdinal*)endConn;
      stk::mesh::ConnectivityOrdinal* endOrds = newOrds+newCapacity;
      stk::mesh::Permutation* newPerms = addPermutations ? (stk::mesh::Permutation*)endOrds : nullptr;

      const size_t unusedAllocSize = 0;
      poolAlloc.deallocate(m_data, unusedAllocSize);

      m_conn = newConn;
      m_ord = newOrds;
      m_perm = newPerms;
      m_data = newData;
      m_capacity = newCapacity;
    }

    for(uint16_t i=0; i<static_cast<uint16_t>(stk::topology::NUM_RANKS+1); ++i) {
      m_offsets[i] = 0;
    }
  }

  template<typename PoolAllocator>
  void grow_capacity(PoolAllocator& poolAlloc,
                     unsigned numConnectivity,
                     bool addPermutations)
  {
    uint16_t newCapacity = m_capacity + numConnectivity;
    if (newCapacity == 0) {
      return;
    }
    uint16_t sizeofPerm = addPermutations ? sizeof(stk::mesh::Permutation) : 0;
    size_t newNumBytes = newCapacity*(sizeof(stk::mesh::Entity) + sizeof(stk::mesh::ConnectivityOrdinal) + sizeofPerm);

    if(newNumBytes > poolAlloc.max_block_size()) {
      throw_err_if_nullptr(nullptr, newNumBytes, poolAlloc);
    }

    void* newData = poolAlloc.allocate(newNumBytes);
    throw_err_if_nullptr(newData, newNumBytes, poolAlloc);

    stk::mesh::Entity* newConn = (stk::mesh::Entity*)newData;
    stk::mesh::Entity* endConn = newConn+newCapacity;
    stk::mesh::ConnectivityOrdinal* newOrds = (stk::mesh::ConnectivityOrdinal*)endConn;
    stk::mesh::ConnectivityOrdinal* endOrds = newOrds+newCapacity;
    stk::mesh::Permutation* newPerms = addPermutations ? (stk::mesh::Permutation*)endOrds : nullptr;

    if (m_data != nullptr) {
      uint16_t curSize = m_offsets[stk::topology::NUM_RANKS];
      for(uint16_t i=0; i<curSize; ++i) {
        newConn[i] = m_conn[i];
        newOrds[i] = m_ord[i];
        if (addPermutations) {
          newPerms[i] = m_perm[i];
        }
      }

      const size_t unusedAllocSize = 0;
      poolAlloc.deallocate(m_data, unusedAllocSize);
    }

    m_conn = newConn;
    m_ord = newOrds;
    m_perm = newPerms;
    m_data = newData;
    m_capacity = newCapacity;
  }

  template<typename PoolAllocator>
  void add_permutations(PoolAllocator& poolAlloc)
  {
    const bool alreadyHavePermutations = m_perm != nullptr;
    if (alreadyHavePermutations) {
      return;
    }

    uint16_t sizeofPerm = sizeof(stk::mesh::Permutation);
    size_t newNumBytes = m_capacity*(sizeof(stk::mesh::Entity) + sizeof(stk::mesh::ConnectivityOrdinal) + sizeofPerm);
    void* newData = poolAlloc.allocate(newNumBytes);
    throw_err_if_nullptr(newData, newNumBytes, poolAlloc);

    stk::mesh::Entity* newConn = (stk::mesh::Entity*)newData;
    stk::mesh::Entity* endConn = newConn+m_capacity;
    stk::mesh::ConnectivityOrdinal* newOrds = (stk::mesh::ConnectivityOrdinal*)endConn;
    stk::mesh::ConnectivityOrdinal* endOrds = newOrds+m_capacity;
    stk::mesh::Permutation* newPerms = (stk::mesh::Permutation*)endOrds;

    uint16_t curSize = m_offsets[stk::topology::NUM_RANKS];
    for(uint16_t i=0; i<curSize; ++i) {
      newConn[i] = m_conn[i];
      newOrds[i] = m_ord[i];
    }

    const size_t unusedAllocSize = 0;
    poolAlloc.deallocate(m_data, unusedAllocSize);

    m_conn = newConn;
    m_ord = newOrds;
    m_perm = newPerms;
    m_data = newData;
    for(uint16_t i=0; i<m_capacity; ++i) {
      m_perm[i] = INVALID_PERMUTATION;
    }
  }

  bool is_empty() const { return m_offsets[stk::topology::NUM_RANKS] == 0; }

private:
  template<typename PoolAllocator>
  void throw_err_if_nullptr(void* ptr, size_t numBytes, PoolAllocator& poolAlloc)
  {
    if (ptr == nullptr && numBytes > 0) {
      std::ostringstream oss;
      oss<<"poolAlloc returned nullptr, failed to allocate "<<numBytes<<" bytes\n";
      poolAlloc.print_state(oss);
      typename PoolAllocator::usage_statistics usageStats;
      poolAlloc.get_usage_statistics(usageStats);
      oss<<"min-block-size: "<<poolAlloc.min_block_size()<<" max-block-size: "<<poolAlloc.max_block_size()<<std::endl;
      oss<<"reserved blocks: "<<usageStats.reserved_blocks<<" consumed blocks: "<<usageStats.consumed_blocks<<std::endl;
      oss<<"reserved bytes: "<<usageStats.reserved_bytes<<" consumed bytes: "<<usageStats.consumed_bytes<<std::endl;
      oss<<"superblock bytes: "<<usageStats.superblock_bytes<<" capacity superblocks: "<<usageStats.capacity_superblocks;
      oss<<" consumed superblocks: "<<usageStats.consumed_superblocks<<std::endl;
      STK_ThrowErrorMsg(oss.str());
    }
  }

  stk::mesh::Entity* m_conn = nullptr;
  stk::mesh::ConnectivityOrdinal* m_ord = nullptr;
  stk::mesh::Permutation* m_perm = nullptr;
  void* m_data = nullptr;
  uint16_t m_capacity = 0;
  uint16_t m_offsets[stk::topology::NUM_RANKS+1] = {0};
};

template<typename NgpSpace>
class MeshConnectivity
{
public:
  using MemSpace = typename NgpSpace::mem_space;
  using PoolAllocatorType = Kokkos::MemoryPool<MemSpace>;

//There are some magic numbers here for setting up the memory pool.
//The pool is organized into superblocks, and it's a bit of a balancing
//act to set the pool up with the right number of superblocks.
//The number of superblocks basically needs to be at least as big as the
//number of powers of two which are less than or contain the max
//individual allocation size.
//
//Each individual allocation actually consumes the smallest power-of-two
//that contains the allocation size. That dictates which superblock the
//allocation is taken from.
//
//So the memory pool needs to be larger than the actual number of bytes that
//we think we need. For sufficiently large problems, a ratio of 2.0 seems
//safe. But for small problems it is much less predictable and 8.0 seems to
//be necessary in some cases. But for small problems it doesn't matter since
//the size of of the pool is much smaller than available memory anyway.
//
//Yes this is unsatisfyingly hand-wavy.
//
  static constexpr size_t bigPoolHeuristic = 32768;
  static constexpr float paddingFactorForLargePool = 2.0;
  static constexpr float paddingFactorForSmallPool = 8.0;
  static constexpr size_t minBlksPerSuperblock = 4;
  static constexpr size_t minNumSuperblocks = 32;

  MeshConnectivity(size_t numEntities, size_t numTotalAlloc, size_t minAllocBlk, size_t maxAllocBlk)
   : m_pool(MemSpace(),
            static_cast<size_t>((numTotalAlloc>bigPoolHeuristic ? paddingFactorForLargePool:paddingFactorForSmallPool)*numTotalAlloc),
            minAllocBlk, maxAllocBlk,
            std::max(static_cast<size_t>(minBlksPerSuperblock*maxAllocBlk),
                     static_cast<size_t>(numTotalAlloc/minNumSuperblocks))),
     m_entityConn("EntityConn",static_cast<size_t>(numEntities*1.4)),
     m_numEntities(numEntities),
     m_updateCounter(0)
  {
  }

  KOKKOS_DEFAULTED_FUNCTION MeshConnectivity() = default;

  KOKKOS_DEFAULTED_FUNCTION MeshConnectivity(const MeshConnectivity&) = default;
  KOKKOS_DEFAULTED_FUNCTION MeshConnectivity(MeshConnectivity&&) = default;
  KOKKOS_DEFAULTED_FUNCTION MeshConnectivity& operator=(const MeshConnectivity&) = default;
  KOKKOS_DEFAULTED_FUNCTION MeshConnectivity& operator=(MeshConnectivity&&) = default;

  KOKKOS_FUNCTION
  ~MeshConnectivity() {
  }

  size_t get_pool_capacity() const
  {
    return m_pool.capacity();
  }

  size_t get_pool_max_block_size() const
  {
    return m_pool.max_block_size();
  }

  void print_state() const
  {
    m_pool.print_state(std::cout);
    typename PoolAllocatorType::usage_statistics usageStats;
    m_pool.get_usage_statistics(usageStats);
    std::cout<<"min-block-size: "<<m_pool.min_block_size()<<" max-block-size: "<<m_pool.max_block_size()<<std::endl;
    std::cout<<"reserved blocks: "<<usageStats.reserved_blocks<<" consumed blocks: "<<usageStats.consumed_blocks<<std::endl;
    std::cout<<"reserved bytes: "<<usageStats.reserved_bytes<<" consumed bytes: "<<usageStats.consumed_bytes<<std::endl;
    std::cout<<"superblock bytes: "<<usageStats.superblock_bytes<<" capacity superblocks: "<<usageStats.capacity_superblocks;
    std::cout<<" consumed superblocks: "<<usageStats.consumed_superblocks<<std::endl;
    std::cout<<" update-count: "<<get_update_count()<<std::endl;
  }

  size_t get_num_entities() const
  { return m_numEntities; }

  void update_num_entities(size_t newNumEntities)
  {
    if (newNumEntities > m_entityConn.extent(0)) {
      size_t newAllocSize = static_cast<size_t>(newNumEntities*1.4);
      Kokkos::resize(m_entityConn, newAllocSize);
    }
    m_numEntities = newNumEntities;
  }

  void reset(stk::mesh::Entity entity, uint16_t newCapacity, bool addPermutations)
  {
    EntityConnectivity& entConn = m_entityConn(entity.local_offset());
    entConn.reset(m_pool, newCapacity, addPermutations);
  }

  size_t get_update_count() const { return m_updateCounter; }

  void set_update_count(size_t newUpdateCount) { m_updateCounter = newUpdateCount; }

  bool is_empty(stk::mesh::Entity entity) const
  {
    STK_NGP_ThrowAssert(entity.local_offset() < m_numEntities);
    return m_entityConn(entity.local_offset()).is_empty();
  }

  void set_capacity(stk::mesh::Entity entity, uint16_t capacity, bool addPermutations)
  {
    EntityConnectivity& entConn = m_entityConn(entity.local_offset());
    entConn.grow_capacity(m_pool, capacity, addPermutations);
  }

  void add_connectivity(stk::mesh::Entity entity,
                        stk::mesh::EntityRank rank,
                        unsigned numConnectivity,
                        const stk::mesh::Entity* entities,
                        const stk::mesh::ConnectivityOrdinal* ords,
                        const stk::mesh::Permutation* perms,
                        bool addPermutations)
  {
    EntityConnectivity& entConn = m_entityConn(entity.local_offset());
    entConn.add_connectivity(m_pool, rank, numConnectivity, entities, ords, perms, addPermutations);
  }

  KOKKOS_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::Entity entity, stk::mesh::EntityRank connectedRank) const
  {
    STK_NGP_ThrowAssert(entity.local_offset() < m_numEntities);
    return m_entityConn[entity.local_offset()].get_connected_entities(connectedRank);
  }

  KOKKOS_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::Entity entity, stk::mesh::EntityRank connectedRank) const
  {
    STK_NGP_ThrowAssert(entity.local_offset() < m_numEntities);
    return m_entityConn[entity.local_offset()].get_connected_ordinals(connectedRank);
  }

  KOKKOS_INLINE_FUNCTION
  ConnectedPermutations get_connected_permutations(stk::mesh::Entity entity, stk::mesh::EntityRank connectedRank) const
  {
    STK_NGP_ThrowAssert(entity.local_offset() < m_numEntities);
    return m_entityConn[entity.local_offset()].get_connected_permutations(connectedRank);
  }

  const PoolAllocatorType& get_memory_pool() const
  {
    return m_pool;
  }

private:
  PoolAllocatorType m_pool;
  Kokkos::View<EntityConnectivity*,MemSpace> m_entityConn;
  size_t m_numEntities = 0;
  size_t m_updateCounter = 0;
};

}
}
}

#endif
