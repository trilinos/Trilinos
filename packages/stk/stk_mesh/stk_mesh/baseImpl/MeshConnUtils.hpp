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
// 

#ifndef stk_mesh_base_impl_MeshConnUtils_hpp
#define stk_mesh_base_impl_MeshConnUtils_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/MeshConnectivity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>

#include <vector>
#include <algorithm>
#include <tuple>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

std::tuple<size_t,size_t,size_t> connectivity_memory_pool_sizes(const stk::mesh::BulkData& bulk);

size_t count_entity_connectivity(const stk::mesh::BulkData& bulk,
                                 stk::mesh::EntityRank endRank, stk::mesh::Entity entity);

template<typename MeshConnType>
void fill_mesh_connectivity(const BulkData& bulk, MeshConnType& meshConn)
{
  const size_t numEntities = bulk.get_size_of_entity_index_space();
  auto [numTotalAlloc, minAllocBlk, maxAllocBlk] = stk::mesh::impl::connectivity_memory_pool_sizes(bulk);

  if (meshConn.get_num_entities() > 0) {
    //if not empty, destroy existing connectivity to keep high-water-mark down.
    meshConn = MeshConnType(meshConn.get_num_entities(), 512, 16, 64);
  }

  //Set some minimums since the memory-pool in the MeshConnectivity object
  //gets weird if zeros get passed in. All "real" meshes should have values
  //larger than these minimums, and if we do use these minimums then the
  //memory amounts being used are very small anyway.
  numTotalAlloc = std::max(std::size_t{2048},numTotalAlloc);
  minAllocBlk = std::max(std::size_t{16}, minAllocBlk);
  maxAllocBlk = std::max(std::size_t{64}, maxAllocBlk);

  meshConn = MeshConnType(numEntities, numTotalAlloc, minAllocBlk, maxAllocBlk);
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());

  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    const bool addPermutations = stk::mesh::does_rank_have_valid_permutations(rank);
    stk::mesh::for_each_entity_run(bulk, rank, meta.universal_part(),
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity ent)
      {
        size_t num = stk::mesh::impl::count_entity_connectivity(mesh, endRank, ent);
        meshConn.set_capacity(ent, num, addPermutations);
        for(stk::mesh::EntityRank connRank = stk::topology::NODE_RANK; connRank<endRank; ++connRank) {
          auto conn = mesh.get_connected_entities(ent, connRank);
          if (conn.size() > 0) {
            const auto* ords = mesh.begin_ordinals(ent, connRank);
            const auto* perms = mesh.begin_permutations(ent, connRank);
            const auto* endPerms = mesh.end_permutations(ent, connRank);
            const auto* permsPtr = (endPerms-perms) > 0 ? perms : nullptr;
            meshConn.add_connectivity(ent, connRank, conn.size(), conn.data(), ords, permsPtr, addPermutations);
          }
        }
      });
  }
  meshConn.set_update_count(bulk.synchronized_count());
}

template<typename MeshConnType>
void reset_deleted_entities(const BulkData& mesh, MeshConnType& meshConn)
{
  const size_t numMeshConnEntities = meshConn.get_num_entities();
  for(size_t entOffset=1; entOffset<numMeshConnEntities; ++entOffset) {
    Entity ent(entOffset);
    if (!mesh.is_valid(ent)) {
      meshConn.reset(ent, 0, true);
    }
  }
}

constexpr inline size_t sizeofEntityOrdinalPerm = sizeof(stk::mesh::Entity)
                           + sizeof(stk::mesh::ConnectivityOrdinal)
                           + sizeof(stk::mesh::Permutation);

inline
unsigned int_power_of_two_that_contains(unsigned N)
{
//borrowed this from Kokkos, because the Kokkos one is private/inaccessible
  return N ? Kokkos::bit_width(N-1) : 0;
}

inline
size_t get_alloc_block_size(size_t numBytes)
{
  return (1LU << int_power_of_two_that_contains(numBytes));
}

void add_or_increment_alloc_blk(size_t allocSize, size_t numBlks,
                                std::vector<std::pair<size_t,size_t>>& allocBlockSizes);
template<typename MeshConnType>
std::pair<bool,size_t>
is_connectivity_different(const stk::mesh::BulkData& bulk,
                          stk::mesh::Entity entity,
                          EntityRank endRank,
                          const MeshConnType& meshConn)
{
  STK_ThrowAssertMsg(entity.local_offset() < meshConn.get_num_entities(),"entity out of range for MeshConnectivity");
  if (meshConn.is_empty(entity)) {
    return std::make_pair(true, count_entity_connectivity(bulk, endRank, entity));
  }

  const auto& meshIndex = bulk.mesh_index(entity);
  bool connectivityIsDifferent = false;
  size_t numConnectivity = 0;
  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    auto conn = meshIndex.bucket->get_connected_entities(meshIndex.bucket_ordinal, rank);
    const unsigned connSize = conn.size();
    numConnectivity += connSize;
    auto mc_conn = meshConn.get_connected_entities(entity, rank);
    if (connSize != mc_conn.size()) {
      connectivityIsDifferent = true;
      continue;
    }
    if (!connectivityIsDifferent && connSize > 0) {
      const auto* ords = meshIndex.bucket->begin_ordinals(meshIndex.bucket_ordinal, rank);
      const auto* perms = meshIndex.bucket->begin_permutations(meshIndex.bucket_ordinal, rank);
      const auto* endPerms = meshIndex.bucket->end_permutations(meshIndex.bucket_ordinal, rank);
      const auto* permsPtr = (endPerms-perms) > 0 ? perms : nullptr;

      auto mc_ords = meshConn.get_connected_ordinals(entity, rank);
      auto mc_perms = meshConn.get_connected_permutations(entity, rank);
      if ((permsPtr != nullptr) && (mc_perms.size() == 0)) {
        connectivityIsDifferent = true;
        break;
      }
      for(unsigned i=0; i<connSize; ++i) {
        if (conn[i] != mc_conn[i] || ords[i] != mc_ords[i] ||
            (permsPtr != nullptr && perms[i] != mc_perms[i])) {
          connectivityIsDifferent = true;
          break;
        }
      }
    }
  }

  return std::make_pair(connectivityIsDifferent,numConnectivity);
}

template<typename MeshConnType>
void get_alloc_block_if_needed(const stk::mesh::BulkData& mesh,
                               MeshConnType& meshConn,
                               size_t prevNumMeshEntities,
                               stk::mesh::Entity ent,
                               stk::mesh::EntityRank endRank,
                               std::vector<std::pair<size_t,size_t>>& allocBlockSizes,
                               std::vector<stk::mesh::Entity>& modifiedEntities)
{
  STK_ThrowAssertMsg(mesh.is_valid(ent),"get_alloc_block_..., invalid entity "<<mesh.entity_key(ent)<<", offset="<<ent.local_offset());
  bool connectivityIsDifferent = false;
  size_t numConn = 0;
  if (ent.local_offset() >= prevNumMeshEntities) {
    connectivityIsDifferent = true;
    numConn = count_entity_connectivity(mesh, endRank, ent);
  }
  else {
    std::tie(connectivityIsDifferent, numConn) = is_connectivity_different(mesh, ent, endRank, meshConn);
  }

  if (connectivityIsDifferent) {
    modifiedEntities.push_back(ent);
    meshConn.reset(ent, 0, true);

    if (numConn > 0) {
      const size_t allocSize = get_alloc_block_size(numConn*sizeofEntityOrdinalPerm);
      add_or_increment_alloc_blk(allocSize, 1, allocBlockSizes);
    }
  }
}

template<typename MeshConnType>
std::tuple<std::vector<std::pair<size_t,size_t>>,std::vector<stk::mesh::Entity>>
dealloc_and_get_needed_allocs_full_mesh(const stk::mesh::BulkData& bulk,
                                        MeshConnType& meshConn)
{
  reset_deleted_entities(bulk, meshConn);

  const size_t numBulkEntities = bulk.get_size_of_entity_index_space();
  const size_t oldNumMeshConnEntities = meshConn.get_num_entities();
  if (numBulkEntities != oldNumMeshConnEntities) {
    meshConn.update_num_entities(numBulkEntities);
  }

  std::vector<std::pair<size_t,size_t>> allocBlockSizes;
  std::vector<stk::mesh::Entity> modifiedEntities;

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());

  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    stk::mesh::for_each_entity_run(bulk, rank, meta.universal_part(),
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity ent)
      {
        get_alloc_block_if_needed(mesh, meshConn, oldNumMeshConnEntities, ent,
                                  endRank, allocBlockSizes, modifiedEntities);
      });
  }

  return std::make_tuple(allocBlockSizes, modifiedEntities);
}

template<typename MeshConnType>
std::tuple<std::vector<std::pair<size_t,size_t>>,std::vector<stk::mesh::Entity>>
dealloc_and_get_needed_allocs_using_entity_states(const stk::mesh::BulkData& bulk,
                                                  MeshConnType& meshConn)
{
  reset_deleted_entities(bulk, meshConn);

  const size_t numBulkEntities = bulk.get_size_of_entity_index_space();
  const size_t oldNumMeshConnEntities = meshConn.get_num_entities();
  if (numBulkEntities != oldNumMeshConnEntities) {
    meshConn.update_num_entities(numBulkEntities);
  }

  std::vector<std::pair<size_t,size_t>> allocBlockSizes;
  std::vector<stk::mesh::Entity> modifiedEntities;

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());

  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    stk::mesh::for_each_entity_run(bulk, rank, meta.universal_part(),
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity ent)
      {
        if (mesh.state(ent) != Unchanged) {
          get_alloc_block_if_needed(mesh, meshConn, oldNumMeshConnEntities, ent,
                                    endRank, allocBlockSizes, modifiedEntities);
        }
      });
  }

  return std::make_tuple(allocBlockSizes, modifiedEntities);
}

template<typename MeshConnType>
void update_mesh_connectivity(const BulkData& mesh, MeshConnType& meshConn,
                              const std::vector<Entity>& modifiedEntities)
{
  const size_t numEntities = mesh.get_size_of_entity_index_space();
  if (numEntities > meshConn.get_num_entities()) {
    meshConn.update_num_entities(numEntities);
  }

  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());

  for(Entity ent : modifiedEntities) {
    if (mesh.state(ent) == Deleted) {
      continue;
    }

    STK_ThrowRequireMsg(mesh.is_valid(ent), "Error in update_mesh_connectivity, invalid entity "<<mesh.entity_key(ent));

    EntityRank rank = mesh.entity_rank(ent);
    const bool addPermutations = stk::mesh::does_rank_have_valid_permutations(rank);
    size_t num = stk::mesh::impl::count_entity_connectivity(mesh, endRank, ent);
    meshConn.set_capacity(ent, num, addPermutations);

    for(stk::mesh::EntityRank connRank = stk::topology::NODE_RANK; connRank<endRank; ++connRank) {
      auto conn = mesh.get_connected_entities(ent, connRank);
      if (conn.size() > 0) {
        const auto* ords = mesh.begin_ordinals(ent, connRank);
        const auto* perms = mesh.begin_permutations(ent, connRank);
        const auto* endPerms = mesh.end_permutations(ent, connRank);
        const auto* permsPtr = (endPerms-perms) > 0 ? perms : nullptr;
        meshConn.add_connectivity(ent, connRank, conn.size(), conn.data(), ords, permsPtr, addPermutations);
      }
    }
  }

  meshConn.set_update_count(mesh.synchronized_count());
}

template<typename MemPoolType>
std::vector<std::pair<size_t,size_t>>
get_available_blocks(const MemPoolType& memPool)
{
  std::vector<std::pair<size_t,size_t>> availableBlocks;
  const int numSuperBlocks = memPool.number_of_superblocks();
  for(int n=0; n<numSuperBlocks; ++n) {
    int block_size = 0, block_count_capacity = 0, block_count_used = 0;
    memPool.superblock_state(n, block_size, block_count_capacity, block_count_used);
    if (block_count_capacity > block_count_used) {
      size_t numAvail = static_cast<size_t>(block_count_capacity - block_count_used);
      add_or_increment_alloc_blk(block_size, numAvail, availableBlocks);
    }
  }

  return availableBlocks;
}

bool blocks_are_available(std::vector<std::pair<size_t, size_t>>& allocBlocksNeeded,
                          std::vector<std::pair<size_t, size_t>>& allocBlocksAvailable);

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_base_impl_MeshConnUtils_hpp

