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

#include <stk_mesh/baseImpl/MeshConnUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <vector>
#include <utility>
#include <algorithm>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

size_t count_entity_connectivity(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank endRank, stk::mesh::Entity entity)
{
  const MeshIndex& mi = bulk.mesh_index(entity);
  size_t numConnectivity = 0;
  for(stk::mesh::EntityRank connRank = stk::topology::NODE_RANK; connRank<endRank; ++connRank) {
    numConnectivity += mi.bucket->num_connectivity(mi.bucket_ordinal, connRank);
  }
  return numConnectivity;
}

std::tuple<size_t,size_t,size_t> connectivity_memory_pool_sizes(const stk::mesh::BulkData& bulk)
{
  size_t minAllocBlk = std::numeric_limits<size_t>::max();
  size_t maxAllocBlk = 0;
  size_t numTotalAlloc = 0;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());
  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    stk::mesh::for_each_entity_run(bulk, rank, meta.universal_part(),
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity ent)
      {
        const size_t numConn = count_entity_connectivity(mesh, endRank, ent);
        const size_t numBytes = numConn*sizeofEntityOrdinalPerm;
        const size_t pwr2 = int_power_of_two_that_contains(numBytes);
        const size_t allocSize = 1LU << pwr2;
        numTotalAlloc += allocSize;
        minAllocBlk = std::min(minAllocBlk, allocSize);
        maxAllocBlk = std::max(maxAllocBlk, allocSize);
      }
    );
  }

  if (minAllocBlk == std::numeric_limits<size_t>::max()) {
    minAllocBlk = 0;
  }

  return std::make_tuple(numTotalAlloc, minAllocBlk, maxAllocBlk);
}

void add_or_increment_alloc_blk(size_t allocSize, size_t numBlks,
                                std::vector<std::pair<size_t,size_t>>& allocBlockSizes)
{
  bool found = false;

  for(auto& sizeAndCount : allocBlockSizes) {
    if (sizeAndCount.first == allocSize) {
      sizeAndCount.second += numBlks;
      found = true;
      break;
    }
  }
  
  if (!found) {
    stk::util::insert_keep_sorted_and_unique(std::make_pair(allocSize,numBlks), allocBlockSizes);
  }
}

bool blocks_are_available(std::vector<std::pair<size_t, size_t>>& allocBlocksNeeded,
                          std::vector<std::pair<size_t, size_t>>& allocBlocksAvailable)
{

  for(auto& needSizeAndCount : allocBlocksNeeded) {
    const size_t neededBlkSize = needSizeAndCount.first;
    size_t& numNeeded = needSizeAndCount.second;

    bool foundEnough = false;
    for(auto& availSizeAndCount : allocBlocksAvailable) {
      const size_t availBlkSize = availSizeAndCount.first;

      if (availBlkSize >= neededBlkSize) {
        size_t& numAvail = availSizeAndCount.second;

        if (numAvail > 0) {
          size_t numToSubtract = std::min(numAvail, numNeeded);
          numAvail -= numToSubtract;
          numNeeded -= numToSubtract;

          if (numNeeded == 0) {
            foundEnough = true;
            break;
          }
        }
      }
    }

    if (!foundEnough) {
      return false;
    }
  }

  return true;
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

