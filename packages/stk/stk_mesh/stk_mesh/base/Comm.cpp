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
// 

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_sum
#include "stk_mesh/base/Bucket.hpp"     // for has_superset, Bucket
#include "stk_mesh/base/Types.hpp"      // for EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }



namespace stk {
namespace mesh {

void fillNumEntitiesPerRankOnThisProc(const BulkData & M, std::vector<size_t>&local, const Selector *selector)
{
    const MetaData & S = MetaData::get(M);
    const EntityRank entity_rank_count = static_cast<EntityRank>(S.entity_rank_count());
    size_t numEntityRanks = entity_rank_count;
    local.clear();
    local.resize(numEntityRanks,0);

    Selector owns = S.locally_owned_part();
    for(EntityRank i = stk::topology::NODE_RANK; i < numEntityRanks; ++i)
    {
        const BucketVector & ks = M.buckets(i);
        BucketVector::const_iterator ik;

        for(ik = ks.begin(); ik != ks.end(); ++ik)
        {
            if(selector && !(*selector)(**ik))
                continue;
            if(owns(**ik))
            {
                local[i] += (*ik)->size();
            }
        }
    }
}

void comm_mesh_counts( const BulkData & M ,
                       std::vector<size_t> & globalCounts ,
                       const Selector *selector)
{
    std::vector<size_t> localCounts;
    fillNumEntitiesPerRankOnThisProc(M, localCounts, selector);

    size_t numEntityRanks = localCounts.size();
    globalCounts.resize(numEntityRanks, 0);

    all_reduce_sum(M.parallel(), &localCounts[0], &globalCounts[0], numEntityRanks);

    return;
}

void comm_mesh_counts( const BulkData & M ,
                       std::vector<size_t> & globalCounts,
                       std::vector<size_t> & min_counts,
                       std::vector<size_t> & max_counts,
                       const Selector *selector)
{
    std::vector<size_t> localEntityCounts;
    fillNumEntitiesPerRankOnThisProc(M, localEntityCounts, selector);

    size_t numEntityRanks = localEntityCounts.size();
    globalCounts.resize(numEntityRanks, 0);
    min_counts.resize(numEntityRanks, 0);
    max_counts.resize(numEntityRanks, 0);

    all_reduce_sum(M.parallel(), &localEntityCounts[0], &globalCounts[0], numEntityRanks);

    all_reduce_min(M.parallel(), &localEntityCounts[0], &min_counts[0], numEntityRanks);
    all_reduce_max(M.parallel(), &localEntityCounts[0], &max_counts[0], numEntityRanks);

    return;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

