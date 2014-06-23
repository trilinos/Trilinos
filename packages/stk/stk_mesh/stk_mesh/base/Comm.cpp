/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

