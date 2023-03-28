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

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_sum
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for has_superset, Bucket
#include "stk_mesh/base/Types.hpp"      // for EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }



namespace stk {
namespace mesh {

void fillNumEntitiesPerRankOnThisProc(const BulkData & M, std::vector<size_t>&local, const Selector *selector)
{
    const MetaData & meta = M.mesh_meta_data();
    const EntityRank numEntityRanks = meta.entity_rank_count();
    local.clear();
    local.resize(numEntityRanks,0);

    Selector owns = meta.locally_owned_part();
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

    all_reduce_sum(M.parallel(), localCounts.data(), globalCounts.data(), numEntityRanks);

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

    all_reduce_sum(M.parallel(), localEntityCounts.data(), globalCounts.data(), numEntityRanks);

    all_reduce_min(M.parallel(), localEntityCounts.data(), min_counts.data(), numEntityRanks);
    all_reduce_max(M.parallel(), localEntityCounts.data(), max_counts.data(), numEntityRanks);

    return;
}

//----------------------------------------------------------------------

struct ghostObj {
  stk::topology::rank_t entityRank;
  stk::mesh::EntityId entityId;
  int destProc;
};

stk::mesh::EntityProcVec
send_non_owned_entities_to_owner(stk::mesh::BulkData& stkBulk, const stk::mesh::EntityProcVec& entities)
{
  int numProcs = stkBulk.parallel_size();

  std::vector<std::vector<ghostObj> > ghostRequestSend(numProcs);

  stk::mesh::EntityProcVec owned_entities;

  for(const stk::mesh::EntityProc& entityAndProc : entities) {
    stk::mesh::Entity entity = entityAndProc.first;
    int ownerRank = stkBulk.parallel_owner_rank(entity);
    if(ownerRank != stkBulk.parallel_rank()) {
      ghostObj curGhost; 
      curGhost.entityRank = stkBulk.entity_rank(entity);
      curGhost.entityId = stkBulk.identifier(entity);
      curGhost.destProc = entityAndProc.second;
      ghostRequestSend[ownerRank].push_back(curGhost);
    }      
    else {
      owned_entities.push_back(entityAndProc);
    }      
  }

  std::vector<std::vector<ghostObj> > ghostRequestRecv(numProcs);
  MPI_Comm mpi_comm(stkBulk.parallel()); 
  stk::parallel_data_exchange_t(ghostRequestSend, ghostRequestRecv, mpi_comm);

  for(int iproc = 0; iproc < numProcs; ++iproc) {
    for(unsigned j = 0; j < ghostRequestRecv[iproc].size(); ++j) {
      ghostObj& curGhost = ghostRequestRecv[iproc][j];
      stk::mesh::Entity curEntity = stkBulk.get_entity(curGhost.entityRank,curGhost.entityId);
      STK_ThrowAssert(stkBulk.is_valid(curEntity));
      STK_ThrowAssert(stkBulk.bucket(curEntity).owned());
      owned_entities.emplace_back(curEntity, curGhost.destProc);
    }      
  }
  return owned_entities;
}

} // namespace mesh
} // namespace stk

