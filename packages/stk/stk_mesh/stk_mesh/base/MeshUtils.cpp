// Copyright (c) 2014, Sandia Corporation.
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

#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_topology/topology.hpp>

namespace stk {
namespace mesh {

void find_ghosted_nodes_that_need_to_be_shared(const stk::mesh::BulkData & bulk, stk::mesh::EntityVector& ghosted_nodes_that_are_now_shared)
{
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
    if (endRank >= stk::topology::END_RANK)
    {
        endRank = stk::topology::END_RANK;
    }

    for (stk::mesh::EntityRank rank=stk::topology::EDGE_RANK; rank<endRank; ++rank)
    {
        const stk::mesh::BucketVector& entity_buckets = bulk.buckets(rank);
        for(size_t i=0; i<entity_buckets.size(); ++i)
        {
            const stk::mesh::Bucket& bucket = *entity_buckets[i];
            if ( bucket.owned() )
            {
                for(size_t n=0; n<bucket.size(); ++n)
                {
                    const stk::mesh::Entity * nodes = bulk.begin_nodes(bucket[n]);
                    unsigned num_nodes = bulk.num_nodes(bucket[n]);
                    for (unsigned j=0;j<num_nodes;++j)
                    {
                        if (bulk.in_receive_ghost(bulk.entity_key(nodes[j])))
                        {
                            ghosted_nodes_that_are_now_shared.push_back(nodes[j]);
                        }
                    }
                }
            }
        }
    }

    stk::util::sort_and_unique(ghosted_nodes_that_are_now_shared);
}

//----------------------------------------------------------------------

void fixup_ghosted_to_shared_nodes(stk::mesh::BulkData & bulk)
{
    stk::mesh::EntityVector ghosted_nodes_that_are_now_shared;
    find_ghosted_nodes_that_need_to_be_shared(bulk, ghosted_nodes_that_are_now_shared);

    stk::CommSparse comm(bulk.parallel());

    for (int phase=0;phase<2;++phase)
    {
        for (size_t i = 0; i < ghosted_nodes_that_are_now_shared.size(); ++i)
        {
            stk::mesh::Entity node = ghosted_nodes_that_are_now_shared[i];
            int proc = bulk.parallel_owner_rank(node);
            comm.send_buffer(proc).pack<stk::mesh::EntityKey>(bulk.entity_key(node));
        }
        if (phase == 0 )
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    stk::mesh::EntityVector sharedNodes;
    for (int process=0;process<bulk.parallel_size();++process)
    {
        while(comm.recv_buffer(process).remaining())
        {
            stk::mesh::EntityKey key;
            comm.recv_buffer(process).unpack<stk::mesh::EntityKey>(key);

            stk::mesh::Entity entity = bulk.get_entity(key);
            if ( bulk.state(entity) != stk::mesh::Deleted && bulk.is_valid(entity) )
            {
                bulk.add_node_sharing(entity, process);
                sharedNodes.push_back(entity);
            }
        }
    }
/////////////////////////

    stk::CommSparse commSecondStage(bulk.parallel());
    for (int phase=0;phase<2;++phase)
    {
        for (size_t i=0;i<sharedNodes.size();++i)
        {
            std::vector<int> procs;
            stk::mesh::EntityKey key = bulk.entity_key(sharedNodes[i]);
            bulk.comm_shared_procs(key, procs);
            for (size_t j=0;j<procs.size();++j)
            {
                if ( procs[j] != bulk.parallel_rank() )
                {
                    commSecondStage.send_buffer(procs[j]).pack<int>(bulk.parallel_rank()).pack<stk::mesh::EntityKey>(key);
                    for (size_t k=0;k<procs.size();++k)
                    {
                        commSecondStage.send_buffer(procs[j]).pack<int>(procs[k]).pack<stk::mesh::EntityKey>(key);
                    }
                }
            }
        }
        if (phase == 0 )
        {
            commSecondStage.allocate_buffers();
        }
        else
        {
            commSecondStage.communicate();
        }
    }

    for (int proc_that_sent_message=0;proc_that_sent_message<bulk.parallel_size();++proc_that_sent_message)
    {
        if ( proc_that_sent_message == bulk.parallel_rank() ) continue;
        while(commSecondStage.recv_buffer(proc_that_sent_message).remaining())
        {
            stk::mesh::EntityKey key;
            int sharingProc;
            commSecondStage.recv_buffer(proc_that_sent_message).unpack<int>(sharingProc).unpack<stk::mesh::EntityKey>(key);
            if ( sharingProc != bulk.parallel_rank() )
            {
                stk::mesh::Entity entity = bulk.get_entity(key);
                if ( bulk.state(entity) != stk::mesh::Deleted && bulk.is_valid(entity) && !bulk.in_shared(key, sharingProc) )
                {
                    bulk.add_node_sharing(entity, sharingProc);
                }
            }
        }
    }
}

} // namespace mesh
} // namespace stk
