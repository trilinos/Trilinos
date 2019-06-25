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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/baseImpl/check_comm_list.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>  // for EntityCommDatabase
#include <stk_mesh/base/EntityCommListInfo.hpp>
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityCommInfo
#include <stk_util/parallel/Parallel.hpp>

namespace {

std::vector<int> all_procs_except_local_proc(int num_procs, int local_proc)
{
  std::vector<int> all_other_procs;
  for(int p=0; p<local_proc; ++p) all_other_procs.push_back(p);
  for(int p=local_proc+1; p<num_procs; ++p) all_other_procs.push_back(p);
  return all_other_procs;
}

TEST(CheckCommList, test_shared)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int num_procs = stk::parallel_machine_size(communicator);
    if (num_procs == 3) {

        int local_proc = stk::parallel_machine_rank(communicator);
        int owner = 1;
        stk::mesh::EntityCommListInfoVector comm_list;
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, 99);
        stk::mesh::Entity entity;//don't need valid entity for this test
        std::vector<int> sharing_procs = all_procs_except_local_proc(num_procs, local_proc);
        stk::mesh::EntityComm entity_comm;
        entity_comm.owner_rank = owner;
        entity_comm.comm_map.push_back(stk::mesh::EntityCommInfo(0, sharing_procs[0]));
        entity_comm.comm_map.push_back(stk::mesh::EntityCommInfo(0, sharing_procs[1]));
        stk::mesh::EntityCommListInfo comm_info = {key, entity, nullptr, 0, owner, &entity_comm};
        comm_list.push_back(comm_info);

        EXPECT_TRUE(stk::mesh::impl::is_comm_list_globally_consistent(communicator, comm_list));
    }
}

TEST(CheckCommList, test_ghosted)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int num_procs = stk::parallel_machine_size(communicator);
    if (num_procs == 3) {

        int local_proc = stk::parallel_machine_rank(communicator);
        int owner = 1;
        stk::mesh::EntityCommListInfoVector comm_list;
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, 99);
        stk::mesh::Entity entity;//don't need valid entity for this test
        stk::mesh::EntityComm entity_comm;
        entity_comm.owner_rank = owner;
        if (local_proc == owner) {
            entity_comm.comm_map.push_back(stk::mesh::EntityCommInfo(3, 0));
            entity_comm.comm_map.push_back(stk::mesh::EntityCommInfo(3, 2));
        }
        else {
            entity_comm.comm_map.push_back(stk::mesh::EntityCommInfo(3, owner));
        }
        stk::mesh::EntityCommListInfo comm_info = {key, entity, nullptr, 0, owner, &entity_comm};
        comm_list.push_back(comm_info);

        EXPECT_TRUE(stk::mesh::impl::is_comm_list_globally_consistent(communicator, comm_list));
    }
}

}

