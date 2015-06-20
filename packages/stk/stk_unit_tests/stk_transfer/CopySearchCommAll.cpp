// Copyright (c) 2015, Sandia Corporation.
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

#include "CopySearchCommAll.hpp"
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <sstream>

namespace stk {
namespace transfer {

void CopySearchCommAll::do_search(const CopyTransferMeshBase & mesha,
                                  const CopyTransferMeshBase & meshb,
                                  KeyToTargetProcessor & key_to_target_processor
                                  )
{
  key_to_target_processor.clear();
  m_remote_keys.clear();
  const CopyTransferMeshBase::MeshIDVector & meshb_ids = meshb.get_mesh_ids();
  const ParallelMachine comm = meshb.comm();
  const int my_proc = parallel_machine_rank(comm);
  const int num_proc = parallel_machine_size(comm);
  stk::CommAll commAll(comm);
  for (int phase=0;phase<2;++phase)
  {
    for (size_t id_index=0 ; id_index<meshb_ids.size() ; ++id_index)
    {
      const Mesh_ID key = meshb_ids[id_index];
      if (mesha.is_locally_owned(key))
      {
        key_to_target_processor[key]=my_proc;
        continue;
      }
      m_remote_keys.insert(key);
      for (int send_proc = 0 ; send_proc < num_proc ; ++send_proc)
      {
        if (my_proc == send_proc) { continue; }
        commAll.send_buffer(send_proc).pack<Mesh_ID>(key);
      }
    }

    if (phase == 0 )
    {
      commAll.allocate_buffers(num_proc/4);
    }
    else
    {
      commAll.communicate();
    }
  }

  for (int recv_proc=0;recv_proc<num_proc;++recv_proc)
  {
    if ( my_proc != recv_proc )
    {
      while(commAll.recv_buffer(recv_proc).remaining())
      {
        Mesh_ID key;
        commAll.recv_buffer(recv_proc).unpack<Mesh_ID>(key);
        if (mesha.is_locally_owned(key))
        {
          key_to_target_processor[key] = recv_proc;
        }
      }
    }
  }
}


}  } // namespace transfer stk
