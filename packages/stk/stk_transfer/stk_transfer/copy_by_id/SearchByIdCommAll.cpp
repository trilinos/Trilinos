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

#include "SearchByIdCommAll.hpp"
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <sstream>

namespace stk {
namespace transfer {

void SearchByIdCommAll::do_search(const TransferCopyByIdMeshAdapter & mesha,
                                  const TransferCopyByIdMeshAdapter & meshb,
                                  KeyToTargetProcessor & key_to_target_processor)
{
  key_to_target_processor.clear();
  m_remote_keys.clear();
  const TransferCopyByIdMeshAdapter::MeshIDVector & meshb_ids = meshb.get_mesh_ids();
  const ParallelMachine comm = meshb.comm();
  const int my_proc = parallel_machine_rank(comm);
  const int num_proc = parallel_machine_size(comm);
  stk::CommSparse commSparse(comm);

  for (int phase = 0; phase < 2; ++phase) {
    for (size_t id_index = 0; id_index < meshb_ids.size(); ++id_index) {
      const Mesh_ID key = meshb_ids[id_index];
      if (mesha.is_locally_owned(key)) {
        key_to_target_processor.emplace_back(key, my_proc);
        continue;
      }
      m_remote_keys.insert(key);
      for (int send_proc = 0; send_proc < num_proc; ++send_proc) {
        if (my_proc == send_proc) { continue; }
        commSparse.send_buffer(send_proc).pack<Mesh_ID>(key);
      }
    }

    if (phase == 0) {
      commSparse.allocate_buffers();
    }
    else {
      commSparse.communicate();
    }
  }

  for (int recv_proc = 0; recv_proc < num_proc; ++recv_proc) {
    if (my_proc != recv_proc) {
      while (commSparse.recv_buffer(recv_proc).remaining()) {
        Mesh_ID key;
        commSparse.recv_buffer(recv_proc).unpack<Mesh_ID>(key);
        if (mesha.is_locally_owned(key)) {
          key_to_target_processor.emplace_back(key, recv_proc);
        }
      }
    }
  }
}


}  } // namespace transfer stk
