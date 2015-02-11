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

#include "CopyTransfer.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <sstream>

namespace stk {
namespace transfer {


void CopyTransfer::do_transfer(const KeyToTargetProcessor & key_to_target_processor, const CopyTransferMeshBase & mesha, CopyTransferMeshBase & meshb)
{
  const unsigned numValsb = meshb.num_fields();
  const ParallelMachine comm = mesha.comm();
  const int my_proc = parallel_machine_rank(comm);
  const int num_proc = parallel_machine_size(comm);
  stk::CommSparse commSparse(comm);
  MeshIDSet remote_keys = m_search.get_remote_keys();
  for (int phase=0;phase<2;++phase)
  {
    KeyToTargetProcessor::const_iterator map_it = key_to_target_processor.begin();
    for ( ; map_it!=key_to_target_processor.end() ;  ++map_it) {
      const Mesh_ID key = map_it->first;
      const int target_proc = map_it->second;
      if (target_proc == my_proc)
      {
        if (0 == phase) {
          // copy directly from A to B:
          for (unsigned f=0 ; f<numValsb ; ++f) {
            const double * f_dataA = mesha.field_data(key,f);
            double * f_dataB = meshb.field_data(key,f);
            const unsigned this_field_size = mesha.field_data_size(key,f);
            for (unsigned index=0 ; index<this_field_size ; ++index) {
              f_dataB[index] = f_dataA[index];
            }
          }
        }
        continue;
      }
      commSparse.send_buffer(target_proc).pack<Mesh_ID>(key);
      for (unsigned f=0; f<numValsb; ++f)  {
        const unsigned this_field_size = mesha.field_data_size(key,f);
        // FIX to allow copying of int
        const double * f_data = mesha.field_data(key,f);
        for (unsigned index=0 ; index<this_field_size ; ++index) {
          commSparse.send_buffer(target_proc).pack<double>(f_data[index]);
        }
      }
    }

    if (phase == 0 )
    {
      commSparse.allocate_buffers();
    }
    else
    {
      commSparse.communicate();
    }
  }

  std::ostringstream error_msg;
  unsigned error_count = 0;
  for (int recv_proc=0;recv_proc<num_proc;++recv_proc)
  {
    if ( my_proc != recv_proc )
    {
      while(commSparse.recv_buffer(recv_proc).remaining())
      {
        Mesh_ID key;
        commSparse.recv_buffer(recv_proc).unpack<Mesh_ID>(key);
        if (remote_keys.find(key) == remote_keys.end()) {
          ++error_count;
          error_msg << "P" << my_proc << " Error, proc = " << recv_proc << " sent unrequested data for key = " << key << std::endl;
        } else {
          remote_keys.erase(key);
        }
        for (unsigned f=0 ; f<numValsb ; ++f) {
          double * f_data = meshb.field_data(key,f);
          const unsigned this_field_size = meshb.field_data_size(key,f);
          for (unsigned index=0 ; index<this_field_size ; ++index) {
            commSparse.recv_buffer(recv_proc).unpack<double>(f_data[index]);
          }
        }
      }
    }
  }
  if (!remote_keys.empty()) {
    error_msg << "P" << my_proc << " Error, Did not receive all the requested keys from other processors, unfulfilled keys = [";
    for (MeshIDSet::const_iterator set_it=remote_keys.begin() ; set_it != remote_keys.end() ; ++set_it) {
      error_msg << meshb.print_mesh_id(*set_it) << " ";
    }
    error_msg << "]" << std::endl;
    ++error_count;
  }
  all_reduce( comm , ReduceSum<1>( & error_count ) );
  if ( error_count ) {
    all_write_string( comm , std::cerr , error_msg.str() );

    ThrowErrorMsg("Error in communication during CopyTransfer!\n");
  }
}

}  } // namespace transfer stk
