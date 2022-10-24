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

#include "stk_util/parallel/BroadcastArg.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for MPI_Bcast, parallel_machine_rank, MPI_BYTE
#include "stk_util/stk_config.h"           // for STK_HAS_MPI
#include <cstddef>                         // for size_t
#include <algorithm>                       // for copy
#include <string>                          // for string, basic_string


namespace stk {

BroadcastArg::BroadcastArg(
  stk::ParallelMachine  parallel_machine,
  int                   argc,
  char **               argv)
{
  const int root_rank = 0;  
  pack_strings_on_root(argc, argv, parallel_machine, root_rank);
  broadcast_strings(argc, parallel_machine, root_rank);
  split_strings();

}

void BroadcastArg::pack_strings_on_root(int argc, char** argv, ParallelMachine parallel_machine, int root_rank)
{
  int rank = stk::parallel_machine_rank(parallel_machine);

  if (rank == root_rank) {
    for (int i = 0; i < argc; ++i) {
      m_string_storage += argv[i];
      m_string_storage += '\0';
    }
  }  
}

void BroadcastArg::broadcast_strings(int argc, ParallelMachine parallel_machine, int root_rank)
{
#ifdef STK_HAS_MPI
  int rank = stk::parallel_machine_rank(parallel_machine);

  int lengths_buffer[2];
  if (rank == root_rank) {
    lengths_buffer[0] = argc;
    lengths_buffer[1] = m_string_storage.size();
  }
  MPI_Bcast(lengths_buffer, 2, MPI_INT, root_rank, parallel_machine);
  m_argc = lengths_buffer[0];
  int buffer_length = lengths_buffer[1];

  if (rank != root_rank) {
    m_string_storage.resize(buffer_length);
  }
  MPI_Bcast(&(m_string_storage[0]), buffer_length, MPI_BYTE, root_rank, parallel_machine);
#endif
}

void BroadcastArg::split_strings()
{
  if(m_argc > 0) {
      m_argv_storage.resize(m_argc);
      int idx = 0;
      for (int i = 0; i < m_argc; ++i) {
          m_argv_storage[i] = &(m_string_storage[idx]);
          while (m_string_storage[idx] != '\0') {
              ++idx;
          }
          ++idx;
      }
  }  

  m_argv = m_argv_storage.data();
}

} // namespace stk
