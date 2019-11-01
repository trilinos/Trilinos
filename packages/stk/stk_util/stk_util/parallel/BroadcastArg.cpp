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

#include <stk_util/parallel/BroadcastArg.hpp>
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for copy
#include <string>                       // for basic_string, string
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc


namespace stk {

BroadcastArg::BroadcastArg(
  stk::ParallelMachine  parallel_machine,
  int                   argc,
  char **               argv)
{
    m_argc = 0;
    m_argv = nullptr;


  int rank = stk::parallel_machine_rank(parallel_machine);

  size_t buffer_length = 0;
  char * buffer = nullptr;
  
// Populate m_argc, m_buffer and buffer_length on rank 0 or !STK_HAS_MPI
  if (rank == 0) {
    std::string s;
    for (int i = 0; i < argc; ++i) {
      s += argv[i];
      s += '\0';
    }
    
    buffer_length = s.size();
    if(buffer_length > 0) {
        buffer = new char[buffer_length];
    }
    
    std::copy(s.begin(), s.end(), buffer);
  }

// if STK_HAS_MPI, broadcast m_argc, buffer and buffer_length to processors
#ifdef STK_HAS_MPI

  int lengths_buffer[2];
  if (rank == 0) {
    lengths_buffer[0] = argc;
    lengths_buffer[1] = buffer_length;
    MPI_Bcast(lengths_buffer, 2, MPI_INT, 0, parallel_machine);
  } else {
    MPI_Bcast(lengths_buffer, 2, MPI_INT, 0, parallel_machine);
  }
  m_argc = lengths_buffer[0];
  buffer_length = lengths_buffer[1];
  if(buffer_length > 0) {
      if (rank == 0) {
          MPI_Bcast(buffer, buffer_length, MPI_BYTE, 0, parallel_machine);
      }
      else {
          buffer = new char[buffer_length];
          MPI_Bcast(buffer, buffer_length, MPI_BYTE, 0, parallel_machine); 
      }    
  }
#endif
  if(m_argc > 0) {
      m_argv = new char *[m_argc];
      char *c = buffer;
      for (int i = 0; i < m_argc; ++i) {
          m_argv[i] = c;
          while (*c != '\0') {
              ++c;
          }
          ++c;
      }
  }
}


BroadcastArg::~BroadcastArg()
{
    if(m_argc > 0) {
        delete [] m_argv[0];
    }
    delete [] m_argv;
}

} // namespace stk
