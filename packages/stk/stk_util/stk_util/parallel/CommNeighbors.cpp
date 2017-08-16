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

#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>

namespace stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

void CommNeighbors::rank_error( const char * method , int p ) const
{
  std::ostringstream os ;
  os << "stk::CommNeighbors::" << method
     << "(" << p << ") ERROR: Not in [0:" << m_size << ")" ;
  throw std::range_error( os.str() );
}

//----------------------------------------------------------------------

CommNeighbors::~CommNeighbors()
{
  m_comm = parallel_machine_null();
  m_size = 0 ;
  m_rank = 0 ;
}

CommNeighbors::CommNeighbors( ParallelMachine comm, const std::vector<int>& neighbor_procs)
  : m_comm( comm ),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(neighbor_procs),
    m_recv_procs(neighbor_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
  std::sort(m_neighbor_procs.begin(), m_neighbor_procs.end());
}

CommNeighbors::CommNeighbors( ParallelMachine comm, const std::vector<int>& send_procs, const std::vector<int>& recv_procs)
  : m_comm( comm ),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(send_procs),
    m_recv_procs(recv_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
  stk::util::sort_and_unique(m_send_procs);
  stk::util::sort_and_unique(m_recv_procs);
}

void CommNeighbors::communicate()
{
  const int mpitag = 10101, mpitag2 = 10102;
  int maxRecvProcs = m_recv_procs.size();
  int maxSendProcs = m_send_procs.size();
  int max = std::max(maxSendProcs, maxRecvProcs);
  std::vector<int> recv_sizes(maxRecvProcs, 0);
  std::vector<MPI_Request> requests(max, MPI_REQUEST_NULL);
  std::vector<MPI_Request> requests2(max, MPI_REQUEST_NULL);
  std::vector<MPI_Request> requests3(max, MPI_REQUEST_NULL);
  std::vector<MPI_Status> statuses(max);
  for(int i=0; i<maxRecvProcs; ++i) {
      int p = m_recv_procs[i];
      MPI_Irecv(&recv_sizes[i], 1, MPI_INT, p, mpitag, m_comm, &requests[i]);
  }

  int numSends = 0;
  for(int p : m_send_procs) {
      int send_size = m_send[p].size_in_bytes();
      MPI_Ssend(&send_size, 1, MPI_INT, p, mpitag, m_comm);
      if (send_size > 0) {
          MPI_Issend(m_send[p].raw_buffer(), m_send[p].size_in_bytes(), MPI_BYTE, p, mpitag2, m_comm, &requests2[numSends++]);
      }
  }

  MPI_Status status;
  int numRecvProcs = 0;
  for(int i=0; i<maxRecvProcs; ++i) {
      int idx = 0;
      MPI_Waitany(maxRecvProcs, &requests[0], &idx, &status);
      int p = status.MPI_SOURCE;
      if (recv_sizes[idx] > 0) {
          m_recv[p].resize(recv_sizes[idx]);
          MPI_Irecv(m_recv[p].raw_buffer(), m_recv[p].size_in_bytes(), MPI_BYTE, p, mpitag2, m_comm, &requests3[numRecvProcs]);
          numRecvProcs++;
      }   
  }

  if (numRecvProcs > 0) {
      MPI_Waitall(numRecvProcs, requests3.data(), statuses.data());
  }
  if (numSends > 0) {
      MPI_Waitall(numSends, requests2.data(), statuses.data());
  }
}

#endif

}

