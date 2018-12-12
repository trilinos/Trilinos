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
#include <cstring>
#include <vector>
#include <algorithm>

#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>

namespace stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#if MPI_VERSION >= 3
#define STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#ifdef OMPI_MAJOR_VERSION
//OpenMPI 3.1.x seems to have a bug in the MPI_Neighbor* functions.
#if OMPI_MAJOR_VERSION == 3 && OMPI_MINOR_VERSION == 1
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif
//OpenMPI 2.x.y doesn't seem to support MPI_Neighbor* functions either...
#if OMPI_MAJOR_VERSION == 2
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#endif

//the MPI_Neighbor functions seem to be unacceptably slow with intel mpi
#ifdef I_MPI_VERSION
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

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
  if (m_created_dist_graph) {
    MPI_Comm_free(&m_comm);
  }
  m_comm = stk::parallel_machine_null();
  m_created_dist_graph = false;
  m_size = 0 ;
  m_rank = 0 ;
}

stk::ParallelMachine CommNeighbors::setup_neighbor_comm(stk::ParallelMachine fullComm,
                                                        const std::vector<int>& sendProcs,
                                                        const std::vector<int>& recvProcs)
{
  stk::ParallelMachine neighborComm = fullComm;
#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM
  MPI_Info info;
  MPI_Info_create(&info);
  int reorder = 0;
  const int* weights = (int*)MPI_UNWEIGHTED;
  MPI_Dist_graph_create_adjacent(fullComm,
              recvProcs.size(), recvProcs.data(), weights,
              sendProcs.size(), sendProcs.data(), weights,
              info, reorder, &neighborComm);
  m_created_dist_graph = true;
  MPI_Info_free(&info);
#endif
  return neighborComm;
}

CommNeighbors::CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& neighbor_procs)
  : m_comm( comm ),
    m_created_dist_graph( false ),
    m_size( stk::parallel_machine_size( comm ) ),
    m_rank( stk::parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(neighbor_procs),
    m_recv_procs(neighbor_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
#ifdef OMPI_MAJOR_VERSION
#if OMPI_MAJOR_VERSION < 2
//if open-mpi version 1.10, MPI_Neighbor_* functions can't handle
//empty send/recv lists.
  if (neighbor_procs.empty()) {
    m_send_procs.push_back(m_rank);
    m_recv_procs.push_back(m_rank);
  }
#endif
#endif
  stk::util::sort_and_unique(m_send_procs);
  stk::util::sort_and_unique(m_recv_procs);

  if (comm != MPI_COMM_NULL && m_size > 1) {
    m_comm = setup_neighbor_comm(comm, m_send_procs, m_recv_procs);
  }
}

CommNeighbors::CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& send_procs, const std::vector<int>& recv_procs)
  : m_comm( comm ),
    m_created_dist_graph( false ),
    m_size( stk::parallel_machine_size( comm ) ),
    m_rank( stk::parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(send_procs),
    m_recv_procs(recv_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);

  //combine send-procs and recv-procs to make them be equal. This shouldn't be
  //necessary, but if one of the lists is empty it triggers a crash in
  //Open-MPI 1.10. When we upgrade to a newer version of Open-MPI then
  //presumably we won't need to do this.
  std::vector<int> symmNeighbors = m_send_procs;
  symmNeighbors.insert(symmNeighbors.end(), m_recv_procs.begin(), m_recv_procs.end());
  stk::util::sort_and_unique(symmNeighbors);
#ifdef OMPI_MAJOR_VERSION
#if OMPI_MAJOR_VERSION < 2
  if (symmNeighbors.empty()) {
    symmNeighbors.push_back(m_rank);
  }
#endif
#endif
  m_send_procs = symmNeighbors;
  m_recv_procs = symmNeighbors;

  if (comm != MPI_COMM_NULL && m_size > 1) {
    m_comm = setup_neighbor_comm(comm, m_send_procs, m_recv_procs);
  }
}

void setup_buffers(const std::vector<int>& procs,
                   const std::vector<CommBufferV>& data,
                   std::vector<int>& counts,
                   std::vector<int>& displs,
                   std::vector<unsigned char>& buf, int m_rank)
{
  counts.resize(procs.size());
  displs.resize(procs.size());

  int totalBytes = 0;
  for(size_t i=0; i<procs.size(); ++i) {
    int p = procs[i];
    counts[i] = data[p].size_in_bytes();
    displs[i] = totalBytes;
    totalBytes += counts[i];
  }

  buf.resize(totalBytes, 0);
  buf.clear();

  for(size_t i=0; i<procs.size(); ++i) {
    int p = procs[i];
    const unsigned char* rawbuf = data[p].raw_buffer();
    size_t len = counts[i];
    buf.insert(buf.end(), rawbuf, rawbuf+len);
  }
}

void store_recvd_data(const std::vector<unsigned char>& recvBuf,
                      const std::vector<int>& recvCounts,
                      const std::vector<int>& recvDispls,
                      const std::vector<int>& recvProcs,
                      std::vector<CommBufferV>& recvBuffers)
{
  for(size_t i=0; i<recvProcs.size(); ++i) {
    int p = recvProcs[i];
    const unsigned char* buf = recvBuf.data() + recvDispls[i];
    int len = recvCounts[i];
    recvBuffers[p].resize(len);
    std::memcpy(recvBuffers[p].raw_buffer(), buf, len);
  }
}

void CommNeighbors::perform_neighbor_communication(MPI_Comm neighborComm,
                                                   const std::vector<unsigned char>& sendBuf,
                                                   const std::vector<int>& sendCounts,
                                                   const std::vector<int>& sendDispls,
                                                         std::vector<unsigned char>& recvBuf,
                                                         std::vector<int>& recvCounts,
                                                         std::vector<int>& recvDispls)
{
  ThrowAssertMsg(sendCounts.size()==m_send_procs.size(), "Error, sendCounts should be same size as m_send_procs.");
  ThrowAssertMsg(sendDispls.size()==m_send_procs.size(), "Error, sendDispls should be same size as m_send_procs.");
  ThrowAssertMsg(recvCounts.size()==m_recv_procs.size(), "Error, recvCounts should be same size as m_recv_procs.");
  ThrowAssertMsg(recvDispls.size()==m_recv_procs.size(), "Error, recvDispls should be same size as m_recv_procs.");

#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM
  const int* sendCountsPtr = sendCounts.size() > 0 ? sendCounts.data() : nullptr;
  const int* recvCountsPtr = recvCounts.size() > 0 ? recvCounts.data() : nullptr;

  MPI_Neighbor_alltoall((void*)sendCountsPtr, 1, MPI_INT,
                        (void*)recvCountsPtr, 1, MPI_INT, neighborComm);
#endif

  int totalRecv = 0;
  for(size_t i=0; i<recvCounts.size(); ++i) {
    recvDispls[i] = totalRecv;
    totalRecv += recvCounts[i];
  }
  recvBuf.resize(totalRecv);

#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM 
  const unsigned char* sendBufPtr = sendBuf.size() > 0 ? sendBuf.data() : nullptr;
  const unsigned char* recvBufPtr = recvBuf.size() > 0 ? recvBuf.data() : nullptr;
  const int* sendDisplsPtr = sendDispls.size() > 0 ? sendDispls.data() : nullptr;
  const int* recvDisplsPtr = recvDispls.size() > 0 ? recvDispls.data() : nullptr;
  MPI_Neighbor_alltoallv(
      (void*)sendBufPtr, sendCountsPtr, sendDisplsPtr, MPI_BYTE,
      (void*)recvBufPtr, recvCountsPtr, recvDisplsPtr, MPI_BYTE, neighborComm);
#endif
}

#ifndef STK_MPI_SUPPORTS_NEIGHBOR_COMM

void old_communicate(MPI_Comm comm,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     const std::vector<CommBufferV>& send_buffers,
                     std::vector<CommBufferV>& recv_buffers)
{
  const int mpitag = 10101, mpitag2 = 10102;
  int maxRecvProcs = recv_procs.size();
  int maxSendProcs = send_procs.size();
  int max = std::max(maxSendProcs, maxRecvProcs);
  std::vector<int> recv_sizes(maxRecvProcs, 0); 
  std::vector<MPI_Request> requests(max, MPI_REQUEST_NULL);
  std::vector<MPI_Request> requests2(max, MPI_REQUEST_NULL);
  std::vector<MPI_Request> requests3(max, MPI_REQUEST_NULL);
  std::vector<MPI_Status> statuses(max);
  for(int i=0; i<maxRecvProcs; ++i) {
      int p = recv_procs[i];
      MPI_Irecv(&recv_sizes[i], 1, MPI_INT, p, mpitag, comm, &requests[i]);
  }

  int numSends = 0;
  for(int p : send_procs) {
      int send_size = send_buffers[p].size_in_bytes();
      MPI_Ssend(&send_size, 1, MPI_INT, p, mpitag, comm);
      if (send_size > 0) {
          MPI_Issend(send_buffers[p].raw_buffer(), send_buffers[p].size_in_bytes(), MPI_BYTE, p, mpitag2, comm, &requests2[numSends++]);
      }   
  }

  MPI_Status status;
  int numRecvProcs = 0;
  for(int i=0; i<maxRecvProcs; ++i) {
      int idx = 0;
      MPI_Waitany(maxRecvProcs, &requests[0], &idx, &status);
      int p = status.MPI_SOURCE;
      if (recv_sizes[idx] > 0) {
          recv_buffers[p].resize(recv_sizes[idx]);
          MPI_Irecv(recv_buffers[p].raw_buffer(), recv_buffers[p].size_in_bytes(), MPI_BYTE, p, mpitag2, comm, &requests3[numRecvProcs]);
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

void CommNeighbors::communicate()
{
  if (m_size == 1) {
    int len = m_send[0].size_in_bytes();
    const unsigned char* buf = m_send[0].raw_buffer();
    m_recv[0].resize(len);
    std::memcpy(m_recv[0].raw_buffer(), buf, len);
    return;
  }

#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM
  std::vector<int> sendCounts, recvCounts(m_recv_procs.size(), 0);
  std::vector<int> sendDispls, recvDispls(m_recv_procs.size(), 0);
  std::vector<unsigned char> sendBuf, recvBuf;

  setup_buffers(m_send_procs, m_send, sendCounts, sendDispls, sendBuf, m_rank);

  perform_neighbor_communication(m_comm, sendBuf, sendCounts, sendDispls,
                                         recvBuf, recvCounts, recvDispls);

  store_recvd_data(recvBuf, recvCounts, recvDispls, m_recv_procs, m_recv);

#else

  old_communicate(m_comm, m_send_procs, m_recv_procs, m_send, m_recv);

#endif
}

#endif

}

