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

#include "stk_util/parallel/CommNeighbors.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg
#include "stk_util/util/SortAndUnique.hpp"  // for sort_and_unique
#include <cstddef>                          // for size_t
#include <cstring>                          // for memcpy
#include <iostream>                         // for operator<<, basic_ostream::operator<<, basic_...
#include <memory>                           // for allocator_traits<>::value_type
#include <stdexcept>                        // for range_error
#include <vector>                           // for vector


namespace stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

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
  if (neighbor_procs.empty()) {
    //at least some MPI_Neighbor_* implementations can't handle empty neighbor lists.
    m_send_procs.push_back(m_rank);
    m_recv_procs.push_back(m_rank);
  }
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
  if (symmNeighbors.empty()) {
    symmNeighbors.push_back(m_rank);
  }
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
    const int p = recvProcs[i];
    const unsigned char* buf = recvBuf.data() + recvDispls[i];
    const int len = recvCounts[i];
    recvBuffers[p].resize(len);
    if (len > 0) {
      std::memcpy(recvBuffers[p].raw_buffer(), buf, len);
    }
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
  STK_ThrowAssertMsg(sendCounts.size()==m_send_procs.size(), "Error, sendCounts should be same size as m_send_procs.");
  STK_ThrowAssertMsg(sendDispls.size()==m_send_procs.size(), "Error, sendDispls should be same size as m_send_procs.");
  STK_ThrowAssertMsg(recvCounts.size()==m_recv_procs.size(), "Error, recvCounts should be same size as m_recv_procs.");
  STK_ThrowAssertMsg(recvDispls.size()==m_recv_procs.size(), "Error, recvDispls should be same size as m_recv_procs.");

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
      MPI_Waitany(maxRecvProcs, requests.data(), &idx, &status);
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
    const int len = m_send[0].size_in_bytes();
    const unsigned char* buf = m_send[0].raw_buffer();
    m_recv[0].resize(len);
    if (len > 0) {
      std::memcpy(m_recv[0].raw_buffer(), buf, len);
    }
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

void CommNeighbors::reset_buffers() {
  for(auto&& s : m_send) {
    s.resize(0);
  }
  for(auto&& r : m_recv) {
    r.resize(0);
  }
}

void CommNeighbors::swap_send_recv()
{
  m_send.swap(m_recv);
  m_send_data.swap(m_recv_data);
  m_send_procs.swap(m_recv_procs);
}

#endif

}

