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

#ifndef stk_util_parallel_DataExchangeKnownPatternNonBlocking_hpp
#define stk_util_parallel_DataExchangeKnownPatternNonBlocking_hpp

#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)

#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg

namespace stk {

class DataExchangeKnownPatternNonBlocking
{
  public:
    DataExchangeKnownPatternNonBlocking(MPI_Comm comm, int tag_hint=11176);

    MPI_Comm get_comm() const { return m_comm; }

    ~DataExchangeKnownPatternNonBlocking();

    // Start non-blocking sends and receives.  The recvLists must be sized
    // correctly to indicate how many items are expected to be received from
    // each process
    template <typename T>
    void start_nonblocking(std::vector< std::vector<T> > &sendLists,
                           std::vector< std::vector<T> > &recvLists);

    // Wait for the receives to finish, calling func on each receive
    // buffer as it becomes ready.
    // func must be callable as func(int rank, std::vector<T> recv_buf)
    // (although taking the recv_buf by reference is usually better)
    template <typename T, typename Tfunc>
    void complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func);

    // Complete sends.  After this function is called, the send buffers can be
    // safely modified.  This function is only needed of more than one round of
    // communication will be done with the same object (the destructor will
    // complete the sends if needed).
    void complete_sends();


  private:
    template <typename T>
    void start_recvs(std::vector< std::vector<T> > &recvLists);
    
    template <typename T>
    void start_sends(std::vector< std::vector<T> > &sendLists);

    MPI_Comm m_comm;
    MPITag m_tag;

    std::vector<int> m_recvRankMap;
    std::vector<int> m_sendRankMap;
    std::vector<MPI_Request> m_recvReqs;
    std::vector<MPI_Request> m_sendReqs;
    bool m_areSendsInProgress = false;
    bool m_areRecvsInProgress = false;
};


template <typename T>
void DataExchangeKnownPatternNonBlocking::start_nonblocking(std::vector< std::vector<T> > &sendLists,
                                                            std::vector< std::vector<T> > &recvLists)
{
  start_recvs(recvLists);
  start_sends(sendLists);
}


template <typename T>
void DataExchangeKnownPatternNonBlocking::start_recvs(std::vector< std::vector<T> > &recvLists)
{
  STK_ThrowRequireMsg(recvLists.size() == size_t(stk::parallel_machine_size(m_comm)), "recvLists must have length comm_size()");
  STK_ThrowRequireMsg(!m_areRecvsInProgress,
                  "cannot start new round of communication until the recvs from the previous round are complete");

  m_recvRankMap.resize(0);
  m_recvReqs.resize(0);
  for (size_t rank=0; rank < recvLists.size(); ++rank)
  {
    auto& recvBuf = recvLists[rank];
    if (recvBuf.size() > 0)
    {
      m_recvRankMap.push_back(rank);
      m_recvReqs.emplace_back();

      MPI_Irecv(recvBuf.data(), recvBuf.size()*sizeof(T), MPI_BYTE, rank, 
                m_tag, m_comm, &(m_recvReqs.back()));
    }
  }

  m_areRecvsInProgress = true;
}

template <typename T>
void DataExchangeKnownPatternNonBlocking::start_sends(std::vector< std::vector<T> > &sendLists)
{
  STK_ThrowRequireMsg(sendLists.size() == size_t(stk::parallel_machine_size(m_comm)), "sendLists must have length comm_size()");
  STK_ThrowRequireMsg(!m_areSendsInProgress, 
                  "cannot start new round of communication until the sends from the previous round are complete");
  
  m_sendRankMap.resize(0);
  m_sendReqs.resize(0);
  for (size_t rank=0; rank < sendLists.size(); ++rank)
  {
    auto& sendBuf = sendLists[rank];
    if (sendBuf.size() > 0)
    {
      m_sendRankMap.push_back(rank);
      m_sendReqs.emplace_back();

      MPI_Isend(sendBuf.data(), sendBuf.size()*sizeof(T), MPI_BYTE, rank,
                m_tag, m_comm, &(m_sendReqs.back()));
    }
  }

  m_areSendsInProgress = true;
}

template <typename T, typename Tfunc>
void DataExchangeKnownPatternNonBlocking::complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func)
{
  STK_ThrowRequireMsg(recvLists.size() == size_t(stk::parallel_machine_size(m_comm)), "recvLists must have length comm_size()");
  STK_ThrowRequireMsg(m_areRecvsInProgress, "Receives must have been started before they can be completed");

  for (size_t i=0; i < m_recvReqs.size(); ++i)
  {
    int idx;
    MPI_Waitany(m_recvReqs.size(), m_recvReqs.data(), &idx, MPI_STATUS_IGNORE);

    int rank = m_recvRankMap[idx];
    func(rank, recvLists[rank]);
  }

  m_areRecvsInProgress = false;
}

}  // namespace

#endif
#endif
