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


#ifndef stk_util_parallel_DataExchangeUnknownPatternNonBlocking_hpp
#define stk_util_parallel_DataExchangeUnknownPatternNonBlocking_hpp

#include <stddef.h>
#include <vector>
#include <deque>
#include <cassert>
#include <iostream>
#include <set>

#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)

#include "stk_util/parallel/ReceiveCounter.hpp"
#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg

namespace stk {

// General data exchange where communication pattern is not known a priori
class DataExchangeUnknownPatternNonBlocking
{
  public:
    static constexpr int Unknown = -1;

    DataExchangeUnknownPatternNonBlocking(MPI_Comm comm, int tag_hint=11173) :
      m_comm(comm),
      m_tagHint(tag_hint),
      m_tag(stk::get_mpi_tag_manager().get_tag(comm, tag_hint)),
      m_recvCounter(comm),
      m_recvcount(0),
      m_numRecvsExpected(Unknown)
    {
      // the MPI standard is a little fuzzy on whether MPI_Requests need
      // to have stable memory addresses.  Reserve enough space to make
      // sure push_back/emplace_back never have to reallocate
      int commSize;
      MPI_Comm_size(comm, &commSize);
      m_sendReqs.reserve(commSize);
      m_recvReqs.reserve(commSize);
    }

    ~DataExchangeUnknownPatternNonBlocking()
    {
      if (m_areSendsInProgress) {
        MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
      }
    }

    MPI_Comm get_comm() const { return m_comm; }


    template <typename T>
    void start_nonblocking(std::vector< std::vector<T> > &sendLists,
                           std::vector< std::vector<T> > &recvLists,
                           int numRecvsExpected=Unknown);

    // this function guarantees that all the receives have been posted
    // It does *not* guarantee the receives have been waited on (ie. that
    // the data is actually in recv_lists)
    template <typename T>
    void post_nonblocking_receives(std::vector< std::vector<T> > &recvLists)
    {
      post_nonblocking_receives_impl(recvLists);
      m_areRecvsInProgress = true;
    }

    // Wait for the receives to finish, calling func on each receive
    // buffer as it becomes ready.
    // func must be callable as func(int rank, std::vector<T> recv_buf)
    // (although taking the recv_buf by reference is usually better)
    template <typename T, typename Tfunc>
    void complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func);

    // wait for sends to complete, after which time the caller can overwrite the send
    // buffers
    void complete_sends()
    {
      if (!m_areSendsInProgress) {
        return;
      }

      MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
      m_areSendsInProgress = false;
    }

    bool are_sends_in_progress() const { return m_areSendsInProgress; }

    bool are_recvs_in_progress() const { return m_areRecvsInProgress; }

  private:
    // maps an index in m_recvReqs to MPI rank
    int get_receive_request_rank(const int idx) { return m_recvRankMap.at(idx); }

    // maps an index in m_sendReqs to MPI rank
    int get_sends_request_rank(const int idx) { return m_sendRankMap.at(idx); }

    void reset();

    template <typename T>
    void start_sends(std::vector< std::vector<T> >& sendLists);

    // do a probe loop until the number of iterations since a probe has been received is equal to
    // nmisses, or all messages have been received.
    // By default, the first criteria is disabled, and the probe loop continues untill all messages have bee
    // received
    template <typename T>
    void post_nonblocking_receives_impl(std::vector< std::vector<T> > &recvLists, int nmisses=-1);

    void yield();

    MPI_Comm m_comm;
    int m_tagHint;
    MPITag m_tag;

    ReceiveCounter m_recvCounter;
    std::vector<int> m_recvRankMap;
    std::vector<int> m_sendRankMap;
    std::vector<MPI_Request> m_recvReqs;
    std::vector<MPI_Request> m_sendReqs;
    bool m_areSendsInProgress = false;
    bool m_areRecvsInProgress = false;
    int m_recvcount;
    int m_numRecvsExpected;
};


template <typename T>
void DataExchangeUnknownPatternNonBlocking::start_sends(std::vector< std::vector<T> >& sendLists)
{
    int commSize;
    MPI_Comm_size(m_comm, &commSize);
    STK_ThrowRequireMsg(sendLists.size() == static_cast<size_t>(commSize), "send list must have length equal to the number of procs");

    for (unsigned int i=0; i < sendLists.size(); ++i) {
      if (sendLists[i].size() > 0) {
          m_sendReqs.emplace_back();
          MPI_Isend(sendLists[i].data(), sendLists[i].size()*sizeof(T), MPI_BYTE, i, m_tag, m_comm, &(m_sendReqs.back()));
          m_sendRankMap.push_back(i);
      }
    }

    m_areSendsInProgress = true;
}


template <typename T>
void DataExchangeUnknownPatternNonBlocking::start_nonblocking(std::vector< std::vector<T> > &sendLists,
                                                              std::vector< std::vector<T> > &recvLists,
                                                              int numRecvsExpected)
{
  STK_ThrowRequireMsg(sendLists.size() == recvLists.size(), "send and receive lists must be same size");
  STK_ThrowRequireMsg(numRecvsExpected >= -1, "num_recvs_expected must be a positive value or -1");
  reset();

  m_numRecvsExpected = numRecvsExpected;
  if (m_numRecvsExpected == Unknown) {
    m_recvCounter.start_receive_count(get_send_counts(sendLists));
  }
  start_sends(sendLists);

  post_nonblocking_receives_impl(recvLists, 100);
  m_areRecvsInProgress = true;
}


template <typename T>
void DataExchangeUnknownPatternNonBlocking::post_nonblocking_receives_impl(std::vector< std::vector<T> > &recvLists, int nmisses)
{
  int nItersSinceReceive = 0;

  while (true)
  {
    int flag = false;
    MPI_Status status;
#ifdef STK_IMPROBE_WORKING
    MPI_Message message;
    MPI_Improbe(MPI_ANY_SOURCE, m_tag, m_comm, &flag, &message, &status);
#else
    MPI_Iprobe(MPI_ANY_SOURCE, m_tag, m_comm, &flag, &status);
#endif

    if (flag) {
      int count = 0;
      int src = status.MPI_SOURCE;
      MPI_Get_count(&status, MPI_BYTE, &count);

      m_recvReqs.emplace_back();
      m_recvRankMap.push_back(src);
      assert(count % sizeof(T) == 0);
      recvLists[src].resize(count / sizeof(T));
#ifdef STK_IMPROBE_WORKING
      MPI_Imrecv(recvLists[src].data(), count, MPI_BYTE, &message, &(m_recvReqs.back()));
#else
      MPI_Irecv(recvLists[src].data(), count, MPI_BYTE, src, m_tag, m_comm, &(m_recvReqs.back()));
#endif

      m_recvcount++;
      nItersSinceReceive = 0;
    }

    if (m_numRecvsExpected == Unknown && m_recvCounter.is_complete()) {
      m_numRecvsExpected = m_recvCounter.get_receive_count();
      assert(m_recvReqs.size() <= static_cast<size_t>(m_numRecvsExpected));
    }

    if (m_recvcount == m_numRecvsExpected) {
      break;
    }

    if (nItersSinceReceive == nmisses) {
      break;
    }

    nItersSinceReceive++;
    yield();
  }
}


template <typename T, typename Tfunc>
void DataExchangeUnknownPatternNonBlocking::complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func)
{
  if (!m_areRecvsInProgress) {
    return;
  }

  for (size_t i=0; i < m_recvReqs.size(); ++i)
  {
    int idx;
    MPI_Waitany(m_recvReqs.size(), m_recvReqs.data(), &idx, MPI_STATUS_IGNORE);
    int rank = get_receive_request_rank(idx);
    func(rank, recvLists[rank]);
  }

  m_areRecvsInProgress = false;
}

} // namespace
#endif
#endif
