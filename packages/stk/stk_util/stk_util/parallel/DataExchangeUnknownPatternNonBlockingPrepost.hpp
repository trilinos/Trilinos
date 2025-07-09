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


#ifndef stk_util_parallel_DataExchangeUnknownPatternNonBlockingPrepost_hpp
#define stk_util_parallel_DataExchangeUnknownPatternNonBlockingPrepost_hpp

#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)

#include "stk_util/parallel/ReceiveSizeCounter.hpp"
#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include "stk_util/parallel/DataExchangeKnownPatternNonBlocking.hpp"

namespace stk {

// General data exchange where communication pattern is not known a priori
// Does a collective operation first to figure out how much data will be received,
// then posts MPI_Irecvs and MPI_Isends
class DataExchangeUnknownPatternNonBlockingPrepost
{
  public:
    static constexpr int Unknown = -1;
    static constexpr size_t MAX_MESSAGE_SIZE = std::numeric_limits<int>::max();

    DataExchangeUnknownPatternNonBlockingPrepost(MPI_Comm comm, [[maybe_unused]] int /*tag_hint*/=11173) :
      m_recvSizeCounter(comm),
      m_exchanger(comm)
    {
    }

    ~DataExchangeUnknownPatternNonBlockingPrepost() {}

    MPI_Comm get_comm() const { return m_exchanger.get_comm(); }


    template <typename T>
    void start_nonblocking(std::vector< std::vector<T> > &sendLists,
                           std::vector< std::vector<T> > &recvLists,
                           [[maybe_unused]] int /*numRecvsExpected*/=Unknown)
    {
      STK_ThrowRequireMsg(!m_areSendsInProgress, "must complete sends before starting new round of communication");      
      STK_ThrowRequireMsg(!m_areRecvsInProgress, "must complete receives before starting new round of communication");      
      
      m_recvSizeCounter.start(sendLists, true);
      const auto& recvSizes = m_recvSizeCounter.get_receive_sizes();
      for (size_t i=0; i < recvLists.size(); ++i)
      {
        recvLists[i].resize(recvSizes[i]);
      }
      
      m_exchanger.start_nonblocking(sendLists, recvLists);
      m_areRecvsInProgress = true;
      m_areSendsInProgress = true;
    }

    template <typename T>
    void post_nonblocking_receives(std::vector< std::vector<T> > & /*recvLists*/) {}

    // Wait for the receives to finish, calling func on each receive
    // buffer as it becomes ready.
    // func must be callable as func(int rank, std::vector<T> recv_buf)
    // (although taking the recv_buf by reference is usually better)
    template <typename T, typename Tfunc>
    void complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func)
    {
      m_exchanger.complete_receives(recvLists, func);
      m_areRecvsInProgress = false;
    }

    // wait for sends to complete, after which time the caller can overwrite the send
    // buffers
    void complete_sends()
    {
      m_exchanger.complete_sends();
      m_areSendsInProgress = false;
    }

    bool are_sends_in_progress() const { return m_areSendsInProgress; }

    bool are_recvs_in_progress() const { return m_areRecvsInProgress; }

  private:

    ReceiveSizeCounter m_recvSizeCounter;
    stk::DataExchangeKnownPatternNonBlocking m_exchanger;
    bool m_areSendsInProgress = false;
    bool m_areRecvsInProgress = false;
};

}

#endif
#endif