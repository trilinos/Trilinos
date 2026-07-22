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

#include "stk_util/parallel/DataExchangeKnownPatternUserDataNonBlocking.hpp"
#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)

#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include "ReceiveCounter.hpp"

namespace stk {

class DataExchangeKnownPatternNonBlocking
{
  public:
    static constexpr size_t MAX_MESSAGE_SIZE = std::numeric_limits<int>::max();
  
    DataExchangeKnownPatternNonBlocking(MPI_Comm comm, int tag_hint=11176);

    MPI_Comm get_comm() const { return m_exchanger.get_comm(); }

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
    
    struct RecvCounts
    {
      int numPosted = 0;
      int numCompleted = 0;
    };
    
    DataExchangeKnownPatternUserDataNonBlocking m_exchanger;
    std::vector<int> m_sendRanks;
    std::vector<int> m_recvRanks;
    std::vector<PointerAndSize> m_sendData;
    std::vector<PointerAndSize> m_recvData;
};


template <typename T>
void DataExchangeKnownPatternNonBlocking::start_nonblocking(std::vector< std::vector<T> > &sendLists,
                                                            std::vector< std::vector<T> > &recvLists)
{
  STK_ThrowRequireMsg(sendLists.size() == size_t(stk::parallel_machine_size(get_comm())), "sendLists must have length comm_size()");  
  STK_ThrowRequireMsg(recvLists.size() == size_t(stk::parallel_machine_size(get_comm())), "recvLists must have length comm_size()");
  
  m_sendRanks.resize(0);
  m_sendData.resize(0);
  m_recvRanks.resize(0);
  m_recvData.resize(0);
  for (size_t i=0; i < sendLists.size(); ++i)
  {
    if (sendLists[i].size() > 0)
    {
      m_sendRanks.push_back(i);
      m_sendData.emplace_back(reinterpret_cast<unsigned char*>(sendLists[i].data()),
                              sendLists[i].size() * sizeof(T));
    }
    
    if (recvLists[i].size() > 0)
    {
      m_recvRanks.push_back(i);
      m_recvData.emplace_back(reinterpret_cast<unsigned char*>(recvLists[i].data()),
                              recvLists[i].size() * sizeof(T));
    }    
  }
  
  m_exchanger.start_nonblocking(m_sendData, m_sendRanks, m_recvData, m_recvRanks);
}


template <typename T, typename Tfunc>
void DataExchangeKnownPatternNonBlocking::complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func)
{
  STK_ThrowRequireMsg(recvLists.size() == size_t(stk::parallel_machine_size(get_comm())), "recvLists must have length comm_size()");
  
  auto funcWrapper = [&](int rank, const PointerAndSize& /*buf*/)
  {
    func(rank, recvLists[rank]);
  };
  
  m_exchanger.complete_receives(m_recvData, m_recvRanks, funcWrapper);

}

}  // namespace

#endif
#endif
