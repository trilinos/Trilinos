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

#ifndef stk_util_parallel_ReceiveCounter_hpp
#define stk_util_parallel_ReceiveCounter_hpp

#include <cstddef>
#include <vector>

#include "stk_util/parallel/Parallel.hpp"   // for MPI

namespace stk {

// given the number of values each process is going to send to each other process,
// computes the number of other processes each process will receive from.
// Does this in a non-blocking manner
class ReceiveCounter
{
  public:
    explicit ReceiveCounter(MPI_Comm comm) :
      m_comm(comm)
    {}

    ReceiveCounter(const ReceiveCounter& rhs) = delete;

    ReceiveCounter& operator=(const ReceiveCounter& rhs) = delete;

    void start_receive_count(const std::vector<int> &send_counts);

    bool is_complete();

    int get_receive_count();

  private:
    MPI_Comm m_comm;
    std::vector<int> m_sendCount;
    bool m_recvFinished = true;
    int m_nrecv = -1;
    MPI_Request m_recvReq;
};

template <typename T>
constexpr size_t max_elements_per_message(size_t maxMsgSizeInBytes)
{
  return maxMsgSizeInBytes / sizeof(T);
}

template <typename T>
std::vector<int>  get_send_counts(std::vector< std::vector<T> > sendLists, size_t maxSendSizeInBytes)
{
  std::vector<int> sendCounts(sendLists.size());
  for (unsigned int i=0; i < sendLists.size(); ++i)
  {
    //size_t bytesToSend = sendLists[i].size() * sizeof(T);
    //size_t elementsPerMessage = max_
    int numSends = sendLists[i].size() / max_elements_per_message<T>(maxSendSizeInBytes);
    numSends    += sendLists[i].size() %  max_elements_per_message<T>(maxSendSizeInBytes) == 0 ? 0 : 1;
    sendCounts[i] = numSends;
  }

  return sendCounts;
}

}  // namespace
#endif
