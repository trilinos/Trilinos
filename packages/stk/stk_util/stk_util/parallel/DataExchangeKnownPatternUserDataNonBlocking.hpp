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

#ifndef stk_util_parallel_DataExchangeKnownPatternUserDataNonBlocking_hpp
#define stk_util_parallel_DataExchangeKnownPatternUserDataNonBlocking_hpp

#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)

#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg

namespace stk {
  
struct PointerAndSize
{
  PointerAndSize(unsigned char* ptr_, size_t size_) :
    ptr(ptr_),
    size(size_)
  {}
  
  PointerAndSize() = default;
  
  unsigned char* ptr = nullptr;
  size_t size = 0;  
};

class DataExchangeKnownPatternUserDataNonBlocking
{
  public:
    DataExchangeKnownPatternUserDataNonBlocking(MPI_Comm comm, int tag_hint=11176):
      m_comm(comm),
      m_tag(stk::get_mpi_tag_manager().get_tag(comm, tag_hint))
    {}
    
    MPI_Comm get_comm() const { return m_comm; }
    
    ~DataExchangeKnownPatternUserDataNonBlocking();
    
    // Start non-blocking sends and receives.  The recvData must be sized
    // correctly to indicate how many items are expected to be received from
    // each process
    void start_nonblocking(const std::vector<PointerAndSize>& sendData, const std::vector<int>& sendRanks,
                           std::vector<PointerAndSize>& recvData, const std::vector<int>& recvRanks);

    // Wait for the receives to finish, calling func on each receive
    // buffer as it becomes ready.
    // func must be callable as func(int rank, PointerAndSize recv_buf)
    // (although taking the recv_buf by reference is usually better)
    template <typename Tfunc>
    void complete_receives(std::vector<PointerAndSize>& recvData, const std::vector<int>& recvRanks, Tfunc func);

    // Complete sends.  After this function is called, the send buffers can be
    // safely modified.  This function is only needed of more than one round of
    // communication will be done with the same object (the destructor will
    // complete the sends if needed).
    void complete_sends();
  
  private:  
  
    static constexpr size_t MAX_MESSAGE_SIZE = std::numeric_limits<int>::max();
    
    int computeNumSegments(size_t totalSize, size_t maxSegmentSize)
    {
      int numSegments = totalSize / maxSegmentSize;
      numSegments += (totalSize > numSegments * maxSegmentSize) ? 1 : 0;      
      return numSegments;
    }

    void start_sends(const std::vector<PointerAndSize>& sendData, 
                     const std::vector<int>& sendRanks);
        
    void start_recvs(std::vector<PointerAndSize> &recvData,
                     const std::vector<int>& recvRanks);
    

    struct RecvCounts
    {
      int numPosted = 0;
      int numCompleted = 0;
    };
      
    MPI_Comm m_comm;
    MPITag m_tag;

    std::vector<int> m_recvIndexMap;
    std::vector<RecvCounts> m_recvCounts;
    std::vector<MPI_Request> m_recvReqs;
    std::vector<MPI_Request> m_sendReqs;
    bool m_areSendsInProgress = false;
    bool m_areRecvsInProgress = false;  
}; 


template <typename Tfunc>
void DataExchangeKnownPatternUserDataNonBlocking::complete_receives(std::vector<PointerAndSize> &recvData, const std::vector<int>& recvRanks, Tfunc func)
{
  STK_ThrowRequireMsg(recvData.size() == recvRanks.size(), "recvData and recvRanks must have same length");
  STK_ThrowRequireMsg(m_areRecvsInProgress, "Receives must have been started before they can be completed");

  std::vector<int> rankRecvCounts(recvRanks.size(), 0);
  for (size_t i=0; i < m_recvReqs.size(); ++i)
  {
    int idxInReqs;
    MPI_Waitany(m_recvReqs.size(), m_recvReqs.data(), &idxInReqs, MPI_STATUS_IGNORE);

    int idxInRecvRanks = m_recvIndexMap[idxInReqs];
    m_recvCounts[idxInRecvRanks].numCompleted++;
    if (m_recvCounts[idxInRecvRanks].numCompleted == m_recvCounts[idxInRecvRanks].numPosted)
    {
      func(recvRanks[idxInRecvRanks], recvData[idxInRecvRanks]);
    }
  }

  m_areRecvsInProgress = false;
}


}

#endif
#endif