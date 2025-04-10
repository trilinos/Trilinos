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

#ifndef stk_util_parallel_ReceiveSizeCounter_hpp
#define stk_util_parallel_ReceiveSizeCounter_hpp

#include "stk_util/parallel/Parallel.hpp"

#ifdef STK_HAS_MPI

#include <cstddef>
#include <cassert>
#include <vector>
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/MPITagManager.hpp"

namespace stk {

class ReceiveSizeCounter
{
  public:
    using ULL = unsigned long long;
    
    explicit ReceiveSizeCounter(MPI_Comm comm) :
      m_comm(comm),
      m_sendStatus(stk::parallel_machine_size(comm)),
      m_recvCounts(stk::parallel_machine_size(comm)),      
      m_tag(get_mpi_tag_manager().get_tag(comm))
    {}
    
    ~ReceiveSizeCounter()
    {
      MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);      
    }
    
    template <typename T>
    void start(const std::vector<std::vector<T>>& send_bufs, bool wait=false)
    {
      int commSize = stk::parallel_machine_size(m_comm);
      STK_ThrowRequireMsg(send_bufs.size() == size_t(commSize), "send_bufs must have length equal to number of procs");
      STK_ThrowRequireMsg(m_numReceived == m_nrecvs, "must complete previous receive before starting new one");
      
      MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
      m_numReceived = 0;
      m_nrecvs = 0;
      std::fill(m_recvCounts.begin(), m_recvCounts.end(), 0);
      m_sendStatus.resize(0);
      m_sendCounts.resize(0);
      m_tag = get_mpi_tag_manager().get_tag(m_comm);      
      for (const std::vector<T>& buf : send_bufs)
      {
        m_sendStatus.push_back(buf.size() > 0 ? 1 : 0);
        if (buf.size() > 0)
        {
          m_sendCounts.push_back(buf.size());
        }
      }
      
      MPI_Reduce_scatter_block(m_sendStatus.data(), &m_nrecvs, 1, MPI_INT, MPI_SUM, m_comm);
      
      m_recvCountsCompressed.resize(m_nrecvs);
      m_recvReqs.resize(m_nrecvs);
      for (int i=0; i < m_nrecvs; ++i)
      {
        MPI_Irecv(&(m_recvCountsCompressed[i]), 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, m_tag, m_comm, &(m_recvReqs[i]));
      }
      
      m_sendReqs.resize(0);
      int idx = 0;
      for (int i=0; i < commSize; ++i)
      {
        if (m_sendStatus[i] > 0)
        {
          m_sendReqs.emplace_back();
          MPI_Isend(&(m_sendCounts[idx++]), 1, MPI_UNSIGNED_LONG_LONG, i, m_tag, m_comm, &(m_sendReqs.back()));
        }
      }
      
            
      if (wait)
      {
        finish_receives();
      }
    }
    
    bool is_complete();

    const std::vector<ULL>& get_receive_sizes() const;
  
  private:
    void finish_receives();
    
    MPI_Comm m_comm;
    std::vector<int> m_sendStatus;
    std::vector<ULL> m_sendCounts;
    std::vector<ULL> m_recvCounts;
    int m_nrecvs = 0;
    int m_numReceived = 0;
    
    MPITag m_tag;
    std::vector<MPI_Request> m_sendReqs;
    std::vector<MPI_Request> m_recvReqs;
    std::vector<ULL> m_recvCountsCompressed;
};

}

#endif
#endif
