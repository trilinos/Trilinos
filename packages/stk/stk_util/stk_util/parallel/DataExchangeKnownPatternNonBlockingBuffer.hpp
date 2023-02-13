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

#ifndef stk_util_parallel_DataExchangeKnownPatternNonBlockingBuffer_hpp
#define stk_util_parallel_DataExchangeKnownPatternNonBlockingBuffer_hpp

#include "stk_util/parallel/Parallel.hpp"   // for MPI

#if defined(STK_HAS_MPI)


#include <vector>
#include "stk_util/parallel/ManagedBufferBase.hpp"   // for ManagedBufferBase
#include "stk_util/parallel/DataExchangeKnownPatternNonBlocking.hpp"

namespace stk {

// General data exchange where communication pattern is known a priori,
// and the object manages the send and receive buffer memory.
// See DataExchangeUnknownPatternNonBlocking for the API
template <typename T>
class DataExchangeKnownPatternNonBlockingBuffer : public ManagedBufferBase<T>
{
  using ManagedBufferBase<T>::m_sendBufs;
  using ManagedBufferBase<T>::m_recvBufs;

  public:
    DataExchangeKnownPatternNonBlockingBuffer(MPI_Comm comm, int tag_hint=11176) :
      ManagedBufferBase<T>(comm),
      m_exchanger(comm, tag_hint)
    {}

    MPI_Comm get_comm() const { return m_exchanger.get_comm(); }

    void start_nonblocking()
    {
      m_exchanger.start_nonblocking(m_sendBufs, m_recvBufs);
      this->set_sends_in_progress(true);
      this->set_recvs_in_progress(true);
    }


    template <typename Tfunc>
    void complete_receives(Tfunc func)
    {
      m_exchanger.complete_receives(m_recvBufs, func);
      this->set_recvs_in_progress(false);
    }

    void complete_sends()
    {
      m_exchanger.complete_sends();
      this->set_sends_in_progress(false);
    }

  private:
    DataExchangeKnownPatternNonBlocking m_exchanger;
};


class DataExchangeKnownPatternNonBlockingCommBuffer : public ManagedCommBufferBase
{
  public:
    DataExchangeKnownPatternNonBlockingCommBuffer(MPI_Comm comm, int tag_hint=11176);

    MPI_Comm get_comm() const { return m_exchanger.get_comm(); }

    // sets the size (in bytes) that are expected to be received from the given rank
    void set_recv_buffer_size(int rank, size_t bytes);

    void allocate_recv_buffers();

    void start_nonblocking();

    template <typename Tfunc>
    void complete_receives(Tfunc func)
    {
      auto func_inner = [&](int rank, std::vector<unsigned char>&)
      {
        func(rank, set_recv_buffer_storage(rank));
      };

      m_exchanger.complete_receives(m_recvBufStorage, func_inner);
      set_recvs_in_progress(false);
      m_recvBufsAllocated = false;
    }

    void complete_sends();

  private:
    DataExchangeKnownPatternNonBlocking m_exchanger;
    bool m_recvBufsAllocated;
};


}  // namespace

#endif
#endif