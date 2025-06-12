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

#ifndef stk_util_parallel_DataExchangeUnknownPatternNonBlockingBuffer_hpp
#define stk_util_parallel_DataExchangeUnknownPatternNonBlockingBuffer_hpp

#include <vector>
#include "stk_util/parallel/ManagedBufferBase.hpp"   // for ManagedBufferBase
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"

namespace stk {

// General data exchange where communication pattern is not known a priori, where
// buffers are managed by the class
// See DataExchangeUnknownPatternNonBlocking for the API
template <typename T>
class DataExchangeUnknownPatternNonBlockingBuffer : public ManagedBufferBase<T>
{
  using ManagedBufferBase<T>::m_sendBufs;
  using ManagedBufferBase<T>::m_recvBufs;

  public:
    DataExchangeUnknownPatternNonBlockingBuffer(MPI_Comm comm, int tag_hint=11173) :
      ManagedBufferBase<T>(comm),
      m_exchanger(comm, tag_hint)
    {}

    void start_nonblocking(int numRecvsExpected=DataExchangeUnknownPatternNonBlocking::Unknown)
    {
      m_exchanger.start_nonblocking(m_sendBufs, m_recvBufs, numRecvsExpected);
      update_buffer_flags();
    }

    void post_nonblocking_receives()
    {
      m_exchanger.post_nonblocking_receives(m_recvBufs);
      update_buffer_flags();
    }

    template <typename Tfunc>
    void complete_receives(Tfunc func)
    {
      m_exchanger.complete_receives(m_recvBufs, func);
      update_buffer_flags();
    }

    void complete_sends()
    {
      m_exchanger.complete_sends();
      update_buffer_flags();
    }

    bool are_sends_in_progress() const { return m_exchanger.are_sends_in_progress(); }

    bool are_recvs_in_progress() const { return m_exchanger.are_recvs_in_progress(); }

  private:
    void update_buffer_flags()
    {
      ManagedBufferBase<T>::set_sends_in_progress(are_sends_in_progress());
      ManagedBufferBase<T>::set_recvs_in_progress(are_recvs_in_progress());
    }

    DataExchangeUnknownPatternNonBlocking m_exchanger;
};


class DataExchangeUnknownPatternNonBlockingCommBuffer : public ManagedCommBufferBase
{
  public:
    DataExchangeUnknownPatternNonBlockingCommBuffer(MPI_Comm comm, int tag_hint=11173) :
      ManagedCommBufferBase(comm),
      m_exchanger(comm, tag_hint)
    {}

    void start_nonblocking(int numRecvsExpected=DataExchangeUnknownPatternNonBlocking::Unknown)
    {
      if (!get_send_buffers_allocated())
        throw std::runtime_error("Cannot start sends before allocating buffers (CommBuffer requires 2 pass approach to size and then populate buffer");

      m_exchanger.start_nonblocking(m_sendBufStorage, m_recvBufStorage, numRecvsExpected);
      update_buffer_flags();
    }

    void post_nonblocking_receives()
    {
      m_exchanger.post_nonblocking_receives(m_recvBufStorage);
      update_buffer_flags();
    }

    template <typename Tfunc>
    void complete_receives(Tfunc func)
    {
      auto func_vector = [&](int rank, std::vector<unsigned char>& /*buf_vec*/)
      {
        auto& commBuf = set_recv_buffer_storage(rank);
        func(rank, commBuf);
      };

      m_exchanger.complete_receives(m_recvBufStorage, func_vector);
      update_buffer_flags();
    }

    void complete_sends()
    {
      m_exchanger.complete_sends();
      update_buffer_flags();
    }

    bool are_sends_in_progress() const { return m_exchanger.are_sends_in_progress(); }

    bool are_recvs_in_progress() const { return m_exchanger.are_recvs_in_progress(); }

  private:
    void update_buffer_flags()
    {
      set_sends_in_progress(are_sends_in_progress());
      set_recvs_in_progress(are_recvs_in_progress());
    }

    DataExchangeUnknownPatternNonBlocking m_exchanger;
};

}  // namespace

#endif