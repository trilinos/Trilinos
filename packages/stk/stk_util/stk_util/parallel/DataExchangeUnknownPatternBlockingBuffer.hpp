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

#ifndef stk_util_parallel_DataExchangeUnknownPatternBlockingBuffer_hpp
#define stk_util_parallel_DataExchangeUnknownPatternBlockingBuffer_hpp

#include <vector>
#include "stk_util/parallel/ManagedBufferBase.hpp"   // for ManagedBufferBase
#include "stk_util/parallel/DataExchangeUnknownPatternBlocking.hpp"

namespace stk {

// General data exchange where communication pattern is not known a priori, where
// buffers are managed by the class
template <typename T>
class DataExchangeUnknownPatternBlockingBuffer : public ManagedBufferBase<T>
{
  using ManagedBufferBase<T>::m_sendBufs;
  using ManagedBufferBase<T>::m_recvBufs;

  public:
    DataExchangeUnknownPatternBlockingBuffer(MPI_Comm comm, int tag_hint=11173) :
      ManagedBufferBase<T>(comm),
      m_exchanger(comm, tag_hint)
    {}

    void execute(int numRecvsExpected=DataExchangeUnknownPatternNonBlocking::Unknown)
    {
      m_exchanger.execute(m_sendBufs, m_recvBufs, numRecvsExpected);
    }


  private:
    DataExchangeUnknownPatternBlocking m_exchanger;
};


class DataExchangeUnknownPatternBlockingCommBuffer : public ManagedCommBufferBase
{
  public:
    DataExchangeUnknownPatternBlockingCommBuffer(MPI_Comm comm, int tag_hint=11173) :
      ManagedCommBufferBase(comm),
      m_exchanger(comm, tag_hint)
    {}

    void execute(int numRecvsExpected=DataExchangeUnknownPatternNonBlocking::Unknown)
    {
      clear_recv_bufs();
      m_exchanger.execute(m_sendBufStorage, m_recvBufStorage, numRecvsExpected);
      set_recv_buffer_storage();
      reset_send_commbufs();
    }


  private:
    DataExchangeUnknownPatternBlocking m_exchanger;
};

} // namespace

#endif