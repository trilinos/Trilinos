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

#include "stk_util/stk_config.h"               // for STK_HAS_MPI
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer, CommBufferAlign
#include <algorithm>                           // for copy, max
#include <iostream>                            // for operator<<, basic_ostream, ostringstream
#include <memory>                              // for allocator_traits<>::value_type
#include <stdexcept>                           // for runtime_error, logic_error, range_error
#include <string>                              // for char_traits, string
#include <vector>                              // for vector


namespace stk {

//-----------------------------------------------------------------------

#if STK_MIN_COUPLING_VERSION < 6
namespace {

inline
size_t align_quad( size_t n )
{
  enum { Size = 4 * sizeof(int) };
  return n + CommBufferAlign<Size>::align(n);
}

}

#endif



void CommSparse::reset_buffers()
{
  if (m_exchanger) {
    for (int p=0; p < m_size; ++p)
    {
      m_exchanger->get_send_buf(p).reset();
      m_exchanger->get_recv_buf(p).reset();
    }
  } else {
    m_null_comm_send_buffer.reset();
    m_null_comm_recv_buffer.reset();
  }

  m_num_recvs = DataExchangeUnknownPatternNonBlocking::Unknown;
}

#if STK_MIN_COUPLING_VERSION < 6
void CommSparse::allocate_data(std::vector<CommBuffer>& bufs, std::vector<unsigned char>& data)
{
  size_t n_size = 0;
  for ( size_t i = 0 ; i < bufs.size() ; ++i ) {
    n_size += align_quad( bufs[i].size() );
  }

  // Allocate space for buffers

  data.reserve(n_size);
  unsigned char * p_data = data.data();

  for ( unsigned i = 0 ; i < bufs.size() ; ++i ) {
    CommBuffer & b = bufs[i] ;
    size_t sz = b.size();
    b.set_buffer_ptrs(p_data, p_data, p_data + sz);
    p_data += align_quad( sz );
  }
}
#endif

bool CommSparse::allocate_buffers()
{
  if (m_exchanger) {
    m_exchanger->allocate_send_buffers();
  } else {
    size_t size = m_null_comm_send_buffer.size();
    m_null_comm_storage.resize(size);
    auto* ptr = m_null_comm_storage.data();
    m_null_comm_send_buffer.set_buffer_ptrs(ptr, ptr, ptr + size);
    m_null_comm_recv_buffer.set_buffer_ptrs(ptr, ptr, ptr + size);
  }

  return false;
}

void CommSparse::allocate_buffers(const std::vector<int>& /*send_procs*/, const std::vector<int>& recv_procs)
{
  allocate_buffers();
  m_num_recvs = recv_procs.size();
}

void CommSparse::verify_send_buffers_filled()
{
#ifndef NDEBUG
  for ( int i = 0 ; i < m_size ; ++i ) {
    // Verify the send buffers have been filled
    if ( send_buffer(i).remaining() ) {
      std::ostringstream msg ;
      msg << "stk::CommSparse::communicate LOCAL[" << m_rank << "] ERROR: Send[" << i
          << "] Buffer not filled." ;
      throw std::underflow_error( msg.str() );
    }
  }
#endif
}

bool CommSparse::communicate(bool deallocateSendBuffers)
{
  bool returnValue = false;
#ifdef STK_HAS_MPI
  if (m_exchanger) {
    auto noExtraWork = [](){};
    auto noUnpacker = [](int /*rank*/, stk::CommBuffer& /*buf*/) {};
    communicate_with_extra_work_and_unpacker(noExtraWork, noUnpacker, false);
  }

  for (int i=0; i < parallel_size(); ++i) {
    if (send_buffer(i).capacity() > 0 || recv_buffer(i).capacity() > 0) {
      returnValue = true;
      break;
    }
  }

  if (deallocateSendBuffers) {
    if (m_exchanger) {
      m_exchanger->deallocate_send_bufs();
    }
  }
#endif
  return returnValue;
}

void CommSparse::communicate_with_extra_work_and_unpacker(
                    const std::function<void()>& workFunctor,
                    const std::function<void(int fromProc, CommBuffer& buf)>& unpackFunctor,
                    bool deallocateSendBuffers)
{
#ifdef STK_HAS_MPI
  if (m_exchanger) {
    verify_send_buffers_filled();
  
    m_exchanger->start_nonblocking(m_num_recvs);
    workFunctor();
    m_exchanger->post_nonblocking_receives();
    m_exchanger->complete_receives(unpackFunctor);
    m_exchanger->complete_sends();
  }

  if (deallocateSendBuffers) {
    if (m_exchanger) {
      m_exchanger->deallocate_send_bufs();
    }
  }
#endif
}

}

