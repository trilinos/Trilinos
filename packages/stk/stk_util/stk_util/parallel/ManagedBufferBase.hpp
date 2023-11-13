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

#ifndef stk_util_parallel_ManagedBufferBase_hpp
#define stk_util_parallel_ManagedBufferBase_hpp

#include <stddef.h>
#include <vector>
#include "stk_util/parallel/Parallel.hpp"   // for MPI
#include "stk_util/parallel/CommBuffer.hpp"

namespace stk {

template <typename T>
class ManagedBufferBase
{
  public:
    explicit ManagedBufferBase(MPI_Comm comm) :
    m_comm(comm)
    {
      int commSize;
      MPI_Comm_size(comm, &commSize);
      m_sendBufs.resize(commSize);
      m_recvBufs.resize(commSize);
    }

    virtual ~ManagedBufferBase() {}

    ManagedBufferBase(const ManagedBufferBase&) = delete;

    ManagedBufferBase operator=(const ManagedBufferBase&) = delete;

    using value_type = T;

    MPI_Comm get_comm() const { return m_comm; }

    std::vector<T>& get_send_buf(int rank)
    { 
      if (m_sendsInProgress)
        throw std::runtime_error("cannot access send buffers while sends are in progress (did you forget to complete the sends?)");

      return m_sendBufs[rank];
    }

    const std::vector<T>& get_send_buf(int rank) const
    { 
      return m_sendBufs[rank];
    }

    std::vector<T>& get_recv_buf(int rank)
    { 
      if (m_recvsInProgress)
        throw std::runtime_error("cannot access receive buffers while receives are in progress");

      return m_recvBufs[rank];
    }

    const std::vector<T>& get_recv_buf(int rank) const
    {
      if (m_recvsInProgress)
        throw std::runtime_error("cannot access receive buffers while receives are in progress");

      return m_recvBufs[rank];
    }

    void clear_send_bufs()
    {
      for (size_t i=0; i < m_sendBufs.size(); ++i) {
        get_send_buf(i).clear();
      }
    }

    void clear_recv_bufs()
    {
      for (size_t i=0; i < m_recvBufs.size(); ++i) {
        get_recv_buf(i).clear();
      }
    }

  protected:
    void set_sends_in_progress(bool flag) { m_sendsInProgress = flag; }

    void set_recvs_in_progress(bool flag) { m_recvsInProgress = flag; }

    std::vector< std::vector<T> > m_sendBufs;
    std::vector< std::vector<T> > m_recvBufs;
    bool m_sendsInProgress = false;
    bool m_recvsInProgress = false;
    MPI_Comm m_comm;
};

class ManagedCommBufferBase
{
  public:
    explicit ManagedCommBufferBase(MPI_Comm comm) :
    m_comm(comm)
    {
      int commSize = parallel_machine_size(comm);
      m_sendBufs.resize(commSize);
      m_sendBufStorage.resize(commSize);
      m_recvBufs.resize(commSize);
      m_recvBufStorage.resize(commSize);
    }

    ManagedCommBufferBase(const ManagedCommBufferBase&) = delete;

    ManagedCommBufferBase operator=(const ManagedCommBufferBase&) = delete;

    virtual ~ManagedCommBufferBase() {}

    MPI_Comm get_comm() const { return m_comm; }

    stk::CommBuffer& get_send_buf(int rank)
    {
      if (m_sendBuffersDeallocated)
        throw std::runtime_error("send buffers deallocated!");
      if (m_sendsInProgress)
        throw std::runtime_error("cannot access send buffers while sends are in progress (did you forget to complete the sends?)");

      return m_sendBufs[rank];
    }

    const stk::CommBuffer& get_send_buf(int rank) const
    { 
      if (m_sendBuffersDeallocated)
        throw std::runtime_error("send buffers deallocated!");
      return m_sendBufs[rank];
    }

    stk::CommBuffer& get_recv_buf(int rank)
    { 
      if (m_recvsInProgress)
        throw std::runtime_error("cannot access receive buffers while receives are in progress");

      return m_recvBufs[rank]; 
    }

    const stk::CommBuffer& get_recv_buf(int rank) const 
    {
      if (m_recvsInProgress)
        throw std::runtime_error("cannot access receive buffers while receives are in progress");

      return m_recvBufs[rank];
    }

    void allocate_send_buffers();

    void clear_send_bufs();

    void clear_recv_bufs();

    void deallocate_send_bufs();

  protected:

    void set_recv_buffer_storage();

    CommBuffer& set_recv_buffer_storage(int rank);

    void reset_send_commbufs();

    void set_sends_in_progress(bool flag) { m_sendsInProgress = flag; }

    void set_recvs_in_progress(bool flag) { m_recvsInProgress = flag; }

    bool get_send_buffers_allocated() const;

    std::vector< stk::CommBuffer > m_sendBufs;
    std::vector< stk::CommBuffer > m_recvBufs;
    std::vector< std::vector<unsigned char> > m_sendBufStorage;
    std::vector< std::vector<unsigned char> > m_recvBufStorage;
    bool m_sendsInProgress = false;
    bool m_recvsInProgress = false;
    bool m_sendBuffersAllocated = false;
    bool m_sendBuffersDeallocated = false;
    MPI_Comm m_comm;
};

} // namespace


#endif
