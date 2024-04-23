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

#ifndef stk_util_parallel_CommSparse_hpp
#define stk_util_parallel_CommSparse_hpp

#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/Parallel.hpp"      // for ParallelMachine, parallel_machine_null
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include <stddef.h>
#include <cstddef>                             // for size_t
#include <vector>                              // for vector

//------------------------------------------------------------------------

namespace stk {

/** Given the send sizes determine the receive sizes and also
 * set output vectors of send-procs and recv-procs.
 * Send and receive size arrays are dimensioned to
 * the size of the parallel machine. (i.e., the number
 * of MPI processor ranks.)
 * Output vectors for send-procs and recv-procs will have
 * length num-send-procs and num-recv-procs respectively.
 */
// delete coupling version 5 is deprecated
void comm_recv_procs_and_msg_sizes(ParallelMachine comm,
                     const unsigned * const send_size,
                     unsigned * const recv_size,
                     std::vector<int>& output_send_procs,
                     std::vector<int>& output_recv_procs);

void comm_recv_procs_and_msg_sizes(ParallelMachine comm ,
                                   const std::vector<CommBuffer>& send_bufs ,
                                         std::vector<CommBuffer>& recv_bufs,
                                   std::vector<int>& send_procs,
                                   std::vector<int>& recv_procs);

/** Given send sizes (of length number-of-MPI-processor-ranks) and
 * send-procs and recv-procs (of length number-of-procs-to-send/recv-with),
 * set recv sizes (recv_size array has length number-of-MPI-processor-ranks).
 */

void comm_recv_msg_sizes(ParallelMachine comm ,
                     const unsigned * const send_size ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     unsigned * const recv_size);

void comm_recv_msg_sizes(ParallelMachine comm ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     const std::vector<CommBuffer>& send_bufs,
                     std::vector<CommBuffer>& recv_bufs);

class CommSparse {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for a given processor */
  CommBuffer & send_buffer( int p )
  {
    STK_ThrowAssertMsg(p < m_size,"CommSparse::send_buffer: "<<p<<" out of range [0:"<<m_size<<")");
    stk::util::print_unsupported_version_warning(5, __LINE__, __FILE__);

    if (stk::util::get_common_coupling_version() >= 6) {
      if (m_exchanger) {
        return m_exchanger->get_send_buf(p);
      } else {
        return m_null_comm_send_buffer;
      }
    } else {
      return m_send[p] ;
    }
  }

  const CommBuffer & send_buffer( int p ) const
  {
    STK_ThrowAssertMsg(p < m_size,"CommSparse::send_buffer: "<<p<<" out of range [0:"<<m_size<<")");
    stk::util::print_unsupported_version_warning(5, __LINE__, __FILE__);

    if (stk::util::get_common_coupling_version() >= 6) {
      if (m_exchanger) {
        return m_exchanger->get_send_buf(p);
      } else {
        return m_null_comm_send_buffer;
      }
    } else {
      return m_send[p] ;
    }
  }

  /** Obtain the message buffer for a given processor */
  CommBuffer & recv_buffer( int p )
  {
    STK_ThrowAssertMsg(p < m_size,"CommSparse::recv_buffer: "<<p<<" out of range [0:"<<m_size<<")");
    stk::util::print_unsupported_version_warning(5, __LINE__, __FILE__);

    if (stk::util::get_common_coupling_version() >= 6) {
      if (m_exchanger) {
        return m_exchanger->get_recv_buf(p);
      } else {
        return m_null_comm_recv_buffer;
      }
    } else {
      return m_recv[p] ;
    }
  }

  /** Obtain the message buffer for a given processor */
  const CommBuffer & recv_buffer( int p ) const
  {
    STK_ThrowAssertMsg(p < m_size,"CommSparse::recv_buffer: "<<p<<" out of range [0:"<<m_size<<")");
    stk::util::print_unsupported_version_warning(5, __LINE__, __FILE__);

    if (stk::util::get_common_coupling_version() >= 6) {
      if (m_exchanger) {
        return m_exchanger->get_recv_buf(p);
      } else {
        return m_null_comm_recv_buffer;
      }
    } else {
      return m_recv[p] ;
    }
  }

  //----------------------------------------
  /** Construct for a to-be-sized communication.
   *  Allocate surrogate send buffers to enable
   *  no-op packing for the purpose of send sizing.
   *  Surrogate send scenario:
   *  1) Surrogate send buffers are "packed" for sizing where
   *     packing sizes are recorded but no data is copied.
   *  2) 'allocate_buffers()' is called to allocate
   *     buffers.
   *  3) Send buffers are identically packed; however, this
   *     packing copies data into the send buffers.
   */
  explicit CommSparse( ParallelMachine comm)
    : m_comm( comm ),
      m_size( parallel_machine_size( comm ) ),
      m_rank( parallel_machine_rank( comm ) ),
#if STK_MIN_COUPLING_VERSION < 6
      m_send(m_size),
      m_recv(m_size),
      m_send_data(),
      m_recv_data(),
      m_send_procs(),
      m_recv_procs(),
#endif
      m_exchanger(nullptr)
  {
    if (comm != MPI_COMM_NULL  && stk::util::get_common_coupling_version() >= 6) {
      m_exchanger = std::make_shared<DataExchangeUnknownPatternNonBlockingCommBuffer>(comm);
    }
  }

  CommSparse(const CommSparse&) = delete;

  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing.
   */
  bool allocate_buffers();

  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing, with user-specified
   *  vectors of procs to send and receive with. Allowing the send/recv
   *  procs to be specified allows the CommSparse class to avoid a
   *  round of communication.
   */
  void allocate_buffers(const std::vector<int>& send_procs, const std::vector<int>& recv_procs);

  /** Communicate send buffers to receive buffers.  */
  bool communicate(bool deallocateSendBuffers = true);

  /** Communicate send buffers to receive buffers, interleave unpacking with
   *    caller-provided functor.  */
  template<typename UNPACK_ALGORITHM>
  void communicate_with_unpack(const UNPACK_ALGORITHM & unpacker)
  {
    auto noExtraWork = [](){};
    communicate_with_extra_work_and_unpacker(noExtraWork, unpacker, true);
  }

  /** Communicate send buffers to receive buffers, interleave extra work.
   */
  template<typename EXTRA_WORK>
  void communicate_with_extra_work(const EXTRA_WORK & extraWork)
  {
    auto noUnpacker = [](int, CommBuffer&){};
    communicate_with_extra_work_and_unpacker(extraWork, noUnpacker, true);
  }

  /** Communicate send buffers to receive buffers, interleave extra work and
   *  also unpacking with caller-provided functor.
   */
  template<typename EXTRA_WORK, typename UNPACK_ALGORITHM>
  void communicate_with_extra_work_and_unpack(const EXTRA_WORK & extraWork, const UNPACK_ALGORITHM & unpacker)
  {
    communicate_with_extra_work_and_unpacker(extraWork, unpacker, true);
  }

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

private:

#if STK_MIN_COUPLING_VERSION < 6
  void allocate_data(std::vector<CommBuffer>& bufs, std::vector<unsigned char>& data);
#endif
  void verify_send_buffers_filled();
  void communicate_with_extra_work_and_unpacker(const std::function<void()>& workFunctor,
                                                const std::function<void(int fromProc, CommBuffer& buf)>& unpackFunctor,
                                                bool deallocateSendBuffers = false);

  ParallelMachine m_comm ;
  int             m_size ;
  int             m_rank ;
#if STK_MIN_COUPLING_VERSION < 6
  std::vector<CommBuffer> m_send;
  std::vector<CommBuffer> m_recv;
  std::vector<unsigned char> m_send_data;
  std::vector<unsigned char> m_recv_data;
  std::vector<int> m_send_procs;
  std::vector<int> m_recv_procs;
#endif

  int             m_num_recvs = DataExchangeUnknownPatternNonBlocking::Unknown;
  std::shared_ptr<DataExchangeUnknownPatternNonBlockingCommBuffer> m_exchanger;

  stk::CommBuffer m_null_comm_send_buffer;
  stk::CommBuffer m_null_comm_recv_buffer;
  std::vector<unsigned char> m_null_comm_storage;
};

template<typename PACK_ALGORITHM>
bool pack_and_communicate(stk::CommSparse & comm, const PACK_ALGORITHM & algorithm, bool deallocateSendBuffers = true)
{
  stk::util::print_unsupported_version_warning(5, __LINE__, __FILE__);

  if (stk::util::get_common_coupling_version() >= 6) {
    algorithm();
    comm.allocate_buffers();
    algorithm();
    return comm.communicate(deallocateSendBuffers);
  } else {
    algorithm();
    const bool actuallySendingOrReceiving = comm.allocate_buffers();
    if (actuallySendingOrReceiving) {
        algorithm();
        comm.communicate(deallocateSendBuffers);
    }
    return actuallySendingOrReceiving; 
  }
}

template<typename UNPACK_ALGORITHM>
void unpack_communications(stk::CommSparse & comm, const UNPACK_ALGORITHM & algorithm)
{
    for(int proc_id=0; proc_id<comm.parallel_size(); ++proc_id)
    {
        if (proc_id != comm.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                algorithm(proc_id);
            }
        }
    }
}

template <typename T>
void pack_vector_to_proc(stk::CommSparse& comm, const std::vector<T>& data, int otherProc)
{
  comm.send_buffer(otherProc).pack(data);
}

template <typename T>
void unpack_vector_from_proc(stk::CommSparse& comm, std::vector<T>& data, int fromProc)
{
  comm.recv_buffer(fromProc).unpack(data);
}

}

//----------------------------------------------------------------------

#endif

