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

void CommSparse::allocate_buffers(const std::vector<int>& send_procs, const std::vector<int>& recv_procs)
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

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Aug 2025
//----------------------------------------------------------------------

#if defined(STK_HAS_MPI)

static const int STK_COMMSPARSE_MPI_TAG_MSG_SIZING  = 10101;
static const int STK_COMMSPARSE_MPI_TAG_PROC_SIZING = 10111;

STK_DEPRECATED
void comm_recv_procs_and_msg_sizes(ParallelMachine comm ,
                                   const std::vector<CommBuffer>& send_bufs ,
                                         std::vector<CommBuffer>& recv_bufs,
                                   std::vector<int>& send_procs,
                                   std::vector<int>& recv_procs)
{
  static const char method[] = "stk::comm_procs_and_msg_recv_sizes" ;

  const int p_size = parallel_machine_size( comm );

  int result = MPI_SUCCESS ;

  MPI_Datatype uint_type = MPI_LONG_LONG;
  if (sizeof(int) == sizeof(unsigned))
    uint_type = MPI_INT;
  else if (sizeof(long) == sizeof(unsigned))
    uint_type = MPI_LONG;
  else if (sizeof(long long) == sizeof(unsigned))
    uint_type = MPI_LONG_LONG;
  else {
    std::ostringstream msg ;
    msg << method << " ERROR: No matching MPI type found for size_t argument";
    throw std::runtime_error(msg.str());
  }

  std::vector<unsigned> buf;
  buf.reserve(p_size*2);
  int* recvcounts = reinterpret_cast<int*>(buf.data());
  unsigned * tmp = &buf[p_size];
  send_procs.clear();
  send_procs.reserve(16);
  for ( int i = 0 ; i < p_size ; ++i ) {
    recvcounts[i] = 1;
    tmp[i] = 0;
    if ( send_bufs[i].size() > 0 ) {
      tmp[i] = 1 ;
      send_procs.push_back(i);
    }
  }

  unsigned num_recv = 0;

  result = MPI_Reduce(tmp,recvcounts,p_size,uint_type,MPI_SUM,0,comm);
  STK_ThrowRequireMsg(result == MPI_SUCCESS, method << " ERROR: " << result << " == MPI_Reduce");

  result = MPI_Scatter(recvcounts,1,uint_type,&num_recv,1,uint_type,0,comm);
  STK_ThrowRequireMsg(result == MPI_SUCCESS, method << " ERROR: " << result << " == MPI_Scatter");

  // do point-to-point send/recvs
  const int mpi_tag = STK_COMMSPARSE_MPI_TAG_PROC_SIZING;

  MPI_Request request_null = MPI_REQUEST_NULL ;
  MPI_Status init_status;
  std::vector<MPI_Request> request( num_recv , request_null );
  std::vector<MPI_Status>  status(  num_recv , init_status );

  // Post receives for point-to-point message sizes

  for ( unsigned i = 0 ; i < num_recv ; ++i ) {
    unsigned    * const p_buf     = & buf[i] ;
    MPI_Request * const p_request = & request[i] ;
    result = MPI_Irecv( p_buf , 1 , uint_type,
                        MPI_ANY_SOURCE , mpi_tag , comm , p_request );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Irecv");
  }

  // Send the point-to-point message sizes,

  for ( size_t i = 0 ; i < send_procs.size() ; ++i ) {
    int      dst = send_procs[i];
    unsigned value = send_bufs[dst].size();
    result = MPI_Send( & value , 1 , uint_type, dst , mpi_tag , comm );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Send");
  }

  // Wait for all receives

  {
    MPI_Request * const p_request = request.data();
    MPI_Status  * const p_status  = status.data();
    result = MPI_Waitall( num_recv , p_request , p_status );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Waitall");
  }

  recv_procs.resize(num_recv);

  // Set the receive message sizes

  for ( unsigned i = 0 ; i < num_recv ; ++i ) {
    MPI_Status * const recv_status = & status[i] ;
    const int recv_proc = recv_status->MPI_SOURCE ;

#ifndef NDEBUG
    //debug-mode-only error check
    const int recv_tag  = recv_status->MPI_TAG ;
    int recv_count  = 0 ;

    MPI_Get_count( recv_status , uint_type , & recv_count );

    if ( recv_tag != mpi_tag || recv_count != 1 ) {
      std::ostringstream msg ;
      const int p_rank = parallel_machine_rank( comm );
      msg << method << " ERROR: Received buffer mismatch " ;
      msg << "P" << p_rank << " <- P" << recv_proc ;
      msg << "  " << 1 << " != " << recv_count ;
      throw std::runtime_error( msg.str() );
    }
#endif

    recv_bufs[ recv_proc ].set_size(buf[i]);
    recv_procs[i] = recv_proc;
  }
}

STK_DEPRECATED
void comm_recv_msg_sizes(ParallelMachine comm ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     const std::vector<CommBuffer>& send_bufs,
                     std::vector<CommBuffer>& recv_bufs)
{
  static const char method[] = "stk::comm_recv_msg_sizes" ;

  int result = MPI_SUCCESS ;

  MPI_Datatype uint_type = MPI_LONG_LONG;
  if (sizeof(int) == sizeof(unsigned))
    uint_type = MPI_INT;
  else if (sizeof(long) == sizeof(unsigned))
    uint_type = MPI_LONG;
  else if (sizeof(long long) == sizeof(unsigned))
    uint_type = MPI_LONG_LONG;
  else {
    std::ostringstream msg ;
    msg << method << " ERROR: No matching MPI type found for unsigned";
    throw std::runtime_error(msg.str());
  }

  // do point-to-point send/recvs

  const int mpi_tag = STK_COMMSPARSE_MPI_TAG_MSG_SIZING ;

  MPI_Request request_null = MPI_REQUEST_NULL ;
  const unsigned num_recv = recv_procs.size();
  const unsigned num_send = send_procs.size();
  std::vector<MPI_Request> request( num_recv , request_null );
  std::vector<MPI_Status>  status(  num_recv );

  std::vector<unsigned> recv_sizes(num_recv);

  // Post receives for point-to-point message sizes

  for ( unsigned i = 0 ; i < num_recv ; ++i ) {
    unsigned    * const p_buf     = & recv_sizes[i] ;
    MPI_Request * const p_request = & request[i] ;
    result = MPI_Irecv( p_buf , 1 , uint_type,
                        recv_procs[i] , mpi_tag , comm , p_request );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Irecv");
  }

  // Send the point-to-point message sizes,

  for ( unsigned i = 0 ; i < num_send ; ++i ) {
    int      dst = send_procs[i];
    unsigned value = send_bufs[dst].size() ;
    result = MPI_Send( & value , 1 , uint_type, dst , mpi_tag , comm );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Send");
  }

  // Wait for all receives

  {
    MPI_Request * const p_request = request.data();
    result = MPI_Waitall( num_recv , p_request , MPI_STATUSES_IGNORE );
    STK_ThrowRequireMsg(MPI_SUCCESS == result, method << " ERROR: " << result << " == MPI_Waitall");
  }

  for(unsigned i=0; i<num_recv; ++i) {
    recv_bufs[recv_procs[i]].set_size(recv_sizes[i]);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

STK_DEPRECATED
void comm_recv_procs_and_msg_sizes(ParallelMachine comm ,
                                   const std::vector<CommBuffer>& send_bufs ,
                                         std::vector<CommBuffer>& recv_bufs,
                                   std::vector<int>& send_procs,
                                   std::vector<int>& recv_procs)
{
  recv_bufs = send_bufs;
  if (send_bufs[0].size() > 0) {
    send_procs.resize(1);
    send_procs[0] = 0;
    recv_procs = send_procs;
  }
}

STK_DEPRECATED
void comm_recv_msg_sizes(ParallelMachine comm ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     const std::vector<CommBuffer>& send_bufs,
                     std::vector<CommBuffer>& recv_bufs)
{
  recv_bufs = send_bufs;
}

//----------------------------------------------------------------------

#endif

#endif // STK_HIDE_DEPRECATED_CODE

}

