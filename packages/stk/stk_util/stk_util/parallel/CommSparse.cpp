// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <vector>

#include <boost/static_assert.hpp> 
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

static const int STK_MPI_TAG_MSG_SIZING      = 01010101;
static const int STK_MPI_TAG_PROC_SIZING = 01110111;
static const int STK_MPI_TAG_DATA        = 11011011;

namespace {

void communicate_any( ParallelMachine p_comm ,
                        const std::vector<CommBuffer > & send ,
                        std::vector<CommBuffer > & recv,
                        const std::vector<int>& send_procs,
                        const std::vector<int>& recv_procs )
{
  static const int mpi_tag = STK_MPI_TAG_DATA ;

  //------------------------------
  // Receive count

  const unsigned num_recv = recv_procs.size();
  const unsigned num_send = send_procs.size();

  //------------------------------
  // Post receives for specific processors with specific sizes

  MPI_Request request_null = MPI_REQUEST_NULL ;
  std::vector<MPI_Request> request(num_recv, request_null);

  for ( unsigned i = 0 ; i < num_recv ; ++i ) {
    int proc = recv_procs[i];
    recv[proc].reset();
    const unsigned recv_size = recv[proc].capacity();
    void * const   recv_buf  = recv[proc].buffer();
    if (MPI_SUCCESS != MPI_Irecv( recv_buf , recv_size , MPI_BYTE , proc , mpi_tag , p_comm , & request[i] )) {
      std::cerr<<"stk::communicate_any ERROR in MPI_Irecv."<<std::endl;
    }
  }

  // This sync is necessary to ensure the IRecvs happen before the Sends.
  MPI_Barrier( p_comm );

  for ( unsigned i = 0 ; i < num_send ; ++i ) {
    int proc = send_procs[i];
    const unsigned send_size = send[proc].capacity();
    void * const   send_buf  = send[proc].buffer();
    MPI_Send( send_buf , send_size , MPI_BYTE , proc , mpi_tag , p_comm );
  }

  std::vector<MPI_Status>  status(  num_recv );
  if (MPI_SUCCESS != MPI_Waitall( num_recv , &request[0] , &status[0] )) {
    std::cerr<<"stk::communicate_any ERROR in MPI_Waitall."<<std::endl;
  }

#ifndef NDEBUG
  const int p_rank = parallel_machine_rank( p_comm );
  for ( unsigned i = 0 ; i < num_recv ; ++i ) {
    MPI_Status * const recv_status = & status[i] ;
    const int recv_proc = recv_status->MPI_SOURCE ;
    const int recv_tag  = recv_status->MPI_TAG ;
    const int recv_plan = recv[recv_proc].capacity();
    int recv_count = 0 ;

    MPI_Get_count( recv_status , MPI_BYTE , & recv_count );

    if ( recv_tag != mpi_tag || recv_count != recv_plan ) {
      std::cerr << "stk::communicate_any LOCAL[" << p_rank << "] ERROR: Recv["
            << recv_proc << "] Size( "
            << recv_count << " != expected " << recv_plan << " ) , " ;
    }
  }
#endif
}

}

#else

// Not parallel


#endif

//----------------------------------------------------------------------

namespace {

inline
size_t align_quad( size_t n )
{
  enum { Size = 4 * sizeof(int) };
  return n + CommBufferAlign<Size>::align(n);
}

}

void CommSparse::rank_error( const char * method , int p ) const
{
  std::ostringstream os ;
  os << "stk::CommSparse::" << method
     << "(" << p << ") ERROR: Not in [0:" << m_size << ")" ;
  throw std::range_error( os.str() );
}

//----------------------------------------------------------------------

CommSparse::~CommSparse()
{
  m_comm = parallel_machine_null();
  m_size = 0 ;
  m_rank = 0 ;
  m_send.clear();
  m_recv.clear();
}

CommSparse::CommSparse()
  : m_comm( parallel_machine_null() ),
    m_size( 0 ), m_rank( 0 ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(),
    m_recv_procs()
{}

CommSparse::CommSparse( ParallelMachine comm)
  : m_comm( comm ),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(),
    m_recv_procs()
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
}

//----------------------------------------------------------------------

void CommSparse::reset_buffers()
{
  for (size_t i=0 ; i<m_send.size(); ++i) {
    m_send[i].reset();
  }
  for (size_t i=0 ; i<m_recv.size(); ++i) {
    m_recv[i].reset();
  }
}

//----------------------------------------------------------------------

void CommSparse::swap_send_recv()
{
  if ( m_recv.empty() ) {
    // ERROR
    std::string
      msg("stk::CommSparse::swap_send_recv(){ NULL recv buffers }" );
    throw std::logic_error( msg );
  }

  m_send.swap(m_recv);
}

//----------------------------------------------------------------------

void CommSparse::allocate_data(std::vector<CommBuffer>& bufs, std::vector<unsigned char>& data)
{
  size_t n_size = 0;
  for ( size_t i = 0 ; i < bufs.size() ; ++i ) {
    n_size += align_quad( bufs[i].size() );
  }

  // Allocate space for buffers

  data.reserve(n_size);
  unsigned char * p_data = &data[0];

  for ( unsigned i = 0 ; i < bufs.size() ; ++i ) {
    CommBuffer & b = bufs[i] ;
    size_t sz = b.size();
    b.m_beg = p_data ;
    b.m_ptr = p_data ;
    b.m_end = p_data + sz ;
    p_data += align_quad( sz );
  }
}

void CommSparse::allocate_buffers()
{
  m_send.resize(m_size);
  m_recv.resize(m_size);

  if (m_size > 1) {
    comm_recv_procs_and_msg_sizes(m_comm, m_send, m_recv, m_send_procs, m_recv_procs);
    allocate_data(m_send, m_send_data);
    allocate_data(m_recv, m_recv_data);
  }
  else {
    allocate_data(m_send, m_send_data);
    m_recv = m_send;
    m_recv_data = m_send_data;
    if (m_send[0].capacity() > 0) {
      m_send_procs.resize(1);
      m_send_procs[0] = 0;
      m_recv_procs = m_send_procs;
    }
  }
}

void CommSparse::allocate_buffers(const std::vector<int>& send_procs, const std::vector<int>& recv_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
  
  m_send_procs = send_procs;
  m_recv_procs = recv_procs;

  if (m_size > 1) {
    comm_recv_msg_sizes(m_comm , send_procs, recv_procs, m_send, m_recv);
    allocate_data(m_send, m_send_data);
    allocate_data(m_recv, m_recv_data);
  }
  else {
    m_recv = m_send;
    m_recv_data = m_send_data;
  }
}

void CommSparse::communicate()
{

#ifndef NDEBUG
  for ( int i = 0 ; i < m_size ; ++i ) {
    // Verify the send buffers have been filled, reset the buffer pointers
    if ( m_send[i].remaining() ) {
      std::ostringstream msg ;
      msg << "stk::CommSparse::communicate LOCAL[" << m_rank << "] ERROR: Send[" << i
          << "] Buffer not filled." ;
      throw std::underflow_error( msg.str() );
    }
  }
#endif

  if ( 1 < m_size ) {
    // Do the communication to exchange the send/recv buffers
    communicate_any( m_comm , m_send , m_recv, m_send_procs, m_recv_procs );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined(STK_HAS_MPI)

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
  int* recvcounts = reinterpret_cast<int*>(&buf[0]);
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

  result = MPI_Reduce_scatter(tmp,&num_recv,recvcounts,uint_type,MPI_SUM,comm);

  if ( result != MPI_SUCCESS ) {
    // PARALLEL ERROR
    std::ostringstream msg ;
    msg << method << " ERROR: " << result << " == MPI_Reduce_scatter" ;
    throw std::runtime_error( msg.str() );
  }

  // do point-to-point send/recvs

  const int mpi_tag = STK_MPI_TAG_PROC_SIZING;

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
    if ( MPI_SUCCESS != result ) {
      // LOCAL ERROR
      std::ostringstream msg ;
      msg << method << " ERROR: " << result << " == MPI_Irecv" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // Send the point-to-point message sizes,

  for ( size_t i = 0 ; i < send_procs.size() ; ++i ) {
    int      dst = send_procs[i];
    unsigned value = send_bufs[dst].size();
    result = MPI_Send( & value , 1 , uint_type, dst , mpi_tag , comm );
    if ( MPI_SUCCESS != result ) {
      // LOCAL ERROR
      std::ostringstream msg ;
      msg << method << " ERROR: " << result << " == MPI_Send" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // Wait for all receives

  {
    MPI_Request * const p_request = (request.empty() ? NULL : & request[0]) ;
    MPI_Status  * const p_status  = (status.empty() ? NULL : & status[0]) ;
    result = MPI_Waitall( num_recv , p_request , p_status );
  }
  if ( MPI_SUCCESS != result ) {
    // LOCAL ERROR ?
    std::ostringstream msg ;
    msg << method << " ERROR: " << result << " == MPI_Waitall" ;
    throw std::runtime_error( msg.str() );
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

  const int mpi_tag = STK_MPI_TAG_MSG_SIZING ;

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
    if ( MPI_SUCCESS != result ) {
      // LOCAL ERROR
      std::ostringstream msg ;
      msg << method << " ERROR: " << result << " == MPI_Irecv" ;
      throw std::runtime_error( msg.str() );
    }
  }

  //barrier to make sure recvs have been posted before sends are launched:
  MPI_Barrier( comm );

  // Send the point-to-point message sizes,

  for ( unsigned i = 0 ; i < num_send ; ++i ) {
    int      dst = send_procs[i];
    unsigned value = send_bufs[dst].size() ;
    result = MPI_Send( & value , 1 , uint_type, dst , mpi_tag , comm );
    if ( MPI_SUCCESS != result ) {
      // LOCAL ERROR
      std::ostringstream msg ;
      msg << method << " ERROR: " << result << " == MPI_Send" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // Wait for all receives

  {
    MPI_Request * const p_request = (request.empty() ? NULL : & request[0]) ;
    MPI_Status  * const p_status  = (status.empty() ? NULL : & status[0]) ;
    result = MPI_Waitall( num_recv , p_request , p_status );
  }
  if ( MPI_SUCCESS != result ) {
    // LOCAL ERROR ?
    std::ostringstream msg ;
    msg << method << " ERROR: " << result << " == MPI_Waitall" ;
    throw std::runtime_error( msg.str() );
  }

  for(unsigned i=0; i<num_recv; ++i) {
    recv_bufs[recv_procs[i]].set_size(recv_sizes[i]);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

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

}

