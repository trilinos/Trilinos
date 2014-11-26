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
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {

  bool CommAll::sm_verbose            = false;
  bool CommAll::sm_verbose_proc0_only = true;

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

enum { STK_MPI_TAG_SIZING = 0 , STK_MPI_TAG_DATA = 1 };

// Communicate in sparse or dense mode, as directed during allocation

namespace {

bool all_to_all_dense( ParallelMachine p_comm ,
                       const CommBuffer * const send ,
                       const CommBuffer * const recv ,
                       std::ostream & msg )
{
  typedef unsigned char * ucharp ;

  static const char method[] = "stk::CommAll::communicate" ;

  int result = MPI_SUCCESS ;

  {
    const int p_size = parallel_machine_size( p_comm );

    std::vector<int> tmp( p_size * 4 );

    int * const send_counts = (tmp.empty() ? NULL : & tmp[0]) ;
    int * const send_displs = send_counts + p_size ;
    int * const recv_counts = send_displs + p_size ;
    int * const recv_displs = recv_counts + p_size ;

    unsigned char * const ps = static_cast<ucharp>(send[0].buffer());
    unsigned char * const pr = static_cast<ucharp>(recv[0].buffer());

    for ( int i = 0 ; i < p_size ; ++i ) {
      const CommBuffer & send_buf = send[i] ;
      const CommBuffer & recv_buf = recv[i] ;

      send_counts[i] = send_buf.capacity();
      recv_counts[i] = recv_buf.capacity();

      send_displs[i] = static_cast<ucharp>(send_buf.buffer()) - ps ;
      recv_displs[i] = static_cast<ucharp>(recv_buf.buffer()) - pr ;
    }

    result = MPI_Alltoallv( ps , send_counts , send_displs , MPI_BYTE ,
                            pr , recv_counts , recv_displs , MPI_BYTE ,
                            p_comm );

    if ( MPI_SUCCESS != result ) {
      msg << method << " GLOBAL ERROR: " << result << " == MPI_Alltoallv" ;
    }
  }

  return MPI_SUCCESS == result ;
}

bool all_to_all_sparse( ParallelMachine p_comm ,
                        const CommBuffer * const send ,
                        const CommBuffer * const recv ,
                        std::ostream & msg )
{
  static const char method[] = "stk::CommAll::communicate" ;
  static const int mpi_tag = STK_MPI_TAG_DATA ;

  int result = MPI_SUCCESS ;

  {
    const int p_size = parallel_machine_size( p_comm );
    const int p_rank = parallel_machine_rank( p_comm );

    //------------------------------
    // Receive count

    unsigned num_recv = 0 ;

    for ( int i = 0 ; i < p_size ; ++i ) {
      if ( recv[i].capacity() ) { ++num_recv ; }
    }

    //------------------------------
    // Post receives for specific processors with specific sizes

    MPI_Request request_null = MPI_REQUEST_NULL ;
    std::vector<MPI_Request> request( num_recv , request_null );
    std::vector<MPI_Status>  status(  num_recv );

    unsigned count = 0 ;

    for ( int i = 0 ; result == MPI_SUCCESS && i < p_size ; ++i ) {
      const unsigned recv_size = recv[i].capacity();
      void * const   recv_buf  = recv[i].buffer();
      if ( recv_size ) {
        result = MPI_Irecv( recv_buf , recv_size , MPI_BYTE ,
                            i , mpi_tag , p_comm , & request[count] );
        ++count ;
      }
    }

    if ( MPI_SUCCESS != result ) {
      msg << method << " LOCAL[" << p_rank << "] ERROR: "
          << result << " == MPI_Irecv , " ;
    }

    //------------------------------
    // Sync to allow ready sends and for a potential error

    int local_error = MPI_SUCCESS == result ? 0 : 1 ;
    int global_error = 0;

    // This sync is necessary to ensure the IRecvs happen before the Rsends. It
    // is an error to call Rsend before the corresponding Irecv.
    result = MPI_Allreduce( & local_error , & global_error ,
                            1 , MPI_INT , MPI_SUM , p_comm );

    if ( MPI_SUCCESS != result ) {
      msg << method << " GLOBAL ERROR: " << result << " == MPI_Allreduce" ;
    }
    else if ( global_error ) {
      result = MPI_ERR_UNKNOWN ;
    }
    else {
      // Everything is local from here on out, no more syncs

      //------------------------------
      // Ready-send the buffers, rotate the send processor
      // in a simple attempt to smooth out the communication traffic.

      for ( int i = 0 ; MPI_SUCCESS == result && i < p_size ; ++i ) {
        const int dst = ( i + p_rank ) % p_size ;
        const unsigned send_size = send[dst].capacity();
        void * const   send_buf  = send[dst].buffer();
        if ( send_size ) {
          result = MPI_Rsend( send_buf , send_size , MPI_BYTE ,
                              dst , mpi_tag , p_comm );
        }
      }

      if ( MPI_SUCCESS != result ) {
        msg << method << " LOCAL ERROR: " << result << " == MPI_Rsend , " ;
      }
      else {
        MPI_Request * const p_request = (request.empty() ? NULL : & request[0]) ;
        MPI_Status  * const p_status  = (status.empty() ? NULL : & status[0]) ;

        result = MPI_Waitall( num_recv , p_request , p_status );
      }

      if ( MPI_SUCCESS != result ) {
        msg << method << " LOCAL[" << p_rank << "] ERROR: "
            << result << " == MPI_Waitall , " ;
      }
      else {

        for ( unsigned i = 0 ; i < num_recv ; ++i ) {
          MPI_Status * const recv_status = & status[i] ;
          const int recv_proc = recv_status->MPI_SOURCE ;
          const int recv_tag  = recv_status->MPI_TAG ;
          const int recv_plan = recv[recv_proc].capacity();
          int recv_count = 0 ;

          MPI_Get_count( recv_status , MPI_BYTE , & recv_count );

          if ( recv_tag != mpi_tag || recv_count != recv_plan ) {
            msg << method << " LOCAL[" << p_rank << "] ERROR: Recv["
                << recv_proc << "] Size( "
                << recv_count << " != " << recv_plan << " ) , " ;
            result = MPI_ERR_UNKNOWN ;
          }
        }
      }
    }
  }

  return MPI_SUCCESS == result ;
}

#ifdef NDEBUG
void communicate_sparse( ParallelMachine p_comm ,
                        const CommBuffer * const send ,
                        const CommBuffer * const recv )
{
  static const int mpi_tag = STK_MPI_TAG_DATA ;

  int result = MPI_SUCCESS ;

  const int p_size = parallel_machine_size( p_comm );

  //------------------------------
  // Receive count

  unsigned num_recv = 0 ;

  for ( int i = 0 ; i < p_size ; ++i ) {
    if ( recv[i].capacity() ) { ++num_recv ; }
  }

  //------------------------------
  // Post receives for specific processors with specific sizes

  MPI_Request request_null = MPI_REQUEST_NULL ;
  std::vector<MPI_Request> request( num_recv , request_null );
  std::vector<MPI_Status>  status(  num_recv );

  unsigned count = 0 ;

  for ( int i = 0 ; result == MPI_SUCCESS && i < p_size ; ++i ) {
    const unsigned recv_size = recv[i].capacity();
    if ( recv_size ) {
      void * const   recv_buf  = recv[i].buffer();
      result = MPI_Irecv( recv_buf , recv_size , MPI_BYTE , i , mpi_tag , p_comm , & request[count] );
      ++count ;
    }
  }

  // This sync is necessary to ensure the IRecvs happen before the Sends.
  MPI_Barrier( p_comm );

  for ( int i = 0 ; i < p_size ; ++i ) {
    const unsigned send_size = send[i].capacity();
    if ( send_size ) {
      void * const   send_buf  = send[i].buffer();
      MPI_Send( send_buf , send_size , MPI_BYTE , i , mpi_tag , p_comm );
    }
  }

  result = MPI_Waitall( num_recv , &request[0] , &status[0] );

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
      std::cerr << "stk::communicate_sparse LOCAL[" << p_rank << "] ERROR: Recv["
            << recv_proc << "] Size( "
            << recv_count << " != expected " << recv_plan << " ) , " ;
    }
  }
#endif
}
#endif

}

#else

// Not parallel

namespace {

bool all_to_all_dense( ParallelMachine ,
                       const CommBuffer * const send ,
                       const CommBuffer * const recv ,
                       std::ostream & )
{ return send == recv ; }

bool all_to_all_sparse( ParallelMachine ,
                        const CommBuffer * const send ,
                        const CommBuffer * const recv ,
                        std::ostream & )
{ return send == recv ; }

}

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

//----------------------------------------------------------------------

void CommBuffer::pack_overflow() const
{
  std::ostringstream os ;
  os << "stk::CommBuffer::pack<T>(...){ overflow by " ;
  os << remaining() ;
  os << " bytes. }" ;
  throw std::overflow_error( os.str() );
}

void CommBuffer::unpack_overflow() const
{
  std::ostringstream os ;
  os << "stk::CommBuffer::unpack<T>(...){ overflow by " ;
  os << remaining();
  os << " bytes. }" ;
  throw std::overflow_error( os.str() );
}

void CommAll::rank_error( const char * method , int p ) const
{
  std::ostringstream os ;
  os << "stk::CommAll::" << method
     << "(" << p << ") ERROR: Not in [0:" << m_size << ")" ;
  throw std::range_error( os.str() );
}

//----------------------------------------------------------------------

CommBuffer::CommBuffer()
  : m_beg(NULL), m_ptr(NULL), m_end(NULL)
{ }

CommBuffer::~CommBuffer()
{ }

void CommBuffer::deallocate( const unsigned number , CommBuffer * buffers )
{
  if ( NULL != buffers ) {
    for ( unsigned i = 0 ; i < number ; ++i ) {
      ( buffers + i )->~CommBuffer();
    }
    free( buffers );
  }
}

CommBuffer * CommBuffer::allocate(
  const unsigned number , const unsigned * const size )
{
  const size_t n_base = align_quad( number * sizeof(CommBuffer) );
  size_t n_size = n_base ;

  if ( NULL != size ) {
    for ( unsigned i = 0 ; i < number ; ++i ) {
      n_size += align_quad( size[i] );
    }
  }

  // Allocate space for buffers

  void * const p_malloc = malloc( n_size );

  CommBuffer * const b_base =
    p_malloc != NULL ? reinterpret_cast<CommBuffer*>(p_malloc)
                     : reinterpret_cast<CommBuffer*>( NULL );

  if ( p_malloc != NULL ) {

    for ( unsigned i = 0 ; i < number ; ++i ) {
      new( b_base + i ) CommBuffer();
    }

    if ( NULL != size ) {

      ucharp ptr = reinterpret_cast<ucharp>( p_malloc );

      ptr += n_base ;

      for ( unsigned i = 0 ; i < number ; ++i ) {
        CommBuffer & b = b_base[i] ;
        b.m_beg = ptr ;
        b.m_ptr = ptr ;
        b.m_end = ptr + size[i] ;
        ptr += align_quad( size[i] );
      }
    }
  }

  return b_base ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

CommAll::~CommAll()
{
  try {
    if ( m_recv == m_send ) {
      CommBuffer::deallocate( m_size , m_send );
    } else {
      CommBuffer::deallocate( m_size , m_send );
      CommBuffer::deallocate( m_size , m_recv );
    }
  } catch(...){}
  m_comm = parallel_machine_null();
  m_size = 0 ;
  m_rank = 0 ;
  m_send = NULL ;
  m_recv = NULL ;
}

CommAll::CommAll(bool propagate_local_error_flags)
  : m_comm( parallel_machine_null() ),
    m_propagate_local_error_flags(propagate_local_error_flags),
    m_size( 0 ), m_rank( 0 ),
    m_bound( 0 ),
    m_max( 0 ),
    m_send(NULL),
    m_recv(NULL)
{}

CommAll::CommAll( ParallelMachine comm, bool propagate_local_error_flags)
  : m_comm( comm ),
    m_propagate_local_error_flags(propagate_local_error_flags),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_bound( 0 ),
    m_max( 0 ),
    m_send(NULL),
    m_recv(NULL)
{
  m_send = CommBuffer::allocate( m_size , NULL );

  if ( NULL == m_send ) {
    std::string msg("stk::CommAll::CommAll FAILED malloc");
    throw std::runtime_error(msg);
  }
}

bool CommAll::allocate_buffers( const unsigned num_msg_bounds ,
                                const bool symmetric ,
                                const bool local_flag )
{
  const unsigned zero = 0 ;
  std::vector<unsigned> tmp( m_size , zero );

  for ( int i = 0 ; i < m_size ; ++i ) {
    tmp[i] = m_send[i].size();
  }

  const unsigned * const send_size = (tmp.empty() ? NULL : & tmp[0]) ;
  const unsigned * const recv_size = symmetric ? (tmp.empty() ? NULL : & tmp[0]) : NULL ;

  return allocate_buffers( m_comm, num_msg_bounds,
                           send_size, recv_size, local_flag );
}

//----------------------------------------------------------------------

void CommAll::reset_buffers()
{
  if ( m_send ) {
    CommBuffer * m = m_send ;
    CommBuffer * const me = m + m_size ;
    for ( ; m != me ; ++m ) { m->reset(); }
  }
  if ( m_recv && 1 < m_size ) {
    CommBuffer * m = m_recv ;
    CommBuffer * const me = m + m_size ;
    for ( ; m != me ; ++m ) { m->reset(); }
  }
}

//----------------------------------------------------------------------

void CommAll::swap_send_recv()
{
  if ( m_recv == NULL ) {
    // ERROR
    std::string
      msg("stk::CommAll::swap_send_recv(){ NULL recv buffers }" );
    throw std::logic_error( msg );
  }

  CommBuffer * tmp_msg = m_send ;
  m_send = m_recv ;
  m_recv = tmp_msg ;
}

//----------------------------------------------------------------------

bool CommAll::allocate_buffers( ParallelMachine comm ,
                                const unsigned num_msg_bounds ,
			        const unsigned * const send_size ,
			        const unsigned * const recv_size ,
			        const bool local_flag )
{

  static const char method[] = "stk::CommAll::allocate_buffers" ;
  const unsigned uzero = 0 ;

  if (send_size == recv_size) {
    return allocate_symmetric_buffers(comm, send_size);
  }

  CommBuffer::deallocate( m_size , m_send );
  CommBuffer::deallocate( m_size , m_recv );

  m_comm = comm ;
  m_size = parallel_machine_size( comm );
  m_rank = parallel_machine_rank( comm );
  m_bound = num_msg_bounds ;

  std::ostringstream msg ;

  //--------------------------------
  // Buffer allocation

  {
    const bool send_none = NULL == send_size ;

    std::vector<unsigned> tmp_send ;

    if ( send_none ) { tmp_send.resize( m_size , uzero ); }

    const unsigned * const send = send_none ? (tmp_send.empty() ? NULL : & tmp_send[0]) : send_size ;

    m_send = CommBuffer::allocate( m_size , send );

    if ( 1 < m_size ) {

      std::vector<unsigned> tmp_recv ;

      const bool recv_tbd = NULL == recv_size ;

      if ( recv_tbd ) { // Had better be globally consistent.

        tmp_recv.resize( m_size , uzero );

        unsigned * const r = (tmp_recv.empty() ? NULL : & tmp_recv[0]) ;

        comm_sizes( m_comm , m_bound , m_max , send , r );
      }

      const unsigned * const recv = recv_tbd  ? (tmp_recv.empty() ? NULL : & tmp_recv[0]) : recv_size ;

      m_recv = CommBuffer::allocate( m_size , recv );
    }
    else {
      m_recv = m_send ;
    }
  }

  if (m_propagate_local_error_flags) {

      bool error_alloc = m_send == NULL || m_recv == NULL ;

      //--------------------------------
      // Propagation of error flag, input flag, and quick/cheap/approximate
      // verification of send and receive messages.
      // Is the number and total size of messages consistent?
      // Sum message counts and sizes for grouped processors.
      // Sent are positive and received are negative.
      // Should finish with all total counts of zero.

      enum { NPSum  = 7 };
      enum { Length = 2 + 2 * NPSum };

      int64_t local_result[ Length ];
      int64_t global_result[ Length ];

      std::fill( local_result , local_result+Length, 0 );

      local_result[ Length - 2 ] = error_alloc ;
      local_result[ Length - 1 ] = local_flag ;

      if ( ! error_alloc ) {

        const unsigned r = 2 * ( m_rank % NPSum );

        for ( int i = 0 ; i < m_size ; ++i ) {
          const unsigned n_send = m_send[i].capacity();
          const unsigned n_recv = m_recv[i].capacity();

          const unsigned s = 2 * ( i % NPSum );

          local_result[s]   += n_send ? 1 : 0 ;
          local_result[s+1] += n_send ;

          local_result[r]   -= n_recv ? 1 : 0 ;
          local_result[r+1] -= n_recv ;
        }
      }

      if (m_size > 1) {
        all_reduce_sum( m_comm , local_result , global_result , Length );
      }
      else {
        std::copy(local_result, local_result+Length, global_result);
      }

      error_alloc   = global_result[ Length - 2 ] ;
      bool global_flag   = global_result[ Length - 1 ] ;

      bool ok = true ;

      for ( unsigned i = 0 ; ok && i < 2 * NPSum ; ++i ) {
        ok = 0 == global_result[i] ;
      }

      if ( error_alloc || ! ok ) {
        msg << method << " ERROR:" ;
        if ( error_alloc   ) { msg << " Failed memory allocation ," ; }
        if ( ! ok          ) { msg << " Parallel inconsistent send/receive ," ; }
        throw std::runtime_error( msg.str() );
      }

      return global_flag ;
  }

  return true;
}

bool CommAll::allocate_buffers( ParallelMachine comm ,
                                unsigned const* const send_sizes,
                                unsigned const* const recv_sizes )
{
  CommBuffer::deallocate( m_size , m_send );
  CommBuffer::deallocate( m_size , m_recv );

  m_comm = comm ;
  m_size = parallel_machine_size( comm );
  m_rank = parallel_machine_rank( comm );
  m_bound = ~0u; // force sparse

  m_send = CommBuffer::allocate( m_size, send_sizes );
  m_recv = CommBuffer::allocate( m_size, recv_sizes );

  return true;
}

bool CommAll::allocate_symmetric_buffers( ParallelMachine comm ,
                                          unsigned const* const buf_sizes )
{
  CommBuffer::deallocate( m_size , m_send );
  CommBuffer::deallocate( m_size , m_recv );

  m_comm = comm ;
  m_size = parallel_machine_size( comm );
  m_rank = parallel_machine_rank( comm );
  m_bound = ~0u; // force sparse

  m_send = CommBuffer::allocate( m_size, buf_sizes );
  m_recv = CommBuffer::allocate( m_size, buf_sizes );

  return true;
}

//----------------------------------------------------------------------

void CommAll::communicate()
{

  for ( int i = 0 ; i < m_size ; ++i ) {
    // Verify the send buffers have been filled, reset the buffer pointers

    if ( m_send[i].remaining() ) {
      std::ostringstream msg ;
      msg << "stk::CommAll::communicate LOCAL[" << m_rank << "] ERROR: Send[" << i
          << "] Buffer not filled." ;
      throw std::underflow_error( msg.str() );
    }
    m_recv[i].reset();
  }

  if ( 1 < m_size ) {
#ifdef NDEBUG
    // Do the communication to exchange the send/recv buffers
    communicate_sparse( m_comm , m_send , m_recv );
#else
    bool ok = false;
    std::ostringstream msg;
    if ( m_bound < m_max ) {
      ok = all_to_all_dense( m_comm , m_send , m_recv , msg );
    }
    else {
      ok = all_to_all_sparse( m_comm , m_send , m_recv , msg );
    }

    if ( ! ok ) { throw std::runtime_error( msg.str() ); }
#endif
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

CommBroadcast::CommBroadcast( ParallelMachine comm , int root_rank )
  : m_comm( comm ),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_root_rank( root_rank ),
    m_buffer()
{}

bool CommBroadcast::allocate_buffer( const bool local_flag )
{
  static const char method[] = "stk::CommBroadcast::allocate_buffer" ;

  int root_rank_min = m_root_rank ;
  int root_rank_max = m_root_rank ;
  unsigned root_send_size = m_root_rank == m_rank ? m_buffer.size() : 0 ;
  unsigned flag = local_flag ;

  all_reduce( m_comm , ReduceMin<1>( & root_rank_min ) &
                       ReduceMax<1>( & root_rank_max ) &
                       ReduceMax<1>( & root_send_size ) &
                       ReduceBitOr<1>( & flag ) );

  if ( root_rank_min != root_rank_max ) {
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED: inconsistent root processor" );
    throw std::runtime_error( msg );
  }

  m_buffer.m_beg = static_cast<CommBuffer::ucharp>( malloc( root_send_size ) );
  m_buffer.m_ptr = m_buffer.m_beg ;
  m_buffer.m_end = m_buffer.m_beg + root_send_size ;

  return flag ;
}

CommBroadcast::~CommBroadcast()
{
  try {
    if ( m_buffer.m_beg ) { free( static_cast<void*>( m_buffer.m_beg ) ); }
  } catch(...) {}
  m_buffer.m_beg = NULL ;
  m_buffer.m_ptr = NULL ;
  m_buffer.m_end = NULL ;
}

CommBuffer & CommBroadcast::recv_buffer()
{
  return m_buffer ;
}

CommBuffer & CommBroadcast::send_buffer()
{
  static const char method[] = "stk::CommBroadcast::send_buffer" ;

  if ( m_root_rank != m_rank ) {
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED: is not root processor" );
    throw std::runtime_error( msg );
  }

  return m_buffer ;
}

void CommBroadcast::communicate()
{
#if defined( STK_HAS_MPI )
  {
    const int count = m_buffer.capacity();
    void * const buf = m_buffer.buffer();

    const int result = MPI_Bcast( buf, count, MPI_BYTE, m_root_rank, m_comm);

    if ( MPI_SUCCESS != result ) {
      std::ostringstream msg ;
      msg << "stk::CommBroadcast::communicate ERROR : "
          << result << " == MPI_Bcast" ;
      throw std::runtime_error( msg.str() );
    }
  }
#endif

  m_buffer.reset();
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

CommGather::~CommGather()
{
  try {
    free( static_cast<void*>( m_send.m_beg ) );

    if ( NULL != m_recv_count ) { free( static_cast<void*>( m_recv_count ) ); }

    if ( NULL != m_recv ) { CommBuffer::deallocate( m_size , m_recv ); }
  } catch(...){}
}

void CommGather::reset()
{
  m_send.reset();

  if ( NULL != m_recv ) {
    for ( int i = 0 ; i < m_size ; ++i ) { m_recv[i].reset(); }
  }
}

CommBuffer & CommGather::recv_buffer( int p )
{
  static CommBuffer empty ;

  return m_size <= p ? empty : (
         m_size <= 1 ? m_send : m_recv[p] );
}

//----------------------------------------------------------------------

CommGather::CommGather( ParallelMachine comm ,
                        int root_rank , unsigned send_size )
  : m_comm( comm ),
    m_size( parallel_machine_size( comm ) ),
    m_rank( parallel_machine_rank( comm ) ),
    m_root_rank( root_rank ),
    m_send(),
    m_recv(NULL),
    m_recv_count(NULL),
    m_recv_displ(NULL)
{
  m_send.m_beg = static_cast<CommBuffer::ucharp>( malloc( send_size ) );
  m_send.m_ptr = m_send.m_beg ;
  m_send.m_end = m_send.m_beg + send_size ;

#if defined( STK_HAS_MPI )

  if ( 1 < m_size ) {

    const bool is_root = m_rank == m_root_rank ;

    if ( is_root ) {
      m_recv_count = static_cast<int*>( malloc(2*m_size*sizeof(int)) );
      m_recv_displ = m_recv_count + m_size ;
    }

    MPI_Gather( & send_size ,    1 , MPI_INT ,
                  m_recv_count , 1 , MPI_INT ,
                  m_root_rank , m_comm );

    if ( is_root ) {
      m_recv = CommBuffer::allocate( m_size ,
                 reinterpret_cast<unsigned*>( m_recv_count ) );

      for ( int i = 0 ; i < m_size ; ++i ) {
        m_recv_displ[i] = m_recv[i].m_beg - m_recv[0].m_beg ;
      }
    }
  }

#endif

}


void CommGather::communicate()
{
#if defined( STK_HAS_MPI )

  if ( 1 < m_size ) {

    const int send_count = m_send.capacity();

    void * const send_buf = m_send.buffer();
    void * const recv_buf = m_rank == m_root_rank ? m_recv->buffer() : NULL ;

    MPI_Gatherv( send_buf , send_count , MPI_BYTE ,
                 recv_buf , m_recv_count , m_recv_displ , MPI_BYTE ,
                 m_root_rank , m_comm );
  }

#endif

  reset();
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

bool comm_dense_sizes( ParallelMachine comm ,
                       const unsigned * const send_size ,
                             unsigned * const recv_size ,
                       bool local_flag )
{
  static const char method[] = "stk::comm_dense_sizes" ;

  const unsigned zero = 0 ;
  const int p_size = parallel_machine_size( comm );

  std::vector<unsigned> send_buf( p_size * 2 , zero );
  std::vector<unsigned> recv_buf( p_size * 2 , zero );

  for ( int i = 0 ; i < p_size ; ++i ) {
    const int i2 = i * 2 ;
    send_buf[i2]   = send_size[i] ;
    send_buf[i2+1] = local_flag ;
  }

  {
    unsigned * const ps = (send_buf.empty() ? NULL : & send_buf[0]) ;
    unsigned * const pr = (recv_buf.empty() ? NULL : & recv_buf[0]) ;

    const int result =
       MPI_Alltoall( ps , 2 , MPI_UNSIGNED , pr , 2 , MPI_UNSIGNED , comm );

    if ( MPI_SUCCESS != result ) {
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED: MPI_SUCCESS != MPI_Alltoall" );
      throw std::runtime_error( msg );
    }
  }

  bool global_flag = false ;

  for ( int i = 0 ; i < p_size ; ++i ) {
    const int i2 = i * 2 ;
    recv_size[i] = recv_buf[i2] ;
    if ( recv_buf[i2+1] ) { global_flag = true ; }
  }

  return global_flag ;
}

//----------------------------------------------------------------------

namespace {

extern "C" {

void sum_np_max_2_op(
  void * inv , void * outv , int * len , ParallelDatatype * )
{
  const int np = *len - 2 ;
  size_t * ind  = static_cast<size_t *>(inv);
  size_t * outd = static_cast<size_t *>(outv);

  // Sum all but the last two
  // the last two are maximum

  for ( int i = 0 ; i < np ; ++i ) {
    *outd += *ind ;
    ++outd ;
    ++ind ;
  }
  if ( outd[0] < ind[0] ) { outd[0] = ind[0] ; }
  if ( outd[1] < ind[1] ) { outd[1] = ind[1] ; }
}

}

}

bool comm_sizes( ParallelMachine comm ,
                 const unsigned   num_msg_bound ,
                       unsigned & num_msg_maximum ,
                 const unsigned * const send_size ,
                       unsigned * const recv_size ,
                 bool local_flag )
{
  static const char method[] = "stk::comm_unknown_sizes" ;
  const size_t uzero = 0 ;

  static MPI_Op mpi_op = MPI_OP_NULL ;

  if ( mpi_op == MPI_OP_NULL ) {
    // Is fully commutative
    MPI_Op_create( sum_np_max_2_op , 1 , & mpi_op );
  }

  const int p_size = parallel_machine_size( comm );
  const int p_rank = parallel_machine_rank( comm );

  int result = MPI_SUCCESS ;

  std::ostringstream msg ;

  num_msg_maximum = 0 ;

  size_t num_recv = 0 ;
  size_t max_msg  = 0 ;
  bool     global_flag = false ;

  {
    std::vector<size_t> send_buf( p_size + 2 , uzero );
    std::vector<size_t> recv_buf( p_size + 2 , uzero );

    size_t * const p_send = (send_buf.empty() ? NULL : & send_buf[0]) ;
    size_t * const p_recv = (recv_buf.empty() ? NULL : & recv_buf[0]) ;

    for ( int i = 0 ; i < p_size ; ++i ) {
      recv_size[i] = 0 ; // Zero output
      if ( send_size[i] ) {
        send_buf[i] = 1 ;
        ++max_msg ;
      }
    }
    send_buf[p_size]   = max_msg ;
    send_buf[p_size+1] = local_flag ;

    BOOST_STATIC_ASSERT(sizeof(long long) == sizeof(size_t));
    result = MPI_Allreduce(p_send,p_recv,p_size,MPI_LONG_LONG,MPI_SUM,comm);

    if ( result != MPI_SUCCESS ) {
      // PARALLEL ERROR
      msg << method << " ERROR: " << result << " == MPI_AllReduce" ;
      throw std::runtime_error( msg.str() );
    }

    result = MPI_Allreduce(p_send+p_size,p_recv+p_size,2,MPI_LONG_LONG,MPI_MAX,comm);

    if ( result != MPI_SUCCESS ) {
      // PARALLEL ERROR
      msg << method << " ERROR: " << result << " == 2nd MPI_AllReduce" ;
      throw std::runtime_error( msg.str() );
    }
    num_recv    = recv_buf[ p_rank ] ;
    max_msg     = recv_buf[ p_size ] ;
    global_flag = recv_buf[ p_size + 1 ] ;

    // max_msg is now the maximum send count,
    // Loop over receive counts to determine
    // if a receive count is larger.

    for ( int i = 0 ; i < p_size ; ++i ) {
      if ( max_msg < recv_buf[i] ) { max_msg = recv_buf[i] ; }
    }
  }

  num_msg_maximum = max_msg ;
  if ( false /*num_msg_bound < max_msg && p_size < 1024*/ ) {
    // Dense, pay for an all-to-all

    result =
      MPI_Alltoall( const_cast<unsigned*>(send_size) , 1 , MPI_LONG_LONG ,
                     recv_size , 1 , MPI_UNSIGNED , comm );

    if ( MPI_SUCCESS != result ) {
      // LOCAL ERROR ?
      msg << method << " ERROR: " << result << " == MPI_Alltoall" ;
      throw std::runtime_error( msg.str() );
    }
  }
  else if ( max_msg ) {
    // Sparse, just do point-to-point

    const int mpi_tag = STK_MPI_TAG_SIZING ;

    MPI_Request request_null = MPI_REQUEST_NULL ;
    std::vector<MPI_Request> request( num_recv , request_null );
    std::vector<MPI_Status>  status(  num_recv );
    std::vector<size_t>    buf( num_recv );

    // Post receives for point-to-point message sizes

    for ( unsigned i = 0 ; i < num_recv ; ++i ) {
      size_t    * const p_buf     = & buf[i] ;
      MPI_Request * const p_request = & request[i] ;
      result = MPI_Irecv( p_buf , 1 , MPI_LONG_LONG ,
                          MPI_ANY_SOURCE , mpi_tag , comm , p_request );
      if ( MPI_SUCCESS != result ) {
        // LOCAL ERROR
        msg << method << " ERROR: " << result << " == MPI_Irecv" ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Send the point-to-point message sizes,
    // rotate the sends in an attempt to balance the message traffic.

    for ( int i = 0 ; i < p_size ; ++i ) {
      int      dst = ( i + p_rank ) % p_size ;
      size_t value = send_size[dst] ;
      if ( value ) {
        result = MPI_Send( & value , 1 , MPI_LONG_LONG , dst , mpi_tag , comm );
        if ( MPI_SUCCESS != result ) {
          // LOCAL ERROR
          msg << method << " ERROR: " << result << " == MPI_Send" ;
          throw std::runtime_error( msg.str() );
        }
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
      msg << method << " ERROR: " << result << " == MPI_Waitall" ;
      throw std::runtime_error( msg.str() );
    }

    // Set the receive message sizes

    for ( unsigned i = 0 ; i < num_recv ; ++i ) {
      MPI_Status * const recv_status = & status[i] ;
      const int recv_proc = recv_status->MPI_SOURCE ;
      const int recv_tag  = recv_status->MPI_TAG ;
      int recv_count  = 0 ;

      MPI_Get_count( recv_status , MPI_LONG_LONG , & recv_count );

      if ( recv_tag != mpi_tag || recv_count != 1 ) {
        msg << method << " ERROR: Received buffer mismatch " ;
        msg << "P" << p_rank << " <- P" << recv_proc ;
        msg << "  " << 1 << " != " << recv_count ;
        throw std::runtime_error( msg.str() );
      }

      const size_t r_size = buf[i] ;
      if (r_size > 2140000000) {
	throw std::overflow_error("r_size > 2.1B");
      }
      recv_size[ recv_proc ] = r_size ;
    }
  }


  return global_flag ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

bool comm_sizes( ParallelMachine ,
                 const unsigned ,
                       unsigned & num_msg_maximum ,
                 const unsigned * const send_size ,
                       unsigned * const recv_size ,
                 bool local_flag )
{
  num_msg_maximum = send_size[0] ? 1 : 0 ;

  recv_size[0] = send_size[0] ;

  return local_flag ;
}

bool comm_dense_sizes( ParallelMachine ,
                       const unsigned * const send_size ,
                             unsigned * const recv_size ,
                       bool local_flag )
{
  recv_size[0] = send_size[0] ;

  return local_flag ;
}

//----------------------------------------------------------------------

#endif

}

