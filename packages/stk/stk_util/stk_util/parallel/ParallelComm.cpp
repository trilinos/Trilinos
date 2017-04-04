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

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

enum { STK_MPI_TAG_SIZING = 0 , STK_MPI_TAG_DATA = 1 };

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

//----------------------------------------------------------------------

CommBuffer::CommBuffer()
  : m_beg(nullptr), m_ptr(nullptr), m_end(nullptr)
{ }

CommBuffer::~CommBuffer()
{ }

void CommBuffer::deallocate( const unsigned number , CommBuffer * buffers )
{
  if ( nullptr != buffers ) {
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

  if ( nullptr != size ) {
    for ( unsigned i = 0 ; i < number ; ++i ) {
      n_size += align_quad( size[i] );
    }
  }

  // Allocate space for buffers

  void * const p_malloc = malloc( n_size );

  CommBuffer * const b_base =
    p_malloc != nullptr ? reinterpret_cast<CommBuffer*>(p_malloc)
                        : reinterpret_cast<CommBuffer*>( NULL );

  if ( p_malloc != nullptr ) {

    for ( unsigned i = 0 ; i < number ; ++i ) {
      new( b_base + i ) CommBuffer();
    }

    if ( nullptr != size ) {

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
  m_buffer.m_beg = nullptr ;
  m_buffer.m_ptr = nullptr ;
  m_buffer.m_end = nullptr ;
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

#if defined( STK_HAS_MPI )

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

//
//  Determine the number of items each other process will send to the current processor
//
std::vector<int> ComputeReceiveList(std::vector<int>& sendSizeArray, MPI_Comm &mpi_communicator) {
  const int msg_tag = 10240;
  int num_procs = sendSizeArray.size();
  int my_proc;
  MPI_Comm_rank(mpi_communicator, &my_proc);
  std::vector<int> receiveSizeArray(num_procs, 0);
  //
  //  Determine the total number of messages every processor will receive
  //
  std::vector<int> local_number_to_receive(num_procs, 0);
  std::vector<int> global_number_to_receive(num_procs, 0);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(sendSizeArray[iproc] > 0) local_number_to_receive[iproc] = 1;
  }
  MPI_Allreduce(local_number_to_receive.data(), global_number_to_receive.data(), num_procs, MPI_INT, MPI_SUM, mpi_communicator);
  //
  //  Now each processor knows how many messages it will recive, but does not know the message lengths or where
  //  the messages will be recived from.  Need to extract this information.
  //  Post a recieve for each expected message.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  int num_to_recv = global_number_to_receive[my_proc];
  std::vector<int> recv_size_buffers(num_to_recv);
  for(int imsg = 0; imsg < num_to_recv; ++imsg) {
    int *recv_buffer = &(recv_size_buffers[imsg]);
    MPI_Irecv(recv_buffer, 1, MPI_INT, MPI_ANY_SOURCE,
              msg_tag, mpi_communicator, &recv_handles[imsg]);
  }
  MPI_Barrier(mpi_communicator);
  //
  //  Send message lengths
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(sendSizeArray[iproc] > 0) {
      int send_length = sendSizeArray[iproc];
      MPI_Send(&send_length, 1, MPI_INT, iproc, msg_tag, mpi_communicator);
    }
  }
  //
  //  Get each message and place the length in the proper place in the length array
  //
  for(int imsg = 0; imsg < num_to_recv; ++imsg) {
    MPI_Status status;
    MPI_Wait(&recv_handles[imsg], &status);
    receiveSizeArray[status.MPI_SOURCE] = recv_size_buffers[imsg];
  }

  return receiveSizeArray;
}

#endif

}

