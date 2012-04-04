/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshFixture.hpp>

#include <ParallelDistributedComm.hpp>

#include <Kokkos_Host.hpp>

//----------------------------------------------------------------------------

#ifdef HAVE_MPI

namespace {

template< class ValueType , class Device >
class AsyncExchange {
public:
  typedef Kokkos::CrsArray< void ,        Device >  recv_part_type ;
  typedef Kokkos::CrsArray< unsigned ,    Device >  send_map_type ;
  typedef Kokkos::Impl::MemoryView< ValueType , Device >  buffer_type ;

  static const int mpi_tag = 11 ;


  MPI_Comm                    mpi_comm ;
  const size_t                n_proc ;
  const size_t                my_proc ;
  std::vector< MPI_Request >  recv_request ;
  recv_part_type              recv_part ;
  send_map_type               send_map ;
  buffer_type                 recv_buffer ;
  buffer_type                 send_buffer ;
  size_t                      chunk_size ;

  AsyncExchange( comm::Machine         arg_machine ,
                 const recv_part_type & arg_recv_part ,
                 const send_map_type & arg_send_map ,
                 const size_t          arg_chunk )
  : mpi_comm( arg_machine.mpi_comm )
  , n_proc( comm::size( arg_machine ) )
  , my_proc( comm::rank( arg_machine ) )
  , recv_request()
  , recv_part( arg_recv_part )
  , send_map( arg_send_map )
  , recv_buffer()
  , send_buffer()
  , chunk_size( arg_chunk )
  {
    const size_t recv_msg_base = recv_part.row_entry_end(0);
    const size_t recv_msg_size = recv_part.row_entry_end(n_proc-1) -
                                 recv_msg_base ;
    size_t recv_msg_count = 0 ;
    size_t send_max_size = 0 ;

    for ( size_t i = 1 ; i < n_proc ; ++i ) {

      if ( recv_part.row_entry_end(i) > recv_part.row_entry_begin(i) ) {
        ++recv_msg_count ;
      }

      const size_t send_count = send_map.row_entry_end(i) -
                                send_map.row_entry_begin(i);

      send_max_size = std::max( send_max_size , send_count );
    }

    recv_request.assign( recv_msg_count , MPI_REQUEST_NULL );
    recv_buffer .allocate( recv_msg_size * chunk_size , std::string() );
    send_buffer .allocate( send_max_size * chunk_size , std::string() );

    for ( size_t i = 1 , j = 0 ; i < n_proc ; ++i ) {
      const int proc = ( i + my_proc ) % n_proc ;

      const size_t begin =
        chunk_size * ( recv_part.row_entry_begin(i) - recv_msg_base );

      const size_t count =
        chunk_size * ( recv_part.row_entry_end(i) -
                       recv_part.row_entry_begin(i) );

      if ( count ) {
        MPI_Irecv( recv_buffer.ptr_on_device() + begin ,
                   count * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm , & recv_request[j] );
        ++j ;
      }
    }

    MPI_Barrier( mpi_comm );
  }

  void wait()
  { 
    std::vector< MPI_Status > recv_status( recv_request.size() );

    MPI_Waitall( recv_request.size() , & recv_request[0] , & recv_status[0] );

    for ( size_t i = 1 , j = 0 ; i < n_proc ; ++i ) {
      const int proc = ( i + my_proc ) % n_proc ;
      const size_t recv_count = chunk_size * ( recv_part.row_entry_end(i) -
                                               recv_part.row_entry_begin(i) );

      if ( recv_count ) {
        int recv_size = 0 ;

        MPI_Get_count( & recv_status[j] , MPI_BYTE , & recv_size );

        if ( ( proc != (int) recv_status[j].MPI_SOURCE ) ||
             ( recv_size != (int)( recv_count * sizeof(ValueType) ) ) ) {
          std::ostringstream msg ;
          msg << "AsyncExchange error:"
              << " P" << my_proc << " received from P"
              << recv_status[j].MPI_SOURCE
              << " size " << recv_size
              << " expected " << recv_count * sizeof(ValueType)
              << " from P" << proc ;
          throw std::runtime_error( msg.str() ); 
        }

        ++j ;
      }
    }
  }
};


//----------------------------------------------------------------------------


template< typename coordinate_scalar_type >
void test_box_fixture_verify_parallel(
  comm::Machine machine ,
  const BoxMeshFixture< coordinate_scalar_type , Kokkos::Host > & fixture )
{
  const size_t proc_count = comm::size( machine );
  const size_t proc_local = comm::rank( machine ) ;
  const size_t chunk_size = 3 ;

  // Communicate node coordinates to verify communication and setup.

  AsyncExchange< coordinate_scalar_type , Kokkos::Host >
    exchange( machine , fixture.node_part , fixture.node_send , chunk_size );

  for ( size_t i = 1 ; i < proc_count ; ++i ) {
    const int proc = ( proc_local + i ) % proc_count ;
    const size_t count = chunk_size * ( fixture.node_send.row_entry_end(i) -
                                        fixture.node_send.row_entry_begin(i) );
    if ( count ) {
      // Pack send buffer for one process
      for ( size_t k = 0 ,
                   j = fixture.node_send.row_entry_begin(i) ;
                   j < fixture.node_send.row_entry_end(i) ; ++j ) {
        const size_t node_id = fixture.node_send(j);
        exchange.send_buffer[k++] = fixture.node_coords(node_id,0);
        exchange.send_buffer[k++] = fixture.node_coords(node_id,1);
        exchange.send_buffer[k++] = fixture.node_coords(node_id,2);
      }
      // Ready-send to that process
      MPI_Rsend( exchange.send_buffer.ptr_on_device() ,
                 count * sizeof(coordinate_scalar_type) ,
                 MPI_BYTE , proc , exchange.mpi_tag , machine.mpi_comm );
    }
  }

  // Could do something else here ...

  // Wait for incoming data

  exchange.wait();

  // Unpack recv buffer
   
  unsigned local_error = 0 ;
  unsigned global_error = 0 ;

  for ( size_t k = 0 ,
               j = fixture.node_part.row_entry_begin(1) ;
               j < fixture.node_part.row_entry_end(proc_count) ; ++j ) {

    const coordinate_scalar_type x = exchange.recv_buffer[k++];
    const coordinate_scalar_type y = exchange.recv_buffer[k++];
    const coordinate_scalar_type z = exchange.recv_buffer[k++];

    if ( x != fixture.node_coords(j,0) ||
         y != fixture.node_coords(j,1) ||
         z != fixture.node_coords(j,2) ) {
      ++local_error ;
    }
  }

  MPI_Allreduce( & local_error , & global_error ,
                 1 , MPI_UNSIGNED , MPI_SUM , machine.mpi_comm );

  if ( global_error ) {
    throw std::runtime_error( std::string("coordinate exchange failed") );
  }
}

}

#endif /* #ifdef HAVE_MPI */

//----------------------------------------------------------------------------

void test_box_fixture( comm::Machine machine )
{
  typedef int coordinate_scalar_type ;

  const size_t proc_count = comm::size( machine );
  const size_t proc_local = comm::rank( machine ) ;
  const size_t nodes_nx = 100 ;
  const size_t nodes_ny = 200 ;
  const size_t nodes_nz = 300 ;

  BoxMeshFixture< coordinate_scalar_type , Kokkos::Host >
   fixture( proc_count , proc_local , nodes_nx , nodes_ny , nodes_nz );

  test_box_fixture_verify_parallel( machine , fixture );
}

