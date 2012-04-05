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

#ifndef KOKKOS_FEMESHFIXTURE_HPP
#define KOKKOS_FEMESHFIXTURE_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_CrsArray.hpp>
#include <Kokkos_MDArray.hpp>
#include <Kokkos_MultiVector.hpp>

#include <ParallelDistributedComm.hpp>

//----------------------------------------------------------------------------

template< typename IntType ,
          typename FloatType ,
          unsigned ElemNodeCount ,
          class Device >
struct FEMeshFixture {

  static const IntType element_node_count = ElemNodeCount ;

  // typedef Array< FloatType[3] , Device > node_coords_type ;
  // typedef Array< IntType[ElemNodeCount] , Device > elem_node_ids_type ;

  typedef Kokkos::MDArray<  FloatType ,  Device >  node_coords_type ;
  typedef Kokkos::MDArray<  IntType ,    Device >  elem_node_ids_type ;
  typedef Kokkos::CrsArray< IntType[2] , Device >  node_elem_ids_type ;
  typedef Kokkos::CrsArray< void ,       Device >  node_part_type ;
  typedef Kokkos::CrsArray< unsigned ,   Device >  node_send_type ;

  node_coords_type    node_coords ;
  elem_node_ids_type  elem_node_ids ;
  node_elem_ids_type  node_elem_ids ;
  node_part_type      node_part ;
  node_send_type      node_send ;
};

//----------------------------------------------------------------------------

#ifdef HAVE_MPI

template< class ValueType , class Device >
class AsyncExchange {
public:
  typedef Kokkos::CrsArray< void ,     Device >  recv_part_type ;
  typedef Kokkos::CrsArray< unsigned , Device >  send_map_type ;
  typedef Kokkos::Impl::MemoryView< ValueType , Device >  buffer_type ;

  typedef typename buffer_type::HostMirror    host_buffer_type ;
  typedef typename recv_part_type::HostMirror host_recv_part_type ;
  typedef typename send_map_type ::HostMirror host_send_map_type ;

  static const int mpi_tag = 11 ;

  MPI_Comm                    mpi_comm ;
  host_recv_part_type         host_recv_part ;
  host_send_map_type          host_send_map ;
  const size_t                n_proc ;
  const size_t                my_proc ;
  const size_t                send_count ;
  const size_t                recv_count ;
  const size_t                chunk_size ;
  buffer_type                 recv_buffer ;
  buffer_type                 send_buffer ;
  host_buffer_type            host_recv_buffer ;
  host_buffer_type            host_send_buffer ;
  std::vector< MPI_Request >  recv_request ;

  AsyncExchange( comm::Machine          arg_machine ,
                 const recv_part_type & arg_recv_part ,
                 const send_map_type  & arg_send_map ,
                 const size_t           arg_chunk )
  : mpi_comm( arg_machine.mpi_comm )
  , host_recv_part( Kokkos::create_mirror( arg_recv_part ,
                                           Kokkos::Impl::MirrorUseView() ) )
  , host_send_map(  Kokkos::create_mirror( arg_send_map ,
                                           Kokkos::Impl::MirrorUseView() ) )
  , n_proc(  comm::size( arg_machine ) )
  , my_proc( comm::rank( arg_machine ) )
  , send_count( host_send_map.row_entry_end( n_proc - 1 ) -
                host_send_map.row_entry_end( 0 ) )
  , recv_count( host_recv_part.row_entry_end( n_proc - 1 ) -
                host_recv_part.row_entry_end( 0 ) )
  , chunk_size( arg_chunk )
  , recv_buffer()
  , send_buffer()
  , host_recv_buffer()
  , host_send_buffer()
  , recv_request()
  {
    const size_t recv_msg_base = host_recv_part.row_entry_end(0);

    size_t recv_msg_count = 0 ;

    for ( size_t i = 1 ; i < n_proc ; ++i ) {

      if ( host_recv_part.row_entry_end(i) > host_recv_part.row_entry_begin(i) ) {
        ++recv_msg_count ;
      }
    }

    recv_request.assign( recv_msg_count , MPI_REQUEST_NULL );
    recv_buffer .allocate( recv_count * chunk_size , std::string() );
    send_buffer .allocate( send_count * chunk_size , std::string() );


    typedef Kokkos::Impl::Factory< host_buffer_type ,
                                   Kokkos::Impl::MirrorUseView >
      buffer_mirror_factory ;

    host_recv_buffer = buffer_mirror_factory::create( recv_buffer , recv_count * chunk_size );
    host_send_buffer = buffer_mirror_factory::create( send_buffer , send_count * chunk_size );

    for ( size_t i = 1 , j = 0 ; i < n_proc ; ++i ) {
      const int proc = ( i + my_proc ) % n_proc ;

      const size_t begin =
        chunk_size * ( host_recv_part.row_entry_begin(i) - recv_msg_base );

      const size_t count =
        chunk_size * ( host_recv_part.row_entry_end(i) -
                       host_recv_part.row_entry_begin(i) );

      if ( count ) {
        MPI_Irecv( host_recv_buffer.ptr_on_device() + begin ,
                   count * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm , & recv_request[j] );
        ++j ;
      }
    }
  }

  void send()
  {
    // Copy from the device send buffer to host mirror send buffer
    // and then ready-send the data:

    typedef Kokkos::Impl::Factory< host_buffer_type , buffer_type >
      buffer_copy_factory ;

    // Make sure packing is complete before copying
    Device::fence();

    buffer_copy_factory::deep_copy( host_send_buffer , send_buffer ,
                                    send_count * chunk_size );

    // Wait for all receives to be posted before ready-sending
    MPI_Barrier( mpi_comm );

    for ( size_t i = 1 ; i < n_proc ; ++i ) {
      const int proc = ( i + my_proc) % n_proc ;
      const size_t begin = chunk_size * host_send_map.row_entry_begin(i);
      const size_t count = chunk_size * host_send_map.row_entry_end(i) - begin ;

      if ( count ) { // Ready-send to that process
        MPI_Rsend( host_send_buffer.ptr_on_device() + begin ,
                   count * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm );
      }
    }
  }

  void wait_receive()
  { 
    // Wait for data to be received into the host mirror receive buffer
    // and then deep copy to the device receive buffer.

    std::vector< MPI_Status > recv_status( recv_request.size() );

    MPI_Waitall( recv_request.size() , & recv_request[0] , & recv_status[0] );

    for ( size_t i = 1 , j = 0 ; i < n_proc ; ++i ) {
      const int proc = ( i + my_proc ) % n_proc ;
      const size_t recv_count = chunk_size * ( host_recv_part.row_entry_end(i) -
                                               host_recv_part.row_entry_begin(i) );

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

    typedef Kokkos::Impl::Factory< buffer_type , host_buffer_type >
      buffer_copy_factory ;

    buffer_copy_factory::deep_copy( recv_buffer , host_recv_buffer , 
                                    recv_count * chunk_size );
  }
};

#else /* ! #ifdef HAVE_MPI */


#endif /* ! #ifdef HAVE_MPI */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_FEMESHFIXTURE_HPP */

