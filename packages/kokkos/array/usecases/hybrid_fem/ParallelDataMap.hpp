/*
//@HEADER
// ************************************************************************
// 
//          KokkosArray: Node API and Parallel Node Kernels
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

#ifndef KOKKOS_PARALLELDATAMAP_HPP
#define KOKKOS_PARALLELDATAMAP_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <KokkosArray_Array.hpp>
#include <KokkosArray_Host.hpp>
#include <ParallelComm.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Parallel distributed data mapping
 *
 *  ordering { interior : { owned items not sent elsewhere }
 *             send     : { owned items sent }
 *             receive  : { not-owned items received } }
 *
 *  recv { { N ghosted items from process P : ( P , N ) } }
 *
 *  send { { N send items to process P : ( P , N ) } }
 *
 *  send_item { send item offsets within 'send' range }
 */
struct ParallelDataMap {
  comm::Machine               machine ;
  Array< unsigned[2], Host >  host_recv ;
  Array< unsigned[2], Host >  host_send ;
  Array< unsigned ,   Host >  host_send_item ;
  unsigned                    count_interior ;
  unsigned                    count_send ;
  unsigned                    count_owned ; // = count_interior + count_send
  unsigned                    count_receive ;
};

template< class ArrayType , class Rank = void > struct PackArray ;
template< class ArrayType , class Rank = void > struct UnpackArray ;

template< class ValueType , class Device , class DataMap >
class AsyncExchange ;

} // namespace KokkosArray

//----------------------------------------------------------------------------
// Application call procedure:
//
// construct: AsyncExchange object 
// * pack send buffer on device
// initiate: copy send buffer from device to host
// * dispatch asynchronous local work
// complete: send/receive on host, copy receive buffer to device
// * unpack receive buffer on device
// destroy: AsyncExchange object
//
//----------------------------------------------------------------------------

#ifdef HAVE_MPI

namespace KokkosArray {

template< class ValueType , class Device >
class AsyncExchange< ValueType, Device , KokkosArray::ParallelDataMap > {
public:

  typedef Device                                              device_type ;
  typedef KokkosArray::ParallelDataMap                             data_map_type ;
  typedef KokkosArray::Impl::MemoryView< ValueType , device_type > buffer_dev_type ;
  typedef KokkosArray::Impl::MemoryView< ValueType , Host >        buffer_host_type ;

private:

  static const int mpi_tag = 11 ;

  const data_map_type  data_map ;
  unsigned             chunk_size ;
  unsigned             send_count_max ;
  buffer_host_type     host_send_buffer ;
  buffer_dev_type      dev_buffer ;

public:

  const buffer_dev_type & buffer() const { return dev_buffer ; }

  AsyncExchange( const data_map_type & arg_data_map ,
                 const size_t          arg_chunk )
  : data_map( arg_data_map )
  , chunk_size( arg_chunk )
  , send_count_max( 0 )
  , host_send_buffer()
  , dev_buffer()
  {
    dev_buffer.allocate( arg_chunk * std::max( arg_data_map.count_send ,
                                               arg_data_map.count_receive ),
                         std::string() );

    host_send_buffer.allocate( arg_chunk * data_map.count_send ,
                               std::string() );

    const size_t send_msg_count = arg_data_map.host_send.dimension(0);

    for ( size_t i = 0 ; i < send_msg_count ; ++i ) {
      send_count_max = std::max( send_count_max ,
                                 (unsigned) arg_data_map.host_send(i,1) );
    }
  }

  //------------------------------------------------------------------------

  void send()
  {
    // Copy send buffer from the device to host memory for sending

    KokkosArray::Impl::Factory< buffer_host_type , buffer_dev_type >
      ::deep_copy( host_send_buffer , dev_buffer ,
                   data_map.count_send * chunk_size );

    // Done with the device until communication is complete.
    // Application can dispatch asynchronous work on the device.
  }

  // Application can dispatch local work to device ...
  // No communication progress until main thread calls 'complete'

  void receive()
  {
    const size_t recv_msg_count = data_map.host_recv.dimension(0);
    const size_t send_msg_count = data_map.host_send.dimension(0);

    buffer_host_type host_recv_buffer ;

    host_recv_buffer.allocate( data_map.count_receive * chunk_size ,
                               std::string() );

    std::vector< MPI_Request > recv_request(recv_msg_count,MPI_REQUEST_NULL);

    {
      ValueType * ptr = host_recv_buffer.ptr_on_device();

      for ( size_t i = 0 ; i < recv_msg_count ; ++i ) {
        const int proc  = data_map.host_recv(i,0);
        const int count = data_map.host_recv(i,1) * chunk_size ;

        MPI_Irecv( ptr , count * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , data_map.machine.mpi_comm ,
                   & recv_request[i] );

        ptr += count ;
      }
    }

    // Wait for all receives to be posted before ready-sending

    MPI_Barrier( data_map.machine.mpi_comm );

    {
      buffer_host_type send_msg_buffer ;

      send_msg_buffer.allocate( send_count_max * chunk_size , std::string() );

      for ( size_t i = 0 , j = 0 ; i < send_msg_count ; ++i ) {
        const int proc  = data_map.host_send(i,0);
        const int count = data_map.host_send(i,1);

        for ( int k = 0 , km = 0 ; k < count ; ++k , ++j ) {
          const int km_end = km + chunk_size ;
          int ki = chunk_size * data_map.host_send_item(j);

          for ( ; km < km_end ; ++km , ++ki ) {
            send_msg_buffer[km] = host_send_buffer[ki];
          }
        }

        MPI_Rsend( send_msg_buffer.ptr_on_device(),
                   count * chunk_size * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , data_map.machine.mpi_comm );
      }
    }

    for ( size_t i = 0 ; i < recv_msg_count ; ++i ) {
      MPI_Status recv_status ;
      int recv_which = 0 ;
      int recv_size  = 0 ;

      MPI_Waitany( recv_msg_count , & recv_request[0] ,
                   & recv_which , & recv_status );

      const int recv_proc = recv_status.MPI_SOURCE ;

      MPI_Get_count( & recv_status , MPI_BYTE , & recv_size );

      // Verify message properly received:

      const int  expected_proc = data_map.host_recv(recv_which,0);
      const int  expected_size = data_map.host_recv(recv_which,1) *
                                 chunk_size * sizeof(ValueType);

      if ( ( expected_proc != recv_proc ) ||
           ( expected_size != recv_size ) ) {
        std::ostringstream msg ;
        msg << "AsyncExchange error:"
            << " P" << comm::rank( data_map.machine )
            << " received from P" << recv_proc
            << " size "     << recv_size
            << " expected " << expected_size
            << " from P"    << expected_proc ;
        throw std::runtime_error( msg.str() ); 
      }
    }

    // Touching device:

    KokkosArray::Impl::Factory< buffer_dev_type , buffer_host_type >
      ::deep_copy( dev_buffer , host_recv_buffer , 
                   data_map.count_receive * chunk_size );
  }
};

}

#else /* ! #ifdef HAVE_MPI */


#endif /* ! #ifdef HAVE_MPI */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_PARALLELDATAMAP_HPP */


