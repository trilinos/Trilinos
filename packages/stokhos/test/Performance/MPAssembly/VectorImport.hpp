//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_VECTORIMPORT_HPP
#define KOKKOS_VECTORIMPORT_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include <Teuchos_CommHelpers.hpp>

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
class VectorImport ;

} // namespace Example
} // namespace Kokkos

#if ! defined( KOKKOS_ENABLE_MPI )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
struct VectorImport {

  const Teuchos::RCP<const Teuchos::Comm<int> > comm ;
  const unsigned count_owned ;
  const unsigned count_receive ;

  VectorImport( const Teuchos::RCP<const Teuchos::Comm<int> > arg_comm ,
                const CommMessageType & ,
                const CommMessageType & ,
                const CommIdentType   & ,
                const unsigned arg_count_owned ,
                const unsigned arg_count_receive )
    : comm( arg_comm )
    , count_owned( arg_count_owned )
    , count_receive( arg_count_receive )
    {}

  inline
  void operator()( const VectorType & ) const {}
};


} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else /* defined( KOKKOS_ENABLE_MPI ) */

#include <Teuchos_DefaultMpiComm.hpp>

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
class VectorImport {
private:

  // rank == 1 or array_layout == LayoutRight
  static_assert(
             ( VectorType::rank == 1 ) ||
             std::is_same< typename VectorType::array_layout , Kokkos::LayoutRight >::value,
             "Kokkos::Example::VectorImport Assert Fail: rank != 1 or array_layout != LayoutRight" );

  typedef typename VectorType::HostMirror HostVectorType ;
  typedef typename CommMessageType::HostMirror HostCommMessageType;

  enum { ReceiveInPlace =
    std::is_same< typename VectorType::memory_space ,
                           typename HostVectorType::memory_space >::value };

  const CommMessageType  recv_msg ;
  const CommMessageType  send_msg ;
  const CommIdentType    send_nodeid ;
  HostCommMessageType    host_recv_msg ;
  HostCommMessageType    host_send_msg ;
  VectorType             send_buffer ;
  HostVectorType         host_send_buffer ;
  HostVectorType         host_recv_buffer ;
  unsigned               chunk ;

public:

  const Teuchos::RCP<const Teuchos::Comm<int> >  comm ;
  const unsigned         count_owned ;
  const unsigned         count_receive ;

  struct Pack {
    typedef typename VectorType::execution_space execution_space ;
    const CommIdentType  index ;
    const VectorType     source ;
    const VectorType     buffer ;

    KOKKOS_INLINE_FUNCTION
    void operator()( const unsigned i ) const
      { buffer( i ) = source( index(i) ); }

    Pack( const CommIdentType  & arg_index ,
          const VectorType     & arg_source ,
          const VectorType     & arg_buffer )
      : index( arg_index )
      , source( arg_source )
      , buffer( arg_buffer )
    {
      Kokkos::parallel_for( index.extent(0) , *this );
      execution_space().fence();
    }
  };

  VectorImport( const Teuchos::RCP<const Teuchos::Comm<int> > & arg_comm ,
                const CommMessageType & arg_recv_msg ,
                const CommMessageType & arg_send_msg ,
                const CommIdentType   & arg_send_nodeid ,
                const unsigned arg_count_owned ,
                const unsigned arg_count_receive )
    : recv_msg( arg_recv_msg )
    , send_msg( arg_send_msg )
    , send_nodeid( arg_send_nodeid )
    , host_recv_msg()
    , host_send_msg()
    , send_buffer()
    , host_send_buffer()
    , host_recv_buffer()
    , comm( arg_comm )
    , count_owned( arg_count_owned )
    , count_receive( arg_count_receive )
    {
      host_recv_msg = Kokkos::create_mirror_view( recv_msg );
      host_send_msg = Kokkos::create_mirror_view( send_msg );
      Kokkos::deep_copy( host_recv_msg , recv_msg );
      Kokkos::deep_copy( host_send_msg , send_msg );
      if ( ! ReceiveInPlace ) {
        host_recv_buffer = HostVectorType("recv_buffer",count_receive);
      }

      unsigned send_count = 0 ;
      for ( unsigned i = 0 ; i < send_msg.extent(0) ; ++i ) { send_count += host_send_msg(i,1); }
      send_buffer      = VectorType("send_buffer",send_count);
      host_send_buffer = Kokkos::create_mirror_view( send_buffer );
    }

  inline
  void operator()( const VectorType & v ) const
  {
    typedef typename VectorType::value_type  scalar_type ;
    typedef typename HostVectorType::value_type  host_scalar_type ;

    const Teuchos::MpiComm<int> & teuchos_mpi_comm = dynamic_cast< const Teuchos::MpiComm<int> & >( *comm );

    MPI_Comm mpi_comm = * teuchos_mpi_comm.getRawMpiComm();

    const int mpi_tag = 42 ;
    const unsigned chunk = v.extent(1);

    // Subvector for receives
    const std::pair<unsigned,unsigned> recv_range( count_owned , count_owned + count_receive );
    const VectorType recv_vector = Kokkos::subview( v , recv_range );

    std::vector< MPI_Request > recv_request( recv_msg.extent(0) , MPI_REQUEST_NULL );

    // Post receives
    if (ReceiveInPlace) {
      scalar_type * ptr = recv_vector.data();

      for ( size_t i = 0 ; i < recv_msg.extent(0) ; ++i ) {
        const int proc  = host_recv_msg(i,0);
        const int count = host_recv_msg(i,1) * chunk ;

        MPI_Irecv( ptr , count * sizeof(scalar_type) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm , & recv_request[i] );

        ptr += count ;
      }
    }
    else {
      host_scalar_type * ptr = host_recv_buffer.data();

      for ( size_t i = 0 ; i < recv_msg.extent(0) ; ++i ) {
        const int proc  = host_recv_msg(i,0);
        const int count = host_recv_msg(i,1) * chunk ;

        MPI_Irecv( ptr , count * sizeof(host_scalar_type) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm , & recv_request[i] );

        ptr += count ;
      }

    }

    MPI_Barrier( mpi_comm );

    { // Pack and send
      const Pack pack( send_nodeid , v , send_buffer );

      Kokkos::deep_copy( host_send_buffer , send_buffer );

      host_scalar_type * ptr = host_send_buffer.data();

      for ( size_t i = 0 ; i < send_msg.extent(0) ; ++i ) {
        const int proc  = host_send_msg(i,0);
        const int count = host_send_msg(i,1) * chunk ;

        // MPI_Ssend blocks until
        // (1) a receive is matched for the message and
        // (2) the send buffer can be re-used.
        //
        // It is suggested that MPI_Ssend will have the best performance:
        // http://www.mcs.anl.gov/research/projects/mpi/sendmode.html .

        MPI_Ssend( ptr ,
                   count * sizeof(host_scalar_type) , MPI_BYTE ,
                   proc , mpi_tag , mpi_comm );

        ptr += count ;
      }
    }

    // Wait for receives and verify:

    for ( size_t i = 0 ; i < recv_msg.extent(0) ; ++i ) {
      MPI_Status recv_status ;
      int recv_which = 0 ;
      int recv_size  = 0 ;

      MPI_Waitany( recv_msg.extent(0) , & recv_request[0] , & recv_which , & recv_status );

      const int recv_proc = recv_status.MPI_SOURCE ;

      MPI_Get_count( & recv_status , MPI_BYTE , & recv_size );

      // Verify message properly received:

      const int  expected_proc = host_recv_msg(recv_which,0);
      const int  expected_size = host_recv_msg(recv_which,1) * chunk * sizeof(scalar_type);

      if ( ( expected_proc != recv_proc ) ||
           ( expected_size != recv_size ) ) {

        int local_rank  = 0 ;

        MPI_Comm_rank( mpi_comm , & local_rank );

        std::ostringstream msg ;
        msg << "VectorImport error:"
            << " P" << local_rank
            << " received from P" << recv_proc
            << " size "     << recv_size
            << " expected " << expected_size
            << " from P"    << expected_proc ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Copy received data to device memory.

    if ( ! ReceiveInPlace ) { Kokkos::deep_copy( recv_vector , host_recv_buffer ); }
  }
};

/*
template< class CommMessageType , class CommIdentType ,
          class S, class L, class D, class M >
class VectorImport< CommMessageType, CommIdentType,
                    Kokkos::View<S,L,D,M,Kokkos::Impl::ViewMPVectorContiguous> >
{
public:

  typedef Kokkos::Impl::ViewMPVectorContiguous Specialize;
  typedef Kokkos::View<S,L,D,M,Specialize> VectorType;

private:

  typedef typename VectorType::flat_array_type FlatVectorType;
  typedef VectorImport<CommMessageType, CommIdentType, FlatVectorType> FlatVectorImportType;

  FlatVectorImportType flat_import;

public:

  VectorImport( const Teuchos::RCP<const Teuchos::Comm<int> > & arg_comm ,
                const CommMessageType & arg_recv_msg ,
                const CommMessageType & arg_send_msg ,
                const CommIdentType   & arg_send_nodeid ,
                const unsigned arg_count_owned ,
                const unsigned arg_count_receive ) :
    flat_import( arg_comm,
                 arg_recv_msg,
                 arg_send_msg,
                 arg_send_nodeid,
                 arg_count_owned,
                 arg_count_receive ) {}

  inline void operator()( const VectorType & v ) const
  {
    FlatVectorType flat_v = v;
    flat_import(flat_v);
  }

};
*/

} // namespace Example
} // namespace Kokkos

#endif

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VECTORIMPORT_HPP */
