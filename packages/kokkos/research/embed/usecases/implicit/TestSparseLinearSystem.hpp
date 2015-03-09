/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef SPARSELINEARSYSTEM_HPP
#define SPARSELINEARSYSTEM_HPP

#include <cmath>
#include <impl/Kokkos_Timer.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_CrsArray.hpp>

#include <TestBlas1.hpp>
#include <TestCrsMatrix.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_MPI )

template< typename AScalarType ,
          typename VScalarType ,
          class Device >
class Operator {
private:

  typedef CrsMatrix<AScalarType,Device>          matrix_type ;
  typedef View<VScalarType*,LayoutRight,Device>  vector_type ;
  typedef typename vector_type::scalar_type      scalar_type ;
  typedef View<scalar_type*,Device>              buffer_type ;
  typedef typename buffer_type::HostMirror       host_buffer_type ;

  const CrsMatrix<AScalarType,Device> m_A ;
  const ParallelDataMap               m_map ;
  host_buffer_type                    m_host_send_message ;
  host_buffer_type                    m_host_send_buffer ;
  host_buffer_type                    m_host_recv_buffer ;
  std::vector< MPI_Request >          m_recv_request ;
  const int                           m_row_count ;
  const int                           m_col_count ;
  const int                           m_chunk ;

  static const int mpi_tag = 11 ;

  static
  inline
  unsigned send_count_max( const ParallelDataMap & map )
  {
    unsigned max = 0 ;
    for ( unsigned i = 0 ; i < map.host_send.dimension_0() ; ++i ) {
      max = std::max( max , (unsigned) map.host_send(i,1) );
    }
    return max ;
  }

  void post_recv()
  {
    // Post receives:
    const unsigned recv_msg_count = m_map.host_recv.dimension_0();

    scalar_type * ptr = m_host_recv_buffer.ptr_on_device();

    for ( unsigned i = 0 ; i < recv_msg_count ; ++i ) {
      const int proc  = m_map.host_recv(i,0);
      const int count = m_map.host_recv(i,1) * m_chunk ;

      MPI_Irecv( ptr , count * sizeof(scalar_type) , MPI_BYTE ,
                 proc , mpi_tag , m_map.machine.mpi_comm ,
                 & m_recv_request[i] );

      ptr += count ;
    }
  }

  void send( const vector_type & v )
  {
    // Input vector: [ owned-interior , owned-send , receive ]

    const std::pair<unsigned,unsigned> send_range( m_map.count_interior , m_map.count_interior + m_map.count_send );

    vector_type vsend = subview( v , send_range );

    Impl::DeepCopy<HostSpace,typename Device::memory_space>( m_host_send_buffer.ptr_on_device() ,
                                                             vsend.ptr_on_device() ,
                                                             m_map.count_send * m_chunk * sizeof(scalar_type) );

    for ( unsigned i = 0 , j = 0 ; i < m_map.host_send.dimension_0() ; ++i ) {
      const int proc  = m_map.host_send(i,0);
      const int count = m_map.host_send(i,1);

      // Gather send data to contiguous buffer:

      for ( int k = 0 , km = 0 ; k < count ; ++k , ++j ) {
        const int km_end = km + m_chunk ;
        for ( int ki = m_chunk * m_map.host_send_item(j) ; km < km_end ; ++km , ++ki ) {
          m_host_send_message[km] = m_host_send_buffer[ki];
        }
      }

      // MPI_Ssend blocks until
      // (1) a receive is matched for the message and
      // (2) the send buffer can be re-used.
      //
      // It is suggested that MPI_Ssend will have the best performance:
      // http://www.mcs.anl.gov/research/projects/mpi/sendmode.html .

      MPI_Ssend( m_host_send_message.ptr_on_device(),
                 count * m_chunk * sizeof(scalar_type) , MPI_BYTE ,
                 proc , mpi_tag , m_map.machine.mpi_comm );
    }
  }

  void recv( const vector_type & v )
  {
    const size_t recv_msg_count = m_recv_request.size();
    const std::pair<unsigned,unsigned> recv_range( m_map.count_owned , m_map.count_owned + m_map.count_receive );

    const vector_type vrecv = subview( v , recv_range );

    // Wait for receives and verify:

    for ( size_t i = 0 ; i < recv_msg_count ; ++i ) {
      MPI_Status recv_status ;
      int recv_which = 0 ;
      int recv_size  = 0 ;

      MPI_Waitany( recv_msg_count , & m_recv_request[0] , & recv_which , & recv_status );

      const int recv_proc = recv_status.MPI_SOURCE ;

      MPI_Get_count( & recv_status , MPI_BYTE , & recv_size );

      // Verify message properly received:

      const int  expected_proc = m_map.host_recv(recv_which,0);
      const int  expected_size = m_map.host_recv(recv_which,1) *
                                 m_chunk * sizeof(scalar_type);

      if ( ( expected_proc != recv_proc ) ||
           ( expected_size != recv_size ) ) {
        std::ostringstream msg ;
        msg << "MatrixMultiply communication error:"
            << " P" << comm::rank( m_map.machine )
            << " received from P" << recv_proc
            << " size "     << recv_size
            << " expected " << expected_size
            << " from P"    << expected_proc ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Copy received data to device memory.

    Impl::DeepCopy<typename Device::memory_space,HostSpace>( vrecv.ptr_on_device() ,
                                                             m_host_recv_buffer.ptr_on_device() ,
                                                             m_map.count_receive * m_chunk * sizeof(scalar_type) );
  }

public:

  Operator( const ParallelDataMap                & arg_data_map ,
            const CrsMatrix<AScalarType,Device>  & arg_A ,
            const unsigned                         arg_chunk )
    : m_A( arg_A )
    , m_map( arg_data_map )
    , m_host_send_message("MultiplySendMessage", send_count_max( arg_data_map ) * arg_chunk )
    , m_host_send_buffer( "MultiplySendBuffer" , arg_data_map.count_send * arg_chunk )
    , m_host_recv_buffer( "MultiplyRecvBuffer" , arg_data_map.count_receive * arg_chunk )
    , m_recv_request( arg_data_map.host_recv.dimension_0() , MPI_REQUEST_NULL )
    , m_row_count( arg_data_map.count_owned )
    , m_col_count( arg_data_map.count_owned + arg_data_map.count_receive )
    , m_chunk( arg_chunk )
    { }

  void apply( const vector_type  & x ,
              const vector_type  & y )
  {
    // Gather off-processor data for 'x'

    post_recv();

    send( x );

    recv( x );

    Impl::Multiply<matrix_type,vector_type,vector_type>( m_A, m_row_count, m_col_count, x, y);
  }
};

#else

template< typename AScalarType ,
          typename VScalarType ,
          class Device >
class Operator {
private:

  typedef CrsMatrix<AScalarType,Device>          matrix_type ;
  typedef View<VScalarType*,LayoutRight,Device>  vector_type ;
  typedef typename vector_type::scalar_type      scalar_type ;

  const CrsMatrix<AScalarType,Device> m_A ;
  const int                           m_row_count ;
  const int                           m_col_count ;

public:

  Operator( const ParallelDataMap                & arg_data_map ,
            const CrsMatrix<AScalarType,Device>  & arg_A ,
            const unsigned                         /* arg_chunk */ )
    : m_A( arg_A )
    , m_row_count( arg_data_map.count_owned )
    , m_col_count( arg_data_map.count_owned + arg_data_map.count_receive )
    { }

  void apply( const vector_type  & x ,
              const vector_type  & y )
  {
    Impl::Multiply<matrix_type,vector_type,vector_type>( m_A, m_row_count, m_col_count, x, y );
  }
};

#endif

//----------------------------------------------------------------------------

template< typename AScalarType , typename VScalarType , class Device >
void cgsolve(
  const ParallelDataMap                 data_map ,
  const CrsMatrix<AScalarType,Device>   A ,
  const View<VScalarType*,LayoutRight,Device> b ,
  const View<VScalarType*,LayoutRight,Device> x ,
  size_t & iteration ,
  double & normr ,
  double & iter_time ,
  const size_t maximum_iteration = 200 ,
  const double tolerance = 1.0e-12 )
{
  typedef View<VScalarType*,LayoutRight,Device> vector_type ;

  const size_t count_owned = data_map.count_owned ;
  const size_t count_total = data_map.count_owned + data_map.count_receive ;

  Operator<AScalarType,VScalarType,Device> matrix_operator( data_map , A , x.dimension_1() );

  // Need input vector to matvec to be owned + received
  vector_type pAll ( "cg::p" , count_total );

  vector_type p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_owned) );
  vector_type r ( "cg::r" , count_owned );
  vector_type Ap( "cg::Ap", count_owned );

  /* r = b - A * x ; */

  /* p  = x      */ deep_copy( p , x );
  /* Ap = A * p  */ matrix_operator.apply( pAll , Ap );
  /* r  = b - Ap */ waxpby( data_map , 1.0 , b , -1.0 , Ap , r );
  /* p  = r      */ deep_copy( p , r );

  double old_rdot = dot( data_map , r );

  normr     = std::sqrt( old_rdot );
  iteration = 0 ;

  Kokkos::Impl::Timer wall_clock ;

  while ( tolerance < normr && iteration < maximum_iteration ) {

    /* pAp_dot = dot( p , Ap = A * p ) */

    /* Ap = A * p  */ matrix_operator.apply( pAll , Ap );

    const double pAp_dot = dot( data_map , p , Ap );
    const double alpha   = old_rdot / pAp_dot ;

    /* x += alpha * p ;  */ axpy( data_map,  alpha, p , x );
    /* r -= alpha * Ap ; */ axpy( data_map, -alpha, Ap, r );

    const double r_dot = dot( data_map , r );
    const double beta  = r_dot / old_rdot ;

    /* p = r + beta * p ; */ xpby( data_map , r , beta , p );

    normr = std::sqrt( old_rdot = r_dot );
    ++iteration ;
  }

  iter_time = iteration ? wall_clock.seconds() / double(iteration) : 0 ;
}

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef SPARSELINEARSYSTEM_HPP */

