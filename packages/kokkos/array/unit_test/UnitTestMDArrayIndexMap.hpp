/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

template< class > class UnitTestMDArrayIndexMap ;

template<>
class UnitTestMDArrayIndexMap< KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE        device_type ;
  typedef device_type::memory_space  memory_space ;

  typedef Kokkos::MDArrayView< int , device_type > array_type ;
  typedef Kokkos::Impl::MDArrayIndexMap< memory_space , Kokkos::Impl::MDArrayIndexMapRight > map_right_type ;
  typedef Kokkos::Impl::MDArrayIndexMap< memory_space , Kokkos::Impl::MDArrayIndexMapLeft > map_left_type ;

  enum { NP = 1000 , N1 = 10 , N2 = 20 };

  static std::string name()
  {
    std::string tmp ;
    tmp.append( "UnitTestMDArrayIndexMap< Kokkos::" );
    tmp.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
    tmp.append( " >" );
    return tmp ;
  }

  void error( const char * msg ) const
  {
    std::string tmp = name();
    tmp.append( msg );
    throw std::runtime_error( tmp );
  }

  array_type     m_left ;
  array_type     m_right ;
  map_left_type  m_map_left ;
  map_right_type m_map_right ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int iwork ) const
  {
    for ( int i = 0 ; i < N1 ; ++i ) {
      for ( int j = 0 ; j < N2 ; ++j ) {
        m_left( iwork,i,j) = m_map_left.offset( iwork,i,j);
        m_right(iwork,i,j) = m_map_right.offset(iwork,i,j);
      }
    }
  }

  UnitTestMDArrayIndexMap()
    : m_left(  Kokkos::create_mdarray< array_type >( NP , N1 , N2 ) )
    , m_right( Kokkos::create_mdarray< array_type >( NP , N1 , N2 ) )
    , m_map_left(  NP , N1 , N2 , 0 , 0 , 0 , 0 , 0 )
    , m_map_right( NP , N1 , N2 , 0 , 0 , 0 , 0 , 0 )
  {
    typedef Kokkos::Serial< Kokkos::HostMemory ,
                            Kokkos::Impl::MDArrayIndexMapLeft >  host_left ;
    typedef Kokkos::Serial< Kokkos::HostMemory ,
                            Kokkos::Impl::MDArrayIndexMapRight > host_right ;

    typedef Kokkos::MDArrayView< int , host_left >  h_array_left_type ;
    typedef Kokkos::MDArrayView< int , host_right > h_array_right_type ;

    h_array_left_type  h_left  = Kokkos::create_mdarray< h_array_left_type >(  NP , N1 , N2 );
    h_array_right_type h_right = Kokkos::create_mdarray< h_array_right_type >( NP , N1 , N2 );

    Kokkos::parallel_for( NP , *this );

    Kokkos::deep_copy( h_left ,  m_left );
    Kokkos::deep_copy( h_right , m_right );

    int verify = 0 ;
    for ( int j = 0 ; j < N2 ; ++j ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int iwork = 0 ; iwork < NP ; ++iwork , ++verify ) {
          if ( verify != h_left(iwork,i,j) ) {
            std::ostringstream msg ;
            msg << "  left( " << iwork << "," << i << "," << j << ") = "
                << h_left(iwork,i,j) << " != " << verify << std::endl ;
            error( msg.str().c_str() );
          }
        }
      }
    }

    verify = 0 ;
    for ( int iwork = 0 ; iwork < NP ; ++iwork ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int j = 0 ; j < N2 ; ++j , ++verify ) {
          if ( verify != h_right(iwork,i,j) ) {
            std::ostringstream msg ;
            msg << "  right( " << iwork << "," << i << "," << j << ") = "
                << h_right(iwork,i,j) << " != " << verify << std::endl ;
            error( msg.str().c_str() );
          }
        }
      }
    }
  }
};

}

/*--------------------------------------------------------------------------*/

