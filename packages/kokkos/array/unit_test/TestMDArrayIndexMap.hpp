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

#include <gtest/gtest.h>

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

template< class > class TestMDArrayIndexMap ;

template<>
class TestMDArrayIndexMap< KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE        device_type ;
  typedef device_type::memory_space  memory_space ;
  typedef Kokkos::Impl::MemoryManager<memory_space>  memory_manager ;

  typedef Kokkos::MDArray< int , device_type > array_type ;
  typedef Kokkos::Impl::MDArrayIndexMapRight< memory_space >  map_right_type ;
  typedef Kokkos::Impl::MDArrayIndexMapLeft<  memory_space >  map_left_type ;

  enum { NP = 1000 , N1 = 10 , N2 = 20 };

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

  TestMDArrayIndexMap()
    : m_left()
    , m_right()
    , m_map_left()
    , m_map_right()
  { run_test(); }

  void run_test()
  {
    typedef Kokkos::Impl::HostMapped< map_left_type >   host_left ;
    typedef Kokkos::Impl::HostMapped< map_right_type >  host_right ;

    const size_t left_alignment_jump =
      memory_manager::preferred_alignment<int>( NP ) - NP ;

    m_left =  Kokkos::create_mdarray< array_type >( NP , N1 , N2 );
    m_right = Kokkos::create_mdarray< array_type >( NP , N1 , N2 );
    m_map_left .assign<int>( NP , N1 , N2 , 0 , 0 , 0 , 0 , 0 );
    m_map_right.assign<int>( NP , N1 , N2 , 0 , 0 , 0 , 0 , 0 );

    int verify = 0 ;
    for ( int j = 0 ; j < N2 ; ++j ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int iwork = 0 ; iwork < NP ; ++iwork , ++verify ) {
          ASSERT_EQ( (unsigned)verify, (unsigned)m_map_left.offset(iwork,i,j));
        }
        verify += left_alignment_jump ;
      }
    }

    verify = 0 ;
    for ( int iwork = 0 ; iwork < NP ; ++iwork ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int j = 0 ; j < N2 ; ++j , ++verify ) {
          ASSERT_EQ( (unsigned)verify, (unsigned)m_map_right.offset(iwork,i,j));
        }
      }
    }

    typedef array_type::HostMirror h_array_type ;

    h_array_type  h_left  = Kokkos::create_mdarray< h_array_type >(  NP , N1 , N2 );
    h_array_type h_right = Kokkos::create_mdarray< h_array_type >( NP , N1 , N2 );

    Kokkos::parallel_for( NP , *this );

    Kokkos::deep_copy( h_left ,  m_left );
    Kokkos::deep_copy( h_right , m_right );

    verify = 0 ;
    for ( int j = 0 ; j < N2 ; ++j ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int iwork = 0 ; iwork < NP ; ++iwork , ++verify ) {
          ASSERT_EQ( (unsigned)verify, (unsigned)h_left(iwork,i,j));
        }
        verify += left_alignment_jump ;
      }
    }

    verify = 0 ;
    for ( int iwork = 0 ; iwork < NP ; ++iwork ) {
      for ( int i = 0 ; i < N1 ; ++i ) {
        for ( int j = 0 ; j < N2 ; ++j , ++verify ) {
          ASSERT_EQ( (unsigned)verify, (unsigned)h_right(iwork,i,j));
        }
      }
    }
  }
};

}

/*--------------------------------------------------------------------------*/

