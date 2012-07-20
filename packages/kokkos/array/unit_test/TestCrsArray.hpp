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

#include <impl/KokkosArray_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< class > class TestCrsArray ;

template<>
class TestCrsArray< KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE device ;


  TestCrsArray()
  {
    run_test_void();
    run_test_graph();
    run_test_graph2();
  }

  void run_test_void()
  {
    typedef KokkosArray::PrefixSum< int , device > dView ;
    typedef dView::HostMirror hView ;

    enum { LENGTH = 1000 };
    dView dx , dy , dz ;
    hView hx , hy , hz ;
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
    ASSERT_FALSE(hx);
    ASSERT_FALSE(hy);
    ASSERT_FALSE(hz);

    std::vector<long> x_size( LENGTH );
    std::vector<size_t> y_size( LENGTH );
    size_t sum = 0 ;
    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const int size = ( i + 1 ) % 16 ;
      x_size[i] = size ;
      y_size[i] = size ;
      sum += size ;
    }

    dx = KokkosArray::create_prefixsum<dView>( "dx" , x_size );
    dy = KokkosArray::create_prefixsum<dView>( "dy" , y_size );

    ASSERT_TRUE(dx);
    ASSERT_TRUE(dy);
    ASSERT_NE( dx , dy );

    ASSERT_EQ( dx.length() , LENGTH );
    ASSERT_EQ( dy.length() , LENGTH );
    ASSERT_EQ( (size_t) dx.sum() , sum );
    ASSERT_EQ( (size_t) dy.sum() , sum );

    hx = KokkosArray::create_mirror( dx );
    hy = KokkosArray::create_mirror( dy );

    ASSERT_EQ( hx.length() , LENGTH );
    ASSERT_EQ( hy.length() , LENGTH );
    ASSERT_EQ( (size_t) hx.sum() , sum );
    ASSERT_EQ( (size_t) hy.sum() , sum );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t x_begin = hx[i];
      const size_t x_end   = hx[i+1];
      const size_t y_begin = hy[i];
      const size_t y_end   = hy[i+1];
      ASSERT_EQ( (long) ( x_end - x_begin ) , x_size[i] );
      ASSERT_EQ( (size_t) ( y_end - y_begin ) , y_size[i] );
    }

    dz = dx ; ASSERT_EQ( dx, dz); ASSERT_NE( dy, dz);
    dz = dy ; ASSERT_EQ( dy, dz); ASSERT_NE( dx, dz);

    dx = dView();
    ASSERT_FALSE(dx);
    ASSERT_TRUE( dy);
    ASSERT_TRUE( dz);
    dy = dView();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_TRUE( dz);
    dz = dView();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
  }

  void run_test_graph()
  {
    typedef KokkosArray::CrsArray< unsigned , device > dView ;
    typedef dView::HostMirror hView ;

    enum { LENGTH = 1000 };
    dView dx ;
    hView hx ;

    std::vector< std::vector< int > > graph( LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      graph[i].reserve(8);
      for ( size_t j = 0 ; j < 8 ; ++j ) {
        graph[i].push_back( i + j * 3 );
      }
    }

    dx = KokkosArray::create_crsarray<dView>( "dx" , graph );
    hx = KokkosArray::create_mirror( dx );
   
    ASSERT_EQ( hx.row_map.length() , LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t begin = hx.row_map[i];
      const size_t n = hx.row_map[i+1] - begin ;
      ASSERT_EQ( n , graph[i].size() );
      for ( size_t j = 0 ; j < n ; ++j ) {
        ASSERT_EQ( (int) hx.entries( j + begin ) , graph[i][j] );
      }
    }
  }

  void run_test_graph2()
  {
    typedef KokkosArray::CrsArray< unsigned[3] , device > dView ;
    typedef dView::HostMirror hView ;

    enum { LENGTH = 10 };

    std::vector< size_t > sizes( LENGTH );

    size_t total_length = 0 ;

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      total_length += ( sizes[i] = 6 + i % 4 );
    }

    dView dx = KokkosArray::create_crsarray<dView>( sizes );
    hView hx = KokkosArray::create_mirror( dx );
    hView mx = KokkosArray::create_mirror( dx );

    ASSERT_EQ( (size_t) dx.row_map.length() , (size_t) LENGTH );
    ASSERT_EQ( (size_t) hx.row_map.length() , (size_t) LENGTH );
    ASSERT_EQ( (size_t) mx.row_map.length() , (size_t) LENGTH );

    ASSERT_EQ( (size_t) dx.entries.dimension(0) , (size_t) total_length );
    ASSERT_EQ( (size_t) hx.entries.dimension(0) , (size_t) total_length );
    ASSERT_EQ( (size_t) mx.entries.dimension(0) , (size_t) total_length );

    ASSERT_EQ( (size_t) dx.entries.dimension(1) , (size_t) 3 );
    ASSERT_EQ( (size_t) hx.entries.dimension(1) , (size_t) 3 );
    ASSERT_EQ( (size_t) mx.entries.dimension(1) , (size_t) 3 );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t entry_begin = hx.row_map[i];
      const size_t entry_end   = hx.row_map[i+1];
      for ( size_t j = entry_begin ; j < entry_end ; ++j ) {
        hx.entries(j,0) = j + 1 ;
        hx.entries(j,1) = j + 2 ;
        hx.entries(j,2) = j + 3 ;
      }
    }

    KokkosArray::deep_copy( dx , hx );

    KokkosArray::deep_copy( mx , dx );
   
    ASSERT_EQ( mx.row_map.length() , LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t entry_begin = mx.row_map[i];
      const size_t entry_end   = mx.row_map[i+1];
      ASSERT_EQ( ( entry_end - entry_begin ) , sizes[i] );
      for ( size_t j = entry_begin ; j < entry_end ; ++j ) {
        ASSERT_EQ( (size_t) mx.entries( j , 0 ) , ( j + 1 ) );
        ASSERT_EQ( (size_t) mx.entries( j , 1 ) , ( j + 2 ) );
        ASSERT_EQ( (size_t) mx.entries( j , 2 ) , ( j + 3 ) );
      }
    }
  }

};

}

/*--------------------------------------------------------------------------*/

