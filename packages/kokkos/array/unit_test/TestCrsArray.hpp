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
  }

  void run_test_void()
  {
    typedef Kokkos::CrsArray< void , device > dView ;
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

    std::vector<int> x_row_size( LENGTH );
    std::vector<size_t> y_row_size( LENGTH );
    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const int size = ( i + 1 ) % 16 ;
      x_row_size[i] = size ;
      y_row_size[i] = size ;
    }

    dx = Kokkos::create_crsarray<dView>( "dx" , x_row_size );
    dy = Kokkos::create_crsarray<dView>( "dy" , y_row_size );

    ASSERT_TRUE(dx);
    ASSERT_TRUE(dy);
    ASSERT_NE( dx , dy );

    ASSERT_EQ( dx.row_count() , LENGTH );
    ASSERT_EQ( dy.row_count() , LENGTH );

    hx = Kokkos::create_mirror( dx );
    hy = Kokkos::create_mirror( dy );

    ASSERT_EQ( hx.row_count() , LENGTH );
    ASSERT_EQ( hy.row_count() , LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t x_range_begin = hx.row_entry_begin(i);
      const size_t x_range_end   = hx.row_entry_end(i);
      const size_t y_range_begin = hy.row_entry_begin(i);
      const size_t y_range_end   = hy.row_entry_end(i);
      ASSERT_EQ( (int) ( x_range_end - x_range_begin ) , x_row_size[i] );
      ASSERT_EQ( (size_t) ( y_range_end - y_range_begin ) , y_row_size[i] );
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
#if 1
    typedef Kokkos::CrsArray< unsigned , device > dView ;
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

    dx = Kokkos::create_crsarray<dView>( "dx" , graph );
    hx = Kokkos::create_mirror( dx );
   
    ASSERT_EQ( hx.row_count() , LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t begin = hx.row_entry_begin(i);
      const size_t n = hx.row_entry_end(i) - begin ;
      ASSERT_EQ( n , graph[i].size() );
      for ( size_t j = 0 ; j < n ; ++j ) {
        ASSERT_EQ( (int) hx( j + begin ) , graph[i][j] );
      }
    }
#endif
  }

  void run_test_graph2()
  {
#if 0
    typedef Kokkos::CrsArray< unsigned[2] , device > dView ;
    typedef dView::HostMirror hView ;

    enum { LENGTH = 1000 };
    dView dx ;
    hView hx ;

    std::vector< std::vector< int[2] > > graph( LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      graph[i].reserve(8);
      for ( size_t j = 0 ; j < 8 ; ++j ) {
        const int data[2] = { i + j * 3 , 1 + i + j * 3 };
        graph[i].push_back( data );
      }
    }

    dx = Kokkos::create_crsarray<dView>( "dx" , graph );
    hx = Kokkos::create_mirror( dx );
   
    ASSERT_EQ( hx.row_count() , LENGTH );

    for ( size_t i = 0 ; i < LENGTH ; ++i ) {
      const size_t begin = hx.row_entry_begin(i);
      const size_t n = hx.row_entry_end(i) - begin ;
      ASSERT_EQ( n , graph[i].size() );
      for ( size_t j = 0 ; j < n ; ++j ) {
        ASSERT_EQ( (int) hx( j + begin , 0 ) , graph[i][j][0] );
        ASSERT_EQ( (int) hx( j + begin , 1 ) , graph[i][j][1] );
      }
    }
#endif
  }

};

}

/*--------------------------------------------------------------------------*/

