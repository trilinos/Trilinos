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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

/*--------------------------------------------------------------------------*/

namespace TestViewSubview {

template< class Space >
void test_left_0()
{
  typedef Kokkos::View< int [2][3][4][5][2][3][4][5] , Kokkos::LayoutLeft , Space >
    view_static_8_type ;

  view_static_8_type  x_static_8("x_static_left_8");

  Kokkos::View<int,Kokkos::LayoutLeft,Space> x0 = Kokkos::subview( x_static_8 , 0, 0, 0, 0, 0, 0, 0, 0 );

  ASSERT_TRUE( & x0() == & x_static_8(0,0,0,0,0,0,0,0) );

  Kokkos::View<int*,Kokkos::LayoutLeft,Space> x1 =
    Kokkos::subview( x_static_8, Kokkos::pair<int,int>(0,2), 1, 2, 3, 0, 1, 2, 3 );

  ASSERT_TRUE( & x1(0) == & x_static_8(0,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x1(1) == & x_static_8(1,1,2,3,0,1,2,3) );

  Kokkos::View<int**,Kokkos::LayoutLeft,Space> x2 =
    Kokkos::subview( x_static_8, Kokkos::pair<int,int>(0,2), 1, 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & x2(0,0) == & x_static_8(0,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x2(0,1) == & x_static_8(0,1,2,3,1,1,2,3) );
  ASSERT_TRUE( & x2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  // Kokkos::View<int**,Kokkos::LayoutLeft,Space> error_2 =
  Kokkos::View<int**,Kokkos::LayoutStride,Space> sx2 =
    Kokkos::subview( x_static_8, 1, Kokkos::pair<int,int>(0,2), 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & sx2(0,0) == & x_static_8(1,0,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(0,1) == & x_static_8(1,0,2,3,1,1,2,3) );
  ASSERT_TRUE( & sx2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  Kokkos::View<int****,Kokkos::LayoutStride,Space> sx4 =
    Kokkos::subview( x_static_8, 0, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 1, Kokkos::pair<int,int>(1,3) /* of [5] */
                               , 2, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 3, Kokkos::pair<int,int>(2,4) /* of [5] */
                   );

  for ( int i0 = 0 ; i0 < (int) sx4.dimension_0() ; ++i0 )
  for ( int i1 = 0 ; i1 < (int) sx4.dimension_1() ; ++i1 )
  for ( int i2 = 0 ; i2 < (int) sx4.dimension_2() ; ++i2 )
  for ( int i3 = 0 ; i3 < (int) sx4.dimension_3() ; ++i3 ) {
    ASSERT_TRUE( & sx4(i0,i1,i2,i3) == & x_static_8(0,0+i0, 1,1+i1, 2,0+i2, 3,2+i3) );
  }
}

template< class Space >
void test_left_1()
{
  typedef Kokkos::View< int ****[2][3][4][5] , Kokkos::LayoutLeft , Space >
    view_static_8_type ;

  view_static_8_type  x_static_8("x_static_left_8",2,3,4,5);

  Kokkos::View<int,Kokkos::LayoutLeft,Space> x0 = Kokkos::subview( x_static_8 , 0, 0, 0, 0, 0, 0, 0, 0 );

  ASSERT_TRUE( & x0() == & x_static_8(0,0,0,0,0,0,0,0) );

  Kokkos::View<int*,Kokkos::LayoutLeft,Space> x1 =
    Kokkos::subview( x_static_8, Kokkos::pair<int,int>(0,2), 1, 2, 3, 0, 1, 2, 3 );

  ASSERT_TRUE( & x1(0) == & x_static_8(0,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x1(1) == & x_static_8(1,1,2,3,0,1,2,3) );

  Kokkos::View<int**,Kokkos::LayoutLeft,Space> x2 =
    Kokkos::subview( x_static_8, Kokkos::pair<int,int>(0,2), 1, 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & x2(0,0) == & x_static_8(0,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & x2(0,1) == & x_static_8(0,1,2,3,1,1,2,3) );
  ASSERT_TRUE( & x2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  // Kokkos::View<int**,Kokkos::LayoutLeft,Space> error_2 =
  Kokkos::View<int**,Kokkos::LayoutStride,Space> sx2 =
    Kokkos::subview( x_static_8, 1, Kokkos::pair<int,int>(0,2), 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & sx2(0,0) == & x_static_8(1,0,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(0,1) == & x_static_8(1,0,2,3,1,1,2,3) );
  ASSERT_TRUE( & sx2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  Kokkos::View<int****,Kokkos::LayoutStride,Space> sx4 =
    Kokkos::subview( x_static_8, 0, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 1, Kokkos::pair<int,int>(1,3) /* of [5] */
                               , 2, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 3, Kokkos::pair<int,int>(2,4) /* of [5] */
                   );

  for ( int i0 = 0 ; i0 < (int) sx4.dimension_0() ; ++i0 )
  for ( int i1 = 0 ; i1 < (int) sx4.dimension_1() ; ++i1 )
  for ( int i2 = 0 ; i2 < (int) sx4.dimension_2() ; ++i2 )
  for ( int i3 = 0 ; i3 < (int) sx4.dimension_3() ; ++i3 ) {
    ASSERT_TRUE( & sx4(i0,i1,i2,i3) == & x_static_8(0,0+i0, 1,1+i1, 2,0+i2, 3,2+i3) );
  }
}

template< class Space >
void test_right_0()
{
  typedef Kokkos::View< int [2][3][4][5][2][3][4][5] , Kokkos::LayoutRight , Space >
    view_static_8_type ;

  view_static_8_type  x_static_8("x_static_right_8");

  Kokkos::View<int,Kokkos::LayoutRight,Space> x0 = Kokkos::subview( x_static_8 , 0, 0, 0, 0, 0, 0, 0, 0 );

  ASSERT_TRUE( & x0() == & x_static_8(0,0,0,0,0,0,0,0) );

  Kokkos::View<int*,Kokkos::LayoutRight,Space> x1 =
    Kokkos::subview( x_static_8, 0, 1, 2, 3, 0, 1, 2, Kokkos::pair<int,int>(1,3) );

  ASSERT_TRUE( & x1(0) == & x_static_8(0,1,2,3,0,1,2,1) );
  ASSERT_TRUE( & x1(1) == & x_static_8(0,1,2,3,0,1,2,2) );

  Kokkos::View<int**,Kokkos::LayoutRight,Space> x2 =
    Kokkos::subview( x_static_8, 0, 1, 2, Kokkos::pair<int,int>(1,3)
                               , 0, 1, 2, Kokkos::pair<int,int>(1,3) );

  ASSERT_TRUE( & x2(0,0) == & x_static_8(0,1,2,1,0,1,2,1) );
  ASSERT_TRUE( & x2(1,0) == & x_static_8(0,1,2,2,0,1,2,1) );
  ASSERT_TRUE( & x2(0,1) == & x_static_8(0,1,2,1,0,1,2,2) );
  ASSERT_TRUE( & x2(1,1) == & x_static_8(0,1,2,2,0,1,2,2) );

  // Kokkos::View<int**,Kokkos::LayoutRight,Space> error_2 =
  Kokkos::View<int**,Kokkos::LayoutStride,Space> sx2 =
    Kokkos::subview( x_static_8, 1, Kokkos::pair<int,int>(0,2), 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & sx2(0,0) == & x_static_8(1,0,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(0,1) == & x_static_8(1,0,2,3,1,1,2,3) );
  ASSERT_TRUE( & sx2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  Kokkos::View<int****,Kokkos::LayoutStride,Space> sx4 =
    Kokkos::subview( x_static_8, 0, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 1, Kokkos::pair<int,int>(1,3) /* of [5] */
                               , 2, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 3, Kokkos::pair<int,int>(2,4) /* of [5] */
                   );

  for ( int i0 = 0 ; i0 < (int) sx4.dimension_0() ; ++i0 )
  for ( int i1 = 0 ; i1 < (int) sx4.dimension_1() ; ++i1 )
  for ( int i2 = 0 ; i2 < (int) sx4.dimension_2() ; ++i2 )
  for ( int i3 = 0 ; i3 < (int) sx4.dimension_3() ; ++i3 ) {
    ASSERT_TRUE( & sx4(i0,i1,i2,i3) == & x_static_8(0,0+i0, 1,1+i1, 2,0+i2, 3,2+i3) );
  }
}

template< class Space >
void test_right_1()
{
  typedef Kokkos::View< int ****[2][3][4][5] , Kokkos::LayoutRight , Space >
    view_static_8_type ;

  view_static_8_type  x_static_8("x_static_right_8",2,3,4,5);

  Kokkos::View<int,Kokkos::LayoutRight,Space> x0 = Kokkos::subview( x_static_8 , 0, 0, 0, 0, 0, 0, 0, 0 );

  ASSERT_TRUE( & x0() == & x_static_8(0,0,0,0,0,0,0,0) );

  Kokkos::View<int*,Kokkos::LayoutRight,Space> x1 =
    Kokkos::subview( x_static_8, 0, 1, 2, 3, 0, 1, 2, Kokkos::pair<int,int>(1,3) );

  ASSERT_TRUE( & x1(0) == & x_static_8(0,1,2,3,0,1,2,1) );
  ASSERT_TRUE( & x1(1) == & x_static_8(0,1,2,3,0,1,2,2) );

  Kokkos::View<int**,Kokkos::LayoutRight,Space> x2 =
    Kokkos::subview( x_static_8, 0, 1, 2, Kokkos::pair<int,int>(1,3)
                               , 0, 1, 2, Kokkos::pair<int,int>(1,3) );

  ASSERT_TRUE( & x2(0,0) == & x_static_8(0,1,2,1,0,1,2,1) );
  ASSERT_TRUE( & x2(1,0) == & x_static_8(0,1,2,2,0,1,2,1) );
  ASSERT_TRUE( & x2(0,1) == & x_static_8(0,1,2,1,0,1,2,2) );
  ASSERT_TRUE( & x2(1,1) == & x_static_8(0,1,2,2,0,1,2,2) );

  // Kokkos::View<int**,Kokkos::LayoutRight,Space> error_2 =
  Kokkos::View<int**,Kokkos::LayoutStride,Space> sx2 =
    Kokkos::subview( x_static_8, 1, Kokkos::pair<int,int>(0,2), 2, 3
                               , Kokkos::pair<int,int>(0,2), 1, 2, 3 );

  ASSERT_TRUE( & sx2(0,0) == & x_static_8(1,0,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(1,0) == & x_static_8(1,1,2,3,0,1,2,3) );
  ASSERT_TRUE( & sx2(0,1) == & x_static_8(1,0,2,3,1,1,2,3) );
  ASSERT_TRUE( & sx2(1,1) == & x_static_8(1,1,2,3,1,1,2,3) );

  Kokkos::View<int****,Kokkos::LayoutStride,Space> sx4 =
    Kokkos::subview( x_static_8, 0, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 1, Kokkos::pair<int,int>(1,3) /* of [5] */
                               , 2, Kokkos::pair<int,int>(0,2) /* of [3] */
                               , 3, Kokkos::pair<int,int>(2,4) /* of [5] */
                   );

  for ( int i0 = 0 ; i0 < (int) sx4.dimension_0() ; ++i0 )
  for ( int i1 = 0 ; i1 < (int) sx4.dimension_1() ; ++i1 )
  for ( int i2 = 0 ; i2 < (int) sx4.dimension_2() ; ++i2 )
  for ( int i3 = 0 ; i3 < (int) sx4.dimension_3() ; ++i3 ) {
    ASSERT_TRUE( & sx4(i0,i1,i2,i3) == & x_static_8(0,0+i0, 1,1+i1, 2,0+i2, 3,2+i3) );
  }
}

}

