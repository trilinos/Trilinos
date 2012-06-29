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

template< typename T, class > class TestArray ;

template<typename T>
class TestArray< T, KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE device ;
  typedef KokkosArray::Host        host ;

  enum { NP = 1000 ,
         N1 = 3 ,
         N2 = 5 ,
         N3 = 7 };

  typedef KokkosArray::Array< T , device > dView1 ;
  typedef KokkosArray::Array< T[N1] , device > dView2 ;
  typedef KokkosArray::Array< T[N1][N2] , device > dView3 ;
  typedef KokkosArray::Array< T[N1][N2][N3] , device > dView4 ;

  typedef typename dView1::HostMirror  hView1 ;
  typedef typename dView2::HostMirror  hView2 ;
  typedef typename dView3::HostMirror  hView3 ;
  typedef typename dView4::HostMirror  hView4 ;

  TestArray() { run_test(); }

  void run_test()
  {
    dView4 dx , dy , dz ;
    hView4 hx , hy , hz ;
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
    ASSERT_FALSE(hx);
    ASSERT_FALSE(hy);
    ASSERT_FALSE(hz);

    dx = KokkosArray::create_array<dView4>( "dx" , NP );
    dy = KokkosArray::create_array<dView4>( "dy" , NP );

    ASSERT_TRUE(dx);
    ASSERT_TRUE(dy);
    ASSERT_NE( dx , dy );

    ASSERT_EQ( dx.dimension(0) , NP );
    ASSERT_EQ( dx.dimension(1) , N1 );
    ASSERT_EQ( dx.dimension(2) , N2 );
    ASSERT_EQ( dx.dimension(3) , N3 );

    ASSERT_EQ( dy.dimension(0) , NP );
    ASSERT_EQ( dy.dimension(1) , N1 );
    ASSERT_EQ( dy.dimension(2) , N2 );
    ASSERT_EQ( dy.dimension(3) , N3 );

    hx = KokkosArray::create_mirror( dx );
    hy = KokkosArray::create_mirror( dy );

    size_t count = 0 ;
    for ( size_t ip = 0 ; ip < NP ; ++ip ) {
    for ( size_t i1 = 0 ; i1 < dView4::template Dimension<1>::value ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < dView4::template Dimension<2>::value ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < dView4::template Dimension<3>::value ; ++i3 ) {
      hx(ip,i1,i2,i3) = ++count ;
    }}}}

    KokkosArray::deep_copy( dx , hx );
    KokkosArray::deep_copy( dy , dx );
    KokkosArray::deep_copy( hy , dy );

#if 0
    for ( size_t ip = 0 ; ip < NP ; ++ip ) {
    for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
      { ASSERT_EQ( hx(ip,i1,i2,i3) , hy(ip,i1,i2,i3) ); }
    }}}}
#endif

    dz = dx ; // ASSERT_EQ( dx, dz); ASSERT_NE( dy, dz);
    dz = dy ; // ASSERT_EQ( dy, dz); ASSERT_NE( dx, dz);

    dx = dView4();
    ASSERT_FALSE(dx);
    ASSERT_TRUE( dy);
    ASSERT_TRUE( dz);
    dy = dView4();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_TRUE( dz);
    dz = dView4();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
  }
};

}

/*--------------------------------------------------------------------------*/

