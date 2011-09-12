/*
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

*/

#ifndef util_ArrayHelper_hpp
#define util_ArrayHelper_hpp

#include <cstddef>
#include <string>
#include <vector>


namespace phdmesh {

//----------------------------------------------------------------------
/** \cond */

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayDimTagList ;

template< class , unsigned > struct ArrayDimTagListAt ;

class ArrayStride {};

template< unsigned R > struct ArrayHelper ;

template< unsigned , unsigned > struct array_bounds_check_rank_is_equal ;
template< unsigned , unsigned > struct array_bounds_check_ordinal_is_less ;

//----------------------------------------------------------------------

void array_stride_natural_tag(
  const int rank , const ArrayDimTag ** const dst ,
                   const ArrayDimTag * const * const src );

void array_stride_fortran_tag(
  const int rank , const ArrayDimTag ** const dst ,
                   const ArrayDimTag * const * const src );

size_t array_stride_size( const int rank , const size_t * const stride );

void array_stride_from_natural_sizes(
  const int rank , size_t * const stride , const size_t * const size );

void array_stride_from_fortran_sizes(
  const int rank , size_t * const stride , const size_t * const size );

void array_stride_to_natural_sizes(
  const int rank , const size_t * const stride , size_t * const size );

void array_stride_to_fortran_sizes(
  const int rank , const size_t * const stride , size_t * const size );

void array_stride_to_natural_index(
  const int rank , const size_t * const stride ,
  const int offset , int * const indices );

void array_stride_to_fortran_index(
  const int rank , const size_t * const stride ,
  const int offset , int * const indices );

//----------------------------------------------------------------------

#define ARRAY_INDEX_NATURAL_8( STRIDE, I1, I2, I3, I4, I5, I6, I7, I8 ) \
  ( I8             + I7 * STRIDE[0] + I6 * STRIDE[1] + I5 * STRIDE[2] + \
    I4 * STRIDE[3] + I3 * STRIDE[4] + I2 * STRIDE[5] + I1 * STRIDE[6] )

#define ARRAY_INDEX_NATURAL_7( STRIDE, I1, I2, I3, I4, I5, I6, I7 ) \
  ( I7             + I6 * STRIDE[0] + I5 * STRIDE[1] + I4 * STRIDE[2] + \
    I3 * STRIDE[3] + I2 * STRIDE[4] + I1 * STRIDE[5] )

#define ARRAY_INDEX_NATURAL_6( STRIDE, I1, I2, I3, I4, I5, I6 ) \
  ( I6             + I5 * STRIDE[0] + I4 * STRIDE[1] + I3 * STRIDE[2] + \
    I2 * STRIDE[3] + I1 * STRIDE[4] )

#define ARRAY_INDEX_NATURAL_5( STRIDE, I1, I2, I3, I4, I5 ) \
  ( I5 + I4 * STRIDE[0] + I3 * STRIDE[1] + I2 * STRIDE[2] + I1 * STRIDE[3] )

#define ARRAY_INDEX_NATURAL_4( STRIDE, I1, I2, I3, I4 ) \
  ( I4 + I3 * STRIDE[0] + I2 * STRIDE[1] + I1 * STRIDE[2] )

#define ARRAY_INDEX_NATURAL_3( STRIDE, I1, I2, I3 ) \
  ( I3 + I2 * STRIDE[0] + I1 * STRIDE[1] )

#define ARRAY_INDEX_NATURAL_2( STRIDE, I1, I2 ) \
  ( I2 + I1 * STRIDE[0] )

#define ARRAY_STRIDE_NATURAL_8( STRIDE, N1, N2, N3, N4, N5, N6, N7, N8 ) \
 STRIDE[7] = N1 * ( STRIDE[6] = N2 * ( STRIDE[5] = N3 * ( STRIDE[4] = N4 * ( \
 STRIDE[3] = N5 * ( STRIDE[2] = N6 * ( STRIDE[1] = N7 * ( STRIDE[0] = N8 )))))))

#define ARRAY_STRIDE_NATURAL_7( STRIDE, N1, N2, N3, N4, N5, N6, N7 ) \
 STRIDE[6] = N1 * ( STRIDE[5] = N2 * ( STRIDE[4] = N3 * ( STRIDE[3] = N4 * ( \
 STRIDE[2] = N5 * ( STRIDE[1] = N6 * ( STRIDE[0] = N7 ))))))

#define ARRAY_STRIDE_NATURAL_6( STRIDE, N1, N2, N3, N4, N5, N6 ) \
 STRIDE[5] = N1 * ( STRIDE[4] = N2 * ( STRIDE[3] = N3 * ( STRIDE[2] = N4 * ( \
 STRIDE[1] = N5 * ( STRIDE[0] = N6 )))))

#define ARRAY_STRIDE_NATURAL_5( STRIDE, N1, N2, N3, N4, N5 ) \
 STRIDE[4] = N1 * ( STRIDE[3] = N2 * ( \
 STRIDE[2] = N3 * ( STRIDE[1] = N4 * ( STRIDE[0] = N5 ))))

#define ARRAY_STRIDE_NATURAL_4( STRIDE, N1, N2, N3, N4 ) \
 STRIDE[3] = N1 * ( STRIDE[2] = N2 * ( STRIDE[1] = N3 * ( STRIDE[0] = N4 )))

#define ARRAY_STRIDE_NATURAL_3( STRIDE, N1, N2, N3 ) \
 STRIDE[2] = N1 * ( STRIDE[1] = N2 * ( STRIDE[0] = N3 ))

#define ARRAY_STRIDE_NATURAL_2( STRIDE, N1, N2 ) \
 STRIDE[1] = N1 * ( STRIDE[0] = N2 )

//----------------------------------------------------------------------

#define ARRAY_INDEX_FORTRAN_8( STRIDE, I1, I2, I3, I4, I5, I6, I7, I8 ) \
  ( I1             + I2 * STRIDE[0] + I3 * STRIDE[1] + I4 * STRIDE[2] + \
    I5 * STRIDE[3] + I6 * STRIDE[4] + I7 * STRIDE[5] + I8 * STRIDE[6] )

#define ARRAY_INDEX_FORTRAN_7( STRIDE, I1, I2, I3, I4, I5, I6, I7 ) \
  ( I1             + I2 * STRIDE[0] + I3 * STRIDE[1] + I4 * STRIDE[2] + \
    I5 * STRIDE[3] + I6 * STRIDE[4] + I7 * STRIDE[5] )

#define ARRAY_INDEX_FORTRAN_6( STRIDE, I1, I2, I3, I4, I5, I6 ) \
  ( I1             + I2 * STRIDE[0] + I3 * STRIDE[1] + I4 * STRIDE[2] + \
    I5 * STRIDE[3] + I6 * STRIDE[4] )

#define ARRAY_INDEX_FORTRAN_5( STRIDE, I1, I2, I3, I4, I5 ) \
  ( I1 + I2 * STRIDE[0] + I3 * STRIDE[1] + I4 * STRIDE[2] + I5 * STRIDE[3] )

#define ARRAY_INDEX_FORTRAN_4( STRIDE, I1, I2, I3, I4 ) \
  ( I1 + I2 * STRIDE[0] + I3 * STRIDE[1] + I4 * STRIDE[2] )

#define ARRAY_INDEX_FORTRAN_3( STRIDE, I1, I2, I3 ) \
  ( I1 + I2 * STRIDE[0] + I3 * STRIDE[1] )

#define ARRAY_INDEX_FORTRAN_2( STRIDE, I1, I2 ) \
  ( I1 + I2 * STRIDE[0] )

#define ARRAY_STRIDE_FORTRAN_8( STRIDE, N1, N2, N3, N4, N5, N6, N7, N8 ) \
 STRIDE[7] = N8 * ( STRIDE[6] = N7 * ( STRIDE[5] = N6 * ( STRIDE[4] = N5 * ( \
 STRIDE[3] = N4 * ( STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 )))))))

#define ARRAY_STRIDE_FORTRAN_7( STRIDE, N1, N2, N3, N4, N5, N6, N7 ) \
 STRIDE[6] = N7 * ( STRIDE[5] = N6 * ( STRIDE[4] = N5 * ( \
 STRIDE[3] = N4 * ( STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 ))))))

#define ARRAY_STRIDE_FORTRAN_6( STRIDE, N1, N2, N3, N4, N5, N6 ) \
 STRIDE[5] = N6 * ( STRIDE[4] = N5 * ( \
 STRIDE[3] = N4 * ( STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 )))))

#define ARRAY_STRIDE_FORTRAN_5( STRIDE, N1, N2, N3, N4, N5 ) \
 STRIDE[4] = N5 * ( \
 STRIDE[3] = N4 * ( STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 ))))

#define ARRAY_STRIDE_FORTRAN_4( STRIDE, N1, N2, N3, N4 ) \
 STRIDE[3] = N4 * ( STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 )))

#define ARRAY_STRIDE_FORTRAN_3( STRIDE, N1, N2, N3 ) \
 STRIDE[2] = N3 * ( STRIDE[1] = N2 * ( STRIDE[0] = N1 ))

#define ARRAY_STRIDE_FORTRAN_2( STRIDE, N1, N2 ) \
 STRIDE[1] = N2 * ( STRIDE[0] = N1 )

//----------------------------------------------------------------------

void array_bounds_checking( int , int );

//----------------------------------------------------------------------

#if defined( ARRAY_BOUNDS_CHECKING )

void array_bounds_checking( const bool ,
                            const int ,
                            const size_t * const ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 ,
                            const int = 0 );

#define ARRAY_BOUNDS_CHECKING_ORDINAL( SIZE , ORD ) \
	array_bounds_checking(SIZE,ORD);

#define ARRAY_BOUNDS_CHECKING_8(I1,I2,I3,I4,I5,I6,I7,I8) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3,I4,I5,I6,I7,I8);

#define ARRAY_BOUNDS_CHECKING_7(I1,I2,I3,I4,I5,I6,I7) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3,I4,I5,I6,I7);

#define ARRAY_BOUNDS_CHECKING_6(I1,I2,I3,I4,I5,I6) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3,I4,I5,I6);

#define ARRAY_BOUNDS_CHECKING_5(I1,I2,I3,I4,I5) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3,I4,I5);

#define ARRAY_BOUNDS_CHECKING_4(I1,I2,I3,I4) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3,I4);

#define ARRAY_BOUNDS_CHECKING_3(I1,I2,I3) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2,I3);

#define ARRAY_BOUNDS_CHECKING_2(I1,I2) \
        array_bounds_checking(Natural,rank(),m_stride,I1,I2);

#define ARRAY_BOUNDS_CHECKING_1(I1) \
        array_bounds_checking(Natural,rank(),m_stride,I1);

#else

#define ARRAY_BOUNDS_CHECKING_ORDINAL( SIZE , ORD ) /**/
#define ARRAY_BOUNDS_CHECKING_8(I1,I2,I3,I4,I5,I6,I7,I8) /**/
#define ARRAY_BOUNDS_CHECKING_7(I1,I2,I3,I4,I5,I6,I7) /**/
#define ARRAY_BOUNDS_CHECKING_6(I1,I2,I3,I4,I5,I6) /**/
#define ARRAY_BOUNDS_CHECKING_5(I1,I2,I3,I4,I5) /**/
#define ARRAY_BOUNDS_CHECKING_4(I1,I2,I3,I4) /**/
#define ARRAY_BOUNDS_CHECKING_3(I1,I2,I3) /**/
#define ARRAY_BOUNDS_CHECKING_2(I1,I2) /**/
#define ARRAY_BOUNDS_CHECKING_1(I1) /**/

#endif

//----------------------------------------------------------------------

template<> struct ArrayHelper<1>
{
  template< typename T >
  static void zero( T * const dst ) { *dst = 0 ; }

  template< typename T >
  static void copy( T * const dst , const T * const src ) { *dst = *src ; }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        *dst = a.dimension(r);
      }
      else {
        *dst = a.dimension(0);
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * , const size_t n )
    { dst[0] = n ; }
};

template<> struct ArrayHelper<2>
{
  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) );
      }
      else {
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) );
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = n * dst[0] ;
    }
};

template<> struct ArrayHelper<3>
{
  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[2] = a.dimension(r-2) * (
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) ));
      }
      else {
        dst[2] = a.dimension(2) * (
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) ));
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = n * dst[1] ;
    }
};

template<> struct ArrayHelper<4>
{

  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
      dst[3] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[3] = a.dimension(r-3) * (
        dst[2] = a.dimension(r-2) * (
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) )));
      }
      else {
        dst[3] = a.dimension(3) * (
        dst[2] = a.dimension(2) * (
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) )));
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = n * dst[2] ;
    }
};

template<> struct ArrayHelper<5> 
{

  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
      dst[3] = 0 ;
      dst[4] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[4] = a.dimension(r-4) * (
        dst[3] = a.dimension(r-3) * (
        dst[2] = a.dimension(r-2) * (
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) ))));
      }
      else {
        dst[4] = a.dimension(4) * (
        dst[3] = a.dimension(3) * (
        dst[2] = a.dimension(2) * (
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) ))));
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = n * dst[3] ;
    }
};

template<> struct ArrayHelper<6>
{
  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
      dst[3] = 0 ;
      dst[4] = 0 ;
      dst[5] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
      dst[5] = src[5] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[5] = a.dimension(r-5) * (
        dst[4] = a.dimension(r-4) * (
        dst[3] = a.dimension(r-3) * (
        dst[2] = a.dimension(r-2) * (
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) )))));
      }
      else {
        dst[5] = a.dimension(5) * (
        dst[4] = a.dimension(4) * (
        dst[3] = a.dimension(3) * (
        dst[2] = a.dimension(2) * (
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) )))));
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
      dst[5] = n * dst[4] ;
    }
};

template<> struct ArrayHelper<7>
{
  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
      dst[3] = 0 ;
      dst[4] = 0 ;
      dst[5] = 0 ;
      dst[6] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
      dst[5] = src[5] ;
      dst[6] = src[6] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      if ( a.natural() ) {
        const int r = a.rank() - 1 ;
        dst[6] = a.dimension(r-6) * (
        dst[5] = a.dimension(r-5) * (
        dst[4] = a.dimension(r-4) * (
        dst[3] = a.dimension(r-3) * (
        dst[2] = a.dimension(r-2) * (
        dst[1] = a.dimension(r-1) * (
        dst[0] = a.dimension(r) ))))));
      }
      else {
        dst[6] = a.dimension(6) * (
        dst[5] = a.dimension(5) * (
        dst[4] = a.dimension(4) * (
        dst[3] = a.dimension(3) * (
        dst[2] = a.dimension(2) * (
        dst[1] = a.dimension(1) * (
        dst[0] = a.dimension(0) ))))));
      }
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
      dst[5] = src[5] ;
      dst[6] = n * dst[5] ;
    }
};

template<> struct ArrayHelper<8>
{
  template< typename T >
  static void zero( T * const dst )
    {
      dst[0] = 0 ;
      dst[1] = 0 ;
      dst[2] = 0 ;
      dst[3] = 0 ;
      dst[4] = 0 ;
      dst[5] = 0 ;
      dst[6] = 0 ;
      dst[7] = 0 ;
    }

  template< typename T >
  static void copy( T * const dst , const T * const src )
    {
      dst[0] = src[0] ;
      dst[1] = src[1] ;
      dst[2] = src[2] ;
      dst[3] = src[3] ;
      dst[4] = src[4] ;
      dst[5] = src[5] ;
      dst[6] = src[6] ;
      dst[7] = src[7] ;
    }

  template< class A >
  static void stride( size_t * const dst , const A & a )
    {
      const int r = a.rank();
      size_t n = 1 ;
      int i ;
      if ( a.natural() ) {
        for ( i = 0 ; i < r ; ++i ) {
          dst[i] = n *= a.dimension( ( r - 1 ) - i );
        }
      }
      else {
        for ( i = 0 ; i < r ; ++i ) {
          dst[i] = n *= a.dimension(i);
        }
      }
      for ( ; i < 8 ; ++i ) { dst[i] = 0 ; } 
    }

  static void stride_append( size_t * const dst ,
                             const size_t * src , const size_t n )
    {
      for ( int i = 0 ; i < 7 ; ++i ) { dst[i] = src[i] ; }
      dst[7] = n * dst[6] ;
    }

  template< class A >
  static void tags( const ArrayDimTag ** const dst , const A & a )
    {
      const int r = a.rank();
      int i ;
      if ( a.natural() ) {
        for ( i = 0 ; i < r ; ++i ) {
          dst[i] = a.tag( ( r - 1 ) - i );
        }
      }
      else {
        for ( i = 0 ; i < r ; ++i ) {
          dst[i] = a.tag(i);
        }
      }
      for ( ; i < 8 ; ++i ) { dst[i] = NULL ; } 
    }
};

//----------------------------------------------------------------------

template< typename Scalar , class T1 , class T2 , class T3 , class T4 ,
                            class T5 , class T6 , class T7 , class T8 >
struct ArrayReverse ;

template< typename Scalar , class T1 >
struct ArrayReverse<Scalar,T1,void,void,void,void,void,void,void> {
  typedef ArrayNatural<Scalar,T1> natural_type ;
  typedef ArrayFortran<Scalar,T1> fortran_type ;
};

template< typename Scalar , class T1 , class T2 >
struct ArrayReverse<Scalar,T1,T2,void,void,void,void,void,void> {
  typedef ArrayNatural<Scalar,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T2,T1> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class T3 >
struct ArrayReverse<Scalar,T1,T2,T3,void,void,void,void,void> {
  typedef ArrayNatural<Scalar,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T3,T2,T1> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class T3 , class T4 >
struct ArrayReverse<Scalar,T1,T2,T3,T4,void,void,void,void> {
  typedef ArrayNatural<Scalar,T4,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T4,T3,T2,T1> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 , class T5 >
struct ArrayReverse<Scalar,T1,T2,T3,T4,T5,void,void,void> {
  typedef ArrayNatural<Scalar,T5,T4,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T5,T4,T3,T2,T1> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 >
struct ArrayReverse<Scalar,T1,T2,T3,T4,T5,T6,void,void> {
  typedef ArrayNatural<Scalar,T6,T5,T4,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T6,T5,T4,T3,T2,T1> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 >
struct ArrayReverse<Scalar,T1,T2,T3,T4,T5,T6,T7,void> {
  typedef ArrayNatural<Scalar,T7,T6,T5,T4,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T7,T6,T5,T4,T3,T2,T1> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayReverse {
  typedef ArrayNatural<Scalar,T8,T7,T6,T5,T4,T3,T2,T1> natural_type ;
  typedef ArrayFortran<Scalar,T8,T7,T6,T5,T4,T3,T2,T1> fortran_type ;
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayTruncate ;

template< typename Scalar >
struct ArrayTruncate<Scalar,void,void,void,void,void,void,void,void> {
  typedef void natural_type ;
  typedef void fortran_type ;
};

template< typename Scalar , class T1 >
struct ArrayTruncate<Scalar,T1,void,void,void,void,void,void,void> {
  typedef void natural_type ;
  typedef void fortran_type ;
};

template< typename Scalar , class T1 , class T2 >
struct ArrayTruncate<Scalar,T1,T2,void,void,void,void,void,void> {
  typedef ArrayNatural<Scalar,T2> natural_type ;
  typedef ArrayFortran<Scalar,T1> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class T3 >
struct ArrayTruncate<Scalar,T1,T2,T3,void,void,void,void,void> {
  typedef ArrayNatural<Scalar,T2,T3> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class T3 , class T4 >
struct ArrayTruncate<Scalar,T1,T2,T3,T4,void,void,void,void> {
  typedef ArrayNatural<Scalar,T2,T3,T4> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2,T3> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 , class T5 >
struct ArrayTruncate<Scalar,T1,T2,T3,T4,T5,void,void,void> {
  typedef ArrayNatural<Scalar,T2,T3,T4,T5> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2,T3,T4> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 >
struct ArrayTruncate<Scalar,T1,T2,T3,T4,T5,T6,void,void> {
  typedef ArrayNatural<Scalar,T2,T3,T4,T5,T6> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2,T3,T4,T5> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 >
struct ArrayTruncate<Scalar,T1,T2,T3,T4,T5,T6,T7,void> {
  typedef ArrayNatural<Scalar,T2,T3,T4,T5,T6,T7> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2,T3,T4,T5,T6> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayTruncate {
  typedef ArrayNatural<Scalar,T2,T3,T4,T5,T6,T7,T8> natural_type ;
  typedef ArrayFortran<Scalar,T1,T2,T3,T4,T5,T6,T7> fortran_type ;
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class TA >
struct ArrayAppend ;

template< typename Scalar , class TA >
struct ArrayAppend<Scalar,void,void,void,void,void,void,void,TA> {
  typedef ArrayNatural<Scalar,TA> natural_type ;
  typedef ArrayFortran<Scalar,TA> fortran_type ;
};

template< typename Scalar , class T1 , class TA >
struct ArrayAppend<Scalar,T1,void,void,void,void,void,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,TA> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class TA >
struct ArrayAppend<Scalar,T1,T2,void,void,void,void,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1,T2>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,TA> fortran_type ;
};

template< typename Scalar , class T1 , class T2 , class T3 , class TA >
struct ArrayAppend<Scalar,T1,T2,T3,void,void,void,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1,T2,T3>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,T3,TA> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 , class TA >
struct ArrayAppend<Scalar,T1,T2,T3,T4,void,void,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1,T2,T3,T4>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,T3,T4,TA> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 , class T5 , class TA >
struct ArrayAppend<Scalar,T1,T2,T3,T4,T5,void,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1,T2,T3,T4,T5>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,T3,T4,T5,TA> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class TA >
struct ArrayAppend<Scalar,T1,T2,T3,T4,T5,T6,void,TA> {
  typedef ArrayNatural<Scalar,TA,T1,T2,T3,T4,T5,T6>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,T3,T4,T5,T6,TA> fortran_type ;
};

template< typename Scalar ,
          class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class TA >
struct ArrayAppend {
  typedef ArrayNatural<Scalar,TA,T1,T2,T3,T4,T5,T6,T7>    natural_type ;
  typedef ArrayFortran<Scalar,   T1,T2,T3,T4,T5,T6,T7,TA> fortran_type ;
};

//----------------------------------------------------------------------

template<> struct array_bounds_check_rank_is_equal<1,1> {};
template<> struct array_bounds_check_rank_is_equal<2,2> {};
template<> struct array_bounds_check_rank_is_equal<3,3> {};
template<> struct array_bounds_check_rank_is_equal<4,4> {};
template<> struct array_bounds_check_rank_is_equal<5,5> {};
template<> struct array_bounds_check_rank_is_equal<6,6> {};
template<> struct array_bounds_check_rank_is_equal<7,7> {};
template<> struct array_bounds_check_rank_is_equal<8,8> {};

//----------------------------------------------------------------------

template<> struct array_bounds_check_ordinal_is_less<8,0> {};
template<> struct array_bounds_check_ordinal_is_less<8,1> {};
template<> struct array_bounds_check_ordinal_is_less<8,2> {};
template<> struct array_bounds_check_ordinal_is_less<8,3> {};
template<> struct array_bounds_check_ordinal_is_less<8,4> {};
template<> struct array_bounds_check_ordinal_is_less<8,5> {};
template<> struct array_bounds_check_ordinal_is_less<8,6> {};
template<> struct array_bounds_check_ordinal_is_less<8,7> {};

template<> struct array_bounds_check_ordinal_is_less<7,0> {};
template<> struct array_bounds_check_ordinal_is_less<7,1> {};
template<> struct array_bounds_check_ordinal_is_less<7,2> {};
template<> struct array_bounds_check_ordinal_is_less<7,3> {};
template<> struct array_bounds_check_ordinal_is_less<7,4> {};
template<> struct array_bounds_check_ordinal_is_less<7,5> {};
template<> struct array_bounds_check_ordinal_is_less<7,6> {};

template<> struct array_bounds_check_ordinal_is_less<6,0> {};
template<> struct array_bounds_check_ordinal_is_less<6,1> {};
template<> struct array_bounds_check_ordinal_is_less<6,2> {};
template<> struct array_bounds_check_ordinal_is_less<6,3> {};
template<> struct array_bounds_check_ordinal_is_less<6,4> {};
template<> struct array_bounds_check_ordinal_is_less<6,5> {};

template<> struct array_bounds_check_ordinal_is_less<5,0> {};
template<> struct array_bounds_check_ordinal_is_less<5,1> {};
template<> struct array_bounds_check_ordinal_is_less<5,2> {};
template<> struct array_bounds_check_ordinal_is_less<5,3> {};
template<> struct array_bounds_check_ordinal_is_less<5,4> {};

template<> struct array_bounds_check_ordinal_is_less<4,0> {};
template<> struct array_bounds_check_ordinal_is_less<4,1> {};
template<> struct array_bounds_check_ordinal_is_less<4,2> {};
template<> struct array_bounds_check_ordinal_is_less<4,3> {};

template<> struct array_bounds_check_ordinal_is_less<3,0> {};
template<> struct array_bounds_check_ordinal_is_less<3,1> {};
template<> struct array_bounds_check_ordinal_is_less<3,2> {};

template<> struct array_bounds_check_ordinal_is_less<2,0> {};
template<> struct array_bounds_check_ordinal_is_less<2,1> {};

template<> struct array_bounds_check_ordinal_is_less<1,0> {};

//----------------------------------------------------------------------

template< class T1 >
struct ArrayDimTagList<T1,void,void,void,void,void,void,void> {
  enum { Rank = 1 };
};

template< class T1 , class T2 >
struct ArrayDimTagList<T1,T2,void,void,void,void,void,void> {
  enum { Rank = 2 };
};

template< class T1 , class T2 , class T3 >
struct ArrayDimTagList<T1,T2,T3,void,void,void,void,void> {
  enum { Rank = 3 };
};

template< class T1 , class T2 , class T3 , class T4 >
struct ArrayDimTagList<T1,T2,T3,T4,void,void,void,void> {
  enum { Rank = 4 };
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 >
struct ArrayDimTagList<T1,T2,T3,T4,T5,void,void,void> {
  enum { Rank = 5 };
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 >
struct ArrayDimTagList<T1,T2,T3,T4,T5,T6,void,void> {
  enum { Rank = 6 };
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 >
struct ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,void> {
  enum { Rank = 7 };
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagList {
  enum { Rank = 8 };
};

//----------------------------------------------------------------------

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 0 > {
  typedef T1 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 1 > {
  typedef T2 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 2 > {
  typedef T3 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 3 > {
  typedef T4 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 4 > {
  typedef T5 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 5 > {
  typedef T6 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 6 > {
  typedef T7 type ;
};

template< class T1 , class T2 , class T3 , class T4 ,
          class T5 , class T6 , class T7 , class T8 >
struct ArrayDimTagListAt< ArrayDimTagList<T1,T2,T3,T4,T5,T6,T7,T8> , 7 > {
  typedef T8 type ;
};

//----------------------------------------------------------------------

namespace {

template< class T >
inline
const ArrayDimTag * array_tag_descriptor() { return & T::descriptor() ; }

template<>
inline
const ArrayDimTag * array_tag_descriptor<void>() { return NULL ; }

template< class TagList >
const ArrayDimTag ** array_tags()
{
  static const ArrayDimTag * t[8] =
    {
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,0>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,1>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,2>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,3>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,4>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,5>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,6>::type >() ,
      array_tag_descriptor< typename ArrayDimTagListAt<TagList,7>::type >()
    };

  return t ;
}

}

/** \endcond */
//----------------------------------------------------------------------

} // namespace phdmesh

#endif

