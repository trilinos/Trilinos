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

#ifndef KOKKOS_VIEWOFFSET_HPP
#define KOKKOS_VIEWOFFSET_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Shape.hpp>

namespace Kokkos { namespace Impl {

template < class ShapeType , class LayoutType , typename Enable = void>
struct ViewOffset ;

//----------------------------------------------------------------------------

template< class T , unsigned R , typename Enable = void >
struct AssignViewOffsetDimension ; 

template< class T >
struct AssignViewOffsetDimension< T , 0 , typename enable_if<( 0 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N0 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 1 , typename enable_if<( 1 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N1 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 2 , typename enable_if<( 2 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N2 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 3 , typename enable_if<( 3 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N3 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 4 , typename enable_if<( 4 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N4 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 5 , typename enable_if<( 5 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N5 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 6 , typename enable_if<( 6 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N6 = n ; } };

template< class T >
struct AssignViewOffsetDimension< T , 7 , typename enable_if<( 7 < T::rank_dynamic )>::type >
{ AssignViewOffsetDimension( T & s , unsigned n ) { s.N7 = n ; } };

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 >= rank OR 0 == rank_dynamic ) : has padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutLeft
                 , typename enable_if<( 1 >= ShapeType::rank
                                        ||
                                        0 == ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t     size_type ;
  typedef ShapeType  shape_type ;
  typedef LayoutLeft array_layout ;

  enum { has_padding = false };

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { AssignViewOffsetDimension< ViewOffset , R >( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( 1 == int(ShapeRHS::rank)
                                    &&
                                    1 == int(shape_type::rank)
                                    &&
                                    1 == int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  KOKKOS_INLINE_FUNCTION
  void set_padding() {}

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      if ( 0 < shape_type::rank ) { s[0] = 1 ; }
      if ( 1 < shape_type::rank ) { s[1] = shape_type::N0 ; }
      if ( 2 < shape_type::rank ) { s[2] = s[1] * shape_type::N1 ; }
      if ( 3 < shape_type::rank ) { s[3] = s[2] * shape_type::N2 ; }
      if ( 4 < shape_type::rank ) { s[4] = s[3] * shape_type::N3 ; }
      if ( 5 < shape_type::rank ) { s[5] = s[4] * shape_type::N4 ; }
      if ( 6 < shape_type::rank ) { s[6] = s[5] * shape_type::N5 ; }
      if ( 7 < shape_type::rank ) { s[7] = s[6] * shape_type::N6 ; }
    }

  KOKKOS_INLINE_FUNCTION size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_1() const { return shape_type::N0 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_2() const { return shape_type::N0 * shape_type::N1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_3() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N7 ; }

  // rank 1
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    { return i0 + shape_type::N0 * i1 ; }

  //rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0
                      , I1 const& i1
                      , I2 const& i2
                      ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * i2 );
    }

  //rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * i3 ));
    }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * i4 )));
    }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * (
             i4 + shape_type::N4 * i5 ))));
    }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6) const
  {
    return i0 + shape_type::N0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7) const
  {
    return i0 + shape_type::N0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * (
           i6 + shape_type::N6 * i7 ))))));
  }
};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutLeft
                 , typename enable_if<( 1 < ShapeType::rank
                                        &&
                                        0 < ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t     size_type ;
  typedef ShapeType  shape_type ;
  typedef LayoutLeft array_layout ;

  enum { has_padding = true };

  unsigned S0 ;

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { AssignViewOffsetDimension< ViewOffset , R >( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); S0 = shape_type::N0 ; }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      S0 = shape_type::N0 ;
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) > 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      S0 = rhs.S0 ;
    }

  KOKKOS_INLINE_FUNCTION
  void set_padding()
    {
      enum { div   = MEMORY_ALIGNMENT / shape_type::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % shape_type::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < S0 ) {

        const unsigned count_mod = S0 % ( div ? div : 1 );

        if ( count_mod ) { S0 += align - count_mod ; }
      }
    }

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      if ( 0 < shape_type::rank ) { s[0] = 1 ; }
      if ( 1 < shape_type::rank ) { s[1] = S0 ; }
      if ( 2 < shape_type::rank ) { s[2] = s[1] * shape_type::N1 ; }
      if ( 3 < shape_type::rank ) { s[3] = s[2] * shape_type::N2 ; }
      if ( 4 < shape_type::rank ) { s[4] = s[3] * shape_type::N3 ; }
      if ( 5 < shape_type::rank ) { s[5] = s[4] * shape_type::N4 ; }
      if ( 6 < shape_type::rank ) { s[6] = s[5] * shape_type::N5 ; }
      if ( 7 < shape_type::rank ) { s[7] = s[6] * shape_type::N6 ; }
    }

  KOKKOS_INLINE_FUNCTION size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_1() const { return S0 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_2() const { return S0 * shape_type::N1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_3() const { return S0 * shape_type::N1 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N7 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 , I1 const & i1) const
    { return i0 + S0 * i1 ; }

  //rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * i2 );
    }

  //rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * i3 ));
    }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * i4 )));
    }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * (
             i4 + shape_type::N4 * i5 ))));
    }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i0 + S0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i0 + S0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * (
           i6 + shape_type::N6 * i7 ))))));
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 >= rank OR 1 >= rank_dynamic ) : no padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutRight
                 , typename enable_if<( 1 >= ShapeType::rank
                                        ||
                                        1 >= ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t       size_type;
  typedef ShapeType    shape_type;
  typedef LayoutRight  array_layout ;

  enum { has_padding = false };

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { AssignViewOffsetDimension< ViewOffset , R >( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( 1 == int(ShapeRHS::rank)
                                    &&
                                    1 == int(shape_type::rank)
                                    &&
                                    1 == int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  KOKKOS_INLINE_FUNCTION
  void set_padding() {}

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  size_type stride_R() const
    {
      return shape_type::N1 * shape_type::N2 * shape_type::N3 *
             shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    };

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < shape_type::rank ) { s[7] = n ; n *= shape_type::N7 ; }
      if ( 6 < shape_type::rank ) { s[6] = n ; n *= shape_type::N6 ; }
      if ( 5 < shape_type::rank ) { s[5] = n ; n *= shape_type::N5 ; }
      if ( 4 < shape_type::rank ) { s[4] = n ; n *= shape_type::N4 ; }
      if ( 3 < shape_type::rank ) { s[3] = n ; n *= shape_type::N3 ; }
      if ( 2 < shape_type::rank ) { s[2] = n ; n *= shape_type::N2 ; }
      if ( 1 < shape_type::rank ) { s[1] = n ; n *= shape_type::N1 ; }
      if ( 0 < shape_type::rank ) { s[0] = n ; }
    }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const { return shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const { return shape_type::N7 * shape_type::N6 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_3() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_2() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_1() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_0() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 * shape_type::N1 ; }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i1 + shape_type::N1 * i0 ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 ));
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i3 + shape_type::N3 * (
             i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 )));
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i4 + shape_type::N4 * (
             i3 + shape_type::N3 * (
             i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 ))));
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
  {
    return i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 )))));
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 ))))));
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i7 + shape_type::N7 * (
           i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 )))))));
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 < rank AND 1 < rank_dynamic ) : has padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutRight
                 , typename enable_if<( 1 < ShapeType::rank
                                        &&
                                        1 < ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t       size_type;
  typedef ShapeType    shape_type;
  typedef LayoutRight  array_layout ;

  enum { has_padding = true };

  unsigned SR ;

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { AssignViewOffsetDimension< ViewOffset , R >( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
    {
      shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 );
      SR = shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      SR = shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) > 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      SR = rhs.SR ;
    }

  KOKKOS_INLINE_FUNCTION
  void set_padding()
    {
      enum { div   = MEMORY_ALIGNMENT / shape_type::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % shape_type::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < SR ) {

        const unsigned count_mod = SR % ( div ? div : 1 );

        if ( count_mod ) { SR += align - count_mod ; }
      }
    }

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const { return shape_type::N0 * SR ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < shape_type::rank ) { s[7] = n ; n *= shape_type::N7 ; }
      if ( 6 < shape_type::rank ) { s[6] = n ; n *= shape_type::N6 ; }
      if ( 5 < shape_type::rank ) { s[5] = n ; n *= shape_type::N5 ; }
      if ( 4 < shape_type::rank ) { s[4] = n ; n *= shape_type::N4 ; }
      if ( 3 < shape_type::rank ) { s[3] = n ; n *= shape_type::N3 ; }
      if ( 2 < shape_type::rank ) { s[2] = n ; n *= shape_type::N2 ; }
      if ( 1 < shape_type::rank ) { s[1] = n ; n *= shape_type::N1 ; }
      if ( 0 < shape_type::rank ) { s[0] = SR ; }
    }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const { return shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const { return shape_type::N7 * shape_type::N6 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_3() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_2() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_1() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_0() const { return SR ; }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i1 + i0 * SR ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i2 + shape_type::N2 * ( i1 ) +
             i0 * SR ;
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i3 + shape_type::N3 * (
             i2 + shape_type::N2 * ( i1 )) +
             i0 * SR ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i4 + shape_type::N4 * (
             i3 + shape_type::N3 * (
             i2 + shape_type::N2 * ( i1 ))) +
             i0 * SR ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
  {
    return i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 )))) +
           i0 * SR ;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 ))))) +
           i0 * SR ;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i7 + shape_type::N7 * (
           i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 )))))) +
           i0 * SR ;
  }
};

//----------------------------------------------------------------------------
// LayoutStride : 
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutStride
                 , typename enable_if<( 0 < ShapeType::rank )>::type >
  : public ShapeType
{
  typedef size_t        size_type;
  typedef ShapeType     shape_type;
  typedef LayoutStride  array_layout ;

  size_type S[ shape_type::rank ];

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { AssignViewOffsetDimension< ViewOffset , R >( *this , n ); }

  template< class ShapeRHS , class Layout >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset<ShapeRHS,Layout> & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank) )>::type * = 0 )
    { rhs.stride(S); }

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    {
      // Maximum stride * corresponding dimension
      int j = 0 ;
      size_type m = 0 ;
      for ( int i = 0 ; i < shape_type::rank ; ++i ) {
        if ( m < S[i] ) { m = S[j=i]; }
      }
      return m * dimension( *this , j );
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    { for ( int i = 0 ; i < shape_type::rank ; ++i ) { s[i] = S[i] ; } }

  // rank 1
  template <typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0 ) const
    {
      return i0 * S[0] ;
    }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i0 * S[0] + i1 * S[1] ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] ;
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] + i6 * S[6] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] + i6 * S[6] + i7 * S[7] ;
    }
};

}} // namespace Kokkos::Impl

#endif //KOKKOS_VIEWOFFSET_HPP

