/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOSARRAY_SHAPELEFT_HPP
#define KOKKOSARRAY_SHAPELEFT_HPP

#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template < unsigned ValueSize >
struct ShapeMap< Shape<ValueSize,0> , LayoutLeft >
{
  typedef Shape<ValueSize,0> shape_type ;

  static inline
  size_t offset( const shape_type & , const unsigned ) { return 0 ; }

  template< class MemorySpace >
  static inline
  size_t stride( const shape_type & ) { return 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  size_t allocation_count( const shape_type & , unsigned ) { return 1 ; }
};

template < unsigned ValueSize , unsigned s0 >
struct ShapeMap< Shape<ValueSize,1,s0> , LayoutLeft >
{
  typedef Shape<ValueSize,1,s0> shape_type ;

  static inline
  size_t offset( const shape_type & shape , const unsigned , const size_t i0 )
  { assert_shape_bounds( shape, i0 ); return i0 ; }

  template< class MemorySpace >
  static inline
  size_t stride( const shape_type & ) { return 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  size_t allocation_count( const shape_type & shape , unsigned ) { return shape.N0 ; }
};

template < unsigned ValueSize , unsigned Rank ,
           unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
           unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
struct ShapeMap< Shape<ValueSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> , LayoutLeft >
{
  typedef Shape<ValueSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> shape_type ;

  static inline
  size_t offset(
    const shape_type & shape , const unsigned stride ,
    const size_t i0     , const size_t i1     ,
    const size_t i2 = 0 , const size_t i3 = 0 ,
    const size_t i4 = 0 , const size_t i5 = 0 ,
    const size_t i6 = 0 , const size_t i7 = 0 )
  {
    assert_shape_bounds( shape, i0, i1, i2, i3, i4, i5, i6, i7 );

    return i0 + stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * (
           i5 + shape.N5 * (
           i6 + shape.N6 * i7 ))))));
  }

  template< class MemorySpace >
  static inline
  size_t stride( const shape_type & shape )
  {
    return MemorySpace::preferred_alignment( shape.scalar_size , shape.N0 );
  }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  size_t allocation_count( const shape_type & shape , unsigned stride )
  {
    return stride * shape.N1 * shape.N2 * shape.N3 *
           shape.N4 * shape.N5 * shape.N6 * shape.N7 ;
  }
};

//----------------------------------------------------------------------------
/** \brief  Subset of a multivector */

template< unsigned ValueSize , unsigned s1 >
struct SubShape< LayoutLeft , Shape< ValueSize, 1, 0 > ,
                 LayoutLeft , Shape< ValueSize, 2, 0,s1 > >
{
  typedef Shape< ValueSize, 1, 0 >    DstShape ;
  typedef Shape< ValueSize, 2, 0,s1 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride , const size_t i1 )
  {
    assert_shape_bounds( src , 0 , i1 );
    offset = ShapeMap<SrcShape,LayoutLeft>::offset( src , src_stride , 0 , i1 );
    shape.N0     = src.N0 ;
    stride       = 1 ;
  }
};

template< unsigned ValueSize , unsigned s1 >
struct SubShape< LayoutRight , Shape< ValueSize, 1, 0 > ,
                 LayoutLeft ,  Shape< ValueSize, 2, 0,s1 > >
{
  typedef Shape< ValueSize, 1, 0 >    DstShape ;
  typedef Shape< ValueSize, 2, 0,s1 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride , const size_t i1 )
  {
    assert_shape_bounds( src , 0 , i1 );
    offset = ShapeMap<SrcShape,LayoutLeft>::offset( src , src_stride , 0 , i1 );
    shape.N0     = src.N0 ;
    stride       = 1 ;
  }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_SHAPELEFT_HPP */

