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

#ifndef KOKKOSARRAY_SHAPERIGHT_HPP
#define KOKKOSARRAY_SHAPERIGHT_HPP

#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class T , unsigned RankDynamic , unsigned Rank >
inline
size_t allocation_count( const Shape<LayoutRight,T,RankDynamic,Rank> & shape )
{
  return Rank < 1 ? 1 : shape.N0 * shape.Stride ;
}

template< class T , class MemorySpace >
struct ShapeMap< Shape<LayoutRight,T,0,0> , MemorySpace > {};

template< class T , unsigned RankDynamic , class MemorySpace >
struct ShapeMap< Shape<LayoutRight,T,RankDynamic,1> , MemorySpace >
{
  inline static
  size_t stride( const Shape<LayoutRight,T,RankDynamic,1> shape ) 
  { return 1 ; }
};

template< class T , unsigned RankDynamic , unsigned Rank , class MemorySpace >
struct ShapeMap< Shape<LayoutRight,T,RankDynamic,Rank> , MemorySpace >
{
  inline static
  size_t stride( const Shape<LayoutRight,T,RankDynamic,Rank> shape ) 
  {
    const size_t block_count = 
                      shape.N1 * (
      2 == Rank ? 1 : shape.N2 * (
      3 == Rank ? 1 : shape.N3 * (
      4 == Rank ? 1 : shape.N4 * (
      5 == Rank ? 1 : shape.N5 * (
      6 == Rank ? 1 : shape.N6 * (
      7 == Rank ? 1 : shape.N7 ))))));

    return MemorySpace::preferred_alignment( shape.value_size , block_count );
  }
};

//----------------------------------------------------------------------------

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,1> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,1> shape ,
                const size_t i0 )
    {
      assert_shape_bounds( shape, i0 );

      return i0 ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,2> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,2> shape ,
                const size_t i0 , const size_t i1 )
    {
      assert_shape_bounds( shape, i0, i1 );

      return i1 + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,3> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,3> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 )
    {
      assert_shape_bounds( shape, i0, i1, i2 );

      return i2 + shape.N2 * (
             i1 ) + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,4> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,4> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 , const size_t i3 )
    {
      assert_shape_bounds( shape, i0, i1, i2, i3 );

      return i3 + shape.N3 * (
             i2 + shape.N2 * (
             i1 )) + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,5> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,5> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 , const size_t i3 ,
                const size_t i4 )
    {
      assert_shape_bounds( shape, i0, i1, i2, i3, i4 );

      return i4 + shape.N4 * (
             i3 + shape.N3 * (
             i2 + shape.N2 * (
             i1 ))) + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,6> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,6> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 , const size_t i3 ,
                const size_t i4 , const size_t i5 )
    {
      assert_shape_bounds( shape, i0, i1, i2, i3, i4, i5 );

      return i5 + shape.N5 * (
             i4 + shape.N4 * (
             i3 + shape.N3 * (
             i2 + shape.N2 * (
             i1 )))) + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,7> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,7> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 , const size_t i3 ,
                const size_t i4 , const size_t i5 ,
                const size_t i6 )
    {
      assert_shape_bounds( shape, i0, i1, i2, i3, i4, i5, i6 );

      return i6 + shape.N6 * (
             i5 + shape.N5 * (
             i4 + shape.N4 * (
             i3 + shape.N3 * (
             i2 + shape.N2 * (
             i1 ))))) + i0 * shape.Stride ;
    }
};

template< class T , unsigned RankDynamic >
struct ShapeOffset< Shape<LayoutRight,T,RankDynamic,8> >
{
  inline static
  size_t apply( const Shape<LayoutRight,T,RankDynamic,8> shape ,
                const size_t i0 , const size_t i1 ,
                const size_t i2 , const size_t i3 ,
                const size_t i4 , const size_t i5 ,
                const size_t i6 , const size_t i7 )
    {
      assert_shape_bounds( shape, i0, i1, i2, i3, i4, i5, i6, i7 );

      return i7 + shape.N7 * (
             i6 + shape.N6 * (
             i5 + shape.N5 * (
             i4 + shape.N4 * (
             i3 + shape.N3 * (
             i2 + shape.N2 * (
             i1 )))))) + i0 * shape.Stride ;
    }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_SHAPERIGHT_HPP */

