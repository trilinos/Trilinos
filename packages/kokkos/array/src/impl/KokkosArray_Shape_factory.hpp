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

#ifndef KOKKOSARRAY_SHAPE_FACTORY_HPP
#define KOKKOSARRAY_SHAPE_FACTORY_HPP

#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class T , unsigned RankDynamic , unsigned Rank >
inline
size_t allocation_count( const Shape<LayoutLeft,T,RankDynamic,Rank> & shape )
{
  return ( 0 == Rank ? 1 : (
           1 == Rank ? shape.N0 : shape.Stride * shape.N1 * (
           2 == Rank ? 1 : shape.N2 * (
           3 == Rank ? 1 : shape.N3 * (
           4 == Rank ? 1 : shape.N4 * (
           5 == Rank ? 1 : shape.N5 * (
           6 == Rank ? 1 : shape.N6 * (
           7 == Rank ? 1 : shape.N7 ))))))));
}

template< class T , unsigned RankDynamic , unsigned Rank >
inline
size_t allocation_count( const Shape<LayoutRight,T,RankDynamic,Rank> & shape )
{
  return ( 0 == Rank ? 1 : shape.N0 * (
           1 == Rank ? 1 : shape.Stride ));
}

template< class ShapeType >
inline
size_t shape_right_block_count( const ShapeType & shape )
{
  enum { Rank = ShapeType::rank };

  return shape.N1 * ( 2 == Rank ? 1 : shape.N2 * (
                      3 == Rank ? 1 : shape.N3 * (
                      4 == Rank ? 1 : shape.N4 * (
                      5 == Rank ? 1 : shape.N5 * (
                      6 == Rank ? 1 : shape.N6 * (
                      7 == Rank ? 1 : shape.N7 ))))));
}

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,0>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,0> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;
    return shape ;
  }
};

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,1>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,1> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;
    return shape ;
  }
};

template < class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,Rank> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,1,1>, MemorySpace> {

  typedef Shape<LayoutLeft,T,1,1> output_type ;

  inline static
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    return shape ;
  }
};

template < class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,1,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,1,Rank> output_type ;

  inline static
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

//----------------------------------------------------------------------------

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,2,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,2,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,3,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,3,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,4,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,4,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,5,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,5,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,6,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,6,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,7,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,7,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,8,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,8,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;
    shape.N7 = n7 ;

    const size_t right_block_count = shape_right_block_count( shape );

    shape.Stride = 1 == right_block_count ? shape.N0 :
      MemorySpace::preferred_alignment( shape.value_size , shape.N0 );

    return shape ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,0,0>,MemorySpace> {

  typedef Shape<LayoutRight,T,0,0> output_type ;

  inline static
  output_type create()
  {
    return output_type();
  }
};

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,0,1>,MemorySpace> {

  typedef Shape<LayoutRight,T,0,1> output_type ;

  inline static
  output_type create()
  {
    return output_type();
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,0,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,0,Rank> output_type ;

  inline static 
  output_type create()
  {
    output_type shape ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,1,1>,MemorySpace> {

  typedef Shape<LayoutRight,T,1,1> output_type ;

  inline static 
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,1,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,1,Rank> output_type ;

  inline static 
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,2,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,2,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,3,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,3,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,4,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,4,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,5,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,5,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,6,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,6,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,7,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,7,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,8,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,8,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;
    shape.N7 = n7 ;

    const size_t block_count = shape_right_block_count( shape );

    shape.Stride = 1 == block_count ? 1 :
      MemorySpace::preferred_alignment( shape.value_size , block_count );

    return shape ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

