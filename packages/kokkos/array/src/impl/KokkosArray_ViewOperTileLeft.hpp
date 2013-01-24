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

#ifndef KOKKOSARRAY_VIEWOPERTILELEFT_HPP
#define KOKKOSARRAY_VIEWOPERTILELEFT_HPP

#include <impl/KokkosArray_ShapeTileLeft.hpp>

namespace KokkosArray {
namespace Impl {

#if defined( KOKKOSARRAY_EXPRESSION_CHECK )

#else

#endif
//----------------------------------------------------------------------------

template< class MemorySpace ,
          typename ValueType , unsigned ValueSize , unsigned M, unsigned N, unsigned s0 , unsigned s1 >
class ViewOper< MemorySpace ,
                ValueType , Shape< LayoutTileLeft<M,N,false/*not power of 2*/>, ValueSize, 2, s0, s1 > >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                         * m_ptr_on_device ;
  Shape<LayoutTileLeft<M,N,false>,ValueSize,2,s0,s1> m_shape ;

public:

  enum { TILE_DIMENSION_0 = M, TILE_DIMENSION_1 = N };

  typedef KokkosArray::View< ValueType[M][N], LayoutLeft, MemorySpace, MemoryUnmanaged > tile_type;

  template< typename iType0, typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  tile_type tile( const iType0 & itile , const iType1 & jtile ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );

      //TODO: assert within tile bounds

      // Generate this tile 'subview' with the proper shape (stride == M).
      // A default constructed Left-shape will attempt to align the leading dimension
      // on a cache-line length.

      typedef Shape<LayoutLeft,ValueSize,2,M,N> tile_shape_type ;

      return tile_type( m_ptr_on_device + M*N * (itile + jtile * tiles_in_dimension_0()),
                        tile_shape_type::create_unpadded() );
    }


  KOKKOSARRAY_INLINE_FUNCTION
    size_t tiles_in_dimension_0() const
    {
      return (m_shape.N0 + M-1)/M;
    }

  KOKKOSARRAY_INLINE_FUNCTION
    size_t tiles_in_dimension_1() const
    {
      return (m_shape.N1 + N-1)/N;
    }


  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_tile_index_0(iType i) const
    {
      return i/M;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_tile_index_1(iType i) const
    {
      return i/N;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_local_tile_index_0(iType i) const
    {
      return i%M;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_local_tile_index_1(iType i) const
    {
      return i%N;
    }



  template< typename iType0, typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i , const iType1 & j ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i,j );

      return  *  (m_ptr_on_device                             //base ptr
                 + (M*N*((i/M) + (j/N)*((m_shape.N0+M-1)/M))  //tile
                 + (i%M) + ((j%N)*M)))                        //offset
              ;
    }
};


//----------------------------------------------------------------------------
template< class MemorySpace ,
          typename ValueType , unsigned ValueSize , unsigned M, unsigned N, unsigned s0 , unsigned s1 >
class ViewOper< MemorySpace ,
                ValueType , Shape< LayoutTileLeft<M,N,true/*power of 2*/>, ValueSize, 2, s0, s1 > >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                         * m_ptr_on_device ;
  Shape<LayoutTileLeft<M,N,true>,ValueSize,2,s0,s1> m_shape ;

  enum {
      M_SHIFT = Impl::power_of_2<M>::value
    , N_SHIFT = Impl::power_of_2<N>::value
    , M_MASK =  M-1
    , N_MASK =  N-1
  };

public:

  enum { TILE_DIMENSION_0 = M, TILE_DIMENSION_1 = N };

  typedef KokkosArray::View< ValueType[M][N], LayoutLeft, MemorySpace, MemoryUnmanaged > tile_type;

  template< typename iType0, typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  tile_type tile( const iType0 & itile , const iType1 & jtile ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );

      //TODO: assert within tile bounds

      // Generate this tile 'subview' with the proper shape (stride == M).
      // A default constructed Left-shape will attempt to align the leading dimension
      // on a cache-line length.

      typedef Shape<LayoutLeft,ValueSize,2,M,N> tile_shape_type ;

      return tile_type( m_ptr_on_device
                        + ( (itile + jtile * tiles_in_dimension_0() ) << ( M_SHIFT + N_SHIFT )),
                        tile_shape_type::create_unpadded() );
    }

  KOKKOSARRAY_INLINE_FUNCTION
    size_t tiles_in_dimension_0() const
    {
      return (m_shape.N0 + M_MASK) >> M_SHIFT;
    }

  KOKKOSARRAY_INLINE_FUNCTION
    size_t tiles_in_dimension_1() const
    {
      return (m_shape.N1 + N_MASK) >> N_SHIFT;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_tile_index_0(iType i) const
    {
      return i>>M_SHIFT;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_tile_index_1(iType i) const
    {
      return i>>N_SHIFT;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_local_tile_index_0(iType i) const
    {
      return i & M_MASK;
    }

  template <typename iType>
  KOKKOSARRAY_INLINE_FUNCTION
    size_t global_to_local_tile_index_1(iType i) const
    {
      return i & N_MASK;
    }

  template< typename iType0, typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i , const iType1 & j ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i,j );

      // Use care to insert necessary parentheses as the
      // shift operators have lower precedence than the arithmatic operators.

      return  *  (m_ptr_on_device                                                                       //base ptr
                 + ((((i>>M_SHIFT) + (j>>N_SHIFT)*((m_shape.N0+M_MASK)>>M_SHIFT)) << (M_SHIFT + N_SHIFT))  //tile
                 + (i & M_MASK) + ((j & N_MASK)<<M_SHIFT)))                                             //offset

              ;
    }
};
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWOPERTILELEFT_HPP */


