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

#ifndef KOKKOSARRAY_VIEWLEFT_HPP
#define KOKKOSARRAY_VIEWLEFT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstViewType >
struct ViewAssignment<
  DstViewType ,
  typename DstViewType::memory_space ,
  typename enable_if< (
    ( is_same< typename DstViewType::array_layout , LayoutLeft >::value )
    &&
    ( is_same< typename DstViewType::memory_traits , MemoryManaged >::value )
  ) >::type >
{
  typedef typename DstViewType::shape_type shape_type ;

private:

  typedef typename DstViewType::memory_space  memory_space ;

  static inline
  void allocate( DstViewType & dst , const std::string & label )
  {
    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    const size_t allocation_count =
      dst.m_stride   * dst.m_shape.N1 * dst.m_shape.N2 * dst.m_shape.N3 *
      dst.m_shape.N4 * dst.m_shape.N5 * dst.m_shape.N6 * dst.m_shape.N7 ;

    dst.m_ptr_on_device = (typename DstViewType::scalar_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::scalar_type) ,
                              sizeof(typename DstViewType::scalar_type) ,
                              allocation_count );
  }

public:

  // Same data type, same layout, different device; used to create a mirror.
  template< class D , class M >
  ViewAssignment( DstViewType & dst , const View< typename DstViewType::data_type ,
                                                  typename DstViewType::layout_type ,
                                                  D , M > & src )
  {
    dst.m_shape = src.m_shape ;
    dst.m_stride = src.m_stride ;
    allocate( dst , "mirror" );
  }

  ViewAssignment( DstViewType & dst , const std::string & label , const shape_type shape )
  {
    dst.m_shape = shape ;
    dst.m_stride =
      memory_space::preferred_alignment( dst.m_shape.scalar_size , dst.m_shape.N0 );

    allocate( dst , label );
  }

  ViewAssignment( DstViewType & dst , const std::string & label ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
    dst.m_stride =
      memory_space::preferred_alignment( dst.m_shape.scalar_size , dst.m_shape.N0 );

    allocate( dst , label );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Subview assignments:

namespace KokkosArray {
namespace Impl {

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 1 )
  ) , unsigned_<1> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 )
  {
    assert_shape_bounds( src.shape() , i0 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 2 )
  ) , unsigned_<2> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 + src.m_stride * i1 ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 3 )
  ) , unsigned_<3> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * i2 );

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 4 )
  ) , unsigned_<4> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * i3 ));

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 5 )
  ) , unsigned_<5> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * i4 )));

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 6 )
  ) , unsigned_<6> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * i5 ))));

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 7 )
  ) , unsigned_<7> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 , const unsigned i6 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 , i6 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * (
        i5 + src.m_shape.N5 * i6 )))));

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 8 )
  ) , unsigned_<8> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 , const unsigned i6 , const unsigned i7 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * (
        i5 + src.m_shape.N5 * (
        i6 + src.m_shape.N6 * i7 ))))));

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_VIEWLEFT_HPP */

