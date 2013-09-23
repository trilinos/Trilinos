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

#ifndef KOKKOS_VIEWDEFAULT_HPP
#define KOKKOS_VIEWDEFAULT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct ViewAssignment< LayoutDefault , LayoutDefault , void >
{
  typedef LayoutDefault Specialize ;

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-1 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                  ), unsigned >::type i0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;

    assert_shape_bounds( src.m_shape , 1 , i0 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-2 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape , 2 , i0 , i1 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device = src.m_ptr_on_device + i0 + src.m_stride.value * i1 ;
    }
    else {
      dst.m_ptr_on_device = src.m_ptr_on_device + i1 + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-3 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 3, i0, i1, i2 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 ));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i2 + src.m_shape.N2 * (
          i1 ) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-4 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 4, i0, i1, i2, i3 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 )));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-5 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 5, i0, i1, i2, i3, i4);

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 ))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-6 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 6, i0, i1, i2, i3, i4, i5);

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 )))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )))) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-7 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 7, i0, i1, i2, i3, i4, i5, i6 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 + src.m_shape.N5 * (
          i6 ))))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i6 + src.m_shape.N6 * (
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))))) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-8 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 ,
                  const unsigned i7 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    enum { is_left = is_same< typename src_traits::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 8, i0, i1, i2, i3, i4, i5, i6, i7 );

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride.value * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 + src.m_shape.N5 * (
          i6 + src.m_shape.N6 * i7 ))))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i7 + src.m_shape.N7 * (
          i6 + src.m_shape.N6 * (
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )))))) + i0 * src.m_stride.value ;
    }

    ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from range of Rank-1 array, either layout */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;
    typedef typename traits_type::shape_type shape_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = 0 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_shape , 1 , range.first );
      assert_shape_bounds( src.m_shape , 1 , range.second - 1 );

      dst.m_shape.N0 = range.second - range.first ;
      dst.m_ptr_on_device = src.m_ptr_on_device + range.first ;

      ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ), unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = src.m_shape.N0 ;
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_stride.value * i1 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutRight Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ), ALL >::type & )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = src.m_shape.N1 ;
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_stride.value * i0 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract LayoutRight Rank-N array from range of LayoutRight Rank-N array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank > 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic > 0 )
                  )>::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;
    typedef typename traits_type::shape_type shape_type ;
    typedef typename View<DT,DL,DD,DM,Specialize>::stride_type stride_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    shape_type ::assign( dst.m_shape, 0, 0, 0, 0, 0, 0, 0, 0 );
    stride_type::assign( dst.m_stride , 0 );
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_shape , 8 , range.first ,      0,0,0,0,0,0,0);
      assert_shape_bounds( src.m_shape , 8 , range.second - 1 , 0,0,0,0,0,0,0);

      shape_type::assign( dst.m_shape, range.second - range.first ,
                          src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                          src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

      stride_type::assign( dst.m_stride , src.m_stride.value );

      dst.m_ptr_on_device = src.m_ptr_on_device + range.first * src.m_stride.value ;

      ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType0,iType0> & range0 ,
                  const std::pair<iType1,iType1> & range1 ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;
    typedef typename traits_type::shape_type shape_type ;
    enum { left = is_same< typename traits_type::array_layout , LayoutLeft >::value };

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = 0 ;
    dst.m_shape.N1      = 0 ;
    dst.m_stride.value  = 0 ;
    dst.m_ptr_on_device = 0 ;

    if ( range0.first < range0.second && range1.first < range1.second ) {
      assert_shape_bounds( src.m_shape , 2 , range0.first , range1.first );
      assert_shape_bounds( src.m_shape , 2 , range0.second - 1 , range1.second - 1 );

      dst.m_shape.N0 = range0.second - range0.first ;
      dst.m_shape.N1 = range1.second - range1.first ;
      dst.m_stride   = src.m_stride ;

      if ( left ) {
        // operator: dst.m_ptr_on_device[ i0 + dst.m_stride * i1 ]
        dst.m_ptr_on_device = src.m_ptr_on_device + range0.first + dst.m_stride.value * range1.first ;
      }
      else {
        // operator: dst.m_ptr_on_device[ i0 * dst.m_stride + i1 ]
        dst.m_ptr_on_device = src.m_ptr_on_device + range0.first * dst.m_stride.value + range1.first ;
      }

      ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Deep copy data from compatible value type, layout, rank, and specialization.
   *          Check the dimensions and allocation lengths at runtime.
   */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline static
  void deep_copy( const View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename Impl::enable_if<(
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::scalar_type ,
                                   typename ViewTraits<ST,SL,SD,SM>::non_const_scalar_type >::value
                    &&
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                                   typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                  )>::type * = 0 )
  {
    typedef typename ViewTraits<DT,DL,DD,DM>::memory_space dst_memory_space ;
    typedef typename ViewTraits<ST,SL,SD,SM>::memory_space src_memory_space ;

    if ( dst.m_ptr_on_device != src.m_ptr_on_device ) {

      Impl::assert_shapes_are_equal( dst.m_shape , src.m_shape );

      const size_t nbytes = dst.m_shape.scalar_size * capacity( dst.m_shape , dst.m_stride );

      DeepCopy< dst_memory_space , src_memory_space >( dst.m_ptr_on_device , src.m_ptr_on_device , nbytes );
    }
  }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWDEFAULT_HPP */

