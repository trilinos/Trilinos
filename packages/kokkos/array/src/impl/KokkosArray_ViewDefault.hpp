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

#ifndef KOKKOSARRAY_VIEWDEFAULT_HPP
#define KOKKOSARRAY_VIEWDEFAULT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

struct PhysicalLayout;
struct LayoutDefault ;

template< typename ScalarType , class Rank , class RankDynamic , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutLeft , Rank , RankDynamic ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutDefault type ; };

template< typename ScalarType , class Rank , class RankDynamic , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutRight , Rank , RankDynamic ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutDefault type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutDefault , void , void >
{
  typedef LayoutDefault Specialize ;

  // Allocation from dimensions:

  template< class T , class L , class D , class M >
  inline
  ViewAssignment( View<T,L,D,M,Specialize> & dst ,
                  const typename enable_if< ViewTraits<T,L,D,M>::is_managed , std::string >::type & label ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    typedef ViewTraits<T,L,D,M> traits ;
    typedef typename traits::shape_type   shape_type ;
    typedef typename traits::memory_space memory_space ;
    typedef typename traits::scalar_type  scalar_type ;

    enum { is_left = is_same< typename traits::array_layout , LayoutLeft >::value };

    ViewTracking< traits >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );

    dst.m_stride = traits::rank <= 1 ? ( 0 /* unused */ ) :
       memory_space::preferred_alignment( dst.m_shape.scalar_size ,
         ( is_left ? dst.m_shape.N0
                   : dst.m_shape.N1 * dst.m_shape.N2 * dst.m_shape.N3 * dst.m_shape.N4 *
                     dst.m_shape.N5 * dst.m_shape.N6 * dst.m_shape.N7 ) );

    dst.m_ptr_on_device = (scalar_type *)
      memory_space::allocate( label , typeid(scalar_type) , sizeof(scalar_type) , dst.capacity() );

    ViewInitialize< View<T,L,D,M,Specialize> >::apply( dst );
  }

  // Allocation of a mirror:

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,Specialize> ,
                             typename View<ST,SL,SD,SM,Specialize>::HostMirror >::value
                  >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits ;
    typedef typename traits::memory_space memory_space ;
    typedef typename traits::scalar_type  scalar_type ;

    dst.m_shape  = src.m_shape ;
    dst.m_stride = src.m_stride ;

    dst.m_ptr_on_device = (scalar_type *)
      memory_space::allocate( "mirror" , typeid(scalar_type) , sizeof(scalar_type) , dst.capacity() );

    ViewInitialize< View<DT,DL,DD,DM,Specialize> >::apply( dst );
  }
};

template<>
struct ViewAssignment< LayoutDefault , LayoutDefault , void >
{
  typedef LayoutDefault Specialize ;

  /** \brief  Assign view with compatible shape and same layout */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                             typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                    &&
                    ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                     typename ViewTraits<ST,SL,SD,SM>::shape_type >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,Specialize> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< DstViewType >::increment( dst.m_ptr_on_device );
  }

  /** \brief  Assign view with different layout and rank=1 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ( ! is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                                 typename ViewTraits<ST,SL,SD,SM>::array_layout >::value )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,Specialize> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape, src.m_shape.N0 );

    dst.m_stride        = 0 ; // Unused for rank=1
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< DstViewType >::increment( dst.m_ptr_on_device );
  }

  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    && 
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                  ), unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = src.m_shape.N0 ;
    dst.m_stride        = 0 ;
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_stride * i1 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }

  /** \brief  Extract Rank-1 array from Rank-1 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
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
    dst.m_stride        = 0 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_shape , range.first );
      assert_shape_bounds( src.m_shape , range.second - 1 );

      dst.m_shape.N0 = range.second - range.first ;
      dst.m_ptr_on_device = src.m_ptr_on_device + range.first ;

      ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
    }
  }
};

// Assign scalar from Array

template<>
struct ViewAssignment< LayoutScalar , LayoutDefault , void >
{
  //------------------------------------
  // Scalar from Rank=1

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                  ), unsigned >::type i0 )
  {
    typedef View<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=2

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape , i0 , i1 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device = src.m_ptr_on_device + i0 + src.m_stride * i1 ;
    }
    else {
      dst.m_ptr_on_device = src.m_ptr_on_device + i1 + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=3

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i2 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * i2 );
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
        i2 + src.m_shape.N2 * (
        i1 ) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=4

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i3 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * i3 ));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
        i3 + src.m_shape.N3 * (
        i2 + src.m_shape.N2 * (
        i1 )) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=5

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i4 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * i4 )));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=6

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i5 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * i5 ))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=7

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i6 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 + src.m_shape.N5 * i6 )))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i6 + src.m_shape.N6 * (
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=8

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i7 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6, i7 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
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
          i1 )))))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

template< class Traits , class OperLayout , unsigned OperRank ,
          typename iType0 = int , typename iType1 = int ,
          typename iType2 = int , typename iType3 = int ,
          typename iType4 = int , typename iType5 = int ,
          typename iType6 = int , typename iType7 = int ,
          bool Match = ( is_same< typename Traits::array_layout , OperLayout >::value &&
                         ( unsigned(Traits::rank) == OperRank ) &&
                         iType0(0) == 0 &&
                         iType1(0) == 0 &&
                         iType2(0) == 0 &&
                         iType3(0) == 0 &&
                         iType4(0) == 0 &&
                         iType5(0) == 0 &&
                         iType6(0) == 0 &&
                         iType7(0) == 0 ) >
struct EnableViewOper ;

template< class Traits ,
          class OperLayout , unsigned OperRank ,
          typename iType0 , typename iType1 ,
          typename iType2 , typename iType3 ,
          typename iType4 , typename iType5 ,
          typename iType6 , typename iType7 >
struct EnableViewOper< Traits , OperLayout , OperRank ,
                       iType0 , iType1 , iType2 , iType3 ,
                       iType4 , iType5 , iType6 , iType7 ,
                       true >
{ typedef typename Traits::value_type type ; };

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType , class DeviceType , class MemoryTraits >
class View< DataType , LayoutType , DeviceType , MemoryTraits , Impl::LayoutDefault >
  : public ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits >
{
public:

  typedef ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits > traits ;

private:

  template< class , class , class > friend class Impl::ViewAssignment ;

  typedef Impl::ViewAssignment<Impl::LayoutDefault> alloc ;
  typedef Impl::ViewAssignment<Impl::LayoutDefault,Impl::LayoutDefault> assign ;

  typename traits::value_type * m_ptr_on_device ;
  typename traits::shape_type   m_shape ;
  unsigned                      m_stride ;

public:

  typedef View< typename traits::const_data_type ,
                typename traits::layout_type ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::layout_type ,
                Host > HostMirror ;

  //------------------------------------
  // Shape

  enum { Rank = traits::rank };

  KOKKOSARRAY_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_shape , i ); }

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOSARRAY_INLINE_FUNCTION
  ~View() { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View() : m_ptr_on_device(0), m_stride(0) { traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0); }

  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0) { assign( *this , rhs ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign( *this , rhs ); return *this ; }

  //------------------------------------
  // Copy or assign compatible array:

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    { assign( *this , rhs ); return *this ; }

  //------------------------------------
  // Allocation.

  explicit
  View( const std::string & label ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    { alloc( *this, label, n0, n1, n2, n3, n4, n5, n6, n7 ); }

  //------------------------------------
  // Is not allocated

  KOKKOSARRAY_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  // Array member access operators enabled if
  // (1) a zero value of all argument types are compile-time comparable to zero
  // (2) the rank matches the number of arguments
  // (3) the memory space is valid for the access
  //------------------------------------
  // LayoutLeft:

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 1, iType0 >::type & operator[]
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 1, iType0 >::type & operator()
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * i1 ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * i6 ))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * (
                              i6 + m_shape.N6 * i7 )))))) ];
    }

  //------------------------------------
  // LayoutRight:

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 1, iType0 >::type & operator[]
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 1, iType0 >::type & operator()
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * (
                              i1 ) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 4, iType0, iType1, iType2, iType3
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i7 + m_shape.N7 * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))))) + i0 * m_stride ];
    }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to this view specialization.

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::size_type stride( iType * const s ) const
  {
    enum { is_left = Impl::is_same< typename traits::array_layout , LayoutLeft >::value };

    if ( 1 == Rank ) {
      s[0] = 1 ;
    }
    else if ( is_left ) {
      s[0] = 1 ;
      s[1] = m_stride ;
      for ( int i = 2 ; i < Rank ; ++i ) { s[i] = s[i-1] * dimension(i-1); }
    }
    else {
      s[0] = m_stride ;
      s[Rank-1] = 1 ;
      for ( int i = Rank - 2 ; 0 < i ; --i ) { s[i] = s[i+1] * dimension(i+1); }
    }
  }

  // Count of contiguously allocated data members including padding.
  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::size_type capacity() const
  {
    enum { is_left = Impl::is_same< typename traits::array_layout , LayoutLeft >::value };

    return Rank == 1 ? typename traits::size_type(m_shape.N0) : (
           is_left ? ( m_stride   * m_shape.N1 * m_shape.N2 * m_shape.N3 *
                       m_shape.N4 * m_shape.N5 * m_shape.N6 * m_shape.N7 )
                   : ( m_stride * m_shape.N0 ) );
  }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWDEFAULT_HPP */

