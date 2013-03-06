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

template< class ViewTraits , unsigned R , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , LayoutLeft , R , MemorySpace , MemoryTraits , 
 typename enable_if<( R > 1 )>::type >
{ typedef LayoutLeft type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutLeft , void , void >
{

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  size_t allocation_count( const View<T,L,D,M,LayoutLeft> & dst )
  {
    return
      dst.m_stride   * dst.m_shape.N1 * dst.m_shape.N2 * dst.m_shape.N3 *
      dst.m_shape.N4 * dst.m_shape.N5 * dst.m_shape.N6 * dst.m_shape.N7 ;
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void decrement( View<T,L,D,M,LayoutLeft> & dst )
  {
    typedef View<T,L,D,M,LayoutLeft> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::decrement( dst.m_ptr_on_device );
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void increment( View<T,L,D,M,LayoutLeft> & dst )
  {
    typedef View<T,L,D,M,LayoutLeft> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::increment( dst.m_ptr_on_device );
  }

private:

  template< class T , class L , class D , class M >
  inline
  void allocate( View<T,L,D,M,LayoutLeft> & dst , const std::string & label )
  {
    typedef View<T,L,D,M,LayoutLeft> DstViewType ;
    typedef typename DstViewType::scalar_type   scalar_type ;
    typedef typename DstViewType::memory_space  memory_space ;

    const size_t count = allocation_count( dst );

    dst.m_ptr_on_device = (scalar_type *)
      memory_space::allocate( label , typeid(scalar_type) , sizeof(scalar_type) , count );

    ViewInitialize< DstViewType >::apply( dst );
  }

public:

  template< class T , class L , class D , class M >
  inline
  ViewAssignment( View<T,L,D,M,LayoutLeft> & dst ,
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
    typedef View<T,L,D,M,LayoutLeft> DstViewType ;
    typedef typename DstViewType::scalar_type   scalar_type ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;

    decrement( dst );

    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );

    dst.m_stride =
      memory_space::preferred_alignment( dst.m_shape.scalar_size , dst.m_shape.N0 );

    allocate( dst , label );
  }

  template< class T , class L , class D , class M >
  inline
  ViewAssignment( View<T,L,D,M,LayoutLeft> & dst ,
                  const typename enable_if< ViewTraits<T,L,D,M>::is_managed , std::string >::type & label ,
                  const typename ViewTraits<T,L,D,M>::shape_type shape )
  {
    ViewAssignment( dst, label, shape.N0, shape.N1, shape.N2, shape.N3,
                                shape.N4, shape.N5, shape.N6, shape.N7 );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,LayoutLeft> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,LayoutLeft> ,
                             typename View<ST,SL,SD,SM,LayoutLeft>::HostMirror >::value
                  >::type * = 0 )
  {
    dst.m_shape  = src.m_shape ;
    dst.m_stride = src.m_stride ;
    allocate( dst , "mirror" );
  }
};

template<>
struct ViewAssignment< LayoutLeft , LayoutLeft , void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutLeft> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft> & src ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                     typename ViewTraits<ST,SL,SD,SM>::shape_type >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,LayoutLeft> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewAssignment< LayoutLeft >::decrement( dst );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewAssignment< LayoutLeft >::increment( dst );
  }
};

template<>
struct ViewAssignment< LayoutScalar , LayoutLeft , void >
{

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                  ) , unsigned >::type i1 )
  {
    assert_shape_bounds( src.m_shape , i0 , i1 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 + src.m_stride * i1 ;

    ViewAssignment< LayoutScalar >::increment( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                  ) , unsigned >::type i2 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * i2 );

    ViewAssignment< LayoutScalar >::increment( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                  ) , unsigned >::type i3 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2, i3 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * i3 ));

    ViewAssignment< LayoutScalar >::increment( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                  ) , unsigned >::type i4 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * i4 )));

    ViewAssignment< LayoutScalar >::increment( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                  ) , unsigned >::type i5 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * i5 ))));

    ViewAssignment< LayoutScalar >::increment( dst );
  }


  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                  ) , unsigned >::type i6 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * (
        i5 + src.m_shape.N5 * i6 )))));

    ViewAssignment< LayoutScalar >::increment( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutLeft>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 ,
                  const typename enable_if< (
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                  ) , unsigned >::type i7 )
  {
    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6, i7 );

    ViewAssignment< LayoutScalar >::decrement( dst );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
        i0 + src.m_stride * (
        i1 + src.m_shape.N1 * (
        i2 + src.m_shape.N2 * (
        i3 + src.m_shape.N3 * (
        i4 + src.m_shape.N4 * (
        i5 + src.m_shape.N5 * (
        i6 + src.m_shape.N6 * i7 ))))));

    ViewAssignment< LayoutScalar >::increment( dst );
  }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWLEFT_HPP */

