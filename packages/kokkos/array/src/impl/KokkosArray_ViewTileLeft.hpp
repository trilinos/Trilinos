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

#ifndef KOKKOSARRAY_VIEWTILELEFT_HPP
#define KOKKOSARRAY_VIEWTILELEFT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

struct LayoutTileLeftFast ;
struct LayoutTileLeftSlow ;

template< class ViewTraits , unsigned N0, unsigned N1, class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , LayoutTileLeft<N0,N1,true> , 2 , MemorySpace , MemoryTraits , void >
{ typedef LayoutTileLeftFast type ; };

template< class ViewTraits , unsigned N0, unsigned N1, class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , LayoutTileLeft<N0,N1,false> , 2 , MemorySpace , MemoryTraits , void >
{ typedef LayoutTileLeftSlow type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutTileLeftFast , void , void >
{
  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  size_t allocation_count( const View<T,L,D,M,LayoutTileLeftFast> & dst )
  {
    typedef typename ViewTraits<T,L,D,M>::array_layout layout ;

    return
      layout::N0 * layout::N1 * 
         ( ( dst.m_shape.N0 + layout::N0 - 1 ) / layout::N0 ) *
         ( ( dst.m_shape.N1 + layout::N1 - 1 ) / layout::N1 );
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void decrement( View<T,L,D,M,LayoutTileLeftFast> & dst )
  {
    typedef View<T,L,D,M,LayoutTileLeftFast> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::decrement( dst.m_ptr_on_device );
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void increment( View<T,L,D,M,LayoutTileLeftFast> & dst )
  {
    typedef View<T,L,D,M,LayoutTileLeftFast> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::increment( dst.m_ptr_on_device );
  }

private:

  template< class DT , class DL , class DD , class DM >
  inline
  void allocate( View<DT,DL,DD,DM,LayoutTileLeftFast> & dst , const std::string label )
  {
    typedef View<DT,DL,DD,DM,LayoutTileLeftFast>  DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;

    decrement( dst );

    const size_t count = allocation_count( dst );

    dst.m_ptr_on_device = (typename DstViewType::scalar_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::scalar_type) ,
                              sizeof(typename DstViewType::scalar_type) ,
                              count );

    ViewInitialize< DstViewType >::apply( dst );
  }

public:

  template< class DT , class DL , class DD , class DM >
  inline
  ViewAssignment( View<DT,DL,DD,DM,LayoutTileLeftFast> & dst ,
                  const typename enable_if< ViewTraits<DT,DL,DD,DM>::is_managed , std::string >::type & label ,
                  const size_t n0 ,
                  const size_t n1 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 )
  {
    typedef View<DT,DL,DD,DM,LayoutTileLeftFast>  DstViewType ;

    dst.m_stride   = 1 ;
    dst.m_shape.N0 = n0 ;
    dst.m_shape.N1 = n1 ;

    allocate( dst , label );
  }


  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,LayoutTileLeftFast> & dst ,
                  const View<ST,SL,SD,SM,LayoutTileLeftFast> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,LayoutTileLeftFast> ,
                             typename View<ST,SL,SD,SM,LayoutTileLeftFast>::HostMirror >::value
                  >::type * = 0 )
  {
    dst.m_shape  = src.m_shape ;
    dst.m_stride = src.m_stride ;
    allocate( dst , "mirror" );
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutTileLeftFast , LayoutTileLeftFast, void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutTileLeftFast> & dst ,
                  const View<ST,SL,SD,SM,LayoutTileLeftFast> & src ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                     typename ViewTraits<ST,SL,SD,SM>::shape_type >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,LayoutTileLeftFast> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewAssignment< LayoutTileLeftFast >::decrement( dst );

    shape_type::assign( dst.m_shape, src.m_shape.N0 , src.m_shape.N1 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewAssignment< LayoutTileLeftFast >::increment( dst );
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutLeft , LayoutTileLeftFast, void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutLeft> & dst ,
                  const View<ST,SL,SD,SM,LayoutTileLeftFast> & src ,
                  const unsigned i0 ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 0 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::shape_type::N0) == unsigned(SL::N0) )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::shape_type::N1) == unsigned(SL::N1) )
                  ), unsigned >::type i1 )
  {
    typedef View<DT,DL,DD,DM,LayoutLeft> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewAssignment< LayoutLeft >::decrement( dst );

    enum { N0 = SL::N0 };
    enum { N1 = SL::N1 };
    enum { SHIFT_0 = power_of_two<N0>::value };
    enum { MASK_0 = N0 - 1 };
    enum { SHIFT_1 = power_of_two<N1>::value };

    const unsigned NT0 = ( src.dimension_0() + MASK_0 ) >> SHIFT_0 ;

    dst.m_ptr_on_device = src.m_ptr_on_device + (( i0 + i1 * NT0 ) << ( SHIFT_0 + SHIFT_1 ));
    dst.m_stride        = N0 ;

    ViewAssignment< LayoutLeft >::increment( dst );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class DT , unsigned N0 , unsigned N1 , class DD , class DM >
struct ViewAssignment<
  View<DT,LayoutTileLeft<N0,N1,true>,DD,DM,LayoutTileLeftFast> ,
  typename DD::memory_space ,
  typename enable_if< (
    ViewTraits<DT,LayoutTileLeft<N0,N1,true>,DD,DM>::is_managed
  ) >::type >
{
  typedef View<DT,LayoutTileLeft<N0,N1,true>,DD,DM,LayoutTileLeftFast> DstViewType ;
  typedef typename DstViewType::shape_type shape_type ;

  KOKKOSARRAY_INLINE_FUNCTION static
  size_t allocation_count( const DstViewType & dst )
  { return ViewAssignment< LayoutTileLeftFast >::allocation_count( dst ); }

  // Same data type, same layout, different device; used to create a mirror.
  template< class D , class M >
  ViewAssignment( DstViewType & dst , const View< typename DstViewType::data_type ,
                                                  typename DstViewType::layout_type ,
                                                  D , M > & src )
  { ViewAssignment< LayoutTileLeftFast >( dst , src ); }

  ViewAssignment( DstViewType & dst , const std::string & label ,
                  const size_t n0 , const size_t n1 ,
                  const size_t ,
                  const size_t ,
                  const size_t ,
                  const size_t ,
                  const size_t ,
                  const size_t )
  { ViewAssignment< LayoutTileLeftFast >( dst , label , n0 , n1 ); }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWTILELEFT_HPP */

