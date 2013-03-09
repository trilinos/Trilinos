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

template< class ViewTraits , class ValueType , unsigned N0, unsigned N1, class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , ValueType , LayoutTileLeft<N0,N1,true> , 2 , MemorySpace , MemoryTraits , void >
{ typedef LayoutTileLeftFast type ; };

template< class ViewTraits , class ValueType , unsigned N0, unsigned N1, class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , ValueType , LayoutTileLeft<N0,N1,false> , 2 , MemorySpace , MemoryTraits , void >
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

private:

  template< class DT , class DL , class DD , class DM >
  inline
  void allocate( View<DT,DL,DD,DM,LayoutTileLeftFast> & dst , const std::string label )
  {
    typedef View<DT,DL,DD,DM,LayoutTileLeftFast>  DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    const size_t count = allocation_count( dst );

    dst.m_ptr_on_device = (typename DstViewType::value_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::value_type) ,
                              sizeof(typename DstViewType::value_type) ,
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

    dst.m_shape.N0 = n0 ;
    dst.m_shape.N1 = n1 ;
    dst.m_tile_N0  = ( n0 + DstViewType::MASK_0 ) >> DstViewType::SHIFT_0 ;

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
    dst.m_shape   = src.m_shape ;
    dst.m_tile_N0 = src.m_tile_N0 ;
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

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape, src.m_shape.N0 , src.m_shape.N1 );

    dst.m_tile_N0       = src.m_tile_N0 ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutDefault , LayoutTileLeftFast, void >
{
  /** \brief Extracting a single tile from a tiled view */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutTileLeftFast> & src ,
                  const unsigned i0 ,
                  const typename enable_if<(
                    is_same< View<DT,DL,DD,DM,LayoutDefault> ,
                             typename View<ST,SL,SD,SM,LayoutTileLeftFast>::tile_type >::value
                  ), unsigned >::type i1 )
  {
    typedef View<DT,DL,DD,DM,LayoutDefault> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    enum { N0 = SL::N0 };
    enum { N1 = SL::N1 };
    enum { SHIFT_0 = power_of_two<N0>::value };
    enum { MASK_0 = N0 - 1 };
    enum { SHIFT_1 = power_of_two<N1>::value };

    const unsigned NT0 = ( src.dimension_0() + MASK_0 ) >> SHIFT_0 ;

    dst.m_ptr_on_device = src.m_ptr_on_device + (( i0 + i1 * NT0 ) << ( SHIFT_0 + SHIFT_1 ));
    dst.m_stride        = N0 ;

    ViewTracking< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType , class DeviceType , class MemoryTraits >
class View< DataType , LayoutType , DeviceType , MemoryTraits , Impl::LayoutTileLeftFast >
  : public ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits >
{
private:
  template< class , class , class > friend class Impl::ViewAssignment ;

  typedef ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits > traits ;

  typedef Impl::ViewAssignment<Impl::LayoutTileLeftFast> alloc ;

  typedef Impl::ViewAssignment<Impl::LayoutTileLeftFast,
                               Impl::LayoutTileLeftFast> assign ;

  typename traits::value_type * m_ptr_on_device ;
  typename traits::shape_type   m_shape ;
  unsigned                      m_tile_N0 ;

  typedef typename traits::array_layout layout ;

  enum { SHIFT_0 = Impl::power_of_two<layout::N0>::value };
  enum { SHIFT_1 = Impl::power_of_two<layout::N1>::value };
  enum { MASK_0  = layout::N0 - 1 };
  enum { MASK_1  = layout::N1 - 1 };

public:

  typedef Impl::LayoutTileLeftFast specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::layout_type ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::layout_type ,
                Host > HostMirror ;

  enum { Rank = 2 };

  KOKKOSARRAY_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_2() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_3() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_4() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_5() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_6() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_7() const { return 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  View() : m_ptr_on_device(0) {}

  KOKKOSARRAY_INLINE_FUNCTION
  ~View() { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0) { assign( *this , rhs ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign( *this , rhs ); return *this ; }

  //------------------------------------
  // Array allocator and member access operator:

  View( const std::string & label , const size_t n0 , const size_t n1 )
    : m_ptr_on_device(0) { alloc( *this , label , n0 , n1 ); }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );

      // Use care to insert necessary parentheses as the
      // shift operators have lower precedence than the arithmatic operators.

      return m_ptr_on_device[
        // ( ( Tile offset                               ) *  ( Tile size       ) )
         + ( ( (i0>>SHIFT_0) + m_tile_N0 * (i1>>SHIFT_1) ) << (SHIFT_0 + SHIFT_1) )
        // ( Offset within tile                       )
         + ( (i0 & MASK_0) + ((i1 & MASK_1)<<SHIFT_0) ) ] ;
    }

  //------------------------------------
  // Tile specialization specific declarations and functions:

  typedef View< typename traits::value_type [ layout::N0 ][ layout::N1 ] ,
                LayoutLeft ,
                typename traits::device_type ,
                MemoryUnmanaged >
    tile_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_t tiles_in_dimension_0() const { return m_tile_N0 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_t tiles_in_dimension_1() const { return ( m_shape.N1 + MASK_1 ) >> SHIFT_1 ; }


  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_t global_to_tile_index_0( const iType & global_i0 ) const
    { return global_i0 >> SHIFT_0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_t global_to_tile_index_1( const iType & global_i1 ) const
    { return global_i1 >> SHIFT_1 ; }


  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_t global_to_local_tile_index_0( const iType & global_i0 ) const
    { return global_i0 & MASK_0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_t global_to_local_tile_index_1( const iType & global_i1 ) const
    { return global_i1 & MASK_1 ; }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWTILELEFT_HPP */

