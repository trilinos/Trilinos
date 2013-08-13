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
//----------------------------------------------------------------------------

template< class ShapeType , class LayoutType , class Enable = void >
class LayoutStride ;

/* Arrays with no stride
 *   (1) rank <= 1
 *   (2a) LayoutLeft  OR
 *   (2b) LayoutRight
 */
template< class ShapeType , class LayoutType >
class LayoutStride< ShapeType , LayoutType ,
                    typename enable_if<(
                      ( ShapeType::rank <= 1 ) &&
                      ( is_same<LayoutType,LayoutLeft >::value ||
                        is_same<LayoutType,LayoutRight>::value )
                    )>::type >
{
public:

  enum { dynamic = false };
  enum { value = 0 };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const ShapeType & , const bool = true ) {}
};

/* Arrays with compile-time stride that is not padded:
 *   (1) rank > 1
 *   (2) rank_dynamic == 0
 *   (3) LayoutLeft OR LayoutRight
 */
template< class ShapeType , class LayoutType >
class LayoutStride< ShapeType , LayoutType ,
                    typename enable_if<(
                      ( is_same<LayoutType,LayoutLeft >::value ||
                        is_same<LayoutType,LayoutRight>::value ) &&
                      ( 1 <  ShapeType::rank ) &&
                      ( 0 == ShapeType::rank_dynamic )
                    )>::type >
{
private:

  enum { left = is_same< LayoutType , LayoutLeft >::value };

public:

  enum { dynamic = false };

  // Left layout arrays are aligned on the first dimension.
  // Right layout arrays are aligned on blocks of the 2-7th dimensions.

  enum { value = left ? ShapeType::N0
                      : ShapeType::N1 * ShapeType::N2 * ShapeType::N3 *
                        ShapeType::N4 * ShapeType::N5 * ShapeType::N6 * ShapeType::N7 };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const ShapeType & , const bool = true ) {}
};

/* Arrays with compile-time stride that is padded:
 *   (1) rank-N with 1 < N
 *   (2) LayoutRight and rank_dynamic == 1 
 */
template< class ShapeType, class LayoutType >
class LayoutStride< ShapeType , LayoutType ,
                    typename enable_if<(
                      is_same< LayoutType, LayoutRight>::value &&
                      ( 1 <  ShapeType::rank ) && 
                      ( 1 == ShapeType::dynamic_rank )
                    )>::type >
{
private:

  enum { div   = MEMORY_ALIGNMENT / ShapeType::scalar_size };
  enum { mod   = MEMORY_ALIGNMENT % ShapeType::scalar_size };
  enum { align = ( ShapeType::rank_dynamic && 0 == mod ) ? div : 0 };

  // Left layout arrays are aligned on the first dimension.
  // Right layout arrays are aligned on blocks of the 2-7th dimensions.
  enum { count = ShapeType::N1 * ShapeType::N2 * ShapeType::N3 *
                 ShapeType::N4 * ShapeType::N5 * ShapeType::N6 * ShapeType::N7 };

  enum { count_mod = count % ( div ? div : 1 ) };

public:

  enum { dynamic = false };
  enum { value = ( align && count_mod && ( MEMORY_ALIGNMENT_THRESHOLD * align < count ) )
               ? ( count + align - count_mod )
               : count };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned , const bool = true ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const ShapeType & ) {}
};

/* Otherwise array has runtime stride that is padded. */
template< class ShapeType , class LayoutType , class Enable >
class LayoutStride
{
public:

  enum { dynamic = true };

  unsigned value ;

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & stride , const unsigned n ) { stride.value = n ; }

  KOKKOS_INLINE_FUNCTION static
  void assign_no_padding( LayoutStride & vs , const ShapeType & sh )
    {
      enum { left = is_same< LayoutType , LayoutLeft >::value };

      // Left  layout arrays are aligned on the first dimension.
      // Right layout arrays are aligned on blocks of the 2-8th dimensions.
      vs.value = ShapeType::rank <= 1 ? 0 : (
                 left ? sh.N0
                      : sh.N1 * sh.N2 * sh.N3 * sh.N4 * sh.N5 * sh.N6 * sh.N7 );
    }

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & vs , const ShapeType & sh )
    {
      enum { div   = MEMORY_ALIGNMENT / ShapeType::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % ShapeType::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      assign_no_padding( vs , sh );

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < vs.value ) {

        const unsigned count_mod = vs.value % ( div ? div : 1 );

        if ( count_mod ) { vs.value += align - count_mod ; }
      }
    }
};

template< class ShapeType , class LayoutType >
KOKKOS_INLINE_FUNCTION
size_t capacity( const ShapeType & shape ,
                 const LayoutStride< ShapeType , LayoutType > & stride )
{
  enum { left = is_same< LayoutType , LayoutLeft >::value };

  return ShapeType::rank <= 1 ? size_t(shape.N0) : (
         left ? size_t( stride.value * shape.N1 * shape.N2 * shape.N3 * shape.N4 * shape.N5 * shape.N6 * shape.N7 )
              : size_t( stride.value * shape.N0 ));
}

template< typename iType , class ShapeType , class LayoutType >
KOKKOS_INLINE_FUNCTION
void stride( iType * const s , const ShapeType & shape ,
                               const LayoutStride< ShapeType , LayoutType > & stride )
{
  enum { rank = ShapeType::rank };
  enum { left = is_same< LayoutType , LayoutLeft >::value };

  if ( 0 < rank ) {
    if ( 1 == rank ) {
      s[0] = 1 ;
    }
    else if ( left ) {
      s[0] = 1 ;
      s[1] = stride.value ;
      if ( 2 < rank ) { s[2] = s[1] * shape.N1 ; }
      if ( 3 < rank ) { s[3] = s[2] * shape.N2 ; }
      if ( 4 < rank ) { s[4] = s[3] * shape.N3 ; }
      if ( 5 < rank ) { s[5] = s[4] * shape.N4 ; }
      if ( 6 < rank ) { s[6] = s[5] * shape.N5 ; }
      if ( 7 < rank ) { s[7] = s[6] * shape.N6 ; }
    }
    else {
      s[rank-1] = 1 ;
      if ( 7 < rank ) { s[6] = s[7] * shape.N7 ; }
      if ( 6 < rank ) { s[5] = s[6] * shape.N6 ; }
      if ( 5 < rank ) { s[4] = s[5] * shape.N5 ; }
      if ( 4 < rank ) { s[3] = s[4] * shape.N4 ; }
      if ( 3 < rank ) { s[2] = s[3] * shape.N3 ; }
      if ( 2 < rank ) { s[1] = s[2] * shape.N2 ; }
      s[0] = stride.value ;
    }
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  View tracking increment/decrement only happens when
 *          view memory is managed and executing in the host space.
 */
template< class ViewTraits ,
          class MemorySpace  = typename ViewTraits::memory_space ,
          class MemoryTraits = typename ViewTraits::memory_traits ,
          class ExecSpec     = Kokkos::ExecutionSpace >
struct ViewTracking {
  KOKKOS_INLINE_FUNCTION static void increment( const void * ) {}
  KOKKOS_INLINE_FUNCTION static void decrement( const void * ) {}
};

template< class ViewTraits , class MemorySpace , class MemoryTraits >
struct ViewTracking< ViewTraits , MemorySpace , MemoryTraits ,
          typename enable_if< MemoryTraits::Unmanaged == 0 , HostSpace >::type >
{
  KOKKOS_INLINE_FUNCTION static void increment( const void * ptr )
    { MemorySpace::increment( ptr ); }

  KOKKOS_INLINE_FUNCTION static void decrement( const void * ptr )
    { MemorySpace::decrement( ptr ); }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutDefault , void , void >
{
  typedef LayoutDefault Specialize ;

  // Allocate a managed view from dimensions:
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
    typedef View<T,L,D,M,Specialize> view_type ;
    typedef typename view_type::device_type   device_type ;
    typedef typename view_type::memory_space  memory_space ;
    typedef typename view_type::scalar_type   scalar_type ;
    typedef typename view_type::shape_type    shape_type ;
    typedef typename view_type::stride_type   stride_type ;

    ViewTracking< view_type >::decrement( dst.m_ptr_on_device );

    shape_type ::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
    stride_type::assign( dst.m_stride , dst.m_shape );

    const size_t cap = capacity( dst.m_shape , dst.m_stride );

    dst.m_ptr_on_device = (scalar_type *)
      memory_space::allocate( label , typeid(scalar_type) , sizeof(scalar_type) , cap );

    ViewInitialize< device_type > init( dst );
  }

  // Allocate using the dimensions and stride of the input array.
  // The specialization, layout, and rank must be compatible.
  // Primarily used for creating mirrors.
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  typename enable_if< (
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                             typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank_dynamic) >= unsigned(ViewTraits<ST,SL,SD,SM>::rank_dynamic) )
                  ) >::type * = 0 )
  {
    typedef View<DT,DL,DD,DM>  dst_view_type ;
    typedef View<ST,SL,SD,SM>  src_view_type ;

    typedef typename src_view_type::memory_space  src_memory_space ;
    typedef typename dst_view_type::device_type   dst_device_type ;
    typedef typename dst_view_type::memory_space  dst_memory_space ;
    typedef typename dst_view_type::scalar_type   dst_scalar_type ;
    typedef typename dst_view_type::shape_type    dst_shape_type ;
    typedef typename dst_view_type::stride_type   dst_stride_type ;

    ViewTracking< dst_view_type >::decrement( dst.m_ptr_on_device );

    if ( src.m_ptr_on_device ) {
      std::string label = src_memory_space::query_label( src.m_ptr_on_device );
      label.append("_mirror");

      dst_shape_type::assign( dst.m_shape, src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                                           src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

      dst_stride_type::assign( dst.m_stride , src.m_stride.value );

      const size_t cap = capacity( dst.m_shape , dst.m_stride );

      dst.m_ptr_on_device = (dst_scalar_type *)
        dst_memory_space::allocate( label , typeid(dst_scalar_type) , sizeof(dst_scalar_type) , cap );

      ViewInitialize< dst_device_type > init( dst );
    }
    else {
      dst_shape_type ::assign( dst.m_shape, 0, 0, 0, 0, 0, 0, 0, 0 );
      dst_stride_type::assign( dst.m_stride , 0 );
      dst.m_ptr_on_device = 0 ;
    }
  }

  // Creating unmanaged view from ptr and dimensions.
  // Must have a dynamic stride that can be set.
  template< class T , class L , class D, class M >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( View<T,L,D,M,Specialize> & dst ,
                  const typename enable_if<
                    (
                      ! ViewTraits<T,L,D,M>::is_managed
                      &&
                      LayoutStride< typename ViewTraits<T,L,D,M>::shape_type ,
                                    typename ViewTraits<T,L,D,M>::array_layout >::dynamic
                    ) , typename ViewTraits<T,L,D,M>::scalar_type >::type * ptr ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    typedef View<T,L,D,M> view_type ;
    typedef typename view_type::scalar_type   scalar_type ;
    typedef typename view_type::shape_type    shape_type ;
    typedef typename view_type::stride_type   stride_type ;

    shape_type ::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
    stride_type::assign_no_padding( dst.m_stride , dst.m_shape );

    dst.m_ptr_on_device = (scalar_type *) ptr;
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutDefault , LayoutDefault , void >
{
  typedef LayoutDefault Specialize ;

  /** \brief  Assign compatible view */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value ) >::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,Specialize> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;
    typedef typename View<DT,DL,DD,DM,Specialize>::stride_type stride_type ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< DstViewType >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from range of Rank-1 array, either layout */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
                  ), unsigned >::type i0 ,
                  const ALL & )
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
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
    typedef typename View<DT,DL,DD,DM,LayoutDefault>::stride_type stride_type ;

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
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
  void deep_copy( const View<DT,DL,DD,DM,LayoutDefault> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
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
/** \brief  Enable view parentheses operator for a given layout, rank,
 *          and arguments of an integral type.
 */
template< class Traits , class OperLayout , unsigned OperRank ,
          typename iType0 = int , typename iType1 = int ,
          typename iType2 = int , typename iType3 = int ,
          typename iType4 = int , typename iType5 = int ,
          typename iType6 = int , typename iType7 = int ,
          class Enable = void >
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
                       typename enable_if<
                         is_same< typename Traits::array_layout , OperLayout >::value &&
                         ( unsigned(Traits::rank) == OperRank ) &&
                         iType0(0) == 0 &&
                         iType1(0) == 0 &&
                         iType2(0) == 0 &&
                         iType3(0) == 0 &&
                         iType4(0) == 0 &&
                         iType5(0) == 0 &&
                         iType6(0) == 0 &&
                         iType7(0) == 0 >::type >
{ typedef typename Traits::value_type type ; };

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::LayoutDefault >
  : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:

  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

private:

  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef Impl::ViewAssignment<Impl::LayoutDefault> alloc ;

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  typename traits::value_type * m_ptr_on_device ;
  typename traits::shape_type   m_shape ;
  stride_type                   m_stride ;
  
public:

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                Host ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  enum { Rank = traits::rank };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_shape.N0
           * m_shape.N1
           * m_shape.N2
           * m_shape.N3
           * m_shape.N4
           * m_shape.N5
           * m_shape.N6
           * m_shape.N7
           ;
  }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_shape , i ); }

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0)
    {
      traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0);
      stride_type::assign(m_stride,0);
    }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0)
    {
      (void) Impl::ViewAssignment< typename traits::specialize ,
                                   typename traits::specialize >( *this , rhs );
    }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs )
    {
      (void) Impl::ViewAssignment< typename traits::specialize ,
                                   typename traits::specialize >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Copy or assign compatible array:

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    : m_ptr_on_device(0)
    {
      typedef View<RT,RL,RD,RM,typename traits::specialize> ViewRHS ;

      enum { OK = Impl::ViewAssignable< View , ViewRHS >::value };

      (void) Impl::ViewAssignment< typename traits::specialize ,
                                   typename traits::specialize ,
                                   typename Impl::enable_if< OK >::type >( *this , rhs );
    }

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    {
      typedef View<RT,RL,RD,RM,typename traits::specialize> ViewRHS ;

      enum { OK = Impl::ViewAssignable< View , ViewRHS >::value };

      (void) Impl::ViewAssignment< typename traits::specialize ,
                                   typename traits::specialize ,
                                   typename Impl::enable_if< OK >::type >( *this , rhs );

      return *this ;
    }

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
    { (void)alloc( *this, label, n0, n1, n2, n3, n4, n5, n6, n7 ); }

  // Create unmanaged View from ptr, can be called in functors, NOTE: no alignment padding is performed
  KOKKOS_INLINE_FUNCTION
  View( const typename traits::scalar_type* ptr,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    {
      (void)alloc( *this, ptr, n0, n1, n2, n3, n4, n5, n6, n7 );
    }
  //------------------------------------
  // Is not allocated

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  // Array member access operators enabled if
  // (1) a zero value of all argument types are compile-time comparable to zero
  // (2) the rank matches the number of arguments
  // (3) the memory space is valid for the access
  //------------------------------------
  // LayoutLeft, rank 1:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 1, iType0 >::type & operator[]
    ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 1, iType0 >::type & operator()
    ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 1, iType0 >::type & operator()
    ( const iType0 & i0 , const int , const int = 0 , const int = 0 ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  // LayoutLeft, rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * i1 ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const int , const int = 0 ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * i1 ];
    }

  // LayoutLeft, rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  // LayoutLeft, rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const int , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  // LayoutLeft, rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const int , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

  // LayoutLeft, rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const int , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

  // LayoutLeft, rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * i6 ))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * i6 ))))) ];
    }

  // LayoutLeft, rank 8:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutLeft, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * (
                              i6 + m_shape.N6 * i7 )))))) ];
    }

  //------------------------------------
  // LayoutRight, rank 1:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 1, iType0 >::type & operator[]
    ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 1, iType0 >::type & operator()
    ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 1, iType0 >::type & operator()
    ( const iType0 & i0 , const int , const int = 0 , const int = 0 ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  // LayoutRight, rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 2, iType0, iType1 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const int , const int = 0 ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride.value ];
    }

  // LayoutRight, rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
      const int = 0 , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * ( i1 ) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper< traits, LayoutRight, 3, iType0, iType1, iType2 >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * ( i1 ) + i0 * m_stride.value ];
    }

  // LayoutRight, rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 4, iType0, iType1, iType2, iType3
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 4, iType0, iType1, iType2, iType3
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const int , const int = 0 , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride.value ];
    }

  // LayoutRight, rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 5, iType0, iType1, iType2, iType3 , iType4
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const int , const int = 0 , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride.value ];
    }

  // LayoutRight, rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const int , const int = 0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 6, iType0, iType1, iType2, iType3 , iType4, iType5
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride.value ];
    }

  // LayoutRight, rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_stride.value ];
    }

  // LayoutRight, rank 8:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::EnableViewOper<
    traits, LayoutRight, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7
  >::type & operator()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i7 + m_shape.N7 * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))))) + i0 * m_stride.value ];
    }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::capacity( m_shape , m_stride ); }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWDEFAULT_HPP */

