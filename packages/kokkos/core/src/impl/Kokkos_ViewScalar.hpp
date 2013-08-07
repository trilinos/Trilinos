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

#ifndef KOKKOS_VIEWSCALAR_HPP
#define KOKKOS_VIEWSCALAR_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Specialization for a Rank = 0 array */
struct LayoutScalar ;

template< typename ScalarType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutLeft , unsigned_<0> , unsigned_<0> ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutScalar type ; };

template< typename ScalarType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutRight , unsigned_<0> , unsigned_<0> ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutScalar type ; };

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Deep copy of scalar value in Host space */

template< typename ValueType , class Arg1 , class Arg2 , class Arg3 >
inline
void deep_copy( ValueType & dst ,
                const View< ValueType , Arg1 , Arg2 , Arg3 , Impl::LayoutScalar > & src ,
                typename Impl::enable_if< (
                  Impl::is_same< typename ViewTraits<ValueType,Arg1,Arg2,Arg3>::memory_space , HostSpace >::value 
                ) >::type * = 0 )
{ dst = src ; }

template< typename ValueType , class Arg1 , class Arg2 , class Arg3 >
inline
void deep_copy( const View< ValueType , Arg1 , Arg2 , Arg3 , Impl::LayoutScalar > & dst ,
                const ValueType & src ,
                typename Impl::enable_if< (
                  Impl::is_same< typename ViewTraits<ValueType,Arg1,Arg2,Arg3>::memory_space , HostSpace >::value
                ) >::type * = 0 )
{ dst = src ; }

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct ViewAssignment< LayoutScalar , void , void >
{
  template< class T , class L , class D , class M >
  ViewAssignment( View<T,L,D,M,LayoutScalar> & dst ,
                  typename enable_if< (
                    is_same< typename ViewTraits<T,L,D,M>::memory_traits ,
                             MemoryManaged >::value
                  ) , const std::string >::type & label )
  {
    typedef View<T,L,D,M,LayoutScalar> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = (typename DstViewType::value_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::value_type) ,
                              sizeof(typename DstViewType::value_type) ,
                              1 );

    ViewInitialize< typename DstViewType::device_type > init( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutScalar> & ,
                  typename enable_if< (
                    is_same< View<DT,DL,DD,DM,LayoutScalar> ,
                             typename View<ST,SL,SD,SM,LayoutScalar>::HostMirror >::value
                  ) >::type * = 0 )
  {
    (void)ViewAssignment( dst , "mirror" );
  }
};

template<>
struct ViewAssignment< LayoutScalar , LayoutScalar , void >
{
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutScalar> & src ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits ;
    typedef typename traits::memory_space  memory_space ;
    typedef typename traits::memory_traits memory_traits ;

    ViewTracking< traits >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< traits >::increment( dst.m_ptr_on_device );
  }

  /** \brief  Deep copy data from compatible value type, layout, rank, and specialization.  */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline static
  void deep_copy( const View<DT,DL,DD,DM,Impl::LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,Impl::LayoutScalar> & src ,
                  const typename Impl::enable_if<(
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::value_type ,
                                   typename ViewTraits<ST,SL,SD,SM>::non_const_value_type >::value
                  )>::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    if ( dst.m_ptr_on_device != src.m_ptr_on_device ) {

      DeepCopy< typename dst_traits::memory_space ,
                typename src_traits::memory_space >( dst.m_ptr_on_device , src.m_ptr_on_device ,
                                                     sizeof(typename dst_traits::value_type) );
    }
  }

  //------------------------------------
};

//----------------------------------------------------------------------------
// Assign scalar from default layout array

template<>
struct ViewAssignment< LayoutScalar , LayoutDefault , void >
{
  //------------------------------------
  // Scalar from Rank=1

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                  ) , unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape , 2 , i0 , i1 );

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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                  ) , unsigned >::type i2 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 3 , i0, i1, i2 );

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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                  ) , unsigned >::type i3 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 4 , i0, i1, i2, i3 );

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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                  ) , unsigned >::type i4 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 5 , i0, i1, i2, i3, i4 );

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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                  ) , unsigned >::type i5 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 6 , i0, i1, i2, i3, i4, i5 );

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
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                  ) , unsigned >::type i6 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 7 , i0, i1, i2, i3, i4, i5, i6 );

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
  KOKKOS_INLINE_FUNCTION
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
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                  ) , unsigned >::type i7 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, 8 , i0, i1, i2, i3, i4, i5, i6, i7 );

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

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

template< class DataType , class Arg1Type , class Arg2Type , class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::LayoutScalar >
  : public ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type >
{
private:

  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type > traits ;

  typedef Impl::ViewAssignment< Impl::LayoutScalar > alloc ;
  typedef Impl::ViewAssignment< Impl::LayoutScalar ,
                                Impl::LayoutScalar > assign ;

  typename traits::value_type * m_ptr_on_device ;

public:

  typedef View< typename traits::const_data_type,
                typename traits::array_layout,
                typename traits::device_type,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                Host ,
                void > HostMirror ;

  enum { Rank = 0 };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return typename traits::shape_type(); }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return 1 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & ) const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0) {}

  KOKKOS_INLINE_FUNCTION
  ~View()
    { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs )
    : m_ptr_on_device(0) { (void)assign( *this , rhs ); }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { (void)assign( *this , rhs ); return *this ; }

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,Impl::LayoutScalar> & rhs )
    : m_ptr_on_device(0) { (void)assign( *this , rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,Impl::LayoutScalar> & rhs )
    { (void)assign( *this , rhs ); return *this ; }

  //------------------------------------
  /** \brief  Allocation constructor */
  explicit
  View( const std::string & label ) : m_ptr_on_device(0) { (void)alloc( *this , label ); }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  const View & operator = ( const typename traits::value_type & rhs ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      *m_ptr_on_device = rhs ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const typename traits::value_type & rhs )
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      *m_ptr_on_device = rhs ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  operator typename traits::value_type & () const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type & operator()() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type & operator*() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const { return 1 ; }
};

//----------------------------------------------------------------------------

} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWSCALAR_HPP */

