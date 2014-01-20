// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_MP_VECTOR_HPP
#define KOKKOS_VIEW_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Stokhos_ViewStorage.hpp"
#include <Kokkos_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct ViewSpecializeSacadoMPVector {};

template< class StorageType ,
          class Rank , class RankDynamic ,
          class MemoryTraits >
struct ViewSpecialize<
  typename StorageType::value_type ,
  Sacado::MP::Vector< StorageType > , // View::value_type
  LayoutLeft ,
  Rank , RankDynamic ,
  typename StorageType::device_type::memory_space ,
  MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class StorageType ,
          class Rank , class RankDynamic ,
          class MemoryTraits >
struct ViewSpecialize<
  const typename StorageType::value_type ,
  const Sacado::MP::Vector< StorageType > , // View::value_type
  LayoutLeft ,
  Rank , RankDynamic ,
  typename StorageType::device_type::memory_space ,
  MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class StorageType ,
          class Rank , class RankDynamic ,
          class MemoryTraits >
struct ViewSpecialize<
  typename StorageType::value_type ,
  Sacado::MP::Vector< StorageType > , // View::value_type
  LayoutRight ,
  Rank , RankDynamic ,
  typename StorageType::device_type::memory_space ,
  MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class StorageType ,
          class Rank , class RankDynamic ,
          class MemoryTraits >
struct ViewSpecialize<
  const typename StorageType::value_type ,
  const Sacado::MP::Vector< StorageType > , // View::value_type
  LayoutRight ,
  Rank , RankDynamic ,
  typename StorageType::device_type::memory_space ,
  MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class T , class Device > struct RebindStokhosStorageDevice ;

//----------------------------------------------------------------------------
/** \brief  Enable view parentheses operator for
 *          match of layout and integral arguments.
 *          If correct rank define type from traits,
 *          otherwise define type as an error message.
 *
 *          We have to create our own here because the one in Kokkos_View.hpp
 *          was recently changed in an incompatible way.
 */
template< class ReturnType , class Traits , class Layout , unsigned Rank ,
          typename iType0 = int , typename iType1 = int ,
          typename iType2 = int , typename iType3 = int ,
          typename iType4 = int , typename iType5 = int ,
          typename iType6 = int , typename iType7 = int ,
          class Enable = void >
struct ViewEnableArrayOperLayout ;

template< class ReturnType , class Traits , class Layout , unsigned Rank ,
          typename iType0 , typename iType1 ,
          typename iType2 , typename iType3 ,
          typename iType4 , typename iType5 ,
          typename iType6 , typename iType7 >
struct ViewEnableArrayOperLayout<
   ReturnType , Traits , Layout , Rank ,
   iType0 , iType1 , iType2 , iType3 ,
   iType4 , iType5 , iType6 , iType7 ,
   typename enable_if<
     iType0(0) == 0 && iType1(0) == 0 && iType2(0) == 0 && iType3(0) == 0 &&
     iType4(0) == 0 && iType5(0) == 0 && iType6(0) == 0 && iType7(0) == 0 &&
     is_same< typename Traits::array_layout , Layout >::value &&
     ( unsigned(Traits::rank) == Rank )
   >::type >
{
  typedef ReturnType type ;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace ViewError {

struct sacado_mp_vector_partition_constructor_requires_unmanaged_view {};

} // namespace ViewError
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/**\brief  View::value_type  == Sacado::MP::Vector< Stokhos::StorageType<...> >
 *         View::scalar_type == StorageType<...>::value_type
 */
template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::ViewSpecializeSacadoMPVector >
  : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:

  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef typename traits::value_type                   sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type  stokhos_storage_type ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::is_static ? stokhos_storage_type::static_size : 0 };

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  typename traits::scalar_type  * m_ptr_on_device ;
  typename traits::shape_type     m_shape ;
  stride_type                     m_stride ;
  typename traits::device_type::size_type m_storage_size ;

  typedef Stokhos::ViewStorage<
    typename stokhos_storage_type::ordinal_type ,
    typename traits::scalar_type ,
    StokhosStorageStaticDimension ,
      /* LayoutRight has stride-one stokhos storage */
    ( Impl::is_same< typename traits::array_layout , LayoutRight >::value ? 1 : 0 ) ,
    typename traits::device_type >  stokhos_view_storage_type ;

public:

  // This needs to be public so that we know what the return type of () is
  typedef Sacado::MP::Vector< stokhos_view_storage_type >  sacado_mp_vector_view_type ;

  // Whether the storage type is statically sized
  static const bool is_static = stokhos_storage_type::is_static;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > non_const_type ;

  typedef View< typename traits::array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > array_type ;

  typedef View< typename Impl::RebindStokhosStorageDevice<
                  typename traits::data_type ,
                  typename traits::device_type::host_mirror_device_type >::type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  // Rank for multidimensional array of the Sacado::MP::Vector value_type
  // is one less than the rank of the array of intrinsic scalar_type defined by the shape.
  enum { Rank = traits::rank - 1 };

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

private:

  // Restrict allocation to 'StokhosStorageStaticDimension'
  inline
  void verify_dimension_storage_static_size() const
  {
    if ( dimension( unsigned(Rank) ) % ( StokhosStorageStaticDimension ? StokhosStorageStaticDimension : 1 ) ) {
      std::ostringstream msg ;
      msg << "Kokkos::View< Sacado::MP::Vector<StorageType , ... > allocation dimension ("
          << dimension( unsigned(Rank) )
          << ") must be a multiple of StorageType::static_size ("
          << StokhosStorageStaticDimension
          << ")" ;
      Impl::throw_runtime_exception( msg.str() );
    }
  }

#if defined( KOKKOS_EXPRESSION_CHECK )
  KOKKOS_INLINE_FUNCTION
  void verify_dimension_storage_size( const typename traits::device_type & dev ) const
  {
    const int length = dimension( Rank );

    const Impl::integral_nonzero_constant< int , StokhosStorageStaticDimension >
      per_thread( ! StokhosStorageStaticDimension ? length / dev.team_size() : 0 );

    if ( per_thread.value * dev.team_size() != length ) {
      const char msg[] = "Kokkos::View< Sacado::MP::Vector ... > incompatible vector-size : team-size" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort(msg);
#else
      throw std::runtime_error(msg);
#endif
    }
  }
#else
  KOKKOS_INLINE_FUNCTION
  void verify_dimension_storage_size( const typename traits::device_type & ) const {}
#endif

public:

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
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
    }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Construct or assign compatible view:

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    : m_ptr_on_device(0)
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
    }

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
      return *this ;
    }

  //------------------------------------

  struct Partition {
    typename traits::size_type rank ;
    typename traits::size_type size ;
    KOKKOS_INLINE_FUNCTION
    Partition( typename traits::size_type r ,
               typename traits::size_type s )
      : rank( r ), size( s ) {}
  };

  template< class ViewRHS >
  KOKKOS_INLINE_FUNCTION
  View( const ViewRHS & rhs ,
        typename Impl::enable_if< is_view< ViewRHS >::value &&
                                  Impl::is_same< typename ViewRHS::traits::specialize ,
                                                 typename traits::specialize >::value &&
                                  Impl::is_same< typename ViewRHS::traits::array_type ,
                                                 typename traits::array_type >::value &&
                                  Impl::is_same< typename ViewRHS::traits::device_type ,
                                                 typename traits::device_type >::value &&
                                  ( ! traits::is_managed ) , Partition >::type & part )
    : m_ptr_on_device(0)
    {
      typedef typename traits::shape_type   shape_type ;
      typedef typename traits::scalar_type  scalar_type ;

      if ( rhs.m_storage_size % part.size ) {
        const char msg[] = "Kokkos::View< Sacado::MP::Vector ... > unbalanced partitioning" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
        cuda_abort(msg);
#else
        throw std::runtime_error(msg);
#endif
      }

      shape_type::assign( m_shape ,
                          ( Rank == 0 ? rhs.m_shape.N0 / part.size : rhs.m_shape.N0 ) ,
                          ( Rank == 1 ? rhs.m_shape.N1 / part.size : rhs.m_shape.N1 ) ,
                          ( Rank == 2 ? rhs.m_shape.N2 / part.size : rhs.m_shape.N2 ) ,
                          ( Rank == 3 ? rhs.m_shape.N3 / part.size : rhs.m_shape.N3 ) ,
                          ( Rank == 4 ? rhs.m_shape.N4 / part.size : rhs.m_shape.N4 ) ,
                          ( Rank == 5 ? rhs.m_shape.N5 / part.size : rhs.m_shape.N5 ) ,
                          ( Rank == 6 ? rhs.m_shape.N6 / part.size : rhs.m_shape.N6 ) ,
                          ( Rank == 7 ? rhs.m_shape.N7 / part.size : rhs.m_shape.N7 ) );

      stride_type::assign( m_stride , rhs.m_stride.value );

      // Original Sacado::MP::Vector length
      m_storage_size = rhs.m_storage_size ;

      if ( Impl::is_same< typename traits::array_layout , LayoutLeft >::value ) {
        m_ptr_on_device = rhs.m_ptr_on_device + part.rank *
                        ( 0 == Rank ? m_shape.N0 : m_stride.value * m_shape.N1 * (
                        ( 1 == Rank ? 1 : m_shape.N2 * (
                        ( 2 == Rank ? 1 : m_shape.N3 * (
                        ( 3 == Rank ? 1 : m_shape.N4 * (
                        ( 4 == Rank ? 1 : m_shape.N5 * (
                        ( 5 == Rank ? 1 : m_shape.N6 * (
                        ( 6 == Rank ? 1 : m_shape.N7 )))))))))))));
      }
      else { // if ( Impl::is_same< typename traits::array_layout , LayoutRight >::value )
        m_ptr_on_device = rhs.m_ptr_on_device + part.rank * Impl::dimension( m_shape , unsigned(Rank) );
      }
    }

  //------------------------------------
  // Allocation of a managed view with possible alignment padding.

  typedef Impl::if_c< traits::is_managed ,
                      std::string ,
                      Impl::ViewError::allocation_constructor_requires_managed >
   if_allocation_constructor ;

  explicit inline
  View( const typename if_allocation_constructor::type & label ,
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
      typedef typename traits::memory_space  memory_space ;
      typedef typename traits::shape_type    shape_type ;
      typedef typename traits::scalar_type   scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_shape );

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_shape , unsigned(Rank) );
      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );

      (void) Impl::ViewFill< array_type >( *this , typename traits::scalar_type() );
    }

  explicit inline
  View( const AllocateWithoutInitializing & ,
        const typename if_allocation_constructor::type & label ,
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
      typedef typename traits::memory_space  memory_space ;
      typedef typename traits::shape_type    shape_type ;
      typedef typename traits::scalar_type   scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_shape );

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_shape , unsigned(Rank) );
      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );
    }

  //------------------------------------
  // Assign an unmanaged View from pointer, can be called in functors.
  // No alignment padding is performed.

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::scalar_type * ,
                      Impl::ViewError::user_pointer_constructor_requires_unmanaged >
    if_user_pointer_constructor ;

  View( typename if_user_pointer_constructor::type ptr ,
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
      typedef typename traits::shape_type   shape_type ;
      typedef typename traits::scalar_type  scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_no_padding( m_stride , m_shape );

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_shape , unsigned(Rank) );
      m_ptr_on_device = if_user_pointer_constructor::select( ptr );
    }


  //------------------------------------
  // Is not allocated

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // Scalar operator on traits::rank == 1

  typedef Impl::if_c< ( traits::rank == 1 ),
                      sacado_mp_vector_view_type ,
                      Impl::ViewError::scalar_operator_called_from_non_scalar_view >
    if_scalar_operator ;

  KOKKOS_INLINE_FUNCTION
  typename if_scalar_operator::type
    operator()() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device ,
        m_shape.N0 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 2:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type , traits, LayoutLeft, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + i0 ,
        m_shape.N1 ,
        m_stride.value ) );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type , traits, LayoutLeft, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + i0 ,
        m_shape.N1 ,
        m_stride.value ) );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 2, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + i0 ,
        m_shape.N1 ,
        m_stride.value ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * ( i1 )),
        m_shape.N2 ,
        m_stride.value * m_shape.N1 ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 3, iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * ( i1 )),
        m_shape.N2 ,
        m_stride.value * m_shape.N1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 ))),
        m_shape.N3 ,
        m_stride.value * m_shape.N1 * m_shape.N2 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 ))),
        m_shape.N3 ,
        m_stride.value * m_shape.N1 * m_shape.N2 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 )))),
        m_shape.N4 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 )))),
        m_shape.N4 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 ))))),
        m_shape.N5 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 ))))),
        m_shape.N5 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
                 const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 )))))),
        m_shape.N6 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
        const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ,
        const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 )))))),
        m_shape.N6 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, traits::rank 8:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 + m_shape.N5 * (
                            i6 ))))))),
        m_shape.N7 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 * m_shape.N6 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const int ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 + m_shape.N5 * (
                            i6 ))))))),
        m_shape.N7 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 * m_shape.N6 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 2:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_stride.value * i0 ) ,
        m_shape.N1 , 1 ) );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 2, iType0 >::type
    at( const iType0 & i0 , const int, const int, const int,
        const int, const int, const int, const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_stride.value * i0 ) ,
        m_shape.N1 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0);
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i1 ) + m_stride.value * i0 ) ,
        m_shape.N2 , 1 ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 3, iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0);
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i1 ) + m_stride.value * i0 ) ,
        m_shape.N2 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i2 +
                            m_shape.N2 * ( i1 )) +
                            m_stride.value * i0 ) ,
        m_shape.N3 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i2 +
                            m_shape.N2 * ( i1 )) +
                            m_stride.value * i0 ) ,
        m_shape.N3 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))) +
                            m_stride.value * i0 ) ,
        m_shape.N4 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))) +
                            m_stride.value * i0 ) ,
        m_shape.N4 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))) +
                            m_stride.value * i0 ) ,
        m_shape.N5 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))) +
                            m_stride.value * i0 ) ,
        m_shape.N5 , 1 ) );
    }


  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 7:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))))) +
                            m_stride.value * i0 ) ,
        m_shape.N6 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))))) +
                            m_stride.value * i0 ) ,
        m_shape.N6 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, traits::rank 8:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i6 +
                            m_shape.N6 * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))))) +
                            m_stride.value * i0 ) ,
        m_shape.N7 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOperLayout< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i6 +
                            m_shape.N6 * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))))) +
                            m_stride.value * i0 ) ,
        m_shape.N7 , 1 ) );
    }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_INLINE_FUNCTION
  typename traits::scalar_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::capacity( m_shape , m_stride ); }

  // Static storage size
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type static_storage_size() const
  { return StokhosStorageStaticDimension; }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Analyze the array shape of a Sacado::MP::Vector.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::MP::Vector< StorageType > , ... >
 *  can be determined at compile-time.
 *
 *  The dimension associated with the MP::Vector is always dynamic.
 */
template< class StorageType >
struct AnalyzeShape< Sacado::MP::Vector< StorageType > >
  : public ShapeInsert< typename AnalyzeShape< typename StorageType::value_type >::shape , 0 >::type
{
private:
  typedef AnalyzeShape< typename StorageType::value_type > nested ;
public:

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::scalar_type            scalar_type ;
  typedef typename nested::const_scalar_type      const_scalar_type ;
  typedef typename nested::non_const_scalar_type  non_const_scalar_type ;

  typedef typename nested::array_type           * array_type ;
  typedef typename nested::const_array_type     * const_array_type ;
  typedef typename nested::non_const_array_type * non_const_array_type ;

  typedef       Sacado::MP::Vector< StorageType >  type ;
  typedef const Sacado::MP::Vector< StorageType >  const_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_type ;

  typedef       Sacado::MP::Vector< StorageType >  value_type ;
  typedef const Sacado::MP::Vector< StorageType >  const_value_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_value_type ;
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewSpecializeSacadoMPVector , ViewSpecializeSacadoMPVector , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector> & dst
                , const View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef typename View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector>::shape_type   shape_type ;
    typedef typename View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector>::stride_type  stride_type ;

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_storage_size  = src.m_storage_size ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }
};

template<>
struct ViewAssignment< LayoutDefault , ViewSpecializeSacadoMPVector , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutDefault> & dst
                , const View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef typename View<DT,DL,DD,DM,LayoutDefault>::shape_type   shape_type ;
    typedef typename View<DT,DL,DD,DM,LayoutDefault>::stride_type  stride_type ;

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

template< class T , class Device >
struct RebindStokhosStorageDevice< T * , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , class Device >
struct RebindStokhosStorageDevice< T [] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , unsigned N , class Device >
struct RebindStokhosStorageDevice< T[N] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type type[N] ;
};

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< Sacado::MP::Vector< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::MP::Vector< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef typename NewVectorApply::type type ;
};

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< const Sacado::MP::Vector< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::MP::Vector< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef const typename NewVectorApply::type type ;
};

//----------------------------------------------------------------------------

} // namespace Impl

// Type name for a local, unmanaged view with possibly a different static size
template <typename ViewType,
          unsigned LocalSize,
          unsigned Rank = ViewType::Rank,
          bool isStatic = ViewType::is_static>
struct LocalMPVectorView {};

template <typename ViewType,
          unsigned LocalSize>
struct LocalMPVectorView< ViewType, LocalSize, 1, true > {
  typedef typename ViewType::value_type vector_type;
  typedef typename ViewType::array_layout array_layout;
  typedef typename ViewType::device_type device_type;
  typedef typename vector_type::storage_type storage_type;
  typedef typename storage_type::template apply_N<LocalSize> StorageApply;
  typedef typename StorageApply::type local_storage_type;
  typedef Sacado::MP::Vector< local_storage_type > local_value_type;

  typedef Kokkos::View< local_value_type*,
                        array_layout,
                        device_type,
                        Kokkos::MemoryUnmanaged > type;
};

template <typename ViewType,
          unsigned LocalSize>
struct LocalMPVectorView<ViewType, LocalSize, 1, false> {
  typedef typename ViewType::value_type vector_type;
  typedef typename ViewType::array_layout array_layout;
  typedef typename ViewType::device_type device_type;

  typedef Kokkos::View< vector_type*,
                        array_layout,
                        device_type,
                        Kokkos::MemoryUnmanaged > type;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEW_MP_VECTOR_HPP */
