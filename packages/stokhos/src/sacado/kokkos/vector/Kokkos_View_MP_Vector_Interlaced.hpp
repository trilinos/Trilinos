// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_MP_VECTOR_INTERLACED_HPP
#define KOKKOS_VIEW_MP_VECTOR_INTERLACED_HPP

#include "Sacado_MP_Vector.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Stokhos_ViewStorage.hpp"
#include <Kokkos_Core.hpp>

#include "Kokkos_View_Utils.hpp"
#include "Kokkos_View_MP_Vector_Utils.hpp"

/*
 * Specialization for Kokkos::View<Sacado::MP::Vector<Storage>...>
 * where the Sacado dimension is interlaced for LayoutLeft.
 *
 * Currently it can't be used at the same time as other such View
 * specializations due to conflicting specializations of AnalyzeShape.
 */

namespace Kokkos {
namespace Impl {

struct ViewMPVectorInterlaced {};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewMPVectorInterlaced
  , LayoutLeft
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewMPVectorInterlaced type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewMPVectorInterlaced
  , LayoutRight
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewMPVectorInterlaced type ;
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

/**\brief  View::value_type  == Sacado::MP::Vector< Stokhos::StorageType<...> > */
template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::ViewMPVectorInterlaced >
  : public ViewTraits< DataType
                     , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::array_layout
                     , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::device_type
                     , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::memory_traits
                     >
{
public:

  typedef ViewTraits< DataType
                    , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::array_layout
                    , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::device_type
                    , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::memory_traits
                    > traits ;

  // Type of const views with same value type
  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  // Type of non-const views with same value type
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > non_const_type ;

  // Host mirror
  typedef View< typename Impl::RebindStokhosStorageDevice<
                  typename traits::data_type ,
                  typename traits::host_mirror_space::memory_space >::type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                void > HostMirror ;

  // Equivalent array type for this view.
  typedef View< typename traits::array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > array_type ;

  // Equivalent const array type for this view.
  typedef View< typename traits::const_array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_array_type ;

  // Equivalent host array type for this view.
  typedef View< typename traits::array_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                typename traits::memory_traits > host_array_type ;

  // Equivalent const host array type for this view.
  typedef View< typename traits::const_array_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                typename traits::memory_traits > host_const_array_type ;

  typedef typename traits::value_type                   sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type  stokhos_storage_type ;
  typedef typename stokhos_storage_type::value_type     intrinsic_scalar_type ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::static_size };
  typedef integral_nonzero_constant< unsigned , StokhosStorageStaticDimension > sacado_size_type;

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;
  typedef typename array_type::traits::shape_type array_shape_type;

  typename stokhos_storage_type::value_type  * m_ptr_on_device ;
  typename traits::shape_type                  m_shape ;
  array_shape_type                             m_array_shape ; // Shape of intrinsic array
  stride_type                                  m_stride ;
  typename traits::execution_space::size_type      m_storage_size ; // Storage size of sacado dimension
  sacado_size_type                             m_sacado_size ; // Size of sacado dimension
  Impl::ViewDataManagement< traits >           m_management ;
  Impl::AllocationTracker                      m_tracker ;
  // Note:  if the view is partitioned, m_sacado_size != m_storage_size.
  // We always have m_storage_size >= m_sacado_size

  typedef Stokhos::ViewStorage<
    typename stokhos_storage_type::ordinal_type ,
    typename stokhos_storage_type::value_type ,
    StokhosStorageStaticDimension ,
      /* LayoutRight has stride-one stokhos storage */
    ( Impl::is_same< typename traits::array_layout , LayoutRight >::value ? 1 : 0 ) ,
    typename traits::device_type >  stokhos_view_storage_type ;

public:

  // This needs to be public so that we know what the return type of () is
  typedef Sacado::MP::Vector< stokhos_view_storage_type > reference_type ;

  // Whether the storage type is statically sized
  static const bool is_static = stokhos_storage_type::is_static;

  // Whether sacado dimension is contiguous
  static const bool is_contiguous =
    Impl::is_same< typename traits::array_layout , LayoutRight >::value;

  //------------------------------------
  // Shape for the Sacado::MP::Vector value_type ignores the internal static array length.
  enum { Rank = traits::rank };

  // Rank corresponding to the sacado dimension
  enum { Sacado_Rank = Rank+1 };

  KOKKOS_FORCEINLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type size() const
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
  KOKKOS_FORCEINLINE_FUNCTION
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
  void verify_dimension_storage_size( const typename traits::execution_space & dev ) const
  {
    const int length = dimension( Rank );

    const Impl::integral_nonzero_constant< int , StokhosStorageStaticDimension >
      per_thread( ! StokhosStorageStaticDimension ? length / dev.team_size() : 0 );

    if ( per_thread.value * dev.team_size() != length ) {
      Kokkos::abort("Kokkos::View< Sacado::MP::Vector ... > incompatible vector-size : team-size");
    }
  }
#else
  KOKKOS_INLINE_FUNCTION
  void verify_dimension_storage_size( const typename traits::execution_space & ) const {}
#endif

public:

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() { }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0), m_storage_size(0), m_sacado_size(0)
    {
      traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0);
      array_shape_type::assign(m_array_shape,0,0,0,0,0,0,0,0);
      stride_type::assign(m_stride,0);
    }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0), m_storage_size(0), m_sacado_size(0)
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
  // Allocation of a managed view with possible alignment padding.

  template< class AllocationProperties >
  explicit inline
  View( const AllocationProperties & prop ,
        // Impl::ViewAllocProp::size_type exists when the traits and allocation properties
        // are valid for allocating viewed memory.
        const typename Impl::ViewAllocProp< traits , AllocationProperties >::size_type n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    {
      typedef Impl::ViewAllocProp< traits , AllocationProperties > Alloc ;

      typedef typename traits::memory_space              memory_space ;
      typedef typename traits::shape_type                shape_type ;
      typedef typename stokhos_storage_type::value_type  scalar_type ;

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      array_shape_type::assign( m_array_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_array_shape );
      m_storage_size  = Impl::dimension( m_array_shape , unsigned(Rank) );
      m_sacado_size = m_storage_size;

      verify_dimension_storage_static_size();

      m_tracker = memory_space::allocate_and_track( Alloc::label( prop ) , sizeof(scalar_type) * Impl::capacity( m_array_shape , m_stride ) );

      m_ptr_on_device = (scalar_type *) m_tracker.alloc_ptr();

      (void) Kokkos::Impl::ViewDefaultConstruct< typename traits::execution_space , scalar_type , Alloc::Initialize >
          ( m_ptr_on_device , Impl::capacity( m_array_shape , m_stride ) );
    }

  //------------------------------------
  // Assign an unmanaged View from pointer, can be called in functors.
  // No alignment padding is performed.

  template< typename T >
  View( T * ptr ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        typename Impl::enable_if<(
          ( Impl::is_same<T,typename traits::value_type>::value ||
            Impl::is_same<T,typename traits::const_value_type>::value ) &&
          ! traits::is_managed ),
        const size_t >::type n7 = 0 )
    : m_ptr_on_device(ptr)
    {
      typedef typename traits::shape_type  shape_type ;

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      array_shape_type::assign( m_array_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_no_padding( m_stride , m_shape );
      m_storage_size  = Impl::dimension( m_array_shape , unsigned(Rank) );
      m_sacado_size = m_storage_size;
      m_management.set_unmanaged();

      verify_dimension_storage_static_size();
    }

  //------------------------------------
  // Is not allocated

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // Scalar operator on traits::rank == 1

  typedef std::conditional< ( traits::rank == 1 ),
                      reference_type ,
                      Impl::ViewError::scalar_operator_called_from_non_scalar_view >
    if_scalar_operator ;

  KOKKOS_FORCEINLINE_FUNCTION
  typename if_scalar_operator::type
    operator()() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device ,
        m_shape.N0 , 1 ) );
    }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 2:

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type , traits, LayoutLeft, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + i0 ,
        m_shape.N1 ,
        m_stride.value ) );
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_stride.value * i0 ) ,
        m_shape.N1 , 1 ) );
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type , traits, typename traits::array_layout, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 2,
                                      iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * ( i1 )),
        m_shape.N2 ,
        m_stride.value * m_shape.N1 ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0);
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i1 ) + m_stride.value * i0 ) ,
        m_shape.N2 , 1 ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 3,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 ))),
        m_shape.N3 ,
        m_stride.value * m_shape.N1 * m_shape.N2 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i2 +
                            m_shape.N2 * ( i1 )) +
                            m_stride.value * i0 ) ,
        m_shape.N3 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 4,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 )))),
        m_shape.N4 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))) +
                            m_stride.value * i0 ) ,
        m_shape.N4 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 5,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 ))))),
        m_shape.N5 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( m_storage_size * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))) +
                            m_stride.value * i0 ) ,
        m_shape.N5 , 1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 6,
                                      iType0, iType1, iType2, iType3, iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , int , int , int ) const
    { return operator()(i0,i1,i2,i3,i4); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
                 const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 )))))),
        m_shape.N6 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 7,
                                      iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , int , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 8:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );

      // Strided storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
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

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return reference_type( stokhos_view_storage_type(
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type ,
                                      traits, typename traits::array_layout, 8,
                                      iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5,i6); }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::value_type::storage_type::value_type *
    data() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_array_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::capacity( m_array_shape , m_stride ); }

  // Static storage size
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type sacado_size() const
  { return m_sacado_size.value; }
};

/** \brief  A deep copy between views of the same specialization, compatible type,
 *          same rank, same layout are handled by that specialization.
 */
template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
inline
void deep_copy( const View<DT,DL,DD,DM,Impl::ViewMPVectorInterlaced> & dst ,
                const View<ST,SL,SD,SM,Impl::ViewMPVectorInterlaced> & src ,
                typename Impl::enable_if<(
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewMPVectorInterlaced>::intinsic_scalar_type ,
                                 typename View<ST,SL,SD,SM,Impl::ViewMPVectorInterlaced>::intinsic_scalar_type >::value
                  &&
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewMPVectorInterlaced>::array_layout ,
                                 typename View<ST,SL,SD,SM,Impl::ViewMPVectorInterlaced>::array_layout >::value
                  &&
                  ( unsigned(View<DT,DL,DD,DM,Impl::ViewMPVectorInterlaced>::rank) ==
                    unsigned(View<ST,SL,SD,SM,Impl::ViewMPVectorInterlaced>::rank) )
                )>::type * = 0 )
{
  typedef  View<DT,DL,DD,DM,Impl::ViewMPVectorInterlaced>  dst_type ;
  typedef  View<ST,SL,SD,SM,Impl::ViewMPVectorInterlaced>  src_type ;

  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( dst.data() != src.data() ) {

    Impl::assert_shapes_are_equal( dst.shape() , src.shape() );

    const size_t nbytes = sizeof(typename dst_type::value_type::storage_type::value_type) * dst.span();

    Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
  }
}

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
 */
template< class StorageType >
struct AnalyzeShape< Sacado::MP::Vector< StorageType > >
  : Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > // Treat as a scalar
{
private:

  typedef AnalyzeShape< typename StorageType::value_type > nested ;

public:

  typedef typename ViewMPVectorInterlaced specialize ;

  typedef Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > shape ;

  // If ( ! StorageType::is_static ) then 0 == StorageType::static_size and the first array declaration is not used.
  // However, the compiler will still generate this type declaration and it must not have a zero length.
  typedef typename
    std::conditional< StorageType::is_static
        , typename nested::array_intrinsic_type [ StorageType::is_static ? StorageType::static_size : 1 ]
        , typename nested::array_intrinsic_type *
        >::type array_intrinsic_type ;

  typedef typename
    std::conditional< StorageType::is_static
        , typename nested::const_array_intrinsic_type [ StorageType::is_static ? StorageType::static_size : 1 ]
        , typename nested::const_array_intrinsic_type *
        >::type const_array_intrinsic_type ;

  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::MP::Vector< StorageType >  type ;
  typedef const Sacado::MP::Vector< StorageType >  const_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_type ;

  typedef       Sacado::MP::Vector< StorageType >  value_type ;
  typedef const Sacado::MP::Vector< StorageType >  const_value_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_value_type ;
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewMPVectorInterlaced , ViewMPVectorInterlaced , void >
{
  typedef ViewMPVectorContiguous specialize ;

  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst
                , const View<ST,SL,SD,SM,specialize> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM>                   dst_traits ;
    typedef View<DT,DL,DD,DM,specialize>              dst_type ;
    typedef typename dst_type::shape_type             shape_type ;
    typedef typename dst_type::array_shape_type       array_shape_type ;
    typedef typename dst_type::stride_type            stride_type ;

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );
    array_shape_type::assign( dst.m_array_shape,
                              src.m_array_shape.N0 , src.m_array_shape.N1 , src.m_array_shape.N2 , src.m_array_shape.N3 ,
                              src.m_array_shape.N4 , src.m_array_shape.N5 , src.m_array_shape.N6 , src.m_array_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );
    dst.m_ptr_on_device = src.m_ptr_on_device ;
    dst.m_storage_size  = src.m_storage_size ;
    dst.m_sacado_size  = src.m_sacado_size ;
    dst.m_tracker = src.m_tracker ;
  }

  //------------------------------------
  /** \brief  Partition of compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst
                , const View<ST,SL,SD,SM,specialize> & src
                , typename enable_if<(
                    // Same intrinsic scalar type
                    is_same< typename View<DT,DL,DD,DM,specialize>::intrinsic_scalar_type ,
                             typename View<ST,SL,SD,SM,specialize>::intrinsic_scalar_type >::value
                    &&
                    // Same memory space
                    is_same< typename View<DT,DL,DD,DM,specialize>::memory_space ,
                             typename View<ST,SL,SD,SM,specialize>::memory_space >::value
                    &&
                    // Same layout
                    is_same< typename View<DT,DL,DD,DM,specialize>::array_layout ,
                             typename View<ST,SL,SD,SM,specialize>::array_layout >::value
                    &&
                    // Same rank
                    ( unsigned(View<DT,DL,DD,DM,specialize>::rank) ==
                      unsigned(View<ST,SL,SD,SM,specialize>::rank) )
                    &&
                    // Destination is not managed
                    ! View<DT,DL,DD,DM,specialize>::is_managed
                  ), const Sacado::MP::VectorPartition & >::type part )
  {
    typedef ViewTraits<DT,DL,DD,DM>                           dst_traits ;
    typedef View<DT,DL,DD,DM,specialize>                      dst_type ;
    typedef typename dst_type::shape_type                     dst_shape_type ;
    typedef typename dst_type::array_shape_type               dst_array_shape_type ;
    typedef typename dst_type::stride_type                    dst_stride_type ;
    typedef typename dst_traits::value_type                   dst_sacado_mp_vector_type ;
    typedef typename dst_sacado_mp_vector_type::storage_type  dst_stokhos_storage_type ;

    enum { DstRank         = dst_type::rank };
    enum { DstStaticLength = dst_stokhos_storage_type::static_size };

    const int length = part.end - part.begin ;

    if ( DstStaticLength && DstStaticLength != length ) {
      Kokkos::abort("Kokkos::View< Sacado::MP::Vector ... > incompatible partitioning");
    }

    unsigned dims[8];
    dims[0] = src.m_array_shape.N0;
    dims[1] = src.m_array_shape.N1;
    dims[2] = src.m_array_shape.N2;
    dims[3] = src.m_array_shape.N3;
    dims[4] = src.m_array_shape.N4;
    dims[5] = src.m_array_shape.N5;
    dims[6] = src.m_array_shape.N6;
    dims[7] = src.m_array_shape.N7;
    unsigned rank = dst_type::rank;

    dst_shape_type::assign( dst.m_shape,
                            dims[0] , dims[1] , dims[2] , dims[3] ,
                            dims[4] , dims[5] , dims[6] , dims[7] );

    dims[rank] = length;
    dst_array_shape_type::assign( dst.m_array_shape,
                                  dims[0] , dims[1] , dims[2] , dims[3] ,
                                  dims[4] , dims[5] , dims[6] , dims[7] );

    dst_stride_type::assign( dst.m_stride , src.m_stride.value );

    // Original Sacado::MP::Vector length
    dst.m_storage_size = src.m_storage_size ;
    dst.m_sacado_size = length;

    if ( Impl::is_same< typename dst_traits::array_layout , LayoutLeft >::value ) {
      dst.m_ptr_on_device = src.m_ptr_on_device + part.begin *
                      ( 0 == DstRank ? 1 : dst.m_stride.value * (
                      ( 1 == DstRank ? 1 : dst.m_shape.N1 * (
                      ( 2 == DstRank ? 1 : dst.m_shape.N2 * (
                      ( 3 == DstRank ? 1 : dst.m_shape.N3 * (
                      ( 4 == DstRank ? 1 : dst.m_shape.N4 * (
                      ( 5 == DstRank ? 1 : dst.m_shape.N5 * (
                      ( 6 == DstRank ? 1 : dst.m_shape.N6 )))))))))))));
    }
    else { // if ( Impl::is_same< typename traits::array_layout , LayoutRight >::value )
      dst.m_ptr_on_device = src.m_ptr_on_device + part.begin ;
    }
    dst.m_tracker = src.m_tracker ;
  }
};

template<>
struct ViewAssignment< ViewDefault , ViewMPVectorInterlaced , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,ViewMPVectorInterlaced>::array_type & dst
                , const    View<ST,SL,SD,SM,ViewMPVectorInterlaced> & src )
  {
    typedef View<ST,SL,SD,SM,ViewMPVectorInterlaced>          src_type ;
    typedef typename src_type::value_type                     src_sacado_mp_vector_type ;
    typedef typename src_sacado_mp_vector_type::storage_type  src_stokhos_storage_type ;

    typedef typename src_type::array_type   dst_type ;
    typedef typename dst_type::shape_type   dst_shape_type ;
    typedef typename dst_type::stride_type  dst_stride_type ;

    dst_shape_type::assign( dst.m_shape,
                            src.m_array_shape.N0 , src.m_array_shape.N1 , src.m_array_shape.N2 , src.m_arrat_shape.N3 ,
                            src.m_array_shape.N4 , src.m_array_shape.N5 , src.m_arrat_shape.N6 , src.m_arrat_shape.N7 );

    dst_stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = reinterpret_cast< typename dst_type::value_type *>( src.m_ptr_on_device );

    dst.m_tracker = src.m_tracker ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEW_MP_VECTOR_HPP */
