// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_FAD_HPP
#define KOKKOS_VIEW_FAD_HPP

// Make sure the user really wants these View specializations
#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Sacado_Traits.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_AnalyzeSacadoShape.hpp"
#include "impl/Kokkos_Error.hpp"
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
#include "Cuda/Kokkos_Cuda_abort.hpp"
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Sacado {
namespace Fad {

/**\brief  Define a partition of a View of Sacado::MP::Vector type */
struct VectorPartition {
  unsigned begin ;
  unsigned end ;

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  VectorPartition( const iType0 & i0 , const iType1 & i1 ) : begin(i0), end(i1) {}
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct ViewSpecializeSacadoFad {};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewSpecializeSacadoFad
  , LayoutLeft
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoFad type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewSpecializeSacadoFad
  , LayoutRight
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoFad type ;
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace ViewError {

struct sacado_fad_partition_constructor_requires_unmanaged_view {};

} // namespace ViewError
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/*
 * \brief Specialization of Kokkos::View for any Sacado Fad value type.
 *
 * Usage is to append another dimension to the view of size fad_size+1 where
 * fad_size is the size of the desired derivative array.
*/
template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type ,
            Impl::ViewSpecializeSacadoFad >
  : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:

  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

  typedef typename traits::value_type fad_type ;
  typedef typename Sacado::ValueType<fad_type>::type fad_value_type ;
  typedef typename Kokkos::Impl::add_const<fad_value_type>::type const_fad_value_type ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  enum { FadStaticDimension = Sacado::StaticSize<fad_type>::value };

  /* LayoutRight has stride-one storage */
  enum { FadStaticStride = ( Impl::is_same< typename traits::array_layout , LayoutRight >::value ? 1 : 0 ) };

  typedef Impl::AnalyzeSacadoShape< typename traits::data_type,
                                    typename traits::array_layout > analyze_sacado_shape;
  typedef Impl::ViewOffset< typename analyze_sacado_shape::shape ,
                            typename traits::array_layout > offset_map_type ;

  fad_value_type                             * m_ptr_on_device ;
  offset_map_type                              m_offset_map ;
  typename traits::device_type::size_type      m_storage_size ;
  Impl::ViewDataManagement< traits >           m_management ;

public:

  // This needs to be public so that we know what the return type of () is
  typedef typename Sacado::ViewFadType<fad_type, FadStaticDimension, FadStaticStride>::type fad_view_type ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > non_const_type ;

  typedef View< typename analyze_sacado_shape::array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > array_type ;

  typedef View< typename traits::data_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  // Rank for multidimensional array of the Fad value_type
  // is one less than the rank of the array of intrinsic fad_value_type defined by the shape.
  enum { Rank = traits::rank };

  KOKKOS_FORCEINLINE_FUNCTION typename traits::shape_type shape() const {
    typedef typename traits::shape_type shape_type;
    shape_type s;
    shape_type::assign(
      s, m_offset_map.N0, m_offset_map.N1, m_offset_map.N2, m_offset_map.N3,
         m_offset_map.N4, m_offset_map.N5, m_offset_map.N6, m_offset_map.N7 );
    return s ;
  }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_0() const {
    return unsigned(Rank) >= 1 ? m_offset_map.N0 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_1() const {
    return unsigned(Rank) >= 2 ? m_offset_map.N1 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_2() const {
    return unsigned(Rank) >= 3 ? m_offset_map.N2 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_3() const {
    return unsigned(Rank) >= 4 ? m_offset_map.N3 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_4() const {
    return unsigned(Rank) >= 5 ? m_offset_map.N4 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_5() const {
    return unsigned(Rank) >= 6 ? m_offset_map.N5 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_6() const {
    return unsigned(Rank) >= 7 ? m_offset_map.N6 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_7() const {
    return unsigned(Rank) >= 8 ? m_offset_map.N7 : 1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type size() const {
    return   dimension_0()
           * dimension_1()
           * dimension_2()
           * dimension_3()
           * dimension_4()
           * dimension_5()
           * dimension_6()
           * dimension_7()
           ;
  }

  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return i < iType(Rank) ? Impl::dimension( m_offset_map , i ) : 1 ; }

  //------------------------------------

private:

  // Restrict allocation to 'FadStaticDimension'
  inline
  void verify_dimension_storage_static_size() const
  {
    if ( Impl::dimension( m_offset_map , unsigned(Rank) ) % ( FadStaticDimension ? FadStaticDimension+1 : 1 ) ) {
      std::ostringstream msg ;
      msg << "Kokkos::View< FadType , ... > allocation dimension ("
          << Impl::dimension( m_offset_map , unsigned(Rank) )
          << ") must be a multiple of StorageType::static_size ("
          << FadStaticDimension+1
          << ")" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort( msg.str().c_str() );
#else
      Impl::throw_runtime_exception( msg.str() );
#endif
    }
  }

public:

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() { m_management.decrement( m_ptr_on_device ); }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0)
    { m_offset_map.assign(0,0,0,0,0,0,0,0); }

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
 
      typedef typename traits::memory_space  memory_space ;

      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
      m_offset_map.set_padding();

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_offset_map , unsigned(Rank) );

      m_ptr_on_device = (fad_value_type *)
        memory_space::allocate( Alloc::label( prop ),
                                typeid(fad_value_type) ,
                                sizeof(fad_value_type) ,
                                m_offset_map.capacity() );

      if ( Alloc::initialize() ) {
        (void) Kokkos::Impl::DefaultConstruct
          < typename traits::execution_space , fad_value_type >
            ( m_ptr_on_device , m_offset_map.capacity() );
      }
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
            Impl::is_same<T,fad_value_type>::value ||
            Impl::is_same<T,const_fad_value_type>::value
          ),
        const size_t >::type n7 = 0 )
    : m_ptr_on_device(ptr)
    {
      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );

      verify_dimension_storage_static_size();

      m_storage_size = Impl::dimension( m_offset_map , unsigned(Rank) );

      m_management.set_unmanaged();
    }

  //------------------------------------
  // Assign unmanaged View to portion of Device shared memory

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::device_type ,
                      Impl::ViewError::device_shmem_constructor_requires_unmanaged >
      if_device_shmem_constructor ;

  explicit KOKKOS_INLINE_FUNCTION
  View( typename if_device_shmem_constructor::type & dev ,
        const unsigned n0 = 0 ,
        const unsigned n1 = 0 ,
        const unsigned n2 = 0 ,
        const unsigned n3 = 0 ,
        const unsigned n4 = 0 ,
        const unsigned n5 = 0 ,
        const unsigned n6 = 0 ,
        const unsigned n7 = 0 )
    : m_ptr_on_device(0)
    {
      enum { align = 8 };
      enum { mask  = align - 1 };

      m_offset_map.assign(  n0, n1, n2, n3, n4, n5, n6, n7 );

      typedef Impl::if_c< ! traits::is_managed ,
                          fad_value_type * ,
                          Impl::ViewError::device_shmem_constructor_requires_unmanaged >
        if_device_shmem_pointer ;

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_offset_map , unsigned(Rank) );

      // Select the first argument:
      m_ptr_on_device = if_device_shmem_pointer::select(
        (fad_value_type *) dev.get_shmem( shmem_size(n0,n1,n2,n3,n4,n5,n6,n7) ) );
    }

  static KOKKOS_INLINE_FUNCTION
  unsigned shmem_size( const unsigned n0 = 0 ,
                       const unsigned n1 = 0 ,
                       const unsigned n2 = 0 ,
                       const unsigned n3 = 0 ,
                       const unsigned n4 = 0 ,
                       const unsigned n5 = 0 ,
                       const unsigned n6 = 0 ,
                       const unsigned n7 = 0 )
  {
    enum { align = 8 };
    enum { mask  = align - 1 };

    offset_map_type offset_map ;

    offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );

    return unsigned( sizeof(fad_value_type) * offset_map.capacity() + unsigned(mask) ) & ~unsigned(mask) ;
  }

  //------------------------------------
  // Is not allocated

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // Scalar operator on traits::rank == 0

  typedef Impl::if_c< ( traits::rank == 0 ),
                      fad_view_type ,
                      Impl::ViewError::scalar_operator_called_from_non_scalar_view >
    if_scalar_operator ;

  KOKKOS_FORCEINLINE_FUNCTION
  typename if_scalar_operator::type
    operator()() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return fad_view_type( m_ptr_on_device , m_storage_size-1 , 1 );
    }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 1:

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type , traits, LayoutLeft, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as fad dimension
      return fad_view_type( m_ptr_on_device + m_offset_map(i0,0) ,
                            m_storage_size-1 ,
                            m_offset_map.stride_1() );
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type( m_ptr_on_device + m_offset_map(i0,0),
                            m_storage_size-1 , 1 );
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type , traits, typename traits::array_layout, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 1,
                                      iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_offset_map, i0, i1, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type( m_ptr_on_device + m_offset_map(i0,i1,0) ,
                            m_storage_size-1 ,
                            m_offset_map.stride_2() );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_offset_map, i0, i1, 0);
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 2,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_offset_map, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,0) ,
        m_storage_size-1 ,
        m_offset_map.stride_3() );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_offset_map, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 3,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_offset_map, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,0) ,
        m_storage_size-1 ,
        m_offset_map.stride_4() );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_offset_map, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 4,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_offset_map, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,0) ,
        m_storage_size-1 ,
        m_offset_map.stride_5() );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_offset_map, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 5,
                                      iType0, iType1, iType2, iType3, iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , int , int , int ) const
    { return operator()(i0,i1,i2,i3,i4); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
                 const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_offset_map, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,i5,0) ,
        m_storage_size-1 ,
        m_offset_map.stride_6() );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_offset_map, i0, i1, i2, i3, i4, i5, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,i5,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 6,
                                      iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , int , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_offset_map, i0, i1, i2, i3, i4, i5, i6, 0 );

      // Strided storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,i5,i6,0) ,
        m_storage_size-1 ,
        m_offset_map.stride_7() );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_offset_map, i0, i1, i2, i3, i4, i5, i6, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the fad dimension
      return fad_view_type(
        m_ptr_on_device + m_offset_map(i0,i1,i2,i3,i4,i5,i6,0) ,
        m_storage_size-1 , 1 );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< fad_view_type ,
                                      traits, typename traits::array_layout, 7,
                                      iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5,i6); }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_FORCEINLINE_FUNCTION
  fad_value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  void stride( iType * const s ) const
    { m_offset_map.stride( s ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type capacity() const
    { return size(); }

  // Size (in bytes) of allocatable data
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type data_capacity() const
    { return sizeof(fad_value_type) * m_offset_map.capacity(); }

  // Static storage size
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type storage_size() const
    { return m_storage_size; }
};

/**
 * \brief A deep copy between views of the same specialization, compatible
 * type, same rank, same layout are handled by that specialization.
 *
 * We compare the nested fad_value_type instead of the view value_type as it
 * allows deep_copy to work with views of different, but compatible Fad types.
 */
template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
inline
void deep_copy( const View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad> & dst ,
                const View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad> & src ,
                typename Impl::enable_if<(
                  ( Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>::fad_value_type ,
                                   typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>::fad_value_type >::value ||
                    Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>::const_fad_value_type ,
                                   typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>::fad_value_type >::value ||
                    Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>::fad_value_type ,
                                   typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>::const_fad_value_type >::value )
                  &&
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>::array_layout ,
                                 typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>::array_layout >::value
                  &&
                  ( unsigned(View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>::rank) ==
                    unsigned(View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>::rank) )
                )>::type * = 0 )
{
  typedef  View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoFad>  dst_type ;
  typedef  View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoFad>  src_type ;

  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( dst.ptr_on_device() != src.ptr_on_device() ) {

    Impl::assert_shapes_are_equal( dst.shape() , src.shape() );

    const size_t nbytes = dst.data_capacity();

    Impl::DeepCopy< dst_memory_space , src_memory_space >(
      dst.ptr_on_device() , src.ptr_on_device() , nbytes );
  }
}

 // Overload of deep_copy for Fad views intializing to a constant scalar
template< typename T, typename L, typename D, typename M >
void deep_copy(
  const View<T,L,D,M,Impl::ViewSpecializeSacadoFad>& view ,
  const typename View<T,L,D,M,Impl::ViewSpecializeSacadoFad>::fad_value_type& value )
{
  typedef View<T,L,D,M,Impl::ViewSpecializeSacadoFad> ViewType;
  typedef typename ViewType::fad_value_type ScalarType;
  if (value == ScalarType(0))
    Impl::ViewFill< typename ViewType::array_type >( view , value );
  else
    Impl::ViewFill< ViewType >( view , value );
}

template< class T , class L , class D , class M >
typename Impl::enable_if<(
    View<T,L,D,M,Impl::ViewSpecializeSacadoFad>::is_managed
  ), typename View<T,L,D,M,Impl::ViewSpecializeSacadoFad>::HostMirror >::type
inline
create_mirror( const View<T,L,D,M,Impl::ViewSpecializeSacadoFad> & src )
{
  typedef View<T,L,D,M,Impl::ViewSpecializeSacadoFad>  view_type ;
  typedef typename view_type::HostMirror    host_view_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef typename view_type::size_type     size_type ;

  // 'view' is managed therefore we can allocate a
  // compatible host_view through the ordinary constructor.

  std::string label = memory_space::query_label( src.ptr_on_device() );
  label.append("_mirror");

  size_type dims[8];
  for (size_type i=0; i<8; ++i)
    dims[i] = src.dimension(i);
  dims[unsigned(view_type::Rank)] = src.storage_size();

  return host_view_type( label ,
                         dims[0] ,
                         dims[1] ,
                         dims[2] ,
                         dims[3] ,
                         dims[4] ,
                         dims[5] ,
                         dims[6] ,
                         dims[7] );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct ViewAssignment< ViewSpecializeSacadoFad , ViewSpecializeSacadoFad , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewSpecializeSacadoFad> & dst
                , const View<ST,SL,SD,SM,ViewSpecializeSacadoFad> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_storage_size  = src.m_storage_size ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;
    dst.m_management      = src.m_management ;

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Partition of compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewSpecializeSacadoFad> & dst
                , const View<ST,SL,SD,SM,ViewSpecializeSacadoFad> & src
                , const Sacado::Fad::VectorPartition & part
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ! ViewTraits<DT,DL,DD,DM>::is_managed
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM>                           dst_traits ;
    typedef View<DT,DL,DD,DM,ViewSpecializeSacadoFad>         dst_type ;
    typedef typename dst_traits::value_type                   dst_fad_type ;

    enum { DstRank         = dst_type::Rank };
    enum { DstStaticLength = Sacado::StaticSize<dst_fad_type>::value };

    const int length = part.end - part.begin ;

    if ( DstStaticLength && DstStaticLength != length ) {
      const char msg[] = "Kokkos::View< Fad ... > incompatible partitioning" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort(msg);
#else
      throw std::runtime_error(msg);
#endif
    }

    // Copy the offset map:
    dst.m_offset_map.assign( src.m_offset_map );

    // Override the last dimension of the offset map:
    dst.m_offset_map.assign<DstRank>( length );

    dst.m_storage_size = src.m_storage_size ;

    dst.m_ptr_on_device = src.m_ptr_on_device + part.begin * (
      ( 0 == DstRank ? dst.m_offset_map.stride_0() :
      ( 1 == DstRank ? dst.m_offset_map.stride_1() :
      ( 2 == DstRank ? dst.m_offset_map.stride_2() :
      ( 3 == DstRank ? dst.m_offset_map.stride_3() :
      ( 4 == DstRank ? dst.m_offset_map.stride_4() :
      ( 5 == DstRank ? dst.m_offset_map.stride_5() :
      ( 6 == DstRank ? dst.m_offset_map.stride_6() :
      ( 7 == DstRank ? dst.m_offset_map.stride_7() : 0 )))))))));
  }
};


template<>
struct ViewAssignment< ViewDefault , ViewSpecializeSacadoFad , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,ViewSpecializeSacadoFad>::array_type & dst
                , const    View<ST,SL,SD,SM,ViewSpecializeSacadoFad> & src )
  {
    typedef View<ST,SL,SD,SM,ViewSpecializeSacadoFad>  src_type ;
    typedef typename src_type::array_type  dst_type ;

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_ptr_on_device = reinterpret_cast< typename dst_type::value_type *>( src.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_management.increment( dst.m_ptr_on_device );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

#endif /* #ifndef KOKKOS_VIEW_FAD_HPP */
