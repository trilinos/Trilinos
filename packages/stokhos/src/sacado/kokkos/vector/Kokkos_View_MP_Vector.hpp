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

namespace Sacado {
namespace MP {

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

struct ViewSpecializeSacadoMPVector {};
struct ViewSpecializeSacadoMPVectorStatic {};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType 
  , ViewSpecializeSacadoMPVector
  , LayoutLeft
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType 
  , ViewSpecializeSacadoMPVector
  , LayoutRight
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType 
  , ViewSpecializeSacadoMPVectorStatic
  , LayoutLeft
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVectorStatic type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType 
  , ViewSpecializeSacadoMPVectorStatic
  , LayoutRight
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVectorStatic type ;
};


template< class T , class Device > struct RebindStokhosStorageDevice ;

//----------------------------------------------------------------------------

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

  typename stokhos_storage_type::value_type  * m_ptr_on_device ;
  typename traits::shape_type                  m_shape ;
  stride_type                                  m_stride ;
  typename traits::device_type::size_type      m_storage_size ;

  typedef Stokhos::ViewStorage<
    typename stokhos_storage_type::ordinal_type ,
    typename stokhos_storage_type::value_type ,
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
      typedef typename traits::memory_space              memory_space ;
      typedef typename traits::shape_type                shape_type ;
      typedef typename stokhos_storage_type::value_type  scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_shape );

      verify_dimension_storage_static_size();

      m_storage_size  = Impl::dimension( m_shape , unsigned(Rank) );
      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );

      (void) Impl::ViewFill< array_type >( *this , typename array_type::value_type() );
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
      typedef typename traits::memory_space              memory_space ;
      typedef typename traits::shape_type                shape_type ;
      typedef typename stokhos_storage_type::value_type  scalar_type ;

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

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_no_padding( m_stride , m_shape );

      verify_dimension_storage_static_size();

      m_storage_size = Impl::dimension( m_shape , unsigned(Rank) );
    }

  //------------------------------------
  // Is not allocated

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // Scalar operator on traits::rank == 1

  typedef Impl::if_c< ( traits::rank == 1 ),
                      sacado_mp_vector_view_type ,
                      Impl::ViewError::scalar_operator_called_from_non_scalar_view >
    if_scalar_operator ;

  KOKKOS_FORCEINLINE_FUNCTION
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
  // Array operators, traits::rank 2:

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type , traits, LayoutLeft, 2, iType0 >::type
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type , traits, typename traits::array_layout, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, typename traits::array_layout, 2,
                                      iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, typename traits::array_layout, 3,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, typename traits::array_layout, 4,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, typename traits::array_layout, 5,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
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
    ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::capacity( m_shape , m_stride ); }

  // Static storage size
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type static_storage_size() const
  { return StokhosStorageStaticDimension; }
};

/** \brief  A deep copy between views of the same specialization, compatible type,
 *          same rank, same layout are handled by that specialization.
 */
template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
inline
void deep_copy( const View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVector> & dst ,
                const View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVector> & src ,
                typename Impl::enable_if<(
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVector>::value_type::storage_type::value_type ,
                                 typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVector>::value_type::storage_type::value_type >::value
                  &&
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVector>::array_layout ,
                                 typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVector>::array_layout >::value
                  &&
                  ( unsigned(View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVector>::rank) ==
                    unsigned(View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVector>::rank) )
                )>::type * = 0 )
{
  typedef  View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVector>  dst_type ;
  typedef  View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVector>  src_type ;

  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( dst.ptr_on_device() != src.ptr_on_device() ) {

    Impl::assert_shapes_are_equal( dst.shape() , src.shape() );

    const size_t nbytes = sizeof(typename dst_type::value_type::storage_type::value_type) * dst.capacity();

    Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.ptr_on_device() , src.ptr_on_device() , nbytes );
  }
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/**\brief  View::value_type  == Sacado::MP::Vector< Stokhos::FixedStaticStorage<...> > */
template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::ViewSpecializeSacadoMPVectorStatic >
  : public ViewTraits< DataType
                     , LayoutRight
                     , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::device_type
                     , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::memory_traits
                     >
{
public:

  typedef ViewTraits< DataType
                    , LayoutRight
                    , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::device_type
                    , typename ViewTraits< DataType , Arg1Type, Arg2Type, Arg3Type >::memory_traits
                    > traits ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef typename traits::value_type                   sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type  stokhos_storage_type ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::static_size };

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  typename traits::value_type  * m_ptr_on_device ;
  typename traits::shape_type    m_shape ;
  unsigned                       m_stride ;

  // original_static_dimension = m_stride * StokhosStorageStaticDimension
  // m_stride = 1 for original allocation.

public:

  // This needs to be public so that we know what the return type of () is
  typedef Sacado::MP::Vector< stokhos_storage_type >  sacado_mp_vector_view_type ;

  // Whether the storage type is statically sized
  static const bool is_static = true ;

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
  // Shape for the Sacado::MP::Vector value_type ignores the internal static array length.
  enum { Rank = traits::rank };

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
      m_stride = 0 ;
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
      typedef typename traits::value_type    value_type ;

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
      m_ptr_on_device = (value_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(value_type) ,
                                sizeof(value_type) ,
                                Impl::cardinality_count( m_shape ) );

      (void) Impl::ViewFill< View >( *this , typename traits::value_type() );
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
      typedef typename traits::value_type    value_type ;

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
      m_ptr_on_device = (value_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(value_type) ,
                                sizeof(value_type) ,
                                Impl::cardinality_count( m_shape ) );
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
      typedef typename traits::shape_type   shape_type ;
      typedef typename traits::value_type   value_type ;

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
    }

  //------------------------------------
  // Is not allocated

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // Scalar operator on traits::rank == 0

  typedef Impl::if_c< ( traits::rank == 0 ) ,
                      typename traits::value_type ,
                      Impl::ViewError::scalar_operator_called_from_non_scalar_view >
    if_scalar_operator ;

  KOKKOS_FORCEINLINE_FUNCTION
  typename if_scalar_operator::type &
    operator()() const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return if_scalar_operator::select( *m_ptr_on_device );
    }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 1:

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, LayoutRight, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      // May have partitioned
      return m_ptr_on_device[ m_stride * i0 ];
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, LayoutRight, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, LayoutRight, 1, iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * ( i1 + m_shape.N1 * i0 ) ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 2,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * ( i2 + m_shape.N2 * ( i1 + m_shape.N1 * i0 ) ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 3,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, i3 );

      return m_ptr_on_device[ m_stride * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 + m_shape.N1 * ( i0 )))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 4,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, i4 );

      return m_ptr_on_device[ m_stride * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 + m_shape.N1 * ( i0 ))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 5,
                                      iType0, iType1, iType2, iType3, iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , int , int , int ) const
    { return operator()(i0,i1,i2,i3,i4); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 6:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, i5 );

      return m_ptr_on_device[ m_stride * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 + m_shape.N1 * ( i0 )))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 6,
                                      iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const int , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 7:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, i6 );

      return m_ptr_on_device[ m_stride * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 + m_shape.N1 * ( i0 ))))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & ,
                                      traits, LayoutRight, 7,
                                      iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5,i6); }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::cardinality_count( m_shape ) * m_stride ; }

  // Static storage size
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type static_storage_size() const
  { return StokhosStorageStaticDimension; }
};

/** \brief  A deep copy between views of the same specialization, compatible type,
 *          same rank, same layout are handled by that specialization.
 */
template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
inline
void deep_copy( const View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVectorStatic> & dst ,
                const View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVectorStatic> & src ,
                typename Impl::enable_if<(
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVectorStatic>::value_type::storage_type::value_type ,
                                 typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVectorStatic>::value_type::storage_type::value_type >::value
                  &&
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVectorStatic>::array_layout ,
                                 typename View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVectorStatic>::array_layout >::value
                  &&
                  ( unsigned(View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVectorStatic>::rank) ==
                    unsigned(View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVectorStatic>::rank) )
                )>::type * = 0 )
{
  typedef  View<DT,DL,DD,DM,Impl::ViewSpecializeSacadoMPVectorStatic>  dst_type ;
  typedef  View<ST,SL,SD,SM,Impl::ViewSpecializeSacadoMPVectorStatic>  src_type ;

  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( (void*) dst.ptr_on_device() != (void*) src.ptr_on_device() ) {

    Impl::assert_shapes_are_equal( dst.shape() , src.shape() );

    const size_t nbytes = sizeof(typename dst_type::value_type) * dst.capacity();

    Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.ptr_on_device() , src.ptr_on_device() , nbytes );
  }
}

} // namespace Kokkos

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
  : if_c< StorageType::is_static
        , Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > // Treat as a scalar
        , typename ShapeInsert< typename AnalyzeShape< typename StorageType::value_type >::shape , 0 >::type
        >::type
{
private:

  typedef AnalyzeShape< typename StorageType::value_type > nested ;

public:

  typedef typename
    if_c< StorageType::is_static 
        , ViewSpecializeSacadoMPVectorStatic
        , ViewSpecializeSacadoMPVector
        >::type specialize ;

  typedef typename
    if_c< StorageType::is_static
        , Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 >
        , typename ShapeInsert< typename nested::shape , 0 >::type 
        >::type shape ;

  // If ( ! StorageType::is_static ) then 0 == StorageType::static_size and the first array declaration is not used.
  // However, the compiler will still generate this type declaration and it must not have a zero length.
  typedef typename
    if_c< StorageType::is_static
        , typename nested::array_type [ StorageType::is_static ? StorageType::static_size : 1 ]
        , typename nested::array_type *
        >::type array_type ;

  typedef typename
    if_c< StorageType::is_static
        , typename nested::const_array_type [ StorageType::is_static ? StorageType::static_size : 1 ]
        , typename nested::const_array_type *
        >::type const_array_type ;

  typedef array_type non_const_array_type ;

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
    typedef ViewTraits<DT,DL,DD,DM>                         dst_traits ;
    typedef View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector>  dst_type ;
    typedef typename dst_type::shape_type                   shape_type ;
    typedef typename dst_type::stride_type                  stride_type ;

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_storage_size  = src.m_storage_size ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Partition of compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector> & dst
                , const View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector> & src
                , const Sacado::MP::VectorPartition & part
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ! ViewTraits<DT,DL,DD,DM>::is_managed
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM>                           dst_traits ;
    typedef View<DT,DL,DD,DM,ViewSpecializeSacadoMPVector>    dst_type ;
    typedef typename dst_type::shape_type                     dst_shape_type ;
    typedef typename dst_type::stride_type                    dst_stride_type ;
    typedef typename dst_traits::value_type                   dst_sacado_mp_vector_type ;
    typedef typename dst_sacado_mp_vector_type::storage_type  dst_stokhos_storage_type ;

    enum { DstRank         = dst_type::Rank };
    enum { DstStaticLength = dst_stokhos_storage_type::is_static ? dst_stokhos_storage_type::static_size : 0 };

    const int length = part.end - part.begin ;

    if ( DstStaticLength && DstStaticLength != length ) {
      const char msg[] = "Kokkos::View< Sacado::MP::Vector ... > incompatible partitioning" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort(msg);
#else
      throw std::runtime_error(msg);
#endif
    }

    dst_shape_type::assign( dst.m_shape ,
                            ( DstRank == 0 ? length : src.m_shape.N0 ) ,
                            ( DstRank == 1 ? length : src.m_shape.N1 ) ,
                            ( DstRank == 2 ? length : src.m_shape.N2 ) ,
                            ( DstRank == 3 ? length : src.m_shape.N3 ) ,
                            ( DstRank == 4 ? length : src.m_shape.N4 ) ,
                            ( DstRank == 5 ? length : src.m_shape.N5 ) ,
                            ( DstRank == 6 ? length : src.m_shape.N6 ) ,
                            ( DstRank == 7 ? length : src.m_shape.N7 ) );

    dst_stride_type::assign( dst.m_stride , src.m_stride.value );

    // Original Sacado::MP::Vector length
    dst.m_storage_size = src.m_storage_size ;

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
  }
};

template<>
struct ViewAssignment< ViewSpecializeSacadoMPVectorStatic , ViewSpecializeSacadoMPVectorStatic , void >
{
  typedef ViewSpecializeSacadoMPVectorStatic specialize ;

  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst
                , const View<ST,SL,SD,SM,specialize> & src
                , const typename enable_if<(
                    ViewAssignable< View<DT,DL,DD,DM,specialize> ,
                                    View<ST,SL,SD,SM,specialize> >::value
                    )>::type * = 0
                  )
  {
    typedef View<DT,DL,DD,DM,specialize>   dst_type ;
    typedef typename dst_type::shape_type  shape_type ;

    ViewTracking< dst_type >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_type >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Partition of compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst
                , const View<ST,SL,SD,SM,specialize> & src
                , typename enable_if<(
                    // Same intrinsic value type
                    is_same< typename View<DT,DL,DD,DM,specialize>::value_type::storage_type::value_type ,
                             typename View<ST,SL,SD,SM,specialize>::value_type::storage_type::value_type >::value
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
                    ( unsigned(View<DT,DL,DD,DM,specialize>::Rank) ==
                      unsigned(View<ST,SL,SD,SM,specialize>::Rank) )
                    &&
                    // Destination is not managed
                    ! View<DT,DL,DD,DM,specialize>::is_managed
                  ), const Sacado::MP::VectorPartition & >::type part )
  {
    typedef View<ST,SL,SD,SM,specialize>   src_type ;
    typedef View<DT,DL,DD,DM,specialize>   dst_type ;
    typedef typename dst_type::shape_type  dst_shape_type ;
    typedef typename dst_type::value_type  dst_sacado_mp_vector_type ;
    typedef typename src_type::value_type  src_sacado_mp_vector_type ;

    typedef typename dst_sacado_mp_vector_type::storage_type  dst_stokhos_storage_type ;
    typedef typename src_sacado_mp_vector_type::storage_type  src_stokhos_storage_type ;

    // Must have: begin = i * dst_stokhos_storage_type::static_size
    //            end   = begin + dst_stokhos_storage_type::static_size
    //
    if ( part.begin % dst_stokhos_storage_type::static_size ||
         part.begin + dst_stokhos_storage_type::static_size != part.end ) {
      const char msg[] = "Kokkos::View< Sacado::MP::Vector ... > incompatible partitioning" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort(msg);
#else
      throw std::runtime_error(msg);
#endif
    }

    dst_shape_type::assign( dst.m_shape ,
                            src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                            src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst.m_stride = src.m_stride * src_stokhos_storage_type::static_size
                                / dst_stokhos_storage_type::static_size ;

    dst.m_ptr_on_device = reinterpret_cast<typename dst_type::value_type*>( src.m_ptr_on_device )
                        + part.begin / dst_stokhos_storage_type::static_size ;
  }
};


template<>
struct ViewAssignment< ViewDefault , ViewSpecializeSacadoMPVector , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector>::array_type & dst
                , const    View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector> & src )
  {
    typedef View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector>    src_type ;
    typedef typename src_type::value_type                     src_sacado_mp_vector_type ;
    typedef typename src_sacado_mp_vector_type::storage_type  src_stokhos_storage_type ;

    typedef typename src_type::array_type   dst_type ;
    typedef typename dst_type::shape_type   dst_shape_type ;
    typedef typename dst_type::stride_type  dst_stride_type ;

    ViewTracking< dst_type >::decrement( dst.m_ptr_on_device );

    dst_shape_type::assign( dst.m_shape,
                            src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                            src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst_stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = reinterpret_cast< typename dst_type::value_type *>( src.m_ptr_on_device );

    Impl::ViewTracking< dst_type >::increment( dst.m_ptr_on_device );
  }
};

template<>
struct ViewAssignment< ViewDefault , ViewSpecializeSacadoMPVectorStatic , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,ViewSpecializeSacadoMPVectorStatic>::array_type & dst
                , const    View<ST,SL,SD,SM,ViewSpecializeSacadoMPVectorStatic> & src )
  {
    typedef View<ST,SL,SD,SM,ViewSpecializeSacadoMPVector>    src_type ;
    typedef typename src_type::value_type                     src_sacado_mp_vector_type ;
    typedef typename src_sacado_mp_vector_type::storage_type  src_stokhos_storage_type ;

    typedef typename src_type::array_type   dst_type ;
    typedef typename dst_type::shape_type   dst_shape_type ;
    typedef typename dst_type::stride_type  dst_stride_type ;

    if ( src.m_stride != 1 ) {
      const char msg[] = "Kokkos::View< Sacado::MP::Vector ... > incompatible assignment" ;
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      cuda_abort(msg);
#else
      throw std::runtime_error(msg);
#endif
    }

    ViewTracking< dst_type >::decrement( dst.m_ptr_on_device );

    dst_shape_type::assign( dst.m_shape,
                            src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                            src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst_stride_type::assign_no_padding( dst.m_stride , dst.m_shape );

    dst.m_ptr_on_device = reinterpret_cast< typename dst_type::value_type *>( src.m_ptr_on_device );

    Impl::ViewTracking< dst_type >::increment( dst.m_ptr_on_device );
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
