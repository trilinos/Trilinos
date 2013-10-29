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

#ifndef SACADO_VEIW_MP_VECTOR_HPP
#define SACADO_VEIW_MP_VECTOR_HPP

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
  typename enable_if< StorageType::is_static , typename StorageType::device_type::memory_space >::type ,
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
  typename enable_if< StorageType::is_static , typename StorageType::device_type::memory_space >::type ,
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
  typename enable_if< StorageType::is_static , typename StorageType::device_type::memory_space >::type ,
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
  typename enable_if< StorageType::is_static , typename StorageType::device_type::memory_space >::type ,
  MemoryTraits >
{
  typedef ViewSpecializeSacadoMPVector type ;
};

template< class T , class Device > struct RebindStokhosStorageDevice ;

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/**\brief  View::value_type  == Sacado::MP::Vector< StorageType >
 *         View::scalar_type == StorageType::value_type
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

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  typename traits::scalar_type * m_ptr_on_device ;
  typename traits::shape_type    m_shape ;
  stride_type                    m_stride ;



  typedef typename traits::value_type                   sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type  stokhos_storage_type ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::static_size };

  typedef Stokhos::ViewStorage<
    typename stokhos_storage_type::ordinal_type ,
    typename traits::scalar_type ,
    StokhosStorageStaticDimension ,
      /* LayoutRight has stride-one stokhos storage */
    ( Impl::is_same< typename traits::array_layout , LayoutRight >::value ? 1 : 0 ) ,
    typename traits::device_type >  stokhos_view_storage_type ;

  typedef Sacado::MP::Vector< stokhos_view_storage_type >  sacado_mp_vector_view_type ;

public:

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename Impl::RebindStokhosStorageDevice<
                  typename traits::data_type ,
                  typename traits::device_type::host_mirror_device_type >::type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  // Rank for multidimensional array of the intrinsic scalar_type,
  // not the Sacado::MP::Vector value_type.
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

private:

  template< class ViewRHS >
  KOKKOS_INLINE_FUNCTION
  void assign_compatible_view( const ViewRHS & rhs ,
                               typename Impl::enable_if< Impl::ViewAssignable< View , ViewRHS >::value >::type * = 0 )
  {
    typedef typename traits::shape_type    shape_type ;
    typedef typename traits::memory_space  memory_space ;
    typedef typename traits::memory_traits memory_traits ;

    Impl::ViewTracking< traits >::decrement( m_ptr_on_device );

    shape_type::assign( m_shape,
                        rhs.m_shape.N0 , rhs.m_shape.N1 , rhs.m_shape.N2 , rhs.m_shape.N3 ,
                        rhs.m_shape.N4 , rhs.m_shape.N5 , rhs.m_shape.N6 , rhs.m_shape.N7 );

    stride_type::assign( m_stride , rhs.m_stride.value );

    m_ptr_on_device = rhs.m_ptr_on_device ;

    Impl::ViewTracking< traits >::increment( m_ptr_on_device );
  }

  // Restrict allocation to a multiple of 'StokhosStorageStaticDimension'
  inline
  void verify_dimension_storage_static_size() const
  {
    if ( dimension( Rank - 1 ) % StokhosStorageStaticDimension ) {
      Impl::throw_runtime_exception( std::string("Kokkos::View< Sacado::MP::Vector<StorageType , ... > allocation dimension must be multple of StorageType::static_size" ) );
    }
  }

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
  View( const View & rhs ) : m_ptr_on_device(0) { assign_compatible_view( rhs ); }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign_compatible_view( rhs ); return *this ; }

  //------------------------------------
  // Construct or assign compatible view:

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    : m_ptr_on_device(0) { assign_compatible_view( rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,typename traits::specialize> & rhs )
    { assign_compatible_view( rhs ); return *this ; }

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

      // Restrict allocation to a multiple of 'StokhosStorageStaticDimension'

      verify_dimension_storage_static_size();

      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );

      Impl::ViewInitialize< typename traits::device_type > init( *this );
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

      // Restrict allocation to a multiple of 'StokhosStorageStaticDimension'

      verify_dimension_storage_static_size();

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

      // Restrict dimension to a multiple of 'StokhosStorageStaticDimension'

      verify_dimension_storage_static_size();

      m_ptr_on_device = if_user_pointer_constructor::select( ptr );
    }

  //------------------------------------
  // Is not allocated

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 1:

  template< typename T >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( Impl::is_same< T , typename traits::device_type >::value &&
                             Impl::is_same< LayoutLeft, typename traits::array_layout >::value 
                           ), sacado_mp_vector_view_type >::type
    operator() ( const T & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + dev.team_rank() * StokhosStorageStaticDimension ,
        m_shape.N0 , 1 ) );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutLeft, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutLeft, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutLeft, 1, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 2:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type , traits, LayoutLeft, 2, iType0 >::type
    operator() ( const iType0 & i0 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + i0 + m_stride.value * ( dev.team_rank() * StokhosStorageStaticDimension ) ,
        m_shape.N1 , m_stride.value ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * i1 ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 2, iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * i1 ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_stride.N1 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ),
        m_shape.N2 ,
        m_stride.value * m_shape.N1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 3, iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ) ),
        m_shape.N3 ,
        m_stride.value * m_shape.N1 * m_shape.N2 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ) ) ),
        m_shape.N4 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 , iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3 , iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 6:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ) ) ) ),
        m_shape.N5 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3 , iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 + m_stride.value * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 7:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
                 const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 + m_shape.N5 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ) ) ) ) ),
        m_shape.N6 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3 , iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

  //------------------------------------
  //------------------------------------
  // LayoutLeft, array operators, rank 8:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Strided storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( i0 + m_stride.value * (
                            i1 + m_shape.N1 * (
                            i2 + m_shape.N2 * (
                            i3 + m_shape.N3 * (
                            i4 + m_shape.N4 * (
                            i5 + m_shape.N5 * (
                            i6 + m_shape.N6 * (
                            dev.team_rank() * StokhosStorageStaticDimension ) ) ) ) ) ) ) ),
        m_shape.N7 ,
        m_stride.value * m_shape.N1 * m_shape.N2 * m_shape.N3 * m_shape.N4 * m_shape.N5 * m_shape.N6 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3 , iType4, iType5, iType6, iType7 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  //------------------------------------
  // LayoutRight, array operators, rank 1:

  template< typename T >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( Impl::is_same< T , typename traits::device_type >::value &&
                             Impl::is_same< LayoutRight , typename traits::array_layout >::value 
                           ), sacado_mp_vector_view_type >::type
    operator() ( const T & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension ) ,
        m_shape.N0 , 1 ) );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutRight, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutRight, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & , traits, LayoutRight, 1, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 2:

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 2, iType0 >::type
    operator() ( const iType0 & i0 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension + m_stride.value * i0 ) ,
        m_shape.N1 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride.value ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 3:

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 3, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension + 
                            m_shape.N2 * ( i1 ) + m_stride.value * i0 ) ,
        m_shape.N2 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * ( i1 ) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * ( i1 ) + i0 * m_stride.value ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 4:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )) +
                            m_stride.value * i0 ) ,
        m_shape.N3 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride.value ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))) +
                            m_stride.value * i0 ) ,
        m_shape.N4 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride.value ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 6:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))) +
                            m_stride.value * i0 ) ,
        m_shape.N5 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride.value ];
    }

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 7:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension +
                            m_shape.N6 * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 ))))) +
                            m_stride.value * i0 ) ,
        m_shape.N6 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

  //------------------------------------
  //------------------------------------
  // LayoutRight, array operators, rank 8:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< sacado_mp_vector_view_type ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ,
                 const typename traits::device_type & dev ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, dev.team_rank() * StokhosStorageStaticDimension );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      // Contiguous storage with right-most index as the stokhos dimension
      return sacado_mp_vector_view_type( stokhos_view_storage_type(
        m_ptr_on_device + ( dev.team_rank() * StokhosStorageStaticDimension +
                            m_shape.N7 * ( i6 +
                            m_shape.N6 * ( i5 +
                            m_shape.N5 * ( i4 +
                            m_shape.N4 * ( i3 +
                            m_shape.N3 * ( i2 +
                            m_shape.N2 * ( i1 )))))) +
                            m_stride.value * i0 ) ,
        m_shape.N7 ,
        1 ) );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type & ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
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

template< class T , class Device >
struct RebindStokhosStorageDevice< T * , Device >
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


} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef SACADO_VIEW_MP_VECTOR_HPP */

