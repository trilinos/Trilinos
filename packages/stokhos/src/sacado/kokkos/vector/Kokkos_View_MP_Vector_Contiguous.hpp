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

#ifndef KOKKOS_VIEW_MP_VECTOR_CONTIGUOUS_HPP
#define KOKKOS_VIEW_MP_VECTOR_CONTIGUOUS_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_MP_Vector.hpp"
#include "Sacado_MP_VectorTraits.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_AnalyzeSacadoShape.hpp"
#include "Kokkos_View_Utils.hpp"
#include "Kokkos_View_MP_Vector_Utils.hpp"

/*
 * Specialization for Kokkos::View<Sacado::MP::Vector<Storage>...>
 * where the Sacado dimension is ALWAYS kept contiguous regardless of
 * layout.
 *
 * Currently it can't be used at the same time as other such View
 * specializations due to conflicting specializations of AnalyzeShape.
 */

namespace Kokkos {
namespace Impl {

struct ViewMPVectorContiguous {};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewMPVectorContiguous
  , LayoutLeft
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewMPVectorContiguous type ;
};

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize
  < ValueType
  , ViewMPVectorContiguous
  , LayoutRight
  , MemorySpace
  , MemoryTraits >
{
  typedef ViewMPVectorContiguous type ;
};

//----------------------------------------------------------------------------

template < class Storage, bool is_static = Storage::is_static >
struct MPVectorAllocation;

template <class Storage>
struct MPVectorAllocation<Storage, true> {
  typedef Sacado::MP::Vector<Storage> value_type;
  typedef typename Storage::value_type scalar_type;

  scalar_type * m_scalar_ptr_on_device;
  Kokkos::Impl::AllocationTracker m_tracker;

  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation() : m_scalar_ptr_on_device(0), m_tracker() {}

  // Allocate scalar_type and value_type arrays
  template <class ExecSpace, class MemSpace, bool Initialize,
            class LabelType, class ShapeType>
  inline
  value_type*
  allocate(const LabelType& label,
           const ShapeType& shape,
           const unsigned vector_size) {
    typedef MemSpace memory_space;
    typedef ExecSpace execution_space;

    const size_t count = Impl::cardinality_count( shape ) * vector_size;
    m_tracker = memory_space::allocate_and_track( label ,
                                                  sizeof(scalar_type)*count );
    m_scalar_ptr_on_device =
      reinterpret_cast<scalar_type *>(m_tracker.alloc_ptr());

    // Initialize data
    if (count > 0) {
      (void) Impl::ViewDefaultConstruct< execution_space , scalar_type , Initialize >( m_scalar_ptr_on_device , count );
    }

    return reinterpret_cast<value_type*>(m_scalar_ptr_on_device);
  }

  // Assign scalar_type pointer to give ptr
  template <typename T>
  void assign(T * ptr) {
    m_scalar_ptr_on_device = (scalar_type*)ptr;
    m_tracker.clear();
  }

};

template <class Storage>
struct MPVectorAllocation<Storage, false> {
  typedef Sacado::MP::Vector<Storage> value_type;
  typedef typename Storage::value_type scalar_type;

  scalar_type * m_scalar_ptr_on_device;
  Kokkos::Impl::AllocationTracker m_tracker;

  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation() : m_scalar_ptr_on_device(), m_tracker() {}

  // Allocate scalar_type and value_type arrays
  template <class ExecSpace, class MemSpace, bool Initialize,
            class LabelType, class ShapeType>
  inline
  value_type*
  allocate(const LabelType& label,
           const ShapeType& shape,
           const unsigned vector_size) {
    typedef MemSpace memory_space;
    typedef ExecSpace execution_space;

    // Allocate space for contiguous MP::Vector values
    // and for MP::Vector itself.  We do this in one
    // chunk so that Kokkos' memory tracking works
    // properly.  However by doing it this way
    // we could run into alignment issues.  Not sure if
    // this is the best choice from a locality perspective
    // either.
    const size_t num_vec = Impl::cardinality_count( shape );
    const size_t count_scalars = num_vec * vector_size;
    const size_t size_scalars = count_scalars * sizeof(scalar_type);
    const size_t size_values = num_vec * sizeof(value_type);
    const size_t size = size_scalars + size_values;
    m_tracker = memory_space::allocate_and_track( label, size );
    char *data = reinterpret_cast<char*>(m_tracker.alloc_ptr());
    m_scalar_ptr_on_device = (scalar_type *) data;
    value_type * ptr = (value_type *) (data + size_scalars);

     // Initialize data
    (void) Impl::ViewDefaultConstruct< execution_space , scalar_type , Initialize >( m_scalar_ptr_on_device , count_scalars );

    // Construct each MP::Vector using memory in ptr array,
    // setting pointer to MP::Vector values from values array
    // Equivalent to:
    // value_type* p = ptr;
    // scalar_type* sp = m_scalar_ptr_on_device;
    // for (size_t i=0; i<num_vec; ++i) {
    //   new (p++) value_type(vector_size, sp, false);
    //   sp += vector_size;
    // }
    parallel_for( num_vec,
                  VectorInit<execution_space>( ptr, m_scalar_ptr_on_device,
                                               vector_size ) );

    return ptr;
  }

  // Assign scalar_type pointer to give ptr
  // This makes BIG assumption on how the data was allocated
  template <typename T>
  void assign(T * ptr) {
    if (ptr != 0) {
      m_scalar_ptr_on_device = ((value_type *)ptr)->coeff();
      m_tracker.clear();
    }
    else {
      m_scalar_ptr_on_device = 0;
      m_tracker.clear();
    }
  }

  template <class ExecSpace>
  struct VectorInit {
    typedef ExecSpace execution_space;
    value_type* p;
    scalar_type* sp;
    const unsigned vector_size;

    KOKKOS_INLINE_FUNCTION
    VectorInit(value_type* p_, scalar_type* sp_, const unsigned vector_size_) :
      p(p_), sp(sp_), vector_size(vector_size_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_t i) const {
      new (p+i) value_type(vector_size, sp+i*vector_size, false);
    }
  };
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

// Overload of deep_copy for MP::Vector views intializing to a constant scalar
template< typename T, typename L, typename D, typename M >
void deep_copy(
  const View<T,L,D,M,Impl::ViewMPVectorContiguous>& view ,
  const typename View<T,L,D,M,Impl::ViewMPVectorContiguous>::intrinsic_scalar_type& value )
{
  typedef View<T,L,D,M,Impl::ViewMPVectorContiguous> ViewType;
  Impl::ViewFill< typename ViewType::flat_array_type >( view , value );
}

/**\brief View::value_type  == Sacado::MP::Vector< Storage<...> > */
template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::ViewMPVectorContiguous >
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

  typedef typename traits::value_type                   sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type  stokhos_storage_type ;
  typedef typename stokhos_storage_type::value_type     intrinsic_scalar_type ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::static_size };
  typedef Sacado::integral_nonzero< unsigned , StokhosStorageStaticDimension > sacado_size_type;

  typedef Impl::ViewOffset< typename traits::shape_type ,
                            typename traits::array_layout > offset_map_type ;

  typedef Impl::AnalyzeSacadoShape< typename traits::data_type,
                                    typename traits::array_layout > analyze_sacado_shape;

  typedef Impl::MPVectorAllocation<stokhos_storage_type> allocation_type;

  typename traits::value_type           * m_ptr_on_device ;
  allocation_type                         m_allocation;
  offset_map_type                         m_offset_map ;
  unsigned                                m_stride ;
  typename traits::execution_space::size_type m_storage_size ; // Storage size of sacado dimension
  sacado_size_type                        m_sacado_size ; // Size of sacado dimension
  Impl::ViewDataManagement< traits >      m_management ;

  // Note:  if the view is partitioned, m_sacado_size != m_storage_size.
  // We always have m_storage_size >= m_sacado_size

  // original_sacado_size = m_stride * m_sacado_size
  // m_stride = 1 for original allocation.

  // Stride here has a slightly different meaning than in the standard
  // View implementation.  For the moment we are assuming no padding within
  // the view array itself and stride is to allow for partitioning the view
  // by dividing up the scalar type.
  //
  // I suspect we could combine this with the way the stride is managed in
  // the default view, in which case, I don't think we even need a
  // specialization
  //
  // For reshaping by folding the sacado dimension into its next adjacent
  // dimension, padding wouldn't generally work.  So unless there becomes
  // a way to turn padding off in the default view, a specialization
  // will be necessary.

public:

  // The return type of operator()
  typedef typename traits::value_type & reference_type ;

  // Whether the storage type is statically sized
  static const bool is_static = stokhos_storage_type::is_static ;

  // Whether sacado dimension is contiguous
  static const bool is_contiguous = true;

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
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                void > HostMirror ;

  // Equivalent array type for this view.
  typedef View< typename analyze_sacado_shape::array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > array_type ;

  // Equivalent const array type for this view.
  typedef View< typename analyze_sacado_shape::const_array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_array_type ;

  // Equivalent host array type for this view.
  typedef View< typename analyze_sacado_shape::array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                typename traits::memory_traits > host_array_type ;

  // Equivalent const host array type for this view.
  typedef View< typename analyze_sacado_shape::const_array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                typename traits::memory_traits > host_const_array_type ;

  // Equivalent flattened array type for this view.
  typedef View< typename analyze_sacado_shape::flat_array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > flat_array_type ;

  // Equivalent const flattened array type for this view.
  typedef View< typename analyze_sacado_shape::const_flat_array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_flat_array_type ;

  //------------------------------------
  // Shape for the Sacado::MP::Vector value_type ignores the internal static array length.
  enum { Rank = traits::rank };

  // Rank corresponding to the sacado dimension
  enum { Sacado_Rank = Impl::is_same< typename traits::array_layout, LayoutLeft >::value ? 0 : Rank+1 };

  KOKKOS_FORCEINLINE_FUNCTION typename traits::shape_type shape() const { return m_offset_map ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_0() const { return m_offset_map.N0 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_1() const { return m_offset_map.N1 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_2() const { return m_offset_map.N2 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_3() const { return m_offset_map.N3 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_4() const { return m_offset_map.N4 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_5() const { return m_offset_map.N5 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_6() const { return m_offset_map.N6 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type dimension_7() const { return m_offset_map.N7 ; }
  KOKKOS_FORCEINLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_offset_map.N0
           * m_offset_map.N1
           * m_offset_map.N2
           * m_offset_map.N3
           * m_offset_map.N4
           * m_offset_map.N5
           * m_offset_map.N6
           * m_offset_map.N7
           ;
  }

  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_offset_map , i ); }

  // Dimensions of view, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  void dimensions( iType * const dims ) const
  {
    if (Rank >= 1) dims[0] = m_offset_map.N0;
    if (Rank >= 2) dims[1] = m_offset_map.N1;
    if (Rank >= 3) dims[2] = m_offset_map.N2;
    if (Rank >= 4) dims[3] = m_offset_map.N3;
    if (Rank >= 5) dims[4] = m_offset_map.N4;
    if (Rank >= 6) dims[5] = m_offset_map.N5;
    if (Rank >= 7) dims[6] = m_offset_map.N6;
    if (Rank >= 8) dims[7] = m_offset_map.N7;
  }

  //------------------------------------

public:

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() {}

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0), m_storage_size(0), m_sacado_size(0)
    {
      m_offset_map.assign(0,0,0,0,0,0,0,0);
      m_stride = 1 ; // need to initialize to 1 as there are checks for
                     // m_stride != 1 that need to work for empty views
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
      typedef typename traits::execution_space execution_space;
      typedef typename traits::memory_space memory_space;

      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
      m_storage_size =
        Impl::GetSacadoSize<unsigned(Rank)>::eval(n0,n1,n2,n3,n4,n5,n6,n7);
      if (m_storage_size == 0)
        m_storage_size = global_sacado_mp_vector_size;
      m_sacado_size = m_storage_size;
      m_ptr_on_device =
        m_allocation.template allocate<execution_space,memory_space,Alloc::Initialize>( Alloc::label( prop ), m_offset_map, m_sacado_size.value );
    }

  template< class AllocationProperties , typename iType >
  explicit inline
  View( const AllocationProperties & prop ,
        const iType * const n ,
        // Impl::ViewAllocProp::size_type exists when the traits and allocation properties
        // are valid for allocating viewed memory.
        const typename Impl::ViewAllocProp< traits , AllocationProperties >::size_type = 0 )
    : m_ptr_on_device(0)
    {
      typedef Impl::ViewAllocProp< traits , AllocationProperties > Alloc ;
      typedef typename traits::execution_space execution_space;
      typedef typename traits::memory_space memory_space;

      const size_t n0 = Rank >= 0 ? n[0] : 0 ;
      const size_t n1 = Rank >= 1 ? n[1] : 0 ;
      const size_t n2 = Rank >= 2 ? n[2] : 0 ;
      const size_t n3 = Rank >= 3 ? n[3] : 0 ;
      const size_t n4 = Rank >= 4 ? n[4] : 0 ;
      const size_t n5 = Rank >= 5 ? n[5] : 0 ;
      const size_t n6 = Rank >= 6 ? n[6] : 0 ;
      const size_t n7 = Rank >= 7 ? n[7] : 0 ;
      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
      m_storage_size =
        Impl::GetSacadoSize<unsigned(Rank)>::eval(n0,n1,n2,n3,n4,n5,n6,n7);
      if (m_storage_size == 0)
        m_storage_size = global_sacado_mp_vector_size;
      m_sacado_size = m_storage_size;
      m_ptr_on_device =
        m_allocation.template allocate<execution_space,memory_space,Alloc::Initialize>( Alloc::label( prop ), m_offset_map, m_sacado_size.value );
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
            Impl::is_same<T,typename traits::value_type>::value ||
            Impl::is_same<T,typename traits::non_const_value_type>::value
          ),
        const size_t >::type n7 = 0 )
    : m_ptr_on_device(ptr)
    {
      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
      m_stride = 1 ;
      m_storage_size = Impl::GetSacadoSize<unsigned(Rank)>::eval(n0,n1,n2,n3,n4,n5,n6,n7);
      if (m_storage_size == 0)
        m_storage_size = global_sacado_mp_vector_size;
      m_sacado_size = m_storage_size;
      m_allocation.assign(ptr);
      m_management.set_unmanaged();
    }

  //------------------------------------
  /** \brief  Constructors for subviews requires following
   *          type-compatibility condition, enforce via StaticAssert.
   *
   *  Impl::is_same< View ,
   *                 typename Impl::ViewSubview< View<D,A1,A2,A3,Impl::ViewDefault>
   *                                           , ArgType0 , ArgType1 , ArgType2 , ArgType3
   *                                           , ArgType4 , ArgType5 , ArgType6 , ArgType7
   *                 >::type >::value
   */
  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
          , class SubArg4_type , class SubArg5_type , class SubArg6_type , class SubArg7_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2 , const SubArg3_type & arg3
      , const SubArg4_type & arg4 , const SubArg5_type & arg5
      , const SubArg6_type & arg6 , const SubArg7_type & arg7
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
          , class SubArg4_type , class SubArg5_type , class SubArg6_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2 , const SubArg3_type & arg3
      , const SubArg4_type & arg4 , const SubArg5_type & arg5
      , const SubArg6_type & arg6
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
          , class SubArg4_type , class SubArg5_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2 , const SubArg3_type & arg3
      , const SubArg4_type & arg4 , const SubArg5_type & arg5
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
          , class SubArg4_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2 , const SubArg3_type & arg3
      , const SubArg4_type & arg4
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2 , const SubArg3_type & arg3
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type , class SubArg2_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      , const SubArg2_type & arg2
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type , class SubArg1_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0 , const SubArg1_type & arg1
      );

  template< class D , class A1 , class A2 , class A3
          , class SubArg0_type
          >
  KOKKOS_INLINE_FUNCTION
  View( const View<D,A1,A2,A3,Impl::ViewMPVectorContiguous> & src
      , const SubArg0_type & arg0
      );

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
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutRight, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_offset_map, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      // May have partitioned
      return m_ptr_on_device[ m_stride * i0 ];
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutRight, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutRight, 1, iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1) ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 2,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_offset_map, i0, i1, i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 3,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_offset_map, i0, i1, i2, i3 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 4,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_offset_map, i0, i1, i2, i3, i4 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
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
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_offset_map, i0, i1, i2, i3, i4, i5 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4,i5) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
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
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_offset_map, i0, i1, i2, i3, i4, i5, i6 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4,i5,i6) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutRight, 7,
                                      iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
        const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , int ) const
    { return operator()(i0,i1,i2,i3,i4,i5,i6); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 1:

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutLeft, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_offset_map, i0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      // May have partitioned
      return m_ptr_on_device[ m_stride * i0 ];
    }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutLeft, 1, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    { return operator()( i0 ); }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type, traits, LayoutLeft, 1, iType0 >::type
    at( const iType0 & i0 , int , int , int , int , int , int , int ) const
    { return operator()(i0); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 2:

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1) ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 2,
                                      iType0, iType1 >::type
    at( const iType0 & i0 , const iType1 & i1 , int , int , int , int , int , int ) const
    { return operator()(i0,i1); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 3:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_offset_map, i0, i1, i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 3,
                                      iType0, iType1, iType2 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , int , int , int , int , int ) const
    { return operator()(i0,i1,i2); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 4:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_offset_map, i0, i1, i2, i3 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 4,
                                      iType0, iType1, iType2, iType3 >::type
    at( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , int , int , int , int ) const
    { return operator()(i0,i1,i2,i3); }

  //------------------------------------
  //------------------------------------
  // Array operators, traits::rank 5:

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_offset_map, i0, i1, i2, i3, i4 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 5,
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
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_offset_map, i0, i1, i2, i3, i4, i5 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4,i5) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 6,
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
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_offset_map, i0, i1, i2, i3, i4, i5, i6 );

      return m_ptr_on_device[ m_stride * m_offset_map(i0,i1,i2,i3,i4,i5,i6) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5, typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< reference_type,
                                      traits, LayoutLeft, 7,
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
    { m_offset_map.stride( s ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type capacity() const
    { return m_stride * m_offset_map.cardinality(); }

  // Size of sacado dimension
  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::size_type sacado_size() const
    { return m_sacado_size.value; }

  Kokkos::Impl::AllocationTracker const & tracker() const
  { return m_allocation.m_tracker; }
};

/** \brief  A deep copy between views of the same specialization, compatible type,
 *          same rank, same layout are handled by that specialization.
 */
template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
inline
void deep_copy( const View<DT,DL,DD,DM,Impl::ViewMPVectorContiguous> & dst ,
                const View<ST,SL,SD,SM,Impl::ViewMPVectorContiguous> & src ,
                typename Impl::enable_if<(
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewMPVectorContiguous>::intrinsic_scalar_type ,
                                 typename View<ST,SL,SD,SM,Impl::ViewMPVectorContiguous>::intrinsic_scalar_type >::value
                  &&
                  Impl::is_same< typename View<DT,DL,DD,DM,Impl::ViewMPVectorContiguous>::array_layout ,
                                 typename View<ST,SL,SD,SM,Impl::ViewMPVectorContiguous>::array_layout >::value
                  &&
                  ( unsigned(View<DT,DL,DD,DM,Impl::ViewMPVectorContiguous>::rank) ==
                    unsigned(View<ST,SL,SD,SM,Impl::ViewMPVectorContiguous>::rank) )
                )>::type * = 0 )
{
  typedef View<DT,DL,DD,DM,Impl::ViewMPVectorContiguous> dst_type ;
  typedef View<ST,SL,SD,SM,Impl::ViewMPVectorContiguous> src_type ;

  typedef typename dst_type::array_type dst_array_type ;
  typedef typename src_type::array_type src_array_type ;

  dst_array_type dst_array = dst ;
  src_array_type src_array = src ;

  deep_copy( dst_array , src_array );
}

//----------------------------------------------------------------------------

/** \brief  Create compatible view for HostMirror device.
 *
 * Specialized for Sacado::MP::Vector with contiguous layout.
 */
template< class T , class L , class D , class M >
typename Impl::enable_if<(
    View<T,L,D,M,Impl::ViewMPVectorContiguous>::is_managed
  ), typename View<T,L,D,M,Impl::ViewMPVectorContiguous>::HostMirror >::type
inline
create_mirror( const View<T,L,D,M,Impl::ViewMPVectorContiguous> & src )
{
  typedef View<T,L,D,M,Impl::ViewMPVectorContiguous> view_type ;
  typedef typename view_type::HostMirror             host_view_type ;
  //typedef typename view_type::memory_space           memory_space ;

  // 'view' is managed therefore we can allocate a
  // compatible host_view through the ordinary constructor.

  std::string label = src.tracker().label();
  label.append("_mirror");

  unsigned dims[8];
  src.dimensions(dims);
  dims[view_type::Rank] = src.sacado_size();
  return host_view_type( label , dims );
}

template< class DstViewType ,
          class T , class L , class D , class M >
KOKKOS_INLINE_FUNCTION
DstViewType
partition( const View<T,L,D,M,Impl::ViewMPVectorContiguous> & src ,
           const unsigned beg ,
           const unsigned end )
{
  DstViewType dst ;
  const Sacado::MP::VectorPartition part( beg , end );
  Impl::ViewAssignment<typename DstViewType::specialize,
                       Impl::ViewMPVectorContiguous>( dst , src , part );

  return dst ;
}

} // namespace Kokkos

namespace Kokkos {
namespace Impl {

template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type , class SubArg7_type
        >
struct ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                  , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                  , SubArg4_type , SubArg5_type , SubArg6_type , SubArg7_type >
{
private:

  typedef View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >  SrcViewType ;

  enum { V0 = Impl::is_same< SubArg0_type , void >::value ? 1 : 0 };
  enum { V1 = Impl::is_same< SubArg1_type , void >::value ? 1 : 0 };
  enum { V2 = Impl::is_same< SubArg2_type , void >::value ? 1 : 0 };
  enum { V3 = Impl::is_same< SubArg3_type , void >::value ? 1 : 0 };
  enum { V4 = Impl::is_same< SubArg4_type , void >::value ? 1 : 0 };
  enum { V5 = Impl::is_same< SubArg5_type , void >::value ? 1 : 0 };
  enum { V6 = Impl::is_same< SubArg6_type , void >::value ? 1 : 0 };
  enum { V7 = Impl::is_same< SubArg7_type , void >::value ? 1 : 0 };

  // The source view rank must be equal to the input argument rank
  // Once a void argument is encountered all subsequent arguments must be void.
  enum { InputRank =
    Impl::StaticAssert<( SrcViewType::rank ==
                         ( V0 ? 0 : (
                           V1 ? 1 : (
                           V2 ? 2 : (
                           V3 ? 3 : (
                           V4 ? 4 : (
                           V5 ? 5 : (
                           V6 ? 6 : (
                           V7 ? 7 : 8 ))))))) ))
                       &&
                       ( SrcViewType::rank ==
                         ( 8 - ( V0 + V1 + V2 + V3 + V4 + V5 + V6 + V7 ) ) )
    >::value ? SrcViewType::rank : 0 };

  enum { R0 = Impl::ViewOffsetRange< SubArg0_type >::is_range ? 1 : 0 };
  enum { R1 = Impl::ViewOffsetRange< SubArg1_type >::is_range ? 1 : 0 };
  enum { R2 = Impl::ViewOffsetRange< SubArg2_type >::is_range ? 1 : 0 };
  enum { R3 = Impl::ViewOffsetRange< SubArg3_type >::is_range ? 1 : 0 };
  enum { R4 = Impl::ViewOffsetRange< SubArg4_type >::is_range ? 1 : 0 };
  enum { R5 = Impl::ViewOffsetRange< SubArg5_type >::is_range ? 1 : 0 };
  enum { R6 = Impl::ViewOffsetRange< SubArg6_type >::is_range ? 1 : 0 };
  enum { R7 = Impl::ViewOffsetRange< SubArg7_type >::is_range ? 1 : 0 };

  enum { OutputRank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
                    + unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7) };

  // Reverse
  enum { R0_rev = 0 == InputRank ? 0u : (
                  1 == InputRank ? unsigned(R0) : (
                  2 == InputRank ? unsigned(R1) : (
                  3 == InputRank ? unsigned(R2) : (
                  4 == InputRank ? unsigned(R3) : (
                  5 == InputRank ? unsigned(R4) : (
                  6 == InputRank ? unsigned(R5) : (
                  7 == InputRank ? unsigned(R6) : unsigned(R7) ))))))) };

  typedef typename SrcViewType::array_layout  SrcViewLayout ;

  // Choose array layout, attempting to preserve original layout if at all possible.
  typedef typename Impl::if_c<
     ( // Same Layout IF
       // OutputRank 0
       ( OutputRank == 0 )
       ||
       // OutputRank 1 or 2, InputLayout Left, Interval 0
       // because single stride one or second index has a stride.
       ( OutputRank <= 2 && R0 && Impl::is_same<SrcViewLayout,LayoutLeft>::value )
       ||
       // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
       // because single stride one or second index has a stride.
       ( OutputRank <= 2 && R0_rev && Impl::is_same<SrcViewLayout,LayoutRight>::value )
     ), SrcViewLayout , Kokkos::LayoutStride >::type OutputViewLayout ;

  // Choose data type as a purely dynamic rank array to accomodate a runtime range.
  typedef typename Impl::if_c< OutputRank == 0 , typename SrcViewType::value_type ,
          typename Impl::if_c< OutputRank == 1 , typename SrcViewType::value_type *,
          typename Impl::if_c< OutputRank == 2 , typename SrcViewType::value_type **,
          typename Impl::if_c< OutputRank == 3 , typename SrcViewType::value_type ***,
          typename Impl::if_c< OutputRank == 4 , typename SrcViewType::value_type ****,
          typename Impl::if_c< OutputRank == 5 , typename SrcViewType::value_type *****,
          typename Impl::if_c< OutputRank == 6 , typename SrcViewType::value_type ******,
          typename Impl::if_c< OutputRank == 7 , typename SrcViewType::value_type *******,
                                                 typename SrcViewType::value_type ********
  >::type >::type >::type >::type >::type >::type >::type >::type  OutputData ;

  // Choose space.
  // If the source view's template arg1 or arg2 is a space then use it,
  // otherwise use the source view's execution space.

  typedef typename Impl::if_c< Impl::is_space< SrcArg1Type >::value , SrcArg1Type ,
          typename Impl::if_c< Impl::is_space< SrcArg2Type >::value , SrcArg2Type , typename SrcViewType::device_type
  >::type >::type OutputSpace ;

public:

  // If keeping the layout then match non-data type arguments
  // else keep execution space and memory traits.
  typedef typename
    Impl::if_c< Impl::is_same< SrcViewLayout , OutputViewLayout >::value
              , Kokkos::View< OutputData , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
              , Kokkos::View< OutputData , OutputViewLayout , OutputSpace
                            , typename SrcViewType::memory_traits
                            , Impl::ViewMPVectorContiguous >
              >::type  type ;
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Analyze the array shape of a Sacado::MP::Vector.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::MP::Vector< StorageType > , ... >
 *  can be determined at compile-time.
 *
 *  This treats Sacado::MP::Vector as an atomic scalar.
 */
template< class StorageType >
struct AnalyzeShape< Sacado::MP::Vector< StorageType > >
  : Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > // Treat as a scalar
{
private:

  typedef AnalyzeShape< typename StorageType::value_type > nested ;

public:

  typedef ViewMPVectorContiguous specialize ;

  typedef Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > shape ;

  typedef       Sacado::MP::Vector< StorageType >  array_intrinsic_type ;
  typedef const Sacado::MP::Vector< StorageType >  const_array_intrinsic_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_array_intrinsic_type ;

  typedef       Sacado::MP::Vector< StorageType >  type ;
  typedef const Sacado::MP::Vector< StorageType >  const_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_type ;

  typedef       Sacado::MP::Vector< StorageType >  value_type ;
  typedef const Sacado::MP::Vector< StorageType >  const_value_type ;
  typedef       Sacado::MP::Vector< StorageType >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::MP::Vector.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::MP::Vector< StorageType > , ... >
 *  can be determined at compile-time.
 *
 *  This treats Sacado::MP::Vector as an array.
 */
template< class StorageType, class Layout >
struct AnalyzeSacadoShape< Sacado::MP::Vector< StorageType >, Layout >
  : Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > // Treat as a scalar
{
private:

  typedef AnalyzeSacadoShape< typename StorageType::value_type, Layout > nested ;

public:

  typedef ViewMPVectorContiguous specialize ;

  typedef Shape< sizeof(Sacado::MP::Vector< StorageType >) , 0 > shape ;

  // If ( ! StorageType::is_static ) then 0 == StorageType::static_size and the first array declaration is not used.
  // However, the compiler will still generate this type declaration and it must not have a zero length.
  //
  // For LayoutLeft, always use a dynamic dimension, which get's tacked on to the left of the array type
  // I think this means our approach doesn't really work for Sacado::MP::Vector<StaticFixedStorage<...> >[N] ???
  typedef typename
    if_c< StorageType::is_static && is_same<Layout, LayoutRight>::value
        , typename nested::array_intrinsic_type [ StorageType::is_static ? StorageType::static_size : 1 ]
        , typename nested::array_intrinsic_type *
        >::type array_intrinsic_type ;

  typedef typename
    if_c< StorageType::is_static && is_same<Layout, LayoutRight>::value
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

  // Replace MP::Vector<S> with S::value_type in array specification to
  // form the flattened array with the same rank
  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewMPVectorContiguous , ViewMPVectorContiguous , void >
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
                    ViewAssignable< View<DT,DL,DD,DM,specialize> ,
                                    View<ST,SL,SD,SM,specialize> >::value
                    )>::type * = 0
                  )
  {
    dst.m_offset_map.assign( src.m_offset_map );
    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;
    dst.m_allocation    = src.m_allocation ;
    dst.m_storage_size  = src.m_storage_size ;
    dst.m_sacado_size   = src.m_sacado_size;
    dst.m_management    = src.m_management ;
  }

  //------------------------------------
  /** \brief  Partition of compatible value and shape
   *
   * Currently only allowed for static storage types as partitioning is
   * quite expensive for dynamic
   */

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
                    ( unsigned(View<DT,DL,DD,DM,specialize>::Rank) ==
                      unsigned(View<ST,SL,SD,SM,specialize>::Rank) )
                    &&
                    // Destination is not managed
                    ! View<DT,DL,DD,DM,specialize>::is_managed
                    &&
                    // Views have static storage
                    ( View<DT,DL,DD,DM,specialize>::is_static &&
                      View<ST,SL,SD,SM,specialize>::is_static )
                  ), const Sacado::MP::VectorPartition & >::type part )
  {
    typedef View<DT,DL,DD,DM,specialize>   dst_type ;

    // Must have: begin = i * src.m_sacado_size.value
    //            end   = begin + src.m_sacado_size.value
    //
    if ( dst_type::is_static && (part.begin % dst.m_sacado_size.value ||
                                 part.begin + dst.m_sacado_size.value != part.end) ) {
      Kokkos::abort("Kokkos::View< Sacado::MP::Vector ... > incompatible partitioning");
    }

    const int length = part.end - part.begin ;

    dst.m_offset_map.assign( src.m_offset_map );

    // Original Sacado::MP::Vector length
    dst.m_storage_size = src.m_storage_size ;
    dst.m_sacado_size = length;

    dst.m_stride = src.m_stride * src.m_sacado_size.value
                                / dst.m_sacado_size.value ;

    dst.m_ptr_on_device =
      reinterpret_cast<typename dst_type::value_type*>( src.m_ptr_on_device ) +
      part.begin / dst.m_sacado_size.value ;
    dst.m_allocation.m_scalar_ptr_on_device =
      src.m_allocation.m_scalar_ptr_on_device +
      (part.begin / dst.m_sacado_size.value) * src.m_storage_size ;
    dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;

    dst.m_management      = src.m_management ;
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from range of Rank-1 array, either layout */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                          ViewAssignable< View<DT,DL,DD,DM,specialize> ,
                                          View<ST,SL,SD,SM,specialize> >::value
                    &&
                          ( View<ST,SL,SD,SM,specialize>::rank == 1 )
                    &&
                          ( View<DT,DL,DD,DM,specialize>::rank == 1 )
                    &&
                          ( View<DT,DL,DD,DM,specialize>::rank_dynamic == 1 )
                  ) >::type * = 0 )
  {
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_offset_map , 1 , range.first );
      assert_shape_bounds( src.m_offset_map , 1 , range.second - 1 );

      dst.m_offset_map.assign( range.second - range.first , 0,0,0,0,0,0,0 );
      dst.m_ptr_on_device = src.m_ptr_on_device + range.first ;
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + range.first * src.m_storage_size ;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
      dst.m_stride       = src.m_stride ;
      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
      dst.m_management     = src.m_management ;
    }
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
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
    dst.m_offset_map.assign( src.m_offset_map.N0 , 0,0,0,0,0,0,0);
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_offset_map.N0 * i1 ;
    dst.m_allocation.m_scalar_ptr_on_device =
      src.m_allocation.m_scalar_ptr_on_device + src.m_offset_map.N0 * i1 * src.m_storage_size ;
    dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
    dst.m_stride        = src.m_stride ;
    dst.m_storage_size  = src.m_storage_size ;
    dst.m_sacado_size   = src.m_sacado_size;
    dst.m_management      = src.m_management ;
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const std::pair<iType,iType> & range ,
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
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_stride        = 1 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_offset_map , 2 , range.first , i1 );
      assert_shape_bounds( src.m_offset_map , 2 , range.second - 1 , i1 );

      dst.m_management      = src.m_management ;
      dst.m_offset_map.N0 = range.second - range.first ;
      dst.m_ptr_on_device =
        src.m_ptr_on_device + src.m_offset_map(range.first,i1);
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + src.m_offset_map(range.first,i1) * src.m_storage_size ;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
       dst.m_stride      = src.m_stride ;
      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
    }
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), unsigned >::type i1 )
  {
    dst.m_offset_map.assign( src.m_offset_map.N0, 1, 0,0,0,0,0,0);
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_offset_map.N0 * i1 ;
    dst.m_allocation.m_scalar_ptr_on_device =
      src.m_allocation.m_scalar_ptr_on_device + src.m_offset_map.N0 * i1 * src.m_storage_size;
    dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
    dst.m_stride        = src.m_stride ;
    dst.m_storage_size  = src.m_storage_size ;
    dst.m_sacado_size   = src.m_sacado_size;
    dst.m_management      = src.m_management ;
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), unsigned >::type i1 )
  {
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_stride        = 1 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_offset_map , 2 , range.first , i1 );
      assert_shape_bounds( src.m_offset_map , 2 , range.second - 1 , i1 );

      dst.m_management      = src.m_management ;
      dst.m_offset_map.N0 = range.second - range.first ;
      dst.m_offset_map.N1 = 1 ;
      dst.m_offset_map.S0 = range.second - range.first ;
      dst.m_ptr_on_device =
        src.m_ptr_on_device + src.m_offset_map(range.first,i1);
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + src.m_offset_map(range.first,i1) * src.m_storage_size ;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
       dst.m_stride      = src.m_stride ;
      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  ALL ,
                  const std::pair<iType,iType> & range1 ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_stride        = 1 ;
    dst.m_ptr_on_device = 0 ;

    if ( range1.first < range1.second ) {
      assert_shape_bounds( src.m_offset_map , 2 , 0 , range1.first );
      assert_shape_bounds( src.m_offset_map , 2 , src.m_offset_map.N0 - 1 , range1.second - 1 );

      dst.m_offset_map.assign( src.m_offset_map.N0 ,
                               range1.second - range1.first ,
                               0,0,0,0,0,0 );
      dst.m_stride   = src.m_stride ;

      // operator: dst.m_ptr_on_device[ i0 + dst.m_stride * i1 ]
      dst.m_ptr_on_device = src.m_ptr_on_device + dst.m_offset_map.N0 * range1.first ;
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + dst.m_offset_map.N0 * range1.first * src.m_storage_size ;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
      dst.m_management     = src.m_management ;

      // LayoutRight won't work with how we are currently using the stride
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const std::pair<iType,iType> & range0 ,
                  ALL ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_stride        = 1 ;
    dst.m_ptr_on_device = 0 ;

    if ( range0.first < range0.second ) {
      assert_shape_bounds( src.m_offset_map , 2 , range0.first , 0 );
      assert_shape_bounds( src.m_offset_map , 2 , range0.second - 1 , src.m_offset_map.N1 - 1 );

      dst.m_offset_map.assign( range0.second - range0.first ,
                               src.m_offset_map.N1 ,
                               0,0,0,0,0,0 );
      dst.m_stride   = src.m_stride ;

      // operator: dst.m_ptr_on_device[ i0 + dst.m_stride * i1 ]
      dst.m_ptr_on_device = src.m_ptr_on_device + range0.first ;
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + range0.first * src.m_storage_size;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;
      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
      dst.m_management     = src.m_management ;

      // LayoutRight won't work with how we are currently using the stride
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,specialize> & dst ,
                  const View<ST,SL,SD,SM,specialize> & src ,
                  const std::pair<iType0,iType0> & range0 ,
                  const std::pair<iType1,iType1> & range1 ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_offset_map.assign(0,0,0,0,0,0,0,0);
    dst.m_stride        = 1 ;
    dst.m_ptr_on_device = 0 ;

    if ( (range0.first < range0.second && range1.first < range1.second) ) {
      assert_shape_bounds( src.m_offset_map , 2 , range0.first , range1.first );
      assert_shape_bounds( src.m_offset_map , 2 , range0.second - 1 , range1.second - 1 );

      dst.m_offset_map.assign( src.m_offset_map );
      dst.m_offset_map.N0 = range0.second - range0.first ;
      dst.m_offset_map.N1 = range1.second - range1.first ;
      dst.m_stride   = src.m_stride ;

      dst.m_ptr_on_device = src.m_ptr_on_device + src.m_offset_map(range0.first,range1.first);
      dst.m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device + src.m_offset_map(range0.first,range1.first) * src.m_storage_size;
      dst.m_allocation.m_tracker = src.m_allocation.m_tracker ;

      dst.m_storage_size = src.m_storage_size ;
      dst.m_sacado_size  = src.m_sacado_size;
      dst.m_management     = src.m_management ;
    }
  }
};

template<>
struct ViewAssignment< ViewDefault , ViewMPVectorContiguous , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,ViewMPVectorContiguous>::array_type & dst
                , const    View<ST,SL,SD,SM,ViewMPVectorContiguous> & src )
  {
    typedef View<ST,SL,SD,SM,ViewMPVectorContiguous> src_type ;

    if ( src.m_stride != 1 ) {
      Kokkos::abort("Kokkos::View< Sacado::MP::Vector ... > incompatible assignment");
    }

    unsigned dims[8];
    dims[0] = src.m_offset_map.N0;
    dims[1] = src.m_offset_map.N1;
    dims[2] = src.m_offset_map.N2;
    dims[3] = src.m_offset_map.N3;
    dims[4] = src.m_offset_map.N4;
    dims[5] = src.m_offset_map.N5;
    dims[6] = src.m_offset_map.N6;
    dims[7] = src.m_offset_map.N7;
    unsigned rank = src_type::Rank;
    unsigned sacado_size = src.sacado_size();
    if (is_same<typename src_type::array_layout, LayoutLeft>::value) {
      // Move sacado_size to the first dimension, shift all others up one
      for (unsigned i=rank; i>0; --i)
        dims[i] = dims[i-1];
      dims[0] = sacado_size;
    }
    else {
      dims[rank] = sacado_size;
    }
    dst.m_offset_map.assign( dims[0] , dims[1] , dims[2] , dims[3] ,
                             dims[4] , dims[5] , dims[6] , dims[7] );

    dst.m_ptr_on_device = src.m_allocation.m_scalar_ptr_on_device;

    dst.m_tracker = src.m_allocation.m_tracker ;

    dst.m_management      = src.m_management ;
  }

  //------------------------------------
  /**
   * \brief Assign to flattened view where Sacado dimension is combined with
   * most adjacent dimension.  Must have same instrinsic value_type, layout,
   * and rank (add 1 to rank for sacado dimension, remove 1 for flattening).
   *
   * Would like to just use anything that is assignable to the flat_array_type,
   * e.g.,
   *
   * typename enable_if<
   *    (
   *      Impl::ViewAssignable<
   *        View<DT,DL,DD,DM,ViewDefault>,
   *        typename View<ST,SL,SD,SM,ViewMPVectorContiguous>::flat_array_type
   *      >::value
   *    ) >::type * = 0)
   *
   * except this conflicts with the overload above for array_type (since
   * ViewAssignable is loose on the ranks and array_type is actually assignable
   * to flat_array_type).  And we can't use flat_array_type as there are
   * use cases where the view is the same as flat_array_type but with different
   * memory traits.
   */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(
          View<DT,DL,DD,DM,ViewDefault>& dst,
    const View<ST,SL,SD,SM,ViewMPVectorContiguous>& src,
    typename enable_if<
      (
        ( is_same< typename View<DT,DL,DD,DM,ViewDefault>::value_type ,
                   typename View<ST,SL,SD,SM,ViewMPVectorContiguous>::intrinsic_scalar_type >::value ||
          is_same< typename View<DT,DL,DD,DM,ViewDefault>::non_const_value_type ,
                   typename View<ST,SL,SD,SM,ViewMPVectorContiguous>::intrinsic_scalar_type >::value
         ) &&
        is_same< typename View<DT,DL,DD,DM,ViewDefault>::array_layout ,
                 typename View<ST,SL,SD,SM,ViewMPVectorContiguous>::array_layout >::value &&
        ( unsigned(View<DT,DL,DD,DM,ViewDefault>::rank) ==
          unsigned(View<ST,SL,SD,SM,ViewMPVectorContiguous>::rank) )
      ) >::type * = 0)
  {
    typedef View<DT,DL,DD,DM,ViewDefault>   dst_type ;
    typedef typename dst_type::array_layout dst_layout_type ;

    if ( src.m_stride != 1 ) {
      Kokkos::abort("Kokkos::View< Sacado::MP::Vector ... > incompatible assignment");
    }

    // Create flattened shape
    unsigned dims[8];
    dims[0] = src.m_offset_map.N0;
    dims[1] = src.m_offset_map.N1;
    dims[2] = src.m_offset_map.N2;
    dims[3] = src.m_offset_map.N3;
    dims[4] = src.m_offset_map.N4;
    dims[5] = src.m_offset_map.N5;
    dims[6] = src.m_offset_map.N6;
    dims[7] = src.m_offset_map.N7;
    unsigned rank = dst_type::Rank;
    unsigned sacado_size = src.sacado_size();
    if (is_same<dst_layout_type, LayoutLeft>::value) {
      dims[0] = dims[0]*sacado_size;
      dims[rank] = 0;
    }
    else {
      dims[rank-1] = dims[rank-1]*sacado_size;
      dims[rank] = 0;
    }
    dst.m_offset_map.assign( dims[0] , dims[1] , dims[2] , dims[3] ,
                             dims[4] , dims[5] , dims[6] , dims[7] );

    dst.m_ptr_on_device = src.m_allocation.m_scalar_ptr_on_device;

    dst.m_tracker = src.m_allocation.m_tracker ;

    dst.m_management      = src.m_management ;
  }
};

// Specialization for deep_copy( view, view::value_type ) for Cuda
#if defined( KOKKOS_HAVE_CUDA )
template< class OutputView , unsigned Rank >
struct ViewFill< OutputView , Rank ,
                 typename enable_if< is_same< typename OutputView::specialize,
                                              ViewMPVectorContiguous >::value &&
                                     is_same< typename OutputView::execution_space,
                                              Cuda >::value >::type >
{
  typedef typename OutputView::const_value_type   const_value_type ;
  typedef typename OutputView::execution_space    execution_space ;
  typedef typename OutputView::size_type          size_type ;

  template <unsigned VectorLength>
  struct Kernel {
    typedef typename OutputView::execution_space execution_space ;
    const OutputView output;
    const_value_type input;

    Kernel( const OutputView & arg_out , const_value_type & arg_in ) :
      output(arg_out), input(arg_in) {}

    typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;

    KOKKOS_INLINE_FUNCTION
    void operator()( const team_member & dev ) const
    {
      const size_type tidx = dev.team_rank() % VectorLength;
      const size_type tidy = dev.team_rank() / VectorLength;
      const size_type nrow = dev.team_size() / VectorLength;
      const size_type nvec = output.sacado_size();

      const size_type i0 = dev.league_rank() * nrow + tidy;
      if ( i0 >= output.dimension_0() ) return;

      for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
      for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
      for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
      for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
      for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
      for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
      for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      for ( size_type is = tidx ; is < nvec ; is+=VectorLength ) {
        output.at(i0,i1,i2,i3,i4,i5,i6,i7).fastAccessCoeff(is) =
          input.fastAccessCoeff(is) ;
      }}}}}}}}
    }
  };

  ViewFill( const OutputView & output , const_value_type & input )
  {
    if ( Sacado::is_constant(input) ) {
      deep_copy( output , input.fastAccessCoeff(0) );
    }
    else {

      // Coalesced accesses are 128 bytes in size
      typedef typename OutputView::intrinsic_scalar_type scalar_type;
      const unsigned vector_length =
        ( 128 + sizeof(scalar_type)-1 ) / sizeof(scalar_type);

      // 8 warps per block should give good occupancy
      const size_type block_size = 256;

      const size_type rows_per_block = block_size / vector_length;
      const size_type n = output.dimension_0();
      const size_type league_size = ( n + rows_per_block-1 ) / rows_per_block;
      const size_type team_size = rows_per_block * vector_length;
      Kokkos::TeamPolicy< execution_space > config( league_size, team_size );

      parallel_for( config, Kernel<vector_length>(output, input) );
      execution_space::fence();
    }
  }

};
#endif /* #if defined( KOKKOS_HAVE_CUDA ) */

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

// Construct subview of a Rank 8 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type , class SubArg7_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    , const SubArg6_type & arg6
    , const SubArg7_type & arg7
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , SubArg6_type , SubArg7_type >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;
    typedef Impl::ViewOffsetRange< SubArg6_type > R6 ;
    typedef Impl::ViewOffsetRange< SubArg7_type > R7 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , R6::dimension( src.m_offset_map.N6 , arg6 )
                                 , R7::dimension( src.m_offset_map.N7 , arg7 )
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , R6::begin( arg6 )
                                        , R7::begin( arg7 ) );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , R3::begin( arg3 )
                        , R4::begin( arg4 )
                        , R5::begin( arg5 )
                        , R6::begin( arg6 )
                        , R7::begin( arg7 ) ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 7 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    , const SubArg6_type & arg6
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , SubArg6_type , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;
    typedef Impl::ViewOffsetRange< SubArg6_type > R6 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , R6::dimension( src.m_offset_map.N6 , arg6 )
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , R6::begin( arg6 )
                                        , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , R3::begin( arg3 )
                        , R4::begin( arg4 )
                        , R5::begin( arg5 )
                        , R6::begin( arg6 )
                        , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 6 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , 0
                                        , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , R3::begin( arg3 )
                        , R4::begin( arg4 )
                        , R5::begin( arg5 )
                        , 0
                        , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 5 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , 0
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , 0
                                        , 0
                                        , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , R3::begin( arg3 )
                        , R4::begin( arg4 )
                        , 0
                        , 0
                        , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 4 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , 0
                                 , 0
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , 0
                                        , 0
                                        , 0
                                        , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , R3::begin( arg3 )
                        , 0
                        , 0
                        , 0
                        , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 3 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , SubArg2_type , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , 0 , 0 , 0 , 0 , 0);

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , 0 , 0 , 0 , 0 , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , R1::begin( arg1 )
                        , R2::begin( arg2 )
                        , 0 , 0 , 0 , 0 , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 2 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , SubArg1_type , void , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , 0 , 0 , 0 , 0 , 0 , 0 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , 0 , 0 , 0 , 0 , 0 , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 ) , R1::begin( arg1 )
                          , 0 , 0 , 0 , 0 , 0 , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

// Construct subview of a Rank 1 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewMPVectorContiguous >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous > & src
    , const SubArg0_type & arg0
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_allocation()
  , m_offset_map()
  , m_stride(1)
  , m_storage_size(0)
  , m_sacado_size(0)
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewMPVectorContiguous >
                           , SubArg0_type , void , void , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , 0 , 0 , 0 , 0 , 0 , 0 , 0 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , 0 , 0 , 0 , 0 , 0 , 0 , 0 );
      m_allocation = src.m_allocation ;
      m_allocation.m_scalar_ptr_on_device =
        src.m_allocation.m_scalar_ptr_on_device +
        src.m_offset_map( R0::begin( arg0 )
                        , 0 , 0 , 0 , 0 , 0 , 0 , 0 ) *
        src.m_storage_size;
      m_sacado_size = src.m_sacado_size;
      m_storage_size = src.m_storage_size;
      m_stride = src.m_stride;
    }
  }
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEW_MP_VECTOR_CONTIGUOUS_HPP */
