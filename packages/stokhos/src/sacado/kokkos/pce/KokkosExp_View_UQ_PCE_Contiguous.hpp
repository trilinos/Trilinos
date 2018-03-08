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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_UQ_PCE_CONTIGUOUS_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_UQ_PCE_CONTIGUOUS_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_UQ_PCE.hpp"
#include "Sacado_UQ_PCE_Traits.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_AnalyzeStokhosShape.hpp"
#include "Kokkos_View_Utils.hpp"
#include "Kokkos_View_UQ_PCE_Utils.hpp"

#include "impl/KokkosExp_ViewMapping.hpp"

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ViewPCEContiguous {};

template< class ... Args >
struct is_ViewPCEContiguous { enum { value = false }; };

template< class D , class ... P , class ... Args >
struct is_ViewPCEContiguous< Kokkos::View<D,P...> , Args... > {
  enum { value =
    std::is_same< typename Kokkos::ViewTraits<D,P...>::specialize
                , ViewPCEContiguous >::value
    &&
    ( ( sizeof...(Args) == 0 ) ||
      is_ViewPCEContiguous< Args... >::value ) };
};

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

namespace Kokkos {

// Overload of deep_copy for UQ::PCE views intializing to a constant scalar
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  typedef View<DT,DP...> view_type;
  typedef typename view_type::array_type::value_type scalar_type;
  typedef typename FlatArrayType<view_type>::type flat_array_type;
  if (value == scalar_type(0))
    Kokkos::Impl::StokhosViewFill< flat_array_type >( view , value );
  else
    Kokkos::Impl::StokhosViewFill< view_type>( view , value );
}

// Overload of deep_copy for UQ::PCE views intializing to a constant UQ::PCE
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  Kokkos::Impl::StokhosViewFill< View<DT,DP...> >( view , value );
}

// Overload of deep_copy for UQ::PCE views intializing to a constant scalar
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  Kokkos::Impl::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  typedef View<DT,DP...> view_type;
  typedef typename view_type::array_type::value_type scalar_type;
  typedef typename FlatArrayType<view_type>::type flat_array_type;
  if (value == scalar_type(0))
    Kokkos::Impl::StokhosViewFill< flat_array_type >( view , value );
  else
    Kokkos::Impl::StokhosViewFill< view_type>( view , value );
}

// Overload of deep_copy for UQ::PCE views intializing to a constant UQ::PCE
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  Kokkos::Impl::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  Kokkos::Impl::StokhosViewFill< View<DT,DP...> >( view , value );
}

namespace Experimental {
namespace Impl {

// Deep copy between views not assuming contiguous storage of arrays
// Need to use team interface for Cuda
template< class OutputView , class InputView >
struct DeepCopyNonContiguous
{
  typedef typename OutputView::execution_space execution_space ;
  typedef typename execution_space::size_type  size_type ;

  const OutputView output ;
  const InputView  input ;

  DeepCopyNonContiguous( const OutputView & arg_out ,
                         const InputView & arg_in ) :
    output( arg_out ), input( arg_in )
  {
    parallel_for( output.dimension_0() , *this );
    execution_space::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7) ;
    }}}}}}}
  }
};

} // namespace Impl
} // namespace Experimental

/* Specialize for deep copy of UQ::PCE */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Deep copy destination must be non-const" );

  static_assert(
    ( unsigned(ViewTraits<DT,DP...>::rank) ==
      unsigned(ViewTraits<ST,SP...>::rank) )
    , "Deep copy destination and source must have same rank" );

  typedef View<DT,DP...> dst_type ;
  typedef View<ST,SP...> src_type ;
  typedef typename dst_type::array_type dst_array_type ;
  typedef typename src_type::array_type src_array_type ;

  if ( is_allocation_contiguous(dst) && is_allocation_contiguous(src) ) {
    dst_array_type dst_array = dst ;
    src_array_type src_array = src ;
    deep_copy( dst_array , src_array );
  }

  // otherwise, use a custom kernel
  else {

    // If views are in the same memory space, copy component-wise
    if ( Impl::is_same< typename dst_type::memory_space ,
                        typename src_type::memory_space >::value ) {
      Experimental::Impl::DeepCopyNonContiguous< dst_type , src_type >( dst , src );
    }

    else {

      typedef View< typename src_type::non_const_data_type ,
                    typename src_type::array_layout ,
                    typename src_type::execution_space > tmp_src_type;
      typedef typename tmp_src_type::array_type tmp_src_array_type;
      typedef View< typename dst_type::non_const_data_type ,
                    typename dst_type::array_layout ,
                    typename dst_type::execution_space > tmp_dst_type;
      typedef typename tmp_dst_type::array_type tmp_dst_array_type;

      // Copy src into a contiguous view in src's memory space,
      // then copy to dst
      if (  is_allocation_contiguous(dst) &&
           !is_allocation_contiguous(src) ) {
        size_t src_dims[8];
        //src.dimensions(src_dims);
        src_dims[0] = src.dimension_0();
        src_dims[1] = src.dimension_1();
        src_dims[2] = src.dimension_2();
        src_dims[3] = src.dimension_3();
        src_dims[4] = src.dimension_4();
        src_dims[5] = src.dimension_5();
        src_dims[6] = src.dimension_6();
        src_dims[7] = src.dimension_7();
        src_dims[src_type::Rank] = dimension_scalar(src);
        tmp_src_type src_tmp(
          view_alloc("src_tmp" , WithoutInitializing, cijk(src) ) ,
          src_dims[0], src_dims[1], src_dims[2], src_dims[3],
          src_dims[4], src_dims[5], src_dims[6], src_dims[7] );
        Experimental::Impl::DeepCopyNonContiguous< tmp_src_type , src_type >( src_tmp , src );
        dst_array_type dst_array = dst ;
        tmp_src_array_type src_array = src_tmp ;
        deep_copy( dst_array , src_array );
      }

      // Copy src into a contiguous view in dst's memory space,
      // then copy to dst
      else if ( !is_allocation_contiguous(dst) &&
                 is_allocation_contiguous(src) ) {
        size_t dst_dims[8];
        //dst.dimensions(dst_dims);
        dst_dims[0] = dst.dimension_0();
        dst_dims[1] = dst.dimension_1();
        dst_dims[2] = dst.dimension_2();
        dst_dims[3] = dst.dimension_3();
        dst_dims[4] = dst.dimension_4();
        dst_dims[5] = dst.dimension_5();
        dst_dims[6] = dst.dimension_6();
        dst_dims[7] = dst.dimension_7();
        dst_dims[dst_type::Rank] = dimension_scalar(dst);
        tmp_dst_type dst_tmp(
          view_alloc("dst_tmp" , WithoutInitializing, cijk(dst) ) ,
          dst_dims[0], dst_dims[1], dst_dims[2], dst_dims[3],
          dst_dims[4], dst_dims[5], dst_dims[6], dst_dims[7] );
        tmp_dst_array_type dst_array = dst_tmp ;
        src_array_type src_array = src ;
        deep_copy( dst_array , src_array );
        Experimental::Impl::DeepCopyNonContiguous< dst_type , tmp_dst_type >( dst , dst_tmp );
      }

      // Copy src into a contiguous view in src's memory space,
      // copy to a continugous view in dst's memory space, then copy to dst
      else {
        size_t src_dims[8];
        //src.dimensions(src_dims);
        src_dims[0] = src.dimension_0();
        src_dims[1] = src.dimension_1();
        src_dims[2] = src.dimension_2();
        src_dims[3] = src.dimension_3();
        src_dims[4] = src.dimension_4();
        src_dims[5] = src.dimension_5();
        src_dims[6] = src.dimension_6();
        src_dims[7] = src.dimension_7();
        src_dims[src_type::Rank] = dimension_scalar(src);
        tmp_src_type src_tmp(
          view_alloc("src_tmp" , WithoutInitializing, cijk(src) ) ,
          src_dims[0], src_dims[1], src_dims[2], src_dims[3],
          src_dims[4], src_dims[5], src_dims[6], src_dims[7] );
        Experimental::Impl::DeepCopyNonContiguous< tmp_src_type , src_type >( src_tmp , src );
        size_t dst_dims[8];
        //dst.dimensions(dst_dims);
        dst_dims[0] = dst.dimension_0();
        dst_dims[1] = dst.dimension_1();
        dst_dims[2] = dst.dimension_2();
        dst_dims[3] = dst.dimension_3();
        dst_dims[4] = dst.dimension_4();
        dst_dims[5] = dst.dimension_5();
        dst_dims[6] = dst.dimension_6();
        dst_dims[7] = dst.dimension_7();
        dst_dims[dst_type::Rank] = dimension_scalar(dst);
        tmp_dst_type dst_tmp(
          view_alloc("dst_tmp" , WithoutInitializing, cijk(dst) ) ,
          dst_dims[0], dst_dims[1], dst_dims[2], dst_dims[3],
          dst_dims[4], dst_dims[5], dst_dims[6], dst_dims[7] );
        tmp_dst_array_type dst_array = dst_tmp ;
        tmp_src_array_type src_array = src_tmp ;
        deep_copy( dst_array , src_array );
        Experimental::Impl::DeepCopyNonContiguous< dst_type , tmp_dst_type >( dst , dst_tmp );
      }
    }
  }
}

template <typename T, typename ... P>
struct is_view_uq_pce< View<T,P...> > {
  typedef View<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Experimental::Impl::ViewPCEContiguous >::value;
};

template <typename D, typename ... P>
struct FlatArrayType< View<D,P...>,
                      typename std::enable_if< is_view_uq_pce< View<D,P...> >::value >::type > {
  typedef View<D,P...> view_type;
  typedef typename view_type::traits::dimension dimension;
  typedef typename view_type::array_type::value_type flat_value_type;
  typedef typename Kokkos::Impl::ViewDataType< flat_value_type , dimension >::type flat_data_type;
  typedef View<flat_data_type,P...> type;
};

template <typename ViewType>
struct CijkType< ViewType,
                 typename std::enable_if< is_view_uq_pce< ViewType >::value >::type > {
  typedef typename ViewType::non_const_value_type::cijk_type type;
};

template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_view_uq_pce< View<T,P...> >::value, unsigned >::type
dimension_scalar(const View<T,P...>& view) {
  return view.implementation_map().dimension_scalar();
}

template <typename view_type>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_view_uq_pce<view_type>::value,
                typename CijkType<view_type>::type >::type
cijk(const view_type& view) {
  return view.implementation_map().cijk();
}

template <typename view_type>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_view_uq_pce<view_type>::value, bool >::type
is_allocation_contiguous(const view_type& view) {
  return view.implementation_map().is_allocation_contiguous();
}

template <typename ViewType>
ViewType
make_view(const std::string& label,
          const typename CijkType<ViewType>::type& cijk,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(view_alloc(label,cijk),
                  N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const std::string& label,
          const Impl::WithoutInitializing_t& init,
          const typename CijkType<ViewType>::type& cijk,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(view_alloc(label,init,cijk),
                  N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const ViewAllocateWithoutInitializing& init,
          const typename CijkType<ViewType>::type& cijk,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(view_alloc(init.label,
                                           WithoutInitializing,
                                           cijk),
                  N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
typename std::enable_if< is_view_uq_pce<ViewType>::value, ViewType>::type
make_view(typename ViewType::pointer_type ptr,
          const typename CijkType<ViewType>::type& cijk,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  size_t N[8] = { N0, N1, N2, N3, N4, N5, N6, N7 };
  N[ViewType::rank] = cijk.dimension();
  ViewType v(view_wrap(ptr, cijk),
             N[0], N[1], N[2], N[3], N[4], N[5], N[6], N[7]);
  return v;
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
//namespace Experimental {
namespace Impl {

// Allow passing of Cijk tensor through ViewCtorProp
template< typename Value, typename Execution, typename Memory >
struct ViewCtorProp< void , Stokhos::CrsProductTensor<Value, Execution, Memory> >
{
  ViewCtorProp() = default ;
  ViewCtorProp( const ViewCtorProp & ) = default ;
  ViewCtorProp & operator = ( const ViewCtorProp & ) = default ;

  typedef Stokhos::CrsProductTensor<Value, Execution, Memory> type ;

  ViewCtorProp( const type & arg ) : value( arg ) {}
  ViewCtorProp( type && arg ) : value( arg ) {}

  type value ;
};

template <typename AllocProp>
struct ctor_prop_has_cijk
{
  static const bool value = false;
};

template< typename T >
struct ctor_prop_has_cijk< ViewCtorProp<T> >
{
  static const bool value = false;
};

template< typename Value, typename Execution, typename Memory >
struct ctor_prop_has_cijk<
  ViewCtorProp< Stokhos::CrsProductTensor<Value, Execution, Memory> >
  >
{
  static const bool value = true;
};

template< typename T, typename ... P >
struct ctor_prop_has_cijk< ViewCtorProp<T,P...> >
{
  static const bool value =
    ctor_prop_has_cijk< ViewCtorProp<T> >::value ||
    ctor_prop_has_cijk< ViewCtorProp<P...> >::value;
};

} /* namespace Impl */
//} /* namespace Experimental */

template <typename CijkType, typename AllocProp>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< !Impl::ctor_prop_has_cijk<AllocProp>::value,
                         CijkType >::type
extract_cijk(const AllocProp& prop)
{
  return CijkType();
}

template <typename CijkType, typename AllocProp>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< Impl::ctor_prop_has_cijk<AllocProp>::value,
                         CijkType >::type
extract_cijk(const AllocProp& prop)
{
  return ( (const Impl::ViewCtorProp<void,CijkType>&) prop ).value;
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DataType , class ArrayLayout , typename StorageType >
struct ViewDataAnalysis< DataType     /* Original view data type */
                         , ArrayLayout
                         , Sacado::UQ::PCE< StorageType > >
{
private:

  typedef typename StorageType::value_type ScalarType;
  typedef ViewArrayAnalysis< DataType > array_analysis ;

public:

  // Specialized view data mapping:
  typedef Kokkos::Experimental::Impl::ViewPCEContiguous specialize ;

  typedef typename array_analysis::dimension             dimension ;
  typedef typename array_analysis::value_type            value_type ;
  typedef typename array_analysis::const_value_type      const_value_type ;
  typedef typename array_analysis::non_const_value_type  non_const_value_type ;

  // Generate analogous multidimensional array specification type.
  typedef typename
    ViewDataType< value_type , dimension >::type  type ;
  typedef typename
    ViewDataType< const_value_type , dimension >::type  const_type ;
  typedef typename
    ViewDataType< non_const_value_type , dimension >::type  non_const_type ;

private:

  // A const ?
  enum { is_const = std::is_same< value_type , const_value_type >::value };

  // The unwrapped scalar types:
  typedef typename
    std::conditional< is_const , const ScalarType , ScalarType >::type
      scalar_type ;

  typedef ScalarType        non_const_scalar_type ;
  typedef const ScalarType  const_scalar_type ;

  // Prepend or append the pce dimension based on ArrayLayout
  typedef typename array_analysis::dimension::
    template prepend<0>::type
      prepend_scalar_dimension ;
  typedef typename array_analysis::dimension::
    template append<0>::type
      append_scalar_dimension ;
  typedef typename std::conditional<
    std::is_same< ArrayLayout, Kokkos::LayoutLeft>::value,
    prepend_scalar_dimension,
    append_scalar_dimension >::type scalar_dimension;

public:

  // Generate "flattened" multidimensional array specification type.
  typedef typename
    ViewDataType< scalar_type , scalar_dimension >::type scalar_array_type ;

  typedef typename
    ViewDataType< const_scalar_type , scalar_dimension >::type
      const_scalar_array_type ;

  typedef typename
    ViewDataType< non_const_scalar_type , scalar_dimension >::type
      non_const_scalar_array_type ;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

// UQ::PCE allocation for dynamically-sized UQ::PCE types.
// In this case we allocate two chunks of data, the first for the the
// UQ::PCE<Storage> itself and then for the underlying scalar type
// (UQ::PCE<Storage>::value_type).  The memory is laid out with the
// former followed by the latter.
template <class ValueType>
struct PCEAllocation {
  typedef ValueType value_type;
  typedef typename Sacado::ValueType<value_type>::type scalar_type;
  typedef typename value_type::cijk_type cijk_type;

  value_type  * value_ptr;
  scalar_type * scalar_ptr;

  KOKKOS_INLINE_FUNCTION
  static constexpr size_t
  memory_span(const size_t span, const unsigned pce_size) {
    return span * ( pce_size * sizeof(scalar_type) + sizeof(value_type) );
  }

  KOKKOS_INLINE_FUNCTION
  PCEAllocation() : value_ptr(0), scalar_ptr(0) {}

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  PCEAllocation& operator=(const PCEAllocation<T>& a) {
    value_ptr = a.value_ptr;
    scalar_ptr = a.scalar_ptr;
    return *this;
  }

  // We are making an assumption the data is laid out as described above,
  // which in general may not be true if the view is created from memory
  // allocated elsewhere.  We should check for that.
  KOKKOS_INLINE_FUNCTION
  void set(value_type* ptr, const size_t span, const unsigned pce_size) {
    value_ptr = ptr;
    scalar_ptr = reinterpret_cast<scalar_type*>(ptr+span);
  }

  template <class ExecSpace>
  struct PCEConstruct {
    ExecSpace m_space;
    value_type* m_p;
    scalar_type* m_sp;
    size_t m_span;
    unsigned m_pce_size;
    cijk_type m_cijk;

    PCEConstruct() = default;
    PCEConstruct(const PCEConstruct&) = default;
    PCEConstruct& operator=(const PCEConstruct&) = default;

    inline
    PCEConstruct(const ExecSpace& space,
                 value_type* p,
                 scalar_type* sp,
                 const size_t span,
                 const unsigned pce_size,
                 const cijk_type& cijk) :
      m_space(space), m_p(p), m_sp(sp), m_span(span), m_pce_size(pce_size),
      m_cijk(cijk) {}

    inline void execute() {
      if ( ! m_space.in_parallel() ) {
        typedef Kokkos::RangePolicy< ExecSpace > PolicyType ;
        const Kokkos::Impl::ParallelFor< PCEConstruct , PolicyType >
          closure( *this , PolicyType( 0 , m_span ) );
        closure.execute();
        m_space.fence();
      }
      else {
        for ( size_t i = 0 ; i < m_span ; ++i ) operator()(i);
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_t i) const {
      new (m_p+i) value_type(m_cijk, m_pce_size, m_sp+i*m_pce_size, false);
    }
  };

  template <class ExecSpace>
  struct ConstructDestructFunctor {
    typedef Kokkos::Impl::ViewValueFunctor< ExecSpace, scalar_type > ScalarFunctorType ;
    typedef PCEConstruct< ExecSpace > PCEFunctorType ;
    ScalarFunctorType m_scalar_functor;
    PCEFunctorType m_pce_functor;
    bool m_initialize;

    ConstructDestructFunctor() = default;
    ConstructDestructFunctor(const ConstructDestructFunctor&) = default;
    ConstructDestructFunctor& operator=(const ConstructDestructFunctor&) = default;

    ConstructDestructFunctor(const ExecSpace & space,
                             const bool initialize,
                             const size_t span,
                             const unsigned pce_size,
                             const cijk_type& cijk,
                             scalar_type* scalar_ptr,
                             value_type* value_ptr) :
      m_scalar_functor( space , scalar_ptr , span*pce_size ),
      m_pce_functor( space , value_ptr , scalar_ptr , span , pce_size , cijk ),
      m_initialize(initialize) {}

    inline void construct_shared_allocation() {
      // First initialize the scalar_type array
      if (m_initialize)
        m_scalar_functor.construct_shared_allocation();

      // Construct each UQ::PCE using memory in scalar_ptr array,
      // setting pointer to UQ::PCE values from values array
      // Equivalent to:
      // value_type* p = value_ptr;
      // scalar_type* sp = scalar_ptr;
      // for (size_t i=0; i<span; ++i) {
      //   new (p++) value_type(cijk, pce_size, sp, false);
      //   sp += pce_size;
      // }
      // (we always need to do this, regardless of initialization)
      m_pce_functor.execute();
    }

    inline void destroy_shared_allocation() {
      // We only need to (possibly) call the destructor on values in the
      // scalar_type array, since the value_type array is a view into it
      if (m_initialize)
        m_scalar_functor.destroy_shared_allocation();
    }

  };

  template <class ExecSpace>
  inline ConstructDestructFunctor<ExecSpace>
  create_functor(const ExecSpace & space,
                 const bool initialize,
                 const size_t span,
                 const unsigned pce_size,
                 const cijk_type& cijk) const {
    return ConstructDestructFunctor<ExecSpace>(space, initialize, span,
                                               pce_size, cijk, scalar_ptr,
                                               value_ptr);
  }

  // Assign scalar_type pointer to give ptr
  // This makes BIG assumption on how the data was allocated
  template <typename T>
  void assign(T * ptr) {
    value_ptr  = reinterpret_cast<value_type*>(ptr);
    if (ptr != 0)
      scalar_ptr = value_ptr->coeff();
    else
      scalar_ptr = 0;
  }
};

}}} // namespace Kokkos::Experimental::Impl

namespace Kokkos {
namespace Impl {

template< class Traits >
class ViewMapping< Traits , /* View internal mapping */
  typename std::enable_if<
    ( std::is_same< typename Traits::specialize
                  , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
      &&
      ( std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutLeft >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutRight >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    )>::type >
{
private:

  template< class , class ... > friend class ViewMapping ;
  template< class , class ... > friend class Kokkos::View ;

public:
  typedef typename Traits::value_type  sacado_uq_pce_type ;
  typedef typename sacado_uq_pce_type::storage_type stokhos_storage_type ;
  typedef typename stokhos_storage_type::value_type intrinsic_scalar_type ;
  typedef typename
    std::add_const< intrinsic_scalar_type >::type  const_intrinsic_scalar_type ;
  typedef typename sacado_uq_pce_type::cijk_type        cijk_type ;
private:

  typedef Kokkos::Experimental::Impl::PCEAllocation<sacado_uq_pce_type> handle_type;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  // Prepend or append the pce dimension based on array_layout
  typedef ViewArrayAnalysis< typename Traits::data_type > array_analysis ;
  typedef typename array_analysis::dimension array_dimension;
  typedef ViewOffset< typename array_dimension::
                        template append<0>::type,
                      typename Traits::array_layout,
                      void
                      >  append_offset_type ;
  typedef ViewOffset< typename array_dimension::
                        template prepend<0>::type,
                      typename Traits::array_layout,
                      void
                      >  prepend_offset_type ;
  typedef typename std::conditional<
    std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft>::value,
    prepend_offset_type,
    append_offset_type >::type array_offset_type;

  handle_type      m_handle ;
  offset_type      m_offset ;
  unsigned         m_sacado_size ;   // Size of sacado dimension
  cijk_type        m_cijk ;          // Sparse 3 tensor
  bool             m_is_contiguous ; // Is data allocated contiguously

  // Check whether data allocation is contiguous
  // Since View() takes an arbitrary pointer, we can't necessarily assume
  // the data was allocated contiguously
  KOKKOS_INLINE_FUNCTION
  bool is_data_contiguous() const {
    const size_t sz = this->span();
    if (sz == 0)
      return true;
    const intrinsic_scalar_type* last_coeff =
      m_handle.value_ptr[sz-1].coeff();
    const intrinsic_scalar_type* last_coeff_expected =
      m_handle.scalar_ptr + (sz-1)*m_sacado_size;
    return last_coeff == last_coeff_expected;
  }

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  // Rank corresponding to the sacado dimension
  enum { Sacado_Rank = std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft >::value ? 0 : Rank+1 };

  // Using the internal offset mapping so limit to public rank:
  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr size_t extent( const iType & r ) const
    { return m_offset.m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr
  typename Traits::array_layout layout() const
    { return m_offset.layout(); }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const
    { return m_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const
    { return m_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const
    { return m_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const
    { return m_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const
    { return m_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const
    { return m_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const
    { return m_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const
    { return m_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  // Since we all for striding within the data type, we can't guarantee
  // regular striding
  using is_regular = std::false_type ;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const
    { return m_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const
    { return m_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const
    { return m_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const
    { return m_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const
    { return m_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const
    { return m_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const
    { return m_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const
    { return m_offset.stride_7(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const
    { m_offset.stride(s); }

  // Size of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned dimension_scalar() const
    { return m_sacado_size; }

  // Sparse tensor
  KOKKOS_FORCEINLINE_FUNCTION
  cijk_type cijk() const
    { return m_cijk; }

  // Sparse tensor
  KOKKOS_FORCEINLINE_FUNCTION
  void set_cijk(const cijk_type& cijk)
    { m_cijk = cijk; }

  // Is allocation contiguous
  KOKKOS_INLINE_FUNCTION
  bool is_allocation_contiguous() const
    { return m_is_contiguous; }

  // Whether the storage type is statically sized
  static const bool is_static = false ;

  // Whether sacado dimension is contiguous
  static const bool is_contiguous = true;

  //----------------------------------------
  // Range of mapping

  // Return type of reference operators
  typedef sacado_uq_pce_type & reference_type ;

  /** \brief Pointer to underlying memory type */
  typedef sacado_uq_pce_type * pointer_type ;

  /** \brief  Span of the mapped range : [ data() .. data() + span() ) */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const
    { return m_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_offset.span_is_contiguous() ; }

  /** \brief Raw data access */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const
    { return m_handle.value_ptr ; }

  //----------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const
    { return *m_handle.value_ptr; }

  // FIXME:  Check this
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    ! std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const
    { return m_handle.value_ptr[i0]; }

  // FIXME:  Check this
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const
    { return m_handle.value_ptr[ m_offset(i0) ]; }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1) ]; }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2,i3) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2,i3,i4) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2,i3,i4,i5) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2,i3,i4,i5,i6) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 ) const
    { return m_handle.value_ptr[ m_offset(i0,i1,i2,i3,i4,i5,i6,i7) ]; }

  //----------------------------------------

  /** \brief  Span, in bytes, of the required memory */
  KOKKOS_INLINE_FUNCTION
  static size_t memory_span( typename Traits::array_layout const & layout )
    {
      // Do not introduce padding...
      typedef std::integral_constant< unsigned , 0 >  padding ;
      offset_type offset( padding(), layout );
      unsigned sacado_size =
        Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout);
      return handle_type::memory_span( offset.span(), sacado_size );
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION ~ViewMapping() = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping() :
    m_handle(),
    m_offset(),
    m_sacado_size(0),
    m_cijk(),
    m_is_contiguous(true)
    {}

  KOKKOS_INLINE_FUNCTION ViewMapping( const ViewMapping & ) = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( const ViewMapping & ) = default ;

  KOKKOS_INLINE_FUNCTION ViewMapping( ViewMapping && ) = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( ViewMapping && ) = default ;

  template< class ... P >
  KOKKOS_INLINE_FUNCTION
  ViewMapping
    ( ViewCtorProp< P ... > const & prop
    , typename Traits::array_layout const & layout
    )
    : m_handle()
    , m_offset( std::integral_constant< unsigned , 0 >() , layout )
    , m_sacado_size( Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout) )
    {
      m_handle.set( ( (ViewCtorProp<void,pointer_type> const &) prop ).value
                    , m_offset.span(), m_sacado_size );
      m_cijk = extract_cijk<cijk_type>(prop);
#ifndef __CUDA_ARCH__
      if (m_cijk.dimension() == 0)
        m_cijk = getGlobalCijkTensor<cijk_type>();
      // Use 0 or 1 to signal the size wasn't specified in the constructor,
      // since now dimesion_i() == 1 for all i >= rank
      if (m_sacado_size == 0 || m_sacado_size == 1)
        m_sacado_size = m_cijk.dimension();
#endif
      m_is_contiguous = this->is_data_contiguous();
    }

  //----------------------------------------
  /*  Allocate and construct mapped array.
   *  Allocate via shared allocation record and
   *  return that record for allocation tracking.
   */
  template< class ... P >
  SharedAllocationRecord<> *
  allocate_shared( ViewCtorProp< P... > const & prop
                 , typename Traits::array_layout const & layout )
  {
    typedef ViewCtorProp< P... > ctor_prop ;

    typedef typename ctor_prop::execution_space  execution_space ;
    typedef typename Traits::memory_space         memory_space ;
    typedef typename handle_type::template ConstructDestructFunctor<execution_space> functor_type ;
    typedef SharedAllocationRecord< memory_space , functor_type > record_type ;

    // Disallow padding
    typedef std::integral_constant< unsigned , 0 > padding ;

    m_offset = offset_type( padding(), layout );
    m_sacado_size = Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout);
    m_cijk = extract_cijk<cijk_type>(prop);
    if (m_cijk.dimension() == 0)
      m_cijk = getGlobalCijkTensor<cijk_type>();
    // Use 0 or 1 to signal the size wasn't specified in the constructor,
    // since now dimesion_i() == 1 for all i >= rank
    if (m_sacado_size == 0 || m_sacado_size == 1)
      m_sacado_size = m_cijk.dimension();
    m_is_contiguous = true;

    const size_t alloc_size =
      handle_type::memory_span( m_offset.span(), m_sacado_size );

    // Create shared memory tracking record with allocate memory from the memory space
    record_type * const record =
      record_type::allocate( ( (ViewCtorProp<void,memory_space> const &) prop ).value
                           , ( (ViewCtorProp<void,std::string>  const &) prop ).value
                           , alloc_size );

    //  Only set the the pointer and initialize if the allocation is non-zero.
    //  May be zero if one of the dimensions is zero.
    if ( alloc_size ) {

      m_handle.set( reinterpret_cast< pointer_type >( record->data() ),
                    m_offset.span(), m_sacado_size );

      // Assume destruction is only required when construction is requested.
      // The ViewValueFunctor has both value construction and destruction operators.
      record->m_destroy = m_handle.create_functor(
        ( (ViewCtorProp<void,execution_space> const &) prop).value
        , ctor_prop::initialize
        , m_offset.span()
        , m_sacado_size
        , m_cijk );

      // Construct values
      record->m_destroy.construct_shared_allocation();
    }

    return record ;
  }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Assign compatible Sacado::UQ::PCE view mappings.
 *
 *  View<UQ::PCE> = View<UQ::PCE>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::value
    &&
    // Destination view has UQ::PCE
    std::is_same< typename DstTraits::specialize
                , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
    &&
    // Source view has UQ::PCE only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of UQ::PCE requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ||
        ( unsigned(DstTraits::rank) == 1 && unsigned(SrcTraits::rank) == 1 ) ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable
          < typename DstType::offset_type::dimension_type
          , typename SrcType::offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

      dst.m_handle  = src.m_handle ;
      dst.m_offset  = src.m_offset ;
      dst.m_sacado_size = src.m_sacado_size ;
      dst.m_cijk    = src.m_cijk ;
      dst.m_is_contiguous = src.m_is_contiguous ;
    }
};

/**\brief  Assign compatible Sacado::UQ::PCE view mappings.
 *
 *  View<ordinary> = View<UQ::PCE>
 *  where View<ordinay>::Rank = View<UQ::PCE>::Rank+1
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::value
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has UQ::PCE only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
    &&
    // Ranks match
    unsigned(DstTraits::dimension::rank) == unsigned(SrcTraits::dimension::rank)+1
  )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of UQ::PCE requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::scalar_array_type
                    , typename SrcTraits::scalar_array_type >::value ||
        std::is_same< typename DstTraits::scalar_array_type
                    , typename SrcTraits::const_scalar_array_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable<
          typename DstType::offset_type::dimension_type,
          typename SrcType::array_offset_type::dimension_type >::value,
        "View assignment must have compatible dimensions" );

      if ( !src.m_is_contiguous )
        Kokkos::abort("\n\n ****** Kokkos::View< Sacado::UQ::PCE ... >:  can't assign non-contiguous view ******\n\n");

      unsigned dims[8];
      dims[0] = src.m_offset.dimension_0();
      dims[1] = src.m_offset.dimension_1();
      dims[2] = src.m_offset.dimension_2();
      dims[3] = src.m_offset.dimension_3();
      dims[4] = src.m_offset.dimension_4();
      dims[5] = src.m_offset.dimension_5();
      dims[6] = src.m_offset.dimension_6();
      dims[7] = src.m_offset.dimension_7();
      unsigned rank = SrcTraits::dimension::rank;
      unsigned sacado_size = src.m_sacado_size;
      if (std::is_same<typename SrcTraits::array_layout, LayoutLeft>::value) {
        // Move sacado_size to the first dimension, shift all others up one
        for (unsigned i=rank; i>0; --i)
          dims[i] = dims[i-1];
        dims[0] = sacado_size;
      }
      else {
        dims[rank] = sacado_size;
      }
      typedef typename DstType::offset_type dst_offset_type;
      dst.m_offset = dst_offset_type( std::integral_constant< unsigned , 0 >(),
                                      typename DstTraits::array_layout(
                                        dims[0] , dims[1] , dims[2] , dims[3] ,
                                        dims[4] , dims[5] , dims[6] , dims[7] ) );
      dst.m_handle  = src.m_handle.scalar_ptr ;
    }
};

/**\brief  Assign compatible Sacado::UQ::PCE view mappings.
 *
 *  View<ordinary> = View<UQ::PCE>
 *  where View<ordinay>::Rank = View<UQ::PCE>::Rank, i.e., assigning
 *  to the "flattened" view type
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::value
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has UQ::PCE only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
    &&
    // Ranks match
    unsigned(DstTraits::dimension::rank) == unsigned(SrcTraits::dimension::rank)
    )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of UQ::PCE requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::non_const_value_type::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , const typename SrcTraits::non_const_value_type::value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable<
          typename DstType::offset_type::dimension_type,
          typename SrcType::offset_type::dimension_type >::value,
        "View assignment must have compatible dimensions" );

       if ( !src.m_is_contiguous )
        Kokkos::abort("\n\n ****** Kokkos::View< Sacado::UQ::PCE ... >:  can't assign non-contiguous view ******\n\n");

      unsigned dims[8];
      dims[0] = src.m_offset.dimension_0();
      dims[1] = src.m_offset.dimension_1();
      dims[2] = src.m_offset.dimension_2();
      dims[3] = src.m_offset.dimension_3();
      dims[4] = src.m_offset.dimension_4();
      dims[5] = src.m_offset.dimension_5();
      dims[6] = src.m_offset.dimension_6();
      dims[7] = src.m_offset.dimension_7();
      unsigned rank = SrcTraits::dimension::rank;
      unsigned sacado_size = src.m_sacado_size;
      if (std::is_same<typename DstTraits::array_layout, LayoutLeft>::value) {
        dims[0] = dims[0]*sacado_size;
        dims[rank] = 0;
      }
      else {
        dims[rank-1] = dims[rank-1]*sacado_size;
        dims[rank] = 0;
      }
      typedef typename DstType::offset_type dst_offset_type;
      dst.m_offset = dst_offset_type( std::integral_constant< unsigned , 0 >(),
                                      typename DstTraits::array_layout(
                                        dims[0] , dims[1] , dims[2] , dims[3] ,
                                        dims[4] , dims[5] , dims[6] , dims[7] ) );
      dst.m_handle  = src.m_handle.scalar_ptr ;
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Subview mapping

template< class DataType, class ... P , class Arg0, class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      // Source view has UQ::PCE only
      std::is_same< typename Kokkos::ViewTraits<DataType,P...>::specialize
                  , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
      &&
      (
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutStride >::value
      )
    )>::type
  , Kokkos::ViewTraits<DataType,P...>
  , Arg0, Args ... >
{
private:

  typedef Kokkos::ViewTraits<DataType,P...> SrcTraits;

  //static_assert( SrcTraits::rank == sizeof...(Args) , "" );

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Arg0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Arg0,Args...>::value)
    , R2 = bool(is_integral_extent<2,Arg0,Args...>::value)
    , R3 = bool(is_integral_extent<3,Arg0,Args...>::value)
    , R4 = bool(is_integral_extent<4,Arg0,Args...>::value)
    , R5 = bool(is_integral_extent<5,Arg0,Args...>::value)
    , R6 = bool(is_integral_extent<6,Arg0,Args...>::value)
    };

  // Public rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  // Whether right-most non-UQ::PCE rank is a range.
  enum { R0_rev = ( 0 == SrcTraits::rank ? RZ : (
                    1 == SrcTraits::rank ? R0 : (
                    2 == SrcTraits::rank ? R1 : (
                    3 == SrcTraits::rank ? R2 : (
                    4 == SrcTraits::rank ? R3 : (
                    5 == SrcTraits::rank ? R4 : (
                    6 == SrcTraits::rank ? R5 : R6 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1 or 2, InputLayout Left, Interval 0
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0 && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0_rev && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value )
      ), typename SrcTraits::array_layout , Kokkos::LayoutStride
      >::type array_layout ;

  typedef typename SrcTraits::value_type  sacado_uq_pce_type ;

  typedef typename std::conditional< rank == 0 , sacado_uq_pce_type ,
          typename std::conditional< rank == 1 , sacado_uq_pce_type * ,
          typename std::conditional< rank == 2 , sacado_uq_pce_type ** ,
          typename std::conditional< rank == 3 , sacado_uq_pce_type *** ,
          typename std::conditional< rank == 4 , sacado_uq_pce_type **** ,
          typename std::conditional< rank == 5 , sacado_uq_pce_type ***** ,
          typename std::conditional< rank == 6 , sacado_uq_pce_type ****** ,
                                                 sacado_uq_pce_type *******
          >::type >::type >::type >::type >::type >::type >::type
    data_type ;

public:

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  // The presumed type is 'ViewMapping< traits_type , void >'
  // However, a compatible ViewMapping is acceptable.
  template< class DstTraits >
  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< DstTraits , void > & dst
                    , ViewMapping< SrcTraits , void > const & src
                    , Arg0 arg0, Args ... args )
    {
      static_assert(
        ViewMapping< DstTraits , traits_type , void >::is_assignable ,
        "Subview destination type must be compatible with subview derived type" );

      typedef ViewMapping< DstTraits , void > DstType ;
      typedef typename DstType::offset_type  dst_offset_type ;

      const SubviewExtents< SrcTraits::rank , rank >
        extents( src.m_offset.m_dim , arg0 , args... );

      const size_t offset = src.m_offset( extents.domain_offset(0)
                                          , extents.domain_offset(1)
                                          , extents.domain_offset(2)
                                          , extents.domain_offset(3)
                                          , extents.domain_offset(4)
                                          , extents.domain_offset(5)
                                          , extents.domain_offset(6)
                                          , extents.domain_offset(7) );

      dst.m_offset = dst_offset_type( src.m_offset , extents );
      dst.m_handle.value_ptr = src.m_handle.value_ptr + offset;
      dst.m_handle.scalar_ptr =
        src.m_handle.scalar_ptr + offset * src.m_sacado_size;
      dst.m_sacado_size = src.m_sacado_size;
      dst.m_cijk = src.m_cijk;
      dst.m_is_contiguous = src.m_is_contiguous;
    }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Specialization for deep_copy( view, view::value_type ) for Cuda
#if defined( KOKKOS_HAVE_CUDA )
template< class OutputView >
struct StokhosViewFill< OutputView ,
                 typename std::enable_if< std::is_same< typename OutputView::specialize,
                                                        Kokkos::Experimental::Impl::ViewPCEContiguous >::value &&
                                     std::is_same< typename OutputView::execution_space,
                                                   Cuda >::value >::type >
{
  typedef typename OutputView::const_value_type   const_value_type ;
  typedef typename Sacado::ScalarType<const_value_type>::type scalar_type ;
  typedef typename OutputView::execution_space    execution_space ;
  typedef typename OutputView::size_type          size_type ;

  template <unsigned VectorLength>
  struct PCEKernel {
    typedef typename OutputView::execution_space execution_space ;
    const OutputView output;
    const_value_type input;

    PCEKernel( const OutputView & arg_out , const_value_type & arg_in ) :
      output(arg_out), input(arg_in) {}

    typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;

    KOKKOS_INLINE_FUNCTION
    void operator()( const team_member & dev ) const
    {
      const size_type tidx = dev.team_rank() % VectorLength;
      const size_type tidy = dev.team_rank() / VectorLength;
      const size_type nrow = dev.team_size() / VectorLength;
      const size_type nvec = dimension_scalar(output);

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
        output(i0,i1,i2,i3,i4,i5,i6,i7).fastAccessCoeff(is) =
          input.fastAccessCoeff(is) ;
      }}}}}}}}
    }
  };

  template <unsigned VectorLength>
  struct ScalarKernel {
    typedef typename OutputView::execution_space execution_space ;
    const OutputView  output;
    const scalar_type input;

    ScalarKernel( const OutputView & arg_out , const scalar_type & arg_in ) :
      output(arg_out), input(arg_in) {}

    typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;
    KOKKOS_INLINE_FUNCTION
    void operator()( const team_member & dev ) const
    {
      const size_type tidx = dev.team_rank() % VectorLength;
      const size_type tidy = dev.team_rank() / VectorLength;
      const size_type nrow = dev.team_size() / VectorLength;
      const size_type npce = dimension_scalar(output);

      const size_type i0 = dev.league_rank() * nrow + tidy;
      if ( i0 >= output.dimension_0() ) return;

      for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
      for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
      for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
      for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
      for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
      for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
      for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      for ( size_type is = tidx ; is < npce ; is+=VectorLength ) {
        output(i0,i1,i2,i3,i4,i5,i6,i7).fastAccessCoeff(is) =
          is == 0 ? input : scalar_type(0) ;
      }}}}}}}}
    }
  };

  StokhosViewFill( const OutputView & output , const_value_type & input )
  {
    // Coalesced accesses are 128 bytes in size
    typedef typename OutputView::array_type::value_type scalar_type;
    const unsigned vector_length =
      ( 128 + sizeof(scalar_type)-1 ) / sizeof(scalar_type);

    // 8 warps per block should give good occupancy
    const size_type block_size = 256;

    const size_type rows_per_block = block_size / vector_length;
    const size_type n = output.dimension_0();
    const size_type league_size = ( n + rows_per_block-1 ) / rows_per_block;
    const size_type team_size = rows_per_block * vector_length;
    Kokkos::TeamPolicy< execution_space > config( league_size, team_size );

    if (static_cast<unsigned>(input.size()) != dimension_scalar(output) &&
        input.size() != 1)
      Kokkos::abort("StokhosViewFill:  Invalid input value size");

    if (input.size() == 1)
      parallel_for(
        config, ScalarKernel<vector_length>(output, input.fastAccessCoeff(0)) );
    else
      parallel_for( config, PCEKernel<vector_length>(output, input) );
    execution_space::fence();
  }

  StokhosViewFill( const OutputView & output , const scalar_type & input )
  {
    // Coalesced accesses are 128 bytes in size
    typedef typename OutputView::array_type::value_type scalar_type;
    const unsigned vector_length =
      ( 128 + sizeof(scalar_type)-1 ) / sizeof(scalar_type);

    // 8 warps per block should give good occupancy
    const size_type block_size = 256;

    const size_type rows_per_block = block_size / vector_length;
    const size_type n = output.dimension_0();
    const size_type league_size = ( n + rows_per_block-1 ) / rows_per_block;
    const size_type team_size = rows_per_block * vector_length;
    Kokkos::TeamPolicy< execution_space > config( league_size, team_size );

    parallel_for( config, ScalarKernel<vector_length>(output, input) );
    execution_space::fence();
  }

};
#endif /* #if defined( KOKKOS_HAVE_CUDA ) */

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_UQ_PCE_CONTIGUOUS_HPP */
