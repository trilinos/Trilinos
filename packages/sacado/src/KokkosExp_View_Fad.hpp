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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOSCORE)

// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Layout.hpp"
#include "Kokkos_View.hpp"

// Some definition that should exist whether the specializations exist or not

namespace Kokkos {

// Whether a given type is a view with Sacado FAD scalar type
template <typename view_type>
struct is_view_fad { static const bool value = false; };

// Whether a given type is a view with Sacado FAD scalar type with contiguous
// layout
template <typename view_type>
struct is_view_fad_contiguous { static const bool value = false; };

// Template function for extracting sacado dimension
template <typename view_type>
KOKKOS_INLINE_FUNCTION
constexpr unsigned
dimension_scalar(const view_type& view) {
  return 0;
}

// Template function for extracting aligned sacado dimension
template <typename view_type>
KOKKOS_INLINE_FUNCTION
constexpr unsigned
dimension_scalar_aligned(const view_type& view) {
  return dimension_scalar(view);
}

}

// Make sure the user really wants these View specializations
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct ViewSpecializeSacadoFad {};
struct ViewSpecializeSacadoFadContiguous {};

template< class ... Args >
struct is_ViewSpecializeSacadoFad { enum { value = false }; };

template< class D , class ... P , class ... Args >
struct is_ViewSpecializeSacadoFad< Kokkos::View<D,P...> , Args... > {
  enum { value =
    std::is_same< typename Kokkos::ViewTraits<D,P...>::specialize
                , ViewSpecializeSacadoFad >::value
    &&
    ( ( sizeof...(Args) == 0 ) ||
      is_ViewSpecializeSacadoFad< Args... >::value ) };
};

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template <typename T, typename ... P>
struct is_view_fad< View<T,P...> > {
  typedef View<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFadContiguous >::value;
};

template <typename T, typename ... P>
struct is_view_fad_contiguous< View<T,P...> > {
  typedef View<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFadContiguous >::value;
};

}

namespace Kokkos {
namespace Impl {

// Overload view_copy for Fad View's:
//   1.  Should be faster than using Fad directly
//   2.  Fixes issues with hierarchical parallelism since the default
//       implementation uses MDRangePolicy which doesn't work with hierarchical
//       parallelism.
// Needs to go before include of Kokkos_Core.hpp so it is in scope when
// Kokkos_CopyViews.hpp is included by Kokkos_Core.hpp, which internally
// calls view_copy().
template<class DT, class ... DP,
         class ST, class ... SP>
typename std::enable_if< is_view_fad< Kokkos::View<DT,DP...> >::value &&
                         is_view_fad< Kokkos::View<ST,SP...> >::value
                       >::type
view_copy(const Kokkos::View<DT,DP...>& dst, const Kokkos::View<ST,SP...>& src);

template<class Space, class T, class ... P>
struct MirrorType;

} // namespace Impl

// Declare overloads of create_mirror() so they are in scope
// Kokkos_Core.hpp is included below

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror(
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    !std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
        Kokkos::LayoutStride >::value >::type * = 0);


template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror(
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value >::type * = 0);

template<class Space, class T, class ... P>
typename Impl::MirrorType<Space,T,P ...>::view_type
create_mirror(
  const Space&,
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value >::type * = 0);

} // namespace Kokkos

#include "Sacado_Traits.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_LayoutContiguous.hpp"
#include "Kokkos_LayoutNatural.hpp"

namespace Kokkos {
namespace Impl {

// Define our overload of view_copy above.  Needs to happen after including
// Kokkos_Core.hpp since it calls the default implementation
template<class DT, class ... DP,
         class ST, class ... SP>
typename std::enable_if< is_view_fad< Kokkos::View<DT,DP...> >::value &&
                         is_view_fad< Kokkos::View<ST,SP...> >::value
                       >::type
view_copy(const Kokkos::View<DT,DP...>& dst, const Kokkos::View<ST,SP...>& src)
{
  typedef typename Kokkos::View<DT,DP...>::array_type dst_array_type;
  typedef typename Kokkos::View<ST,SP...>::array_type src_array_type;
  view_copy( dst_array_type(dst) , src_array_type(src) );
}

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_view_fad< View<T,P...> >::value, unsigned >::type
dimension_scalar(const View<T,P...>& view) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  return view.implementation_map().dimension_scalar();
#else
  return view.impl_map().dimension_scalar();
#endif
}

template <typename Layout>
struct ApplyNatural {
  typedef LayoutNatural<Layout> type;
};

template <typename Layout>
struct ApplyNatural< LayoutNatural<Layout> > {
  typedef LayoutNatural<Layout> type;
};

template < typename T, typename Enable = void >
struct ArrayScalar;

template < typename T >
struct ArrayScalar< T, typename std::enable_if< !Sacado::IsFad<T>::value >::type > {
  typedef T type;
};

template < typename T >
struct ArrayScalar< T, typename std::enable_if< Sacado::IsFad<T>::value >::type > {
  typedef typename ArrayScalar< typename Sacado::ValueType<T>::type >::type* type;
};


template < typename DataType, int Rank >
struct AppendRankToConvertedFad {
  static_assert( Rank > -1, "Sacado AppendRankToConvertedFad Error: Rank < 0" );
  typedef typename AppendRankToConvertedFad<DataType,Rank-1>::type* type;
};

// terminating specialization
template < typename DataType >
struct AppendRankToConvertedFad< DataType, 0 > {
  typedef DataType type;
};


template < class ArrayLayout, class Enable = void >
struct ViewArrayLayoutSelector;

template < class ArrayLayout >
struct ViewArrayLayoutSelector< ArrayLayout, typename std::enable_if< std::is_same<ArrayLayout, Kokkos::LayoutLeft>::value >::type >
{
  using type = Kokkos::LayoutLeft;
};

template < class ArrayLayout >
struct ViewArrayLayoutSelector< ArrayLayout, typename std::enable_if< std::is_same<ArrayLayout, Kokkos::LayoutRight>::value >::type >
{
  using type = Kokkos::LayoutRight;
};

template < class ArrayLayout >
struct ViewArrayLayoutSelector< ArrayLayout, typename std::enable_if< std::is_same<ArrayLayout, Kokkos::LayoutStride>::value >::type >
{
  using type = Kokkos::LayoutStride;
};

template < typename ViewType, typename Enable = void >
struct PODViewDeepCopyType;

template < typename ViewType >
struct PODViewDeepCopyType< ViewType, typename std::enable_if< is_view_fad<ViewType>::value >::type >
{

  typedef ViewType view_type;
  typedef typename ArrayScalar< typename view_type::value_type >::type fad_converted_type;
  typedef typename AppendRankToConvertedFad< fad_converted_type, view_type::Rank >::type new_data_type;

  typedef typename ViewArrayLayoutSelector<typename view_type::array_layout>::type layout;
  //typedef typename view_type::array_layout layout;
  typedef typename view_type::device_type device;
  typedef typename view_type::memory_traits memory;

  typedef Kokkos::View< new_data_type, layout, device, memory > type;
};

// Not a Fad type
template < typename ViewType >
struct PODViewDeepCopyType< ViewType, typename std::enable_if< !is_view_fad<ViewType>::value >::type > 
{
  typedef ViewType type;
};


template <typename ViewType, typename Enabled = void>
struct NaturalArrayType {
  typedef ViewType type;
};

template <typename D, typename ... P>
struct NaturalArrayType< View<D,P...>,
                            typename std::enable_if< is_view_fad< View<D,P...> >::value >::type > {
  typedef View<D,P...> view_type;
  typedef typename view_type::data_type data_type;
  typedef typename view_type::array_layout layout;
  typedef typename view_type::device_type device;
  typedef typename view_type::memory_traits memory;
  //typedef typename ApplyNatural<layout>::type natural_layout;
  typedef typename ViewArrayLayoutSelector<layout>::type natural_layout;
  typedef View<data_type,natural_layout,device,memory> type;
};

namespace Impl {

template <class OutputView, typename Enabled = void>
struct SacadoViewFill
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::execution_space execution_space ;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    const size_t n1 = output.extent(1);
    const size_t n2 = output.extent(2);
    const size_t n3 = output.extent(3);
    const size_t n4 = output.extent(4);
    const size_t n5 = output.extent(5);
    const size_t n6 = output.extent(6);
    const size_t n7 = output.extent(7);

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_t i7 = 0 ; i7 < n7 ; ++i7 ) {
      output.access(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }

  SacadoViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      const size_t n0 = output.extent(0);
      Kokkos::RangePolicy<execution_space> policy( 0, n0 );
      Kokkos::parallel_for( policy, *this );
      execution_space::fence();
    }
};

}

// Overload of deep_copy for Fad views intializing to a constant scalar
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename Sacado::ScalarType< typename View<DT,DP...>::value_type >::type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  Impl::SacadoViewFill< View<DT,DP...> >( view , value );
}


// Overload of deep_copy for Fad views intializing to a constant Fad
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  Impl::SacadoViewFill< View<DT,DP...> >( view , value );
}


/* Specialize for deep copy of FAD */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  ( std::is_same< typename ViewTraits<DT,DP...>::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value
    ||
    std::is_same< typename ViewTraits<DT,DP...>::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
  &&
  ( std::is_same< typename ViewTraits<ST,SP...>::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value
    ||
    std::is_same< typename ViewTraits<ST,SP...>::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
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

#if 0
  // Current impl
  typedef typename View<DT,DP...>::array_type dst_array_type;
  typedef typename View<ST,SP...>::array_type src_array_type;
  typename NaturalArrayType< dst_array_type >::type dst_array( dst );
  typename NaturalArrayType< src_array_type >::type src_array( src );
#else
  // Copy-assign Views of FadType to Kokkos Views to use Kokkos' deep_copy routine
  typename PODViewDeepCopyType< View<DT,DP...> >::type dst_array( dst );
  typename PODViewDeepCopyType< View<ST,SP...> >::type src_array( src );
#endif
  Kokkos::deep_copy( dst_array , src_array );
}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror( const Kokkos::View<T,P...> & src
             , typename std::enable_if<
                 ( std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
                   std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
               &&
                 ! std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout
                               , Kokkos::LayoutStride >::value
               >::type *
             )
{
  typedef View<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  typename src_type::array_layout layout = src.layout();
  layout.dimension[src_type::rank] = Kokkos::dimension_scalar(src);

  return dst_type(std::string(src.label()).append("_mirror"), layout);
}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror( const Kokkos::View<T,P...> & src
             , typename std::enable_if<
                 ( std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
                   std::is_same< typename ViewTraits<T,P...>::specialize ,
                   Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
                &&
                   std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout
                               , Kokkos::LayoutStride >::value
               >::type *
             )
{
  typedef View<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  Kokkos::LayoutStride layout ;

  layout.dimension[0] = src.extent(0);
  layout.dimension[1] = src.extent(1);
  layout.dimension[2] = src.extent(2);
  layout.dimension[3] = src.extent(3);
  layout.dimension[4] = src.extent(4);
  layout.dimension[5] = src.extent(5);
  layout.dimension[6] = src.extent(6);
  layout.dimension[7] = src.extent(7);

  layout.stride[0] = src.stride_0();
  layout.stride[1] = src.stride_1();
  layout.stride[2] = src.stride_2();
  layout.stride[3] = src.stride_3();
  layout.stride[4] = src.stride_4();
  layout.stride[5] = src.stride_5();
  layout.stride[6] = src.stride_6();
  layout.stride[7] = src.stride_7();

  layout.dimension[src_type::rank] = Kokkos::dimension_scalar(src);

  return dst_type(std::string(src.label()).append("_mirror"), layout);
}

template<class Space, class T, class ... P>
typename Impl::MirrorType<Space,T,P ...>::view_type
create_mirror(const Space& , const Kokkos::View<T,P...> & src
             , typename std::enable_if<
                 ( std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
                   std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
               >::type *) {
  typedef View<T,P...> src_type ;
  typename src_type::array_layout layout = src.layout();
  layout.dimension[src_type::rank] = Kokkos::dimension_scalar(src);
  return typename Impl::MirrorType<Space,T,P ...>::view_type(src.label(),layout);
}

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DataType , class ArrayLayout , class ScalarType , unsigned DimFad >
struct FadViewDataAnalysis
{
private:

  typedef ViewArrayAnalysis< DataType > array_analysis ;

public:

  // Specialized view data mapping:
  typedef ViewSpecializeSacadoFad specialize ;

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

  // Append the FAD static dimension
  typedef typename array_analysis::dimension::
    template append<( DimFad ? DimFad + 1 : 0 )>::type
      scalar_dimension ;

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

// Specialization for LayoutContiguous, where the Fad type is kept contiguous.
// This requires a separate view specialization.
template< class DataType , class ArrayLayout , class ScalarType , unsigned DimFad, unsigned Stride >
struct FadViewDataAnalysis<DataType, LayoutContiguous<ArrayLayout,Stride>, ScalarType, DimFad>
{
private:

  typedef ViewArrayAnalysis< DataType > array_analysis ;

public:

  // For now use the default mapping
  typedef ViewSpecializeSacadoFadContiguous specialize ;

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

  // Prepend/append the FAD dimension
  typedef typename std::conditional<
    std::is_same< ArrayLayout, Kokkos::LayoutLeft >::value,
    typename array_analysis::dimension::
      template prepend<0>::type,
    typename array_analysis::dimension::
      template append<0>::type >::type
    scalar_dimension ;

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

// Specialization for LayoutNatural, where we don't allow striding within
// the FadType.
//
// Currently this is implemented by choosing the default ViewMapping
// specialization.
template< class DataType , class ArrayLayout , class ScalarType , unsigned DimFad >
struct FadViewDataAnalysis<DataType, LayoutNatural<ArrayLayout>, ScalarType, DimFad>
{
private:

  typedef ViewArrayAnalysis< DataType > array_analysis ;

public:

  // For now use the default mapping
  typedef void specialize ;

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

  // Generate "flattened" multidimensional array specification type.
  typedef type            scalar_array_type ;
  typedef const_type      const_scalar_array_type ;
  typedef non_const_type  non_const_scalar_array_type ;

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Sacado {

namespace Fad         { namespace Exp { template< typename > class GeneralFad ; } }

#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
namespace Fad         { template< typename > class DFad ; }
namespace Fad         { template< typename , int > class SFad ; }
namespace Fad         { template< typename , int > class SLFad ; }
#endif

namespace CacheFad    { template< typename > class DFad ; }
namespace ELRFad      { template< typename > class DFad ; }
namespace ELRCacheFad { template< typename > class DFad ; }

namespace CacheFad    { template< typename , int > class SFad ; }
namespace ELRFad      { template< typename , int > class SFad ; }
namespace ELRCacheFad { template< typename , int > class SFad ; }


namespace CacheFad    { template< typename , int > class SLFad ; }
namespace ELRFad      { template< typename , int > class SLFad ; }
namespace ELRCacheFad { template< typename , int > class SLFad ; }
}

namespace Kokkos {
namespace Impl {

#define KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD( NS ) \
template< class DataType , class ArrayLayout , typename ScalarType > \
struct ViewDataAnalysis \
  < DataType     /* Original view data type */ \
  , ArrayLayout \
  , Sacado:: NS ::DFad< ScalarType > \
  > : public FadViewDataAnalysis< DataType, ArrayLayout, ScalarType , 0 > {}; \
\
template< class DataType , class ArrayLayout , typename ScalarType , int N > \
struct ViewDataAnalysis \
  < DataType     /* Original view data type */ \
  , ArrayLayout \
  , Sacado:: NS ::SFad< ScalarType , N > \
  > : public FadViewDataAnalysis< DataType, ArrayLayout, ScalarType , \
       int(Sacado::StaticSize< Sacado:: NS ::SFad< ScalarType , N > >::value) \
       > {}; \
\
template< class DataType , class ArrayLayout , typename ScalarType , int N > \
struct ViewDataAnalysis \
  < DataType     /* Original view data type */ \
  , ArrayLayout \
  , Sacado:: NS ::SLFad< ScalarType , N > \
  > : public FadViewDataAnalysis< DataType, ArrayLayout, ScalarType , \
       int(Sacado::StaticSize< Sacado:: NS ::SLFad< ScalarType , N > >::value) \
       > {}; \

template< class DataType , class ArrayLayout , typename StorageType >
struct ViewDataAnalysis
  < DataType     /* Original view data type */
  , ArrayLayout
  , Sacado::Fad::Exp::GeneralFad< StorageType >
    > : public FadViewDataAnalysis< DataType, ArrayLayout, typename StorageType::value_type , 0 > {};

#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD( Fad )
#endif

KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD( CacheFad )
KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD( ELRFad )
KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD( ELRCacheFad )

#undef KOKKOS_VIEW_DATA_ANALYSIS_SACADO_FAD

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

// Copied from Sacado_ViewFactory
template <class View, class ... ViewPack>
KOKKOS_INLINE_FUNCTION
unsigned dimension_scalar(const View& v, const ViewPack&... views) {
  const unsigned dim0 = dimension_scalar(v);
  const unsigned dim1 = dimension_scalar(views...);
  return dim0 >= dim1 ? dim0 : dim1 ;
}

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos { namespace Impl {

template < typename Specialize, typename A, typename B >
struct CommonViewValueType;

template < typename A, typename B >
struct CommonViewValueType< Kokkos::Impl::ViewSpecializeSacadoFad, A, B >
{
  using value_type = typename Sacado::Promote<A,B>::type ;
};

template < typename A, typename B >
struct CommonViewValueType< Kokkos::Impl::ViewSpecializeSacadoFadContiguous, A, B >
{
  using value_type = typename Sacado::Promote<A,B>::type ;
};


template < class Specialize, class ValueType >
struct CommonViewAllocProp;

template < class ValueType >
struct CommonViewAllocProp< Kokkos::Impl::ViewSpecializeSacadoFad, ValueType >
{
  using value_type = ValueType;
  using scalar_array_type = typename Sacado::ValueType< value_type >::type;
  unsigned fad_dim;
  bool is_view_type;

  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp()
  : fad_dim(0) , is_view_type(false) {}

  // Assume all views are View or DynRankView
  // TODO If assumption is insufficient, better deduction on is_view...
  template < class View >
  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp( const View & view )
  : fad_dim ( dimension_scalar(view) )
  {
    is_view_type = (Kokkos::is_view<View>::value || Kokkos::is_view_fad<View>::value);
  }

  // TODO If assumption is insufficient, better deduction on is_view...
  template < class View, class ... Views >
  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp( const View & view,  const Views & ... views ) 
  : fad_dim ( dimension_scalar(view, views... ) )
  {
    is_view_type = (Kokkos::is_view<View>::value || Kokkos::is_view_fad<View>::value);
  }

};

template < class ValueType >
struct CommonViewAllocProp< Kokkos::Impl::ViewSpecializeSacadoFadContiguous, ValueType >
{
  using value_type = ValueType;
  using scalar_array_type = typename Sacado::ValueType< value_type >::type;
  unsigned fad_dim;
  bool is_view_type;

  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp()
  : fad_dim(0) , is_view_type(false) {}

  // Assume all views are View or DynRankView
  // TODO If assumption is insufficient, better deduction on is_view...
  template < class View >
  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp( const View & view )
  : fad_dim ( dimension_scalar(view) )
  {
    is_view_type = (Kokkos::is_view<View>::value || Kokkos::is_view_fad<View>::value);
  }

  // TODO If assumption is insufficient, better deduction on is_view...
  template < class View, class ... Views >
  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp( const View & view,  const Views & ... views ) 
  : fad_dim ( dimension_scalar(view, views... ) )
  {
    is_view_type = (Kokkos::is_view<View>::value || Kokkos::is_view_fad<View>::value);
  }
};

// Detect if a ViewCtorProp contains a CommonViewAllocProp
template < typename ... P >
struct has_common_view_alloc_prop : public std::false_type {};

template < class Specialize, class ValueType >
struct has_common_view_alloc_prop< CommonViewAllocProp<Specialize, ValueType> > : public std::true_type {};


// Check for CommonViewAllocProp in pack of properties
template < typename ... >
struct check_has_common_view_alloc_prop;

template <>
struct check_has_common_view_alloc_prop<>
{
  enum { value = false };
};

template < typename P >
struct check_has_common_view_alloc_prop<P>
{
  enum { value = has_common_view_alloc_prop< P >::value };
};

template < typename P0, typename ... P >
struct check_has_common_view_alloc_prop<P0, P...>
{
  enum { value = ( (has_common_view_alloc_prop<P0>::value == true) ? true : check_has_common_view_alloc_prop<P...>::value ) };
};

template < typename ... >
struct compute_fad_dim_from_alloc_prop;

template < >
struct compute_fad_dim_from_alloc_prop<> {
  template <typename CtorProp>
  KOKKOS_INLINE_FUNCTION
  static unsigned eval(const CtorProp&) { return 0; }
};

template < typename P >
struct compute_fad_dim_from_alloc_prop<P> {
  template <typename CtorProp>
  KOKKOS_INLINE_FUNCTION
  static unsigned eval(const CtorProp&) { return 0; }
};

template < typename P0, typename ... P >
struct compute_fad_dim_from_alloc_prop<P0,P...> {
  template <typename CtorProp>
  KOKKOS_INLINE_FUNCTION
  static unsigned eval(const CtorProp& prop) {
    unsigned d1 = compute_fad_dim_from_alloc_prop<P0>::eval(prop);
    unsigned d2 = compute_fad_dim_from_alloc_prop<P...>::eval(prop);
    return d1 > d2 ? d1 : d2;
  }
};

template < class ValueType >
struct compute_fad_dim_from_alloc_prop<
  CommonViewAllocProp<ViewSpecializeSacadoFad, ValueType>
  > {
  template <typename CtorProp>
  KOKKOS_INLINE_FUNCTION
  static unsigned eval(const CtorProp& prop) {
    using specialize = ViewSpecializeSacadoFad;
    using CVAP = CommonViewAllocProp< specialize, ValueType >;
    auto cast_prop = ((Kokkos::Impl::ViewCtorProp<void, CVAP> const &)prop).value;
    return cast_prop.fad_dim;
  }
};

template < class ValueType >
struct compute_fad_dim_from_alloc_prop<
  CommonViewAllocProp<ViewSpecializeSacadoFadContiguous, ValueType>
  > {
  template <typename CtorProp>
  KOKKOS_INLINE_FUNCTION
  static unsigned eval(const CtorProp& prop) {
    using specialize = ViewSpecializeSacadoFadContiguous;
    using CVAP = CommonViewAllocProp< specialize, ValueType >;
    auto cast_prop = ((Kokkos::Impl::ViewCtorProp<void, CVAP> const &)prop).value;
    return cast_prop.fad_dim;
  }
};

template <typename Traits, typename ... P >
struct appendFadToLayoutViewAllocHelper
{
  using layout_type = typename Traits::array_layout;
  using specialize = typename Traits::specialize;
  using CtorProp = ViewCtorProp< P... >;

  KOKKOS_INLINE_FUNCTION
  static layout_type returnNewLayoutPlusFad( const CtorProp & arg_prop, const layout_type & arg_layout ) {

    layout_type appended_layout( arg_layout );

    // Static View case - DynRankView layout handled within createLayout calls

    const unsigned fad_dim =
      compute_fad_dim_from_alloc_prop<P...>::eval(arg_prop);
    appended_layout.dimension[ Traits::rank ] = (fad_dim > 0) ? fad_dim : 1;

    return appended_layout;
  }
};

template <typename Layout>
struct prependFadToLayout
{
  using layout_type = Layout;

  template < typename FadSizeType >
  KOKKOS_INLINE_FUNCTION
  static layout_type returnNewLayoutPlusFad( const layout_type & arg_layout, const FadSizeType fad_dim ) {

    layout_type prepended_layout(0,0,0,0,0,0,0,0);

    prepended_layout.dimension[0] = fad_dim;

    for ( int i = 1; i < ARRAY_LAYOUT_MAX_RANK; ++i ) {
      prepended_layout.dimension[i] = arg_layout.dimension[i-1];
    }

    return prepended_layout;
  }
};

} } // namespace Kokkos::Impl


//----------------------------------------------------------------------------


namespace Kokkos {
namespace Impl {

template< class Traits >
class ViewMapping< Traits , /* View internal mapping */
  typename std::enable_if<
    ( std::is_same< typename Traits::specialize
                  , ViewSpecializeSacadoFad >::value
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

  typedef typename Traits::value_type  fad_type ;
  typedef typename Sacado::ValueType< fad_type >::type fad_value_type ;
  typedef typename
    std::add_const< fad_value_type >::type  const_fad_value_type ;

  enum { FadStaticDimension = Sacado::StaticSize< fad_type >::value };
  typedef Sacado::integral_nonzero< unsigned , FadStaticDimension > sacado_size_type;

  // Only LayoutRight has a static stride one
  enum { FadStaticStride =
    std::is_same< typename Traits::array_layout
                , Kokkos::LayoutRight >::value ? 1 : 0 };

  typedef Sacado::integral_nonzero< unsigned , FadStaticStride > sacado_stride_type;

  typedef fad_value_type * handle_type ;

  typedef ViewArrayAnalysis< typename Traits::data_type > array_analysis ;

  // Offset without Fad dimension
  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  // Append the fad dimension for the internal offset mapping.
  typedef ViewOffset
    < typename array_analysis::dimension::
        template append<( unsigned(FadStaticDimension) > 0 ? unsigned(FadStaticDimension) + 1 : 0 )>::type
    , typename Traits::array_layout
    , void
    >  array_offset_type ;

  handle_type  m_impl_handle ;
  offset_type  m_impl_offset ;
  array_offset_type  m_array_offset ;
  sacado_size_type m_fad_size ;
  sacado_stride_type m_fad_stride ;

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  // Using the internal offset mapping so limit to public rank:
  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr size_t extent( const iType & r ) const
    { return m_impl_offset.m_dim.extent(r) ; }

  KOKKOS_INLINE_FUNCTION constexpr
  typename Traits::array_layout layout() const
    { return m_impl_offset.layout(); }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const
    { return m_impl_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const
    { return m_impl_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const
    { return m_impl_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const
    { return m_impl_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const
    { return m_impl_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const
    { return m_impl_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const
    { return m_impl_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const
    { return m_impl_offset.dimension_7(); }

  // Can only be regular layout with uniform striding
  // when LayoutRight with contiguous values so not guaranteed true.
  using is_regular = std::false_type ;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const
    { return m_impl_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const
    { return m_impl_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const
    { return m_impl_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const
    { return m_impl_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const
    { return m_impl_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const
    { return m_impl_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const
    { return m_impl_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const
    { return m_impl_offset.stride_7(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const
    { m_impl_offset.stride(s) ; }

  // Size of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned dimension_scalar() const
    { return m_fad_size.value+1; }

  // trode of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned stride_scalar() const
    { return m_fad_stride.value; }

  //----------------------------------------
  // Range of mapping

  // Return type of reference operators
  typedef typename
    Sacado::ViewFadType< fad_type , FadStaticDimension , FadStaticStride >::type  reference_type ;

  /** \brief Pointer to underlying memory type */
  typedef fad_value_type * pointer_type ;

  /** \brief  Span of the mapped range : [ data() .. data() + span() ) */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const
    { return m_array_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_array_offset.span_is_contiguous() ; }

  /** \brief Raw data access */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const
    { return m_impl_handle ; }

  //----------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const
    { return reference_type( m_impl_handle
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type
  reference( const I0 & i0 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }


  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,i2,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,i2,i3,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,i5,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }


  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return reference_type( m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,i5,i6,0)
                           , m_fad_size.value
                           , m_fad_stride.value ); }

  //----------------------------------------

  /** \brief  Span, in bytes, of the required memory */
  KOKKOS_INLINE_FUNCTION
  static size_t memory_span( typename Traits::array_layout const & layout )
    {
      size_t dims[8];
      for (int i=0; i<8; ++i)
        dims[i] = layout.dimension[i];
      if (unsigned(FadStaticDimension) > 0)
        dims[unsigned(Rank)] = FadStaticDimension+1;

      typename Traits::array_layout alayout(
        dims[0], dims[1], dims[2], dims[3],
        dims[4], dims[5], dims[6], dims[7] );

      // Do not introduce padding...
      typedef std::integral_constant< unsigned , 0 >  padding ;
      return array_offset_type( padding() , alayout ).span() * sizeof(fad_value_type);
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION ~ViewMapping() {}
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_impl_handle(0) , m_impl_offset() , m_array_offset() , m_fad_size(0) , m_fad_stride(0) {}

  KOKKOS_INLINE_FUNCTION ViewMapping( const ViewMapping & ) = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( const ViewMapping & ) = default ;

  KOKKOS_INLINE_FUNCTION ViewMapping( ViewMapping && ) = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( ViewMapping && ) = default ;

  template< class ... P >
  KOKKOS_INLINE_FUNCTION
  ViewMapping
    ( ViewCtorProp< P ... > const & prop
    , typename Traits::array_layout const & local_layout
    )
    : m_impl_handle( ( (ViewCtorProp<void,pointer_type> const &) prop ).value )
    , m_impl_offset( std::integral_constant< unsigned , 0 >()
              , local_layout )
    , m_array_offset( std::integral_constant< unsigned , 0 >()
                    , local_layout )
    // Query m_array_offset, not input, in case of static dimension
    , m_fad_size(
       ( Rank == 0 ? m_array_offset.dimension_0() :
       ( Rank == 1 ? m_array_offset.dimension_1() :
       ( Rank == 2 ? m_array_offset.dimension_2() :
       ( Rank == 3 ? m_array_offset.dimension_3() :
       ( Rank == 4 ? m_array_offset.dimension_4() :
       ( Rank == 5 ? m_array_offset.dimension_5() :
       ( Rank == 6 ? m_array_offset.dimension_6() :
                     m_array_offset.dimension_7() ))))))) - 1 )
    , m_fad_stride(
       ( Rank == 0 ? m_array_offset.stride_0() :
       ( Rank == 1 ? m_array_offset.stride_1() :
       ( Rank == 2 ? m_array_offset.stride_2() :
       ( Rank == 3 ? m_array_offset.stride_3() :
       ( Rank == 4 ? m_array_offset.stride_4() :
       ( Rank == 5 ? m_array_offset.stride_5() :
       ( Rank == 6 ? m_array_offset.stride_6() :
                     m_array_offset.stride_7() ))))))))

    {
      const unsigned fad_dim =
       ( Rank == 0 ? m_array_offset.dimension_0() :
       ( Rank == 1 ? m_array_offset.dimension_1() :
       ( Rank == 2 ? m_array_offset.dimension_2() :
       ( Rank == 3 ? m_array_offset.dimension_3() :
       ( Rank == 4 ? m_array_offset.dimension_4() :
       ( Rank == 5 ? m_array_offset.dimension_5() :
       ( Rank == 6 ? m_array_offset.dimension_6() :
         m_array_offset.dimension_7() )))))));
      if (unsigned(FadStaticDimension) == 0 && fad_dim == 0)
        Kokkos::abort("invalid fad dimension (0) supplied!");
    }

  //----------------------------------------
  /*  Allocate and construct mapped array.
   *  Allocate via shared allocation record and
   *  return that record for allocation tracking.
   */
  template< class ... P >
  SharedAllocationRecord<> *
  allocate_shared( ViewCtorProp< P... > const & prop
                 , typename Traits::array_layout const & local_layout )
  {
    typedef ViewCtorProp< P... > ctor_prop ;

    typedef typename ctor_prop::execution_space  execution_space ;
    typedef typename Traits::memory_space         memory_space ;
    typedef ViewValueFunctor< execution_space , fad_value_type > functor_type ;
    typedef SharedAllocationRecord< memory_space , functor_type > record_type ;

    // Disallow padding
    typedef std::integral_constant< unsigned , 0 > padding ;

    // Check if ViewCtorProp has CommonViewAllocProp - if so, retrieve the fad_size and append to layout
    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };

    m_impl_offset = offset_type( padding(), local_layout );

    typename Traits::array_layout internal_layout =
      (test_traits_check == true)
      ? Kokkos::Impl::appendFadToLayoutViewAllocHelper< Traits, P... >::returnNewLayoutPlusFad(prop, local_layout)
      : local_layout;

    m_array_offset = array_offset_type( padding(), internal_layout );

    const unsigned fad_dim =
      ( Rank == 0 ? m_array_offset.dimension_0() :
      ( Rank == 1 ? m_array_offset.dimension_1() :
      ( Rank == 2 ? m_array_offset.dimension_2() :
      ( Rank == 3 ? m_array_offset.dimension_3() :
      ( Rank == 4 ? m_array_offset.dimension_4() :
      ( Rank == 5 ? m_array_offset.dimension_5() :
      ( Rank == 6 ? m_array_offset.dimension_6() :
        m_array_offset.dimension_7() )))))));
    if (unsigned(FadStaticDimension) == 0 && fad_dim == 0)
      Kokkos::abort("invalid fad dimension (0) supplied!");
    m_fad_size = fad_dim - 1 ;

    m_fad_stride =
       ( Rank == 0 ? m_array_offset.stride_0() :
       ( Rank == 1 ? m_array_offset.stride_1() :
       ( Rank == 2 ? m_array_offset.stride_2() :
       ( Rank == 3 ? m_array_offset.stride_3() :
       ( Rank == 4 ? m_array_offset.stride_4() :
       ( Rank == 5 ? m_array_offset.stride_5() :
       ( Rank == 6 ? m_array_offset.stride_6() :
                     m_array_offset.stride_7() )))))));

    const size_t alloc_size = m_array_offset.span() * sizeof(fad_value_type);

    // Create shared memory tracking record with allocate memory from the memory space
    record_type * const record =
      record_type::allocate( ( (ViewCtorProp<void,memory_space> const &) prop ).value
                           , ( (ViewCtorProp<void,std::string>  const &) prop ).value
                           , alloc_size );

    //  Only set the the pointer and initialize if the allocation is non-zero.
    //  May be zero if one of the dimensions is zero.
    if ( alloc_size ) {

      m_impl_handle = handle_type( reinterpret_cast< pointer_type >( record->data() ) );

      if ( ctor_prop::initialize ) {
        // Assume destruction is only required when construction is requested.
        // The ViewValueFunctor has both value construction and destruction operators.
        record->m_destroy = functor_type( ( (ViewCtorProp<void,execution_space> const &) prop).value
                                        , (fad_value_type *) m_impl_handle
                                        , m_array_offset.span()
                                        );

        // Construct values
        record->m_destroy.construct_shared_allocation();
      }
    }

    return record ;
  }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<FAD>      = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has FAD
    std::is_same< typename DstTraits::specialize
                , ViewSpecializeSacadoFad >::value
    &&
    // Source view has FAD
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFad >::value
  )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcFadType ;

  template< class DstType >
  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcFadType & src
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
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

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
        ViewDimensionAssignable
          < typename DstType::offset_type::dimension_type
          , typename SrcFadType::offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

       static_assert(
        ViewDimensionAssignable
          < typename DstType::array_offset_type::dimension_type
          , typename SrcFadType::array_offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::array_offset_type  dst_array_offset_type ;

      dst.m_impl_handle  = src.m_impl_handle ;
      dst.m_impl_offset  = dst_offset_type( src.m_impl_offset );
      dst.m_array_offset = dst_array_offset_type( src.m_array_offset );
      dst.m_fad_size = src.m_fad_size.value ;
      dst.m_fad_stride = src.m_fad_stride.value ;
    }
};

// Integer argument is the actual rank => ranks 0 to Rank-1 will be assigned
/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<ordinary> = View<FAD>
 *  (TBD)  View<FAD> = View<ordinary>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFad >::value
  )>::type >
{
public:

  enum { is_assignable = true };


  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcFadType ;


  // Helpers to assign, and generate if necessary, ViewOffset to the dst map
  // These are necessary to use Kokkos' deep_copy with nested fads
  template < class DstType, class SrcFadType, class Truth = void >
    struct AssignOffset;

  template < class DstType, class SrcFadType >
    struct AssignOffset< DstType, SrcFadType, typename std::enable_if< ((int)DstType::offset_type::dimension_type::rank != (int)SrcFadType::array_offset_type::dimension_type::rank) >::type >
    {
      // ViewOffset's Dimensions Ranks do not match
      KOKKOS_INLINE_FUNCTION
      static void assign( DstType & dst, const SrcFadType & src )
      {
        typedef typename SrcTraits::value_type TraitsValueType;

        if ( Sacado::IsFad<TraitsValueType>::value
            && Sacado::IsStaticallySized< typename Sacado::ValueType< TraitsValueType >::type >::value
           )
        {
          typedef typename DstType::offset_type::array_layout DstLayoutType;
          //typedef typename ViewArrayLayoutSelector<typename DstType::offset_type::array_layout>::type DstLayoutType;
          typedef typename SrcFadType::array_offset_type::dimension_type SrcViewDimension;

          // This is the static dimension of the inner fad, missing from ViewDimension
          const size_t InnerStaticDim = Sacado::StaticSize< typename Sacado::ValueType< TraitsValueType >::type >::value;

          static constexpr bool is_layout_left =
            std::is_same< DstLayoutType, Kokkos::LayoutLeft>::value;

          typedef typename std::conditional< is_layout_left,
                                             typename SrcViewDimension:: template prepend< InnerStaticDim+1 >::type,
                                             typename SrcViewDimension:: template append < InnerStaticDim+1 >::type
                    >::type SrcViewDimensionAppended;

          typedef std::integral_constant< unsigned , 0 >  padding ;

          typedef ViewOffset< SrcViewDimensionAppended, DstLayoutType > TmpOffsetType;

          auto src_layout = src.m_array_offset.layout();

          if ( is_layout_left ) {
            auto prepend_layout = Kokkos::Impl::prependFadToLayout< DstLayoutType >::returnNewLayoutPlusFad(src_layout, InnerStaticDim+1);
            TmpOffsetType offset_tmp( padding(), prepend_layout );
            dst.m_impl_offset = offset_tmp;
          }
          else {
            TmpOffsetType offset_tmp( padding(), src_layout );
            dst.m_impl_offset = offset_tmp;
          }

        } else {
          Kokkos::abort("Sacado error: Applying AssignOffset for case with nested Fads, but without nested Fads - something went wrong");
        }
      }
    };

  template < class DstType, class SrcFadType >
    struct AssignOffset< DstType, SrcFadType, typename std::enable_if< ((int)DstType::offset_type::dimension_type::rank == (int)SrcFadType::array_offset_type::dimension_type::rank) >::type >
    {
      KOKKOS_INLINE_FUNCTION
      static void assign( DstType & dst, const SrcFadType & src )
      {

        typedef typename DstType::offset_type  dst_offset_type ;
        dst.m_impl_offset  = dst_offset_type( src.m_array_offset );
      }
    };


// If the dst and src mappings are not equal in Rank, the src should come from a View of nested fads
// In the case of two nested fads, the innermost must be an SFad (static Fad)
// The offset_type's are not compatible in the case of nested fads because the ViewDimension's ranks will not agree
// In this case, rather than trying to construct an offset_type from src (which will fail at compile time)
// and assign to dst.m_impl_offset, manually assign the ViewDimension arguments to dst;
// requires appending the missing inner SFad dim + 1 to the Rank-1 ViewDimension
  // DstType and SrcFadType are MAPS...
  template < class DstType >
  KOKKOS_INLINE_FUNCTION static
  void
  assign( DstType & dst
        , const SrcFadType & src
        , const TrackType &
        )
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
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

    static_assert(
        std::is_same< typename DstTraits::array_layout
        , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
        , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );
#if 0
    static_assert(
        std::is_same< typename DstTraits::scalar_array_type
        , typename SrcTraits::scalar_array_type >::value ||
        std::is_same< typename DstTraits::scalar_array_type
        , typename SrcTraits::const_scalar_array_type >::value ,
        "View assignment must have same value type or const = non-const" );
#endif

    AssignOffset< DstType, SrcFadType >::assign( dst, src );

    dst.m_impl_handle  = reinterpret_cast< typename DstType::handle_type >(src.m_impl_handle) ;
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Subview mapping

template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      // Source view has FAD only
      std::is_same< typename SrcTraits::specialize
                  , ViewSpecializeSacadoFad >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    )>::type
  , SrcTraits
  , Args ... >
{
private:

  static_assert( SrcTraits::rank == sizeof...(Args) , "" );

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Args...>::value)
    , R2 = bool(is_integral_extent<2,Args...>::value)
    , R3 = bool(is_integral_extent<3,Args...>::value)
    , R4 = bool(is_integral_extent<4,Args...>::value)
    , R5 = bool(is_integral_extent<5,Args...>::value)
    , R6 = bool(is_integral_extent<6,Args...>::value)
    };

  // Public rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  // Whether right-most non-FAD rank is a range.
  enum { R0_rev = ( 0 == SrcTraits::rank ? RZ : (
                    1 == SrcTraits::rank ? R0 : (
                    2 == SrcTraits::rank ? R1 : (
                    3 == SrcTraits::rank ? R2 : (
                    4 == SrcTraits::rank ? R3 : (
                    5 == SrcTraits::rank ? R4 : (
                    6 == SrcTraits::rank ? R5 : R6 ))))))) };

  // Subview's layout
  // If LayoutRight then FAD is contiguous
  // For LayoutLeft, result is LayoutLeft only if 1st arg is a range,
  // and since last (FAD) dimension is also a range, and these
  // ranges must be consecutive, the input rank must be 1
  typedef typename std::conditional<
    ( /* Same layout IF */
      ( rank == 0 )
      ||
      ( std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value
        &&
        ( rank == 1 ) && R0_rev
      )
      ||
      ( std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value
        &&
        ( rank == 1 ) && (SrcTraits::rank == 1) && R0
      )
    ), typename SrcTraits::array_layout , Kokkos::LayoutStride
    >::type  array_layout ;

  typedef typename SrcTraits::value_type  fad_type ;

  typedef typename std::conditional< rank == 0 , fad_type ,
          typename std::conditional< rank == 1 , fad_type * ,
          typename std::conditional< rank == 2 , fad_type ** ,
          typename std::conditional< rank == 3 , fad_type *** ,
          typename std::conditional< rank == 4 , fad_type **** ,
          typename std::conditional< rank == 5 , fad_type ***** ,
          typename std::conditional< rank == 6 , fad_type ****** ,
                                                 fad_type *******
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


  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< traits_type , void > & dst
                    , ViewMapping< SrcTraits , void > const & src
                    , Args ... args )
    {
      typedef ViewMapping< traits_type , void > DstType ;
      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::array_offset_type  dst_array_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      const SubviewExtents< SrcTraits::rank , rank >
        extents( src.m_impl_offset.m_dim , args... );
      const SubviewExtents< SrcTraits::rank + 1 , rank + 1 >
        array_extents( src.m_array_offset.m_dim , args... , Kokkos::ALL() );

      dst.m_impl_offset = dst_offset_type( src.m_impl_offset , extents );
      dst.m_array_offset = dst_array_offset_type( src.m_array_offset , array_extents );
      dst.m_impl_handle =
        dst_handle_type( src.m_impl_handle +
                         src.m_array_offset( array_extents.domain_offset(0)
                                           , array_extents.domain_offset(1)
                                           , array_extents.domain_offset(2)
                                           , array_extents.domain_offset(3)
                                           , array_extents.domain_offset(4)
                                           , array_extents.domain_offset(5)
                                           , array_extents.domain_offset(6)
                                           , array_extents.domain_offset(7) ) );
      dst.m_fad_size = src.m_fad_size;
      dst.m_fad_stride = src.m_fad_stride.value;
    }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined(HAVE_SACADO_KOKKOSCORE) && \
    defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM) && \
    defined(HAVE_SACADO_VIEW_SPEC) && \
    ! defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_TeuchosCommAdapters.hpp"

namespace Teuchos {

template< typename Ordinal , class SD , class ... SP , class RD , class ... RP >
typename std::enable_if<Kokkos::is_view_fad< Kokkos::View<SD,SP...> >::value &&
                        Kokkos::is_view_fad< Kokkos::View<RD,RP...> >::value
                       >::type
reduceAll
  ( const Comm<Ordinal>& comm,
    const EReductionType reductType ,
    const Ordinal count,
    const Kokkos::View<SD,SP...> & sendBuffer ,
    const Kokkos::View<RD,RP...> & recvBuffer )
{
  // We can't implement reduceAll by extracting the underlying array (since we
  // can't reduce across the derivative dimension) and we can't just extract
  // a pointer due to ViewFad.  In principle we could handle ViewFad in the
  // serializer, but for the time being we just copy the view's into local
  // buffers (on the host).
  typedef Kokkos::View<SD,SP...> SendViewType;
  typedef Kokkos::View<RD,RP...> RecvViewType;
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
    SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
    "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
    "The send View's rank is " << SendViewType::rank << " and the receive "
    "View's rank is " << RecvViewType::rank << ".");

  // Copy send buffer into local array
  Teuchos::Array<send_value_type> localSendBuffer(count);
  typename SendViewType::HostMirror hostSendBuffer =
    Kokkos::create_mirror_view(sendBuffer);
  Kokkos::deep_copy(hostSendBuffer, sendBuffer);
  for (Ordinal i=0; i<count; ++i)
    localSendBuffer[i] = hostSendBuffer(i);

  // Copy receive buffer into local array (necessary to initialize Fad types
  // properly)
  Teuchos::Array<recv_value_type> localRecvBuffer(count);
  typename RecvViewType::HostMirror hostRecvBuffer =
    Kokkos::create_mirror_view(recvBuffer);
  Kokkos::deep_copy(hostRecvBuffer, recvBuffer);
  for (Ordinal i=0; i<count; ++i)
    localRecvBuffer[i] = hostRecvBuffer(i);

  // Do reduce-all
  reduceAll(comm, reductType, count,
            localSendBuffer.getRawPtr(),
            localRecvBuffer.getRawPtr());

  // Copy back into original buffer
  for (Ordinal i=0; i<count; ++i)
    hostRecvBuffer(i) = localRecvBuffer[i];
  Kokkos::deep_copy(recvBuffer, hostRecvBuffer);
}


template< typename Ordinal , typename Serializer ,
          class SD , class ... SP , class RD , class ... RP >
typename std::enable_if<Kokkos::is_view_fad< Kokkos::View<SD,SP...> >::value &&
                        Kokkos::is_view_fad< Kokkos::View<RD,RP...> >::value
                       >::type
reduceAll
  ( const Comm<Ordinal>& comm,
    const Serializer& serializer,
    const EReductionType reductType ,
    const Ordinal count,
    const Kokkos::View<SD,SP...> & sendBuffer ,
    const Kokkos::View<RD,RP...> & recvBuffer )
{
  // We can't implement reduceAll by extracting the underlying array (since we
  // can't reduce across the derivative dimension) and we can't just extract
  // a pointer due to ViewFad.  In principle we could handle ViewFad in the
  // serializer, but for the time being we just copy the view's into local
  // buffers (on the host).
  typedef Kokkos::View<SD,SP...> SendViewType;
  typedef Kokkos::View<RD,RP...> RecvViewType;
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
    SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
    "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
    "The send View's rank is " << SendViewType::rank << " and the receive "    "View's rank is " << RecvViewType::rank << ".");

  // Copy send buffer into local array
  Teuchos::Array<send_value_type> localSendBuffer(count);
  typename SendViewType::HostMirror hostSendBuffer =
    Kokkos::create_mirror_view(sendBuffer);
  Kokkos::deep_copy(hostSendBuffer, sendBuffer);
  for (Ordinal i=0; i<count; ++i)
    localSendBuffer[i] = hostSendBuffer(i);

  // Copy receive buffer into local array (necessary to initialize Fad types
  // properly)
  Teuchos::Array<recv_value_type> localRecvBuffer(count);
  typename RecvViewType::HostMirror hostRecvBuffer =
    Kokkos::create_mirror_view(recvBuffer);
  Kokkos::deep_copy(hostRecvBuffer, recvBuffer);
  for (Ordinal i=0; i<count; ++i)
    localRecvBuffer[i] = hostRecvBuffer(i);

  // Do reduce-all
  reduceAll(comm, serializer, reductType, count,
            localSendBuffer.getRawPtr(),
            localRecvBuffer.getRawPtr());

  // Copy back into original buffer
  for (Ordinal i=0; i<count; ++i)
    hostRecvBuffer(i) = localRecvBuffer[i];
  Kokkos::deep_copy(recvBuffer, hostRecvBuffer);
}


template<typename Ordinal, class D, class ... P  >
typename std::enable_if<Kokkos::is_view_fad< Kokkos::View<D,P...> >::value>::type
broadcast
  ( const Comm<Ordinal>& comm,
    const int rootRank ,
    const Ordinal count,
    const Kokkos::View<D,P...>& buffer)
{
  typedef Kokkos::View<D,P...> view_type;
  typename view_type::array_type array_buffer = buffer;
  Ordinal array_count = count * Kokkos::dimension_scalar(buffer);
  broadcast( comm, rootRank, array_count, array_buffer );
}

template<typename Ordinal,
         typename Serializer ,
         class D, class ... P >
typename std::enable_if<Kokkos::is_view_fad< Kokkos::View<D,P...> >::value>::type
broadcast
  ( const Comm<Ordinal>& comm,
    const Serializer& serializer,
    const int rootRank ,
    const Ordinal count,
    const Kokkos::View<D,P...>& buffer)
{
  typedef Kokkos::View<D,P...> view_type;
  typename view_type::array_type array_buffer = buffer;
  Ordinal array_count = count * Kokkos::dimension_scalar(buffer);
  broadcast( comm, *(serializer.getValueSerializer()), rootRank,
             array_count, array_buffer );
}

} // namespace Teuchos

#endif

//----------------------------------------------------------------------------

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOSCORE)

#include "KokkosExp_View_Fad_Contiguous.hpp"

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_HPP */
