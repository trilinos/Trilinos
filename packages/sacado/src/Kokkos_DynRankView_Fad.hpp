// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP
#define KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP

#include "Sacado_ConfigDefs.h"

// This file is setup to always work even when KokkosContainers (which contains
// Kokkos::DynRankView) isn't enabled.

#if defined(HAVE_SACADO_KOKKOS)

// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Layout.hpp"
//#include "Kokkos_DynRankView.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

namespace Kokkos {

template< class DataType , class ... Properties >
struct ViewTraits ;

template< class DataType , class ... Properties >
class DynRankView ;

namespace Impl {

template<class Space, class T, class ... P>
struct MirrorDRViewType;

}

template <typename view_type>
struct is_dynrankview_fad { static const bool value = false; };

template <typename view_type>
struct is_dynrankview_fad_contiguous { static const bool value = false; };

}

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_View_Fad.hpp"

namespace Kokkos {

// Declare overloads of create_mirror() so they are in scope
// Kokkos_DynRankView.hpp is included below

template< class T , class ... P >
inline
typename Kokkos::DynRankView<T,P...>::HostMirror
create_mirror(
  const Kokkos::DynRankView<T,P...> & src,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    !std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
        Kokkos::LayoutStride >::value >::type * = 0);


template< class T , class ... P >
inline
typename Kokkos::DynRankView<T,P...>::HostMirror
create_mirror(
  const Kokkos::DynRankView<T,P...> & src,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,P...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value >::type * = 0);

template<class Space, class T, class ... P>
typename Impl::MirrorDRViewType<Space,T,P ...>::dest_view_type
create_mirror(
  const Space&,
  const Kokkos::DynRankView<T,P...> & src,
  typename std::enable_if<
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value >::type * = 0);

namespace Impl {

template <unsigned N, typename T, typename... Args>
KOKKOS_FUNCTION auto as_view_of_rank_n(
  DynRankView<T, Args...> v,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,Args...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,Args...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
    >::type * = 0 );

}

} // namespace Kokkos

#include "Kokkos_DynRankView.hpp"

namespace Kokkos {
namespace Impl {

template <>
struct DynRankDimTraits<Kokkos::Impl::ViewSpecializeSacadoFad> {

  enum : size_t{unspecified = ~size_t(0)};

  // Compute the rank of the view from the nonzero dimension arguments.
  // For views of Fad, the rank is one less than the rank determined by the nonzero dimension args
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const size_t N0
                           , const size_t N1
                           , const size_t N2
                           , const size_t N3
                           , const size_t N4
                           , const size_t N5
                           , const size_t N6
                           , const size_t N7 )
  {
    return
      (   (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified && N0 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified) ? 1
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified) ? 2
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified) ? 3
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified) ? 4
      : ( (N7 == unspecified && N6 == unspecified) ? 5
      : ( (N7 == unspecified) ? 6
      : 7 ) ) ) ) ) ) ) );
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const Layout& layout )
  {
    return computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );
  }

  // Compute the rank of the view from the nonzero layout arguments and possible hidden dim
  template <typename Layout, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const ViewCtorProp<P...>& /* prop */, const Layout& layout )
  {
    size_t rank = computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );

    // Check if has_common_view_alloc_prop; if so, return rank+1, else rank
    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    return (test_traits_check == true) ? rank+1 : rank;
  }

  // Create the layout for the rank-7 view.
  // For Fad we have to move the fad dimension to the last (rank 8 since the DynRankView is rank-7)
  // LayoutLeft or LayoutRight
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutRight>::value || std::is_same<Layout , Kokkos::LayoutLeft>::value) , Layout >::type createLayout( const Layout& layout )
  {
    Layout l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
            , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
            , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
            , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
            , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
            , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
            , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
            , layout.dimension[7] != unspecified ? layout.dimension[7] : 1 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  //LayoutStride
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutStride>::value) , Layout>::type createLayout( const Layout& layout )
  {
    Layout      l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
                 , layout.stride[0]
                 , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
                 , layout.stride[1]
                 , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
                 , layout.stride[2]
                 , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
                 , layout.stride[3]
                 , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
                 , layout.stride[4]
                 , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
                 , layout.stride[5]
                 , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
                 , layout.stride[6]
                 , layout.dimension[7] != unspecified ? layout.dimension[7] : 1
                 , layout.stride[7]
                 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;
    const size_t fad_stride = layout.stride[fad_dim];
    l.stride[fad_dim] = 0;
    l.stride[7] = fad_stride;

    return l;
  }

  // If fad_dim is stored in ViewCtorProp
  // LayoutLeft or LayoutRight
  template <typename Traits, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<typename Traits::array_layout , Kokkos::LayoutRight>::value || std::is_same<typename Traits::array_layout , Kokkos::LayoutLeft>::value) , typename Traits::array_layout >::type createLayout( const ViewCtorProp<P...> & arg_prop, const typename Traits::array_layout& layout )
  {
    using Layout = typename Traits::array_layout;

    Layout l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
            , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
            , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
            , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
            , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
            , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
            , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
            , layout.dimension[7] != unspecified ? layout.dimension[7] : 1 );

    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    if (test_traits_check == true) {
      l.dimension[7] = compute_fad_dim_from_alloc_prop<P...>::eval(arg_prop);
    }
    else {
      const unsigned fad_dim = computeRank(layout);
      const size_t fad_size = layout.dimension[fad_dim];
      l.dimension[fad_dim] = 1;
      l.dimension[7] = fad_size;
    }

    return l;
  }

  // If fad_dim is stored in ViewCtorProp
  //LayoutStride
  template <typename Traits, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<typename Traits::array_layout , Kokkos::LayoutStride>::value) , typename Traits::array_layout>::type createLayout( const ViewCtorProp<P...> & arg_prop, const typename Traits::array_layout& layout )
  {
    using Layout = typename Traits::array_layout;

    Layout      l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
                 , layout.stride[0]
                 , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
                 , layout.stride[1]
                 , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
                 , layout.stride[2]
                 , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
                 , layout.stride[3]
                 , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
                 , layout.stride[4]
                 , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
                 , layout.stride[5]
                 , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
                 , layout.stride[6]
                 , layout.dimension[7] != unspecified ? layout.dimension[7] : 1
                 , layout.stride[7]
                 );

    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    const unsigned fad_dim = computeRank(layout);
    if (test_traits_check == true) {
      l.dimension[fad_dim] = 1;
      l.dimension[7] = compute_fad_dim_from_alloc_prop<P...>::eval(arg_prop);
    }
    else {
      const size_t fad_size = layout.dimension[fad_dim];
      l.dimension[fad_dim] = 1;
      l.dimension[7] = fad_size;
    }
    const size_t fad_stride = layout.stride[fad_dim];
    l.stride[fad_dim] = 0;
    l.stride[7] = fad_stride;

    return l;
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  template <typename ViewType, typename ViewArg>
  static ViewType createView( const ViewArg& arg
                            , const size_t N0
                            , const size_t N1
                            , const size_t N2
                            , const size_t N3
                            , const size_t N4
                            , const size_t N5
                            , const size_t N6
                            , const size_t N7 )
  {
    typename ViewType::array_layout l( N0, N1, N2, N3, N4, N5, N6, N7 );
    typename ViewType::array_layout l_fad = createLayout(l);
    return ViewType( arg
                   , l_fad.dimension[0]
                   , l_fad.dimension[1]
                   , l_fad.dimension[2]
                   , l_fad.dimension[3]
                   , l_fad.dimension[4]
                   , l_fad.dimension[5]
                   , l_fad.dimension[6]
                   , l_fad.dimension[7] );
  }

};

}} // namespace Kokkos::Impl

namespace Kokkos {
namespace Impl {

template <unsigned> struct AssignDim7 {
  template <typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const size_t src_dim) {}
};
template <> struct AssignDim7<0u> {
  template <typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const size_t src_dim) {
    dst.N7 = src_dim;
  }
};

// Utility class that handles calculation of the stride in Fad subview
template <unsigned,unsigned> struct AssignFadDimStride {};
template <unsigned StaticDim> struct AssignFadDimStride<0u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = 0 ;
    dst.m_stride.S1 = 0 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S0 ;

    dst.m_dim.N0 = 1 ;
    dst.m_dim.N1 = 1 ;
    dst.m_dim.N2 = 1 ;
    dst.m_dim.N3 = 1 ;
    dst.m_dim.N4 = 1 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N0 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<1u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = 0 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S1 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = 1 ;
    dst.m_dim.N2 = 1 ;
    dst.m_dim.N3 = 1 ;
    dst.m_dim.N4 = 1 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N1 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<2u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S2 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = 1 ;
    dst.m_dim.N3 = 1 ;
    dst.m_dim.N4 = 1 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N2 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<3u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S3 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = src.m_dim.N2 ;
    dst.m_dim.N3 = 1 ;
    dst.m_dim.N4 = 1 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N3 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<4u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S4 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = src.m_dim.N2 ;
    dst.m_dim.N3 = src.m_dim.N3 ;
    dst.m_dim.N4 = 1 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N4 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<5u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S5 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = src.m_dim.N2 ;
    dst.m_dim.N3 = src.m_dim.N3 ;
    dst.m_dim.N4 = src.m_dim.N4 ;
    dst.m_dim.N5 = 1 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N5 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<6u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = src.m_stride.S5 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S6 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = src.m_dim.N2 ;
    dst.m_dim.N3 = src.m_dim.N3 ;
    dst.m_dim.N4 = src.m_dim.N4 ;
    dst.m_dim.N5 = src.m_dim.N5 ;
    dst.m_dim.N6 = 1 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N6 );
  }
};
template <unsigned StaticDim> struct AssignFadDimStride<7u,StaticDim> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = src.m_stride.S5 ;
    dst.m_stride.S6 = src.m_stride.S6 ;
    dst.m_stride.S7 = src.m_stride.S7 ;

    dst.m_dim.N0 = src.m_dim.N0 ;
    dst.m_dim.N1 = src.m_dim.N1 ;
    dst.m_dim.N2 = src.m_dim.N2 ;
    dst.m_dim.N3 = src.m_dim.N3 ;
    dst.m_dim.N4 = src.m_dim.N4 ;
    dst.m_dim.N5 = src.m_dim.N5 ;
    dst.m_dim.N6 = src.m_dim.N6 ;
    AssignDim7<StaticDim>::eval( dst.m_dim, src.m_dim.N7 );
  }
};

// Specializations for subdynrankview
template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize ,
                    Kokkos::Impl::ViewSpecializeSacadoFad >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    ), Kokkos::Impl::DynRankSubviewTag >::type
  , SrcTraits
  , Args ... >
{
private:

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

  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  typedef Kokkos::LayoutStride array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

  typedef value_type******* data_type ;

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


  template< class MemoryTraits >
  struct apply {

    static_assert( Kokkos::is_memory_traits< MemoryTraits >::value , "" );

    typedef Kokkos::ViewTraits
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > traits_type ;

    typedef Kokkos::View
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > type ;
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 7 , rank > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
      return SubviewExtents< 7 , rank >( dim , arg0 , arg1 , arg2 , arg3 ,
                                         arg4 , arg5 , arg6 );
    }
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ArrayExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 8 , rank+1 > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
       return SubviewExtents< 8 , rank+1 >( dim , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 , Kokkos::ALL() );
    }
  };

  typedef DynRankView< value_type , array_layout , typename SrcTraits::device_type , typename SrcTraits::memory_traits >  ret_type;

  template < typename T , class ... P >
  KOKKOS_INLINE_FUNCTION
  static ret_type subview( const unsigned src_rank , Kokkos::DynRankView< T , P...> const & src , Args ... args )
  {

    typedef ViewMapping< traits_type, typename traits_type::specialize>  DstType ;
    typedef ViewMapping< SrcTraits, typename SrcTraits::specialize> SrcType;
    enum { FadStaticDim = SrcType::FadStaticDimension };
    typedef typename std::conditional< (rank==0) , ViewDimension<>
      , typename std::conditional< (rank==1) , ViewDimension<0>
      , typename std::conditional< (rank==2) , ViewDimension<0,0>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0>
      , ViewDimension<0,0,0,0,0,0,0>
      >::type >::type >::type >::type >::type >::type >::type  DstDimType ;
    typedef typename std::conditional< (rank==0) , ViewDimension<FadStaticDim>
      , typename std::conditional< (rank==1) , ViewDimension<0,FadStaticDim>
      , typename std::conditional< (rank==2) , ViewDimension<0,0,FadStaticDim>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0,FadStaticDim>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0,FadStaticDim>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0,FadStaticDim>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0,FadStaticDim>
      , ViewDimension<0,0,0,0,0,0,0,FadStaticDim>
      >::type >::type >::type >::type >::type >::type >::type  DstArrayDimType ;

      typedef ViewOffset< DstDimType , Kokkos::LayoutStride > dst_offset_type ;
      typedef ViewOffset< DstArrayDimType , Kokkos::LayoutStride > dst_array_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      ret_type dst ;

      const SubviewExtents< 7 , rank > extents =
        ExtentGenerator< Args ... >::generator(
          src.m_map.m_impl_offset.m_dim , args... ) ;
      const SubviewExtents< 8 , rank+1 > array_extents =
        ArrayExtentGenerator< Args ... >::generator(
          src.m_map.m_array_offset.m_dim , args... ) ;

      dst_offset_type tempdst( src.m_map.m_impl_offset , extents ) ;
      dst_array_offset_type temparraydst(
        src.m_map.m_array_offset , array_extents ) ;

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_offset.m_dim.N0 = tempdst.m_dim.N0 ;
      dst.m_map.m_impl_offset.m_dim.N1 = tempdst.m_dim.N1 ;
      dst.m_map.m_impl_offset.m_dim.N2 = tempdst.m_dim.N2 ;
      dst.m_map.m_impl_offset.m_dim.N3 = tempdst.m_dim.N3 ;
      dst.m_map.m_impl_offset.m_dim.N4 = tempdst.m_dim.N4 ;
      dst.m_map.m_impl_offset.m_dim.N5 = tempdst.m_dim.N5 ;
      dst.m_map.m_impl_offset.m_dim.N6 = tempdst.m_dim.N6 ;

      dst.m_map.m_impl_offset.m_stride.S0 = tempdst.m_stride.S0;
      dst.m_map.m_impl_offset.m_stride.S1 = tempdst.m_stride.S1;
      dst.m_map.m_impl_offset.m_stride.S2 = tempdst.m_stride.S2;
      dst.m_map.m_impl_offset.m_stride.S3 = tempdst.m_stride.S3;
      dst.m_map.m_impl_offset.m_stride.S4 = tempdst.m_stride.S4;
      dst.m_map.m_impl_offset.m_stride.S5 = tempdst.m_stride.S5;
      dst.m_map.m_impl_offset.m_stride.S6 = tempdst.m_stride.S6;

      // Move last non-unit dim and stride to N7/S7 since subview collapses
      // out all singleton dimensions between the last rank and the fad
      // dimension.  Equivalent to:
      //   dst.m_map.m_impl_offset.m_dim.N* = tempdst.m_dim.N*
      //   dst.m_map.m_impl_offset.m_dim.N7 = tempdst.m_dim.N{rank}
      //   dst.m_map.m_impl_offset.m_stride.S* = tempdst.m_stride.S*
      //   dst.m_map.m_impl_offset.m_stride.S7 = tempdst.m_stride.S{rank}
      AssignFadDimStride<rank,FadStaticDim>::eval( dst.m_map.m_array_offset, temparraydst );

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_handle =
        dst_handle_type(
          src.m_map.m_impl_handle +
          src.m_map.m_array_offset( array_extents.domain_offset(0)
                                  , array_extents.domain_offset(1)
                                  , array_extents.domain_offset(2)
                                  , array_extents.domain_offset(3)
                                  , array_extents.domain_offset(4)
                                  , array_extents.domain_offset(5)
                                  , array_extents.domain_offset(6)
                                  , array_extents.domain_offset(7)
            ) );

      dst.m_map.m_fad_size = src.m_map.m_fad_size;
      dst.m_map.m_fad_stride = src.m_map.m_fad_stride.value;

      dst.m_rank = ( src_rank > 0 ? unsigned(R0) : 0 )
                 + ( src_rank > 1 ? unsigned(R1) : 0 )
                 + ( src_rank > 2 ? unsigned(R2) : 0 )
                 + ( src_rank > 3 ? unsigned(R3) : 0 )
                 + ( src_rank > 4 ? unsigned(R4) : 0 )
                 + ( src_rank > 5 ? unsigned(R5) : 0 )
                 + ( src_rank > 6 ? unsigned(R6) : 0 ) ;

      return dst ;
    }
};

// ViewMapping for copy and copy-assign from View to DynRankView
template <unsigned> struct AssignFadDim7 {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src , const unsigned dim ) {}
};

template <> struct AssignFadDim7<0u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src , const unsigned dim ) {
    dst.m_dim.N7 = src.m_dim.extent(dim);
  }
};

// Copy a layout, moving the Fad dimension to the last
template <typename Layout>
KOKKOS_INLINE_FUNCTION
static Layout
permute_fad_layout(const Layout& src, const unsigned rank) {
  Layout dst = src;
  dst.dimension[rank] = 1;
  dst.dimension[7]    = src.dimension[rank];
  return dst;
}
KOKKOS_INLINE_FUNCTION
static LayoutStride
permute_fad_layout(const LayoutStride& src, const unsigned rank) {
  LayoutStride dst = src;
  dst.dimension[rank] = 1;
  dst.stride[rank]    = 1;
  dst.dimension[7]    = src.dimension[rank];
  dst.stride[7]       = src.stride[rank];
  return dst;
}

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  DynRankView<FAD>      = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has FAD only
    std::is_same< typename DstTraits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value
  ), Kokkos::Impl::ViewToDynRankViewTag >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template < typename DT , typename ... DP , typename ST , typename ... SP >
  KOKKOS_INLINE_FUNCTION static
  void assign( Kokkos::DynRankView< DT , DP... > & dst
             , const Kokkos::View< ST , SP... >& src )
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
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      typedef typename DstType::offset_type dst_offset_type;
      typedef typename DstType::array_offset_type dst_array_offset_type;
      dst.m_map.m_array_offset =
        dst_array_offset_type(std::integral_constant<unsigned,0>(),
                              permute_fad_layout(src.m_map.m_array_offset.layout(),
                                                 SrcTraits::rank) );
      dst.m_map.m_impl_offset =
        dst_offset_type(std::integral_constant<unsigned,0>(),
                        src.m_map.m_impl_offset.layout() );

      dst.m_map.m_impl_handle = src.m_map.m_impl_handle ;
      dst.m_rank = src.rank ;

      dst.m_map.m_fad_size = src.m_map.m_fad_size ;
      dst.m_map.m_fad_stride = src.m_map.m_fad_stride ;
    }
};

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<ordinary> = View<FAD>
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
  ), Kokkos::Impl::ViewToDynRankViewTag >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template < typename DT , typename ... DP , typename ST , typename ... SP >
  KOKKOS_INLINE_FUNCTION static
  void assign( Kokkos::DynRankView< DT , DP... > & dst
             , const Kokkos::View< ST , SP... >& src )
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
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      typedef typename DstType::offset_type dst_offset_type;
      dst.m_map.m_impl_offset =
        dst_offset_type(std::integral_constant<unsigned,0>(),
                        permute_fad_layout(src.m_map.m_array_offset.layout(),
                                           SrcTraits::rank));

      dst.m_map.m_impl_handle  = src.m_map.m_impl_handle ;
      dst.m_rank    = src.Rank ;
    }
};

}} //end Kokkos::Impl

namespace Kokkos {

template <typename T, typename ... P>
struct is_dynrankview_fad< DynRankView<T,P...> > {
  typedef DynRankView<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFadContiguous >::value;
};

template <typename T, typename ... P>
struct is_dynrankview_fad_contiguous< DynRankView<T,P...> > {
  typedef DynRankView<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Impl::ViewSpecializeSacadoFadContiguous >::value;
};

template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_dynrankview_fad< DynRankView<T,P...> >::value, unsigned >::type
dimension_scalar(const DynRankView<T,P...>& view) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  return view.implementation_map().dimension_scalar();
#else
  return view.impl_map().dimension_scalar();
#endif
}


// Overload of deep_copy for Fad views intializing to a constant scalar

template< class DT, class ... DP >
void deep_copy(
  const DynRankView<DT,DP...> & view ,
  const typename Sacado::ScalarType< typename DynRankView<DT,DP...>::value_type >::type & value
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

  Kokkos::fence();
  Kokkos::Impl::DynRankViewFill< DynRankView<DT,DP...> >( view , value );
  Kokkos::fence();
}

// Overload of deep_copy for Fad views intializing to a constant Fad
template< class DT, class ... DP >
void deep_copy(
  const DynRankView<DT,DP...> & view ,
  const typename DynRankView<DT,DP...>::value_type & value
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

  Kokkos::fence();
  Kokkos::Impl::DynRankViewFill< DynRankView<DT,DP...> >( view , value );
  Kokkos::fence();
}

template< class DstType , class SrcType >
inline
void deep_copy
  ( const DstType & dst
  , const SrcType & src
  , typename std::enable_if<(
  ( std::is_same< typename DstType::traits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename DstType::traits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
  &&
  ( std::is_same< typename SrcType::traits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename SrcType::traits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
  &&
  ( Kokkos::is_dyn_rank_view<DstType>::value || Kokkos::is_dyn_rank_view<SrcType
>::value )
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename DstType::traits::value_type ,
                  typename DstType::traits::non_const_value_type >::value
    , "deep_copy requires non-const destination type" );

  typedef DstType dst_type ;
  typedef SrcType src_type ;

  typedef typename dst_type::execution_space  dst_execution_space ;
  typedef typename src_type::execution_space  src_execution_space ;
  typedef typename dst_type::memory_space     dst_memory_space ;
  typedef typename src_type::memory_space     src_memory_space ;

  enum { DstExecCanAccessSrc =
   Kokkos::SpaceAccessibility< typename dst_execution_space::memory_space, src_memory_space >::accessible };

  enum { SrcExecCanAccessDst =
   Kokkos::SpaceAccessibility< typename src_execution_space::memory_space, dst_memory_space >::accessible };

  if ( (void *) dst.data() != (void*) src.data() ) {

    // Concern: If overlapping views then a parallel copy will be erroneous.
    // ...

    // If same type, equal layout, equal dimensions, equal span, and contiguous memory then can byte-wise copy
    if ( rank(src) == 0 && rank(dst) == 0 )
    {
      typedef typename dst_type::value_type::value_type value_type ;
      const size_t nbytes = sizeof(value_type) * dst.span() ;
      Kokkos::fence();
      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
      Kokkos::fence();
    }
    else if ( std::is_same< typename DstType::traits::value_type ,
                       typename SrcType::traits::non_const_value_type >::value &&
         (
           ( std::is_same< typename DstType::traits::array_layout ,
                           typename SrcType::traits::array_layout >::value
             &&
             ( std::is_same< typename DstType::traits::array_layout ,
                             typename Kokkos::LayoutLeft>::value
             ||
               std::is_same< typename DstType::traits::array_layout ,
                             typename Kokkos::LayoutRight>::value
             )
           )
           ||
           (
             rank(dst) == 1
             &&
             rank(src) == 1
           )
         ) &&
         dst.span_is_contiguous() &&
         src.span_is_contiguous() &&
         dst.span() == src.span() &&
         dst.extent(0) == src.extent(0) &&
         dst.extent(1) == src.extent(1) &&
         dst.extent(2) == src.extent(2) &&
         dst.extent(3) == src.extent(3) &&
         dst.extent(4) == src.extent(4) &&
         dst.extent(5) == src.extent(5) &&
         dst.extent(6) == src.extent(6) &&
         dst.extent(7) == src.extent(7) ) {

//      const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span() * dimension_scalar(dst) ;
      //const size_t nbytes = sizeof(typename dst_type::scalar_array_type) * dst.span() ;
      const size_t nbytes = sizeof(typename dst_type::value_type::value_type) * dst.span() ; 
      //dst_type::value_type is outer FAD type, dst_type::value_type::value_type is inner FAD type
      Kokkos::fence();
      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
      Kokkos::fence();
    }
    else if ( std::is_same< typename DstType::traits::value_type ,
                            typename SrcType::traits::non_const_value_type >::value &&
         (
           ( std::is_same< typename DstType::traits::array_layout ,
                           typename SrcType::traits::array_layout >::value
             &&
             std::is_same< typename DstType::traits::array_layout ,
                          typename Kokkos::LayoutStride>::value
           )
           ||
           (
             rank(dst) == 1
             &&
             rank(src) == 1
           )
         ) &&
         dst.span_is_contiguous() &&
         src.span_is_contiguous() &&
         dst.span() == src.span() &&
         dst.extent(0) == src.extent(0) &&
         dst.extent(1) == src.extent(1) &&
         dst.extent(2) == src.extent(2) &&
         dst.extent(3) == src.extent(3) &&
         dst.extent(4) == src.extent(4) &&
         dst.extent(5) == src.extent(5) &&
         dst.extent(6) == src.extent(6) &&
         dst.extent(7) == src.extent(7) &&
         dst.stride_0() == src.stride_0() &&
         dst.stride_1() == src.stride_1() &&
         dst.stride_2() == src.stride_2() &&
         dst.stride_3() == src.stride_3() &&
         dst.stride_4() == src.stride_4() &&
         dst.stride_5() == src.stride_5() &&
         dst.stride_6() == src.stride_6() &&
         dst.stride_7() == src.stride_7()
         ) {

      const size_t nbytes = sizeof(typename dst_type::value_type::value_type) * dst.span() ; 
      Kokkos::fence();
      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
      Kokkos::fence();
    }
    else if ( DstExecCanAccessSrc ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::fence();
      Kokkos::Impl::DynRankViewRemap< dst_type , src_type >( dst , src );
      Kokkos::fence();
    }
    else if ( SrcExecCanAccessDst ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::fence();
      Kokkos::Impl::DynRankViewRemap< dst_type , src_type , src_execution_space >( dst , src );
      Kokkos::fence();
    }
    else {
      Kokkos::Impl::throw_runtime_exception("deep_copy given views that would require a temporary allocation");
    }
  }
  else {
    Kokkos::fence();
  }
}

template< class T , class ... P >
inline
typename Kokkos::DynRankView<T,P...>::HostMirror
create_mirror( const Kokkos::DynRankView<T,P...> & src
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
  typedef DynRankView<T,P...>            src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  typename src_type::array_layout layout = src.layout();
  layout.dimension[src.rank()] = Kokkos::dimension_scalar(src);

  return dst_type(std::string(src.label()).append("_mirror"),
                  Impl::reconstructLayout(layout, src.rank()+1));
}

template< class T , class ... P >
inline
typename Kokkos::DynRankView<T,P...>::HostMirror
create_mirror( const Kokkos::DynRankView<T,P...> & src
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
  typedef DynRankView<T,P...>            src_type ;
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

  layout.dimension[src.rank()] = Kokkos::dimension_scalar(src);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  layout.stride[src.rank()] = src.implementation_map().stride_scalar();
#else
  layout.stride[src.rank()] = src.impl_map().stride_scalar();
#endif

  return dst_type(std::string(src.label()).append("_mirror"),
                  Impl::reconstructLayout(layout, src.rank()+1));
}

template<class Space, class T, class ... P>
typename Impl::MirrorDRViewType<Space,T,P ...>::dest_view_type
create_mirror(const Space& , const Kokkos::DynRankView<T,P...> & src
             , typename std::enable_if<
                 ( std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
                   std::is_same< typename ViewTraits<T,P...>::specialize ,
                     Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
               >::type *) {
  typedef DynRankView<T,P...> src_type ;
  typename src_type::array_layout layout = src.layout();
  layout.dimension[src.rank()] = Kokkos::dimension_scalar(src);
  return typename Impl::MirrorDRViewType<Space,T,P ...>::dest_view_type(
    src.label(),Impl::reconstructLayout(layout, src.rank()+1));
}

namespace Impl {

template <unsigned N, typename T, typename... Args>
KOKKOS_FUNCTION auto as_view_of_rank_n(
  DynRankView<T, Args...> v,
  typename std::enable_if<
    ( std::is_same< typename ViewTraits<T,Args...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
      std::is_same< typename ViewTraits<T,Args...>::specialize ,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
    >::type*)
{
  if (v.rank() != N) {
    KOKKOS_IF_ON_HOST(
        const std::string message =
            "Converting DynRankView of rank " + std::to_string(v.rank()) +
            " to a View of mis-matched rank " + std::to_string(N) + "!";
        Kokkos::abort(message.c_str());)
    KOKKOS_IF_ON_DEVICE(
        Kokkos::abort("Converting DynRankView to a View of mis-matched rank!");)
  }

  auto layout = v.impl_map().layout();
  layout.dimension[v.rank()] = Kokkos::dimension_scalar(v);
  return View<typename RankDataType<T, N>::type, Args...>(
      v.data(), Impl::reconstructLayout(layout, v.rank()+1));
}

}

} // end Kokkos

#endif //defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#include "Kokkos_DynRankView_Fad_Contiguous.hpp"

#endif /* #ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP */
