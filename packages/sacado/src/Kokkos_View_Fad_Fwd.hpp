// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_FAD_FWD_HPP
#define KOKKOS_VIEW_FAD_FWD_HPP

#include "Sacado_ConfigDefs.h"
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
#include "Kokkos_View.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

namespace Kokkos {

// Whether a given type is a view with Sacado FAD scalar type
template <typename view_type>
struct is_view_fad;

}

// Make sure the user really wants these View specializations
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

namespace Kokkos {
namespace Impl {

struct ViewSpecializeSacadoFad;
struct ViewSpecializeSacadoFadContiguous;

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

template<class ExecutionSpace,
         class DT, class ... DP,
         class ST, class ... SP>
typename std::enable_if< is_view_fad< Kokkos::View<DT,DP...> >::value &&
                         is_view_fad< Kokkos::View<ST,SP...> >::value
                       >::type
view_copy(const ExecutionSpace& space,
          const Kokkos::View<DT,DP...>& dst, const Kokkos::View<ST,SP...>& src);

template<class Space, class T, class ... P>
struct MirrorViewType;

} // namespace Impl

// Declare overloads of create_mirror() so they are in scope
// Kokkos_Core.hpp is included later

template< class T , class ... P >
inline
typename std::enable_if<
  ( std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    !std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value,
  typename Kokkos::View<T,P...>::HostMirror>::type
create_mirror(const Kokkos::View<T,P...> & src);

template< class T , class ... P >
inline
typename std::enable_if<
  ( std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value,
  typename Kokkos::View<T,P...>::HostMirror>::type
create_mirror(const Kokkos::View<T,P...> & src);

template<class Space, class T, class ... P,
         typename Enable = std::enable_if_t<Kokkos::is_space<Space>::value>>
typename std::enable_if<
  std::is_same< typename ViewTraits<T,P...>::specialize ,
    Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
  std::is_same< typename ViewTraits<T,P...>::specialize ,
    Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value,
  typename Impl::MirrorViewType<Space,T,P ...>::dest_view_type>::type
create_mirror(const Space&, const Kokkos::View<T,P...> & src);

template< class T , class ... P >
inline
typename std::enable_if<
  ( std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    !std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value,
  typename Kokkos::View<T,P...>::HostMirror>::type
create_mirror(Kokkos::Impl::WithoutInitializing_t wi,
              const Kokkos::View<T,P...> & src);

template< class T , class ... P >
inline
typename std::enable_if<
  ( std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
     std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
    std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value,
  typename Kokkos::View<T,P...>::HostMirror>::type
create_mirror(Kokkos::Impl::WithoutInitializing_t wi,
              const Kokkos::View<T,P...> & src);

template<class Space, class T, class ... P,
          typename Enable = std::enable_if_t<Kokkos::is_space<Space>::value>>
typename std::enable_if<
  ( std::is_same< typename ViewTraits<T,P...>::specialize ,
       Kokkos::Impl::ViewSpecializeSacadoFad >::value ||
     std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ),
  typename Impl::MirrorViewType<Space,T,P ...>::dest_view_type>::type
create_mirror(Kokkos::Impl::WithoutInitializing_t wi,
              const Space&, const Kokkos::View<T,P...> & src);

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        ( std::is_same<typename ViewTraits<T, P...>::specialize,
              Kokkos::Impl::ViewSpecializeSacadoFad>::value ||
          std::is_same< typename ViewTraits<T,P...>::specialize ,
              Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr);

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        ( std::is_same<typename ViewTraits<T, P...>::specialize,
              Kokkos::Impl::ViewSpecializeSacadoFad>::value ||
          std::is_same< typename ViewTraits<T,P...>::specialize ,
              Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value ) &&
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr);

namespace Impl {

template <unsigned N, typename... Args>
KOKKOS_FUNCTION std::enable_if_t<
    N == View<Args...>::Rank &&
    (std::is_same<typename ViewTraits<Args...>::specialize,
                  Kokkos::Impl::ViewSpecializeSacadoFad>::value ||
     std::is_same<typename ViewTraits<Args...>::specialize,
                  Kokkos::Impl::ViewSpecializeSacadoFadContiguous>::value),
    View<Args...>>
as_view_of_rank_n(View<Args...> v);

// Placeholder implementation to compile generic code for DynRankView; should
// never be called
template <unsigned N, typename T, typename... Args>
std::enable_if_t<
    N != View<T, Args...>::Rank &&
        (std::is_same<typename ViewTraits<T, Args...>::specialize,
                      Kokkos::Impl::ViewSpecializeSacadoFad>::value ||
         std::is_same<typename ViewTraits<T, Args...>::specialize,
                      Kokkos::Impl::ViewSpecializeSacadoFadContiguous>::value),
    View<typename RankDataType<typename View<T, Args...>::value_type, N>::type,
         Args...>>
as_view_of_rank_n(View<T, Args...>);

}

namespace Experimental {

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 1))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 2))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 3))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 4))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 5))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 6))>::type* = nullptr);

template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 7))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 1))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 2))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 3))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 4))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 5))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 6))>::type* = nullptr);

template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(
      ( std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFad >::value
        ||
        std::is_same< typename ViewTraits<DT,DP...>::specialize,
        Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value )
      && (unsigned(ViewTraits<DT, DP...>::rank) == 7))>::type* = nullptr);

}

} // namespace Kokkos

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef KOKKOS_VIEW_FAD_FWD_HPP */
