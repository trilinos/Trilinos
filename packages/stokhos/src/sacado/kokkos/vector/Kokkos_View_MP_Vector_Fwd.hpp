// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_MP_VECTOR_FWD_HPP
#define KOKKOS_VIEW_MP_VECTOR_FWD_HPP

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

//----------------------------------------------------------------------------

namespace Sacado {
  namespace MP {
    template <typename Storage >
    class Vector;
  }
}

namespace Kokkos {

  namespace Impl {
    template<class Space, class T, class ... P>
    struct MirrorViewType;
  }

}

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ViewMPVectorContiguous;

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

namespace Kokkos {

// Declare overloads of create_mirror() so they are in scope
// Kokkos_Core.hpp is included below

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(
  const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class T, class... P>
inline auto create_mirror(
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class Space, class T, class... P, typename Enable = std::enable_if_t<Kokkos::is_space<Space>::value>>
inline auto create_mirror(
  const Space& space,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class T, class... P>
inline auto create_mirror(
  Impl::WithoutInitializing_t wi,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class Space, class T, class... P, typename Enable = std::enable_if_t<is_space<Space>::value>>
inline auto create_mirror(
  Impl::WithoutInitializing_t wi,
  const Space& space,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror_view(
  const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
  const Kokkos::View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>* = nullptr
);

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize,
            Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value &&
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr);

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize,
            Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value &&
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr);

// Overload of deep_copy for MP::Vector views intializing to a constant scalar
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for MP::Vector views intializing to a constant MP::Vector
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for MP::Vector views intializing to a constant scalar
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for MP::Vector views intializing to a constant MP::Vector
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
                 )>::type * = 0 );

/* Specialize for deep copy of MP::Vector */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    )>::type * = 0 );

/* Specialize for deep copy of MP::Vector */
template< class ExecSpace, class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const ExecSpace &,
                const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    )>::type * = 0 );

namespace Impl {

template <unsigned N, typename... Args>
KOKKOS_FUNCTION std::enable_if_t<
    N == View<Args...>::rank &&
    std::is_same<typename ViewTraits<Args...>::specialize,
                 Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value,
    View<Args...>>
as_view_of_rank_n(View<Args...> v);

// Placeholder implementation to compile generic code for DynRankView; should
// never be called
template <unsigned N, typename T, typename... Args>
std::enable_if_t<
    N != View<T, Args...>::rank &&
        std::is_same<typename ViewTraits<T, Args...>::specialize,
                     Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value,
    View<typename RankDataType<typename View<T, Args...>::value_type, N>::type,
         Args...>>
as_view_of_rank_n(View<T, Args...>);

}

} // namespace Kokkos

#endif /* #ifndef KOKKOS_VIEW_MP_VECTOR_FWD_HPP */
