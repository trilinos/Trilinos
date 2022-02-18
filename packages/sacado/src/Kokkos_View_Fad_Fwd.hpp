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

#ifndef KOKKOS_VIEW_FAD_FWD_HPP
#define KOKKOS_VIEW_FAD_FWD_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOSCORE)

// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_View.hpp"

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

template<class Space, class T, class ... P>
struct MirrorType;

} // namespace Impl

// Declare overloads of create_mirror() so they are in scope
// Kokkos_Core.hpp is included later

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

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOSCORE)

#endif /* #ifndef KOKKOS_VIEW_FAD_FWD_HPP */
