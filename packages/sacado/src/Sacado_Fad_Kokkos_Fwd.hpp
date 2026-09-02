// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_FWD_HPP
#define SACADO_FAD_KOKKOS_FWD_HPP

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

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

namespace Sacado {

// Whether a given type is a view with Sacado FAD scalar type
template <typename view_type>
struct is_view_fad;

}

#ifdef SACADO_ENABLE_DEPRECATED_CODE
namespace Kokkos {
template<typename T>
using is_view_fad SACADO_DEPRECATED_WITH_COMMENT(
    "Use Sacado::is_view_fad instead of Kokkos::is_view_fad") = Sacado::is_view_fad<T>;
}
#endif

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef SACADO_FAD_KOKKOS_FWD_HPP */
