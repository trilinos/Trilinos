// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_KOKKOS_INNERPRODUCTSPACETRAITS_HPP
#define KOKKOSKERNELS_KOKKOS_INNERPRODUCTSPACETRAITS_HPP

/// \file Kokkos_InnerProductSpaceTraits.hpp
/// \brief Declaration and definition of
/// Kokkos::Details::InnerProductSpaceTraits

#include <KokkosKernels_InnerProductSpaceTraits.hpp>

#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_5
static_assert(false,
              "Header `Kokkos_InnerProductSpaceTraits.hpp` is deprecated, include "
              "`KokkosKernels_InnerProductSpaceTraits.hpp` instead");
#endif

namespace Kokkos {
namespace Details {

template <typename T>
using InnerProductSpaceTraits
    [[deprecated("Use KokkosKernels::Details::InnerProductSpaceTraits from "
                 "KokkosKernels_InnerProductSpaceTraits.hpp header instead")]] =
        ::KokkosKernels::Details::InnerProductSpaceTraits<T>;

template <typename Out, typename In>
using CastPossiblyComplex
    [[deprecated("Use KokkosKernels::Details::CastPossiblyComplex from "
                 "KokkosKernels_InnerProductSpaceTraits.hpp header instead")]] =
        ::KokkosKernels::Details::CastPossiblyComplex<Out, In>;

using ::KokkosKernels::Details::updateDot;

}  // namespace Details
}  // namespace Kokkos

#endif  // KOKKOSKERNELS_KOKKOS_INNERPRODUCTSPACETRAITS_HPP
