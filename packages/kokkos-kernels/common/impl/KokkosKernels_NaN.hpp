//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSKERNELS_NAN_HPP
#define KOKKOSKERNELS_NAN_HPP

#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_NumericTraits.hpp>

namespace KokkosKernels::Impl {

// This could be constexpr if Kokkos::complex ctor was
template <typename T>
KOKKOS_INLINE_FUNCTION T quiet_NaN() {
  if constexpr (std::is_same_v<double, T>) {
    return double(Kokkos::Experimental::quiet_NaN_v<float>);  // Kokkos::Experimetnal::quiet_NaN_v<double>
                                                              // is undefined in
                                                              // device code
  } else if constexpr (Kokkos::ArithTraits<T>::is_complex) {
    using value_type = typename T::value_type;
    return T(quiet_NaN<value_type>(),
             quiet_NaN<value_type>());  // Kokkos::complex ctor is not constexpr
  } else {
    return Kokkos::Experimental::quiet_NaN_v<T>;
  }
}

}  // namespace KokkosKernels::Impl

#endif  // KOKKOSKERNELS_NAN_HPP
