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
#ifndef KOKKOSBLAS1_ROTG_IMPL_HPP_
#define KOKKOSBLAS1_ROTG_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rotg_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class Scalar, class Magnitude,
          typename std::enable_if<!Kokkos::ArithTraits<Scalar>::is_complex, bool>::type = true>
KOKKOS_INLINE_FUNCTION void rotg_impl(Scalar* a, Scalar* b, Magnitude* c, Scalar* s) {
  const Scalar one  = Kokkos::ArithTraits<Scalar>::one();
  const Scalar zero = Kokkos::ArithTraits<Scalar>::zero();

  const Scalar numerical_scaling = Kokkos::abs(*a) + Kokkos::abs(*b);
  if (numerical_scaling == zero) {
    *c = one;
    *s = zero;
    *a = zero;
    *b = zero;
  } else {
    const Scalar scaled_a = *a / numerical_scaling;
    const Scalar scaled_b = *b / numerical_scaling;
    Scalar norm           = Kokkos::sqrt(scaled_a * scaled_a + scaled_b * scaled_b) * numerical_scaling;
    Scalar sign           = Kokkos::abs(*a) > Kokkos::abs(*b) ? *a : *b;
    norm                  = Kokkos::copysign(norm, sign);
    *c                    = *a / norm;
    *s                    = *b / norm;

    Scalar z = one;
    if (Kokkos::abs(*a) > Kokkos::abs(*b)) {
      z = *s;
    }
    if ((Kokkos::abs(*b) >= Kokkos::abs(*a)) && (*c != zero)) {
      z = one / *c;
    }
    *a = norm;
    *b = z;
  }
}

template <class Scalar, class Magnitude,
          typename std::enable_if<Kokkos::ArithTraits<Scalar>::is_complex, bool>::type = true>
KOKKOS_INLINE_FUNCTION void rotg_impl(Scalar* a, Scalar* b, Magnitude* c, Scalar* s) {
  using mag_type = typename Kokkos::ArithTraits<Scalar>::mag_type;

  const Scalar one        = Kokkos::ArithTraits<Scalar>::one();
  const Scalar zero       = Kokkos::ArithTraits<Scalar>::zero();
  const mag_type mag_zero = Kokkos::ArithTraits<mag_type>::zero();

  const mag_type numerical_scaling = Kokkos::abs(*a) + Kokkos::abs(*b);
  if (Kokkos::abs(*a) == zero) {
    *c = mag_zero;
    *s = one;
    *a = *b;
  } else {
    const Scalar scaled_a = Kokkos::abs(*a / numerical_scaling);
    const Scalar scaled_b = Kokkos::abs(*b / numerical_scaling);
    mag_type norm         = Kokkos::abs(Kokkos::sqrt(scaled_a * scaled_a + scaled_b * scaled_b)) * numerical_scaling;
    Scalar unit_a         = *a / Kokkos::abs(*a);
    *c                    = Kokkos::abs(*a) / norm;
    *s                    = unit_a * Kokkos::conj(*b) / norm;
    *a                    = unit_a * norm;
  }
}

template <class SViewType, class MViewType>
struct rotg_functor {
  SViewType a, b;
  MViewType c;
  SViewType s;

  rotg_functor(SViewType const& a_, SViewType const& b_, MViewType const& c_, SViewType const& s_)
      : a(a_), b(b_), c(c_), s(s_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int const) const { rotg_impl(a.data(), b.data(), c.data(), s.data()); }
};

/// \brief Compute Givens rotation coefficients.
template <class ExecutionSpace, class SViewType, class MViewType>
void Rotg_Invoke(ExecutionSpace const& space, SViewType const& a, SViewType const& b, MViewType const& c,
                 SViewType const& s) {
  Kokkos::RangePolicy<ExecutionSpace> rotg_policy(space, 0, 1);
  rotg_functor rotg_func(a, b, c, s);
  Kokkos::parallel_for("KokkosBlas::rotg", rotg_policy, rotg_func);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTG_IMPL_HPP_
