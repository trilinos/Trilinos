// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_ROTG_IMPL_HPP_
#define KOKKOSBLAS1_ROTG_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rotg_spec.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Compute Givens rotation coefficients for real scalars.
/// \tparam Scalar The real scalar type for a and b
/// \tparam Magnitude The real scalar type for c
/// \param[in,out] a The first scalar, overwritten by the first rotation coefficient
/// \param[in,out] b The second scalar, overwritten by the second rotation coefficient
/// \param[out] c The cosine of the rotation
/// \param[out] s The sine of the rotation
template <class Scalar, class Magnitude>
  requires(!KokkosKernels::ArithTraits<Scalar>::is_complex)
KOKKOS_INLINE_FUNCTION void rotg_impl(Scalar* KOKKOS_RESTRICT a, Scalar* KOKKOS_RESTRICT b,
                                      Magnitude* KOKKOS_RESTRICT c, Scalar* KOKKOS_RESTRICT s) {
  const Scalar one   = KokkosKernels::ArithTraits<Scalar>::one();
  const Scalar zero  = KokkosKernels::ArithTraits<Scalar>::zero();
  const Scalar anorm = Kokkos::abs(*a);
  const Scalar bnorm = Kokkos::abs(*b);

  if (bnorm == zero) {
    *c = one;
    *s = zero;
    *b = zero;
    return;
  }

  if (anorm == zero) {
    *c = zero;
    *s = one;
    *a = bnorm;
    *b = one;
    return;
  }

  const Scalar scale = anorm + bnorm;
  const Scalar sign  = anorm > bnorm ? *a : *b;
  Scalar norm        = scale * Kokkos::hypot(*a / scale, *b / scale);
  norm               = Kokkos::copysign(norm, sign);

  *c = *a / norm;
  *s = *b / norm;
  if (anorm > bnorm) {
    *b = *s;
  } else if (*c != zero) {
    *b = one / *c;
  } else {
    *b = one;
  }
  *a = norm;
}

/// \brief Compute Givens rotation coefficients for complex scalars.
/// \tparam Scalar The complex scalar type for a and b
/// \tparam Magnitude The real scalar type for c
/// \param[in,out] a The first scalar, overwritten by the first rotation coefficient
/// \param[in] b The second scalar
/// \param[out] c The cosine of the rotation
/// \param[out] s The sine of the rotation
template <class Scalar, class Magnitude>
  requires KokkosKernels::ArithTraits<Scalar>::is_complex
KOKKOS_INLINE_FUNCTION void rotg_impl(Scalar* KOKKOS_RESTRICT a, const Scalar* KOKKOS_RESTRICT b,
                                      Magnitude* KOKKOS_RESTRICT c, Scalar* KOKKOS_RESTRICT s) {
  using mag_type    = typename KokkosKernels::ArithTraits<Scalar>::mag_type;
  const Scalar zero = KokkosKernels::ArithTraits<Scalar>::zero();

  if (Kokkos::abs(*b) == zero) {
    *c = KokkosKernels::ArithTraits<mag_type>::one();
    *s = zero;
    return;
  }

  if (Kokkos::abs(*a) == zero) {
    *c = KokkosKernels::ArithTraits<mag_type>::zero();
    *s = Kokkos::conj(*b) / Kokkos::abs(*b);
    *a = Scalar(Kokkos::abs(*b));
    return;
  }

  const mag_type scale = Kokkos::abs(*a) + Kokkos::abs(*b);
  const mag_type norm  = scale * Kokkos::hypot(Kokkos::abs(*a) / scale, Kokkos::abs(*b) / scale);
  const Scalar unit_a  = *a / Kokkos::abs(*a);
  *c                   = Kokkos::abs(*a) / norm;
  *s                   = unit_a * Kokkos::conj(*b) / norm;
  *a                   = unit_a * norm;
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
