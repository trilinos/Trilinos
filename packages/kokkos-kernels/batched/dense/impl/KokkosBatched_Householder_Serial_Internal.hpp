// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_HOUSEHOLDER_SERIAL_INTERNAL_HPP
#define KOKKOSBATCHED_HOUSEHOLDER_SERIAL_INTERNAL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialLeftHouseholderInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m_x2,
                                           /* */ ValueType* chi1,
                                           /* */ ValueType* x2, const int x2s,
                                           /* */ ValueType* tau) {
    using value_type = ValueType;
    using KAT        = KokkosKernels::ArithTraits<value_type>;
    using mag_type   = typename KAT::mag_type;
    using KAT_mag    = KokkosKernels::ArithTraits<mag_type>;

    const mag_type zero      = KAT_mag::zero();
    const mag_type one       = KAT_mag::one();
    const mag_type half      = one / (one + one);
    const mag_type minus_one = -one;

    /// compute the 2norm of x2
    mag_type norm_x2_square = zero;
    for (int i = 0; i < m_x2; ++i) {
      const auto x2_at_i = x2[i * x2s];
      norm_x2_square += Kokkos::real(Kokkos::conj(x2_at_i) * x2_at_i);
    }

    /// if norm_x2 is zero, return with trivial values
    if (norm_x2_square == zero) {
      *chi1 = -(*chi1);
      *tau  = half * KAT::one();

      return 0;
    }

    /// compute magnitude of chi1, equal to norm2 of chi1
    const mag_type norm_chi1 = KAT::abs(*chi1);

    /// compute 2 norm of x using norm_chi1 and norm_x2
    const mag_type norm_x = KAT_mag::sqrt(norm_x2_square + norm_chi1 * norm_chi1);

    /// compute alpha
    const mag_type alpha = (Kokkos::real(*chi1) < zero ? one : minus_one) * norm_x;

    /// overwrite x2 with u2
    const value_type chi1_minus_alpha     = *chi1 - alpha;
    const value_type inv_chi1_minus_alpha = one / chi1_minus_alpha;
    for (int i = 0; i < m_x2; ++i) x2[i * x2s] *= inv_chi1_minus_alpha;

    // later consider to use the following
    // SerialScaleInternal::invoke(m_x2, inv_chi1_minus_alpha, x2, x2s);

    /// compute tau
    // Note that in the complex case we have
    // multiple possible expressions for tau
    // we chose the same as LAPACK which
    // guarentees that R is real valued.
    const mag_type chi1_minus_alpha_square = Kokkos::abs(chi1_minus_alpha) * Kokkos::abs(chi1_minus_alpha);
    if constexpr (KAT::is_complex) {
      *tau = alpha / (alpha - *chi1);
    } else {
      *tau = half + half * (norm_x2_square / chi1_minus_alpha_square);
    }

    /// overwrite chi1 with alpha
    *chi1 = alpha;

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
