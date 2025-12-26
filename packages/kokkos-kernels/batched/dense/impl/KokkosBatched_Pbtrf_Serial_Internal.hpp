// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"
#include "KokkosBatched_Syr_Serial_Internal.hpp"
#include "KokkosBatched_Lacgv_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

///
/// Lower
///

template <typename AlgoType>
struct SerialPbtrfInternalLower {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an,
                                           /**/ ValueType *KOKKOS_RESTRICT AB, const int as0, const int as1,
                                           const int kd);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalLower<Algo::Pbtrf::Unblocked>::invoke(const int an,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT AB,
                                                                                    const int as0, const int as1,
                                                                                    const int kd) {
  // Compute the Cholesky factorization A = L*L'.
  for (int j = 0; j < an; ++j) {
    auto a_jj = KokkosKernels::ArithTraits<ValueType>::real(AB[0 * as0 + j * as1]);

    // Check if L (j, j) is positive definite
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (a_jj <= 0) {
      return j + 1;
    }
#endif

    a_jj                  = Kokkos::sqrt(a_jj);
    AB[0 * as0 + j * as1] = a_jj;

    // Compute elements J+1:J+KN of column J and update the
    // trailing submatrix within the band.
    int kn = Kokkos::min(an - j - 1, kd);
    if (kn > 0) {
      // scale to diagonal elements
      const ValueType alpha = 1.0 / a_jj;
      KokkosBlas::Impl::SerialScaleInternal::invoke(kn, alpha, &(AB[1 * as0 + j * as1]), 1);

      // syr or zher (lower) with alpha = -1.0 to diagonal elements
      using op     = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpConj,
                                    KokkosBlas::Impl::OpID>;
      using op_sym = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpReal,
                                        KokkosBlas::Impl::OpID>;
      SerialSyrInternalLower::invoke(op(), op_sym(), kn, -1.0, &(AB[1 * as0 + j * as1]), as0,
                                     &(AB[0 * as0 + (j + 1) * as1]), as0, (as1 - as0));
    }
  }

  return 0;
}

///
/// Upper
///

template <typename AlgoType>
struct SerialPbtrfInternalUpper {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an,
                                           /**/ ValueType *KOKKOS_RESTRICT AB, const int as0, const int as1,
                                           const int kd);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalUpper<Algo::Pbtrf::Unblocked>::invoke(const int an,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT AB,
                                                                                    const int as0, const int as1,
                                                                                    const int kd) {
  // Compute the Cholesky factorization A = U'*U.
  for (int j = 0; j < an; ++j) {
    auto a_jj = KokkosKernels::ArithTraits<ValueType>::real(AB[kd * as0 + j * as1]);

    // Check if U (j,j) is positive definite
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (a_jj <= 0) {
      return j + 1;
    }
#endif
    a_jj                   = Kokkos::sqrt(a_jj);
    AB[kd * as0 + j * as1] = a_jj;

    // Compute elements J+1:J+KN of row J and update the
    // trailing submatrix within the band.
    int kn  = Kokkos::min(kd, an - j - 1);
    int kld = Kokkos::max(1, as0 - 1);
    if (kn > 0) {
      // scale to diagonal elements
      const ValueType alpha = 1.0 / a_jj;
      KokkosBlas::Impl::SerialScaleInternal::invoke(kn, alpha, &(AB[(kd - 1) * as0 + (j + 1) * as1]), kld);

      // zlacgv to diagonal elements (no op for real matrix)
      SerialLacgvInternal::invoke(kn, &(AB[(kd - 1) * as0 + (j + 1) * as1]), (as0 - as1));

      // syr or zher (upper) with alpha = -1.0 to diagonal elements
      using op     = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpConj,
                                    KokkosBlas::Impl::OpID>;
      using op_sym = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpReal,
                                        KokkosBlas::Impl::OpID>;
      SerialSyrInternalUpper::invoke(op(), op_sym(), kn, -1.0, &(AB[(kd - 1) * as0 + (j + 1) * as1]), as0,
                                     &(AB[kd * as0 + (j + 1) * as1]), as0, (as1 - as0));

      // zlacgv to diagonal elements (no op for real matrix)
      SerialLacgvInternal::invoke(kn, &(AB[(kd - 1) * as0 + (j + 1) * as1]), (as0 - as1));
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_
