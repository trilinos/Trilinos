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

#ifndef KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"

namespace KokkosBatched {

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

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an,
                                           /**/ Kokkos::complex<ValueType> *KOKKOS_RESTRICT AB, const int as0,
                                           const int as1, const int kd);
};

///
/// Real matrix
///

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalLower<Algo::Pbtrf::Unblocked>::invoke(const int an,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT AB,
                                                                                    const int as0, const int as1,
                                                                                    const int kd) {
  // Compute the Cholesky factorization A = L*L'.
  for (int j = 0; j < an; ++j) {
    auto a_jj = AB[0 * as0 + j * as1];

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

      // syr (lower) with alpha = -1.0 to diagonal elements
      for (int k = 0; k < kn; ++k) {
        auto x_k = AB[(k + 1) * as0 + j * as1];
        if (x_k != 0) {
          auto temp = -1.0 * x_k;
          for (int i = k; i < kn; ++i) {
            auto x_i = AB[(i + 1) * as0 + j * as1];
            AB[i * as0 + (j + 1 + k - i) * as1] += x_i * temp;
          }
        }
      }
    }
  }

  return 0;
}

///
/// Complex matrix
///
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalLower<Algo::Pbtrf::Unblocked>::invoke(
    const int an,
    /**/ Kokkos::complex<ValueType> *KOKKOS_RESTRICT AB, const int as0, const int as1, const int kd) {
  // Compute the Cholesky factorization A = L*L**H
  for (int j = 0; j < an; ++j) {
    auto a_jj = AB[0 * as0 + j * as1].real();

    // Check if L (j, j) is positive definite
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (a_jj <= 0) {
      AB[0 * as0 + j * as1] = a_jj;
      return j + 1;
    }
#endif

    a_jj                  = Kokkos::sqrt(a_jj);
    AB[0 * as0 + j * as1] = a_jj;

    // Compute elements J+1:J+KN of column J and update the
    // trailing submatrix within the band.
    int kn = Kokkos::min(kd, an - j - 1);
    if (kn > 0) {
      // scale to diagonal elements
      const ValueType alpha = 1.0 / a_jj;
      KokkosBlas::Impl::SerialScaleInternal::invoke(kn, alpha, &(AB[1 * as0 + j * as1]), 1);

      // zher (lower) with alpha = -1.0 to diagonal elements
      for (int k = 0; k < kn; ++k) {
        auto x_k = AB[(k + 1) * as0 + j * as1];
        if (x_k != 0) {
          auto temp                   = -1.0 * Kokkos::conj(x_k);
          AB[k * as0 + (j + 1) * as1] = AB[k * as0 + (j + 1) * as1].real() + (temp * x_k).real();
          for (int i = k + 1; i < kn; ++i) {
            auto x_i = AB[(i + 1) * as0 + j * as1];
            AB[i * as0 + (j + 1 + k - i) * as1] += x_i * temp;
          }
        } else {
          AB[k * as0 + (j + 1) * as1] = AB[k * as0 + (j + 1) * as1].real();
        }
      }
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

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an,
                                           /**/ Kokkos::complex<ValueType> *KOKKOS_RESTRICT AB, const int as0,
                                           const int as1, const int kd);
};

///
/// Real matrix
///
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalUpper<Algo::Pbtrf::Unblocked>::invoke(const int an,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT AB,
                                                                                    const int as0, const int as1,
                                                                                    const int kd) {
  // Compute the Cholesky factorization A = U'*U.
  for (int j = 0; j < an; ++j) {
    auto a_jj = AB[kd * as0 + j * as1];

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

      // syr (upper) with alpha = -1.0 to diagonal elements
      for (int k = 0; k < kn; ++k) {
        auto x_k = AB[(k + kd - 1) * as0 + (j + 1 - k) * as1];
        if (x_k != 0) {
          auto temp = -1.0 * x_k;
          for (int i = 0; i < k + 1; ++i) {
            auto x_i = AB[(i + kd - 1) * as0 + (j + 1 - i) * as1];
            AB[(kd + i) * as0 + (j + 1 + k - i) * as1] += x_i * temp;
          }
        }
      }
    }
  }

  return 0;
}

///
/// Complex matrix
///
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrfInternalUpper<Algo::Pbtrf::Unblocked>::invoke(
    const int an,
    /**/ Kokkos::complex<ValueType> *KOKKOS_RESTRICT AB, const int as0, const int as1, const int kd) {
  // Compute the Cholesky factorization A = U**H * U.
  for (int j = 0; j < an; ++j) {
    auto a_jj = AB[kd * as0 + j * as1].real();

    // Check if U (j,j) is positive definite
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    if (a_jj <= 0) {
      AB[kd * as0 + j * as1] = a_jj;
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

      // zlacgv to diagonal elements
      for (int i = 0; i < kn; ++i) {
        AB[(i + kd - 1) * as0 + (j + 1 - i) * as1] = Kokkos::conj(AB[(i + kd - 1) * as0 + (j + 1 - i) * as1]);
      }

      // zher (upper) with alpha = -1.0 to diagonal elements
      for (int k = 0; k < kn; ++k) {
        auto x_k = AB[(k + kd - 1) * as0 + (j + 1 - k) * as1];
        if (x_k != 0) {
          auto temp = -1.0 * Kokkos::conj(x_k);
          for (int i = 0; i < k + 1; ++i) {
            auto x_i = AB[(i + kd - 1) * as0 + (j + 1 - i) * as1];
            AB[(kd + i) * as0 + (j + 1 + k - i) * as1] += x_i * temp;
          }
        }
      }

      // zlacgv to diagonal elements
      for (int i = 0; i < kn; ++i) {
        AB[(i + kd - 1) * as0 + (j + 1 - i) * as1] = Kokkos::conj(AB[(i + kd - 1) * as0 + (j + 1 - i) * as1]);
      }
    }
  }

  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRF_SERIAL_INTERNAL_HPP_
