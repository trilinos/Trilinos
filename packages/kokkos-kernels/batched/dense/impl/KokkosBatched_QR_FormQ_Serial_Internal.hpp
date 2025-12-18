// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_QR_FORM_Q_SERIAL_INTERNAL_HPP
#define KOKKOSBATCHED_QR_FORM_Q_SERIAL_INTERNAL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"
#include "KokkosBatched_ApplyQ_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialQR_FormQ_Internal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int k,
                                           /* */ ValueType* A, const int as0, const int as1,
                                           /* */ ValueType* t, const int ts,
                                           /* */ ValueType* Q, const int qs0, const int qs1,
                                           /* */ ValueType* w, const bool is_Q_zero = false) {
    using value_type = ValueType;

    /// Given a matrix A that includes QR factorization
    /// it forms a unitary matrix Q
    ///   B = Q = (H0 H1 H2 H3 ... H(k-1)) I
    /// where
    ///   A is m x k (holding H0, H1 ... H(k-1)
    ///   t is k x 1
    ///   B is m x m

    // set identity
    if (is_Q_zero) {
      for (int idx = 0; idx < m; ++idx) {
        Q[(qs0 + qs1) * idx] = value_type(1);
        // SerialSetInternal::invoke(m, value_type(1), Q, qs0 + qs1);
      }
    } else {
      SerialSetIdentityInternal::invoke(m, m, Q, qs0, qs1);
    }

    return SerialApplyQ_LeftForwardInternal::invoke(m, m, k, A, as0, as1, t, ts, Q, qs0, qs1, w);
  }
};

}  // end namespace KokkosBatched

#endif
