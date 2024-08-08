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
#ifndef __KOKKOSBATCHED_QR_FORM_Q_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_QR_FORM_Q_SERIAL_INTERNAL_HPP__

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
    typedef ValueType value_type;

    /// Given a matrix A that includes QR factorization
    /// it forms a unitary matrix Q
    ///   B = Q = (H0 H1 H2 H3 ... H(k-1)) I
    /// where
    ///   A is m x k (holding H0, H1 ... H(k-1)
    ///   t is k x 1
    ///   B is m x m

    // set identity
    if (is_Q_zero)
      SerialSetInternal::invoke(m, value_type(1), Q, qs0 + qs1);
    else
      SerialSetIdentityInternal::invoke(m, Q, qs0, qs1);

    return SerialApplyQ_LeftNoTransForwardInternal ::invoke(m, m, k, A, as0, as1, t, ts, Q, qs0, qs1, w);
  }
};

}  // end namespace KokkosBatched

#endif
