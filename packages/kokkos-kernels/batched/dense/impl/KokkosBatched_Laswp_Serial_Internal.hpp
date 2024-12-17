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

#ifndef KOKKOSBATCHED_LASWP_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_LASWP_SERIAL_INTERNAL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ========================

///
//// Forward pivot apply
///

struct SerialLaswpVectorForwardInternal {
  template <typename IntType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int plen, const IntType *KOKKOS_RESTRICT p, const int ps0,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    for (int i = 0; i < plen; ++i) {
      const int piv = p[i * ps0];
      if (piv != i) {
        const int idx_i = i * as0, idx_p = piv * as0;
        const ValueType tmp = A[idx_i];
        A[idx_i]            = A[idx_p];
        A[idx_p]            = tmp;
      }
    }
    return 0;
  }
};

struct SerialLaswpMatrixForwardInternal {
  template <typename IntType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const int plen, const IntType *KOKKOS_RESTRICT p, const int ps0,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (as0 <= as1) {
      // LayoutLeft like
      for (int j = 0; j < n; j++) {
        ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
        for (int i = 0; i < plen; ++i) {
          const int piv = p[i * ps0];
          if (piv != i) {
            const int idx_i = i * as0, idx_p = piv * as0;
            const ValueType tmp = A_at_j[idx_i];
            A_at_j[idx_i]       = A_at_j[idx_p];
            A_at_j[idx_p]       = tmp;
          }
        }
      }
    } else {
      // LayoutRight like
      for (int i = 0; i < plen; ++i) {
        const int piv = p[i * ps0];
        if (piv != i) {
          const int idx_i = i * as0, idx_p = piv * as0;
          for (int j = 0; j < n; j++) {
            ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
            const ValueType tmp               = A_at_j[idx_i];
            A_at_j[idx_i]                     = A_at_j[idx_p];
            A_at_j[idx_p]                     = tmp;
          }
        }
      }
    }
    return 0;
  }
};

///
/// Backward pivot apply
///

struct SerialLaswpVectorBackwardInternal {
  template <typename IntType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int plen, const IntType *KOKKOS_RESTRICT p, const int ps0,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    for (int i = (plen - 1); i >= 0; --i) {
      const int piv = p[i * ps0];
      if (piv != i) {
        const int idx_i = i * as0, idx_p = piv * as0;
        const ValueType tmp = A[idx_i];
        A[idx_i]            = A[idx_p];
        A[idx_p]            = tmp;
      }
    }
    return 0;
  }
};

struct SerialLaswpMatrixBackwardInternal {
  template <typename IntType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const int plen, const IntType *KOKKOS_RESTRICT p, const int ps0,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (as0 <= as1) {
      // LayoutLeft like
      for (int j = 0; j < n; j++) {
        ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
        for (int i = (plen - 1); i >= 0; --i) {
          const int piv = p[i * ps0];
          if (piv != i) {
            const int idx_i = i * as0, idx_p = piv * as0;
            const ValueType tmp = A_at_j[idx_i];
            A_at_j[idx_i]       = A_at_j[idx_p];
            A_at_j[idx_p]       = tmp;
          }
        }
      }
    } else {
      // LayoutRight like
      for (int i = (plen - 1); i >= 0; --i) {
        const int piv = p[i * ps0];
        if (piv != i) {
          const int idx_i = i * as0, idx_p = piv * as0;
          for (int j = 0; j < n; j++) {
            ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
            const ValueType tmp               = A_at_j[idx_i];
            A_at_j[idx_i]                     = A_at_j[idx_p];
            A_at_j[idx_p]                     = tmp;
          }
        }
      }
    }
    return 0;
  }
};

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_LASWP_SERIAL_INTERNAL_HPP_
