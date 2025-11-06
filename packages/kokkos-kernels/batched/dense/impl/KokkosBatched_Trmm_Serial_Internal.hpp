// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRMM_SERIAL_INTERNAL_HPP
#define KOKKOSBATCHED_TRMM_SERIAL_INTERNAL_HPP

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"

namespace KokkosBatched {

template <typename AlgoType>
struct SerialTrmmInternalLeftLower {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int am, const int an,
                                           const int bm, const int bn, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <typename AlgoType>
struct SerialTrmmInternalLeftUpper {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int am, const int an,
                                           const int bm, const int bn, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <typename AlgoType>
struct SerialTrmmInternalRightLower {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int am, const int an,
                                           const int bm, const int bn, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <typename AlgoType>
struct SerialTrmmInternalRightUpper {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int am, const int an,
                                           const int bm, const int bn, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

// ech-note: use_unit_diag intentionally ignored for now. Compiler can optimize
// it out. Assuming that branching logic (especially on GPU) to handle
// use_unit_diag will use more cycles than simply doing 1.0*B[idx] for the copy
// if use_unit_diag.
template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
    const bool /*use_unit_diag*/, const bool do_conj, const int am, const int an, const int bm, const int bn,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);
  typedef KokkosKernels::ArithTraits<ValueType> AT;
  int left_m  = am;
  int right_n = bn;
  // echo-TODO: See about coniditionally setting conjOp at compile time.
  // auto conjOp = noop;
  // if (do_conj) {
  //  conjOp = AT::conj;
  //}
  // printf("SerialTrmmInternalLeftLower\n");

  auto dotLowerLeftConj = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int left_row,
                              ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1, const int right_col) {
    auto B_elems   = left_row;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += _a[left_row, i] * B[i, right_col]
      sum += AT::conj(_a[left_row * _as0 + i * _as1]) * _b[i * _bs0 + _bs1 * right_col];
    }
    return sum;
  };

  auto dotLowerLeft = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int left_row,
                          ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1, const int right_col) {
    auto B_elems   = left_row;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += _a[left_row, i] * _b[i, right_col]
      sum += _a[left_row * _as0 + i * _as1] * _b[i * _bs0 + _bs1 * right_col];
    }
    return sum;
  };

  if (bm <= 0 || bn <= 0 || am <= 0 || an <= 0) return 0;

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(bm, bn, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(bm, bn, alpha, B, bs0, bs1);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int m = left_m - 1; m >= 0; m--) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int n = 0; n < right_n; n++) {
        if (do_conj) {
          B[m * bs0 + n * bs1] = dotLowerLeftConj(A, as0, as1, m, B, bs0, bs1, n);
        } else {
          B[m * bs0 + n * bs1] = dotLowerLeft(A, as0, as1, m, B, bs0, bs1, n);
        }
      }
    }
  }
  return 0;
}

// ech-note: use_unit_diag intentionally ignored for now. Compiler can optimize
// it out. Assuming that branching logic (especially on GPU) to handle
// use_unit_diag will use more cycles than simply doing 1.0*B[idx] for the copy
// if use_unit_diag.
template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
    const bool /*use_unit_diag*/, const bool do_conj, const int am, const int an, const int bm, const int bn,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);
  typedef KokkosKernels::ArithTraits<ValueType> AT;
  int left_m  = bm;
  int right_n = an;
  // echo-TODO: See about coniditionally setting conjOp at compile time.
  // auto conjOp = noop;
  // if (do_conj) {
  //  conjOp = AT::conj;
  //}

  // Lower triangular matrix is on RHS with the base facing down.
  // Everytime we compute a new output row of B, we must shift over to the
  // right by one in A's column to ensure we skip the 0's.
  auto dotLowerRightConj = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int _am,
                               const int left_row, ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1,
                               const int right_col) {
    auto B_elems   = _am - 1;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = right_col; i <= B_elems; i++) {
      // sum += B[left_row, i] * A[i, right_col]
      sum += _b[_bs0 * left_row + i * _bs1] * AT::conj(_a[i * _as0 + right_col * _as1]);
    }
    return sum;
  };

  auto dotLowerRight = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int _am,
                           const int left_row, ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1,
                           const int right_col) {
    auto B_elems   = _am - 1;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = right_col; i <= B_elems; i++) {
      // sum += B[left_row, i] * A[i, right_col]
      sum += _b[_bs0 * left_row + i * _bs1] * _a[i * _as0 + right_col * _as1];
    }
    return sum;
  };

  if (bm <= 0 || bn <= 0 || am <= 0 || an <= 0) return 0;

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(bm, bn, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(bm, bn, alpha, B, bs0, bs1);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int m = 0; m < left_m; m++) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int n = 0; n < right_n; n++) {
        if (do_conj) {
          B[m * bs0 + n * bs1] = dotLowerRightConj(A, as0, as1, am, m, B, bs0, bs1, n);
        } else {
          B[m * bs0 + n * bs1] = dotLowerRight(A, as0, as1, am, m, B, bs0, bs1, n);
        }
      }
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
    const bool /*use_unit_diag*/, const bool do_conj, const int am, const int an, const int bm, const int bn,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);
  typedef KokkosKernels::ArithTraits<ValueType> AT;
  int left_m  = am;
  int right_n = bn;
  // echo-TODO: See about coniditionally setting conjOp at compile time.
  // auto conjOp = noop;
  // if (do_conj) {
  //  conjOp = AT::conj;
  //}

  auto dotUpperLeftConj = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int _an,
                              const int left_row, ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1,
                              const int right_col) {
    auto B_elems   = _an - left_row - 1;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += A[left_row, i+left_row] * B[i+left_row, right_col]
      sum += AT::conj(_a[left_row * _as0 + (i + left_row) * _as1]) * _b[(i + left_row) * _bs0 + _bs1 * right_col];
    }
    return sum;
  };

  auto dotUpperLeft = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int _an,
                          const int left_row, ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1,
                          const int right_col) {
    auto B_elems   = _an - left_row - 1;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += A[left_row, i+left_row] * B[i+left_row, right_col]
      sum += _a[left_row * _as0 + (i + left_row) * _as1] * _b[(i + left_row) * _bs0 + _bs1 * right_col];
    }
    return sum;
  };

  if (bm <= 0 || bn <= 0 || am <= 0 || an <= 0) return 0;

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(bm, bn, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(bm, bn, alpha, B, bs0, bs1);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int m = 0; m < left_m; ++m) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int n = 0; n < right_n; ++n) {
        if (do_conj) {
          B[m * bs0 + n * bs1] = dotUpperLeftConj(A, as0, as1, an, m, B, bs0, bs1, n);
        } else {
          B[m * bs0 + n * bs1] = dotUpperLeft(A, as0, as1, an, m, B, bs0, bs1, n);
        }
      }
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
    const bool /*use_unit_diag*/, const bool do_conj, const int am, const int an, const int bm, const int bn,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);
  typedef KokkosKernels::ArithTraits<ValueType> AT;
  int left_m  = bm;
  int right_n = an;
  // echo-TODO: See about coniditionally setting conjOp at compile time.
  // auto conjOp = noop;
  // if (do_conj) {
  //  conjOp = AT::conj;
  //}

  auto dotUpperRightConj = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int left_row,
                               ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1, const int right_col) {
    auto B_elems   = right_col;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += B[left_row, i] * A[i, right_col]
      sum += _b[left_row * _bs0 + i * _bs1] * AT::conj(_a[i * _as0 + right_col * _as1]);
    }
    return sum;
  };

  auto dotUpperRight = [&](const ValueType *KOKKOS_RESTRICT _a, const int _as0, const int _as1, const int left_row,
                           ValueType *KOKKOS_RESTRICT _b, const int _bs0, const int _bs1, const int right_col) {
    auto B_elems   = right_col;
    ScalarType sum = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i <= B_elems; i++) {
      // sum += B[left_row, i] * A[i, right_col]
      sum += _b[left_row * _bs0 + i * _bs1] * _a[i * _as0 + right_col * _as1];
    }
    return sum;
  };

  if (bm <= 0 || bn <= 0 || am <= 0 || an <= 0) return 0;

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(bm, bn, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(bm, bn, alpha, B, bs0, bs1);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int m = 0; m < left_m; ++m) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int n = right_n - 1; n >= 0; --n) {
        if (do_conj) {
          B[m * bs0 + n * bs1] = dotUpperRightConj(A, as0, as1, m, B, bs0, bs1, n);
        } else {
          B[m * bs0 + n * bs1] = dotUpperRight(A, as0, as1, m, B, bs0, bs1, n);
        }
      }
    }
  }
  return 0;
}
}  // namespace KokkosBatched
#endif  // KOKKOSBATCHED_TRMM_SERIAL_INTERNAL_HPP
