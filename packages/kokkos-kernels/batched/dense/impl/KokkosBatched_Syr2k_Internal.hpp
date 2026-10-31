// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYR2K_INTERNAL_HPP_
#define KOKKOSBATCHED_SYR2K_INTERNAL_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

/// Lower

struct SerialSyr2kInternalLower {
  template <bool is_trans, typename Op, typename SymOp, typename AlphaType, typename ValueType, typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha,
                                           const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1,
                                           const ValueType *KOKKOS_RESTRICT b, const int bs0, const int bs1,
                                           const BetaType beta, ValueType *KOKKOS_RESTRICT c, const int cs0,
                                           const int cs1);
};

template <bool is_trans, typename Op, typename SymOp, typename AlphaType, typename ValueType, typename BetaType>
KOKKOS_INLINE_FUNCTION int SerialSyr2kInternalLower::invoke(
    Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a,
    const int as0, const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0, const int bs1, const BetaType beta,
    ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    for (int j = 0; j < n; j++) {
      for (int i = j; i < n; i++) {
        c[i * cs0 + j * cs1] = ValueType(0);
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
      for (int i = j + 1; i < n; i++) {
        c[i * cs0 + j * cs1] *= beta;
      }
    }
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha*A*B^T + alpha*B*A^T + beta*C or C: = alpha*A*B^H + conj(alpha)*B*A^H + beta*C
    if constexpr (!is_trans) {
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0) || b[j * bs0 + l * bs1] != ValueType(0)) {
            auto temp1 = alpha * op(b[j * bs0 + l * bs1]);
            auto temp2 = op(alpha * a[j * as0 + l * as1]);
            c[j * cs0 + j * cs1] =
                sym_op(c[j * cs0 + j * cs1] + temp1 * a[j * as0 + l * as1] + temp2 * b[j * bs0 + l * bs1]);
            for (int i = j + 1; i < n; i++) {
              c[i * cs0 + j * cs1] += (temp1 * a[i * as0 + l * as1] + temp2 * b[i * bs0 + l * bs1]);
            }
          }
        }
      }
    } else {
      // C: = alpha*A^T*B + alpha*B^T*A + beta*C or C: = alpha*A^H*B + conj(alpha)*B^H*A + beta*C
      for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        }
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      }
    }
  }

  return 0;
}

/// Upper

struct SerialSyr2kInternalUpper {
  template <bool is_trans, typename Op, typename SymOp, typename AlphaType, typename ValueType, typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha,
                                           const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1,
                                           const ValueType *KOKKOS_RESTRICT b, const int bs0, const int bs1,
                                           const BetaType beta, ValueType *KOKKOS_RESTRICT c, const int cs0,
                                           const int cs1);
};

template <bool is_trans, typename Op, typename SymOp, typename AlphaType, typename ValueType, typename BetaType>
KOKKOS_INLINE_FUNCTION int SerialSyr2kInternalUpper::invoke(
    Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a,
    const int as0, const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0, const int bs1, const BetaType beta,
    ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < j + 1; i++) {
        c[i * cs0 + j * cs1] = ValueType(0);
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < j; i++) {
        c[i * cs0 + j * cs1] *= beta;
      }
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
    }
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha * A * B^T + alpha * B * A^T + beta * C or C: = alpha * A * B^H + conj(alpha) * B * A^H + beta * C
    if constexpr (!is_trans) {
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0) || b[j * bs0 + l * bs1] != ValueType(0)) {
            auto temp1 = alpha * op(b[j * bs0 + l * bs1]);
            auto temp2 = op(alpha * a[j * as0 + l * as1]);
            for (int i = 0; i < j; i++) {
              c[i * cs0 + j * cs1] += (temp1 * a[i * as0 + l * as1] + temp2 * b[i * bs0 + l * bs1]);
            }
            c[j * cs0 + j * cs1] += sym_op(temp1 * a[j * as0 + l * as1] + temp2 * b[j * bs0 + l * bs1]);
          }
        }
      }
    } else {
      // C: = alpha * A^T * B + alpha * B^T * A + beta * C or C: = alpha * A^H * B + conj(alpha) * B^H * A + beta * C
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < j + 1; i++) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        }
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      }
    }
  }

  return 0;
}

///
/// Team Internal Impl
/// ====================

/// Lower

struct TeamSyr2kInternalLower {
  template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
            typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
                                           const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
          typename BetaType>
KOKKOS_INLINE_FUNCTION int TeamSyr2kInternalLower::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                          const int k, const AlphaType alpha,
                                                          const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                          const int as1, const ValueType *KOKKOS_RESTRICT b,
                                                          const int bs0, const int bs1, const BetaType beta,
                                                          ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j)),
                           [&](const int ii) { c[(ii + j) * cs0 + j * cs1] = ValueType(0); });
    });
  } else {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j - 1)),
                           [&](const int ii) { c[(ii + j + 1) * cs0 + j * cs1] *= beta; });
    });
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha*A*B^T + alpha*B*A^T + beta*C or C: = alpha*A*B^H + conj(alpha)*B*A^H + beta*C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j)), [&](const int ii) {
          const int i = ii + j;
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += a[i * as0 + l * as1] * op(b[j * bs0 + l * bs1]);
            temp2 += b[i * bs0 + l * bs1] * op(a[j * as0 + l * as1]);
          }
          c[i * cs0 + j * cs1] += alpha * temp1 + op(alpha) * temp2;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    } else {
      // C: = alpha*A^T*B + alpha*B^T*A + beta*C or C: = alpha*A^H*B + conj(alpha)*B^H*A + beta*C
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j)), [&](const int ii) {
          const int i = ii + j;
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}

/// Upper

struct TeamSyr2kInternalUpper {
  template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
            typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
                                           const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
          typename BetaType>
KOKKOS_INLINE_FUNCTION int TeamSyr2kInternalUpper::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                          const int k, const AlphaType alpha,
                                                          const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                          const int as1, const ValueType *KOKKOS_RESTRICT b,
                                                          const int bs0, const int bs1, const BetaType beta,
                                                          ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1),
                           [&](const int i) { c[i * cs0 + j * cs1] = ValueType(0); });
    });
  } else {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j), [&](const int i) { c[i * cs0 + j * cs1] *= beta; });
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
    });
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha * A * B^T + alpha * B * A^T + beta * C or C: = alpha * A * B^H + conj(alpha) * B * A^H + beta * C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1), [&](const int i) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += a[i * as0 + l * as1] * op(b[j * bs0 + l * bs1]);
            temp2 += b[i * bs0 + l * bs1] * op(a[j * as0 + l * as1]);
          }
          c[i * cs0 + j * cs1] += alpha * temp1 + op(alpha) * temp2;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    } else {
      // C: = alpha * A^T * B + alpha * B^T * A + beta * C or C: = alpha * A^H * B + conj(alpha) * B^H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1), [&](const int i) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}

///
/// TeamVector Internal Impl
/// ========================

/// Lower

struct TeamVectorSyr2kInternalLower {
  template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
            typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
                                           const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
          typename BetaType>
KOKKOS_INLINE_FUNCTION int TeamVectorSyr2kInternalLower::invoke(
    const MemberType &member, Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha,
    const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
    const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      for (int i = j; i < n; i++) {
        c[i * cs0 + j * cs1] = ValueType(0);
      }
    });
  } else {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
      for (int i = j + 1; i < n; i++) {
        c[i * cs0 + j * cs1] *= beta;
      }
    });
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha*A*B^T + alpha*B*A^T + beta*C or C: = alpha*A*B^H + conj(alpha)*B*A^H + beta*C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0) || b[j * bs0 + l * bs1] != ValueType(0)) {
            auto temp1 = alpha * op(b[j * bs0 + l * bs1]);
            auto temp2 = op(alpha * a[j * as0 + l * as1]);
            c[j * cs0 + j * cs1] =
                sym_op(c[j * cs0 + j * cs1] + temp1 * a[j * as0 + l * as1] + temp2 * b[j * bs0 + l * bs1]);
            for (int i = j + 1; i < n; i++) {
              c[i * cs0 + j * cs1] += (temp1 * a[i * as0 + l * as1] + temp2 * b[i * bs0 + l * bs1]);
            }
          }
        }
      });
    } else {
      // C: = alpha*A^T*B + alpha*B^T*A + beta*C or C: = alpha*A^H*B + conj(alpha)*B^H*A + beta*C
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int i = j; i < n; i++) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        }
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}

/// Upper

struct TeamVectorSyr2kInternalUpper {
  template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
            typename BetaType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const AlphaType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
                                           const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename MemberType, typename Op, typename SymOp, typename AlphaType, typename ValueType,
          typename BetaType>
KOKKOS_INLINE_FUNCTION int TeamVectorSyr2kInternalUpper::invoke(
    const MemberType &member, Op op, SymOp sym_op, const int n, const int k, const AlphaType alpha,
    const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1, const ValueType *KOKKOS_RESTRICT b, const int bs0,
    const int bs1, const BetaType beta, ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == BetaType(0)) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      for (int i = 0; i < j + 1; i++) {
        c[i * cs0 + j * cs1] = ValueType(0);
      }
    });
  } else {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      for (int i = 0; i < j; i++) {
        c[i * cs0 + j * cs1] *= beta;
      }
      c[j * cs0 + j * cs1] = beta * sym_op(c[j * cs0 + j * cs1]);
    });
  }

  if (alpha != AlphaType(0)) {
    // C: = alpha * A * B^T + alpha * B * A^T + beta * C or C: = alpha * A * B^H + conj(alpha) * B * A^H + beta * C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0) || b[j * bs0 + l * bs1] != ValueType(0)) {
            auto temp1 = alpha * op(b[j * bs0 + l * bs1]);
            auto temp2 = op(alpha * a[j * as0 + l * as1]);
            for (int i = 0; i < j; i++) {
              c[i * cs0 + j * cs1] += (temp1 * a[i * as0 + l * as1] + temp2 * b[i * bs0 + l * bs1]);
            }
            c[j * cs0 + j * cs1] += sym_op(temp1 * a[j * as0 + l * as1] + temp2 * b[j * bs0 + l * bs1]);
          }
        }
      });
    } else {
      // C: = alpha * A^T * B + alpha * B^T * A + beta * C or C: = alpha * A^H * B + conj(alpha) * B^H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int i = 0; i < j + 1; i++) {
          ValueType temp1(0), temp2(0);
          for (int l = 0; l < k; l++) {
            temp1 += op(a[l * as0 + i * as1]) * b[l * bs0 + j * bs1];
            temp2 += op(b[l * bs0 + i * bs1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += (alpha * temp1 + op(alpha) * temp2);
        }
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}
}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYR2K_INTERNAL_HPP_
