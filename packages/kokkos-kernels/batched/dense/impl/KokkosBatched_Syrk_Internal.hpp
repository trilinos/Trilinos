// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYRK_INTERNAL_HPP_
#define KOKKOSBATCHED_SYRK_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

/// Lower

struct SerialSyrkInternalLower {
  template <bool is_trans, typename ReduceType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int n, const int k, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1,
                                           const ScalarType beta, ValueType *KOKKOS_RESTRICT c, const int cs0,
                                           const int cs1);
};

template <bool is_trans, typename ReduceType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSyrkInternalLower::invoke(Op op, SymOp sym_op, const int n, const int k,
                                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a,
                                                           const int as0, const int as1, const ScalarType beta,
                                                           ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
    if constexpr (!is_trans) {
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0)) {
            auto temp            = alpha * op(a[j * as0 + l * as1]);
            c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1] + temp * a[j * as0 + l * as1]);
            for (int i = j + 1; i < n; i++) {
              c[i * cs0 + j * cs1] += temp * a[i * as0 + l * as1];
            }
          }
        }
      }
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      for (int j = 0; j < n; j++) {
        ReduceType rtemp(0);
        for (int l = 0; l < k; l++) {
          rtemp += sym_op(op(a[l * as0 + j * as1]) * a[l * as0 + j * as1]);
        }
        c[j * cs0 + j * cs1] += alpha * rtemp;
        for (int i = j + 1; i < n; i++) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + i * as1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        }
      }
    }
  }

  return 0;
}

/// Upper

struct SerialSyrkInternalUpper {
  template <bool is_trans, typename ReduceType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int n, const int k, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT a, const int as0, const int as1,
                                           const ScalarType beta, ValueType *KOKKOS_RESTRICT c, const int cs0,
                                           const int cs1);
};

template <bool is_trans, typename ReduceType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSyrkInternalUpper::invoke(Op op, SymOp sym_op, const int n, const int k,
                                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a,
                                                           const int as0, const int as1, const ScalarType beta,
                                                           ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
    if constexpr (!is_trans) {
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0)) {
            auto temp = alpha * op(a[j * as0 + l * as1]);
            for (int i = 0; i < j; i++) {
              c[i * cs0 + j * cs1] += temp * a[i * as0 + l * as1];
            }
            c[j * cs0 + j * cs1] += sym_op(temp * a[j * as0 + l * as1]);
          }
        }
      }
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + i * as1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        }
        ReduceType rtemp(0);
        for (int l = 0; l < k; l++) {
          rtemp += sym_op(op(a[l * as0 + j * as1]) * a[l * as0 + j * as1]);
        }
        c[j * cs0 + j * cs1] += alpha * rtemp;
      }
    }
  }

  return 0;
}

///
/// Team Internal Impl
/// ====================

/// Lower

struct TeamSyrkInternalLower {
  template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ScalarType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamSyrkInternalLower::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                         const int k, const ScalarType alpha,
                                                         const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                         const int as1, const ScalarType beta,
                                                         ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    if constexpr (!is_trans) {
      // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j)), [&](const int ii) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[j * as0 + l * as1]) * a[(ii + j) * as0 + l * as1];
          }
          c[(ii + j) * cs0 + j * cs1] += alpha * temp;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j)), [&](const int ii) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + (ii + j) * as1]) * a[l * as0 + j * as1];
          }
          c[(ii + j) * cs0 + j * cs1] += alpha * temp;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}

/// Upper

struct TeamSyrkInternalUpper {
  template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ScalarType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamSyrkInternalUpper::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                         const int k, const ScalarType alpha,
                                                         const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                         const int as1, const ScalarType beta,
                                                         ValueType *KOKKOS_RESTRICT c, const int cs0, const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1), [&](const int i) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[j * as0 + l * as1]) * a[i * as0 + l * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1), [&](const int i) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + i * as1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        });
        c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1]);
      });
    }
  }

  return 0;
}

///
/// TeamVector Internal Impl
/// ====================

/// Lower

struct TeamVectorSyrkInternalLower {
  template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ScalarType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorSyrkInternalLower::invoke(const MemberType &member, Op op, SymOp sym_op,
                                                               const int n, const int k, const ScalarType alpha,
                                                               const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                               const int as1, const ScalarType beta,
                                                               ValueType *KOKKOS_RESTRICT c, const int cs0,
                                                               const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    if constexpr (!is_trans) {
      // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0)) {
            auto temp            = alpha * op(a[j * as0 + l * as1]);
            c[j * cs0 + j * cs1] = sym_op(c[j * cs0 + j * cs1] + temp * a[j * as0 + l * as1]);
            for (int i = j + 1; i < n; i++) {
              c[i * cs0 + j * cs1] += temp * a[i * as0 + l * as1];
            }
          }
        }
      });
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        ReduceType rtemp(0);
        for (int l = 0; l < k; l++) {
          rtemp += sym_op(op(a[l * as0 + j * as1]) * a[l * as0 + j * as1]);
        }
        c[j * cs0 + j * cs1] += alpha * rtemp;
        for (int i = j + 1; i < n; i++) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + i * as1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        }
      });
    }
  }

  return 0;
}

/// Upper

struct TeamVectorSyrkInternalUpper {
  template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT a, const int as0,
                                           const int as1, const ScalarType beta, ValueType *KOKKOS_RESTRICT c,
                                           const int cs0, const int cs1);
};

template <bool is_trans, typename ReduceType, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorSyrkInternalUpper::invoke(const MemberType &member, Op op, SymOp sym_op,
                                                               const int n, const int k, const ScalarType alpha,
                                                               const ValueType *KOKKOS_RESTRICT a, const int as0,
                                                               const int as1, const ScalarType beta,
                                                               ValueType *KOKKOS_RESTRICT c, const int cs0,
                                                               const int cs1) {
  if (beta == ScalarType(0)) {
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

  if (alpha != ScalarType(0)) {
    // C: = alpha * A * A**T + beta * C or C: = alpha * A * A**H + beta * C
    if constexpr (!is_trans) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int l = 0; l < k; l++) {
          if (a[j * as0 + l * as1] != ValueType(0)) {
            auto temp = alpha * op(a[j * as0 + l * as1]);
            for (int i = 0; i < j; i++) {
              c[i * cs0 + j * cs1] += temp * a[i * as0 + l * as1];
            }
            c[j * cs0 + j * cs1] += sym_op(temp * a[j * as0 + l * as1]);
          }
        }
      });
    } else {
      // C: = alpha * A**T * A + beta * C or C: = alpha * A**H * A + beta * C
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int i = 0; i < j; i++) {
          ValueType temp(0);
          for (int l = 0; l < k; l++) {
            temp += op(a[l * as0 + i * as1]) * a[l * as0 + j * as1];
          }
          c[i * cs0 + j * cs1] += alpha * temp;
        }
        ReduceType rtemp(0);
        for (int l = 0; l < k; l++) {
          rtemp += sym_op(op(a[l * as0 + j * as1]) * a[l * as0 + j * as1]);
        }
        c[j * cs0 + j * cs1] += alpha * rtemp;
      });
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYRK_INTERNAL_HPP_
