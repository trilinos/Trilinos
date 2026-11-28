// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYMM_INTERNAL_HPP_
#define KOKKOSBATCHED_SYMM_INTERNAL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <concepts>
#include <Kokkos_Core.hpp>
#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialSymmInternal {
  template <typename ArgSide, typename ArgUplo, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int m, const int n, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                           const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1);
};

template <typename ArgSide, typename ArgUplo, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSymmInternal::invoke(Op op, SymOp sym_op, const int m, const int n,
                                                      const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                                      const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B,
                                                      const int bs0, const int bs1, const ScalarType beta,
                                                      /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  const ScalarType one(1.0), zero(0.0);

  if (beta == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, C, cs0, cs1);
  else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, beta, C, cs0, cs1);

  if (alpha != zero) {
    if constexpr (std::same_as<ArgSide, Side::Left>) {
      // C := alpha*A*B + beta*C
      if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
        for (int j = 0; j < n; ++j) {
          for (int i = 0; i < m; ++i) {
            ValueType temp1 = alpha * B[i * bs0 + j * bs1];
            ValueType temp2(0);
            for (int k = 0; k < i; ++k) {
              C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
              temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
            }
            C[i * cs0 + j * cs1] += (temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2);
          }
        }
      } else {
        for (int j = 0; j < n; ++j) {
          for (int i = m - 1; i >= 0; --i) {
            ValueType temp1 = alpha * B[i * bs0 + j * bs1];
            ValueType temp2(0);
            for (int k = i + 1; k < m; ++k) {
              C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
              temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
            }
            C[i * cs0 + j * cs1] += (temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2);
          }
        }
      }
    } else {
      // C := alpha*B*A + beta*C
      if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
        for (int j = 0; j < n; ++j) {
          ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + j * bs1];
          }
          for (int k = 0; k < j; ++k) {
            temp1 = alpha * A[k * as0 + j * as1];
            for (int i = 0; i < m; ++i) {
              C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
            }
          }
          for (int k = j + 1; k < n; ++k) {
            temp1 = alpha * op(A[j * as0 + k * as1]);
            for (int i = 0; i < m; ++i) {
              C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
            }
          }
        }
      } else {
        for (int j = 0; j < n; ++j) {
          ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + j * bs1];
          }
          for (int k = 0; k < j; ++k) {
            temp1 = alpha * op(A[j * as0 + k * as1]);
            for (int i = 0; i < m; ++i) {
              C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
            }
          }
          for (int k = j + 1; k < n; ++k) {
            temp1 = alpha * A[k * as0 + j * as1];
            for (int i = 0; i < m; ++i) {
              C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
            }
          }
        }
      }
    }
  }
  return 0;
}

///
/// Team Internal Impl
/// ==================

struct TeamSymmInternal {
  template <typename ArgSide, typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1, const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1);
};

template <typename ArgSide, typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamSymmInternal::invoke(const MemberType &member, Op op, SymOp sym_op, const int m,
                                                    const int n, const ScalarType alpha,
                                                    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                    const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                    const ScalarType beta,
                                                    /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  if constexpr (std::same_as<ArgSide, Side::Left>) {
    // C := alpha*A*B + beta*C
    if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int i) {
          ValueType temp1 = alpha * B[i * bs0 + j * bs1];
          ValueType temp2(0);
          for (int k = 0; k < i; ++k) {
            C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
            temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
          }
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2;
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int ii) {
          const int i     = m - ii - 1;
          ValueType temp1 = alpha * B[i * bs0 + j * bs1];
          ValueType temp2(0);
          for (int k = i + 1; k < m; ++k) {
            C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
            temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
          }
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2;
        });
      });
    }
  } else {
    // C := alpha*B*A + beta*C
    if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int i) {
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * B[i * bs0 + j * bs1];
          for (int k = 0; k < j; ++k) {
            auto temp2 = alpha * A[k * as0 + j * as1];
            C[i * cs0 + j * cs1] += temp2 * B[i * bs0 + k * bs1];
          }
          for (int k = j + 1; k < n; ++k) {
            auto temp2 = alpha * op(A[j * as0 + k * as1]);
            C[i * cs0 + j * cs1] += temp2 * B[i * bs0 + k * bs1];
          }
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
        ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int i) {
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * B[i * bs0 + j * bs1];
          for (int k = 0; k < j; ++k) {
            auto temp2 = alpha * op(A[j * as0 + k * as1]);
            C[i * cs0 + j * cs1] += temp2 * B[i * bs0 + k * bs1];
          }
          for (int k = j + 1; k < n; ++k) {
            auto temp2 = alpha * A[k * as0 + j * as1];
            C[i * cs0 + j * cs1] += temp2 * B[i * bs0 + k * bs1];
          }
        });
      });
    }
  }
  return 0;
}

///
/// TeamVector Internal Impl
/// ========================

struct TeamVectorSymmInternal {
  template <typename ArgSide, typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1, const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1);
};

template <typename ArgSide, typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorSymmInternal::invoke(const MemberType &member, Op op, SymOp sym_op, const int m,
                                                          const int n, const ScalarType alpha,
                                                          const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                          const int as1, const ValueType *KOKKOS_RESTRICT B,
                                                          const int bs0, const int bs1, const ScalarType beta,
                                                          /**/ ValueType *KOKKOS_RESTRICT C, const int cs0,
                                                          const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  if constexpr (std::same_as<ArgSide, Side::Left>) {
    // C := alpha*A*B + beta*C
    if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int i = 0; i < m; ++i) {
          ValueType temp1 = alpha * B[i * bs0 + j * bs1];
          ValueType temp2(0);
          for (int k = 0; k < i; ++k) {
            C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
            temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
          }
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2;
        }
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        for (int i = m - 1; i >= 0; --i) {
          ValueType temp1 = alpha * B[i * bs0 + j * bs1];
          ValueType temp2(0);
          for (int k = i + 1; k < m; ++k) {
            C[k * cs0 + j * cs1] += temp1 * A[k * as0 + i * as1];
            temp2 += B[k * bs0 + j * bs1] * op(A[k * as0 + i * as1]);
          }
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + (temp1 * sym_op(A[i * as0 + i * as1]) + alpha * temp2);
        }
      });
    }
  } else {
    // C := alpha*B*A + beta*C
    if constexpr (std::same_as<ArgUplo, Uplo::Upper>) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
        for (int i = 0; i < m; ++i) {
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * B[i * bs0 + j * bs1];
        }
        for (int k = 0; k < j; ++k) {
          temp1 = alpha * A[k * as0 + j * as1];
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
          }
        }
        for (int k = j + 1; k < n; ++k) {
          temp1 = alpha * op(A[j * as0 + k * as1]);
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
          }
        }
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
        ValueType temp1 = alpha * sym_op(A[j * as0 + j * as1]);
        for (int i = 0; i < m; ++i) {
          C[i * cs0 + j * cs1] = beta * C[i * cs0 + j * cs1] + temp1 * B[i * bs0 + j * bs1];
        }
        for (int k = 0; k < j; ++k) {
          temp1 = alpha * op(A[j * as0 + k * as1]);
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
          }
        }
        for (int k = j + 1; k < n; ++k) {
          temp1 = alpha * A[k * as0 + j * as1];
          for (int i = 0; i < m; ++i) {
            C[i * cs0 + j * cs1] += temp1 * B[i * bs0 + k * bs1];
          }
        }
      });
    }
  }
  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYMM_INTERNAL_HPP_
