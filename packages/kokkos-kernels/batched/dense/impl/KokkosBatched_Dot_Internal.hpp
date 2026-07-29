// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_DOT_INTERNAL_HPP
#define KOKKOSBATCHED_DOT_INTERNAL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialDotInternal {
  // i \in [0,m)
  // C = conj(A(:))*B(:)
  template <typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(Op op, const int m, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    C[0] = ValueType(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      const int idx_a = i * as0, idx_b = i * bs0;
      C[0] += op(A[idx_a]) * B[idx_b];
    }
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, const int m, const int n, const ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B,
                                           const int bs0, const int bs1,
                                           /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    for (int j = 0; j < n; ++j) invoke(op, m, A + j * as1, as0, B + j * bs1, bs0, C + j * cs);
    return 0;
  }
};

///
/// Team Internal Impl
/// ========================

// i \in [0,m)
// C = conj(A(:))*B(:)
struct TeamDotInternal {
  template <typename MemberType, typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    ValueType t(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, m),
        [&](const int &i, ValueType &update) {
          const int idx_a = i * as0, idx_b = i * bs0;
          update += op(A[idx_a]) * B[idx_b];
        },
        t);
    Kokkos::single(Kokkos::PerThread(member), [&]() { C[0] = t; });
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename MemberType, typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
      const ValueType *KOKKOS_RESTRICT B_at_j = B + j * bs1;
      for (int i = 0; i < m; ++i) {
        const int idx_a = i * as0, idx_b = i * bs0;
        t += op(A_at_j[idx_a]) * B_at_j[idx_b];
      }
      Kokkos::single(Kokkos::PerThread(member), [&]() { C[j * cs] = t; });
    });
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================

// i \in [0,m)
// C = conj(A(:))*B(:)
struct TeamVectorDotInternal {
  template <typename MemberType, typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    ValueType t(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, m),
        [&](const int &i, ValueType &update) {
          const int idx_a = i * as0, idx_b = i * bs0;
          update += op(A[idx_a]) * B[idx_b];
        },
        t);
    Kokkos::single(Kokkos::PerThread(member), [&]() { C[0] = t; });
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename MemberType, typename Op, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
      const ValueType *KOKKOS_RESTRICT B_at_j = B + j * bs1;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(member, m),
          [&](const int &i, ValueType &update) {
            const int idx_a = i * as0, idx_b = i * bs0;
            update += op(A_at_j[idx_a]) * B_at_j[idx_b];
          },
          t);
      Kokkos::single(Kokkos::PerThread(member), [&]() { C[j * cs] = t; });
    });
    return 0;
  }
};
}  // namespace Impl

struct [[deprecated("Use KokkosBatched::SerialDot instead")]] SerialDotInternal {
  template <typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int m, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    return Impl::SerialDotInternal::invoke(KokkosBlas::Impl::OpConj(), m, A, as0, B, bs0, C);
  }

  template <typename ValueType, typename MagnitudeType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1,
                                           /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    return Impl::SerialDotInternal::invoke(KokkosBlas::Impl::OpConj(), m, n, A, as0, as1, B, bs0, bs1, C, cs);
  }
};

struct [[deprecated("Use KokkosBatched::TeamDot instead")]] TeamDotInternal {
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    return Impl::TeamDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), m, A, as0, B, bs0, C);
  }

  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    return Impl::TeamDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), m, n, A, as0, as1, B, bs0, bs1, C, cs);
  }
};

struct [[deprecated("Use KokkosBatched::TeamVectorDot instead")]] TeamVectorDotInternal {
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    return Impl::TeamVectorDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), m, A, as0, B, bs0, C);
  }

  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    return Impl::TeamVectorDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), m, n, A, as0, as1, B, bs0, bs1, C,
                                               cs);
  }
};

}  // namespace KokkosBatched

#endif
