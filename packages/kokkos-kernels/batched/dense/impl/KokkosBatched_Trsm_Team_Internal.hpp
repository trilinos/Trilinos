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
#ifndef __KOKKOSBATCHED_TRSM_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSM_TEAM_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"
#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Internal Impl
/// ====================

template <typename AlgoType>
struct TeamTrsmInternalLeftLower {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const bool use_unit_diag, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      // Made this non-const in order to WORKAROUND issue #349
      int iend = m - p - 1;
      int jend = n;

      const ValueType *KOKKOS_RESTRICT a21 = iend ? A + (p + 1) * as0 + p * as1 : NULL;

      ValueType *KOKKOS_RESTRICT b1t = B + p * bs0, *KOKKOS_RESTRICT B2 = iend ? B + (p + 1) * bs0 : NULL;

      member.team_barrier();
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, jend),
                             [&](const int &j) { b1t[j * bs1] = b1t[j * bs1] / alpha11; });
        member.team_barrier();
      }
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, iend * jend), [&](const int &ij) {
        // assume layout right for batched computation
        const int i = ij / jend, j = ij % jend;
        B2[i * bs0 + j * bs1] -= a21[i * as0] * b1t[j * bs1];
      });
    }
  }
  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  constexpr int mbAlgo = Algo::Trsm::Blocked::mb();

  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  // note that parallel range is different ( m*n vs m-1*n);
  if (alpha == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    ///
    /// case host: team size is small and blocksize (mb,nb) is large

    ///
    /// case GPU: team size is large and blocksize (mb,nb) is small
    InnerTrsmLeftLowerUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, bs1);
    InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);

    auto trsm = [&](const int ib, const int jb, const ValueType *KOKKOS_RESTRICT AA,
                    /**/ ValueType *KOKKOS_RESTRICT BB) {
      const int mb    = mbAlgo;
      const int tsize = member.team_size();
      // Made this non-const in order to WORKAROUND issue #349
      int nb = (jb / tsize + jb % tsize > 0);
      int np = jb % nb;
      for (int p = 0; p < ib; p += mb) {
        // Made this non-const in order to WORKAROUND issue #349
        int pb = ((p + mb) > ib ? (ib - p) : mb);

        // trsm update
        const ValueType *KOKKOS_RESTRICT Ap = AA + p * as0 + p * as1;
        /**/ ValueType *KOKKOS_RESTRICT Bp  = BB + p * bs0;

        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, (jb / nb) + (np > 0)), [&](const int jj) {
          // Made this non-const in order to WORKAROUND issue #349
          int j = jj * nb, qb = (j + nb) > jb ? np : nb;
          if (use_unit_diag)
            trsm_u.serial_invoke(Ap, pb, qb, Bp + j * bs1);
          else
            trsm_n.serial_invoke(Ap, pb, qb, Bp + j * bs1);
        });
        member.team_barrier();

        // gemm update
        TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, ib - p - pb, jb, pb, minus_one, Ap + pb * as0, as0, as1,
                                                      Bp, bs0, bs1, one, Bp + pb * bs0, bs0, bs1);
      }
    };

    const bool is_small = true;  //(m*n <= 64*64);
    if (is_small) {
      trsm(m, n, A, B);
    } else {
      // // some cache blocking may need (not priority yet);
      // trsm(m, n, A, B);
    }
  }
  return 0;
}

template <typename AlgoType>
struct TeamTrsmInternalLeftUpper {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const bool use_unit_diag, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  // note that parallel range is different ( m*n vs m-1*n);
  if (alpha == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    ValueType *KOKKOS_RESTRICT B0 = B;
    for (int p = (m - 1); p >= 0; --p) {
      // Made this non-const in order to WORKAROUND issue #349
      int iend = p;
      int jend = n;

      const ValueType *KOKKOS_RESTRICT a01 = A + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT b1t  = B + p * bs0;

      member.team_barrier();
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, jend),
                             [&](const int &j) { b1t[j * bs1] = b1t[j * bs1] / alpha11; });
        member.team_barrier();
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, iend * jend), [&](const int &ij) {
        int i, j;
        if (KokkosKernels::Impl::kk_is_gpu_exec_space<typename MemberType::execution_space>()) {
          i = ij % iend;
          j = ij / iend;
        } else {
          i = ij / jend;
          j = ij % jend;
        }
        B0[i * bs0 + j * bs1] -= a01[i * as0] * b1t[j * bs1];
      });
    }
  }
  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  constexpr int mbAlgo = Algo::Trsm::Blocked::mb();

  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  // note that parallel range is different ( m*n vs m-1*n);
  if (alpha == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    InnerTrsmLeftUpperUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, bs1);
    InnerTrsmLeftUpperNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);

    auto trsm = [&](const int ib, const int jb, const ValueType *KOKKOS_RESTRICT AA,
                    /**/ ValueType *KOKKOS_RESTRICT BB) {
      const int mb    = mbAlgo;  //(ib <=5 ? ib : mbAlgo);
      const int tsize = member.team_size();
      // Made this non-const in order to WORKAROUND issue #349
      int nb = (jb / tsize + jb % tsize > 0);
      int np = jb % nb;
      for (int pp = 0; pp < ib; pp += mb) {
        const int ptmp = (ib - pp - mb), p = (ptmp < 0 ? 0 : ptmp), pb = (mb + (ptmp < 0) * ptmp);

        // trsm update
        const ValueType *KOKKOS_RESTRICT Ap = AA + p * as0 + p * as1;
        /**/ ValueType *KOKKOS_RESTRICT Bp  = BB + p * bs0;

        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, (jb / nb) + (np > 0)), [&](const int &jj) {
          const int j = jj * nb, qb = (j + nb) > jb ? np : nb;
          if (use_unit_diag)
            trsm_u.serial_invoke(Ap, pb, qb, Bp + j * bs1);
          else
            trsm_n.serial_invoke(Ap, pb, qb, Bp + j * bs1);
        });
        member.team_barrier();

        // gemm update
        TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, p, jb, pb, minus_one, Ap - p * as0, as0, as1, Bp, bs0,
                                                      bs1, one, BB, bs0, bs1);
      }
    };

    const bool is_small = true;  //(m*n <= 64*64);
    if (is_small) {
      trsm(m, n, A, B);
    } else {
      // // some cache blocking may need (not priority yet);
      // trsm(m, n, A, B);
    }
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
