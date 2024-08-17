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
#ifndef __KOKKOSBATCHED_TRSV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSV_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"

namespace KokkosBatched {

///
/// Team Internal Impl
/// ====================

///
/// Lower
///

template <typename AlgoType>
struct TeamVectorTrsvInternalLower {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const bool /*use_unit_diag*/, const int /*m*/,
                                           const ScalarType /*alpha*/, const ValueType *KOKKOS_RESTRICT /*A*/,
                                           const int /*as0*/, const int /*as1*/,
                                           /**/ ValueType *KOKKOS_RESTRICT /*b*/, const int /*bs0*/) {
    assert(false && "Error: encounter dummy impl");
    return 0;
  }
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT b, const int bs0) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, alpha, b, bs0);
    if (m <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      const int iend = m - p - 1;

      const ValueType *KOKKOS_RESTRICT a21 = iend ? A + (p + 1) * as0 + p * as1 : NULL;

      ValueType *KOKKOS_RESTRICT beta1 = b + p * bs0, *KOKKOS_RESTRICT b2 = iend ? beta1 + bs0 : NULL;

      member.team_barrier();
      ValueType local_beta1 = *beta1;
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        local_beta1             = local_beta1 / alpha11;

        member.team_barrier();
        Kokkos::single(Kokkos::PerTeam(member), [&]() { *beta1 = local_beta1; });
      }
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, iend),
                           [&](const int &i) { b2[i * bs0] -= a21[i * as0] * local_beta1; });
    }
  }
  return 0;
}

///
/// Upper
///

template <typename AlgoType>
struct TeamVectorTrsvInternalUpper {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const bool /*use_unit_diag*/, const int /*m*/,
                                           const ScalarType /*alpha*/, const ValueType *KOKKOS_RESTRICT /*A*/,
                                           const int /*as0*/, const int /*as1*/,
                                           /**/ ValueType *KOKKOS_RESTRICT /*b*/, const int /*bs0*/) {
    assert(false && "Error: encounter dummy impl");
    return 0;
  }
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT b, const int bs0) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, alpha, b, bs0);
    if (m <= 0) return 0;

    ValueType *KOKKOS_RESTRICT b0 = b;
    for (int p = (m - 1); p >= 0; --p) {
      const int iend = p;

      const ValueType *KOKKOS_RESTRICT a01  = A + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT beta1 = b + p * bs0;

      member.team_barrier();
      ValueType local_beta1 = *beta1;
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        local_beta1             = local_beta1 / alpha11;

        member.team_barrier();
        Kokkos::single(Kokkos::PerTeam(member), [&]() { *beta1 = local_beta1; });
      }
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, iend),
                           [&](const int &i) { b0[i * bs0] -= a01[i * as0] * local_beta1; });
    }
  }
  return 0;
}
}  // namespace KokkosBatched

#endif
