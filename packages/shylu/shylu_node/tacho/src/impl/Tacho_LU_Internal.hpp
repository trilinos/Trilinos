// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_LU_INTERNAL_HPP__
#define __TACHO_LU_INTERNAL_HPP__

/// \file  Tacho_LU_Internal.hpp
/// \brief LU team factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_Team.hpp"

namespace Tacho {

/// LAPACK LU
/// ==========
template <> struct LU<Algo::Internal> {
  template <typename MemberType, typename ViewTypeA, typename ViewTypeP>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    // typedef typename ViewTypeP::non_const_value_type p_value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

    TACHO_TEST_FOR_ABORT(P.extent(0) < 4 * A.extent(0), "P should be 4*A.extent(0) .");

    int r_val(0);
    const ordinal_type m = A.extent(0), n = A.extent(1);
    if (m > 0 && n > 0) {
      /// factorize LU
      LapackTeam<value_type>::getrf(member, m, n, A.data(), A.stride_1(), P.data(), &r_val);
    }
    return r_val;
  }

  template <typename MemberType, typename ViewTypeP>
  KOKKOS_INLINE_FUNCTION static int modify(const MemberType &member, const ordinal_type m, const ViewTypeP &P) {
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

    TACHO_TEST_FOR_ABORT(P.extent(0) < 4 * m, "P should be 4*m .");

    int r_val = 0;
    if (m > 0) {
      ordinal_type *KOKKOS_RESTRICT ipiv = P.data(), *KOKKOS_RESTRICT fpiv = ipiv + m, *KOKKOS_RESTRICT perm = fpiv + m,
                                 *KOKKOS_RESTRICT peri = perm + m;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) {
        perm[i] = i;
        fpiv[i] = ipiv[i] - i - 1;
      });
      member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        for (ordinal_type i = 0; i < m; ++i) {
          /// apply pivots to perm vector
          if (fpiv[i]) {
            const ordinal_type pidx = i + fpiv[i];
            swap(perm[i], perm[pidx]);
          }
        }
      });
      member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { peri[perm[i]] = i; });
    }

    return r_val;
  }
};

} // namespace Tacho

#endif
