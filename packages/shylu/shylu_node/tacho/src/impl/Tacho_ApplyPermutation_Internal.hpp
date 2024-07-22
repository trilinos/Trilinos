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
#ifndef __TACHO_APPLY_PERMUTATION_INTERNAL_HPP__
#define __TACHO_APPLY_PERMUTATION_INTERNAL_HPP__

/// \file  Tacho_ApplyPermutation_Internal.hpp
/// \brief Apply pivots
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

/// row exchange
template <> struct ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal> {
  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeB>
  KOKKOS_INLINE_FUNCTION static int invoke(const ViewTypeA &A, const ViewTypeP &P, const ViewTypeB &B) {
    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        if (n == 1) { /// vector
          for (ordinal_type i = 0; i < m; ++i) {
            const ordinal_type idx = P(i);
            B(i, 0) = A(idx, 0);
          }
        } else { /// matrix
          for (ordinal_type i = 0; i < m; ++i) {
            const ordinal_type idx = P(i);
            for (ordinal_type j = 0; j < n; ++j) {
              B(i, j) = A(idx, j);
            }
          }
        }
      }
    } else {
      Kokkos::printf("Error: ApplyPermutation<Algo::Internal> A extent(0) does not match to P extent(0)\n");
    }
    return 0;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeB>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ViewTypeA &A, const ViewTypeP &P,
                                           const ViewTypeB &B) {
    KOKKOS_IF_ON_DEVICE((
    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        if (n == 1) { /// vector
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const ordinal_type &i) {
            const ordinal_type idx = P(i);
            B(i, 0) = A(idx, 0);
          });
        } else { /// matrix
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m * n), [&](const ordinal_type &ij) {
            const ordinal_type i = ij % m, j = ij / m;
            const ordinal_type idx = P(i);
            B(i, j) = A(idx, j);
          });
        }
      }
    } else {
      Kokkos::printf("Error: ApplyPermutation<Algo::Internal> A extent(0) does not match to P extent(0)\n");
    }))
    KOKKOS_IF_ON_HOST((invoke(A, P, B);))
    return 0;
  }
};

} // namespace Tacho
#endif
