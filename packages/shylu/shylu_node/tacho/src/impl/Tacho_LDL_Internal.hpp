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
#ifndef __TACHO_LDL_INTERNAL_HPP__
#define __TACHO_LDL_INTERNAL_HPP__

/// \file  Tacho_LDL_Internal.hpp
/// \brief LDL team factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_Team.hpp"

namespace Tacho {

/// LAPACK LDL
/// ==========
template <> struct LDL<Uplo::Lower, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P,
                                           const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    // typedef typename ViewTypeP::non_const_value_type p_value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");
    static_assert(ViewTypeW::rank == 1, "W is not rank 1 view.");

    TACHO_TEST_FOR_ABORT(P.extent(0) < 4 * A.extent(0), "P should be 4*A.extent(0) .");

    int r_val(0);
    const ordinal_type m = A.extent(0);
    if (m > 0) {
      /// factorize LDL
      LapackTeam<value_type>::sytrf(member, Uplo::Lower::param, m, A.data(), A.stride_1(), P.data(), W.data(), &r_val);
    }
    return r_val;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  KOKKOS_INLINE_FUNCTION static int modify(MemberType &member, const ViewTypeA &A, const ViewTypeP &P,
                                           const ViewTypeD &D) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");
    static_assert(ViewTypeD::rank == 2, "D is not rank 2 view.");

    TACHO_TEST_FOR_ABORT(D.extent(0) < A.extent(0), "D extent(0) is smaller than A extent(0).");
    TACHO_TEST_FOR_ABORT(D.extent(1) != 2, "D is supposed to store 2x2 blocks .");
    TACHO_TEST_FOR_ABORT(P.extent(0) < 4 * A.extent(0), "P should be 4*A.extent(0) .");

    int r_val = 0;
    const ordinal_type m = A.extent(0);
    if (m > 0) {
      value_type *KOKKOS_RESTRICT Aptr = A.data();
      ordinal_type *KOKKOS_RESTRICT ipiv = P.data(), *KOKKOS_RESTRICT fpiv = ipiv + m, *KOKKOS_RESTRICT perm = fpiv + m,
                                 *KOKKOS_RESTRICT peri = perm + m;
      const value_type one(1);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { perm[i] = i; });
      member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &j) {
        const bool single = (j == (m - 1));
        for (ordinal_type i = 0 /*,cnt=0*/; i < m; ++i) {
          // if (ipiv[i] <= 0) {
          //   if (++cnt%2) {
          //     if (single) {
          //       ipiv[i] = 0; /// invalidate this pivot
          //       fpiv[i] = 0;

          //       D(i,0) = A(i,  i);
          //       D(i,1) = A(i+1,i); /// symmetric
          //       A(i,i) = one;
          //     }
          //   } else {
          //     const ordinal_type fla_pivot = -ipiv[i]-i-1;
          //     if (single) {
          //       fpiv[i] = fla_pivot;
          //     }
          //     if (fla_pivot) {
          //       value_type *KOKKOS_RESTRICT src = Aptr + i;
          //       value_type *KOKKOS_RESTRICT tgt = src + fla_pivot;
          //       if (j<(i-1)) {
          //         const ordinal_type idx = j*m;
          //         swap(src[idx], tgt[idx]);
          //       }
          //     }

          //     if (single) {
          //       D(i,0) = A(i,i-1);
          //       D(i,1) = A(i,i  );
          //       A(i,i-1) = zero; A(i,i) = one;
          //     }
          //   }
          // } else
          {
            const ordinal_type fla_pivot = ipiv[i] - i - 1;
            if (single) {
              fpiv[i] = fla_pivot;
            }
            if (fla_pivot) {
              value_type *src = Aptr + i;
              value_type *tgt = src + fla_pivot;
              if (j < i) {
                const ordinal_type idx = j * m;
                swap(src[idx], tgt[idx]);
              }
            }

            if (single) {
              D(i, 0) = A(i, i);
              A(i, i) = one;
            }
          }

          /// apply pivots to perm vector
          if (single) {
            if (fpiv[i]) {
              const ordinal_type pidx = i + fpiv[i];
              swap(perm[i], perm[pidx]);
            }
          }
        }
      });
      member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { peri[perm[i]] = i; });
    }

    /// no piv version
    // if (m > 0) {
    //   ordinal_type
    //     *KOKKOS_RESTRICT ipiv = P.data(),
    //     *KOKKOS_RESTRICT fpiv = ipiv + m,
    //     *KOKKOS_RESTRICT perm = fpiv + m,
    //     *KOKKOS_RESTRICT peri = perm + m;
    //   const value_type one(1);
    //   Kokkos::parallel_for(Kokkos::TeamVectorRange(member,m),[&](const int &i) {
    //       D(i,0) = A(i,i);
    //       A(i,i) = one;
    //       ipiv[i] = i+1;
    //       fpiv[i] = 0;
    //       perm[i] = i;
    //       peri[i] = i;
    //     });
    // }
    return r_val;
  }
};

} // namespace Tacho

#endif
