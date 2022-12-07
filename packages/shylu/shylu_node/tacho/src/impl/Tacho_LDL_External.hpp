// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_LDL_EXTERNAL_HPP__
#define __TACHO_LDL_EXTERNAL_HPP__

/// \file  Tacho_LDL_External.hpp
/// \brief LAPACK LDL factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_External.hpp"

namespace Tacho {

/// LAPACK LDL
/// ==========
template <> struct LDL<Uplo::Lower, Algo::External> {
  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  inline static int invoke(const ViewTypeA &A, const ViewTypeP &P, const ViewTypeW &W) {
    int r_val = 0;
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");
    static_assert(ViewTypeW::rank == 1, "W is not rank 1 view.");

    TACHO_TEST_FOR_EXCEPTION(P.extent(0) < 4 * A.extent(0), std::runtime_error, "P should be 4*A.extent(0) .");

    const ordinal_type m = A.extent(0);
    if (m > 0) {
      /// factorize LDL
      Lapack<value_type>::sytrf('L', m, A.data(), A.stride_1(), P.data(), W.data(), W.extent(0), &r_val);
      TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, "LAPACK (sytrf) returns non-zero error code.");
    }
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
    return r_val;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P,
                                           const ViewTypeW &W) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    int r_val = 0;
    r_val = invoke(A, P, W);
    return r_val;
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
  }

  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  inline static int modify(const ViewTypeA &A, const ViewTypeP &P, const ViewTypeD &D) {
    int r_val = 0;
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");
    static_assert(ViewTypeD::rank == 2, "D is not rank 2 view.");

    TACHO_TEST_FOR_EXCEPTION(D.extent(0) < A.extent(0), std::runtime_error, "D extent(0) is smaller than A extent(0).");
    TACHO_TEST_FOR_EXCEPTION(D.extent(1) != 2, std::runtime_error, "D is supposed to store 2x2 blocks .");
    TACHO_TEST_FOR_EXCEPTION(P.extent(0) < 4 * A.extent(0), std::runtime_error, "P should be 4*A.extent(0) .");

    const ordinal_type m = A.extent(0);
    if (m > 0) {
      value_type *__restrict__ Aptr = A.data();
      ordinal_type *__restrict__ ipiv = P.data(), *__restrict__ fpiv = ipiv + m, *__restrict__ perm = fpiv + m,
                                 *__restrict__ peri = perm + m;

      const value_type one(1), zero(0);
      for (ordinal_type i = 0; i < m; ++i)
        perm[i] = i;
      for (ordinal_type i = 0; i < m; ++i) {
        if (ipiv[i] < 0) {
          {
            // first pivot
            ipiv[i] = 0; /// invalidate this pivot
            fpiv[i] = 0;

            D(i, 0) = A(i, i);
            D(i, 1) = A(i + 1, i); /// symmetric
            A(i, i) = one;
          }
          {
            // second pivot
            i++;
            const ordinal_type fla_pivot = -ipiv[i] - i - 1;
            fpiv[i] = fla_pivot;
            if (fla_pivot) {
              value_type *__restrict__ src = Aptr + i;
              value_type *__restrict__ tgt = src + fla_pivot;
              for (ordinal_type j = 0; j < (i - 1); ++j) {
                const ordinal_type idx = j * m;
                swap(src[idx], tgt[idx]);
              }
            }

            D(i, 0) = A(i, i - 1);
            D(i, 1) = A(i, i);
            A(i, i - 1) = zero;
            A(i, i) = one;
          }
        } else {
          const ordinal_type fla_pivot = ipiv[i] - i - 1;
          fpiv[i] = fla_pivot;
          if (fla_pivot) {
            value_type *src = Aptr + i;
            value_type *tgt = src + fla_pivot;
            for (ordinal_type j = 0; j < i; ++j) {
              const ordinal_type idx = j * m;
              swap(src[idx], tgt[idx]);
            }
          }

          D(i, 0) = A(i, i);
          A(i, i) = one;
        }

        /// apply pivots to perm vector
        if (fpiv[i]) {
          const ordinal_type pidx = i + fpiv[i];
          swap(perm[i], perm[pidx]);
        }
      }
      for (ordinal_type i = 0; i < m; ++i)
        peri[perm[i]] = i;
    }
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
    return r_val;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  KOKKOS_INLINE_FUNCTION static int modify(MemberType &member, const ViewTypeA &A, const ViewTypeP &P,
                                           const ViewTypeD &D) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    int r_val = 0;
    r_val = modify(A, P, D);
    return r_val;
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
  }
};

} // namespace Tacho

#endif
