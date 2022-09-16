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
#ifndef __TACHO_LU_EXTERNAL_HPP__
#define __TACHO_LU_EXTERNAL_HPP__

/// \file  Tacho_LU_External.hpp
/// \brief LAPACK LU factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_External.hpp"

namespace Tacho {

/// LAPACK LU
/// ==========
template <> struct LU<Algo::External> {
  template <typename ViewTypeA, typename ViewTypeP> inline static int invoke(const ViewTypeA &A, const ViewTypeP &P) {
    int r_val = 0;
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

    TACHO_TEST_FOR_EXCEPTION(P.extent(0) < 4 * A.extent(0), std::runtime_error, "P should be 4*A.extent(0) .");

    const ordinal_type m = A.extent(0), n = A.extent(1);
    if (m > 0 && n > 0) {
      /// factorize LU
      Lapack<value_type>::getrf(m, n, A.data(), A.stride_1(), P.data(), &r_val);
      TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, "LAPACK (getrf) returns non-zero error code.");
    }
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
    return r_val;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    int r_val = 0;
    r_val = invoke(A, P);
    return r_val;
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
  }

  template <typename ViewTypeP> inline static int modify(const ordinal_type m, const ViewTypeP &P) {
    int r_val = 0;
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

    TACHO_TEST_FOR_EXCEPTION(int(P.extent(0)) < 4 * m, std::runtime_error, "P should be 4*m.");

    if (m > 0) {
      ordinal_type *__restrict__ ipiv = P.data();
      ordinal_type *__restrict__ fpiv = ipiv + m;
      ordinal_type *__restrict__ perm = fpiv + m;
      ordinal_type *__restrict__ peri = perm + m;

      for (ordinal_type i = 0; i < m; ++i)
        perm[i] = i;
      for (ordinal_type i = 0; i < m; ++i) {
        const ordinal_type fla_pivot = ipiv[i] - i - 1;
        fpiv[i] = fla_pivot;

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

  template <typename MemberType, typename ViewTypeP>
  KOKKOS_INLINE_FUNCTION static int modify(MemberType &member, ordinal_type m, const ViewTypeP &P) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    int r_val = 0;
    r_val = modify(m, P);
    return r_val;
#else
    TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
#endif
  }
};

} // namespace Tacho

#endif
