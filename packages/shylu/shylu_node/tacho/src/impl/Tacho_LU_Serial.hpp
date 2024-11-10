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
#ifndef __TACHO_LU_SERIAL_HPP__
#define __TACHO_LU_SERIAL_HPP__

/// \file  Tacho_LU_Serial.hpp
/// \brief LAPACK LU factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_Serial.hpp"

namespace Tacho {

/// LAPACK LU
/// ==========
template <> struct LU<Algo::Serial> {
  template <typename ViewTypeA, typename ViewTypeP> inline static int invoke(const ViewTypeA &A, const ViewTypeP &P) {
    int r_val = 0;

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type;

      static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
      static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

      TACHO_TEST_FOR_EXCEPTION(P.extent(0) < 4 * A.extent(0), std::runtime_error, "P should be 4*A.extent(0) .");

      const ordinal_type m = A.extent(0), n = A.extent(1);
      if (m > 0 && n > 0) {
        /// factorize LU
        LapackSerial<value_type>::getrf(m, n, A.data(), A.stride_1(), P.data(), &r_val);
        TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, "LAPACK (getrf) returns non-zero error code.");
      }
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return r_val;
  }

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP>
  inline static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      int r_val = 0;
      r_val = invoke(A, P);
      return r_val;
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
  }

  template <typename ViewTypeP> inline static int modify(const ordinal_type m, const ViewTypeP &P) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeP::execution_space>;

    int r_val = 0;

    if constexpr(runOnHost) {
      static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");

      TACHO_TEST_FOR_EXCEPTION(int(P.extent(0)) < 4 * m, std::runtime_error, "P should be 4*m.");

      if (m > 0) {
        ordinal_type *KOKKOS_RESTRICT ipiv = P.data();
        ordinal_type *KOKKOS_RESTRICT fpiv = ipiv + m;
        ordinal_type *KOKKOS_RESTRICT perm = fpiv + m;
        ordinal_type *KOKKOS_RESTRICT peri = perm + m;

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
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return r_val;
  }

  template <typename MemberType, typename ViewTypeP>
  inline static int modify(MemberType &member, ordinal_type m, const ViewTypeP &P) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeP::execution_space>;

    if constexpr(runOnHost) {
      int r_val = 0;
      r_val = modify(m, P);
      return r_val;
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
  }
};

} // namespace Tacho

#endif
