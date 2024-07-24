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
#ifndef __TACHO_CHOL_EXTERNAL_HPP__
#define __TACHO_CHOL_EXTERNAL_HPP__

/// \file  Tacho_Chol_Upper_ExternalLapack.hpp
/// \brief LAPACK upper Cholesky factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_External.hpp"

namespace Tacho {

/// LAPACK Chol
/// ===========
template <typename ArgUplo> struct Chol<ArgUplo, Algo::External> {
  template <typename ViewTypeA> inline static int invoke(const ViewTypeA &A) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

      int r_val = 0;
      const ordinal_type m = A.extent(0);
      if (m > 0) {
        Lapack<value_type>::potrf(ArgUplo::param, m, A.data(), A.stride_1(), &r_val);
        TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, "LAPACK (potrf) returns non-zero error code.");
      }
      return r_val;
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
      return -1;
    }
  }

  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      int r_val = 0;
      // Kokkos::single(Kokkos::PerTeam(member), [&]() {
      r_val = invoke(A);
      KOKKOS_IF_ON_HOST((TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, "LAPACK (potrf) returns non-zero error code.");))
      //});
      return r_val;
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
      return 0;
    }
  }
};
} // namespace Tacho

#endif
