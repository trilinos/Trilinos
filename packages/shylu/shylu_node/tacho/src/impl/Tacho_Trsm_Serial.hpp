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
#ifndef __TACHO_TRSM_SERIAL_HPP__
#define __TACHO_TRSM_SERIAL_HPP__

/// \file  Tacho_Trsm_Serial.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Serial.hpp"

namespace Tacho {

template <typename ArgSide, typename ArgUplo, typename ArgTransA>
struct Trsm<ArgSide, ArgUplo, ArgTransA, Algo::Serial> {

  template <typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int invoke(const DiagType diagA, const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      typedef typename ViewTypeB::non_const_value_type value_type_b;

      static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
      static_assert(ViewTypeB::rank == 2, "B is not rank 2 view.");

      static_assert(std::is_same<value_type, value_type_b>::value, "A and B do not have the same value type.");

      const ordinal_type m = B.extent(0);
      const ordinal_type n = B.extent(1);

      if (m > 0 && n > 0)
        //Blas<value_type>::trsm(ArgSide::param, ArgUplo::param, ArgTransA::param, diagA.param, m, n, value_type(alpha),
        //                       A.data(), A.stride_1(), B.data(), B.stride_1());
        BlasSerial<value_type>::trsm(ArgSide::param, ArgUplo::param, ArgTransA::param, diagA.param, m, n, value_type(alpha),
                                     A.data(), A.stride_1(), B.data(), B.stride_1());
    } else {
      TACHO_TEST_FOR_ABORT(true, "This function is only allowed in host space.");
    }
    return 0;
  }

  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline  static int invoke(MemberType &member, const DiagType diagA, const ScalarType alpha,
                            const ViewTypeA &A, const ViewTypeB &B) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      invoke(diagA, alpha, A, B);
    } else {
      TACHO_TEST_FOR_ABORT(true, "This function is only allowed in host space.");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
