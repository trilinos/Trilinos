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
#ifndef __TACHO_GEMM_SERIAL_HPP__
#define __TACHO_GEMM_SERIAL_HPP__

/// \file  Tacho_Gemm_Serial.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Serial.hpp"

namespace Tacho {

template <typename ArgTransA, typename ArgTransB> struct Gemm<ArgTransA, ArgTransB, Algo::Serial> {
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int invoke(const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B, const ScalarType beta,
                           const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type_a;
      typedef typename ViewTypeC::non_const_value_type value_type_c;

      const ordinal_type m = C.extent(0), n = C.extent(1),
                         k = (std::is_same<ArgTransB, Trans::NoTranspose>::value ? B.extent(0) : B.extent(1));

      BlasSerial<value_type_a>::gemm(ArgTransA::param, ArgTransB::param, m, n, k,
                                     value_type_a(alpha), A.data(), A.stride_1(),
                                                          B.data(), B.stride_1(),
                                     value_type_c(beta),  C.data(), C.stride_1());
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A,
                           const ViewTypeB &B, const ScalarType beta, const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      invoke(alpha, A, B, beta, C);
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
