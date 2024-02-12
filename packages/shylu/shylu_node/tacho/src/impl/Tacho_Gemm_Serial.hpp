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
