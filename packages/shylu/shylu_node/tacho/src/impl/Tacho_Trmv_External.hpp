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
#ifndef __TACHO_TRMV_EXTERNAL_HPP__
#define __TACHO_TRMV_EXTERNAL_HPP__

/// \file  Tacho_Trmv_External.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTrans> struct Trmv<ArgUplo, ArgTrans, Algo::External> {
  template <typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(const DiagType diag,
                                           const ScalarType alpha, const ViewTypeA &A,
                                                                   const ViewTypeB &B,
                                           const ScalarType beta,  const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      typedef typename ViewTypeB::non_const_value_type value_type_b;
      typedef typename ViewTypeC::non_const_value_type value_type_c;

      static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
      static_assert(ViewTypeB::rank == 2, "B is not rank 2 view.");
      static_assert(ViewTypeC::rank == 2, "C is not rank 2 view.");
      static_assert(std::is_same<value_type, value_type_b>::value && std::is_same<value_type_b, value_type_c>::value,
                    "A, B and C do not have the same value type.");

      const int mC = C.extent(0), nC = C.extent(1);
      if (mC > 0 && nC > 0) {
        const int mA = A.extent(0), nA = A.extent(1);
        if (nC == 1) {
          BlasSerial<value_type>::trmv(ArgUplo::param, ArgTrans::param, diag.param, mA, nA,
                                       value_type(alpha), A.data(), A.stride(1),
                                                          B.data(), B.stride(0),
                                       value_type(beta),  C.data(), C.stride(0));
        } else {
          int k = (ArgTrans::param == 'N' || ArgTrans::param == 'n' ? nA : mA);
          BlasSerial<value_type>::trmm(ArgUplo::param, ArgTrans::param, diag.param, mC, nC, k,
                                       value_type(alpha), A.data(), A.stride(1),
                                                          B.data(), B.stride(1),
                                       value_type(beta),  C.data(), C.stride(1));
        }
      }
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }

  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const DiagType diag,
                                           const ScalarType alpha, const ViewTypeA &A,
                                                                   const ViewTypeB &B,
                                           const ScalarType beta,  const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      // Kokkos::single(Kokkos::PerTeam(member), [&]() {
      invoke(diag, alpha, A, B, beta, C);
      //});
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
