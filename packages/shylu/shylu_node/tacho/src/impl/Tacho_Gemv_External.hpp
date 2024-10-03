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
#ifndef __TACHO_GEMV_EXTERNAL_HPP__
#define __TACHO_GEMV_EXTERNAL_HPP__

/// \file  Tacho_Gemv_External.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgTrans> struct Gemv<ArgTrans, Algo::External> {
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int invoke(const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B, const ScalarType beta,
                           const ViewTypeC &C) {

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

      const ordinal_type m = C.extent(0), n = C.extent(1);

      if (m > 0 && n > 0) {
        if (n == 1) {
          const int mm = A.extent(0), nn = A.extent(1);
          Blas<value_type>::gemv(ArgTrans::param, mm, nn, value_type(alpha), A.data(), A.stride_1(), B.data(),
                                 B.stride_0(), value_type(beta), C.data(), C.stride_0());
        } else {
          const int mm = C.extent(0), nn = C.extent(1), kk = B.extent(0);
          Blas<value_type>::gemm(ArgTrans::param, Trans::NoTranspose::param, mm, nn, kk, value_type(alpha), A.data(),
                                 A.stride_1(), B.data(), B.stride_1(), value_type(beta), C.data(), C.stride_1());
        }
      }
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A,
                                           const ViewTypeB &B, const ScalarType beta, const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      // Kokkos::single(Kokkos::PerTeam(member), [&]() {
      invoke(alpha, A, B, beta, C);
      //});
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
