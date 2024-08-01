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
#ifndef __TACHO_GEMM_TRIANGULAR_INTERNAL_HPP__
#define __TACHO_GEMM_TRIANGULAR_INTERNAL_HPP__

/// \file  Tacho_Gemm_Internal.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Team.hpp"

namespace Tacho {

template <> struct GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::Internal> {
  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A,
                                           const ViewTypeB &B, const ScalarType beta, const ViewTypeC &C) {

    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeB::non_const_value_type value_type_b;
    typedef typename ViewTypeC::non_const_value_type value_type_c;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeB::rank == 2, "B is not rank 2 view.");
    static_assert(ViewTypeC::rank == 2, "C is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_b>::value && std::is_same<value_type_b, value_type_c>::value,
                  "A, B and C do not have the same value type.");

    const ordinal_type m = C.extent(0);
    const ordinal_type n = C.extent(1);
    const ordinal_type k = B.extent(0);

    if (m > 0 && n > 0 && k > 0) {
      if (m == n) {
        BlasTeam<value_type>::gemm_upper(member, Trans::Transpose::param, Trans::NoTranspose::param, m, n, k,
                                         value_type(alpha), A.data(), A.stride_1(), B.data(), B.stride_1(),
                                         value_type(beta), C.data(), C.stride_1());
      } else {
        TACHO_TEST_FOR_ABORT(true, "C is not a square matrix");
      }
    }
    return 0;
  }
};
} // namespace Tacho
#endif
