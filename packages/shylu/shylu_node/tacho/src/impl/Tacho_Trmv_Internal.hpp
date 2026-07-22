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
#ifndef __TACHO_TRMV_INTERNAL_HPP__
#define __TACHO_TRMV_INTERNAL_HPP__

/// \file  Tacho_Trmv_Internal.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Team.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTrans> struct Trmv<ArgUplo, ArgTrans, Algo::Internal> {
  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const DiagType diag,
                                           const ScalarType alpha, const ViewTypeA &A,
                                                                   const ViewTypeB &B,
                                           const ScalarType beta,  const ViewTypeC &C) {

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
      const int mA = A.extent(0), nA = A.extent(1);
      if (n == 1) {
        BlasTeam<value_type>::trmv(member, ArgUplo::param, ArgTrans::param, diag.param,
                                   mA, nA, value_type(alpha), A.data(), A.stride(1),
                                                              B.data(), B.stride(0),
                                           value_type(beta),  C.data(), C.stride(0));
      } else {
        // TODO: need trmm team
        for (ordinal_type j = 0; j < n; j++ ) {
          BlasTeam<value_type>::trmv(member, ArgUplo::param, ArgTrans::param, diag.param,
                                     mA, nA, value_type(alpha), A.data(), A.stride(1),
                                                                &B(0, j), B.stride(0),
                                             value_type(beta),  &C(0, j), C.stride(0));
        }
      }
    }
    return 0;
  }
};

} // namespace Tacho
#endif
