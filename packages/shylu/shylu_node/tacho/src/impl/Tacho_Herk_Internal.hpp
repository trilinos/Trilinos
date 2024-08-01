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
#ifndef __TACHO_HERK_INTERNAL_HPP__
#define __TACHO_HERK_INTERNAL_HPP__

/// \file  Tacho_Herk_Internal.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTrans> struct Herk<ArgUplo, ArgTrans, Algo::Internal> {
  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A,
                                           const ScalarType beta, const ViewTypeC &C) {

    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeC::non_const_value_type value_type_c;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeC::rank == 2, "B is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_c>::value, "A and C do not have the same value type.");

    const ordinal_type n = C.extent(0);
    const ordinal_type k = (std::is_same<ArgTrans, Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));

    if (n > 0 && k > 0)
      BlasTeam<value_type>::herk(member, ArgUplo::param, ArgTrans::param, n, k, value_type(alpha), A.data(),
                                 A.stride_1(), value_type(beta), C.data(), C.stride_1());
    return 0;
  }
};

} // namespace Tacho
#endif
