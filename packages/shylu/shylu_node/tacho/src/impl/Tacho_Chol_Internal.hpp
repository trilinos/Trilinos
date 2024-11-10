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
#ifndef __TACHO_CHOL_INTERNAL_HPP__
#define __TACHO_CHOL_INTERNAL_HPP__

/// \file  Tacho_Chol_Internal.hpp
/// \brief LAPACK upper Cholesky factorization; hand-made team version for cuda
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_Team.hpp"

namespace Tacho {

/// LAPACK Chol
/// ===========
template <typename ArgUplo> struct Chol<ArgUplo, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

    int r_val = 0;
    const ordinal_type m = A.extent(0);
    if (m > 0)
      LapackTeam<value_type>::potrf(member, ArgUplo::param, m, A.data(), A.stride_1(), &r_val);
    return r_val;
  }
};
} // namespace Tacho

#endif
