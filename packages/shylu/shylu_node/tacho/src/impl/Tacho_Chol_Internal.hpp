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
