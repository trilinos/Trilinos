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
#ifndef __TACHO_MODIFY_DIAGONALS_INTERNAL_HPP__
#define __TACHO_MODIFY_DIAGONALS_INTERNAL_HPP__

/// \file  Tacho_ModifyDiagonals_Internal.hpp
/// \brief Modify diagonals
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct ModifyDiagonals<Algo::Internal> {
  template <typename MemberType, typename ViewTypeA, typename ScalarType>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ScalarType &alpha) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

    const ordinal_type m = A.extent(0), n = A.extent(1), min_mn = m > n ? n : m;

    if (min_mn > 0) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, min_mn), [&](const ordinal_type &i) { A(i, i) += alpha; });
    }
    return 0;
  }
};

} // namespace Tacho
#endif
