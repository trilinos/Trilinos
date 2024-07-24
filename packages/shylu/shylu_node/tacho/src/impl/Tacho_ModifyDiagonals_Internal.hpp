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
