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
#ifndef __TACHO_MODIFY_DIAGONALS_ON_DEVICE_HPP__
#define __TACHO_MODIFY_DIAGONALS_ON_DEVICE_HPP__

/// \file  Tacho_ModifyDiagonals_OnDevice.hpp
/// \brief Modify diagonals
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct ModifyDiagonals<Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeA, typename ScalarType>
  inline static int invoke(MemberType &member, const ViewTypeA &A, const ScalarType &alpha) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

    const ordinal_type m = A.extent(0), n = A.extent(1), min_mn = m > n ? n : m;
    ;

    if (min_mn > 0) {
      using exec_space = MemberType;
      using range_policy_type = Kokkos::RangePolicy<exec_space>;

      const auto exec_instance = member;
      const auto policy = range_policy_type(exec_instance, 0, min_mn);
      Kokkos::parallel_for(
          policy, KOKKOS_LAMBDA(const ordinal_type &i) { A(i, i) += alpha; });
    }
    return 0;
  }
};

} // namespace Tacho
#endif
