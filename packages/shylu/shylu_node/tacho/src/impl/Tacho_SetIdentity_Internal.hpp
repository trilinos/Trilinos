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
#ifndef __TACHO_SET_IDENTITY_INTERNAL_HPP__
#define __TACHO_SET_IDENTITY_INTERNAL_HPP__

/// \file  Tacho_SetIdentity_Internal.hpp
/// \brief Set an identity matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct SetIdentity<Algo::Internal> {
  template <typename MemberType, typename ViewTypeA, typename ScalarType>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ScalarType &alpha) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m > 0 && n > 0) {
      const value_type diag(alpha), zero(0);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m),
                             [&](const ordinal_type &i) { A(i, j) = i == j ? diag : zero; });
      });
    }
    return 0;
  }
};

} // namespace Tacho
#endif
