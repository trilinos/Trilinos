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
#ifndef __TACHO_HERK_EXTERNAL_HPP__
#define __TACHO_HERK_EXTERNAL_HPP__

/// \file  Tacho_Herk_External.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTrans> struct Herk<ArgUplo, ArgTrans, Algo::External> {
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  inline static int invoke(const ScalarType alpha, const ViewTypeA &A, const ScalarType beta, const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      typedef typename ViewTypeC::non_const_value_type value_type_c;

      static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
      static_assert(ViewTypeC::rank == 2, "B is not rank 2 view.");

      static_assert(std::is_same<value_type, value_type_c>::value, "A and C do not have the same value type.");

      const ordinal_type n = C.extent(0),
                         k = (std::is_same<ArgTrans, Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));
      if (n > 0 && k > 0) {
        Blas<value_type>::herk(ArgUplo::param, ArgTrans::param, n, k, value_type(alpha), A.data(), A.stride_1(),
                               value_type(beta), C.data(), C.stride_1());
      }
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A,
                                           const ScalarType beta, const ViewTypeC &C) {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename ViewTypeA::execution_space>;

    if constexpr(runOnHost) {
      // Kokkos::single(Kokkos::PerTeam(member), [&]() {
      invoke(alpha, A, beta, C);
      //});
    } else {
      TACHO_TEST_FOR_ABORT(true, ">> This function is only allowed in host space.");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
