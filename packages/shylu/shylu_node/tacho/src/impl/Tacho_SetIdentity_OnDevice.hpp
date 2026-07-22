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
#ifndef __TACHO_SET_IDENTITY_ON_DEVICE_HPP__
#define __TACHO_SET_IDENTITY_ON_DEVICE_HPP__

/// \file  Tacho_SetIdentity_OnDevice.hpp
/// \brief Set an identity matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct SetIdentity<Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeA, typename ScalarType>
  inline static int invoke(MemberType &exec_instance, const ViewTypeA &A, const ScalarType &alpha) {

    typedef typename ViewTypeA::non_const_value_type value_type;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");

    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m > 0 && n > 0) {
      const value_type diag(alpha), zero(0);
      using exec_space = MemberType;
      using team_policy_type = Kokkos::TeamPolicy<exec_space>;

      const auto policy = team_policy_type(exec_instance, n, Kokkos::AUTO);
      Kokkos::parallel_for(
          policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {
            const ordinal_type j = member.league_rank();
            Kokkos::parallel_for(
                Kokkos::TeamVectorRange(member, m),
                [&, diag, zero, A,
                 j](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
                  A(i, j) = i == j ? diag : zero;
                });
          });
    }
    return 0;
  }
};

} // namespace Tacho
#endif
