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
#ifndef __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_ON_DEVICE_HPP
#define __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_ON_DEVICE_HPP

/// \file  Tacho_Scale2x2_BlockInverseDiagonals_OnDevice.hpp
/// \brief inverse scale
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

/// row exchange
template <> struct Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeP, typename ViewTypeD, typename ViewTypeA>
  inline static int invoke(MemberType &exec_instance, const ViewTypeP &P, const ViewTypeD &D, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = A.extent(0), n = A.extent(1);
    if (A.extent(0) == D.extent(0)) {
      if (A.span() > 0) {
        using exec_space = MemberType;

        if (n == 1) {
          Kokkos::RangePolicy<exec_space> policy(exec_instance, 0, m);
          Kokkos::parallel_for(
              policy, KOKKOS_LAMBDA(const ordinal_type &i) {
                const ordinal_type pval = P(i);
                if (pval == 0) {
                  /// do nothing
                } else if (pval < 0) {
                  /// take 2x2 block to D
                  const value_type a00 = D(i - 1, 0), a01 = D(i - 1, 1), a10 = D(i, 0), a11 = D(i, 1);
                  const value_type det = a00 * a11 - a10 * a01;
                  const value_type x0 = A(i - 1, 0), x1 = A(i, 0);

                  A(i - 1, 0) = (a11 * x0 - a10 * x1) / det;
                  A(i, 0) = (-a10 * x0 + a00 * x1) / det;
                } else {
                  const value_type a00 = D(i, 0);
                  A(i, 0) /= a00;
                }
              });
        } else {
          using policy_type = Kokkos::TeamPolicy<exec_space>;
          policy_type policy(exec_instance, m, Kokkos::AUTO);
          Kokkos::parallel_for(
              policy, KOKKOS_LAMBDA(const typename policy_type::member_type &member) {
                const ordinal_type i = member.league_rank();
                const ordinal_type pval = P(i);
                if (pval == 0) {
                  /// do nothing
                } else if (pval < 0) {
                  /// take 2x2 block to D
                  const value_type a00 = D(i - 1, 0), a01 = D(i - 1, 1), a10 = D(i, 0), a11 = D(i, 1);
                  const value_type det = a00 * a11 - a10 * a01;
                  Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [=](const ordinal_type &j) {
                    const value_type x0 = A(i - 1, j), x1 = A(i, j);
                    A(i - 1, j) = (a11 * x0 - a10 * x1) / det;
                    A(i, j) = (-a10 * x0 + a00 * x1) / det;
                  });
                } else {
                  const value_type a00 = D(i, 0);
                  Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n),
                                       [=](const ordinal_type &j) { A(i, j) /= a00; });
                }
              });
        }
      }
    } else {
      printf("Error: Scale2x2_BlockInverseDiagonals<Side::Left,Algo::Internal> A is not square\n");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
