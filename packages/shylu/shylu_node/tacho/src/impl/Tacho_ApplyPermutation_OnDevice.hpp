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
#ifndef __TACHO_APPLY_PERMUTATION_ON_DEVICE_HPP__
#define __TACHO_APPLY_PERMUTATION_ON_DEVICE_HPP__

/// \file  Tacho_ApplyPivots_OnDevice.hpp
/// \brief Apply pivots on device
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeB>
  inline static int invoke(const MemberType &member, const ViewTypeA &A, const ViewTypeP &P, const ViewTypeB &B) {
    // typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = A.extent(0), n = A.extent(1), plen = P.extent(0);

    if (m == plen) {
      if (A.span() > 0) {
        using exec_space = MemberType;
        const auto &exec_instance = member;

        const Kokkos::RangePolicy<exec_space> policy(exec_instance, 0, m * n);
        if (A.span() > 0) {
          if (n == 1) {
            Kokkos::parallel_for(
                policy, KOKKOS_LAMBDA(const ordinal_type &i) {
                  const ordinal_type idx = P(i);
                  B(i, 0) = A(idx, 0);
                });
          } else {
            Kokkos::parallel_for(
                policy, KOKKOS_LAMBDA(const ordinal_type &ij) {
                  const ordinal_type i = ij % m, j = ij / m;
                  const ordinal_type idx = P(i);
                  B(i, j) = A(idx, j);
                });
          }
        }
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
