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
#ifndef __TACHO_NONPIV_LDL_ON_DEVICE_HPP__
#define __TACHO_NONPIV_LDL_ON_DEVICE_HPP__

/// \file  Tacho_NonPivLDL_OnDevice.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <typename ArgUplo> struct LDL_nopiv<ArgUplo, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeA, typename ViewTypeW>
  inline static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeW &W) {

    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    for (ordinal_type i = 0; i < m; i++) {
      ////for (ordinal_type j = i+1; j < m; j++) A(i, j) /= A(i, i);
      const auto &exec_instance = member;
      const auto policy_scale = policy_type(exec_instance, i+1, m);
      Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
          A(i, j) /= A(i, i); });

      ////for (ordinal_type k = i+1; k < m; k++) {
      ////  for (ordinal_type j = k; j < m; j++) A(k, j) -= A(i, k) * A(i, i) * A(i, j);
      ////}
      policy_type policy_update(exec_instance, i+1, m);
      Kokkos::parallel_for(policy_update, KOKKOS_LAMBDA(const ordinal_type &k) {
        for (ordinal_type j = k; j < m; j++) A(k, j) -= A(i, k) * A(i, i) * A(i, j);
      });
    }

    return r_val;
  }
};

} // namespace Tacho
#endif
