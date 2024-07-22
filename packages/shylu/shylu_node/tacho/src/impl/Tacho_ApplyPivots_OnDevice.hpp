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
#ifndef __TACHO_APPLY_PIVOTS_ON_DEVICE_HPP__
#define __TACHO_APPLY_PIVOTS_ON_DEVICE_HPP__

/// \file  Tacho_ApplyPivots_OnDevice.hpp
/// \brief Apply pivots on device
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct ApplyPivots<PivotMode::Lapack, Side::Left, Direct::Forward, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeP, typename ViewTypeA>
  inline static int invoke(MemberType &member, const ViewTypeP &P, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == P.extent(0)) {
      if (A.span() > 0) {
        using exec_space = MemberType;
        using policy_type = Kokkos::RangePolicy<exec_space>;

        const auto &exec_instance = member;
        const auto policy = policy_type(exec_instance, 0, n);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &j) {
              for (ordinal_type i = 0; i < m; ++i) {
                const ordinal_type piv = P(i);
                if (piv == i || piv == 0) {
                  /// no pivot
                } else if (piv > 0) {
                  /// 1x1 pivot
                  const ordinal_type p = piv - 1;
                  const value_type tmp = A(i, j);
                  A(i, j) = A(p, j);
                  A(p, j) = tmp;
                } else {
                  /// 2x2 pivot
                  const int p = -piv - 1;
                  const value_type tmp_a = A(i, j), tmp_b = A(i + 1, j);
                  A(i, j) = A(p, j);
                  A(i + 1, j) = A(p + 1, j);
                  A(p, j) = tmp_a;
                  A(p + 1, j) = tmp_b;
                }
              }
            });
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }
};

template <> struct ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeP, typename ViewTypeA>
  inline static int invoke(MemberType &member, const ViewTypeP &P, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = A.extent(0), n = A.extent(1), plen = P.extent(0);

    if (m == plen) {
      if (A.span() > 0) {
        using exec_space = MemberType;
        using policy_type = Kokkos::RangePolicy<exec_space>;

        const auto exec_instance = member;
        const auto policy = policy_type(exec_instance, 0, n);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &j) {
              for (ordinal_type i = 0; i < m; ++i) {
                const ordinal_type piv = P(i);
                if (piv == 0) {
                  /// no pivot
                } else {
                  /// 1x1 pivot
                  const ordinal_type p = i + piv;
                  const value_type tmp = A(i, j);
                  A(i, j) = A(p, j);
                  A(p, j) = tmp;
                }
              }
            });
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }
};

template <> struct ApplyPivots<PivotMode::Flame, Side::Left, Direct::Backward, Algo::OnDevice> {
  template <typename MemberType, typename ViewTypeP, typename ViewTypeA>
  inline static int invoke(MemberType &member, const ViewTypeP &P, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = A.extent(0), n = A.extent(1), plen = P.extent(0);

    if (m == plen) {
      if (A.span() > 0) {
        using exec_space = MemberType;
        using policy_type = Kokkos::RangePolicy<exec_space>;

        const auto exec_instance = member;
        const auto policy = policy_type(exec_instance, 0, n);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &j) {
              for (ordinal_type i = (m - 1); i >= 0; --i) {
                const ordinal_type piv = P(i);
                if (piv == 0) {
                  /// no pivot
                } else {
                  /// 1x1 pivot
                  const ordinal_type p = i + piv;
                  const value_type tmp = A(i, j);
                  A(i, j) = A(p, j);
                  A(p, j) = tmp;
                }
              }
            });
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A and P dimension does not match");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
