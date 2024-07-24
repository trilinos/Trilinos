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
#ifndef __TACHO_SYMMETRIZE_ON_DEVICE_HPP__
#define __TACHO_SYMMETRIZE_ON_DEVICE_HPP__

/// \file  Tacho_Symmetrize_OnDevice.hpp
/// \brief Symmetrize a matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct Symmetrize<Uplo::Upper, Algo::OnDevice> {
  template <typename ViewTypeA> inline static int invoke(const ViewTypeA &A) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(std::is_same<typename ViewTypeA::memory_space, Kokkos::HostSpace>::value,
                  "A is not accessible from host");

    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(j, i) = A(i, j);
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }

  template <typename ExecSpaceType, typename ViewTypeA>
  inline static int invoke(ExecSpaceType &exec_instance, const ViewTypeA &A) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
        using exec_space = ExecSpaceType;
        const Kokkos::RangePolicy<exec_space> policy(exec_instance, 0, m * m);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &ij) {
              const ordinal_type i = ij % m;
              const ordinal_type j = ij / m;
              if (i < j)
                A(j, i) = A(i, j);
            });
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }
};

template <> struct Symmetrize<Uplo::Lower, Algo::OnDevice> {
  template <typename ViewTypeA> inline static int invoke(const ViewTypeA &A) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(std::is_same<typename ViewTypeA::memory_space, Kokkos::HostSpace>::value,
                  "A is not accessible from host");

    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(i, j) = A(j, i);
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }

  template <typename ExecSpaceType, typename ViewTypeA>
  inline static int invoke(ExecSpaceType &exec_instance, const ViewTypeA &A) {
    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
        using exec_space = ExecSpaceType;
        const Kokkos::RangePolicy<exec_space> policy(exec_instance, 0, m * m);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &ij) {
              const ordinal_type i = ij % m;
              const ordinal_type j = ij / m;
              if (i < j)
                A(i, j) = A(j, i);
            });
      }
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
