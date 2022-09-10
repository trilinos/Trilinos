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
#ifndef __TACHO_SYMMETRIZE_INTERNAL_HPP__
#define __TACHO_SYMMETRIZE_INTERNAL_HPP__

/// \file  Tacho_Symmetrize_Internal.hpp
/// \brief Symmetrize a square block matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct Symmetrize<Uplo::Upper, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A) {
    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j), [&](const ordinal_type &i) { A(j, i) = A(i, j); });
        });
#else
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(j, i) = A(i, j);
#endif
      }
    } else {
      printf("Error: Symmetrize<Algo::Internal> A is not square\n");
    }
    return 0;
  }
};

template <> struct Symmetrize<Uplo::Lower, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A) {
    const ordinal_type m = A.extent(0), n = A.extent(1);

    if (m == n) {
      if (A.span() > 0) {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j), [&](const ordinal_type &i) { A(i, j) = A(j, i); });
        });
#else
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(i, j) = A(j, i);
#endif
      }
    } else {
      printf("Error: Symmetrize<Algo::Internal> A is not square\n");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
