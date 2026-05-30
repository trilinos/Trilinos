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
#ifndef __TACHO_SYMMETRIZE_INTERNAL_HPP__
#define __TACHO_SYMMETRIZE_INTERNAL_HPP__

/// \file  Tacho_Symmetrize_Internal.hpp
/// \brief Symmetrize a square block matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct Symmetrize<Uplo::Upper, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const bool conjugate) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using arith_traits = ArithTraits<value_type>;

    const ordinal_type m = A.extent(0), n = A.extent(1);
    if (m == n) {
      if (A.span() > 0) {
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j), [&](const ordinal_type &i) { A(j, i) = A(i, j); });
        });))
        KOKKOS_IF_ON_HOST((
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(j, i) = (conjugate ? arith_traits::conj(A(i, j)) : A(i, j));))
      }
    } else {
      Kokkos::printf("Error: Symmetrize<Algo::Internal> A is not square\n");
    }
    return 0;
  }
};

template <> struct Symmetrize<Uplo::Lower, Algo::Internal> {
  template <typename MemberType, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const bool conjugate) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using arith_traits = ArithTraits<value_type>;

    const ordinal_type m = A.extent(0), n = A.extent(1);
    if (m == n) {
      if (A.span() > 0) {
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j), [&](const ordinal_type &i) {
            A(i, j) = (conjugate ? arith_traits::conj(A(j, i)) : A(j, i)); });
        });))
        KOKKOS_IF_ON_HOST((
        for (ordinal_type j = 0; j < n; ++j)
          for (ordinal_type i = 0; i < j; ++i)
            A(i, j) = (conjugate ? arith_traits::conj(A(j, i)) : A(j, i));))
      }
    } else {
      Kokkos::printf("Error: Symmetrize<Algo::Internal> A is not square\n");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
