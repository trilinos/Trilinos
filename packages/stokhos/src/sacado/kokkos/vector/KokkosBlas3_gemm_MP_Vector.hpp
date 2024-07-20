// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOSBLAS3_GEMM_MP_VECTOR_HPP
#define KOKKOSBLAS3_GEMM_MP_VECTOR_HPP

#include <type_traits>
#include "Sacado_ConfigDefs.h"

#include "Stokhos_ViewStorage.hpp"
#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "KokkosBlas.hpp"

namespace KokkosBlas
{
template <typename DA, typename... PA,
          typename DB, typename... PB,
          typename DC, typename... PC>
typename std::enable_if<Kokkos::is_view_mp_vector<Kokkos::View<DA, PA...>>::value &&
                        Kokkos::is_view_mp_vector<Kokkos::View<DB, PB...>>::value &&
                        Kokkos::is_view_mp_vector<Kokkos::View<DC, PC...>>::value>::type
gemm(const char transA[],
     const char transB[],
     typename Kokkos::View<DA, PA...>::const_value_type &alpha,
     const Kokkos::View<DA, PA...> &A,
     const Kokkos::View<DB, PB...> &B,
     typename Kokkos::View<DC, PC...>::const_value_type &beta,
     const Kokkos::View<DC, PC...> &C)
{
  // Assert that A, B, and C are in fact matrices
  static_assert(Kokkos::View<DA, PA...>::rank == 2, "GEMM: A must have rank 2 (be a matrix).");
  static_assert(Kokkos::View<DB, PB...>::rank == 2, "GEMM: B must have rank 2 (be a matrix).");
  static_assert(Kokkos::View<DC, PC...>::rank == 2, "GEMM: C must have rank 2 (be a matrix).");

  if (B.extent(1) == 1 && C.extent(1) == 1)
  {
    auto x = Kokkos::subview(B, Kokkos::ALL, 0);
    auto y = Kokkos::subview(C, Kokkos::ALL, 0);
    KokkosBlas::gemv(transA, alpha, A, x, beta, y);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "GEMM: Not implemented for Sacado::MP::Vector scalar type!");
}
} // namespace KokkosBlas

#endif
