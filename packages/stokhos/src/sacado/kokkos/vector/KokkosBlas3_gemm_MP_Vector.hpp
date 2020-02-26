// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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
