//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef __KOKKOSBATCHED_SVD_DECL_HPP__
#define __KOKKOSBATCHED_SVD_DECL_HPP__

/// \author Brian Kelley (bmkelle@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// Given a general matrix A (m x n), compute the full singular value
/// decomposition (SVD): U * diag(s) * V^T = A. U/V are orthogonal and s
/// contains nonnegative values in descending order.
///
/// Currently only supports real-valued matrices.
///
/// Parameters:
///   [in] A
///     General matrix (rank 2 view), m x n.
///     The contents of A are overwritten and undefined after calling this
///     function.
///   [out] U
///     m left singular vectors (in columns). Dimensions m*m.
///   [out] Vt
///     n right singular vectors (in rows). Dimensions n*n.
///   [out] s
///     min(m, n) singular values.
///   [in] W
///     1D contiguous workspace. The required size is max(m, n).
///
/// Preconditions:
///   m == A.extent(0) == U.extent(0) == U.extent(1)
///   n == A.extent(1) == V.extent(0) == V.extent(1)
///   min(m, n) == s.extent(0)
///   W.extent(0) >= max(m, n)
///   W.stride(0) == 1 (contiguous)

struct SVD_USV_Tag {};
struct SVD_S_Tag {};
// Note: Could easily add SV or US tags later if needed

struct SerialSVD {
  // Version to compute full factorization: A == U * diag(s) * Vt
  template <typename AViewType, typename UViewType, typename VtViewType, typename SViewType, typename WViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      SVD_USV_Tag, const AViewType &A, const UViewType &U, const SViewType &s, const VtViewType &Vt, const WViewType &W,
      typename AViewType::const_value_type tol = Kokkos::ArithTraits<typename AViewType::value_type>::zero());

  // Version which computes only singular values
  template <typename AViewType, typename SViewType, typename WViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      SVD_S_Tag, const AViewType &A, const SViewType &s, const WViewType &W,
      typename AViewType::const_value_type tol = Kokkos::ArithTraits<typename AViewType::value_type>::zero());
};

}  // namespace KokkosBatched

#include "KokkosBatched_SVD_Serial_Impl.hpp"

#endif
