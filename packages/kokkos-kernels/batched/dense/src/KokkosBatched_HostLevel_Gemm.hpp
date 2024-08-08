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
#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_DECL_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_DECL_HPP__

// Include explicit specializations of BatchedGemm.
// If ETI_ONLY is disabled, the primary template will
// be inlined into each caller's invocation using non-
// ETI'd template arguments.
#include "KokkosBatched_HostLevel_Gemm_Spec.hpp"

namespace KokkosBatched {
// clang-format off
/// \brief Non-blocking solve of general matrix multiply on a batch of
/// uniform matrices.
///
/// Note: If a TPL is selected, this interface follows the blocking
/// behavior (either blocking or non-blocking) of the TPL vendor's API.
///
/// Note: To leverage SIMD instructions, 4-rank views must be selected via the
/// template parameters documented below.
///
///        C = alpha * op(A) * op(B) + beta * C
///
/// \tparam ArgTransA      Specifies what op does to A:
///
///                        Trans::NoTranspose   for non-transpose
///                        Trans::Transpose     for transpose
///                        Trans::ConjTranspose for conjugate transpose
/// \tparam ArgTransB      Specifies what op does to B:
///
///                        Trans::NoTranspose   for non-transpose
///                        Trans::Transpose     for transpose
///                        Trans::ConjTranspose for conjugate transpose
/// \tparam ArgBatchSzDim  Specifies where the batch dimension is allocated in
///
///                        AViewType, BViewType, and CViewType:
///                        BatchLayout::Left  Batch dimension is leftmost
///                        BatchLayout::Right Batch dimension is rightmost
/// \tparam ScalarType     Specifies the scalar type of alpha and beta
/// \tparam AViewType      Input matrix, as either a 3-rank Kokkos::View or a
///                        4-rank Kokkos::View for SIMD operations.
/// \tparam BViewType      Input matrix, as either a 3-rank Kokkos::View or a
///                        4-rank Kokkos::View for SIMD operations.
/// \tparam CViewType      Input(RHS)/Output(LHS) matrix, as either a 3-rank
///                        Kokkos::View or a 4-rank Kokkos::View for SIMD
///                        operations.
///
/// \param handle [in]     A handle which specifies how to invoke the batched
///                        gemm.
///                        See struct BatchedGemmHandle for details.
/// \param alpha [in]      Input coefficient used for multiplication with A
/// \param A [in]          Input matrix, as a 3-rank Kokkos::View
///
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix A is MxKxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix A is BxMxK
/// \param B [in]          Input matrix, as a 3-rank Kokkos::View
///
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix B is KxNxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix B is BxKxN
/// \param beta [in]       Input coefficient used for multiplication with C
/// \param C [in/out]      Input/Output matrix, as a 3-rank Kokkos::View
///
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix C is MxNxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix C is BxMxN
/// \return 0 upon success, non-zero otherwise
///
/// Usage Example:
///   BatchedGemm<ArgTransA, ArgTransB,
///               ArgBatchSzDim>(handle, alpha, A, B, beta, C);
// clang-format on
template <typename ArgTransA, typename ArgTransB, typename ArgBatchSzDim, typename BatchedGemmHandleType,
          typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
inline int BatchedGemm(BatchedGemmHandleType *const handle, const ScalarType alpha, const AViewType &A,
                       const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Minimize the number of ImplBatchedGemmWrapper instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  using UnifiedAVT = Kokkos::View<typename AViewType::value_type ***, typename AViewType::array_layout,
                                  typename AViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UnifiedBVT = Kokkos::View<typename BViewType::value_type ***, typename BViewType::array_layout,
                                  typename BViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UnifiedCVT = Kokkos::View<typename CViewType::non_const_value_type ***, typename CViewType::array_layout,
                                  typename CViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Go through specialization layer in case ETI'd symbols are available.
  return Impl::BatchedGemmSpec<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType, UnifiedAVT,
                               UnifiedBVT, UnifiedCVT>::run(handle, alpha, A, B, beta, C);
}
}  // namespace KokkosBatched
#endif  // __KOKKOSBATCHED_HOSTLEVEL_GEMM_DECL_HPP__
