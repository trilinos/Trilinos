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
#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_SERIAL_IMPL_HPP__
#include "KokkosBatched_Gemm_Decl.hpp"

namespace KokkosBatched {
namespace Impl {
// clang-format off
/// \brief Non-blocking general matrix multiply on a batch of
/// uniform matrices.
///
///
///        C = alpha * op(A) * op(B) + beta * C
///
/// \tparam ArgTransA           Specifies what op does to A:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose
/// \tparam ArgTransB           Specifies what op does to B:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose
/// \tparam ArgMode             Specifies algorithm mode to use for serial work:
///                             Algo::Gemm::Unblocked  for no register blocking
///                             Algo::Gemm::Blocked    for register blocking
///                             Algo::Gemm::CompactMKL for mkl compact tpl interface
/// \tparam ArgBatchSzDim       Specifies where the batch dimension is allocated in
///                             AViewType, BViewType, and CViewType:
///                             BatchSzDim::Left  Batch dimension is leftmost
///                             BatchSzDim::Right  Batch dimension is rightmost
/// \tparam ArgResultsPerThread Specifies how to divide work among threads. For
///                             this serial interface, each rank specifies how
///                             much work to assign a single thread.
///                             ResultsPerThread::Rank0 Each thread computes a scalar of C
///                             ResultsPerThread::Rank1 Each thread computes a 1-rank chunk of C
///                             ResultsPerThread::Rank2 Each thread computes a 2-rank chunk of C
/// \tparam ScalarType          Specifies the scalar type of alpha and beta
/// \tparam AViewType           Input matrix, as either a 3-rank Kokkos::View or a
///                             4-rank Kokkos::View for SIMD operations.
/// \tparam BViewType           Input matrix, as either a 3-rank Kokkos::View or a
///                             4-rank Kokkos::View for SIMD operations.
/// \tparam CViewType           Input(RHS)/Output(LHS) matrix, as either a 3-rank
///                             Kokkos::View or a 4-rank Kokkos::View for SIMD
///                             operations.
///
///                             See struct BatchedGemmHandle for details.
/// \param alpha [in]           Input coefficient used for multiplication with A
/// \param A [in]               Input matrix, as a 3-rank Kokkos::View
///                             If ArgBatchSzDim == "BatchSzDim::Right", matrix A is MxKxB
///                             If ArgBatchSzDim == "BatchSzDim::Left",  matrix A is BxMxK
/// \param B [in]               Input matrix, as a 3-rank Kokkos::View
///                             If ArgBatchSzDim == "BatchSzDim::Right", matrix B is KxNxB
///                             If ArgBatchSzDim == "BatchSzDim::Left",  matrix B is BxKxN
/// \param beta [in]            Input coefficient used for multiplication with C
/// \param C [in/out]           Input/Output matrix, as a 3-rank Kokkos::View
///                             If ArgBatchSzDim == "BatchSzDim::Right", matrix C is MxNxB
///                             If ArgBatchSzDim == "BatchSzDim::Left",  matrix C is BxMxN
/// \return 0 upon success, non-zero otherwise
///
/// Usage Example:
///   BatchedSerialGemm<ArgTransA, ArgTransB, ArgMode, ArgBatchSzDim,
///                     ArgResultsPerThread, ScalarType, AViewType,
///                     BViewType, CViewType>(alpha, A, B, beta, C).invoke();
// clang-format on
template <class ArgTransA, class ArgTransB, class ArgMode, class ArgBatchSzDim, class ArgResultsPerThread,
          class ScalarType, class AViewType, class BViewType, class CViewType>
class BatchedSerialGemm {
 private:
  AViewType A;
  BViewType B;
  CViewType C;
  ScalarType alpha, beta;
  size_t divisor, c_cols, batch_size;
  ArgBatchSzDim batch_layout_tag;
  ArgTransA transA_tag;
  ArgTransB transB_tag;

  void run() {
    using execution_space = typename CViewType::device_type::execution_space;
    using policy_type     = Kokkos::RangePolicy<ArgResultsPerThread, execution_space>;
    Kokkos::parallel_for("BatchedSerialGemm", policy_type(0, batch_size), *this);
  }

 public:
  int invoke() {
    if (std::is_same<ArgResultsPerThread, ResultsPerThread::Rank0>::value) {
      // Set members for ResultsPerThread::Rank0 operator; these members allow
      // each thread to calculate its C output index
      if (std::is_same<ArgBatchSzDim, BatchLayout::Left>::value) {
        batch_size = C.extent(0);
        divisor    = C.extent(1) * C.extent(2);
        c_cols     = C.extent(2);
      } else {
        batch_size = C.extent(2);
        divisor    = C.extent(0) * C.extent(1);
        c_cols     = C.extent(1);
      }

      // Increase the number of threads by the divisor
      batch_size *= divisor;

      run();
    } else if (std::is_same<ArgResultsPerThread, ResultsPerThread::Rank2>::value) {
      if (std::is_same<ArgBatchSzDim, BatchLayout::Left>::value)
        batch_size = C.extent(0);
      else
        batch_size = C.extent(2);

      run();
    } else {
      std::cerr << "Error: ArgResultsPerThread not supported" << std::endl;
      return -1;
    }
    return 0;
  }

  BatchedSerialGemm(ScalarType _alpha, AViewType _A, BViewType _B, ScalarType _beta, CViewType _C)
      : A(_A), B(_B), C(_C), alpha(_alpha), beta(_beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ResultsPerThread::Rank0 &, const int &i) const {
    // Here, the batch_idx is strided by c_rows * c_cols
    auto batch_idx = i / divisor;
    // For every batch, we need mod in [0, c_rows*c_cols-1]
    auto mod = i % divisor;
    // For every mod, we need a column index in [0, c_cols-1]
    auto col_idx = mod % c_cols;
    // For every mod, we need a row index in [0, c_rows-1]
    auto row_idx = mod / c_cols;

    // Due to taking 1-rank subviews out, we must handle transpose here.
    // Use overloads of subview_wrapper to handle transpose at compile time.
    auto svA_row = subview_wrapper(A, batch_idx, row_idx, Kokkos::ALL(), batch_layout_tag, transA_tag);
    auto svB_col = subview_wrapper(B, batch_idx, Kokkos::ALL(), col_idx, batch_layout_tag, transB_tag);
    auto svC_ele = subview_wrapper(C, batch_idx, row_idx, col_idx, batch_layout_tag);

    // Kokkos::subview(scalar, ALL) or Kokkos::subview(ALL, scalar) always
    // returns a column vector. Since the subviews above handle the
    // matrix transpositions, here we must perform the GEMM on:
    // row_vec x col_vec, which is svA_row' x svB_col to compute the element
    // of C.
    KokkosBatched::SerialGemm<Trans::Transpose, Trans::NoTranspose, ArgMode>::invoke(alpha, svA_row, svB_col, beta,
                                                                                     svC_ele);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ResultsPerThread::Rank2 &, const int &i) const {
    auto svA = subview_wrapper(A, i, Kokkos::ALL(), Kokkos::ALL(), batch_layout_tag);
    auto svB = subview_wrapper(B, i, Kokkos::ALL(), Kokkos::ALL(), batch_layout_tag);
    auto svC = subview_wrapper(C, i, Kokkos::ALL(), Kokkos::ALL(), batch_layout_tag);

    KokkosBatched::SerialGemm<ArgTransA, ArgTransB, ArgMode>::invoke(alpha, svA, svB, beta, svC);
  }
};
}  // namespace Impl
}  // namespace KokkosBatched
#endif