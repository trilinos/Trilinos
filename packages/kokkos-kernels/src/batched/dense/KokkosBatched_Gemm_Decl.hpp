//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
#ifndef __KOKKOSBATCHED_GEMM_DECL_HPP__
#define __KOKKOSBATCHED_GEMM_DECL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

// Includes for non-functor-level routines
#include <KokkosBatched_Gemm_Handle.hpp>
#include <KokkosKernels_ExecSpaceUtils.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBatched {
/********************* BEGIN functor-level routines *********************/
///
/// Serial Gemm
///

template <typename ArgTransA, typename ArgTransB, typename ArgAlgo>
struct SerialGemm {
  template <typename ScalarType, typename AViewType, typename BViewType,
            typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha,
                                           const AViewType &A,
                                           const BViewType &B,
                                           const ScalarType beta,
                                           const CViewType &C);
};

///
/// Team Gemm
///

template <typename MemberType, typename ArgTransA, typename ArgTransB,
          typename ArgAlgo>
struct TeamGemm {
  template <typename ScalarType, typename AViewType, typename BViewType,
            typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const BViewType &B, const ScalarType beta, const CViewType &C);
};

///
/// TeamVector Gemm
///

template <typename MemberType, typename ArgTransA, typename ArgTransB,
          typename ArgAlgo>
struct TeamVectorGemm {
  template <typename ScalarType, typename AViewType, typename BViewType,
            typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const BViewType &B, const ScalarType beta, const CViewType &C);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgTransA, typename ArgTransB,
          typename ArgMode, typename ArgAlgo>
struct Gemm {
  template <typename ScalarType, typename AViewType, typename BViewType,
            typename CViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const BViewType &B, const ScalarType beta, const CViewType &C) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialGemm<ArgTransA, ArgTransB, ArgAlgo>::invoke(alpha, A, B,
                                                                beta, C);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamGemm<MemberType, ArgTransA, ArgTransB, ArgAlgo>::invoke(
          member, alpha, A, B, beta, C);
    }
    return r_val;
  }
};
/********************* END functor-level routines *********************/

/********************* BEGIN non-functor-level routines *********************/

namespace Impl {
/********************* BEGIN forward declarations *********************/
// clang-format off
/// \brief Non-blocking solve of general matrix multiply on a batch of
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
///                             BatchSzDim::Right Batch dimension is rightmost
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
template <class ArgTransA, class ArgTransB, class ArgMode, class ArgBatchSzDim,
          class ArgResultsPerThread, class ScalarType, class AViewType,
          class BViewType, class CViewType>
class BatchedSerialGemm;

// clang-format off
/// \brief Non-blocking solve of general matrix multiply on a batch of
/// uniform matrices with an algorithm based on:
///   B. P. D. J. Kunkel, Julian, “Performance, design, and autotuning of batched gemm for GPUs,”
///   in Lecture Notes in Computer Science, ser. ISC High Performance Computing ’16, vol. 9697, 06 2016.
///
///
///        C = alpha * op(A) * op(B) + beta * C
///
/// \tparam ArgTransA           Specifies what op does to A:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose (unsupported)
/// \tparam ArgTransB           Specifies what op does to B:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose (unsupported)
/// \tparam ArgBatchSzDim       Specifies where the batch dimension is allocated in
///                             AViewType, BViewType, and CViewType:
///                             BatchSzDim::Left  Batch dimension is leftmost
///                             BatchSzDim::Right Batch dimension is rightmost
/// \tparam ArgResultsPerThread Specifies how to divide work among threads. For
///                             this serial interface, each rank specifies how
///                             much work to assign a single thread.
///                             ResultsPerThread::Rank0 Each thread computes a scalar of C
///                             ResultsPerThread::Rank1 Each thread computes a 1-rank chunk of C
///                             ResultsPerThread::Rank2 Each thread computes a 2-rank chunk of C
/// \tparam HandleType          Specifies the handle type of the kernel handle
/// \tparam ScalarType          Specifies the scalar type of alpha and beta
/// \tparam AViewType           Input matrix, as either a 3-rank Kokkos::View or a
///                             4-rank Kokkos::View for SIMD operations.
/// \tparam BViewType           Input matrix, as either a 3-rank Kokkos::View or a
///                             4-rank Kokkos::View for SIMD operations.
/// \tparam CViewType           Input(RHS)/Output(LHS) matrix, as either a 3-rank
///                             Kokkos::View or a 4-rank Kokkos::View for SIMD
///                             operations.
/// \tparam ArgBoundsCheck      Specifies whether to perform global memory access
///                             bounds checks within the functor. Bounds checks
///                             are required when matrix sizes are not evenly divisible
///                             by tile sizes.
///                             BoundsCheck::Yes The functor will     perform bound checks (recommended)
///                             BoundsCheck::No  The functor will NOT perform bound checks
/// \tparam ArgAlphaFmaTag      Specifies whether to apply alpha during fmas.
///                             AlphaFmaTag::Yes alpha will be applied during fma (C = C * alpha + AB).
///                             AlphaFmaTag::No  alpha will be applied during mul (A * B * alpha).
/// \tparam TILE_M              Specifies the number of rows in each tile.
/// \tparam TILE_N              Specifies the number of cols in each tile.
/// \tparam TILE_K              Specifies the number of cols or rows in a tile of A or tile of B, respectively.
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
///                     ScalarType, AViewType, BViewType, CViewType
///                     ArgBoundsCheck, tile_m, tile_n, tile_k>(alpha, A, B, beta, C).invoke();
// clang-format on
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim,
          class HandleType, class ScalarType, class AViewType, class BViewType,
          class CViewType, class ArgBoundsCheck, class ArgAlphaFmaTag,
          int tile_m, int tile_n, int tile_k>
class BatchedDblBufGemm;

// clang-format off
/// \brief Blocking solve of general matrix multiply on a batch of uniform matrices.
///
///
///        C = alpha * op(A) * op(B) + beta * C
///
/// \tparam ArgTransA           Specifies what op does to A:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose (unsupported)
/// \tparam ArgTransB           Specifies what op does to B:
///                             Trans::NoTranspose   for non-transpose
///                             Trans::Transpose     for transpose
///                             Trans::ConjTranspose for conjugate transpose (unsupported)
/// \tparam HandleType          Specifies the handle type of the kernel handle
/// \tparam ScalarType          Specifies the scalar type of alpha and beta
/// \tparam AViewType           Input matrix, as a 3-rank Kokkos::View
/// \tparam BViewType           Input matrix, as a 3-rank Kokkos::View
/// \tparam CViewType           Input(RHS)/Output(LHS) matrix, as a 3-rank
///                             Kokkos::View
///
///                             See struct BatchedGemmHandle for details
/// \param handle [in]          A handle which specifies how to invoke the batched
///                             gemm. handle->get_tpl_params() returns &ninter.
///                             ninter: The number of matrices to interleave.
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
///   BatchedArmplGemm<ArgTransA, ArgTransB, ArgBatchSzDim, HandleType,
///                     ScalarType, AViewType, BViewType, CViewType>
///                     (handle, alpha, A, B, beta, C).invoke();
// clang-format on
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim,
          class HandleType, class ScalarType, class AViewType, class BViewType,
          class CViewType>
class BatchedArmplGemm;
/********************* END forward declarations *********************/
}  // namespace Impl

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
///                        Trans::NoTranspose   for non-transpose
///                        Trans::Transpose     for transpose
///                        Trans::ConjTranspose for conjugate transpose
/// \tparam ArgTransB      Specifies what op does to B:
///                        Trans::NoTranspose   for non-transpose
///                        Trans::Transpose     for transpose
///                        Trans::ConjTranspose for conjugate transpose
/// \tparam ArgBatchSzDim  Specifies where the batch dimension is allocated in
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
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix A is MxKxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix A is BxMxK
/// \param B [in]          Input matrix, as a 3-rank Kokkos::View
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix B is KxNxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix B is BxKxN
/// \param beta [in]       Input coefficient used for multiplication with C
/// \param C [in/out]      Input/Output matrix, as a 3-rank Kokkos::View
///                        If ArgBatchSzDim == "BatchLayout::Right", matrix C is MxNxB
///                        If ArgBatchSzDim == "BatchLayout::Left",  matrix C is BxMxN
/// \return 0 upon success, non-zero otherwise
///
/// Usage Example:
///   BatchedGemm<ArgTransA, ArgTransB,
///               ArgBatchSzDim>(handle, alpha, A, B, beta, C);
// clang-format on
template <typename ArgTransA, typename ArgTransB, typename ArgBatchSzDim,
          typename BatchedGemmHandleType, typename ScalarType,
          typename AViewType, typename BViewType, typename CViewType>
int BatchedGemm(BatchedGemmHandleType *const handle, const ScalarType alpha,
                const AViewType &A, const BViewType &B, const ScalarType beta,
                const CViewType &C) {
  int ret = 0;
  size_t c_m, c_n;
  using ViewValueType = typename CViewType::value_type;
  // Check for valid input views
  static_assert(Kokkos::is_view<AViewType>::value,
                "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value,
                "BViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<CViewType>::value,
                "CViewType must be a Kokkos::View.");
  static_assert(
      std::is_same<ArgTransA, Trans::NoTranspose>::value ||
          std::is_same<ArgTransA, Trans::Transpose>::value,
      "ArgTransA must be either Trans::Transpose or Trans::NoTranspose.");
  static_assert(
      std::is_same<ArgTransB, Trans::NoTranspose>::value ||
          std::is_same<ArgTransB, Trans::Transpose>::value,
      "ArgTransB must be either Trans::Transpose or Trans::NoTranspose.");
  if (is_vector<ViewValueType>::value) {
    // Check ranks of view with underlying SIMD value types
    // For SIMD views, we can have either 3-rank or 4-ranks inputs.
    switch (handle->get_kernel_algo_type()) {
      case BaseKokkosBatchedAlgos::KK_SERIAL:
      case BaseHeuristicAlgos::SQUARE:
      case BaseTplAlgos::ARMPL:
        static_assert(static_cast<int>(AViewType::rank) == 3,
                      "AViewType must have rank 3.");
        static_assert(static_cast<int>(BViewType::rank) == 3,
                      "BViewType must have rank 3.");
        static_assert(static_cast<int>(CViewType::rank) == 3,
                      "CViewType must have rank 3.");
        break;

        // TODO: check this once KK_TEAM is supported
        //        case GemmKokkosBatchedAlgos::KK_TEAM:
        //          static_assert(static_cast<int>(AViewType::rank) == 4,
        //                        "AViewType must have rank 4.");
        //          static_assert(static_cast<int>(BViewType::rank) == 4,
        //                        "BViewType must have rank 4.");
        //          static_assert(static_cast<int>(CViewType::rank) == 4,
        //                        "CViewType must have rank 4.");
        //          break;

      default:
        std::ostringstream os;
        os << "KokkosBatched::BatchedGemm does not support kernelAlgoType = "
           << std::to_string(handle->get_kernel_algo_type())
           << " with SIMD views." << std::endl;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
        break;
    }
  } else {
    // Check ranks of views with underlying scalar value types
    static_assert(static_cast<int>(AViewType::rank) == 3,
                  "AViewType must have rank 3.");
    static_assert(static_cast<int>(BViewType::rank) == 3,
                  "BViewType must have rank 3.");
    static_assert(static_cast<int>(CViewType::rank) == 3,
                  "CViewType must have rank 3.");
  }

  // Check for valid data access patterns
  // Skip checking a_layout == b_layout == c_layout
  // Skip checking for LayoutStride
  using c_layout = typename CViewType::array_layout;
  if (std::is_same<c_layout, Kokkos::LayoutLeft>::value &&
      !std::is_same<ArgBatchSzDim, BatchLayout::Right>::value) {
    throw std::runtime_error(
        "Error: LayoutLeft views require BatchLayout::Right");
  }
  if (std::is_same<c_layout, Kokkos::LayoutRight>::value &&
      !std::is_same<ArgBatchSzDim, BatchLayout::Left>::value) {
    throw std::runtime_error(
        "Error: LayoutRight views require BatchLayout::Left");
  }

  if (std::is_same<ArgBatchSzDim, BatchLayout::Left>::value) {
    // c_b = C.extent(0);
    c_m = C.extent(1);
    c_n = C.extent(2);
  } else {
    // c_b = C.extent(2);
    c_m = C.extent(0);
    c_n = C.extent(1);
  }

  // Begin checking conditions for optimal BatchedGemm invocation.
  using view_scalar_type   = typename CViewType::value_type;
  using layout_type        = typename CViewType::array_layout;
  constexpr bool is_vector = KokkosBatched::is_vector<view_scalar_type>::value;
  constexpr bool on_gpu    = KokkosKernels::Impl::kk_is_gpu_exec_space<
      typename CViewType::execution_space>();
  constexpr bool on_x86_64 = KokkosKernels::Impl::kk_is_x86_64_mem_space<
      typename CViewType::execution_space::memory_space>();
  constexpr bool on_a64fx = KokkosKernels::Impl::kk_is_a64fx_mem_space<
      typename CViewType::execution_space::memory_space>();

  if (handle->enableDebug) {
    std::cout << "view_scalar_type:" << typeid(view_scalar_type).name()
              << std::endl
              << "execution_space:"
              << typeid(typename CViewType::execution_space).name() << std::endl
              << std::endl
              << "is_vector:" << is_vector << std::endl
              << "on_gpu:" << on_gpu << std::endl
              << "on_x86_64:" << on_x86_64 << std::endl
              << "on_a64fx:" << on_a64fx << std::endl;
  }

  switch (handle->get_kernel_algo_type()) {
    ////////////// HEURISTIC ALGOS //////////////
    case BaseHeuristicAlgos::SQUARE:
      if (c_m != c_n) {
        std::ostringstream os;
        os << "KokkosBatched::BatchedGemm does not support kernelAlgoType = "
           << std::to_string(handle->get_kernel_algo_type()) << " when c_m("
           << std::to_string(c_m) << ") != c_n(" << std::to_string(c_n) << ")"
           << std::endl;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }

      // Select optimal resultsPerThread param for BatchedSerialGemm
      using bsgResultsPerThread =
          typename std::conditional<!is_vector && on_gpu,
                                    ResultsPerThread::Rank0,
                                    ResultsPerThread::Rank2>::type;

      // Select optimal mode param for SerialGemm.
      using bsgModeType = typename std::conditional<
          is_vector,
          typename std::conditional<on_gpu || on_x86_64, Algo::Gemm::Blocked,
                                    Algo::Gemm::Unblocked>::type,
          typename std::conditional<
              on_gpu, Algo::Gemm::Unblocked,
              typename std::conditional<on_a64fx, Algo::Gemm::Unblocked,
                                        Algo::Gemm::Blocked>::type>::type>::
          type;

      if (handle->enableDebug) {
        std::cout << "bsgResultsPerThread: "
                  << typeid(bsgResultsPerThread).name() << std::endl
                  << "bsgModeType: " << typeid(bsgModeType).name() << std::endl;
      }

      // if (on_gpu && c_m >= 20 &&
      //     (alpha == 1.0F && beta == 0.0F) ? c_m <= 24 : c_m <= 21) {
      //   // TODO: invoke TeamShmem
      // } else
      if (on_gpu && ((std::is_same<layout_type, Kokkos::LayoutLeft>::value)
                         ? (c_m >= 16)
                         : (c_m >= 24 && c_m <= 32) || c_m >= 40)) {
        handle->teamSz = handle->vecLen = 8;
        constexpr int tile_m = 32, tile_n = 32, tile_k = 8;
#ifdef __CUDACC_RDC__
        constexpr size_t alpha_in_fma_thresh = 24;
#else
        constexpr size_t alpha_in_fma_thresh = 64;
#endif  // __CUDAACC_RDC__

        if (c_m % 32 == 0) {                 // No bounds checking
          if (c_m >= alpha_in_fma_thresh) {  // apply alpha in fma
            ret =
                Impl::BatchedDblBufGemm<
                    ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                    ScalarType, AViewType, BViewType, CViewType,
                    BoundsCheck::No, AlphaTag::Yes, tile_m, tile_n, tile_k>(
                    handle, alpha, A, B, beta, C)
                    .invoke();
          } else {  // apply alpha in mul
            ret =
                Impl::BatchedDblBufGemm<
                    ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                    ScalarType, AViewType, BViewType, CViewType,
                    BoundsCheck::No, AlphaTag::No, tile_m, tile_n, tile_k>(
                    handle, alpha, A, B, beta, C)
                    .invoke();
          }
        } else {                             // bounds checking
          if (c_m >= alpha_in_fma_thresh) {  // apply alpha in fma
            ret =
                Impl::BatchedDblBufGemm<
                    ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                    ScalarType, AViewType, BViewType, CViewType,
                    BoundsCheck::Yes, AlphaTag::Yes, tile_m, tile_n, tile_k>(
                    handle, alpha, A, B, beta, C)
                    .invoke();
          } else {  // apply alpha in mul
            ret =
                Impl::BatchedDblBufGemm<
                    ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                    ScalarType, AViewType, BViewType, CViewType,
                    BoundsCheck::Yes, AlphaTag::No, tile_m, tile_n, tile_k>(
                    handle, alpha, A, B, beta, C)
                    .invoke();
          }
        }
      } else {
        ret = Impl::BatchedSerialGemm<ArgTransA, ArgTransB, bsgModeType,
                                      ArgBatchSzDim, bsgResultsPerThread,
                                      ScalarType, AViewType, BViewType,
                                      CViewType>(alpha, A, B, beta, C)
                  .invoke();
      }
      break;

      //    case BaseHeuristicAlgos::TALL:
      //
      //    case BaseHeuristicAlgos::WIDE:
      ////////////// TPL ALGOS //////////////
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
    case BaseTplAlgos::ARMPL:
      ret = Impl::BatchedArmplGemm<ArgTransA, ArgTransB, ArgBatchSzDim,
                                   BatchedGemmHandleType, ScalarType, AViewType,
                                   BViewType, CViewType>(handle, alpha, A, B,
                                                         beta, C)
                .invoke();
      break;
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
      //    case BaseTplAlgos::MKL:
      //
      //    case GemmTplAlgos::CUBLAS:
      //
      //    case GemmTplAlgos::MAGMA:

    ////////////// KokkosBatched ALGOS //////////////
    case BaseKokkosBatchedAlgos::KK_SERIAL:
      ret =
          Impl::BatchedSerialGemm<ArgTransA, ArgTransB, Algo::Gemm::Unblocked,
                                  ArgBatchSzDim, ResultsPerThread::Rank2,
                                  ScalarType, AViewType, BViewType, CViewType>(
              alpha, A, B, beta, C)
              .invoke();
      break;

      // case GemmKokkosBatchedAlgos::KK_SERIALSIMD:

    case GemmKokkosBatchedAlgos::KK_SERIAL_RANK0:
      ret =
          Impl::BatchedSerialGemm<ArgTransA, ArgTransB, Algo::Gemm::Unblocked,
                                  ArgBatchSzDim, ResultsPerThread::Rank0,
                                  ScalarType, AViewType, BViewType, CViewType>(
              alpha, A, B, beta, C)
              .invoke();
      break;

      //    case GemmKokkosBatchedAlgos::KK_SERIAL_SHMEM:
      //    case GemmKokkosBatchedAlgos::KK_TEAM:
      //    case GemmKokkosBatchedAlgos::KK_TEAMVECTOR:
      //    case GemmKokkosBatchedAlgos::KK_TEAMSIMD:

    case GemmKokkosBatchedAlgos::KK_DBLBUF:
      // Note: The tile sizes of 1x1x1 here will not perform well but must be
      // selected in order to function on all devices since the serial execution
      // space has a max team size of 1. KokkosKernels API users will need to
      // follow an approach similar to KK_SQUARE above for best performance.

      // TODO: Add auto-selection of tile size based on inputs and device type
      ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim,
                                    BatchedGemmHandleType, ScalarType,
                                    AViewType, BViewType, CViewType,
                                    BoundsCheck::Yes, AlphaTag::No, 1, 1, 1>(
                handle, alpha, A, B, beta, C)
                .invoke();
      break;

    default:
      std::ostringstream os;
      os << "KokkosBatched::BatchedGemm does not support kernelAlgoType = "
         << std::to_string(handle->get_kernel_algo_type()) << "." << std::endl;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
      break;
  }
  return ret;
}
/********************* END non-functor-level routines *********************/
}  // namespace KokkosBatched

#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"
#include "KokkosBatched_Gemm_DblBuf_Impl.hpp"
#include "KokkosBatched_Gemm_Armpl_Impl.hpp"

#endif
