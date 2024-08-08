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
#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_IMPL_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_IMPL_HPP__
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Util.hpp>  // Trans, BatchLayout
#include <KokkosKernels_ExecSpaceUtils.hpp>
#include <KokkosKernels_Error.hpp>

#include "KokkosBatched_HostLevel_Gemm_Handle.hpp"  // BatchedGemmHandle
#include "KokkosBatched_HostLevel_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_HostLevel_Gemm_DblBuf_Impl.hpp"
#include "KokkosBatched_HostLevel_Gemm_Armpl_Impl.hpp"

namespace KokkosBatched {
namespace Impl {
//////////////////////////////// tile_m //////////////////////////////////
template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION int kk_gemm_dbl_buf_tile_m() {
  return 32;
}
//////////////////////////////// tile_n //////////////////////////////////
template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION int kk_gemm_dbl_buf_tile_n() {
  return 32;
}
//////////////////////////////// tile_k //////////////////////////////////
template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION int kk_gemm_dbl_buf_tile_k() {
  return 8;
}

// On MI100, batched_scalar_batched_gemm_nt_nt_dcomplex_dcomplex_right fails
// without this. See https://github.com/kokkos/kokkos-kernels/issues/1547.
// This reduces the register allocations (REG_M and REG_N) in the double
// buffering algorithm by a factor of 2.
#if defined(KOKKOS_ENABLE_HIP) && defined(KOKKOS_ARCH_VEGA908)
template <>
constexpr KOKKOS_INLINE_FUNCTION int kk_gemm_dbl_buf_tile_k<Kokkos::HIP>() {
  return 16;
}
#endif
////////////////////////// alpha_in_fma_thresh ////////////////////////////
constexpr KOKKOS_INLINE_FUNCTION size_t kk_gemm_dbl_buf_alpha_in_fma_thresh() {
#ifdef __CUDACC_RDC__
  return 24;
#else
  return 64;
#endif  // __CUDAACC_RDC__
}

template <typename ArgTransA, typename ArgTransB, typename ArgBatchSzDim, typename BatchedGemmHandleType,
          typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
int BatchedGemmImpl(BatchedGemmHandleType *const handle, const ScalarType alpha, const AViewType &A, const BViewType &B,
                    const ScalarType beta, const CViewType &C) {
  int ret = 0;
  size_t c_m, c_n;
  using ViewValueType = typename CViewType::value_type;
  // Check for valid input views
  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value, "BViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<CViewType>::value, "CViewType must be a Kokkos::View.");
  static_assert(std::is_same<ArgTransA, Trans::NoTranspose>::value || std::is_same<ArgTransA, Trans::Transpose>::value,
                "ArgTransA must be either Trans::Transpose or Trans::NoTranspose.");
  static_assert(std::is_same<ArgTransB, Trans::NoTranspose>::value || std::is_same<ArgTransB, Trans::Transpose>::value,
                "ArgTransB must be either Trans::Transpose or Trans::NoTranspose.");
  if constexpr (is_vector<ViewValueType>::value) {
    // Check ranks of view with underlying SIMD value types
    // For SIMD views, we can have either 3-rank or 4-ranks inputs.
    switch (handle->get_kernel_algo_type()) {
      case BaseKokkosBatchedAlgos::KK_SERIAL:
      case BaseHeuristicAlgos::SQUARE:
      case BaseTplAlgos::ARMPL:
        assert(A.rank_dynamic() == 3 && "AViewType must have rank 3.");
        assert(B.rank_dynamic() == 3 && "BViewType must have rank 3.");
        assert(C.rank_dynamic() == 3 && "CViewType must have rank 3.");
        break;
      default:
        std::ostringstream os;
        os << "KokkosBatched::BatchedGemm does not support kernelAlgoType = "
           << std::to_string(handle->get_kernel_algo_type()) << " with SIMD views." << std::endl;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
        break;
    }
  } else {
    // Check ranks of views with underlying scalar value types
    static_assert(static_cast<int>(AViewType::rank) == 3, "AViewType must have rank 3.");
    static_assert(static_cast<int>(BViewType::rank) == 3, "BViewType must have rank 3.");
    static_assert(static_cast<int>(CViewType::rank) == 3, "CViewType must have rank 3.");
  }

  // Check for valid data access patterns
  // Skip checking a_layout == b_layout == c_layout
  // Skip checking for LayoutStride
  using c_layout = typename CViewType::array_layout;
  static_assert(
      !(std::is_same<c_layout, Kokkos::LayoutLeft>::value && !std::is_same<ArgBatchSzDim, BatchLayout::Right>::value),
      "LayoutLeft views require BatchLayout::Right");
  static_assert(
      !(std::is_same<c_layout, Kokkos::LayoutRight>::value && !std::is_same<ArgBatchSzDim, BatchLayout::Left>::value),
      "LayoutRight views require BatchLayout::Left");

  if constexpr (std::is_same<ArgBatchSzDim, BatchLayout::Left>::value) {
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
  using exec_space         = typename CViewType::execution_space;
  constexpr bool is_vector = KokkosBatched::is_vector<view_scalar_type>::value;
  constexpr bool on_gpu    = KokkosKernels::Impl::kk_is_gpu_exec_space<exec_space>();
  constexpr bool on_x86_64 = KokkosKernels::Impl::kk_is_x86_64_mem_space<typename exec_space::memory_space>();
  constexpr bool on_a64fx  = KokkosKernels::Impl::kk_is_a64fx_mem_space<typename exec_space::memory_space>();
  bool out_of_range        = false;

  if (handle->enableDebug) {
    std::cout << "view_scalar_type:" << typeid(view_scalar_type).name() << std::endl
              << "execution_space:" << typeid(exec_space).name() << std::endl
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
           << std::to_string(handle->get_kernel_algo_type()) << " when c_m(" << std::to_string(c_m) << ") != c_n("
           << std::to_string(c_n) << ")" << std::endl;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }

      // Select optimal resultsPerThread param for BatchedSerialGemm
      using bsgResultsPerThread =
          std::conditional_t<!is_vector && on_gpu, ResultsPerThread::Rank0, ResultsPerThread::Rank2>;

      // Select optimal mode param for SerialGemm.
      using bsgModeType = typename std::conditional<
          is_vector, typename std::conditional<on_gpu || on_x86_64, Algo::Gemm::Blocked, Algo::Gemm::Unblocked>::type,
          typename std::conditional<
              on_gpu, Algo::Gemm::Unblocked,
              typename std::conditional<on_a64fx, Algo::Gemm::Unblocked, Algo::Gemm::Blocked>::type>::type>::type;

      if (handle->enableDebug) {
        std::cout << "bsgResultsPerThread: " << typeid(bsgResultsPerThread).name() << std::endl
                  << "bsgModeType: " << typeid(bsgModeType).name() << std::endl;
      }

      if constexpr (on_gpu) {
        if (((std::is_same<layout_type, Kokkos::LayoutLeft>::value) ? (c_m >= 16)
                                                                    : (c_m >= 24 && c_m <= 32) || c_m >= 40)) {
          handle->teamSz = handle->vecLen      = 8;
          constexpr int tile_m                 = Impl::kk_gemm_dbl_buf_tile_m<exec_space>();
          constexpr int tile_n                 = Impl::kk_gemm_dbl_buf_tile_n<exec_space>();
          constexpr int tile_k                 = Impl::kk_gemm_dbl_buf_tile_k<exec_space>();
          constexpr size_t alpha_in_fma_thresh = Impl::kk_gemm_dbl_buf_alpha_in_fma_thresh();

          if (c_m % 32 == 0) {                 // No bounds checking
            if (c_m >= alpha_in_fma_thresh) {  // apply alpha in fma
              ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType,
                                            AViewType, BViewType, CViewType, BoundsCheck::No, AlphaTag::Yes, tile_m,
                                            tile_n, tile_k>(handle, alpha, A, B, beta, C)
                        .invoke();
            } else {  // apply alpha in mul
              ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType,
                                            AViewType, BViewType, CViewType, BoundsCheck::No, AlphaTag::No, tile_m,
                                            tile_n, tile_k>(handle, alpha, A, B, beta, C)
                        .invoke();
            }
          } else {                             // bounds checking
            if (c_m >= alpha_in_fma_thresh) {  // apply alpha in fma
              ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType,
                                            AViewType, BViewType, CViewType, BoundsCheck::Yes, AlphaTag::Yes, tile_m,
                                            tile_n, tile_k>(handle, alpha, A, B, beta, C)
                        .invoke();
            } else {  // apply alpha in mul
              ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType,
                                            AViewType, BViewType, CViewType, BoundsCheck::Yes, AlphaTag::No, tile_m,
                                            tile_n, tile_k>(handle, alpha, A, B, beta, C)
                        .invoke();
            }
          }
        } else {
          out_of_range = true;
        }
      }
      if (!on_gpu || out_of_range) {
        ret = Impl::BatchedSerialGemm<ArgTransA, ArgTransB, bsgModeType, ArgBatchSzDim, bsgResultsPerThread, ScalarType,
                                      AViewType, BViewType, CViewType>(alpha, A, B, beta, C)
                  .invoke();
      }
      break;

      //    case BaseHeuristicAlgos::TALL:
      //
      //    case BaseHeuristicAlgos::WIDE:
      ////////////// TPL ALGOS //////////////
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
    case BaseTplAlgos::ARMPL:
      ret = Impl::BatchedArmplGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType, AViewType,
                                   BViewType, CViewType>(handle, alpha, A, B, beta, C)
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
      ret = Impl::BatchedSerialGemm<ArgTransA, ArgTransB, Algo::Gemm::Unblocked, ArgBatchSzDim, ResultsPerThread::Rank2,
                                    ScalarType, AViewType, BViewType, CViewType>(alpha, A, B, beta, C)
                .invoke();
      break;

      // case GemmKokkosBatchedAlgos::KK_SERIALSIMD:

    case GemmKokkosBatchedAlgos::KK_SERIAL_RANK0:
      ret = Impl::BatchedSerialGemm<ArgTransA, ArgTransB, Algo::Gemm::Unblocked, ArgBatchSzDim, ResultsPerThread::Rank0,
                                    ScalarType, AViewType, BViewType, CViewType>(alpha, A, B, beta, C)
                .invoke();
      break;

      //    case GemmKokkosBatchedAlgos::KK_SERIAL_SHMEM:
      //    case GemmKokkosBatchedAlgos::KK_TEAM:
      //    case GemmKokkosBatchedAlgos::KK_TEAMVECTOR:
      //    case GemmKokkosBatchedAlgos::KK_TEAMSIMD:

    case GemmKokkosBatchedAlgos::KK_DBLBUF:
      // Note: The tile sizes of 1x1x1 here will not perform well but must be
      // selected in order to function on all devices since the serial
      // execution space has a max team size of 1. KokkosKernels API users
      // will need to follow an approach similar to KK_SQUARE above for best
      // performance.

      // TODO: Add auto-selection of tile size based on inputs and device type
      ret = Impl::BatchedDblBufGemm<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType, AViewType,
                                    BViewType, CViewType, BoundsCheck::Yes, AlphaTag::No, 1, 1, 1>(handle, alpha, A, B,
                                                                                                   beta, C)
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
}  // namespace Impl
}  // namespace KokkosBatched
#endif  // __KOKKOSBATCHED_HOSTLEVEL_GEMM_IMPL_HPP__
