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
#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_DBLBUF_IMPL_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_DBLBUF_IMPL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosBatched {
/********************* BEGIN functor-level routines *********************/
/********************* END   functor-level routines *********************/

namespace Impl {
/********************* BEGIN non-functor-level routines *********************/
///
/// Implemented:
/// NT/NT, T/NT, NT/T, T/T
///
/// Not yet immplemented (ConjTranspose):
/// CT/NT, NT/CT, CT/CT
///

struct LayoutLeftTag {};
struct LayoutRightTag {};
template <class>
struct TagFromLayoutHelper;
template <>
struct TagFromLayoutHelper<Kokkos::LayoutLeft> {
  using tag = LayoutLeftTag;
};
template <>
struct TagFromLayoutHelper<Kokkos::LayoutRight> {
  using tag = LayoutRightTag;
};
template <class Layout>
using TagFromLayout = typename TagFromLayoutHelper<Layout>::tag;

// TODO - scaling between (32x32, 64x64)
//   Option 0: Increase number of tiles and figure out how to map kokkos teams
//             into cuda grid. Keep team size and vector lanes constant.
//             TODO: write up small example and ask Christian. [DONE,
//             MdRangePolicy not applicable here]
//   Option 1: Increase register sizes to handle rows/cols past tile size
//   Option 2: Fix league_size and have single team solve full tile followed
//   by same team solving extra rows/cols (without multiplying by the
//   zero rows/cols)

// clang-format off
/// \brief Non-blocking general matrix multiply on a batch of
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
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim, class HandleType, class ScalarType, class AViewType,
          class BViewType, class CViewType, class ArgBoundsCheck, class ArgAlphaFmaTag, int TILE_M, int TILE_N,
          int TILE_K>
class BatchedDblBufGemm {
 private:
  using AlphaMulTag =
      std::conditional_t<std::is_same<ArgAlphaFmaTag, AlphaTag::Yes>::value, AlphaTag::No, AlphaTag::Yes>;

  HandleType *const __handle;
  AViewType __A;
  BViewType __B;
  CViewType __C;
  ScalarType __alpha, __beta;
  ArgTransA __transA_tag;
  ArgTransB __transB_tag;
  ArgBatchSzDim __batch_layout_tag;
  ArgBoundsCheck __bounds_check_tag;
  ArgAlphaFmaTag __alpha_fma_tag;
  AlphaMulTag __alpha_mul_tag;
  int __c_batch_size, __c_m, __c_n;

  using view_value_type      = typename CViewType::value_type;
  using layout_type          = typename CViewType::array_layout;
  using device_type          = typename CViewType::device_type;
  using execution_space_type = typename device_type::execution_space;
  using scratch_space_type   = typename execution_space_type::scratch_memory_space;
  using view_type_2d_scratch = Kokkos::View<view_value_type **, Kokkos::LayoutRight, scratch_space_type>;

 public:
  BatchedDblBufGemm(HandleType *const handle, ScalarType alpha, AViewType A, BViewType B, ScalarType beta, CViewType C)
      : __handle(handle), __A(A), __B(B), __C(C), __alpha(alpha), __beta(beta) {}

  int invoke() {
    __run();
    return 0;
  }

 private:
  void __run() {
    using policy_type = Kokkos::TeamPolicy<TagFromLayout<layout_type>, execution_space_type>;
    using member_type = typename policy_type::member_type;

    // Compile-time expressions required for functor-level register allocations:
    //   Each team uses a shmem buffer and statically allocated register buffer.
    //   Below, we need a 1-1 mapping between GPU threads and register
    //   allocations to ensure that each GPU thread does not step on another
    //   GPU threads' registers. In short, we must map register allocations
    //   to parallel_for loop bounds.
    // TODO: check these expressions for all tile_m, tile_n, tile_k in Z+.
    constexpr int reg_m    = TILE_M / TILE_K;
    constexpr int reg_n    = TILE_N / TILE_K + 2 * !!(TILE_N % TILE_K);
    constexpr int stride_m = TILE_K;
    constexpr int stride_n = TILE_N / reg_n;
    using functor_type     = Functor<member_type, reg_m, reg_n, stride_m, stride_n>;

    functor_type functor(*this, __A, __B, __C);

    if (__handle->enableDebug) {
      std::cout << "algo_type:" << __handle->get_kernel_algo_type() << std::endl
                << "__c_m:" << __c_m << std::endl
                << "__c_n:" << __c_n << std::endl
                << "reg_m:" << reg_m << std::endl
                << "reg_n:" << reg_n << std::endl
                << "stride_m:" << stride_m << std::endl
                << "stride_n:" << stride_n << std::endl;
    }

    // Each team solves a single tile. Within each tile, the team solves
    // all __n_tile_k_tiles one at a time.
    size_t league_size = __c_batch_size * functor.get_n_sub_tiles();
    int team_size      = stride_m;
    int vector_len     = stride_n;

    const int max_team_size =
        policy_type(league_size, Kokkos::AUTO, vector_len).team_size_max(functor, Kokkos::ParallelForTag());
    if (team_size > max_team_size) {
      std::ostringstream os;
      os << "KokkosBatched::BatchedGemm with kernelAlgoType = " << std::to_string(__handle->get_kernel_algo_type())
         << " does not support team_size > " << std::to_string(max_team_size) << "." << std::endl
         << " The tile dimensions must be adjusted." << std::endl;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }

    const int max_vector_len = policy_type(league_size, team_size, Kokkos::AUTO).vector_length_max();
    if (vector_len > max_vector_len) {
      std::ostringstream os;
      os << "KokkosBatched::BatchedGemm with kernelAlgoType = " << std::to_string(__handle->get_kernel_algo_type())
         << " does not support vector_len > " << std::to_string(max_vector_len) << "." << std::endl
         << " The tile dimensions must be adjusted." << std::endl;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }

    if (__handle->enableDebug) {
      std::cout << "max_team_size:" << max_team_size << " team_size:" << team_size << std::endl
                << "max_vector_len:" << max_vector_len << " vector_len:" << vector_len << std::endl
                << "TILE_M:" << TILE_M << std::endl
                << "TILE_N:" << TILE_N << std::endl
                << "TILE_K:" << TILE_K << std::endl;
    }

    // TODO: Use statically allocated shmem
    int shmem_size =
        view_type_2d_scratch::shmem_size(TILE_M, TILE_K) + view_type_2d_scratch::shmem_size(TILE_K, TILE_N);

    // Each member solves a portion of TILE_K in parallel with other members
    policy_type team_policy(league_size, team_size, vector_len);
    team_policy.set_scratch_size(0, Kokkos::PerTeam(shmem_size));

    Kokkos::parallel_for("BatchedDblBufGemm", team_policy, functor);
  }

 public:
  // Make Functor public for cuda 9.
  // See https://github.com/kokkos/kokkos-kernels/issues/1121.
  template <class MemberType, int REG_M, int REG_N, int STRIDE_M, int STRIDE_N>
  class Functor {
   private:
    BatchedDblBufGemm &__ei;
    AViewType __A;
    BViewType __B;
    CViewType __C;
    ScalarType __alpha, __beta;
    int __k;
    size_t __n_sub_tiles;
    unsigned __tiles_per_col, __tiles_per_row;

   public:
    size_t get_n_sub_tiles() { return __n_sub_tiles; }

    // NOTE: We cannot use __ei.{__A,__B,__C,__beta,__alpha,__k} in the operator
    // below. If those are used, we  get an invalid memory error from cuda. I
    // suspect this is due the values not being copied to device and then
    // runtime resolution of the host address &__ei.
    Functor(BatchedDblBufGemm &ei, AViewType A, BViewType B, CViewType C) : __ei(ei), __A(A), __B(B), __C(C) {
      if (std::is_same<ArgBatchSzDim, BatchLayout::Left>::value) {
        ei.__c_batch_size = ei.__C.extent_int(0);
        ei.__c_m          = ei.__C.extent_int(1);
        ei.__c_n          = ei.__C.extent_int(2);
        if (std::is_same<ArgTransA, Trans::Transpose>::value)
          __k = ei.__A.extent_int(1);
        else
          __k = ei.__A.extent_int(2);
      } else {
        ei.__c_batch_size = ei.__C.extent_int(2);
        ei.__c_m          = ei.__C.extent_int(0);
        ei.__c_n          = ei.__C.extent_int(1);
        if (std::is_same<ArgTransA, Trans::Transpose>::value)
          __k = ei.__A.extent_int(0);
        else
          __k = ei.__A.extent_int(1);
      }
      __beta  = ei.__beta;   // Copy to device
      __alpha = ei.__alpha;  // Copy to device
      // To handle truncation of tiles per row/col, round up to one extra tile
      // with '!!'. This extra tile will hang off the edge of the 2-rank matrix.
      // For cases where tiles hang off the edge, we over-compute 0s within
      // registers via a conditional bounds check selected at compile-time.
      __tiles_per_row = ei.__c_m / TILE_M + !!((unsigned)ei.__c_m % TILE_M);
      __tiles_per_col = ei.__c_n / TILE_N + !!((unsigned)ei.__c_n % TILE_N);

      __n_sub_tiles = __tiles_per_row * __tiles_per_col;
    }

    KOKKOS_INLINE_FUNCTION
    void __mul(view_value_type a, view_value_type b, view_value_type &c, const AlphaTag::No &) const { c += a * b; }

    KOKKOS_INLINE_FUNCTION
    void __mul(view_value_type a, view_value_type b, view_value_type &c, const AlphaTag::Yes &) const {
      c += a * b * __alpha;
    }

    KOKKOS_INLINE_FUNCTION
    void __rshmem_and_mul(const int &thread_id, const int &vlane_id, const unsigned &nk, view_value_type reg_a[REG_M],
                          view_value_type reg_b[REG_N], view_value_type reg_c[REG_M][REG_N],
                          view_type_2d_scratch &svA_scr, view_type_2d_scratch &svB_scr) const {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
      for (unsigned k = 0; k < nk; ++k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int m = 0; m < REG_M; ++m) reg_a[m] = svA_scr(thread_id + m * STRIDE_M, k);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int n = 0; n < REG_N; ++n) reg_b[n] = svB_scr(k, vlane_id + n * STRIDE_N);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int m = 0; m < REG_M; ++m) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int n = 0; n < REG_N; ++n) __mul(reg_a[m], reg_b[n], reg_c[m][n], __ei.__alpha_mul_tag);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void __rshmem_and_mul_ll(const int &thread_id, const int &vlane_id, const unsigned &nk,
                             view_value_type reg_a[REG_M], view_value_type reg_b[REG_N],
                             view_value_type reg_c[REG_M][REG_N], view_type_2d_scratch &svA_scr,
                             view_type_2d_scratch &svB_scr) const {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
      for (unsigned k = 0; k < nk; ++k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int m = 0; m < REG_M; ++m) reg_a[m] = svA_scr(k, vlane_id + m * STRIDE_M);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int n = 0; n < REG_N; ++n) reg_b[n] = svB_scr(thread_id + n * STRIDE_N, k);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
        for (int m = 0; m < REG_M; ++m) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int n = 0; n < REG_N; ++n) __mul(reg_a[m], reg_b[n], reg_c[m][n], __ei.__alpha_mul_tag);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(LayoutRightTag, const MemberType &member) const {
      // TODO: use Kokkos view with compile-time size to allocating register??
      //  Then we can use local deep copy for prefetch_reg population.
      // Allocate registers used for prefetching
      view_value_type prefetch_reg_a[REG_M] = {0}, prefetch_reg_b[REG_N] = {0};

      // Allocate registers used for FMAs
      view_value_type reg_a[REG_M] = {0}, reg_b[REG_N] = {0}, reg_c[REG_M][REG_N] = {{0}};
      // TODO: look at local loads and stores via nvprof
      // TODO: look at GPU trace in nvprof to find out how many registers are
      // used.

      unsigned batch_idx = member.league_rank() / __n_sub_tiles;

      // Compute starting tile offsets for each team into svA, svB, svC
      unsigned local_team_idx = member.league_rank() % __n_sub_tiles;
      unsigned start_m        = (local_team_idx / __tiles_per_col) * TILE_M;
      unsigned start_n        = (local_team_idx % __tiles_per_col) * TILE_N;

      int kk;

      // Fetch entire 2-rank sub-matrix
      auto svA =
          subview_wrapper(__A, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag, __ei.__transA_tag);
      auto svB =
          subview_wrapper(__B, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag, __ei.__transB_tag);
      auto svC = subview_wrapper(__C, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag);

      // Allocate scratch memory buffers used for prefetching
      view_type_2d_scratch svA_scr(member.team_scratch(0), TILE_M, TILE_K);
      view_type_2d_scratch svB_scr(member.team_scratch(0), TILE_K, TILE_N);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, STRIDE_M), [&](const int &thread_id) {
        int m_offset = thread_id + start_m;

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, 0, STRIDE_N), [&](const int &vlane_id) {
          int n_offset = vlane_id + start_n;

          // Here we populate scratch memory with one or more "k" tiles for
          // every thread of the team!
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int i = 0; i < REG_N * STRIDE_N; i += STRIDE_N)
            svB_scr(thread_id, vlane_id + i) =
                access_view_bounds_check<view_value_type>(svB, thread_id, n_offset + i, __ei.__bounds_check_tag);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int i = 0; i < REG_M * STRIDE_M; i += STRIDE_M)
            svA_scr(thread_id + i, vlane_id) =
                access_view_bounds_check<view_value_type>(svA, m_offset + i, vlane_id, __ei.__bounds_check_tag);

          // Wait for A, B to reside in scratch memory
          member.team_barrier();

          // Each thread calculates a single dot product in chunks of
          // size TILE_K
          for (kk = 0; kk < __k - TILE_K; kk += TILE_K) {
            int k_tile_offset = kk + TILE_K;

            // Get this threads next TILE_K entries from global memory
            // Each thread has its own copy of prefetch_reg_b.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_N; ++i)
              prefetch_reg_b[i] = access_view_bounds_check<view_value_type>(
                  svB, thread_id + k_tile_offset, n_offset + i * STRIDE_N, __ei.__bounds_check_tag);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_M; ++i)
              prefetch_reg_a[i] = access_view_bounds_check<view_value_type>(
                  svA, m_offset + i * STRIDE_M, vlane_id + k_tile_offset, __ei.__bounds_check_tag);

            __rshmem_and_mul(thread_id, vlane_id, TILE_K, reg_a, reg_b, reg_c, svA_scr, svB_scr);

            // Wait for:
            //   1. prefetch_regs to be populated
            //   2. for shmem to no longer be read from
            member.team_barrier();

            // populate shmem from prefetch registers. Each thread has its own
            // copy of prefetch_reg_b.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_N; ++i) svB_scr(thread_id, vlane_id + i * STRIDE_N) = prefetch_reg_b[i];

              // populate shmem from prefetch registers. Each thread has its own
              // copy of prefetch_reg_a.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_M; ++i) svA_scr(thread_id + i * STRIDE_M, vlane_id) = prefetch_reg_a[i];

            // Wait for shmem stores to land before performing next
            // TILE_K multiply
            member.team_barrier();
          }  // end n_tile_k_tiles loop

          // Multiply last tile, may be a partial tile
          __rshmem_and_mul(thread_id, vlane_id, __k - kk, reg_a, reg_b, reg_c, svA_scr, svB_scr);

          // store results back to global memory
          if (__beta == 0.0F) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int m = 0; m < REG_M; ++m) {
              int cm = m_offset + m * STRIDE_M;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
              for (int n = 0; n < REG_N; ++n) {
                int cn = n_offset + n * STRIDE_N;
                fma_bounds_check(svC, cm, cn, reg_c[m][n], __alpha, __ei.__alpha_fma_tag, __ei.__bounds_check_tag);
              }
            }
          } else {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int m = 0; m < REG_M; ++m) {
              int cm = m_offset + m * STRIDE_M;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
              for (int n = 0; n < REG_N; ++n) {
                int cn = n_offset + n * STRIDE_N;
                fma_bounds_check(svC, cm, cn, reg_c[m][n], __alpha, __beta, __ei.__alpha_fma_tag,
                                 __ei.__bounds_check_tag);
              }
            }
          }
        });
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(LayoutLeftTag, const MemberType &member) const {
      // TODO: use Kokkos view with compile-time size to allocating register??
      //  Then we can use local deep copy for prefetch_reg population.
      // Allocate registers used for prefetching
      view_value_type prefetch_reg_a[REG_M] = {0}, prefetch_reg_b[REG_N] = {0};

      // Allocate registers used for FMAs
      view_value_type reg_a[REG_M] = {0}, reg_b[REG_N] = {0}, reg_c[REG_M][REG_N] = {{0}};
      // TODO: look at local loads and stores via nvprof
      // TODO: look at GPU trace in nvprof to find out how many registers are
      // used.

      unsigned batch_idx = member.league_rank() / __n_sub_tiles;

      // Compute starting tile offsets for each team into svA, svB, svC
      unsigned local_team_idx = member.league_rank() % __n_sub_tiles;
      unsigned start_m        = (local_team_idx % __tiles_per_row) * TILE_M;
      unsigned start_n        = (local_team_idx / __tiles_per_row) * TILE_N;

      int kk;

      // Fetch entire 2-rank sub-matrix
      auto svA =
          subview_wrapper(__A, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag, __ei.__transA_tag);
      auto svB =
          subview_wrapper(__B, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag, __ei.__transB_tag);
      auto svC = subview_wrapper(__C, batch_idx, Kokkos::ALL(), Kokkos::ALL(), __ei.__batch_layout_tag);

      // Allocate scratch memory buffers used for prefetching
      view_type_2d_scratch svA_scr(member.team_scratch(0), TILE_K, TILE_M);
      view_type_2d_scratch svB_scr(member.team_scratch(0), TILE_N, TILE_K);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, STRIDE_N), [&](const int &thread_id) {
        int n_offset = thread_id + start_n;

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, 0, STRIDE_M), [&](const int &vlane_id) {
          int m_offset = vlane_id + start_m;

          // Here we populate scratch memory with one or more "k" tiles for
          // every thread of the team!
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int i = 0; i < REG_N * STRIDE_N; i += STRIDE_N)
            svB_scr(thread_id + i, vlane_id) =
                access_view_bounds_check<view_value_type>(svB, vlane_id, n_offset + i, __ei.__bounds_check_tag);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
          for (int i = 0; i < REG_M * STRIDE_M; i += STRIDE_M)
            svA_scr(thread_id, vlane_id + i) =
                access_view_bounds_check<view_value_type>(svA, m_offset + i, thread_id, __ei.__bounds_check_tag);

          // Wait for A, B to reside in scratch memory
          member.team_barrier();

          // Each thread calculates a single dot product in chunks of
          // size TILE_K
          for (kk = 0; kk < __k - TILE_K; kk += TILE_K) {
            int k_tile_offset = kk + TILE_K;

            // Get this threads next TILE_K entries from global memory
            // Each thread has its own copy of prefetch_reg_b.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_N; ++i)
              prefetch_reg_b[i] = access_view_bounds_check<view_value_type>(
                  svB, vlane_id + k_tile_offset, n_offset + i * STRIDE_N, __ei.__bounds_check_tag);

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_M; ++i)
              prefetch_reg_a[i] = access_view_bounds_check<view_value_type>(
                  svA, m_offset + i * STRIDE_M, thread_id + k_tile_offset, __ei.__bounds_check_tag);

            __rshmem_and_mul_ll(thread_id, vlane_id, TILE_K, reg_a, reg_b, reg_c, svA_scr, svB_scr);

            // Wait for:
            //   1. prefetch_regs to be populated
            //   2. for shmem to no longer be read from
            member.team_barrier();

            // populate shmem from prefetch registers. Each thread has its own
            // copy of prefetch_reg_b.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_N; ++i) svB_scr(thread_id + i * STRIDE_N, vlane_id) = prefetch_reg_b[i];

              // populate shmem from prefetch registers. Each thread has its own
              // copy of prefetch_reg_a.
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int i = 0; i < REG_M; ++i) svA_scr(thread_id, vlane_id + i * STRIDE_M) = prefetch_reg_a[i];

            // Wait for shmem stores to land before performing next
            // TILE_K multiply
            member.team_barrier();
          }  // end n_tile_k_tiles loop

          // Multiply last tile, may be a partial tile
          __rshmem_and_mul_ll(thread_id, vlane_id, __k - kk, reg_a, reg_b, reg_c, svA_scr, svB_scr);

          // store results back to global memory
          if (__beta == 0.0F) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int n = 0; n < REG_N; ++n) {
              int cn = n_offset + n * STRIDE_N;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
              for (int m = 0; m < REG_M; ++m) {
                int cm = m_offset + m * STRIDE_M;
                fma_bounds_check(svC, cm, cn, reg_c[m][n], __alpha, __ei.__alpha_fma_tag, __ei.__bounds_check_tag);
              }
            }
          } else {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
            for (int n = 0; n < REG_N; ++n) {
              int cn = n_offset + n * STRIDE_N;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UNROLL
              for (int m = 0; m < REG_M; ++m) {
                int cm = m_offset + m * STRIDE_M;
                fma_bounds_check(svC, cm, cn, reg_c[m][n], __alpha, __beta, __ei.__alpha_fma_tag,
                                 __ei.__bounds_check_tag);
              }
            }
          }
        });
      });
    }
  };
};
/********************* END non-functor-level routines *********************/
}  // namespace Impl

}  // namespace KokkosBatched

#endif
