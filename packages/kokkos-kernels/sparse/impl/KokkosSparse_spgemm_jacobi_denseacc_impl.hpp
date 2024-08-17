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

#include "KokkosKernels_Utils.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
          typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
          typename c_scalar_view_t, typename dinv_view_t, typename mpool_type>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::JacobiSpGEMMDenseAcc {
  nnz_lno_t numrows;
  nnz_lno_t numcols;

  a_row_view_t row_mapA;
  a_nnz_view_t entriesA;
  a_scalar_view_t valuesA;

  b_row_view_t row_mapB;
  b_nnz_view_t entriesB;
  b_scalar_view_t valuesB;

  c_row_view_t row_mapC;
  c_nnz_view_t entriesC;
  c_scalar_view_t valuesC;

  typename c_scalar_view_t::const_value_type omega;
  dinv_view_t dinv;

  nnz_lno_t *pEntriesC;
  scalar_t *pVals;

  mpool_type memory_space;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;
  const nnz_lno_t team_work_size;

  JacobiSpGEMMDenseAcc(nnz_lno_t m_, nnz_lno_t k_, a_row_view_t row_mapA_, a_nnz_view_t entriesA_,
                       a_scalar_view_t valuesA_, b_row_view_t row_mapB_, b_nnz_view_t entriesB_,
                       b_scalar_view_t valuesB_, c_row_view_t row_mapC_, c_nnz_view_t entriesC_,
                       c_scalar_view_t valuesC_, typename c_scalar_view_t::const_value_type omega_, dinv_view_t dinv_,
                       mpool_type memory_space_, const KokkosKernels::Impl::ExecSpaceType my_exec_space_,
                       nnz_lno_t team_row_chunk_size)
      : numrows(m_),
        numcols(k_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valuesA(valuesA_),
        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valuesB(valuesB_),
        row_mapC(row_mapC_),
        entriesC(entriesC_),
        valuesC(valuesC_),
        omega(omega_),
        dinv(dinv_),
        pEntriesC(entriesC_.data()),
        pVals(valuesC.data()),
        memory_space(memory_space_),
        my_exec_space(my_exec_space_),
        team_work_size(team_row_chunk_size) {}

  KOKKOS_INLINE_FUNCTION
  size_t get_thread_id(const size_t row_index) const {
    switch (my_exec_space) {
      default: return row_index;
#if defined(KOKKOS_ENABLE_SERIAL)
      case KokkosKernels::Impl::Exec_SERIAL: return 0;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
      case KokkosKernels::Impl::Exec_OMP: return Kokkos::OpenMP::impl_hardware_thread_id();
#endif
#if defined(KOKKOS_ENABLE_THREADS)
      case KokkosKernels::Impl::Exec_THREADS: return Kokkos::Threads::impl_hardware_thread_id();
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag &, const team_member_t &teamMember) const {
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    scalar_t *dense_accum = NULL;
    size_t tid            = get_thread_id(team_row_begin + teamMember.team_rank());
    while (dense_accum == NULL) {
      dense_accum = (scalar_t *)(memory_space.allocate_chunk(tid));
    }
    char *marker = (char *)(dense_accum + numcols);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type c_row_begin = row_mapC[row_index];
                           nnz_lno_t *myentries        = pEntriesC + c_row_begin;
                           scalar_t *myvals            = pVals + c_row_begin;

                           nnz_lno_t current_col_index = 0;
                           const size_type col_begin   = row_mapA[row_index];
                           const nnz_lno_t nnza        = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

                           // Insert B
                           size_type rowBegin  = row_mapB(row_index);
                           nnz_lno_t left_work = row_mapB(row_index + 1) - rowBegin;
                           for (int i = 0; i < left_work; ++i) {
                             const size_type adjind = i + rowBegin;
                             nnz_lno_t b_col_ind    = entriesB[adjind];
                             scalar_t b_val         = valuesB[adjind];
                             if (marker[b_col_ind] == 0) {
                               marker[b_col_ind]              = 1;
                               myentries[current_col_index++] = b_col_ind;
                             }
                             dense_accum[b_col_ind] += b_val;
                           }

                           // Insert -omega * dinv * A*B
                           const scalar_t mult = -omega * dinv(row_index, 0);
                           for (nnz_lno_t colind = 0; colind < nnza; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];
                             scalar_t valA   = valuesA[a_col] * mult;

                             rowBegin  = row_mapB(rowB);
                             left_work = row_mapB(rowB + 1) - rowBegin;
                             for (int i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t b_col_ind    = entriesB[adjind];
                               scalar_t b_val         = valuesB[adjind] * valA;
                               if (marker[b_col_ind] == 0) {
                                 marker[b_col_ind]              = 1;
                                 myentries[current_col_index++] = b_col_ind;
                               }
                               dense_accum[b_col_ind] += b_val;
                             }
                           }
                           for (nnz_lno_t i = 0; i < current_col_index; ++i) {
                             nnz_lno_t ind    = myentries[i];
                             myvals[i]        = dense_accum[ind];
                             dense_accum[ind] = 0;
                             marker[ind]      = 0;
                           }
                         });
    memory_space.release_chunk(dense_accum);
  }
};

// ============================================================================
// This is a Jacobi-fused SPGEMM implementation which uses a dense accumulator.
// Copied from KokkosSparse_spgemm_impl_speed.hpp, removed the Cuda code, and
// updated for the jacobi-fused functionality by adding omega and dinv. The
// Cuda code was removed because it is never called (Deveci19 paper does not
// recommend using dense accumulators on GPUs).
// ============================================================================

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t, typename dinv_view_t>
void KokkosSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::KokkosSPGEMM_jacobi_denseacc(c_row_view_t row_mapC_, c_lno_nnz_view_t entriesC_,
                                                        c_scalar_nnz_view_t valuesC_,
                                                        typename c_scalar_nnz_view_t::const_value_type omega,
                                                        dinv_view_t dinv,
                                                        KokkosKernels::Impl::ExecSpaceType my_exec_space_) {
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tDENSE ACC MODE" << std::endl;
  }

  // Initialize the timer
  Kokkos::Timer jacobi_speed_timer;

  // Get suggested vector size, teamsize and row chunk size
  nnz_lno_t brows               = row_mapB.extent(0) - 1;
  size_type bnnz                = valsB.extent(0);
  int suggested_vector_size     = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);

  // Allocate the memory pool
  KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
  int num_chunks                             = concurrency;
  typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, scalar_t> pool_memory_space;
  Kokkos::Timer timer;
  pool_memory_space m_space(num_chunks, this->b_col_cnt + (this->b_col_cnt) / sizeof(scalar_t) + 1, 0, my_pool_type);
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tPool Alloc Time:" << timer.seconds() << std::endl;
    std::cout << "\tPool Size(MB):"
              << sizeof(scalar_t) * (num_chunks * (this->b_col_cnt + (this->b_col_cnt) / sizeof(scalar_t) + 1)) /
                     1024. / 1024.
              << std::endl;
  }

  // Initialize the functor
  JacobiSpGEMMDenseAcc<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
                       const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, c_row_view_t,
                       c_lno_nnz_view_t, c_scalar_nnz_view_t, dinv_view_t, pool_memory_space>
      jacobi(a_row_cnt, b_col_cnt, row_mapA, entriesA, valsA, row_mapB, entriesB, valsB, row_mapC_, entriesC_, valuesC_,
             omega, dinv, m_space, my_exec_space_, team_row_chunk_size);
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tCPU vector_size:" << suggested_vector_size << " team_size:" << suggested_team_size
              << " chunk_size:" << team_row_chunk_size << std::endl;
  }
  timer.reset();

  if (use_dynamic_schedule) {
    Kokkos::parallel_for("KokkosSparse::Jacobi::DenseAcc::Dynamic",
                         dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size,
                                                         suggested_vector_size),
                         jacobi);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::Jacobi::DenseAcc::Static",
        multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
        jacobi);
  }

  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tJacobi COMP TIME:" << timer.seconds() << std::endl;
    std::cout << "\t\tJacobi TOTAL TIME:" << jacobi_speed_timer.seconds() << std::endl;
  }
}
}  // namespace Impl
}  // namespace KokkosSparse
