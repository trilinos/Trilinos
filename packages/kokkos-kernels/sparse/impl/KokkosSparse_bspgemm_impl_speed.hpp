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
          typename c_scalar_view_t, typename mpool_type>
struct KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                     b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::NumericCMEM_CPU {
  using BSPGEMM = KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                                b_lno_nnz_view_t_, b_scalar_nnz_view_t_>;

  nnz_lno_t numrows;
  nnz_lno_t numcols;
  nnz_lno_t block_dim;
  nnz_lno_t block_size;

  a_row_view_t row_mapA;
  a_nnz_view_t entriesA;
  a_scalar_view_t valuesA;

  b_row_view_t row_mapB;
  b_nnz_view_t entriesB;
  b_scalar_view_t valuesB;

  c_row_view_t rowmapC;
  c_nnz_view_t entriesC;
  c_scalar_view_t valuesC;
  mpool_type memory_space;

  nnz_lno_t *pEntriesC;
  scalar_t *pVals;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;
  const nnz_lno_t team_work_size;

  NumericCMEM_CPU(nnz_lno_t m_, nnz_lno_t k_, nnz_lno_t block_dim_, a_row_view_t row_mapA_, a_nnz_view_t entriesA_,
                  a_scalar_view_t valuesA_,

                  b_row_view_t row_mapB_, b_nnz_view_t entriesB_, b_scalar_view_t valuesB_,

                  c_row_view_t rowmapC_, c_nnz_view_t entriesC_, c_scalar_view_t valuesC_, mpool_type memory_space_,
                  const KokkosKernels::Impl::ExecSpaceType my_exec_space_, nnz_lno_t team_row_chunk_size)
      : numrows(m_),
        numcols(k_),
        block_dim(block_dim_),
        block_size(block_dim_ * block_dim_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valuesA(valuesA_),

        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valuesB(valuesB_),

        rowmapC(rowmapC_),
        entriesC(entriesC_),
        valuesC(valuesC_),
        memory_space(memory_space_),
        pEntriesC(entriesC_.data()),
        pVals(valuesC.data()),
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
#if defined(KOKKOS_ENABLE_CUDA)
      case KokkosKernels::Impl::Exec_CUDA: return row_index;
#endif
#if defined(KOKKOS_ENABLE_HIP)
      case KokkosKernels::Impl::Exec_HIP: return row_index;
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
    char *marker = (char *)(dense_accum + numcols * block_size);

    // Performs C[row_index,b_col_ind] += A[row_index,rowB] * B[rowB,b_col_ind]
    // using dense_accum[col] to accumulate scalar values,
    // marker[col] for boolean flags denoting initialized accumulators
    // and col=pEntriesC[i] to index sparse column indices.
    // Note: each CPU thread works on its own row, thus no need for locking.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type c_row_begin = rowmapC[row_index];
                           nnz_lno_t *myentries        = pEntriesC + c_row_begin;
                           scalar_t *myvals            = pVals + c_row_begin * block_size;

                           nnz_lno_t current_col_index = 0;
                           const size_type col_begin   = row_mapA[row_index];
                           const nnz_lno_t row_size    = row_mapA[row_index + 1] - col_begin;

                           for (nnz_lno_t colind = 0; colind < row_size; ++colind) {
                             size_type a_col       = colind + col_begin;
                             nnz_lno_t rowB        = entriesA[a_col];
                             const scalar_t *a_val = &valuesA[a_col * block_size];

                             size_type rowBegin  = row_mapB(rowB);
                             nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;
                             for (int i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t b_col_ind    = entriesB[adjind];
                               const scalar_t *b_val  = &valuesB[adjind * block_size];
                               if (marker[b_col_ind] == 0) {
                                 marker[b_col_ind]              = 1;
                                 myentries[current_col_index++] = b_col_ind;
                               }
                               kk_block_add_mul(block_dim, dense_accum + b_col_ind * block_size, a_val, b_val);
                             }
                           }
                           for (nnz_lno_t i = 0; i < current_col_index; ++i) {
                             nnz_lno_t ind = myentries[i];
                             scalar_t *acc = dense_accum + ind * block_size;
                             kk_block_set(block_dim, myvals + i * block_size, acc);
                             kk_block_init(block_dim, acc);
                             marker[ind] = 0;
                           }
                         });
    memory_space.release_chunk(dense_accum);
  }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__, typename b_row_view_t__,
          typename b_nnz_view_t__, typename b_scalar_view_t__, typename c_row_view_t__, typename c_nnz_view_t__,
          typename c_scalar_view_t__, typename c_nnz_tmp_view_t>

struct KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                     b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::NumericCMEM {
  static constexpr auto scalarAlignPad =
      KokkosBSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::scalarAlignPad;

  nnz_lno_t numrows;
  nnz_lno_t block_dim;
  nnz_lno_t block_size;

  a_row_view_t__ row_mapA;
  a_nnz_view_t__ entriesA;
  a_scalar_view_t__ valuesA;

  b_row_view_t__ row_mapB;
  b_nnz_view_t__ entriesB;
  b_scalar_view_t__ valuesB;

  c_row_view_t__ rowmapC;
  c_nnz_view_t__ entriesC;
  c_scalar_view_t__ valuesC;

  c_nnz_tmp_view_t beginsC;
  c_nnz_tmp_view_t nextsC;

  nnz_lno_t *pbeginsC, *pnextsC, *pEntriesC;
  scalar_t *pvaluesC;

  const size_t shared_memory_size;
  const int vector_size;
  const nnz_lno_t team_work_size;

  const int unit_memory;  // begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size;

  NumericCMEM(nnz_lno_t m_, nnz_lno_t block_dim_, a_row_view_t__ row_mapA_, a_nnz_view_t__ entriesA_,
              a_scalar_view_t__ valuesA_,

              b_row_view_t__ row_mapB_, b_nnz_view_t__ entriesB_, b_scalar_view_t__ valuesB_,

              c_row_view_t__ rowmapC_, c_nnz_view_t__ entriesC_, c_scalar_view_t__ valuesC_,

              c_nnz_tmp_view_t beginsC_, c_nnz_tmp_view_t nextsC_,

              const size_type sharedMemorySize_, const int suggested_vector_size, const nnz_lno_t team_row_chunk_size,
              int suggested_team_size_, bool KOKKOSKERNELS_VERBOSE_)
      : numrows(m_),
        block_dim(block_dim_),
        block_size(block_dim_ * block_dim_),

        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valuesA(valuesA_),

        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valuesB(valuesB_),

        rowmapC(rowmapC_),
        entriesC(entriesC_),
        valuesC(valuesC_),
        beginsC(beginsC_),
        nextsC(nextsC_),
        pbeginsC(beginsC_.data()),
        pnextsC(nextsC_.data()),
        pEntriesC(entriesC_.data()),
        pvaluesC(valuesC_.data()),
        shared_memory_size(sharedMemorySize_),

        vector_size(suggested_vector_size),
        team_work_size(team_row_chunk_size),

        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof(scalar_t) * block_size),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size / 8 / suggested_team_size_) * 8),
        shmem_key_size(),
        shared_memory_hash_func(),
        shmem_hash_size(1) {
    shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 2 - scalarAlignPad) / unit_memory);
    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tNumericCMEM -- thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " initial key size:" << shmem_key_size << std::endl;
    }
    while (shmem_hash_size * 2 <= shmem_key_size) {
      shmem_hash_size = shmem_hash_size * 2;
    }
    shared_memory_hash_func = shmem_hash_size - 1;
    shmem_key_size =
        shmem_key_size + ((shmem_key_size - shmem_hash_size) * sizeof(nnz_lno_t)) / (unit_memory - sizeof(nnz_lno_t));
    shmem_key_size = (shmem_key_size >> 1) << 1;

    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << shmem_hash_size << " shmem_key_size:" << shmem_key_size
                << std::endl;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    // get the beginning and end rows of the team.
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    // used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    // poins to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;

    // holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;
    scalar_t *vals = KokkosKernels::Impl::alignPtrTo<scalar_t>(all_shared_memory);

    KokkosKernels::Experimental::BlockHashmapAccumulator<nnz_lno_t, nnz_lno_t, scalar_t,
                                                         KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm(block_dim, shmem_key_size, shared_memory_hash_func, begins, nexts, keys, vals);

    // issue-508, TODO: understand and re-work below parallel_for loop.
    // Inialize hm2 with correct max_value_size and hashOpRHS
    // global_memory_hash_size is computed, per team of threads -- this is
    // hashOpRHS.

    KokkosKernels::Experimental::BlockHashmapAccumulator<nnz_lno_t, nnz_lno_t, scalar_t,
                                                         KokkosKernels::Experimental::HashOpType::modulo>
        hm2(block_dim, 0, 0, NULL, NULL, NULL, NULL);
    /*
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(global_memory_hash_size, global_memory_hash_size,
        pbeginsC + c_row_begin, pnextsC + c_row_begin, pEntriesC + c_row_begin,
    pvaluesC + c_row_begin);
        */

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          const size_type c_row_begin             = rowmapC[row_index];
          const nnz_lno_t global_memory_hash_size = nnz_lno_t(rowmapC[row_index + 1] - c_row_begin);

          hm2.keys        = pEntriesC + c_row_begin;
          hm2.values      = pvaluesC + c_row_begin * block_size;
          hm2.hash_begins = pbeginsC + c_row_begin;
          hm2.hash_nexts  = pnextsC + c_row_begin;

          // initialize begins.
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, shmem_hash_size), [&](int i) { begins[i] = -1; });

          // initialize hash usage sizes
          Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
            used_hash_sizes[0] = 0;
            used_hash_sizes[1] = 0;
          });

          const size_type col_begin = row_mapA[row_index];
          const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

          for (nnz_lno_t colind = 0; colind < left_work; ++colind) {
            size_type a_col       = colind + col_begin;
            nnz_lno_t rowB        = entriesA[a_col];
            const scalar_t *a_val = &valuesA[a_col * block_size];

            size_type rowBegin   = row_mapB(rowB);
            nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

            while (left_work_) {
              nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work_);
              nnz_lno_t b_col_ind      = -1;
              const scalar_t *b_val    = nullptr;
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, work_to_handle), [&](nnz_lno_t i) {
                const size_type adjind = i + rowBegin;
                b_col_ind              = entriesB[adjind];
                b_val                  = &valuesB[adjind * block_size];
              });

              int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(b_col_ind, a_val, b_val, used_hash_sizes);

              int overall_num_unsuccess = 0;

              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(teamMember, vector_size),
                  [&](const int /* threadid */, int &overall_num_unsuccess_) {
                    overall_num_unsuccess_ += num_unsuccess;
                  },
                  overall_num_unsuccess);

              if (overall_num_unsuccess) {
                nnz_lno_t hash_ = -1;
                if (num_unsuccess) {
                  hash_ = b_col_ind % global_memory_hash_size;
                }

                // int insertion =
                hm2.vector_atomic_insert_into_hash_mergeAdd_with_team_level_list_length(
                    teamMember, vector_size, hash_, b_col_ind, a_val, b_val, used_hash_sizes + 1,
                    global_memory_hash_size);
              }
              left_work_ -= work_to_handle;
              rowBegin += work_to_handle;
            }
          }

          Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
            if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
          });

          size_type num_elements = used_hash_sizes[0];

          size_type written_index = used_hash_sizes[1];
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, num_elements), [&](size_type i) {
            const auto idx = c_row_begin + written_index + i;
            pEntriesC[idx] = keys[i];
            kk_block_set(block_dim, pvaluesC + idx * block_size, vals + i * block_size);
          });
        });
  }

  size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
};

//
// * Notes on KokkosBSPGEMM_numeric_speed *
//
// Prior to this routine, KokkosBSPGEMM_numeric(...) was called
//
//   KokkosBSPGEMM_numeric(...) :
//     if (this->spgemm_algorithm == SPGEMM_KK || SPGEMM_KK_LP ==
//     this->spgemm_algorithm) :
//       call KokkosBSPGEMM_numeric_speed(...)
//     else:
//       call  KokkosBSPGEMM_numeric_hash(...)
//
//
// KokkosBSPGEMM_numeric_speed:
//
// Algorithm selection as follows and matching to kernel Tag:
//
//  Policy typedefs with tags found in: KokkosSparse_spgemm_impl.hpp
//
//  if GPU:
//    "KokkosSparse::NumericCMEM::KKSPEED::GPU" : gpu_team_policy_t,  i.e.
//    GPUTag
//
//  else :
//    "KokkosSparse::NumericCMEM_CPU::DENSE::DYNAMIC" :
//    dynamic_multicore_team_policy_t,  i.e. MultiCoreTag
//    "KokkosSparse::NumericCMEM_CPU::DENSE::STATIC" :  multicore_team_policy_t,
//    i.e. MultiCoreTag
//

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosBSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::KokkosBSPGEMM_numeric_speed(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_,
                                                       c_scalar_nnz_view_t valuesC_,
                                                       KokkosKernels::Impl::ExecSpaceType my_exec_space_) {
  if (Base::KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tSPEED MODE" << std::endl;
  }

  nnz_lno_t brows = this->row_mapB.extent(0) - 1;
  size_type bnnz  = this->valsB.extent(0);

  // get suggested vector size, teamsize and row chunk size.
  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size =
      this->handle->get_team_work_size(suggested_team_size, this->concurrency, Base::a_row_cnt);

  Kokkos::Timer numeric_speed_timer_with_free;

  if (KokkosKernels::Impl::kk_is_gpu_exec_space<typename HandleType::HandleExecSpace>()) {
    // allocate memory for begins and next to be used by the hashmap
    nnz_lno_temp_work_view_t beginsC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C keys"), valuesC_.extent(0));
    nnz_lno_temp_work_view_t nextsC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C nexts"), valuesC_.extent(0));
    Kokkos::deep_copy(beginsC, -1);

    // create the functor.
    NumericCMEM<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t, const_b_lno_row_view_t,
                const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
                nnz_lno_temp_work_view_t>
        sc(Base::a_row_cnt, block_dim, this->row_mapA, this->entriesA, this->valsA, this->row_mapB, this->entriesB,
           this->valsB,

           rowmapC_, entriesC_, valuesC_,

           beginsC, nextsC, this->shmem_size, suggested_vector_size, team_row_chunk_size, suggested_team_size,
           Base::KOKKOSKERNELS_VERBOSE);

    Kokkos::Timer timer1;
    MyExecSpace().fence();

    if (Base::KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tGPU vector_size:" << suggested_vector_size << " team_size:" << suggested_team_size
                << " chunk_size:" << team_row_chunk_size << std::endl;
    }

    timer1.reset();
    // this is basically kkmem without memory pools.
    // only executed for to check the effect of memory pools.
    Kokkos::parallel_for(
        "KokkosSparse::NumericCMEM::KKSPEED::GPU",
        gpu_team_policy_t(Base::a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), sc);
    MyExecSpace().fence();

    if (Base::KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
    }
  } else {
    Kokkos::Timer numeric_speed_timer;
    typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, scalar_t> pool_memory_space;

    KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
    int num_chunks                             = this->concurrency;

    Kokkos::Timer timer1;
    const size_t chunk_size = this->b_col_cnt * block_dim * block_dim + this->b_col_cnt / sizeof(scalar_t) + 1;
    pool_memory_space m_space(num_chunks, chunk_size, 0, my_pool_type);
    MyExecSpace().fence();

    if (Base::KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tPool Alloc Time:" << timer1.seconds() << std::endl;
      std::cout << "\tPool Size(MB):" << sizeof(scalar_t) * (num_chunks * chunk_size) / 1024. / 1024. << std::endl;
    }

    NumericCMEM_CPU<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t, const_b_lno_row_view_t,
                    const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, c_row_view_t, c_lno_nnz_view_t,
                    c_scalar_nnz_view_t, pool_memory_space>
        sc(Base::a_row_cnt, this->b_col_cnt, block_dim, this->row_mapA, this->entriesA, this->valsA, this->row_mapB,
           this->entriesB, this->valsB,

           rowmapC_, entriesC_, valuesC_, m_space, my_exec_space_, team_row_chunk_size);

    MyExecSpace().fence();
    if (Base::KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tCPU vector_size:" << suggested_vector_size << " team_size:" << suggested_team_size
                << " chunk_size:" << team_row_chunk_size << std::endl;
    }
    timer1.reset();

    if (this->use_dynamic_schedule) {
      Kokkos::parallel_for("KokkosSparse::NumericCMEM_CPU::DENSE::DYNAMIC",
                           dynamic_multicore_team_policy_t(Base::a_row_cnt / team_row_chunk_size + 1,
                                                           suggested_team_size, suggested_vector_size),
                           sc);
    } else {
      Kokkos::parallel_for("KokkosSparse::NumericCMEM_CPU::DENSE::STATIC",
                           multicore_team_policy_t(Base::a_row_cnt / team_row_chunk_size + 1, suggested_team_size,
                                                   suggested_vector_size),
                           sc);
    }

    MyExecSpace().fence();

    if (Base::KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
      std::cout << "\t\tNumeric SPEED TIME:" << numeric_speed_timer.seconds() << std::endl;
    }
  }
  if (Base::KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tNumeric SPEED TIME WITH FREE:" << numeric_speed_timer_with_free.seconds() << std::endl;
  }
}
}  // namespace Impl
}  // namespace KokkosSparse
