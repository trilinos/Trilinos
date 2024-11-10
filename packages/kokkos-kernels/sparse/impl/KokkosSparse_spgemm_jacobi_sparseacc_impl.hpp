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

#define HASHSCALAR 107

#include "KokkosKernels_Utils.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
          typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
          typename c_scalar_view_t, typename dinv_view_t, typename pool_memory_type>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::JacobiSpGEMMSparseAcc {
  nnz_lno_t numrows;

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
  scalar_t *pvaluesC;
  const size_t shared_memory_size;
  const int vector_size;
  pool_memory_type memory_space;

  const nnz_lno_t pow2_hash_size;
  const nnz_lno_t max_nnz;
  const nnz_lno_t pow2_hash_func;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;
  const nnz_lno_t team_work_size;

  const int unit_memory;  // begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t thread_shmem_key_size;
  nnz_lno_t thread_shmem_hash_func;
  nnz_lno_t thread_shmem_hash_size;

  nnz_lno_t team_shmem_key_size;
  nnz_lno_t team_shmem_hash_func;
  nnz_lno_t team_shmem_hash_size;
  nnz_lno_t team_cuckoo_key_size;
  nnz_lno_t team_cuckoo_hash_func;

  nnz_lno_t max_first_level_hash_size;
  row_lno_persistent_work_view_t flops_per_row;

  JacobiSpGEMMSparseAcc(nnz_lno_t m_, a_row_view_t row_mapA_, a_nnz_view_t entriesA_, a_scalar_view_t valuesA_,
                        b_row_view_t row_mapB_, b_nnz_view_t entriesB_, b_scalar_view_t valuesB_,
                        c_row_view_t row_mapC_, c_nnz_view_t entriesC_, c_scalar_view_t valuesC_,
                        typename c_scalar_view_t::const_value_type omega_, dinv_view_t dinv_,
                        size_t shared_memory_size_, int vector_size_, pool_memory_type mpool_, nnz_lno_t min_hash_size,
                        nnz_lno_t max_nnz_, int suggested_team_size_,
                        const KokkosKernels::Impl::ExecSpaceType my_exec_space_, nnz_lno_t team_row_chunk_size,
                        double first_level_cut_off, row_lno_persistent_work_view_t flops_per_row_,
                        bool KOKKOSKERNELS_VERBOSE_)
      : numrows(m_),
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
        pvaluesC(valuesC_.data()),
        shared_memory_size(shared_memory_size_),
        vector_size(vector_size_),
        memory_space(mpool_),
        pow2_hash_size(min_hash_size),
        max_nnz(max_nnz_),
        pow2_hash_func(min_hash_size - 1),
        my_exec_space(my_exec_space_),
        team_work_size(team_row_chunk_size),
        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof(scalar_t)),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size / 8 / suggested_team_size_) * 8),
        thread_shmem_key_size(),
        thread_shmem_hash_func(),
        thread_shmem_hash_size(1),
        team_shmem_key_size(),
        team_shmem_hash_func(),
        team_shmem_hash_size(1),
        team_cuckoo_key_size(1),
        team_cuckoo_hash_func(1),
        max_first_level_hash_size(1),
        flops_per_row(flops_per_row_)

  {
    nnz_lno_t tmp_team_cuckoo_key_size =
        ((shared_memory_size - sizeof(nnz_lno_t) * 2) / (sizeof(nnz_lno_t) + sizeof(scalar_t)));

    while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
    team_cuckoo_hash_func = team_cuckoo_key_size - 1;

    // How many extra bytes are needed to align a scalar_t after an array of
    // nnz_lno_t, in the worst case?
    constexpr size_t scalarAlignPad =
        (alignof(scalar_t) > alignof(nnz_lno_t)) ? (alignof(scalar_t) - alignof(nnz_lno_t)) : 0;
    team_shmem_key_size   = ((shared_memory_size - sizeof(nnz_lno_t) * 4 - scalarAlignPad) / unit_memory);
    thread_shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4 - scalarAlignPad) / unit_memory);
    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- sizeof(scalar_t): " << sizeof(scalar_t)
                << "  sizeof(nnz_lno_t): " << sizeof(nnz_lno_t) << "  suggested_team_size: " << suggested_team_size
                << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " initial key size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- team shared_memory:" << shared_memory_size
                << " unit_memory:" << unit_memory << " initial team key size:" << team_shmem_key_size << std::endl;
    }
    while (thread_shmem_hash_size * 2 <= thread_shmem_key_size) {
      thread_shmem_hash_size = thread_shmem_hash_size * 2;
    }
    while (team_shmem_hash_size * 2 <= team_shmem_key_size) {
      team_shmem_hash_size = team_shmem_hash_size * 2;
    }
    team_shmem_hash_func   = team_shmem_hash_size - 1;
    thread_shmem_hash_func = thread_shmem_hash_size - 1;
    team_shmem_key_size    = team_shmem_key_size + ((team_shmem_key_size - team_shmem_hash_size) * sizeof(nnz_lno_t)) /
                                                    (sizeof(nnz_lno_t) * 2 + sizeof(scalar_t));
    team_shmem_key_size = (team_shmem_key_size >> 1) << 1;

    thread_shmem_key_size =
        thread_shmem_key_size + ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t)) /
                                    (sizeof(nnz_lno_t) * 2 + sizeof(scalar_t));
    thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;

    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " resized key size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- team shared_memory:" << shared_memory_size
                << " unit_memory:" << unit_memory << " resized team key size:" << team_shmem_key_size << std::endl;
    }

    max_first_level_hash_size = first_level_cut_off * team_cuckoo_key_size;
    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " initial key size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- team_memory:" << shared_memory_size << " unit_memory:" << unit_memory
                << " initial team key size:" << team_shmem_key_size << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- adjusted hashsize:" << thread_shmem_hash_size
                << " thread_shmem_key_size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tJacobiSpGEMMSparseAcc -- adjusted team hashsize:" << team_shmem_hash_size
                << " team_shmem_key_size:" << team_shmem_key_size << std::endl;
      std::cout << "\t\t  team_cuckoo_key_size:" << team_cuckoo_key_size
                << " team_cuckoo_hash_func:" << team_cuckoo_hash_func
                << " max_first_level_hash_size:" << max_first_level_hash_size << std::endl;
      std::cout << "\t\t  pow2_hash_size:" << pow2_hash_size << " pow2_hash_func:" << pow2_hash_func << std::endl;
    }
  }

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

  // linear probing with tracking.
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag4 &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(memory_space.allocate_chunk(tid));
    }

    nnz_lno_t *used_indices = (nnz_lno_t *)(tmp);
    tmp += max_nnz;
    nnz_lno_t *hash_ids = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;
    scalar_t *hash_values = KokkosKernels::Impl::alignPtrTo<scalar_t>(tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t used_count = 0;

                           // Insert B
                           size_type rowBegin   = row_mapB(row_index);
                           nnz_lno_t left_workB = row_mapB(row_index + 1) - rowBegin;

                           for (nnz_lno_t i = 0; i < left_workB; ++i) {
                             const size_type adjind = i + rowBegin;
                             nnz_lno_t b_col_ind    = entriesB[adjind];
                             scalar_t b_val         = valuesB[adjind];
                             nnz_lno_t hash         = (b_col_ind * HASHSCALAR) & pow2_hash_func;

                             while (true) {
                               if (hash_ids[hash] == -1) {
                                 used_indices[used_count++] = hash;
                                 hash_ids[hash]             = b_col_ind;
                                 hash_values[hash]          = b_val;
                                 break;
                               } else if (hash_ids[hash] == b_col_ind) {
                                 hash_values[hash] += b_val;
                                 break;
                               } else {
                                 hash = (hash + 1) & pow2_hash_func;
                               }
                             }
                           }

                           // Insert -omega * dinv * A*B
                           scalar_t mult             = -omega * dinv(row_index, 0);
                           const size_type col_begin = row_mapA[row_index];
                           const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
                           for (nnz_lno_t ii = 0; ii < left_work; ++ii) {
                             size_type a_col = col_begin + ii;
                             nnz_lno_t rowB  = entriesA[a_col];
                             scalar_t valA   = valuesA[a_col] * mult;

                             rowBegin   = row_mapB(rowB);
                             left_workB = row_mapB(rowB + 1) - rowBegin;

                             for (nnz_lno_t i = 0; i < left_workB; ++i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t b_col_ind    = entriesB[adjind];
                               scalar_t b_val         = valuesB[adjind] * valA;
                               nnz_lno_t hash         = (b_col_ind * HASHSCALAR) & pow2_hash_func;

                               while (true) {
                                 if (hash_ids[hash] == -1) {
                                   used_indices[used_count++] = hash;
                                   hash_ids[hash]             = b_col_ind;
                                   hash_values[hash]          = b_val;
                                   break;
                                 } else if (hash_ids[hash] == b_col_ind) {
                                   hash_values[hash] += b_val;
                                   break;
                                 } else {
                                   hash = (hash + 1) & pow2_hash_func;
                                 }
                               }
                             }
                           }
                           size_type c_row_begin = row_mapC[row_index];
                           for (nnz_lno_t i = 0; i < used_count; ++i) {
                             nnz_lno_t used_index    = used_indices[i];
                             pEntriesC[c_row_begin]  = hash_ids[used_index];
                             pvaluesC[c_row_begin++] = hash_values[used_index];
                             hash_ids[used_index]    = -1;
                           }
                         });
    memory_space.release_chunk(used_indices);
  }

  // assumes that the vector lane is 1, as in cpus
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, scalar_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(pow2_hash_size, pow2_hash_func, NULL, NULL, NULL, NULL);

    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(memory_space.allocate_chunk(tid));
    }
    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *)tmp;
    tmp += pow2_hash_size;

    hm2.hash_begins = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;
    hm2.hash_nexts = (nnz_lno_t *)(tmp);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          nnz_lno_t globally_used_hash_count = 0;
          nnz_lno_t used_hash_sizes          = 0;

          const size_type c_row_begin = row_mapC[row_index];

          hm2.keys   = pEntriesC + c_row_begin;
          hm2.values = pvaluesC + c_row_begin;

          // Insert B
          size_type rowBegin   = row_mapB(row_index);
          nnz_lno_t left_workB = row_mapB(row_index + 1) - rowBegin;

          for (nnz_lno_t i = 0; i < left_workB; ++i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t b_col_ind    = entriesB[adjind];
            scalar_t b_val         = valuesB[adjind];

            hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(b_col_ind, b_val, &used_hash_sizes,
                                                                 &globally_used_hash_count, globally_used_hash_indices);
          }

          // Insert -omega * dinv * A*B
          const scalar_t mult       = -omega * dinv(row_index, 0);
          const size_type col_begin = row_mapA[row_index];
          const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;

          for (nnz_lno_t ii = 0; ii < left_work; ++ii) {
            size_type a_col = col_begin + ii;
            nnz_lno_t rowB  = entriesA[a_col];
            scalar_t valA   = valuesA[a_col] * mult;

            rowBegin   = row_mapB(rowB);
            left_workB = row_mapB(rowB + 1) - rowBegin;

            for (nnz_lno_t i = 0; i < left_workB; ++i) {
              const size_type adjind = i + rowBegin;
              nnz_lno_t b_col_ind    = entriesB[adjind];
              scalar_t b_val         = valuesB[adjind] * valA;

              hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
                  b_col_ind, b_val, &used_hash_sizes, &globally_used_hash_count, globally_used_hash_indices);
            }
          }

          for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i) {
            nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
            hm2.hash_begins[dirty_hash] = -1;
          }
        });
    memory_space.release_chunk(globally_used_hash_indices);
  }

  // This is the default row-based algorithm (thread-sequential)
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // Shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    // Holds the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_hash_size;

    // Points to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_key_size;

    // Holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_key_size;

    // Remainder of shmem allocation for vals
    scalar_t *vals = KokkosKernels::Impl::alignPtrTo<scalar_t>(all_shared_memory);

    // Create the hashmaps
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, scalar_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm(thread_shmem_key_size, thread_shmem_hash_func, begins, nexts, keys, vals);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, scalar_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(pow2_hash_size, pow2_hash_func, NULL, NULL, NULL, NULL);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          scalar_t mult                           = -omega * dinv(row_index, 0);
          const size_type c_row_begin             = row_mapC[row_index];
          const size_type c_row_end               = row_mapC[row_index + 1];
          const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
          bool is_global_alloced                  = false;
          nnz_lno_t *globally_used_hash_indices   = NULL;

          // TODO: HashmapAccumulator should encapsulate growing the linked
          // lists.
          if (global_memory_hash_size > thread_shmem_key_size) {
            volatile nnz_lno_t *tmp = NULL;
            size_t tid              = row_index;

            while (tmp == NULL) {
              Kokkos::single(
                  Kokkos::PerThread(teamMember),
                  [&](volatile nnz_lno_t *&memptr) {
                    memptr = (volatile nnz_lno_t *)(memory_space.allocate_chunk(tid));
                  },
                  tmp);
            }

            is_global_alloced          = true;
            globally_used_hash_indices = (nnz_lno_t *)tmp;
            tmp += pow2_hash_size;
            hm2.hash_begins = (nnz_lno_t *)(tmp);
            tmp += pow2_hash_size;
            hm2.hash_nexts = (nnz_lno_t *)(tmp);
          }
          hm2.keys   = pEntriesC + c_row_begin;
          hm2.values = pvaluesC + c_row_begin;

          // Initialize begins.
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, thread_shmem_hash_size),
                               [&](nnz_lno_t i) { begins[i] = -1; });

          // Initialize hash usage sizes
          Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
            used_hash_sizes[0]          = 0;
            used_hash_sizes[1]          = 0;
            globally_used_hash_count[0] = 0;
          });

          // Insert B
          size_type rowBegin   = row_mapB(row_index);
          nnz_lno_t left_work_ = row_mapB(row_index + 1) - rowBegin;
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work_), [&](nnz_lno_t i) {
            const size_type adjind     = i + rowBegin;
            nnz_lno_t b_col_ind        = entriesB[adjind];
            scalar_t b_val             = valuesB[adjind];
            volatile int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(b_col_ind, b_val, used_hash_sizes);
            if (num_unsuccess) {
              hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                  b_col_ind, b_val, used_hash_sizes + 1, globally_used_hash_count, globally_used_hash_indices);
            }
          });

          // Insert -omega * Dinv * A*B
          const size_type col_begin = row_mapA[row_index];
          const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
          nnz_lno_t ii              = left_work;

          while (ii-- > 0) {
            size_type a_col = col_begin + ii;
            nnz_lno_t rowB  = entriesA[a_col];
            scalar_t valA   = valuesA[a_col] * mult;
            rowBegin        = row_mapB(rowB);
            left_work_      = row_mapB(rowB + 1) - rowBegin;
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, left_work_), [&](nnz_lno_t i) {
              const size_type adjind = i + rowBegin;
              nnz_lno_t b_col_ind    = entriesB[adjind];
              scalar_t b_val         = valuesB[adjind] * valA;
              volatile int num_unsuccess =
                  hm.vector_atomic_insert_into_hash_mergeAdd(b_col_ind, b_val, used_hash_sizes);
              if (num_unsuccess) {
                hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                    b_col_ind, b_val, used_hash_sizes + 1, globally_used_hash_count, globally_used_hash_indices);
              }
            });
          }

          if (is_global_alloced) {
            nnz_lno_t dirty_hashes = globally_used_hash_count[0];
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, dirty_hashes), [&](nnz_lno_t i) {
              nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
              hm2.hash_begins[dirty_hash] = -1;
            });

            Kokkos::single(Kokkos::PerThread(teamMember),
                           [&]() { memory_space.release_chunk(globally_used_hash_indices); });
          }

          Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
            if (used_hash_sizes[0] > thread_shmem_key_size) used_hash_sizes[0] = thread_shmem_key_size;
          });

          nnz_lno_t num_elements  = used_hash_sizes[0];
          nnz_lno_t written_index = used_hash_sizes[1];
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, num_elements), [&](nnz_lno_t i) {
            pEntriesC[c_row_begin + written_index + i] = keys[i];
            pvaluesC[c_row_begin + written_index + i]  = vals[i];
          });
        });
  }

  // The thread-flat-parallel implementation for the case one row fits into
  // shmem
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag4 &, const team_member_t &teamMember) const {
    const nnz_lno_t init_value   = -1;
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    char *all_shared_memory             = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * team_cuckoo_key_size;
    scalar_t *vals = KokkosKernels::Impl::alignPtrTo<scalar_t>(all_shared_memory);

    int thread_rank = teamMember.team_rank();
    int vector_rank = 0;
    typedef typename std::remove_reference<decltype(*used_hash_sizes)>::type atomic_incr_type;
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(teamMember, vector_size),
                          [&](const int /* threadid */, int &update, const bool final) {
                            if (final) {
                              vector_rank = update;
                            }
                            update += 1;
                          });

    int bs           = vector_size * suggested_team_size;
    int vector_shift = thread_rank * vector_size + vector_rank;

    for (nnz_lno_t row_index = team_row_begin; row_index < team_row_end; ++row_index) {
      teamMember.team_barrier();

      scalar_t mult               = -omega * dinv(row_index, 0);
      const size_type c_row_begin = row_mapC[row_index];
      nnz_lno_t *c_row            = entriesC.data() + c_row_begin;
      scalar_t *c_row_vals        = valuesC.data() + c_row_begin;

      // Initialize hashmaps
      nnz_lno_t num_threads = team_cuckoo_key_size / vector_size;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_threads), [&](nnz_lno_t teamind) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, vector_size), [&](nnz_lno_t i) {
          keys[teamind * vector_size + i] = init_value;
          vals[teamind * vector_size + i] = 0;
        });
      });

      Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      // Insert B
      teamMember.team_barrier();
      nnz_lno_t current_b_column_flops = row_mapB[row_index + 1] - row_mapB[row_index];
      for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift < current_b_column_flops;
           vector_read_shift += bs) {
        nnz_lno_t hash = init_value;
        int fail       = 0;

        nnz_lno_t my_b_col = entriesB[row_mapB[row_index] + vector_read_shift];
        scalar_t my_b_val  = valuesB[row_mapB[row_index] + vector_read_shift];

        // now insert it to first level hashmap accumulator.
        hash = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
        fail = 1;
        for (nnz_lno_t trial = hash; trial < team_cuckoo_key_size;) {
          if (keys[trial] == my_b_col) {
            Kokkos::atomic_add(vals + trial, my_b_val);
            fail = 0;
            break;
          } else if (keys[trial] == init_value) {
            if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            }
          } else {
            ++trial;
          }
        }

        if (fail) {
          for (nnz_lno_t trial = 0; trial < hash;) {
            if (keys[trial] == my_b_col) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            } else if (keys[trial] == init_value) {
              if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
                Kokkos::atomic_add(vals + trial, my_b_val);
                fail = 0;
                break;
              }
            } else {
              ++trial;
            }
          }
        }  // end if (fail)
      }    // end for (nnz_lno_t vector_read_shift = vector_shift ...

      // Insert -w Dinv A*B
      const size_type a_col_begin_offset = row_mapA[row_index];

      nnz_lno_t a_col_ind = entriesA[a_col_begin_offset];
      scalar_t a_col_val  = valuesA[a_col_begin_offset] * mult;

      nnz_lno_t current_a_column_offset_inrow = 0;
      nnz_lno_t flops_on_the_left_of_offsett  = 0;
      size_type current_b_read_offsett        = row_mapB[a_col_ind];
      nnz_lno_t current_a_column_flops        = row_mapB[a_col_ind + 1] - current_b_read_offsett;

      nnz_lno_t row_flops = flops_per_row(row_index);

      teamMember.team_barrier();

      for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift < row_flops; vector_read_shift += bs) {
        nnz_lno_t my_b_col_shift = vector_read_shift - flops_on_the_left_of_offsett;
        nnz_lno_t my_b_col       = init_value;
        scalar_t my_b_val        = 0;
        nnz_lno_t hash           = init_value;
        int fail                 = 0;

        if (my_b_col_shift >= current_a_column_flops) {
          do {
            ++current_a_column_offset_inrow;
            my_b_col_shift -= current_a_column_flops;
            flops_on_the_left_of_offsett += current_a_column_flops;
            a_col_ind = entriesA[a_col_begin_offset + current_a_column_offset_inrow];

            current_b_read_offsett = row_mapB[a_col_ind];
            current_a_column_flops = row_mapB[a_col_ind + 1] - current_b_read_offsett;
          } while (my_b_col_shift >= current_a_column_flops);
          a_col_val = valuesA[a_col_begin_offset + current_a_column_offset_inrow];
        }

        my_b_col = entriesB[my_b_col_shift + current_b_read_offsett];

        my_b_val = valuesB[my_b_col_shift + current_b_read_offsett] * a_col_val;

        // now insert it to first level hashmap accumulator.
        hash = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
        fail = 1;
        for (nnz_lno_t trial = hash; trial < team_cuckoo_key_size;) {
          if (keys[trial] == my_b_col) {
            Kokkos::atomic_add(vals + trial, my_b_val);
            fail = 0;
            break;
          } else if (keys[trial] == init_value) {
            if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            }
          } else {
            ++trial;
          }
        }

        if (fail) {
          for (nnz_lno_t trial = 0; trial < hash;) {
            if (keys[trial] == my_b_col) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            } else if (keys[trial] == init_value) {
              if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
                Kokkos::atomic_add(vals + trial, my_b_val);
                fail = 0;
                break;
              }
            } else {
              ++trial;
            }
          }
        }  // end if (fail)
      }    // end for (nnz_lno_t vector_read_shift = vector_shift ...

      teamMember.team_barrier();
      for (nnz_lno_t my_index = vector_shift; my_index < team_cuckoo_key_size; my_index += bs) {
        nnz_lno_t my_key = keys[my_index];
        if (my_key != init_value) {
          scalar_t my_val         = vals[my_index];
          nnz_lno_t write_index   = Kokkos::atomic_fetch_add(used_hash_sizes, atomic_incr_type(1));
          c_row[write_index]      = my_key;
          c_row_vals[write_index] = my_val;
        }
      }
    }  // end for (nnz_lno_t  my_index = vector_shift ...
  }

  // The thread-flat-parallel implementation for the case one row does not fit
  // into shmem
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag6 &, const team_member_t &teamMember) const {
    const nnz_lno_t init_value   = -1;
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    char *all_shared_memory             = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * team_cuckoo_key_size;
    scalar_t *vals = KokkosKernels::Impl::alignPtrTo<scalar_t>(all_shared_memory);

    int thread_rank = teamMember.team_rank();
    int vector_rank = 0;
    typedef typename std::remove_reference<decltype(*used_hash_sizes)>::type atomic_incr_type;
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(teamMember, vector_size),
                          [&](const int /* threadid */, int &update, const bool final) {
                            if (final) {
                              vector_rank = update;
                            }
                            update += 1;
                          });

    int bs           = vector_size * suggested_team_size;
    int vector_shift = thread_rank * vector_size + vector_rank;

    for (nnz_lno_t row_index = team_row_begin; row_index < team_row_end; ++row_index) {
      teamMember.team_barrier();

      scalar_t mult                  = -omega * dinv(row_index, 0);
      const size_type c_row_begin    = row_mapC[row_index];
      const size_type c_row_end      = row_mapC[row_index + 1];
      const nnz_lno_t c_row_size     = c_row_end - c_row_begin;
      nnz_lno_t *c_row               = entriesC.data() + c_row_begin;
      scalar_t *c_row_vals           = valuesC.data() + c_row_begin;
      nnz_lno_t *global_acc_row_keys = c_row;
      scalar_t *global_acc_row_vals  = c_row_vals;
      volatile nnz_lno_t *tmp        = NULL;

      // Initialize hashmaps
      if (c_row_size > max_first_level_hash_size) {
        while (tmp == NULL) {
          Kokkos::single(
              Kokkos::PerTeam(teamMember),
              [&](volatile nnz_lno_t *&memptr) {
                memptr = (volatile nnz_lno_t *)(memory_space.allocate_chunk(row_index));
              },
              tmp);
        }
        global_acc_row_keys = (nnz_lno_t *)(tmp);
        global_acc_row_vals = KokkosKernels::Impl::alignPtrTo<scalar_t>(tmp + pow2_hash_size);

        nnz_lno_t num_threads = pow2_hash_size / vector_size;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_threads), [&](nnz_lno_t teamind) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, vector_size),
                               [&](nnz_lno_t i) { global_acc_row_vals[teamind * vector_size + i] = 0; });
        });
      }

      nnz_lno_t num_threads = team_cuckoo_key_size / vector_size;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_threads), [&](nnz_lno_t teamind) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, vector_size), [&](nnz_lno_t i) {
          keys[teamind * vector_size + i] = init_value;
          vals[teamind * vector_size + i] = 0;
        });
      });

      Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      // Insert B
      teamMember.team_barrier();
      nnz_lno_t current_b_column_flops = row_mapB[row_index + 1] - row_mapB[row_index];
      bool insert_is_on                = true;
      for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift < current_b_column_flops;
           vector_read_shift += bs) {
        nnz_lno_t hash = init_value;
        int fail       = 0;

        nnz_lno_t my_b_col = entriesB[row_mapB[row_index] + vector_read_shift];
        scalar_t my_b_val  = valuesB[row_mapB[row_index] + vector_read_shift];

        // now insert it to first level hashmap accumulator.
        hash                 = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
        fail                 = 1;
        bool try_to_insert   = true;
        nnz_lno_t search_end = team_cuckoo_key_size;
        for (nnz_lno_t trial = hash; trial < search_end;) {
          if (keys[trial] == my_b_col) {
            Kokkos::atomic_add(vals + trial, my_b_val);
            fail = 0;
            break;
          } else if (keys[trial] == init_value) {
            if (!insert_is_on) {
              try_to_insert = false;
              break;
            } else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              Kokkos::atomic_increment(used_hash_sizes);
              if (used_hash_sizes[0] > max_first_level_hash_size) insert_is_on = false;
              fail = 0;
              break;
            }
          } else {
            ++trial;
          }
        }  // end for (nnz_lno_t trial = hash; ...

        if (fail) {
          search_end = hash;

          for (nnz_lno_t trial = 0; try_to_insert && trial < search_end;) {
            if (keys[trial] == my_b_col) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            } else if (keys[trial] == init_value) {
              if (!insert_is_on) {
                break;
              } else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
                Kokkos::atomic_add(vals + trial, my_b_val);
                Kokkos::atomic_increment(used_hash_sizes);
                if (used_hash_sizes[0] > max_first_level_hash_size) insert_is_on = false;
                fail = 0;
                break;
              }
            } else {
              ++trial;
            }
          }

          if (fail) {
            nnz_lno_t new_hash = (my_b_col * HASHSCALAR) & pow2_hash_func;

            for (nnz_lno_t trial = new_hash; trial < pow2_hash_size;) {
              if (global_acc_row_keys[trial] == my_b_col) {
                Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                fail = 0;
                break;
              } else if (global_acc_row_keys[trial] == init_value) {
                if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial, init_value, my_b_col)) {
                  Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                  fail = 0;
                  break;
                }
              } else {
                ++trial;
              }
            }

            if (fail) {
              for (nnz_lno_t trial = 0; trial < new_hash;) {
                if (global_acc_row_keys[trial] == my_b_col) {
                  Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                  break;
                } else if (global_acc_row_keys[trial] == init_value) {
                  if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial, init_value, my_b_col)) {
                    Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                    break;
                  }
                } else {
                  ++trial;
                }
              }
            }
          }
        }  // end if (fail)
      }    // end for (nnz_lno_t vector_read_shift = vector_shift ...

      // Insert - w Dinv A * B
      // insert_is_on = true; // I'm not sure if this is necessary
      const size_type a_col_begin_offset = row_mapA[row_index];

      nnz_lno_t a_col_ind = entriesA[a_col_begin_offset];
      scalar_t a_col_val  = valuesA[a_col_begin_offset] * mult;

      nnz_lno_t current_a_column_offset_inrow = 0;
      nnz_lno_t flops_on_the_left_of_offsett  = 0;
      size_type current_b_read_offsett        = row_mapB[a_col_ind];
      nnz_lno_t current_a_column_flops        = row_mapB[a_col_ind + 1] - current_b_read_offsett;

      nnz_lno_t row_flops = flops_per_row(row_index);

      teamMember.team_barrier();
      for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift < row_flops; vector_read_shift += bs) {
        nnz_lno_t my_b_col_shift = vector_read_shift - flops_on_the_left_of_offsett;
        nnz_lno_t my_b_col       = init_value;
        scalar_t my_b_val        = 0;
        nnz_lno_t hash           = init_value;
        int fail                 = 0;

        if (my_b_col_shift >= current_a_column_flops) {
          do {
            ++current_a_column_offset_inrow;
            my_b_col_shift -= current_a_column_flops;
            flops_on_the_left_of_offsett += current_a_column_flops;
            a_col_ind = entriesA[a_col_begin_offset + current_a_column_offset_inrow];

            current_b_read_offsett = row_mapB[a_col_ind];
            current_a_column_flops = row_mapB[a_col_ind + 1] - current_b_read_offsett;
          } while (my_b_col_shift >= current_a_column_flops);
          a_col_val = valuesA[a_col_begin_offset + current_a_column_offset_inrow];
        }

        my_b_col = entriesB[my_b_col_shift + current_b_read_offsett];
        my_b_val = valuesB[my_b_col_shift + current_b_read_offsett] * a_col_val;

        // now insert it to first level hashmap accumulator.
        hash                 = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
        fail                 = 1;
        bool try_to_insert   = true;
        nnz_lno_t search_end = team_cuckoo_key_size;
        for (nnz_lno_t trial = hash; trial < search_end;) {
          if (keys[trial] == my_b_col) {
            Kokkos::atomic_add(vals + trial, my_b_val);
            fail = 0;
            break;
          } else if (keys[trial] == init_value) {
            if (!insert_is_on) {
              try_to_insert = false;
              break;
            } else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              Kokkos::atomic_increment(used_hash_sizes);
              if (used_hash_sizes[0] > max_first_level_hash_size) insert_is_on = false;
              fail = 0;
              break;
            }
          } else {
            ++trial;
          }
        }  // end for (nnz_lno_t trial = hash ...

        if (fail) {
          search_end = hash;

          for (nnz_lno_t trial = 0; try_to_insert && trial < search_end;) {
            if (keys[trial] == my_b_col) {
              Kokkos::atomic_add(vals + trial, my_b_val);
              fail = 0;
              break;
            } else if (keys[trial] == init_value) {
              if (!insert_is_on) {
                break;
              } else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)) {
                Kokkos::atomic_add(vals + trial, my_b_val);
                Kokkos::atomic_increment(used_hash_sizes);
                if (used_hash_sizes[0] > max_first_level_hash_size) insert_is_on = false;
                fail = 0;
                break;
              }
            } else {
              ++trial;
            }
          }

          if (fail) {
            nnz_lno_t new_hash = (my_b_col * HASHSCALAR) & pow2_hash_func;
            for (nnz_lno_t trial = new_hash; trial < pow2_hash_size;) {
              if (global_acc_row_keys[trial] == my_b_col) {
                Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                fail = 0;
                break;
              } else if (global_acc_row_keys[trial] == init_value) {
                if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial, init_value, my_b_col)) {
                  Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                  fail = 0;
                  break;
                }
              } else {
                ++trial;
              }
            }
            if (fail) {
              for (nnz_lno_t trial = 0; trial < new_hash;) {
                if (global_acc_row_keys[trial] == my_b_col) {
                  Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                  break;
                } else if (global_acc_row_keys[trial] == init_value) {
                  if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial, init_value, my_b_col)) {
                    Kokkos::atomic_add(global_acc_row_vals + trial, my_b_val);
                    break;
                  }
                } else {
                  ++trial;
                }
              }
            }
          }
        }  // end if (fail)

      }  // end for (nnz_lno_t vector_read_shift = vector_shift ...

      teamMember.team_barrier();

      if (tmp != NULL) {
        for (nnz_lno_t my_index = vector_shift; my_index < pow2_hash_size; my_index += bs) {
          nnz_lno_t my_b_col = global_acc_row_keys[my_index];
          if (my_b_col != init_value) {
            scalar_t my_b_val = global_acc_row_vals[my_index];
            int fail          = 1;

            nnz_lno_t hash       = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
            nnz_lno_t search_end = team_cuckoo_key_size;
            for (nnz_lno_t trial = hash; trial < search_end; ++trial) {
              if (keys[trial] == my_b_col) {
                vals[trial] += my_b_val;
                fail = 0;
                break;
              } else if (keys[trial] == init_value) {
                break;
              }
            }
            search_end = hash;

            for (nnz_lno_t trial = 0; trial < search_end; ++trial) {
              if (keys[trial] == my_b_col) {
                vals[trial] += my_b_val;
                fail = 0;
                break;
              } else if (keys[trial] == init_value) {
                break;
              }
            }

            if (fail) {
              nnz_lno_t write_index   = 0;
              write_index             = Kokkos::atomic_fetch_add(used_hash_sizes + 1, atomic_incr_type(1));
              c_row[write_index]      = my_b_col;
              c_row_vals[write_index] = my_b_val;
            }
            global_acc_row_keys[my_index] = init_value;
          }  // end if (my_b_col != init_value)
        }    // end for (nnz_lno_t my_index = vector_shift ...

        teamMember.team_barrier();
        Kokkos::single(Kokkos::PerTeam(teamMember), [&]() { memory_space.release_chunk(global_acc_row_keys); });
      }  // end if (tmp != NULL)

      for (nnz_lno_t my_index = vector_shift; my_index < team_cuckoo_key_size; my_index += bs) {
        nnz_lno_t my_key = keys[my_index];
        if (my_key != init_value) {
          scalar_t my_val         = vals[my_index];
          nnz_lno_t write_index   = 0;
          write_index             = Kokkos::atomic_fetch_add(used_hash_sizes + 1, atomic_incr_type(1));
          c_row[write_index]      = my_key;
          c_row_vals[write_index] = my_val;
        }
      }
    }
  }

  size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
};

// ============================================================================
// This is a Jacobi-fused SPGEMM implementation which uses a sparse accumulator.
// Copied from KokkosSparse_spgemm_impl_kkmem.hpp and updated for the
// jacobi-fused functionality by adding omega and dinv.
// ============================================================================

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t, typename dinv_view_t>
void KokkosSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::KokkosSPGEMM_jacobi_sparseacc(c_row_view_t row_mapC_, c_lno_nnz_view_t entriesC_,
                                                         c_scalar_nnz_view_t valuesC_,
                                                         typename c_scalar_nnz_view_t::const_value_type omega,
                                                         dinv_view_t dinv,
                                                         KokkosKernels::Impl::ExecSpaceType lcl_my_exec_space) {
  using pool_memory_space = KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_lno_t>;
  constexpr bool exec_gpu = KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>();
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tSPARSE ACC MODE" << std::endl;
  }
  // Initialize the variables
  KokkosSparse::SPGEMMAlgorithm algorithm_to_run = this->spgemm_algorithm;
  nnz_lno_t brows                                = row_mapB.extent(0) - 1;
  size_type bnnz                                 = valsB.extent(0);
  int suggested_vector_size                      = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size                        = this->handle->get_suggested_team_size(suggested_vector_size);
  size_t shmem_size_to_use                       = shmem_size;
  row_lno_persistent_work_view_t flops_per_row   = this->handle->get_spgemm_handle()->row_flops;
  size_t original_overall_flops                  = this->handle->get_spgemm_handle()->original_overall_flops;
  nnz_lno_t max_nnz          = this->handle->get_spgemm_handle()->template get_max_result_nnz<c_row_view_t>(row_mapC_);
  size_type overall_nnz      = this->handle->get_spgemm_handle()->get_c_nnz();
  nnz_lno_t min_hash_size    = 1;
  size_t chunksize           = 1;
  double first_level_cut_off = this->handle->get_spgemm_handle()->get_first_level_hash_cut_off();
  int hash_scaler            = this->handle->get_spgemm_handle()->get_min_hash_size_scale();
  nnz_lno_t tmp_max_nnz      = max_nnz;

  // Compute the shared memory variables.
  // These are not actually passed to functors requiring shmem, but used for
  // algorithm selection
  constexpr size_t scalarAlignPad =
      (alignof(scalar_t) > alignof(nnz_lno_t)) ? (alignof(scalar_t) - alignof(nnz_lno_t)) : 0;
  nnz_lno_t unit_memory         = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof(scalar_t);
  nnz_lno_t team_shmem_key_size = ((shmem_size_to_use - sizeof(nnz_lno_t) * 4 - scalarAlignPad) / unit_memory);
  // alignment padding is per-thread for algorithms with per-thread hashmap
  nnz_lno_t thread_memory         = ((shmem_size_to_use / suggested_team_size - scalarAlignPad) / 8) * 8;
  nnz_lno_t thread_shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tinitial JacobiSpGEMMSparseAcc -- thread_memory:" << thread_memory
              << " unit_memory:" << unit_memory << " initial key size:" << thread_shmem_key_size << std::endl;
    std::cout << "\t\tinitial JacobiSpGEMMSparseAcc -- team_memory:" << shmem_size_to_use
              << " unit_memory:" << unit_memory << " initial team key size:" << team_shmem_key_size << std::endl;
  }
  nnz_lno_t thread_shmem_hash_size = 1;
  while (thread_shmem_hash_size * 2 <= thread_shmem_key_size) {
    thread_shmem_hash_size = thread_shmem_hash_size * 2;
  }
  nnz_lno_t team_shmem_hash_size = 1;
  while (team_shmem_hash_size * 2 <= team_shmem_key_size) {
    team_shmem_hash_size = team_shmem_hash_size * 2;
  }

  team_shmem_key_size = team_shmem_key_size + ((team_shmem_key_size - team_shmem_hash_size) * sizeof(nnz_lno_t)) /
                                                  (sizeof(nnz_lno_t) * 2 + sizeof(scalar_t));
  team_shmem_key_size = (team_shmem_key_size >> 1) << 1;
  thread_shmem_key_size =
      thread_shmem_key_size + ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t)) /
                                  (sizeof(nnz_lno_t) * 2 + sizeof(scalar_t));
  thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;

  if (hash_scaler == 0) {
    tmp_max_nnz = KOKKOSKERNELS_MACRO_MAX(max_nnz, nnz_lno_t(this->b_col_cnt / this->concurrency + 1));
  } else {
    tmp_max_nnz *= hash_scaler;
  }

  // Choose the SpGEMM algorithm and corresponding parameters
  if (this->spgemm_algorithm == SPGEMM_KK || this->spgemm_algorithm == SPGEMM_KK_LP) {
    if (exec_gpu) {
      size_type average_row_nnz = overall_nnz / this->a_row_cnt;
      size_t average_row_flops  = original_overall_flops / this->a_row_cnt;

      // If we have very low flops per row, or our maximum number of nnz is
      // pretty small, then we decide on the row-based algorithm
      // (SPGEMM_KK_MEMORY).
      if (SPGEMM_KK_LP != this->spgemm_algorithm && (average_row_nnz < 32 || average_row_flops < 256)) {
        algorithm_to_run = SPGEMM_KK_MEMORY;
        while (average_row_nnz > size_type(thread_shmem_key_size) && suggested_vector_size < 32) {
          suggested_vector_size  = suggested_vector_size * 2;
          suggested_vector_size  = KOKKOSKERNELS_MACRO_MIN(32, suggested_vector_size);
          suggested_team_size    = this->handle->get_suggested_team_size(suggested_vector_size);
          thread_memory          = (shmem_size_to_use / 8 / suggested_team_size) * 8;
          thread_shmem_key_size  = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
          thread_shmem_hash_size = 1;
          while (thread_shmem_hash_size * 2 <= thread_shmem_key_size) {
            thread_shmem_hash_size = thread_shmem_hash_size * 2;
          }
          thread_shmem_key_size =
              thread_shmem_key_size +
              ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t) - scalarAlignPad) /
                  (sizeof(nnz_lno_t) * 2 + sizeof(scalar_t));
          thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;
        }

        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\t\tRunning KKMEM with suggested_vector_size:" << suggested_vector_size
                    << " suggested_team_size:" << suggested_team_size << std::endl;
        }
      }
      // Otherwise, we decide on SPGEMM_KK_MEMORYSPREADTEAM or
      // SPGEMM_KK_MEMORY_BIDSPREADTEAM
      else {
        nnz_lno_t tmp_team_cuckoo_key_size =
            ((shmem_size_to_use - sizeof(nnz_lno_t) * 2 - scalarAlignPad) / (sizeof(nnz_lno_t) + sizeof(scalar_t)));
        int team_cuckoo_key_size = 1;
        while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
        suggested_vector_size = 32;
        suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
        algorithm_to_run      = SPGEMM_KK_MEMORY_BIGSPREADTEAM;
        while (average_row_nnz < team_cuckoo_key_size / 2 * (KOKKOSKERNELS_MACRO_MIN(first_level_cut_off + 0.05, 1))) {
          shmem_size_to_use = shmem_size_to_use / 2;
          tmp_team_cuckoo_key_size =
              ((shmem_size_to_use - sizeof(nnz_lno_t) * 2 - scalarAlignPad) / (sizeof(nnz_lno_t) + sizeof(scalar_t)));
          team_cuckoo_key_size = 1;
          while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;

          suggested_team_size = suggested_team_size / 2;
        }
        if (average_row_flops > size_t(2) * suggested_team_size * suggested_vector_size &&
            average_row_nnz >
                size_type(team_cuckoo_key_size) * (KOKKOSKERNELS_MACRO_MIN(first_level_cut_off + 0.05, 1))) {
          shmem_size_to_use = shmem_size_to_use * 2;
          tmp_team_cuckoo_key_size =
              ((shmem_size_to_use - sizeof(nnz_lno_t) * 2 - scalarAlignPad) / (sizeof(nnz_lno_t) + sizeof(scalar_t)));
          team_cuckoo_key_size = 1;
          while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
          suggested_team_size = suggested_team_size * 2;
        }

        suggested_team_size = KOKKOSKERNELS_MACRO_MAX(2, suggested_team_size);

        if (max_nnz < team_cuckoo_key_size * KOKKOSKERNELS_MACRO_MIN(first_level_cut_off + 0.20, 1)) {
          algorithm_to_run = SPGEMM_KK_MEMORY_SPREADTEAM;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY_SPREADTEAM with "
                         "suggested_vector_size:"
                      << suggested_vector_size << " suggested_team_size:" << suggested_team_size
                      << " shmem_size_to_use:" << shmem_size_to_use << std::endl;
          }
        } else {
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY_BIGSPREADTEAM with "
                         "suggested_vector_size:"
                      << suggested_vector_size << " suggested_team_size:" << suggested_team_size
                      << " shmem_size_to_use:" << shmem_size_to_use << std::endl;
          }
        }
      }
    }
    // If non-GPU, we decide whether we want to use a sparse or a dense
    // acumulator
    else {
      bool run_dense               = false;
      nnz_lno_t max_column_cut_off = this->handle->get_spgemm_handle()->MaxColDenseAcc;
      nnz_lno_t col_size           = this->b_col_cnt;
      if (col_size < max_column_cut_off) {
        run_dense = true;
        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\t\tRunning SPGEMM_KK_DENSE col_size:" << col_size
                    << " max_column_cut_off:" << max_column_cut_off << std::endl;
        }
      } else {
        nnz_lno_t tmp_min_hash_size = 1;
        while (tmp_max_nnz > tmp_min_hash_size) {
          tmp_min_hash_size *= 4;
        }

        size_t kkmem_chunksize = tmp_min_hash_size;  // this is for used hash indices
        kkmem_chunksize += tmp_min_hash_size;        // this is for the hash begins
        kkmem_chunksize += max_nnz;                  // this is for hash nexts
        kkmem_chunksize        = kkmem_chunksize * sizeof(nnz_lno_t) + scalarAlignPad;
        size_t dense_chunksize = (col_size + col_size / sizeof(scalar_t) + 1) * sizeof(scalar_t);

        if (kkmem_chunksize >= dense_chunksize * 0.5) {
          run_dense = true;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_SPEED kkmem_chunksize:" << kkmem_chunksize
                      << " dense_chunksize:" << dense_chunksize << std::endl;
          }
        } else {
          run_dense = false;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY col_size:" << col_size
                      << " max_column_cut_off:" << max_column_cut_off << std::endl;
          }
        }
      }

      if (run_dense) {
        this->KokkosSPGEMM_jacobi_denseacc(row_mapC_, entriesC_, valuesC_, omega, dinv, lcl_my_exec_space);
        return;
      }
    }
  }
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tJacobiSpGEMMSparseAcc -- adjusted hashsize:" << thread_shmem_hash_size
              << " thread_shmem_key_size:" << thread_shmem_key_size << std::endl;
    std::cout << "\t\tJacobiSpGEMMSparseAcc -- adjusted team hashsize:" << team_shmem_hash_size
              << " team_shmem_key_size:" << team_shmem_key_size << std::endl;
  }

  // Compute the memory pool size
  if (exec_gpu) {
    if (algorithm_to_run == SPGEMM_KK_MEMORY_SPREADTEAM) {
      tmp_max_nnz = 1;
    }
  }

  if (algorithm_to_run == SPGEMM_KK_LP) {
    while (tmp_max_nnz > min_hash_size) {
      min_hash_size *= 4;
    }
    chunksize = min_hash_size;                                          // this is for used hash keys
    chunksize += max_nnz;                                               // this is for used hash keys
    chunksize += scalarAlignPad;                                        // for padding betwen keys and values
    chunksize += min_hash_size * sizeof(scalar_t) / sizeof(nnz_lno_t);  // this is for the hash values
  } else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGSPREADTEAM) {
    while (tmp_max_nnz > min_hash_size) {
      min_hash_size *= 2;  // try to keep it as low as possible because hashes
                           // are not tracked.
    }
    chunksize = min_hash_size;                                          // this is for used hash keys
    chunksize += scalarAlignPad;                                        // for padding between keys and values
    chunksize += min_hash_size * sizeof(scalar_t) / sizeof(nnz_lno_t);  // this is for the hash values
  } else {
    while (tmp_max_nnz > min_hash_size) {
      min_hash_size *= 4;
    }
    chunksize = min_hash_size;   // this is for used hash indices
    chunksize += min_hash_size;  // this is for the hash begins
    chunksize += max_nnz;        // this is for hash nexts
  }

  nnz_lno_t num_chunks = this->template compute_num_pool_chunks<pool_memory_space>(chunksize * sizeof(nnz_lno_t),
                                                                                   concurrency / suggested_vector_size);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\t max_nnz: " << max_nnz << " chunk_size:" << chunksize << " min_hash_size:" << min_hash_size
              << " concurrency:" << concurrency << " MyExecSpace().concurrency():" << MyExecSpace().concurrency()
              << " numchunks:" << num_chunks << std::endl;
  }

  // Allocate the memory pool
  KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
  if (exec_gpu) {
    my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;
  }

  Kokkos::Timer timer;
  pool_memory_space m_space(num_chunks, chunksize, -1, my_pool_type);
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    m_space.print_memory_pool();
    std::cout << "\t\tPool Alloc Time:" << timer.seconds() << std::endl;
    std::cout << "\t\tPool Size(MB):" << sizeof(nnz_lno_t) * (num_chunks * chunksize) / 1024. / 1024. << std::endl;
  }

  // Initialize the functor
  JacobiSpGEMMSparseAcc<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
                        const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, c_row_view_t,
                        c_lno_nnz_view_t, c_scalar_nnz_view_t, dinv_view_t, pool_memory_space>
      jacobi(a_row_cnt, row_mapA, entriesA, valsA, row_mapB, entriesB, valsB, row_mapC_, entriesC_, valuesC_, omega,
             dinv, shmem_size_to_use, suggested_vector_size, m_space, min_hash_size, max_nnz, suggested_team_size,
             lcl_my_exec_space, team_row_chunk_size, first_level_cut_off, flops_per_row, KOKKOSKERNELS_VERBOSE);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tvector_size:" << suggested_vector_size << " chunk_size:" << team_row_chunk_size
              << " suggested_team_size:" << suggested_team_size << std::endl;
  }
  timer.reset();

  if (exec_gpu) {
    if (algorithm_to_run == SPGEMM_KK_MEMORY_SPREADTEAM) {
      if (thread_shmem_key_size <= 0) {
        std::cout << "KokkosSPGEMM_jacobi_sparseacc "
                     "SPGEMM_KK_MEMORY_SPREADTEAM: Insufficient shmem "
                     "available for key for hash map accumulator - Terminating"
                  << std::endl;
        std::cout << "    thread_shmem_key_size = " << thread_shmem_key_size << std::endl;
        throw std::runtime_error(
            " KokkosSPGEMM_jacobi_sparseacc SPGEMM_KK_MEMORY_SPREADTEAM: "
            "Insufficient shmem available for key for hash map accumulator ");
      }
      Kokkos::parallel_for(
          "KokkosSparse::Jacobi::SparseAcc::GPU::SPGEMM_KK_MEMORY_SPREADTEAM",
          gpu_team_policy4_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), jacobi);
      MyExecSpace().fence();

    } else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGSPREADTEAM) {
      if (thread_shmem_key_size <= 0) {
        std::cout << "KokkosSPGEMM_jacobi_sparseacc "
                     "SPGEMM_KK_MEMORY_BIGSPREADTEAM: Insufficient shmem "
                     "available for key for hash map accumulator - Terminating"
                  << std::endl;
        std::cout << "    thread_shmem_key_size = " << thread_shmem_key_size << std::endl;
        throw std::runtime_error(
            " KokkosSPGEMM_jacobi_sparseacc SPGEMM_KK_MEMORY_BIGSPREADTEAM: "
            "Insufficient shmem available for key for hash map accumulator ");
      }
      Kokkos::parallel_for(
          "KokkosSparse::Jacobi::SparseAcc::GPU::SPGEMM_KK_MEMORY_"
          "BIGSPREADTEAM",
          gpu_team_policy6_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), jacobi);
    } else {
      if (team_shmem_key_size <= 0) {
        std::cout << "KokkosSPGEMM_jacobi_sparseacc SPGEMM_KK_MEMORY: Insufficient "
                     "shmem available for key for hash map accumulator - Terminating"
                  << std::endl;
        std::cout << "    team_shmem_key_size = " << team_shmem_key_size << std::endl;
        throw std::runtime_error(
            " KokkosSPGEMM_jacobi_sparseacc SPGEMM_KK_MEMORY: Insufficient "
            "shmem available for key for hash map accumulator ");
      }
      Kokkos::parallel_for(
          "KokkosSparse::Jacobi::SparseAcc::GPU::SPGEMM_KK_MEMORY",
          gpu_team_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), jacobi);
    }
    MyExecSpace().fence();
  } else {
    if (algorithm_to_run == SPGEMM_KK_LP) {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for("KokkosSparse::Jacobi::SparseAcc::CPU::LP::Dynamic",
                             dynamic_multicore_team_policy4_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size,
                                                              suggested_vector_size),
                             jacobi);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::Jacobi::SparseAcc::CPU::LP::Static",
            multicore_team_policy4_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            jacobi);
      }
    } else {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for("KokkosSparse::Jacobi::SparseAcc::CPU::Dynamic",
                             dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size,
                                                             suggested_vector_size),
                             jacobi);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::Jacobi::SparseAcc::CPU::Static",
            multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            jacobi);
      }
    }
    MyExecSpace().fence();
  }

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tJacobi COMP TIME:" << timer.seconds() << std::endl;
  }
}
}  // namespace Impl
}  // namespace KokkosSparse
