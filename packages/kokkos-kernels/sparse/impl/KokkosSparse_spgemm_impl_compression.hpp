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
namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t,
          typename pool_memory_space>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::SingleStepZipMatrix {
  const nnz_lno_t numrows;
  const row_view_t row_map;
  const nnz_view_t entries;
  const nnz_lno_t compression_bit_mask;
  const int compression_bit_divide_shift;
  const int vector_size;

  const new_row_view_t new_row_map;

  new_nnz_view_t set_index_begins;
  new_nnz_view_t set_index_nexts;

  new_nnz_view_t set_index_entries;
  new_nnz_view_t set_entries;

  nnz_lno_t *pset_index_begins;
  nnz_lno_t *pset_index_nexts;

  nnz_lno_t *pset_index_entries;
  nnz_lno_t *pset_entries;

  const int shared_memory_size;
  nnz_lno_t team_row_chunk_size;
  pool_memory_space memory_space;
  nnz_lno_t pow2_hash_size;
  nnz_lno_t max_row_size;

  const int unit_memory;  // begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;

  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size, pow2_hash_func;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;

  SingleStepZipMatrix(const row_view_t row_map_, const nnz_view_t entries_, const nnz_lno_t compression_bit_mask_,
                      const int compression_bit_divide_shift_, const int vector_size_, new_row_view_t new_row_map_,
                      /*Global memory hash space. in the size of nnz*/
                      new_nnz_view_t set_index_begins_, new_nnz_view_t set_index_nexts_,
                      new_nnz_view_t set_index_entries_, new_nnz_view_t set_entries_, const int shared_mem,
                      const nnz_lno_t team_row_chunk_size_, int suggested_team_size_, bool KOKKOSKERNELS_VERBOSE_,
                      KokkosKernels::Impl::ExecSpaceType my_exec_space_)
      : numrows(row_map_.extent(0) - 1),
        row_map(row_map_),
        entries(entries_),
        compression_bit_mask(compression_bit_mask_),
        compression_bit_divide_shift(compression_bit_divide_shift_),
        vector_size(vector_size_),
        new_row_map(new_row_map_),
        set_index_begins(set_index_begins_),
        set_index_nexts(set_index_nexts_),
        set_index_entries(set_index_entries_),
        set_entries(set_entries_),
        pset_index_begins(set_index_begins_.data()),
        pset_index_nexts(set_index_nexts_.data()),
        pset_index_entries(set_index_entries_.data()),
        pset_entries(set_entries_.data()),
        shared_memory_size(shared_mem),
        team_row_chunk_size(team_row_chunk_size_),

        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size / 8 / suggested_team_size_) * 8),
        shmem_key_size(),
        shared_memory_hash_func(),
        shmem_hash_size(1),
        my_exec_space(my_exec_space_) {
    shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * (4 + vector_size_ * 2)) / unit_memory);
    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tCOMPRESS -- thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " initial key size:" << shmem_key_size << std::endl;
    }
    while (shmem_hash_size * 2 <= shmem_key_size) {
      shmem_hash_size = shmem_hash_size * 2;
    }
    shared_memory_hash_func = shmem_hash_size - 1;
    shmem_key_size          = shmem_key_size + ((shmem_key_size - shmem_hash_size) * sizeof(nnz_lno_t)) /
                                          (sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t));
    shmem_key_size = (shmem_key_size >> 1) << 1;
    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\t\tCOMPRESS -- adjusted hashsize:" << shmem_hash_size << " shmem_key_size:" << shmem_key_size
                << std::endl;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_ind) {
                           // CPU part only compares it with the previous index.
                           size_type rowBegin  = row_map(row_ind);
                           nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

                           nnz_lno_t used_size = 0;
                           if (left_work > 0) {
                             const nnz_lno_t n       = entries(rowBegin);
                             nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;
                             // nnz_lno_t prev_nset = 1;
                             // prev_nset = prev_nset << (n & compression_bit_mask);

                             for (nnz_lno_t i = 1; i < left_work; ++i) {
                               // nnz_lno_t n_set = 1;
                               const size_type adjind = i + rowBegin;
                               const nnz_lno_t nn     = entries(adjind);
                               nnz_lno_t n_set_index  = nn >> compression_bit_divide_shift;
                               // n_set = n_set << (nn & compression_bit_mask);
                               if (prev_nset_ind != n_set_index) {
                                 ++used_size;
                                 prev_nset_ind = n_set_index;
                               }
                               /*
                               if (prev_nset_ind == n_set_index){
                                 //prev_nset = prev_nset | n_set;
                               } else {
                                 ++used_size;
                                 prev_nset_ind = n_set_index;

                               }
                               */
                             }
                             ++used_size;
                           }
                           new_row_map(row_ind) = used_size;
                         });
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

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag2 &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(max_row_size, pow2_hash_func, NULL, NULL, NULL, NULL);

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
    tmp += max_row_size;
    hm2.keys = (nnz_lno_t *)(tmp);
    tmp += max_row_size;
    hm2.values = (nnz_lno_t *)(tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_ind) {
                           nnz_lno_t globally_used_hash_count = 0;
                           // CPU part only compares it with the previous index.
                           size_type rowBegin  = row_map(row_ind);
                           nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

                           nnz_lno_t used_size = 0;
                           if (left_work > 0) {
                             const nnz_lno_t n       = entries(rowBegin);
                             nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;

                             for (nnz_lno_t i = 1; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;
                               const nnz_lno_t nn     = entries(adjind);
                               nnz_lno_t n_set_index  = nn >> compression_bit_divide_shift;
                               if (prev_nset_ind != n_set_index) {
                                 // std::cout << " pow2_hash_func:" << pow2_hash_func << "
                                 // prev_nset_ind:" << prev_nset_ind << std::endl; insert
                                 // prev_nset_ind to hashmap
                                 hm2.sequential_insert_into_hash_TrackHashes(
                                     prev_nset_ind, &used_size, &globally_used_hash_count, globally_used_hash_indices);
                                 //++used_size;
                                 prev_nset_ind = n_set_index;
                               }
                             }
                             // insert prev_nset_ind to hashmap
                             hm2.sequential_insert_into_hash_TrackHashes(
                                 prev_nset_ind, &used_size, &globally_used_hash_count, globally_used_hash_indices);
                             //++used_size;
                           }
                           // get the size of the hashmap
                           new_row_map(row_ind) = used_size;
                           // std::cout << "globally_used_hash_count:" <<
                           // globally_used_hash_count << std::endl;
                           for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i) {
                             hm2.hash_begins[globally_used_hash_indices[i]] = -1;
                           }
                         });
    memory_space.release_chunk(globally_used_hash_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag2 &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(max_row_size, pow2_hash_func, NULL, NULL, NULL, NULL);

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
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_ind) {
          // CPU part only compares it with the previous index.
          size_type rowBegin  = row_map(row_ind);
          nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

          size_type outrowBegin = new_row_map(row_ind);

          hm2.keys = pset_index_entries + outrowBegin;

          hm2.values                         = pset_entries + outrowBegin;
          nnz_lno_t globally_used_hash_count = 0;

          nnz_lno_t used_size = 0;
          if (left_work > 0) {
            const nnz_lno_t n       = entries(rowBegin);
            nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;
            nnz_lno_t prev_nset     = 1;
            prev_nset               = prev_nset << (n & compression_bit_mask);

            for (nnz_lno_t i = 1; i < left_work; ++i) {
              nnz_lno_t n_set        = 1;
              const size_type adjind = i + rowBegin;
              const nnz_lno_t nn     = entries(adjind);
              nnz_lno_t n_set_index  = nn >> compression_bit_divide_shift;
              n_set                  = n_set << (nn & compression_bit_mask);
              if (prev_nset_ind == n_set_index) {
                prev_nset = prev_nset | n_set;
              } else {
                hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                    prev_nset_ind, prev_nset, &used_size, &globally_used_hash_count, globally_used_hash_indices);

                // pset_index_entries[used_size + outrowBegin] = prev_nset_ind ;
                // pset_entries[used_size + outrowBegin] = prev_nset;
                //++used_size;
                prev_nset_ind = n_set_index;
                prev_nset     = n_set;
              }
            }
            hm2.sequential_insert_into_hash_mergeOr_TrackHashes(prev_nset_ind, prev_nset, &used_size,
                                                                &globally_used_hash_count, globally_used_hash_indices);
            for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i) {
              hm2.hash_begins[globally_used_hash_indices[i]] = -1;
            }
            // pset_index_entries[used_size + outrowBegin] = prev_nset_ind ;
            // pset_entries[used_size + outrowBegin] = prev_nset;
            //++used_size;
          }
        });
    memory_space.release_chunk(globally_used_hash_indices);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_ind) {
                           // CPU part only compares it with the previous index.
                           size_type rowBegin  = row_map(row_ind);
                           nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

                           size_type outrowBegin = new_row_map(row_ind);

                           nnz_lno_t used_size = 0;
                           if (left_work > 0) {
                             const nnz_lno_t n       = entries(rowBegin);
                             nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;
                             nnz_lno_t prev_nset     = 1;
                             prev_nset               = prev_nset << (n & compression_bit_mask);

                             for (nnz_lno_t i = 1; i < left_work; ++i) {
                               nnz_lno_t n_set        = 1;
                               const size_type adjind = i + rowBegin;
                               const nnz_lno_t nn     = entries(adjind);
                               nnz_lno_t n_set_index  = nn >> compression_bit_divide_shift;
                               n_set                  = n_set << (nn & compression_bit_mask);
                               if (prev_nset_ind == n_set_index) {
                                 prev_nset = prev_nset | n_set;
                               } else {
                                 pset_index_entries[used_size + outrowBegin] = prev_nset_ind;
                                 pset_entries[used_size + outrowBegin]       = prev_nset;
                                 ++used_size;
                                 prev_nset_ind = n_set_index;
                                 prev_nset     = n_set;
                               }
                             }
                             pset_index_entries[used_size + outrowBegin] = prev_nset_ind;
                             pset_entries[used_size + outrowBegin]       = prev_nset;
                             ++used_size;
                           }
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_ind) {
                           // CPU part only compares it with the previous index.
                           size_type rowBegin  = row_map(row_ind);
                           nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

                           nnz_lno_t used_size = 0;
                           if (left_work > 0) {
                             const nnz_lno_t n       = entries(rowBegin);
                             nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;
                             nnz_lno_t prev_nset     = 1;
                             prev_nset               = prev_nset << (n & compression_bit_mask);

                             for (nnz_lno_t i = 1; i < left_work; ++i) {
                               nnz_lno_t n_set        = 1;
                               const size_type adjind = i + rowBegin;
                               const nnz_lno_t nn     = entries(adjind);
                               nnz_lno_t n_set_index  = nn >> compression_bit_divide_shift;
                               n_set                  = n_set << (nn & compression_bit_mask);
                               if (prev_nset_ind == n_set_index) {
                                 prev_nset = prev_nset | n_set;
                               } else {
                                 pset_index_entries[used_size + rowBegin] = prev_nset_ind;
                                 pset_entries[used_size + rowBegin]       = prev_nset;
                                 ++used_size;
                                 prev_nset_ind = n_set_index;
                                 prev_nset     = n_set;
                               }
                             }
                             pset_index_entries[used_size + rowBegin] = prev_nset_ind;
                             pset_entries[used_size + rowBegin]       = prev_nset;
                             ++used_size;
                           }
                           new_row_map(row_ind) = rowBegin + used_size;
                         });
  }
  // TODO: Implement the GPU count version.

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    nnz_lno_t row_ind = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();
    // std::cout << "i:" << i << std::endl;
    if (row_ind >= numrows) return;

    // std::cout << "i:" << i << std::endl;
    // how much shared memory a thread will have in team
    // int thread_memory = shared_memory_size / teamMember.team_size();
    // allocate all shared memory
    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    // used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    // allocate memory in the size of vectors
    nnz_lno_t *result_keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * vector_size;
    nnz_lno_t *result_vals = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * vector_size;

    // thread_memory -= vector_size * sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t)
    // * 2;

    // calculate the memory needed for 1 single hash entry.
    // int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
    // //begins, nexts, and keys. No need for vals yet. how many hash entries
    // can be held in shared memory. nnz_lno_t shared_memory_hash_size =
    // thread_memory / unit_memory;

    // points to the beginning of hashes
    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    // poins to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    // holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;
    nnz_lno_t *vals = (nnz_lno_t *)(all_shared_memory);

    // this is a hash for individual row elements. therefore, the size can be at
    // most the number of elements in a row, so the nnz_lno_t can be used
    // instead of size_type here.

    // first level hashmap
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm(shmem_hash_size, shared_memory_hash_func, begins, nexts, keys, vals);

    size_type rowBegin  = row_map(row_ind);
    size_type rowBeginP = rowBegin;

    nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;
#ifdef KOKKOSKERNELSMOREMEM
    // same as second level hash map.
    // second level hashmap.
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::modulo>
        hm2(left_work, left_work, pset_index_begins + rowBegin, pset_index_nexts + rowBegin,
            pset_index_entries + rowBegin, pset_entries + rowBegin);
#else
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(left_work, pow2_hash_func, /*pset_index_begins + rowBegin*/ NULL,
            NULL,  // pset_index_nexts+ rowBegin,
            pset_index_entries + rowBegin, pset_entries + rowBegin);
    bool l2_allocated                     = false;
    nnz_lno_t *globally_used_hash_indices = NULL;
#endif

    // initialize begins.
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, shmem_hash_size), [&](int i) { begins[i] = -1; });
    // initialize begins.
    // these are initialized before the loop.
    // initialize hash usage sizes
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      used_hash_sizes[0]          = 0;
      used_hash_sizes[1]          = 0;
      globally_used_hash_count[0] = 0;
    });

    // lno_t neighbor_set_count = 0;

    while (left_work) {
      // std::cout << "left_work:" << left_work << std::endl;
      nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

      // first get the portion of the work for the vector lane.
      nnz_lno_t n_set_index = -1;
      nnz_lno_t n_set       = 1;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, vector_size), [&](nnz_lno_t i) {
        result_keys[i] = -1;
        // This is required for Nvidia architectures with "independent thread
        // scheduling"
#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
        result_vals[i] = 0;
#endif
      });

      // here we find the column_set indices for each column.
      // if it is integer, we divide the column index by 32 with shifts.
      // result_keys is an array of size vector_size
      // this is used as hashtable.
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, work_to_handle), [&](nnz_lno_t i) {
        const size_type adjind = i + rowBegin;
        const nnz_lno_t n      = entries(adjind);
        n_set_index            = n >> compression_bit_divide_shift;
        n_set                  = n_set << (n & compression_bit_mask);

        size_type new_hash = n_set_index & (vector_size - 1);
        nnz_lno_t r        = -1;
        while (true) {
          if (result_keys[new_hash] == n_set_index) {
            Kokkos::atomic_fetch_or(result_vals + new_hash, n_set);
            break;
          } else if (result_keys[new_hash] == r) {
            if (Kokkos::atomic_compare_exchange_strong(result_keys + new_hash, r, n_set_index)) {
              // MD 4/4/18: on these architectures there can be divergence in
              // the warp. once the keys are set, some other vector lane might
              // be doing a fetch_or before we set with n_set. Therefore it is
              // necessary to do atomic, and set it with zero as above.
#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
              Kokkos::atomic_fetch_or(result_vals + new_hash, n_set);
#else
                  result_vals[new_hash] = n_set;
#endif
              break;
            }
          } else if (++new_hash == decltype(new_hash)(vector_size)) {
            new_hash = 0;
          }
        }
      });
      int num_unsuccess         = 0;
      int overall_num_unsuccess = 0;

      // once the values are written result_keys and result_vals
      // each vector reads their corresponding value and inserts it into level-1
      // hash. if level-1 hash is full, it returns 1 (num_unsuccess=1) a
      // reduction is done to see whether any vector lanes failed.
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, vector_size),
          [&](const int i, int &overall_num_unsuccess_) {
            n_set_index = result_keys[i];
            n_set       = result_vals[i];
            if (n_set_index != -1) {
              num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(n_set_index, n_set, used_hash_sizes);
              overall_num_unsuccess_ += num_unsuccess;
            }
          },
          overall_num_unsuccess);

#ifdef KOKKOSKERNELSMOREMEM
      // if one of the inserts was successfull, which means we run out shared
      // memory
      if (overall_num_unsuccess) {
        if (num_unsuccess) {
          hm2.vector_atomic_insert_into_hash_mergeOr(n_set_index, n_set, used_hash_sizes + 1);
        }
      }
#else
      // if one of the inserts was successfull, which means we run out shared
      // memory
      if (overall_num_unsuccess) {
        // then we allocate second level memory using memory pool.
        if (!l2_allocated) {
          volatile nnz_lno_t *tmp = NULL;
          while (tmp == NULL) {
            Kokkos::single(
                Kokkos::PerThread(teamMember),
                [&](volatile nnz_lno_t *&memptr) {
                  memptr = (volatile nnz_lno_t *)(memory_space.allocate_chunk(row_ind));
                },
                tmp);
          }
          globally_used_hash_indices = (nnz_lno_t *)tmp;
          hm2.hash_begins            = (nnz_lno_t *)(globally_used_hash_indices + pow2_hash_size);
          hm2.hash_nexts             = (nnz_lno_t *)(globally_used_hash_indices + pow2_hash_size * 2);
          l2_allocated               = true;
        }

        // this parallel_for is not really needed.
        // we just need a sync threads at the end of the insertion.
        // Basically, we do not want
        // new_row_map(row_ind) = rowBeginP + used_hash_sizes[0] +
        // used_hash_sizes[1]; to execute before the below insertion finishes.
        // parallel_for will provide this mechanism.
#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, vector_size), [&](nnz_lno_t /*i*/) {
#endif
          // then for those who failed we insert it again to L2-accumulator.
          if (num_unsuccess) {
            hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                n_set_index, n_set, used_hash_sizes + 1, globally_used_hash_count, globally_used_hash_indices);
          }
#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
        });
#endif
      }
#endif

      left_work -= work_to_handle;
      rowBegin += work_to_handle;
    }

    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      if (used_hash_sizes[0] > shmem_hash_size) used_hash_sizes[0] = shmem_hash_size;
      new_row_map(row_ind) = rowBeginP + used_hash_sizes[0] + used_hash_sizes[1];
    });

#ifndef KOKKOSKERNELSMOREMEM
    if (l2_allocated) {
      nnz_lno_t dirty_hashes = globally_used_hash_count[0];
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, dirty_hashes), [&](nnz_lno_t i) {
        nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      });

      Kokkos::single(Kokkos::PerThread(teamMember), [&]() { memory_space.release_chunk(globally_used_hash_indices); });
    }
#endif
    size_type written_index = used_hash_sizes[1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, used_hash_sizes[0]), [&](nnz_lno_t i) {
      pset_index_entries[rowBeginP + written_index + i] = keys[i];
      pset_entries[rowBeginP + written_index + i]       = vals[i];
    });
  }

  size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename in_row_view_t, typename in_nnz_view_t, typename out_rowmap_view_t, typename out_nnz_view_t>
bool KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::compressMatrix(nnz_lno_t n, size_type nnz,
                                                                           in_row_view_t in_row_map,
                                                                           in_nnz_view_t in_entries,

                                                                           out_rowmap_view_t out_row_map,
                                                                           out_nnz_view_t &out_nnz_indices,
                                                                           out_nnz_view_t &out_nnz_sets,
                                                                           bool compress_in_single_step) {
  // get the execution space type.
  KokkosKernels::Impl::ExecSpaceType lcl_my_exec_space = this->handle->get_handle_exec_space();
  constexpr bool exec_gpu                              = KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>();
  // get the suggested vectorlane size based on the execution space, and average
  // number of nnzs per row.
  int suggested_vector_size = this->handle->get_suggested_vector_size(n, nnz);
  // get the suggested team size.
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  // get the chunk size suggested by the handle.
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency, n);
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tn:" << n << " nnz:" << nnz << " vector_size:" << suggested_vector_size
              << " team_size:" << suggested_team_size << " chunk_size::" << team_row_chunk_size
              << " shmem:" << shmem_size << std::endl;
  }

  // size of the lno_t, how many bits it can hold.
  int lnot_size                     = sizeof(nnz_lno_t) * 8;
  int compression_bit_divide_shift_ = 0;
  int val                           = lnot_size;
  while (val > 1) {
    ++compression_bit_divide_shift_;
    val = val >> 1;
  }
  nnz_lno_t compression_bit_mask_ = lnot_size - 1;

  Kokkos::Timer timer1;
  // Allocate memory for the linked list to be used for the hashmap
  out_nnz_view_t set_nexts_;
  out_nnz_view_t set_begins_;
#ifdef KOKKOSKERNELSMOREMEM
  if (exec_gpu) {
    set_nexts_  = out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_nexts_"), nnz);
    set_begins_ = out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_begins_"), nnz);
    Kokkos::deep_copy(set_begins_, -1);
  }
  MyExecSpace().fence();
#endif

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tCompression Allocations:" << timer1.seconds() << std::endl;
  }

  // if compressing in single step, allocate the memory as upperbound.
  // TODO: two step is not there for GPU.

  if (compress_in_single_step || exec_gpu) {
    out_nnz_indices = out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_entries_"), nnz);
    out_nnz_sets    = out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_indices_"), nnz);
  }
  typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_lno_t> pool_memory_space;
  // create functor to compress matrix.
  SingleStepZipMatrix<in_row_view_t, in_nnz_view_t, out_rowmap_view_t, out_nnz_view_t, pool_memory_space>
      sszm_compressMatrix(in_row_map,  // input symbolic matrix
                          in_entries,  // input symbolic matrix

                          compression_bit_mask_,          // modular
                          compression_bit_divide_shift_,  // divisor
                          suggested_vector_size,          // vector size.
                          out_row_map,                    // output row map
                          set_begins_,                    // linked list parallel array begin
                          set_nexts_,                     // linked list parallel array next
                          out_nnz_indices,                // output set bits
                          out_nnz_sets,                   // output set indices
                          shmem_size,                     // shared memory size.
                          team_row_chunk_size             // chunksize.
                          ,
                          suggested_team_size, KOKKOSKERNELS_VERBOSE, lcl_my_exec_space);
  double min_reduction      = this->handle->get_spgemm_handle()->get_compression_cut_off();
  size_t OriginaltotalFlops = this->handle->get_spgemm_handle()->original_overall_flops;

  timer1.reset();
  // bool compression_applied = false;
  if (KokkosKernels::Impl::kk_is_gpu_exec_space<typename HandleType::HandleExecSpace>()) {
#ifndef KOKKOSKERNELSMOREMEM
    size_type max_row_nnz = 0;
    KokkosKernels::Impl::view_reduce_maxsizerow<in_row_view_t, MyExecSpace>(n, in_row_map, max_row_nnz);
    MyExecSpace().fence();
    KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;

    nnz_lno_t min_hash_size = 1;
    while (max_row_nnz > (size_type)(min_hash_size)) {
      min_hash_size *= 2;
    }
    size_t chunksize = min_hash_size;  // this is for used hash indices
    chunksize += min_hash_size;        // this is for the hash begins
    chunksize += max_row_nnz;          // this is for hash nexts

    sszm_compressMatrix.pow2_hash_size = min_hash_size;
    sszm_compressMatrix.pow2_hash_func = min_hash_size - 1;

    nnz_lno_t num_chunks = this->template compute_num_pool_chunks<pool_memory_space>(
        chunksize * sizeof(nnz_lno_t), concurrency / suggested_vector_size);

    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tPOOL chunksize:" << chunksize << " num_chunks:" << num_chunks
                << " min_hash_size:" << min_hash_size << " max_row_nnz:" << max_row_nnz << std::endl;
      std::cout << "\t\tPool Alloc MB:" << (sizeof(nnz_lno_t) * num_chunks * chunksize) / 1024. / 1024. << std::endl;
    }
    nnz_lno_t pool_init_val = -1;
    pool_memory_space m_space(num_chunks, chunksize, pool_init_val, my_pool_type);
    MyExecSpace().fence();
    sszm_compressMatrix.memory_space = m_space;
#endif
    Kokkos::parallel_for("KokkosSparse::SingleStepZipMatrix::GPUEXEC",
                         gpu_team_policy_t(n / suggested_team_size + 1, suggested_team_size, suggested_vector_size),
                         sszm_compressMatrix);
  } else {
    if (!compress_in_single_step) {
      bool use_unordered_compress = true;
      if (use_unordered_compress) {
        size_type max_row_nnz = 0;
        if (n) KokkosKernels::Impl::view_reduce_maxsizerow<in_row_view_t, MyExecSpace>(n, in_row_map, max_row_nnz);
        MyExecSpace().fence();
        KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;

        nnz_lno_t min_hash_size = 1;
        while (max_row_nnz > (size_type)(min_hash_size)) {
          min_hash_size *= 2;
        }
        size_t chunksize = min_hash_size;  // this is for used hash indices
        chunksize += min_hash_size;        // this is for the hash begins
        chunksize += max_row_nnz;          // this is for hash nexts
        chunksize += max_row_nnz;          // this is for hash keys //only for counting.

        sszm_compressMatrix.pow2_hash_size = min_hash_size;
        sszm_compressMatrix.pow2_hash_func = min_hash_size - 1;
        sszm_compressMatrix.max_row_size   = max_row_nnz;
        size_t num_chunks                  = concurrency / suggested_vector_size;

        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\tPOOL chunksize:" << chunksize << " num_chunks:" << num_chunks
                    << " min_hash_size:" << min_hash_size << " max_row_nnz:" << max_row_nnz << std::endl;
          std::cout << "\t\tPool Alloc MB:" << (sizeof(nnz_lno_t) * num_chunks * chunksize) / 1024. / 1024.
                    << std::endl;
        }
        nnz_lno_t pool_init_val = -1;
        pool_memory_space m_space(num_chunks, chunksize, pool_init_val, my_pool_type);
        MyExecSpace().fence();
        sszm_compressMatrix.memory_space = m_space;
      }

      Kokkos::Timer timer_count;
      // HashmapAccumulator is populated here
      if (use_unordered_compress)
        Kokkos::parallel_for(
            "KokkosSparse::TwoStepZipMatrix::use_unordered_compress",
            team_count2_policy_t(n / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sszm_compressMatrix);
      else
        Kokkos::parallel_for(
            "KokkosSparse::TwoStepZipMatrix::use_ordered_compress",
            team_count_policy_t(n / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sszm_compressMatrix);

      MyExecSpace().fence();
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\tCompression Count Kernel:" << timer_count.seconds() << std::endl;
      }
      KokkosKernels::Impl::exclusive_parallel_prefix_sum<out_rowmap_view_t, MyExecSpace>(n + 1, out_row_map);

      {
        nnz_lno_t compressed_maxNumRoughZeros = 0;
        size_t compressedoverall_flops        = 0;
        Kokkos::Timer timer1_t;
        auto new_row_mapB_begin = Kokkos::subview(out_row_map, std::make_pair(nnz_lno_t(0), b_row_cnt));
        auto new_row_mapB_end   = Kokkos::subview(out_row_map, std::make_pair(nnz_lno_t(1), b_row_cnt + 1));
        row_lno_persistent_work_view_t compressed_flops_per_row(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "origianal row flops"), a_row_cnt);

        compressed_maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin,
                                                              new_row_mapB_end, compressed_flops_per_row.data());
        KokkosKernels::Impl::kk_reduce_view2<row_lno_persistent_work_view_t, MyExecSpace>(
            a_row_cnt, compressed_flops_per_row, compressedoverall_flops);
        double ratio = 0;
        if (OriginaltotalFlops) ratio = compressedoverall_flops / ((double)OriginaltotalFlops);
        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\tCompressed Max Row Flops:" << compressed_maxNumRoughZeros << std::endl;
          std::cout << "\t\tCompressed Overall Row Flops:" << compressedoverall_flops << std::endl;
          std::cout << "\t\tCompressed Flops ratio:" << ratio << " min_reduction:" << min_reduction << std::endl;
          std::cout << "\t\tCompressed Max Row Flop Calc Time:" << timer1_t.seconds() << std::endl;
        }

        this->handle->get_spgemm_handle()->compressed_max_row_flops = compressed_maxNumRoughZeros;
        this->handle->get_spgemm_handle()->compressed_overall_flops = compressedoverall_flops;
        if (ratio > min_reduction) {
          return false;
        }
      }

      auto d_c_nnz_size = Kokkos::subview(out_row_map, n);
      auto h_c_nnz_size = Kokkos::create_mirror_view(d_c_nnz_size);
      Kokkos::deep_copy(h_c_nnz_size, d_c_nnz_size);
      typename out_rowmap_view_t::non_const_value_type compressed_b_size = h_c_nnz_size();

      // std::cout << "\tcompressed_b_size:" <<compressed_b_size << std::endl;
      out_nnz_indices =
          out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_entries_"), compressed_b_size);
      out_nnz_sets = out_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "set_indices_"), compressed_b_size);

      sszm_compressMatrix.set_index_entries = out_nnz_indices;
      sszm_compressMatrix.set_entries       = out_nnz_sets;

      sszm_compressMatrix.pset_index_entries = out_nnz_indices.data();
      sszm_compressMatrix.pset_entries       = out_nnz_sets.data();
      if (use_unordered_compress)
        Kokkos::parallel_for(
            "KokkosSparse::TwoStepZipMatrix::fill::use_unordered_compress",
            team_fill2_policy_t(n / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sszm_compressMatrix);
      else
        Kokkos::parallel_for(
            "KokkosSparse::TwoStepZipMatrix::fill::use_unordered_compress",
            team_fill_policy_t(n / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sszm_compressMatrix);
      return true;
    } else {
      // USING DYNAMIC SCHEDULE HERE SLOWS DOWN SIGNIFICANTLY WITH HYPERTHREADS
      Kokkos::parallel_for(
          "KokkosSparse::SingleStepZipMatrix::fill::use_unordered_compress",
          multicore_team_policy_t(n / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
          sszm_compressMatrix);
    }
  }
  MyExecSpace().fence();
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tCompression Kernel time:" << timer1.seconds() << std::endl;
  }
  {
    nnz_lno_t compressed_maxNumRoughZeros = 0;
    size_t compressedoverall_flops        = 0;
    Kokkos::Timer timer1_t;
    auto new_row_mapB_begin = in_row_map;
    auto new_row_mapB_end   = out_row_map;
    row_lno_persistent_work_view_t compressed_flops_per_row(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "origianal row flops"), a_row_cnt);

    compressed_maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin,
                                                          new_row_mapB_end, compressed_flops_per_row.data());
    KokkosKernels::Impl::kk_reduce_view2<row_lno_persistent_work_view_t, MyExecSpace>(
        a_row_cnt, compressed_flops_per_row, compressedoverall_flops);
    double ratio = 0;
    if (OriginaltotalFlops) ratio = compressedoverall_flops / ((double)OriginaltotalFlops);
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tCompressed Max Row Flops:" << compressed_maxNumRoughZeros << std::endl;
      std::cout << "\t\tCompressed Overall Row Flops:" << compressedoverall_flops << std::endl;
      std::cout << "\t\tCompressed Flops ratio:" << ratio << " min_reduction:" << min_reduction << std::endl;

      std::cout << "\t\tCompressed Max Row Flop Calc Time:" << timer1_t.seconds() << std::endl;
    }

    this->handle->get_spgemm_handle()->compressed_max_row_flops = compressed_maxNumRoughZeros;
    this->handle->get_spgemm_handle()->compressed_overall_flops = compressedoverall_flops;
    if (ratio > min_reduction) {
      return false;
    }
  }
  return true;

}  // compressMatrix (end)

}  // namespace Impl
}  // namespace KokkosSparse
