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
#include "KokkosKernels_BitUtils.hpp"
#include <unordered_map>
namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
          typename b_compressed_row_view_t, typename b_nnz_view_t,
          typename c_row_view_t,  // typename nnz_lno_temp_work_view_t,
          typename pool_memory_space>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::StructureC_NC {
  const nnz_lno_t numrows;  // num rows in A

  const a_row_view_t row_mapA;  // A row pointers
  const a_nnz_view_t entriesA;  // A column indices

  const b_original_row_view_t row_pointer_begins_B;
  const b_compressed_row_view_t row_pointer_ends_B;
  b_nnz_view_t b_entries;
  b_nnz_view_t entriesSetsB;

  c_row_view_t rowmapC;
  // nnz_lno_temp_work_view_t entriesSetIndicesC;
  // nnz_lno_temp_work_view_t entriesSetsC;

  const nnz_lno_t pow2_hash_size;
  const nnz_lno_t pow2_hash_func;
  const nnz_lno_t MaxRoughNonZero;

  const size_t shared_memory_size;
  const int vector_size;
  pool_memory_space m_space;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;

  const int unit_memory;  // begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size;
  nnz_lno_t team_row_chunk_size;

  /**
   * \brief constructor
   * \param m_: input row size of A
   * \param row_mapA_: row pointers of A
   * \param entriesA_: col indices of A
   * \param row_ptr_begins_B_: beginning of the rows of B
   * \param row_ptr_ends_B_:end of the rows of B
   * \param entriesSetIndicesB_: column set indices of B
   * \param entriesSetsB_: columns sets of B
   * \param rowmapC_: output rowmap C
   * \param hash_size_: global hashmap hash size.
   * \param MaxRoughNonZero_: max flops for row.
   * \param sharedMemorySize_: shared memory size.
   * \param suggested_team_size_: suggested team size
   * \param team_row_chunk_size_: suggested team chunk size
   * \param my_exec_space_ : execution space.
   */
  StructureC_NC(const nnz_lno_t m_, const a_row_view_t row_mapA_, const a_nnz_view_t entriesA_,
                const b_original_row_view_t row_ptr_begins_B_, const b_compressed_row_view_t row_ptr_ends_B_,
                const b_nnz_view_t entries_b, c_row_view_t rowmapC_, const nnz_lno_t hash_size_,
                const nnz_lno_t MaxRoughNonZero_, const size_t sharedMemorySize_, const int suggested_team_size_,
                const nnz_lno_t team_row_chunk_size_, const int vector_size_, pool_memory_space mpool_,
                const KokkosKernels::Impl::ExecSpaceType my_exec_space_, bool KOKKOSKERNELS_VERBOSE_)
      : numrows(m_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_pointer_begins_B(row_ptr_begins_B_),
        row_pointer_ends_B(row_ptr_ends_B_),
        b_entries(entries_b),
        rowmapC(rowmapC_),
        // entriesSetIndicesC(),
        // entriesSetsC(),
        pow2_hash_size(hash_size_),
        pow2_hash_func(hash_size_ - 1),
        MaxRoughNonZero(MaxRoughNonZero_),
        shared_memory_size(sharedMemorySize_),
        vector_size(vector_size_),
        m_space(mpool_),
        my_exec_space(my_exec_space_),

        // unit memory for a hashmap entry. assuming 1 begin, 1 next, 1 key 1
        // value.
        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size / 8 / suggested_team_size_) * 8),
        shmem_key_size(),
        shared_memory_hash_func(),
        shmem_hash_size(1),
        team_row_chunk_size(team_row_chunk_size_) {
    // how many keys I can hold?
    // thread memory - 3 needed entry for size.
    shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);

    // put the hash size closest power of 2.
    // we round down here, because we want to store more keys,
    // conflicts are cheaper.
    while (shmem_hash_size * 2 <= shmem_key_size) {
      shmem_hash_size = shmem_hash_size * 2;
    }
    // for and opeation we get -1.
    shared_memory_hash_func = shmem_hash_size - 1;

    // increase the key size wit the left over from hash size.
    shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size)) / 3;
    // round it down to 2, because of some alignment issues.
    shmem_key_size = (shmem_key_size >> 1) << 1;

    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\tStructureC "
                << " thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " adjusted hashsize:" << shmem_hash_size << " adjusted shmem_key_size:" << shmem_key_size
                << " using " << (shmem_key_size * 3 + shmem_hash_size) * sizeof(nnz_lno_t) + sizeof(nnz_lno_t) * 3
                << " of thread_memory: " << thread_memory << std::endl;
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

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    // get memory from memory pool.
    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    // set first to globally used hash indices.
    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *)tmp;
    tmp += pow2_hash_size;

    // create hashmap accumulator.
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    // set memory for hash begins.
    hm2.hash_begins = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;

    hm2.hash_nexts = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;

    // holds the keys
    hm2.keys = (nnz_lno_t *)(tmp);
    // tmp += MaxRoughNonZero;
    // hm2.values = (nnz_lno_t *) (tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t globally_used_hash_count = 0;
                           nnz_lno_t used_hash_size           = 0;
                           const size_type col_begin          = row_mapA[row_index];
                           const nnz_lno_t col_size           = row_mapA[row_index + 1] - col_begin;
                           // traverse columns of A.
                           for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];

                             size_type rowBegin  = row_pointer_begins_B(rowB);
                             nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;
                             // traverse columns of B
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;

                               nnz_lno_t b_set_ind = b_entries[adjind];
                               // nnz_lno_t b_set = entriesSetsB[adjind];

                               // insert it to first hash.
                               hm2.sequential_insert_into_hash_TrackHashes(
                                   b_set_ind, &used_hash_size, &globally_used_hash_count, globally_used_hash_indices);
                             }
                           }

                           // when done with all insertions, traverse insertions and get the
                           // size.
                           nnz_lno_t num_el = used_hash_size;

                           // clear the begins.
                           for (int i = 0; i < globally_used_hash_count; ++i) {
                             nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
                             hm2.hash_begins[dirty_hash] = -1;
                           }
                           // set the row size.
                           rowmapC(row_index) = num_el;
                         });

    m_space.release_chunk(globally_used_hash_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreDenseAccumulatorTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    // dense accumulators
    nnz_lno_t *indices      = NULL;
    nnz_lno_t *sets         = NULL;
    volatile nnz_lno_t *tmp = NULL;

    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    // we need as much as column size for sets.
    sets = (nnz_lno_t *)tmp;
    tmp += MaxRoughNonZero;  // this is set as column size before calling dense
                             // accumulators.
    // indices only needs max row size.
    indices = (nnz_lno_t *)tmp;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t index_cnt       = 0;
                           const size_type col_begin = row_mapA[row_index];
                           const nnz_lno_t col_size  = row_mapA[row_index + 1] - col_begin;

                           // traverse columns of A
                           for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
                             size_type a_col = colind + col_begin;

                             nnz_lno_t rowB      = entriesA[a_col];
                             size_type rowBegin  = row_pointer_begins_B(rowB);
                             nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;

                             // traverse columns of B
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t b_set_ind    = b_entries[adjind];
                               // nnz_lno_t b_set = entriesSetsB[adjind];

                               // if sets are not set before, add this to indices.
                               if (sets[b_set_ind] == 0) {
                                 indices[index_cnt++] = b_set_ind;
                                 sets[b_set_ind]      = 1;
                               }
                             }
                           }
                           for (nnz_lno_t ii = 0; ii < index_cnt; ++ii) {
                             nnz_lno_t set_ind = indices[ii];
                             sets[set_ind]     = 0;
                           }
                           rowmapC(row_index) = index_cnt;
                         });

    m_space.release_chunk(sets);
  }
  // this one will be LP on CPUs
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag4 &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    // get memory from memory pool.
    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    nnz_lno_t *used_indices = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;
    nnz_lno_t *hash_ids = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;
    // nnz_lno_t *hash_values = (nnz_lno_t *) (tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t used_count = 0;

                           // nnz_lno_t globally_used_hash_count = 0;
                           // nnz_lno_t used_hash_size = 0;
                           const size_type col_begin = row_mapA[row_index];
                           const nnz_lno_t col_size  = row_mapA[row_index + 1] - col_begin;
                           // traverse columns of A.
                           for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];

                             size_type rowBegin  = row_pointer_begins_B(rowB);
                             nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;
                             // traverse columns of B
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;

                               nnz_lno_t b_set_ind = b_entries[adjind];
                               // nnz_lno_t b_set = entriesSetsB[adjind];
                               nnz_lno_t hash = (b_set_ind * HASHSCALAR) & pow2_hash_func;

                               while (true) {
                                 if (hash_ids[hash] == -1) {
                                   used_indices[used_count++] = hash;
                                   hash_ids[hash]             = b_set_ind;
                                   // hash_values[hash] = b_set;
                                   break;
                                 } else if (hash_ids[hash] == b_set_ind) {
                                   // hash_values[hash] = hash_values[hash] | b_set;
                                   break;
                                 } else {
                                   hash = (hash + 1) & pow2_hash_func;
                                 }
                               }
                             }
                           }

                           for (nnz_lno_t ii = 0; ii < used_count; ++ii) {
                             nnz_lno_t used_index = used_indices[ii];
                             hash_ids[used_index] = -1;
                           }
                           rowmapC(row_index) = used_count;
                         });

    m_space.release_chunk(used_indices);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    nnz_lno_t row_index = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();

    using hashmapType =
        KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                        KokkosKernels::Experimental::HashOpType::bitwiseAnd>;

    if (row_index >= numrows) return;

    // printf("row:%d\n", row_index);

    // int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) *
    // 4;
    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // nnz_lno_t *alloc_global_memory = NULL;
    nnz_lno_t *globally_used_hash_indices = NULL;

    // shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    // used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *)(all_shared_memory);

    all_shared_memory += sizeof(nnz_lno_t);
    // int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
    // nnz_lno_t shmem_key_size = (thread_memory - sizeof(nnz_lno_t) * 3) /
    // unit_memory;

    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    // poins to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;

    // holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    // all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;
    nnz_lno_t *vals = NULL;

    // printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts,
    // keys, vals); return; first level hashmap
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm(shmem_key_size, shared_memory_hash_func, begins, nexts, keys, vals);

    hashmapType hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    // initialize begins.
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, shmem_hash_size), [&](int i) { begins[i] = -1; });

    // initialize hash usage sizes
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      used_hash_sizes[0]          = 0;
      used_hash_sizes[1]          = 0;
      globally_used_hash_count[0] = 0;
    });

    bool is_global_alloced = false;

    const size_type col_end   = row_mapA[row_index + 1];
    const size_type col_begin = row_mapA[row_index];
    const nnz_lno_t col_size  = col_end - col_begin;

    for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
      size_type a_col = colind + col_begin;

      nnz_lno_t rowB     = entriesA[a_col];
      size_type rowBegin = row_pointer_begins_B(rowB);

      nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;

      while (left_work) {
        nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

        nnz_lno_t b_set_ind = -1;  // , b_set = -1;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, work_to_handle), [&](nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          b_set_ind              = b_entries[adjind];
          // b_set = entriesSetsB[adjind];
        });

        int num_unsuccess = hm.vector_atomic_insert_into_hash(b_set_ind, used_hash_sizes);

        int overall_num_unsuccess = 0;

        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&](const int /* threadid */, int &overall_num_unsuccess_) { overall_num_unsuccess_ += num_unsuccess; },
            overall_num_unsuccess);

        if (overall_num_unsuccess) {
          // printf("row:%d\n", row_index);
          if (!is_global_alloced) {
            volatile nnz_lno_t *tmp = NULL;
            size_t tid              = get_thread_id(row_index);
            while (tmp == NULL) {
              Kokkos::single(
                  Kokkos::PerThread(teamMember),
                  [&](volatile nnz_lno_t *&memptr) { memptr = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid)); },
                  tmp);
            }
            is_global_alloced = true;

            globally_used_hash_indices = (nnz_lno_t *)tmp;
            tmp += pow2_hash_size;

            hm2.hash_begins = (nnz_lno_t *)(tmp);
            tmp += pow2_hash_size;

            // poins to the next elements
            hm2.hash_nexts = (nnz_lno_t *)(tmp);
            tmp += MaxRoughNonZero;

            // holds the keys
            hm2.keys = (nnz_lno_t *)(tmp);
            // tmp += MaxRoughNonZero;
            // hm2.values = (nnz_lno_t *) (tmp);
          }

          if (num_unsuccess) {
            // int insertion =
            hm2.vector_atomic_insert_into_hash_TrackHashes(b_set_ind, used_hash_sizes + 1, globally_used_hash_count,
                                                           globally_used_hash_indices);
          }
        }
        left_work -= work_to_handle;
        rowBegin += work_to_handle;
      }
    }

    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
    });

    nnz_lno_t num_elements = used_hash_sizes[0];

    if (is_global_alloced) {
      num_elements += used_hash_sizes[1];

      // now thread leaves the memory as it finds. so there is no need to
      // initialize the hash begins
      nnz_lno_t dirty_hashes = globally_used_hash_count[0];
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, dirty_hashes), [&](nnz_lno_t i) {
        nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      });

      Kokkos::single(Kokkos::PerThread(teamMember), [&]() { m_space.release_chunk(globally_used_hash_indices); });
    }
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() { rowmapC(row_index) = num_elements; });
  }

  ////

  size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
          typename b_compressed_row_view_t, typename b_nnz_view_t,
          typename c_row_view_t,  // typename nnz_lno_temp_work_view_t,
          typename pool_memory_space>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::StructureC {
  const nnz_lno_t numrows;  // num rows in A

  const a_row_view_t row_mapA;  // A row pointers
  const a_nnz_view_t entriesA;  // A column indices

  const b_original_row_view_t row_pointer_begins_B;
  const b_compressed_row_view_t row_pointer_ends_B;
  b_nnz_view_t entriesSetIndicesB;
  b_nnz_view_t entriesSetsB;

  c_row_view_t rowmapC;
  // nnz_lno_temp_work_view_t entriesSetIndicesC;
  // nnz_lno_temp_work_view_t entriesSetsC;

  const nnz_lno_t pow2_hash_size;
  const nnz_lno_t pow2_hash_func;
  const nnz_lno_t MaxRoughNonZero;

  const size_t shared_memory_size;
  const int vector_size;
  pool_memory_space m_space;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;

  const int unit_memory;  // begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size;
  nnz_lno_t team_row_chunk_size;

  /**
   * \brief constructor
   * \param m_: input row size of A
   * \param row_mapA_: row pointers of A
   * \param entriesA_: col indices of A
   * \param row_ptr_begins_B_: beginning of the rows of B
   * \param row_ptr_ends_B_:end of the rows of B
   * \param entriesSetIndicesB_: column set indices of B [CSI]
   * \param entriesSetsB_: columns sets of B             [CS]
   * \param rowmapC_: output rowmap C
   * \param hash_size_: global hashmap hash size.
   * \param MaxRoughNonZero_: max flops for row.         [upper bound on entries
   * per row] \param sharedMemorySize_: shared memory size. \param
   * suggested_team_size_: suggested team size \param team_row_chunk_size_:
   * suggested team chunk size \param my_exec_space_ : execution space.
   */
  StructureC(const nnz_lno_t m_, const a_row_view_t row_mapA_, const a_nnz_view_t entriesA_,
             const b_original_row_view_t row_ptr_begins_B_, const b_compressed_row_view_t row_ptr_ends_B_,
             const b_nnz_view_t entriesSetIndicesB_, const b_nnz_view_t entriesSetsB_, c_row_view_t rowmapC_,
             const nnz_lno_t hash_size_, const nnz_lno_t MaxRoughNonZero_, const size_t sharedMemorySize_,
             const int suggested_team_size_, const nnz_lno_t team_row_chunk_size_, const int vector_size_,
             pool_memory_space mpool_, const KokkosKernels::Impl::ExecSpaceType my_exec_space_,
             bool KOKKOSKERNELS_VERBOSE_)
      : numrows(m_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_pointer_begins_B(row_ptr_begins_B_),
        row_pointer_ends_B(row_ptr_ends_B_),
        entriesSetIndicesB(entriesSetIndicesB_),
        entriesSetsB(entriesSetsB_),
        rowmapC(rowmapC_),
        // entriesSetIndicesC(),
        // entriesSetsC(),
        pow2_hash_size(hash_size_),
        pow2_hash_func(hash_size_ - 1),
        MaxRoughNonZero(MaxRoughNonZero_),
        shared_memory_size(sharedMemorySize_),
        vector_size(vector_size_),
        m_space(mpool_),
        my_exec_space(my_exec_space_),

        // unit memory for a hashmap entry. assuming 1 begin, 1 next, 1 key 1
        // value.
        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size / 8 / suggested_team_size_) * 8),
        shmem_key_size(),
        shared_memory_hash_func(),
        shmem_hash_size(1),
        team_row_chunk_size(team_row_chunk_size_) {
    // how many keys I can hold?
    // thread memory - 3 needed entry for size.
    shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);

    // put the hash size closest power of 2.
    // we round down here, because we want to store more keys,
    // conflicts are cheaper.
    while (shmem_hash_size * 2 <= shmem_key_size) {
      shmem_hash_size = shmem_hash_size * 2;
    }
    // for and opeation we get -1.
    shared_memory_hash_func = shmem_hash_size - 1;

    // increase the key size wit the left over from hash size.
    shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size)) / 3;
    // round it down to 2, because of some alignment issues.
    shmem_key_size = (shmem_key_size >> 1) << 1;

    if (KOKKOSKERNELS_VERBOSE_) {
      std::cout << "\tStructureC "
                << " thread_memory:" << thread_memory << " unit_memory:" << unit_memory
                << " adjusted hashsize:" << shmem_hash_size << " adjusted shmem_key_size:" << shmem_key_size
                << " using " << (shmem_key_size * 3 + shmem_hash_size) * sizeof(nnz_lno_t) + sizeof(nnz_lno_t) * 3
                << " of thread_memory: " << thread_memory << std::endl;
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

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreDenseAccumulatorTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    // dense accumulators
    nnz_lno_t *indices      = NULL;
    nnz_lno_t *sets         = NULL;
    volatile nnz_lno_t *tmp = NULL;

    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    // we need as much as column size for sets.
    sets = (nnz_lno_t *)tmp;
    tmp += MaxRoughNonZero;  // this is set as column size before calling dense
                             // accumulators.
    // indices only needs max row size.
    indices = (nnz_lno_t *)tmp;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t index_cnt       = 0;
                           const size_type col_begin = row_mapA[row_index];
                           const nnz_lno_t col_size  = row_mapA[row_index + 1] - col_begin;

                           // traverse columns of A
                           for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
                             size_type a_col = colind + col_begin;

                             nnz_lno_t rowB      = entriesA[a_col];
                             size_type rowBegin  = row_pointer_begins_B(rowB);
                             nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;

                             // traverse columns of B
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;
                               nnz_lno_t b_set_ind    = entriesSetIndicesB[adjind];
                               nnz_lno_t b_set        = entriesSetsB[adjind];

                               // if sets are not set before, add this to indices.
                               if (sets[b_set_ind] == 0) {
                                 indices[index_cnt++] = b_set_ind;
                               }
                               // make a union.
                               sets[b_set_ind] |= b_set;
                             }
                           }
                           nnz_lno_t num_el = 0;
                           for (nnz_lno_t ii = 0; ii < index_cnt; ++ii) {
                             nnz_lno_t set_ind = indices[ii];
                             nnz_lno_t c_rows  = sets[set_ind];
                             sets[set_ind]     = 0;

                             nnz_lno_t num_el2 = KokkosKernels::Impl::pop_count(c_rows);
                             /*
                                     //count number of set bits
                                     nnz_lno_t num_el2 = 0;
                                     for (; c_rows; num_el2++) {
                                       c_rows = c_rows & (c_rows - 1); // clear the least
                                significant bit set
                                     }
                                     */
                             num_el += num_el2;
                           }
                           rowmapC(row_index) = num_el;
                         });

    m_space.release_chunk(sets);
  }

  // this one will be LP on CPUs
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag4 &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

    // get memory from memory pool.
    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    nnz_lno_t *used_indices = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;
    nnz_lno_t *hash_ids = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;
    nnz_lno_t *hash_values = (nnz_lno_t *)(tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           nnz_lno_t used_count = 0;

                           // nnz_lno_t globally_used_hash_count = 0;
                           // nnz_lno_t used_hash_size = 0;
                           const size_type col_begin = row_mapA[row_index];
                           const nnz_lno_t col_size  = row_mapA[row_index + 1] - col_begin;
                           // traverse columns of A.
                           for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
                             size_type a_col = colind + col_begin;
                             nnz_lno_t rowB  = entriesA[a_col];

                             size_type rowBegin  = row_pointer_begins_B(rowB);
                             nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;
                             // traverse columns of B
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind = i + rowBegin;

                               nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
                               nnz_lno_t b_set     = entriesSetsB[adjind];
                               nnz_lno_t hash      = (b_set_ind * HASHSCALAR) & pow2_hash_func;

                               while (true) {
                                 if (hash_ids[hash] == -1) {
                                   used_indices[used_count++] = hash;
                                   hash_ids[hash]             = b_set_ind;
                                   hash_values[hash]          = b_set;
                                   break;
                                 } else if (hash_ids[hash] == b_set_ind) {
                                   hash_values[hash] = hash_values[hash] | b_set;
                                   break;
                                 } else {
                                   hash = (hash + 1) & pow2_hash_func;
                                 }
                               }
                             }
                           }

                           // when done with all insertions, traverse insertions and get the
                           // size.
                           nnz_lno_t num_el = 0;
                           for (nnz_lno_t ii = 0; ii < used_count; ++ii) {
                             nnz_lno_t used_index = used_indices[ii];
                             nnz_lno_t c_rows     = hash_values[used_index];
                             hash_ids[used_index] = -1;

                             nnz_lno_t num_el2 = KokkosKernels::Impl::pop_count(c_rows);
                             /*
                                       //the number of set bits.
                                       for (; c_rows; num_el2++) {
                                               c_rows = c_rows & (c_rows - 1); // clear the least
                                significant bit set
                                       }
                                       */
                             num_el += num_el2;
                           }
                           rowmapC(row_index) = num_el;
                         });

    m_space.release_chunk(used_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag &, const team_member_t &teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
    using hashmapType =
        KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                        KokkosKernels::Experimental::HashOpType::bitwiseAnd>;

    // get memory from memory pool.
    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    // set first to globally used hash indices.
    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *)tmp;
    tmp += pow2_hash_size;

    // create hashmap accumulator.
    hashmapType hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    // set memory for hash begins.
    hm2.hash_begins = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;

    hm2.hash_nexts = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;

    // holds the keys
    hm2.keys = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;
    hm2.values = (nnz_lno_t *)(tmp);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_index) {
          nnz_lno_t globally_used_hash_count = 0;
          nnz_lno_t used_hash_size           = 0;
          const size_type col_begin          = row_mapA[row_index];
          const nnz_lno_t col_size           = row_mapA[row_index + 1] - col_begin;
          // traverse columns of A.
          for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
            size_type a_col = colind + col_begin;
            nnz_lno_t rowB  = entriesA[a_col];

            size_type rowBegin  = row_pointer_begins_B(rowB);
            nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;
            // traverse columns of B
            for (nnz_lno_t i = 0; i < left_work; ++i) {
              const size_type adjind = i + rowBegin;

              nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
              nnz_lno_t b_set     = entriesSetsB[adjind];

              // insert it to first hash.
              hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                  b_set_ind, b_set, &used_hash_size, &globally_used_hash_count, globally_used_hash_indices);
            }
          }

          // when done with all insertions, traverse insertions and get the
          // size.
          nnz_lno_t num_el = 0;
          for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii) {
            nnz_lno_t c_rows  = hm2.values[ii];
            nnz_lno_t num_el2 = KokkosKernels::Impl::pop_count(c_rows);
            /*
                    //the number of set bits.
                    for (; c_rows; num_el2++) {
                      c_rows = c_rows & (c_rows - 1); // clear the least
               significant bit set
                    }
                    */
            num_el += num_el2;
          }

          // clear the begins.
          for (int i = 0; i < globally_used_hash_count; ++i) {
            nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
            hm2.hash_begins[dirty_hash] = -1;
          }
          // set the row size.
          rowmapC(row_index) = num_el;
        });

    m_space.release_chunk(globally_used_hash_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    nnz_lno_t row_index = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();

    using hashmapType =
        KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                        KokkosKernels::Experimental::HashOpType::bitwiseAnd>;
    if (row_index >= numrows) return;

    // printf("row:%d\n", row_index);

    // int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) *
    // 4;
    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // nnz_lno_t *alloc_global_memory = NULL;
    nnz_lno_t *globally_used_hash_indices = NULL;

    // shift it to the thread private part
    // Each thread gets a partition and each partition consists of 4 arrays:
    // begin, next, keys/ids, values
    all_shared_memory += thread_memory * teamMember.team_rank();

    // used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *)(all_shared_memory);

    all_shared_memory += sizeof(nnz_lno_t);
    // int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
    // nnz_lno_t shmem_key_size = (thread_memory - sizeof(nnz_lno_t) * 3) /
    // unit_memory;

    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    // poins to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;

    // holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;
    nnz_lno_t *vals = (nnz_lno_t *)(all_shared_memory);

    // printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts,
    // keys, vals); return; first level hashmap
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::bitwiseAnd>
        hm(shmem_key_size, shared_memory_hash_func, begins, nexts, keys, vals);

    hashmapType hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    // initialize begins.
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, shmem_hash_size), [&](int i) { begins[i] = -1; });

    // initialize hash usage sizes
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      used_hash_sizes[0]          = 0;
      used_hash_sizes[1]          = 0;
      globally_used_hash_count[0] = 0;
    });

    bool is_global_alloced = false;

    const size_type col_end   = row_mapA[row_index + 1];
    const size_type col_begin = row_mapA[row_index];
    const nnz_lno_t col_size  = col_end - col_begin;

    for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
      size_type a_col = colind + col_begin;

      nnz_lno_t rowB     = entriesA[a_col];
      size_type rowBegin = row_pointer_begins_B(rowB);

      nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;

      while (left_work) {
        nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

        nnz_lno_t b_set_ind = -1, b_set = -1;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, work_to_handle), [&](nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          b_set_ind              = entriesSetIndicesB[adjind];
          b_set                  = entriesSetsB[adjind];
        });

        int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(b_set_ind, b_set, used_hash_sizes);

        int overall_num_unsuccess = 0;

        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&](const int /* threadid */, int &overall_num_unsuccess_) { overall_num_unsuccess_ += num_unsuccess; },
            overall_num_unsuccess);

        if (overall_num_unsuccess) {
          // printf("row:%d\n", row_index);
          if (!is_global_alloced) {
            volatile nnz_lno_t *tmp = NULL;
            size_t tid              = get_thread_id(row_index);
            while (tmp == NULL) {
              Kokkos::single(
                  Kokkos::PerThread(teamMember),
                  [&](volatile nnz_lno_t *&memptr) { memptr = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid)); },
                  tmp);
            }
            is_global_alloced = true;

            globally_used_hash_indices = (nnz_lno_t *)tmp;
            tmp += pow2_hash_size;

            hm2.hash_begins = (nnz_lno_t *)(tmp);
            tmp += pow2_hash_size;

            // poins to the next elements
            hm2.hash_nexts = (nnz_lno_t *)(tmp);
            tmp += MaxRoughNonZero;

            // holds the keys
            hm2.keys = (nnz_lno_t *)(tmp);
            tmp += MaxRoughNonZero;
            hm2.values = (nnz_lno_t *)(tmp);
          }

          if (num_unsuccess) {
            // int insertion =
            hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                b_set_ind, b_set, used_hash_sizes + 1, globally_used_hash_count, globally_used_hash_indices);
          }
        }
        left_work -= work_to_handle;
        rowBegin += work_to_handle;
      }
    }

    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
    });

    /*
    Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
      if (used_hash_sizes[1] > hm2.max_value_size) used_hash_sizes[1] =
    hm2.max_value_size;
    });
    */

    nnz_lno_t num_elements = 0;

    nnz_lno_t num_compressed_elements = used_hash_sizes[0];

    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
        [&](const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
          nnz_lno_t c_rows = hm.values[ii];
          nnz_lno_t num_el = KokkosKernels::Impl::pop_count(c_rows);

          /*
          nnz_lno_t num_el = 0;
          for (; c_rows; num_el++) {
            c_rows &= c_rows - 1; // clear the least significant bit set
          }
          */
          num_nnz_in_row += num_el;
        },
        num_elements);

    if (is_global_alloced) {
      nnz_lno_t num_global_elements      = 0;
      nnz_lno_t num_compressed_elements_ = used_hash_sizes[1];
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, num_compressed_elements_),
          [&](const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
            nnz_lno_t c_rows = hm2.values[ii];
            nnz_lno_t num_el = KokkosKernels::Impl::pop_count(c_rows);

            /*
            nnz_lno_t num_el = 0;
            for (; c_rows; num_el++) {
              c_rows &= c_rows - 1; // clear the least significant bit set
            }
            */
            num_nnz_in_row += num_el;
          },
          num_global_elements);

      // now thread leaves the memory as it finds. so there is no need to
      // initialize the hash begins
      nnz_lno_t dirty_hashes = globally_used_hash_count[0];
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, dirty_hashes), [&](nnz_lno_t i) {
        nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      });

      Kokkos::single(Kokkos::PerThread(teamMember), [&]() { m_space.release_chunk(globally_used_hash_indices); });
      num_elements += num_global_elements;
    }
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() { rowmapC(row_index) = num_elements; });
  }

  size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_r_view_t, typename a_n_view_t, typename b_oldrow_view_t, typename b_row_view_t>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::PredicMaxRowNNZ {
  nnz_lno_t m;          // num rows
  a_r_view_t row_mapA;  // row pointers of a
  a_n_view_t entriesA;  // col
  b_oldrow_view_t row_begins_B;
  b_row_view_t row_end_indices_B;
  const size_type min_val;
  nnz_lno_t team_row_chunk_size;

  size_type *flops_per_row;
  /**
   * \brief Constructor
   * \param m_: num rows in A.
   * \param row_mapA_: row pointers of A
   * \param entriesA_: col indices of A
   * \param row_begins_B_: row begin indices of B
   * \param row_end_indices_B_: row end indices of B
   * \param team_row_chunk_size_: the number of rows assigned to each team.
   */
  PredicMaxRowNNZ(nnz_lno_t m_, a_r_view_t row_mapA_, a_n_view_t entriesA_,

                  b_oldrow_view_t row_begins_B_, b_row_view_t row_end_indices_B_, nnz_lno_t team_row_chunk_size_,
                  size_type *flops_per_row_ = NULL)
      : m(m_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_begins_B(row_begins_B_),
        row_end_indices_B(row_end_indices_B_),
        min_val(((std::numeric_limits<size_type>::lowest()))),
        team_row_chunk_size(team_row_chunk_size_),
        flops_per_row(flops_per_row_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember, size_type &overal_max) const {
    // get the range of rows for team.
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type col_begin = row_mapA[row_index];
                           const size_type col_end   = row_mapA[row_index + 1];
                           const nnz_lno_t left_work = col_end - col_begin;

                           size_type max_num_results_in_row = 0;

                           // get the size of the rows of B, pointed by row of A
                           Kokkos::parallel_reduce(
                               Kokkos::ThreadVectorRange(teamMember, left_work),
                               [&](nnz_lno_t i, size_type &valueToUpdate) {
                                 const size_type adjind   = i + col_begin;
                                 const nnz_lno_t colIndex = entriesA[adjind];
                                 valueToUpdate += row_end_indices_B(colIndex) - row_begins_B(colIndex);
                               },
                               max_num_results_in_row);

                           if (flops_per_row != NULL) {
                             flops_per_row[row_index] = max_num_results_in_row;
                           }

                           // set max.
                           if (overal_max < max_num_results_in_row) {
                             overal_max = max_num_results_in_row;
                           }
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void join(size_type &dst, const size_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_type &dst) const { dst = min_val; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::PredicMaxRowNNZIntersection {
  const nnz_lno_t m, k;       // num rows
  const size_type *row_mapA;  // row pointers of a
  const nnz_lno_t *entriesA;  // col
  const size_type *row_begins_B;
  const size_type *row_end_indices_B;
  const size_type min_val;
  const nnz_lno_t team_row_chunk_size;
  nnz_lno_t *min_result_row_for_each_row;

  /**
   * \brief Constructor
   * \param m_: num rows in A.
   * \param row_mapA_: row pointers of A
   * \param entriesA_: col indices of A
   * \param row_begins_B_: row begin indices of B
   * \param row_end_indices_B_: row end indices of B
   * \param team_row_chunk_size_: the number of rows assigned to each team.
   */
  PredicMaxRowNNZIntersection(const nnz_lno_t m_, const nnz_lno_t k_, const size_type *row_mapA_,
                              const nnz_lno_t *entriesA_,

                              const size_type *row_begins_B_, const size_type *row_end_indices_B_,
                              const nnz_lno_t team_row_chunk_size_, nnz_lno_t *min_result_row_for_each_row_)
      : m(m_),
        k(k_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_begins_B(row_begins_B_),
        row_end_indices_B(row_end_indices_B_),
        min_val(((std::numeric_limits<size_type>::lowest()))),
        team_row_chunk_size(team_row_chunk_size_),
        min_result_row_for_each_row(min_result_row_for_each_row_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember, nnz_lno_t &overal_max) const {
    // get the range of rows for team.
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type col_begin    = row_mapA[row_index];
                           const size_type col_end      = row_mapA[row_index + 1];
                           const nnz_lno_t left_work    = col_end - col_begin;
                           nnz_lno_t min_num_result_row = -1;
                           if (left_work) {
                             nnz_lno_t min_num_results_in_row = this->k;
                             for (nnz_lno_t i = 0; i < left_work; ++i) {
                               const size_type adjind   = i + col_begin;
                               const nnz_lno_t colIndex = entriesA[adjind];
                               nnz_lno_t rowsize        = row_end_indices_B[colIndex] - row_begins_B[colIndex];
                               if (min_num_results_in_row > rowsize) {
                                 min_num_results_in_row = rowsize;
                                 min_num_result_row     = colIndex;
                               }
                             }

                             // set max.
                             if (overal_max < min_num_results_in_row) {
                               overal_max = min_num_results_in_row;
                             }
                           }
                           min_result_row_for_each_row[row_index] = min_num_result_row;
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void join(size_type &dst, const size_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_type &dst) const { dst = min_val; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_r_view_t, typename a_nnz_view_t, typename b_original_row_view_t, typename b_compressed_row_view_t,
          typename b_nnz_view_t, typename c_row_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::symbolic_c_no_compression(nnz_lno_t m, a_r_view_t row_mapA_,
                                                                   a_nnz_view_t entriesA_,

                                                                   b_original_row_view_t b_rowmap_begin,
                                                                   b_compressed_row_view_t b_rowmap_end,
                                                                   b_nnz_view_t entriesb_, c_row_view_t rowmapC,
                                                                   nnz_lno_t maxNumRoughNonzeros) {
  SPGEMMAlgorithm current_spgemm_algorithm             = this->spgemm_algorithm;
  constexpr bool exec_gpu                              = KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>();
  KokkosKernels::Impl::ExecSpaceType lcl_my_exec_space = this->handle->get_handle_exec_space();
  if (exec_gpu) {
    current_spgemm_algorithm = SPGEMM_KK_MEMORY;
  }
  maxNumRoughNonzeros   = KOKKOSKERNELS_MACRO_MIN(this->b_col_cnt, maxNumRoughNonzeros);
  int shmem_size_to_use = shmem_size;

  typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_lno_t> pool_memory_space;

  // get the number of rows and nonzeroes of B.
  nnz_lno_t brows = b_rowmap_end.extent(0) - 1;
  size_type bnnz  = entriesb_.extent(0);

  int max_vector_size       = KokkosKernels::Impl::kk_get_max_vector_size<MyExecSpace>();
  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);

  // this kernel does not really work well if the vector size is less than 4.
  if (suggested_vector_size < 4 && exec_gpu) {
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tsuggested_vector_size:" << suggested_vector_size << " setting it to 4 for Structure kernel"
                << std::endl;
    }
    suggested_vector_size = 4;
  }
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);

  if (this->spgemm_algorithm == SPGEMM_KK || SPGEMM_KK_LP == this->spgemm_algorithm) {
    if (exec_gpu) {
      // then chose the best method and parameters.
      current_spgemm_algorithm = SPGEMM_KK_MEMORY;
      int estimate_compress    = 8;
#ifdef FIRSTPARAMS
      size_t estimate_max_nnz = maxNumRoughNonzeros / estimate_compress;
#else
      // THIS IS BETTER PARAMETER SELECTION.
      size_t original_overall_flops = this->handle->get_spgemm_handle()->original_overall_flops;
      size_t estimate_max_nnz = (sqrt(maxNumRoughNonzeros) * sqrt(original_overall_flops / m)) / estimate_compress;
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\t\testimate_max_nnz:" << estimate_max_nnz << " maxNumRoughNonzeros:" << maxNumRoughNonzeros
                  << " original_overall_flops / m:" << original_overall_flops / m << std::endl;
      }
#endif
      int unit_memory = (sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2);

      int thread_memory = ((shmem_size_to_use / 8 / suggested_team_size) * 8);

      int shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);

      if (estimate_max_nnz / shmem_key_size > 1) {
        int scale = estimate_max_nnz / shmem_key_size;
        while (scale / 2 > 0) {
          scale                 = scale / 2;
          suggested_vector_size = suggested_vector_size * 2;
        }
        suggested_vector_size = KOKKOSKERNELS_MACRO_MIN(max_vector_size, suggested_vector_size);
        suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
      }
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\t\tRunning KKMEM with suggested_vector_size:" << suggested_vector_size
                  << " suggested_team_size:" << suggested_team_size << std::endl;
      }
    } else {
      nnz_lno_t max_column_cut_off = nnz_lno_t(this->handle->get_spgemm_handle()->MaxColDenseAcc);

      nnz_lno_t col_size = this->b_col_cnt;
      if (col_size < max_column_cut_off) {
        current_spgemm_algorithm = SPGEMM_KK_DENSE;
        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\t\tRunning SPGEMM_KK_SPEED col_size:" << col_size
                    << " max_column_cut_off:" << max_column_cut_off << std::endl;
        }
      } else {
        // round up maxNumRoughNonzeros to closest power of 2.
        nnz_lno_t min_hash_size = 1;
        while (maxNumRoughNonzeros > min_hash_size) {
          min_hash_size *= 2;
        }

        size_t kkmem_chunksize = min_hash_size;  // this is for used hash
                                                 // indices
        kkmem_chunksize += min_hash_size;        // this is for the hash begins
        kkmem_chunksize += maxNumRoughNonzeros;  // this is for hash nexts
        kkmem_chunksize += maxNumRoughNonzeros;  // this is for hash keys

        size_t dense_chunksize = col_size + maxNumRoughNonzeros;

        if (kkmem_chunksize >= dense_chunksize * 0.5) {
          current_spgemm_algorithm = SPGEMM_KK_DENSE;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_SPEED kkmem_chunksize:" << kkmem_chunksize
                      << " dense_chunksize:" << dense_chunksize << std::endl;
          }
        } else {
          current_spgemm_algorithm = SPGEMM_KK_MEMORY;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY col_size:" << col_size
                      << " max_column_cut_off:" << max_column_cut_off << std::endl;
          }
        }
      }
    }
  }

  // round up maxNumRoughNonzeros to closest power of 2.
  nnz_lno_t min_hash_size = 1;

  while (maxNumRoughNonzeros > min_hash_size) {
    min_hash_size *= 2;
  }
  min_hash_size *= this->handle->get_spgemm_handle()->get_min_hash_size_scale();

  size_t chunksize = 1;
  // initizalize value for the mem pool
  int pool_init_val = -1;

  // this should not be executed at all.
  if (current_spgemm_algorithm == SPGEMM_KK_LP) {
    pool_init_val = -1;
    chunksize     = min_hash_size;  // this is for keys
    chunksize += min_hash_size;     // this is for the values
  } else if (current_spgemm_algorithm == SPGEMM_KK_DENSE) {
    // nnz_lno_t max_row_size = KOKKOSKERNELS_MACRO_MIN(b_col_cnt,
    // maxNumRoughNonzeros);
    chunksize = b_col_cnt + maxNumRoughNonzeros;
    // if speed is set, and exec space is cpu, then  we use dense accumulators.
    // or if memspeed is set, and concurrency is not high, we use dense
    // accumulators.
    maxNumRoughNonzeros = b_col_cnt;
    pool_init_val       = 0;
  } else {
    pool_init_val = -1;
    // set the chunksize.
    chunksize = min_hash_size;         // this is for used hash indices
    chunksize += min_hash_size;        // this is for the hash begins
    chunksize += maxNumRoughNonzeros;  // this is for hash nexts
    chunksize += maxNumRoughNonzeros;  // this is for hash keys
  }

  // initizalize value for the mem pool
  KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
  if (exec_gpu) {
    my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;
  }

  nnz_lno_t num_chunks = this->template compute_num_pool_chunks<pool_memory_space>(chunksize * sizeof(nnz_lno_t),
                                                                                   concurrency / suggested_vector_size);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tPool Size (MB):" << (num_chunks * chunksize * sizeof(nnz_lno_t)) / 1024. / 1024.
              << " num_chunks:" << num_chunks << " chunksize:" << chunksize << std::endl;
  }
  Kokkos::Timer timer1;
  pool_memory_space m_space(num_chunks, chunksize, pool_init_val, my_pool_type);
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tPool Alloc Time:" << timer1.seconds() << std::endl;
  }

  StructureC_NC<a_r_view_t, a_nnz_view_t, b_original_row_view_t, b_compressed_row_view_t, b_nnz_view_t, c_row_view_t,
                /* nnz_lno_temp_work_view_t,*/ pool_memory_space>
      sc(m, row_mapA_, entriesA_, b_rowmap_begin, b_rowmap_end, entriesb_, rowmapC, min_hash_size, maxNumRoughNonzeros,
         shmem_size_to_use, suggested_team_size, team_row_chunk_size, suggested_vector_size, m_space, lcl_my_exec_space,
         KOKKOSKERNELS_VERBOSE);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tStructureC vector_size:" << suggested_vector_size << " team_size:" << suggested_team_size
              << " chunk_size:" << team_row_chunk_size << " shmem_size:" << shmem_size << std::endl;
  }

  timer1.reset();

  if (exec_gpu) {
    Kokkos::parallel_for("StructureC_NC::GPU_EXEC",
                         gpu_team_policy_t(m / suggested_team_size + 1, suggested_team_size, suggested_vector_size),
                         sc);
  } else {
    if (current_spgemm_algorithm == SPGEMM_KK_DENSE) {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for("KokkosSparse::StructureC_NC::DENSE_DYNAMIC",
                             dynamic_multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1,
                                                                         suggested_team_size, suggested_vector_size),
                             sc);
      } else {
        Kokkos::parallel_for("KokkosSparse::StructureC_NC::DENSE_STATIC",
                             multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1, suggested_team_size,
                                                                 suggested_vector_size),
                             sc);
      }
    } else if (current_spgemm_algorithm == SPGEMM_KK_LP) {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC_NC::LP_DYNAMIC",
            dynamic_multicore_team_policy4_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sc);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC_NC::LP_STATIC",
            multicore_team_policy4_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), sc);
      }
    } else {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC_NC::KKMEM_DYNAMIC",
            dynamic_multicore_team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sc);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC_NC::KKMEM_STATIC",
            multicore_team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), sc);
      }
    }
  }
  MyExecSpace().fence();
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tStructureC Kernel time:" << timer1.seconds() << std::endl << std::endl;
  }
  typename c_row_view_t::non_const_value_type c_nnz_size = 0;
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(m + 1, rowmapC, c_nnz_size);
  this->handle->get_spgemm_handle()->set_c_nnz(c_nnz_size);
  nnz_lno_t c_max_nnz = KokkosSparse::Impl::graph_max_degree<MyExecSpace, size_type, c_row_view_t>(rowmapC);
  this->handle->get_spgemm_handle()->set_max_result_nnz(c_max_nnz);
}  // end: symbolic_c_no_compression

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_r_view_t, typename a_nnz_view_t, typename b_original_row_view_t, typename b_compressed_row_view_t,
          typename b_nnz_view_t, typename c_row_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::symbolic_c(nnz_lno_t m, a_r_view_t row_mapA_,
                                                                       a_nnz_view_t entriesA_,

                                                                       b_original_row_view_t old_row_mapB,
                                                                       b_compressed_row_view_t row_mapB_,
                                                                       b_nnz_view_t entriesSetIndex,
                                                                       b_nnz_view_t entriesSets,

                                                                       c_row_view_t rowmapC,
                                                                       nnz_lno_t maxNumRoughNonzeros) {
  SPGEMMAlgorithm current_spgemm_algorithm = this->spgemm_algorithm;
  constexpr bool exec_gpu = KokkosKernels::Impl::kk_is_gpu_exec_space<typename HandleType::HandleExecSpace>();
  KokkosKernels::Impl::ExecSpaceType lcl_my_exec_space = this->handle->get_handle_exec_space();
  if (exec_gpu) {
    current_spgemm_algorithm = SPGEMM_KK_MEMORY;
  }

  // get the number of rows and nonzeroes of B.
  nnz_lno_t brows             = row_mapB_.extent(0) - 1;
  size_type bnnz              = entriesSetIndex.extent(0);
  size_type compressed_b_size = bnnz;
  if (exec_gpu) {
    KokkosKernels::Impl::kk_reduce_diff_view<b_original_row_view_t, b_compressed_row_view_t, MyExecSpace>(
        brows, old_row_mapB, row_mapB_, compressed_b_size);
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tcompressed_b_size:" << compressed_b_size << " bnnz:" << bnnz << std::endl;
    }
  }
  int max_vector_size       = KokkosKernels::Impl::kk_get_max_vector_size<MyExecSpace>();
  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, compressed_b_size);

  // this kernel does not really work well if the vector size is less than 4.
  if (suggested_vector_size < 4 && exec_gpu) {
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tsuggested_vector_size:" << suggested_vector_size << " setting it to 4 for Structure kernel"
                << std::endl;
    }
    suggested_vector_size = 4;
  }
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  maxNumRoughNonzeros =
      KOKKOSKERNELS_MACRO_MIN(size_t(this->b_col_cnt / sizeof(nnz_lno_t) * 8 + 1), size_t(maxNumRoughNonzeros));
  int shmem_size_to_use = shmem_size;

  if (this->spgemm_algorithm == SPGEMM_KK || SPGEMM_KK_LP == this->spgemm_algorithm) {
    if (exec_gpu) {
      // then chose the best method and parameters.
      current_spgemm_algorithm = SPGEMM_KK_MEMORY;
      int estimate_compress    = 8;
#ifdef FIRSTPARAMS
      size_t estimate_max_nnz = maxNumRoughNonzeros / estimate_compress;
#else
      size_t original_overall_flops = this->handle->get_spgemm_handle()->compressed_overall_flops;
      size_t estimate_max_nnz       = 0;
      if (m > 0) estimate_max_nnz = (sqrt(maxNumRoughNonzeros) * sqrt(original_overall_flops / m)) / estimate_compress;
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\t\testimate_max_nnz:" << estimate_max_nnz << " maxNumRoughNonzeros:" << maxNumRoughNonzeros
                  << " original_overall_flops / m:" << original_overall_flops / m << std::endl;
      }
#endif

      int unit_memory    = (sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2);
      int thread_memory  = ((shmem_size_to_use / 8 / suggested_team_size) * 8);
      int shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);

#ifdef SECONDPARAMS
      if (estimate_max_nnz / shmem_key_size > 1) {
        int scale = estimate_max_nnz / shmem_key_size;
        while (scale / 2 > 0) {
          scale                 = scale / 2;
          suggested_vector_size = suggested_vector_size * 2;
        }
        suggested_vector_size = KOKKOSKERNELS_MACRO_MIN(max_vector_size, suggested_vector_size);
        suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
      }
#else
      // need to adjust the size of shmem based on average row size.
      int thread_shmem_hash_size = 1;
      while (thread_shmem_hash_size * 2 <= shmem_key_size) {
        thread_shmem_hash_size = thread_shmem_hash_size * 2;
      }
      shmem_key_size = shmem_key_size + (shmem_key_size - thread_shmem_hash_size) / 3;
      shmem_key_size = (shmem_key_size >> 1) << 1;
      while (estimate_max_nnz > size_t(shmem_key_size) && suggested_vector_size < max_vector_size) {
        suggested_vector_size  = suggested_vector_size * 2;
        suggested_vector_size  = KOKKOSKERNELS_MACRO_MIN(max_vector_size, suggested_vector_size);
        suggested_team_size    = this->handle->get_suggested_team_size(suggested_vector_size);
        thread_memory          = (shmem_size_to_use / 8 / suggested_team_size) * 8;
        shmem_key_size         = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);
        thread_shmem_hash_size = 1;
        while (thread_shmem_hash_size * 2 <= shmem_key_size) {
          thread_shmem_hash_size = thread_shmem_hash_size * 2;
        }
        shmem_key_size = shmem_key_size + (shmem_key_size - thread_shmem_hash_size) / 3;
        shmem_key_size = (shmem_key_size >> 1) << 1;
        // thread_shmem_key_size * 2;
      }
#endif
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "\t\t\tRunning KKMEM with suggested_vector_size:" << suggested_vector_size
                  << " suggested_team_size:" << suggested_team_size << std::endl;
      }
    } else {
      nnz_lno_t max_column_cut_off = nnz_lno_t(this->handle->get_spgemm_handle()->MaxColDenseAcc);
      nnz_lno_t col_size           = this->b_col_cnt / (sizeof(nnz_lno_t) * 8) + 1;
      if (col_size < max_column_cut_off) {
        current_spgemm_algorithm = SPGEMM_KK_DENSE;
        if (KOKKOSKERNELS_VERBOSE) {
          std::cout << "\t\t\tRunning SPGEMM_KK_SPEED col_size:" << col_size
                    << " max_column_cut_off:" << max_column_cut_off << std::endl;
        }
      } else {
        // round up maxNumRoughNonzeros to closest power of 2.
        nnz_lno_t min_hash_size = 1;
        while (maxNumRoughNonzeros > min_hash_size) {
          min_hash_size *= 2;
        }
        size_t kkmem_chunksize = min_hash_size;  // this is for used hash
                                                 // indices
        kkmem_chunksize += min_hash_size;        // this is for the hash begins
        kkmem_chunksize += maxNumRoughNonzeros;  // this is for hash nexts
        kkmem_chunksize += maxNumRoughNonzeros;  // this is for hash keys
        kkmem_chunksize += maxNumRoughNonzeros;  // this is for hash values

        size_t dense_chunksize = col_size + maxNumRoughNonzeros;

        if (kkmem_chunksize >= dense_chunksize * 0.5) {
          current_spgemm_algorithm = SPGEMM_KK_DENSE;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_SPEED kkmem_chunksize:" << kkmem_chunksize
                      << " dense_chunksize:" << dense_chunksize << std::endl;
          }
        } else {
          current_spgemm_algorithm = SPGEMM_KK_MEMORY;
          if (KOKKOSKERNELS_VERBOSE) {
            std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY col_size:" << col_size
                      << " max_column_cut_off:" << max_column_cut_off << std::endl;
          }
        }
      }
    }
  }
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency, a_row_cnt);

  typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_lno_t> pool_memory_space;

  // round up maxNumRoughNonzeros to closest power of 2.
  nnz_lno_t min_hash_size = 1;
  while (maxNumRoughNonzeros > min_hash_size) {
    min_hash_size *= 2;
  }

  // set the chunksize.
  size_t chunksize = 1;
  // initizalize value for the mem pool
  int pool_init_val = -1;

  if (current_spgemm_algorithm == SPGEMM_KK_LP) {
    // this should not be executed with the new parameter selectio above.
    pool_init_val = -1;
    chunksize     = min_hash_size;     // this is for keys
    chunksize += min_hash_size;        // this is for the values
    chunksize += maxNumRoughNonzeros;  // this is for hash values
  } else {
    pool_init_val = -1;

    chunksize = min_hash_size;         // this is for used hash indices
    chunksize += min_hash_size;        // this is for the hash begins
    chunksize += maxNumRoughNonzeros;  // this is for hash nexts
    chunksize += maxNumRoughNonzeros;  // this is for hash keys
    chunksize += maxNumRoughNonzeros;  // this is for hash values
  }

  if (current_spgemm_algorithm == SPGEMM_KK_DENSE && !exec_gpu) {
    nnz_lno_t col_size     = this->b_col_cnt / (sizeof(nnz_lno_t) * 8) + 1;
    nnz_lno_t max_row_size = KOKKOSKERNELS_MACRO_MIN(col_size, maxNumRoughNonzeros);
    chunksize              = col_size + max_row_size;
    // if speed is set, and exec space is cpu, then  we use dense accumulators.
    // or if memspeed is set, and concurrency is not high, we use dense
    // accumulators.
    maxNumRoughNonzeros = col_size;
    pool_init_val       = 0;
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tDense Acc - COLS:" << col_size << " max_row_size:" << max_row_size << std::endl;
    }
  }
  KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
  if (exec_gpu) {
    my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;
  }

  nnz_lno_t num_chunks = this->template compute_num_pool_chunks<pool_memory_space>(chunksize * sizeof(nnz_lno_t),
                                                                                   concurrency / suggested_vector_size);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tPool Size (MB):" << (num_chunks * chunksize * sizeof(nnz_lno_t)) / 1024. / 1024.
              << " num_chunks:" << num_chunks << " chunksize:" << chunksize << std::endl;
  }
  Kokkos::Timer timer1;
  pool_memory_space m_space(num_chunks, chunksize, pool_init_val, my_pool_type);
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tPool Alloc Time:" << timer1.seconds() << std::endl;
  }

  StructureC<a_r_view_t, a_nnz_view_t, b_original_row_view_t, b_compressed_row_view_t, b_nnz_view_t, c_row_view_t,
             /* nnz_lno_temp_work_view_t,*/ pool_memory_space>
      sc(m, row_mapA_, entriesA_, old_row_mapB, row_mapB_, entriesSetIndex, entriesSets, rowmapC, min_hash_size,
         maxNumRoughNonzeros, shmem_size_to_use, suggested_team_size, team_row_chunk_size, suggested_vector_size,
         m_space, lcl_my_exec_space, KOKKOSKERNELS_VERBOSE);

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tStructureC vector_size:" << suggested_vector_size << " team_size:" << suggested_team_size
              << " chunk_size:" << team_row_chunk_size << " shmem_size:" << shmem_size_to_use << std::endl;
  }

  timer1.reset();

  if (exec_gpu) {
    Kokkos::parallel_for("KokkosSparse::StructureC::GPU_EXEC",
                         gpu_team_policy_t(m / suggested_team_size + 1, suggested_team_size, suggested_vector_size),
                         sc);
  } else {
    if (current_spgemm_algorithm == SPGEMM_KK_DENSE) {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for("KokkosSparse::StructureC::SPGEMM_KK_DENSE::DYNAMIC",
                             dynamic_multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1,
                                                                         suggested_team_size, suggested_vector_size),
                             sc);
      } else {
        Kokkos::parallel_for("KokkosSparse::StructureC::SPGEMM_KK_DENSE::STATIC",
                             multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1, suggested_team_size,
                                                                 suggested_vector_size),
                             sc);
      }
    } else if (current_spgemm_algorithm == SPGEMM_KK_LP) {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC::SPGEMM_KK_LP::DYNAMIC",
            dynamic_multicore_team_policy4_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sc);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC::SPGEMM_KK_LP::STATIC",
            multicore_team_policy4_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), sc);
      }
    } else {
      if (use_dynamic_schedule) {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC::SPGEMM_KK_MEM::DYNAMIC",
            dynamic_multicore_team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
            sc);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::StructureC::SPGEMM_KK_MEM::STATIC",
            multicore_team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size), sc);
      }
    }
  }
  MyExecSpace().fence();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\tStructureC Kernel time:" << timer1.seconds() << std::endl << std::endl;
  }

#if 0
  int read_write_cost = this->handle->get_spgemm_handle()->get_read_write_cost_calc();

  if (read_write_cost && KOKKOSKERNELS_VERBOSE)
  {
	  std::cout << "\t\t!!!!DELETE THIS PART!!!! PRINTING STATS HERE!!!!!" << std::endl;
	  KokkosKernels::Impl::print_1Dview(old_row_mapB);
	  KokkosKernels::Impl::print_1Dview(row_mapB_);
	  KokkosKernels::Impl::kk_reduce_diff_view
	  	  	  <b_original_row_view_t, b_compressed_row_view_t, MyExecSpace>
	  	  	  (brows, old_row_mapB, row_mapB_, compressed_b_size);
	  std::cout << "\tcompressed_b_size:" << compressed_b_size << " bnnz:" << bnnz << std::endl;
	  std::cout << "Given compressed maxNumRoughNonzeros:" << maxNumRoughNonzeros << std::endl;

	  //row_lno_temp_work_view_t flops_per_row("flops_per_row", a_row_cnt);
	  nnz_lno_t r_maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, old_row_mapB,row_mapB_ /*,flops_per_row.data()*/);
	  std::cout << "compressed r_maxNumRoughZeros:" << r_maxNumRoughZeros << std::endl;

	  typename const_a_lno_row_view_t::HostMirror h_row_mapA =  Kokkos::create_mirror_view (row_mapA);
	  Kokkos::deep_copy(h_row_mapA, row_mapA);
	  typename const_a_lno_row_view_t::HostMirror h_entriesA =  Kokkos::create_mirror_view (entriesA);
	  Kokkos::deep_copy(h_entriesA, entriesA);
	  typename const_a_lno_row_view_t::HostMirror h_row_mapB_ =  Kokkos::create_mirror_view (row_mapB_);
	  Kokkos::deep_copy(h_row_mapB_, row_mapB_);
	  typename const_a_lno_row_view_t::HostMirror h_old_row_mapB =  Kokkos::create_mirror_view (old_row_mapB);
	  Kokkos::deep_copy(h_old_row_mapB, old_row_mapB);
	  typename const_a_lno_row_view_t::HostMirror h_row_mapB =  Kokkos::create_mirror_view (row_mapB);
	  Kokkos::deep_copy(h_row_mapB, row_mapB);
	  typename c_row_view_t::HostMirror h_rowmapC = Kokkos::create_mirror_view (rowmapC);
	  Kokkos::deep_copy(h_rowmapC, rowmapC);

	  MyExecSpace().fence();


	  std::unordered_map <size_t, nnz_lno_t> flop_to_row_count;
	  std::unordered_map <size_t, nnz_lno_t> compressed_flop_to_row_count;
	  std::unordered_map <size_t, nnz_lno_t> outputrowsize_to_row_count;
	  std::unordered_map <size_t, nnz_lno_t> output_floponrhs_to_row_count;
	  std::unordered_map <size_t, nnz_lno_t> output_rowcomp_to_row_count;

	  size_t compressed_flops = 0;
	  size_t original_flops = 0;
	  size_t compressd_max_flops= 0;
	  size_t original_max_flops = 0;
	  int GROUPSIZE = 64;

	  for (int i = 0; i < a_row_cnt; ++i){
		  size_type arb = h_row_mapA(i);
		  size_type are = h_row_mapA(i + 1);
		  size_t compressed_row_flops = 0;
		  size_t original_row_flops = 0;
		  for (size_type j = arb; j < are; ++j){
			  nnz_lno_t ae = h_entriesA(j);
			  compressed_row_flops += h_row_mapB_(ae) - h_old_row_mapB(ae);
			  original_row_flops += h_row_mapB(ae + 1) - h_row_mapB(ae);
		  }
		  if (compressed_row_flops > compressd_max_flops) compressd_max_flops = compressed_row_flops;
		  if (original_row_flops > original_max_flops) original_max_flops = original_row_flops;
		  compressed_flops += compressed_row_flops;
		  original_flops += original_row_flops;
		  flop_to_row_count[original_row_flops / GROUPSIZE]++;
		  compressed_flop_to_row_count[compressed_row_flops / GROUPSIZE]++;
		  size_t outputrowsize = h_rowmapC(i);

		  outputrowsize_to_row_count[outputrowsize]++;
		  output_floponrhs_to_row_count[original_row_flops / (are - arb) ]++;
		  output_rowcomp_to_row_count[original_row_flops / (outputrowsize) ]++;

	  }
	  std::cout   << "original_flops:" << original_flops
			  << " compressed_flops:" << compressed_flops
			  << " FLOP_REDUCTION:" << double(compressed_flops) / original_flops
			  << std::endl;
	  std::cout   << "original_max_flops:" << original_max_flops
			  << " compressd_max_flops:" << compressd_max_flops
			  << " MEM_REDUCTION:" << double(compressd_max_flops) / original_max_flops * 2
			  << std::endl;
	  std::cout   << "\tOriginal_B_SIZE:" << bnnz
			  << " Compressed_b_size:" << compressed_b_size
			  << std::endl;
	  std::cout << " AR AC ANNZ BR BC BNNZ original_flops compressed_flops FLOP_REDUCTION original_max_flops compressd_max_flops MEM_REDUCTION riginal_B_SIZE Compressed_b_size B_SIZE_REDUCTION" <<  std::endl;
	  std::cout << " " << a_row_cnt << " " << b_row_cnt << " "
			  	  	   << entriesA.extent(0) << " "
					   << b_row_cnt << " " << b_col_cnt << " "
					   << entriesB.extent(0) << " "
					   <<  original_flops << " "
					   << compressed_flops << " "
					   << double(compressed_flops) / original_flops
					   << " " << original_max_flops << " "
					   << compressd_max_flops << " "
					   << double(compressd_max_flops) / original_max_flops * 2 << " "
					   << bnnz << " " << compressed_b_size
					   <<" "<< double(compressed_b_size) / bnnz  << std::endl;


	  struct sortItem{
		  size_t key;
		  nnz_lno_t val;
		  sortItem(size_t k, nnz_lno_t v): key(k), val(v){}
		  bool operator<(const sortItem & a) const
		  {
		    //return !((this->src < a.src) || (this->src == a.src && this->dst < a.dst));
		    return (this->val > a.val);
		  }
	  };
	  {
		  std::vector<sortItem> flopvec;
		  // Iterate and print keys and values of unordered_map
		  for( const auto& n : flop_to_row_count ) {
			  flopvec.push_back(sortItem(n.first, n.second));
		  }
		  std::sort(flopvec.begin(), flopvec.begin() + flopvec.size());

		  for (size_t i = 0; i < flopvec.size(); ++i){
			  std::cout << "FLOPS:" << (flopvec[i].key) * GROUPSIZE << " - " << (flopvec[i].key + 1) * GROUPSIZE << " NumRow: " <<  flopvec[i].val << std::endl;
		  }
		  std::cout << std::endl << std::endl << std::endl << std::endl;
	  }
	  {
		  std::vector<sortItem> flopvec;
		  // Iterate and print keys and values of unordered_map
		  for( const auto& n : compressed_flop_to_row_count ) {
			  flopvec.push_back(sortItem(n.first, n.second));
		  }
		  std::sort(flopvec.begin(), flopvec.begin() + flopvec.size());

		  for (size_t i = 0; i < flopvec.size(); ++i){
			  std::cout << "COMPRESSED FLOPS:" << (flopvec[i].key) * GROUPSIZE << " - " << (flopvec[i].key  +1) * GROUPSIZE <<  " NumRow: " <<  flopvec[i].val << std::endl;
		  }
		  std::cout << std::endl << std::endl << std::endl << std::endl;
	  }

	  {
		  std::vector<sortItem> flopvec;
		  // Iterate and print keys and values of unordered_map
		  for( const auto& n : outputrowsize_to_row_count ) {
			  flopvec.push_back(sortItem(n.first, n.second));
		  }
		  std::sort(flopvec.begin(), flopvec.begin() + flopvec.size());

		  for (size_t i = 0; i < flopvec.size(); ++i){
			  std::cout << "OUTPUT ROWSIZE:" << (flopvec[i].key)  <<  " NumRow: " <<  flopvec[i].val << std::endl;
		  }
		  std::cout << std::endl << std::endl << std::endl << std::endl;
	  }
	  {
		  std::vector<sortItem> flopvec;
		  // Iterate and print keys and values of unordered_map
		  for( const auto& n : output_floponrhs_to_row_count ) {
			  flopvec.push_back(sortItem(n.first, n.second));
		  }
		  std::sort(flopvec.begin(), flopvec.begin() + flopvec.size());

		  for (size_t i = 0; i < flopvec.size(); ++i){
			  std::cout << "OUTPUT AVG RHS ROWFLOPS:" << (flopvec[i].key) <<  " NumRow: " <<  flopvec[i].val << std::endl;
		  }
		  std::cout << std::endl << std::endl << std::endl << std::endl;
	  }
	  {
		  std::vector<sortItem> flopvec;
		  // Iterate and print keys and values of unordered_map
		  for( const auto& n : output_rowcomp_to_row_count ) {
			  flopvec.push_back(sortItem(n.first, n.second));
		  }
		  std::sort(flopvec.begin(), flopvec.begin() + flopvec.size());

		  for (size_t i = 0; i < flopvec.size(); ++i){
			  std::cout << "OUTPUT AVG COMPRESSION:" << (flopvec[i].key) <<  " NumRow: " <<  flopvec[i].val << std::endl;
		  }
		  std::cout << std::endl << std::endl << std::endl << std::endl;
	  }
  }
#endif
  typename c_row_view_t::non_const_value_type c_nnz_size = 0;
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<MyExecSpace>(m + 1, rowmapC, c_nnz_size);
  this->handle->get_spgemm_handle()->set_c_nnz(c_nnz_size);
  nnz_lno_t c_max_nnz = KokkosSparse::Impl::graph_max_degree<MyExecSpace, size_type, c_row_view_t>(rowmapC);
  this->handle->get_spgemm_handle()->set_max_result_nnz(c_max_nnz);
}  // symbolic_c (end)

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_r_view_t, typename a_n_view_t, typename b_oldrow_view_t, typename b_r_view_t>

size_t KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::getMaxRoughRowNNZ(nnz_lno_t m, a_r_view_t row_mapA_,
                                                                                a_n_view_t entriesA_,

                                                                                b_oldrow_view_t row_pointers_begin_B,
                                                                                b_r_view_t row_pointers_end_B,
                                                                                size_type *flops_per_row /* = NULL*/) {
  // get the execution space type.
  // KokkosKernels::Impl::ExecSpaceType my_exec_space =
  // this->handle->get_handle_exec_space();
  int suggested_vector_size     = this->handle->get_suggested_vector_size(m, entriesA_.extent(0));
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency, m);

  PredicMaxRowNNZ<a_r_view_t, a_n_view_t, b_oldrow_view_t, b_r_view_t> pcnnnz(
      m, row_mapA_, entriesA_, row_pointers_begin_B, row_pointers_end_B, team_row_chunk_size, flops_per_row);

  typename b_oldrow_view_t::non_const_value_type rough_size = 0;
  Kokkos::parallel_reduce("KokkosSparse::PredicMaxRowNNZ::STATIC",
                          team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
                          pcnnnz, rough_size);
  MyExecSpace().fence();

  return rough_size;
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>

struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::PredicMaxRowNNZ_p {
  const nnz_lno_t m;          // num rows
  const size_type *row_mapA;  // row pointers of a
  const nnz_lno_t *entriesA;  // col
  const size_type *row_begins_B;
  const size_type *row_end_indices_B;
  const size_type min_val;
  const nnz_lno_t team_row_chunk_size;

  /**
   * \brief Constructor
   * \param m_: num rows in A.
   * \param row_mapA_: row pointers of A
   * \param entriesA_: col indices of A
   * \param row_begins_B_: row begin indices of B
   * \param row_end_indices_B_: row end indices of B
   * \param team_row_chunk_size_: the number of rows assigned to each team.
   */
  PredicMaxRowNNZ_p(const nnz_lno_t m_, const size_type *row_mapA_, const nnz_lno_t *entriesA_,

                    const size_type *row_begins_B_, const size_type *row_end_indices_B_, nnz_lno_t team_row_chunk_size_)
      : m(m_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        row_begins_B(row_begins_B_),
        row_end_indices_B(row_end_indices_B_),
        min_val(((std::numeric_limits<size_type>::lowest()))),
        team_row_chunk_size(team_row_chunk_size_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember, size_type &overal_max) const {
    // get the range of rows for team.
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
    const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);

    // TODO MD: here do I need a reduce as well?
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
                         [&](const nnz_lno_t &row_index) {
                           const size_type col_begin = row_mapA[row_index];
                           const size_type col_end   = row_mapA[row_index + 1];
                           const nnz_lno_t left_work = col_end - col_begin;

                           size_type max_num_results_in_row = 0;

                           // get the size of the rows of B, pointed by row of A
                           Kokkos::parallel_reduce(
                               Kokkos::ThreadVectorRange(teamMember, left_work),
                               [&](nnz_lno_t i, size_type &valueToUpdate) {
                                 const size_type adjind   = i + col_begin;
                                 const nnz_lno_t colIndex = entriesA[adjind];
                                 valueToUpdate += row_end_indices_B[colIndex] - row_begins_B[colIndex];
                               },
                               max_num_results_in_row);

                           // set max.
                           if (overal_max < max_num_results_in_row) {
                             overal_max = max_num_results_in_row;
                           }
                         });
  }

  KOKKOS_INLINE_FUNCTION
  void join(size_type &dst, const size_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_type &dst) const { dst = min_val; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
size_t KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::getMaxRoughRowNNZ_p(const nnz_lno_t m,
                                                                                  const size_type annz,
                                                                                  const size_type *row_mapA_,
                                                                                  const nnz_lno_t *entriesA_,

                                                                                  const size_type *row_pointers_begin_B,
                                                                                  const size_type *row_pointers_end_B) {
  int suggested_vector_size     = this->handle->get_suggested_vector_size(m, annz);
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency, m);

  PredicMaxRowNNZ_p pcnnnz(m, row_mapA_, entriesA_, row_pointers_begin_B, row_pointers_end_B, team_row_chunk_size);

  size_type rough_size = 0;
  Kokkos::parallel_reduce("KokkosSparse::PredicMaxRowNNZ_P::STATIC",
                          team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
                          pcnnnz, rough_size);
  MyExecSpace().fence();
  return rough_size;
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
size_t
KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
             b_scalar_nnz_view_t_>::getMaxRoughRowNNZIntersection_p(const nnz_lno_t m, const size_type annz,
                                                                    const size_type *row_mapA_,
                                                                    const nnz_lno_t *entriesA_,

                                                                    const size_type *row_pointers_begin_B,
                                                                    const size_type *row_pointers_end_B,
                                                                    nnz_lno_t *min_result_row_for_each_row) {
  // get the execution space type.
  // KokkosKernels::Impl::ExecSpaceType my_exec_space =
  // this->handle->get_handle_exec_space();
  int suggested_vector_size     = this->handle->get_suggested_vector_size(m, annz);
  int suggested_team_size       = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency, m);

  PredicMaxRowNNZIntersection pcnnnz(m, this->b_col_cnt, row_mapA_, entriesA_, row_pointers_begin_B, row_pointers_end_B,
                                     team_row_chunk_size, min_result_row_for_each_row);

  nnz_lno_t rough_size = 0;
  Kokkos::parallel_reduce("KokkosSparse::PredicMaxRowNNZIntersection::STATIC",
                          team_policy_t(m / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
                          pcnnnz, rough_size);
  MyExecSpace().fence();
  return rough_size;
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
          typename b_compressed_row_view_t, typename b_nnz_view_t, typename c_row_view_t,
          typename nnz_lno_temp_work_view_t, typename pool_memory_space>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::NonzeroesC {
  nnz_lno_t numrows;

  a_row_view_t row_mapA;
  a_nnz_view_t entriesA;

  b_original_row_view_t old_row_mapB;
  b_compressed_row_view_t row_mapB;
  b_nnz_view_t entriesSetIndicesB;
  b_nnz_view_t entriesSetsB;

  c_row_view_t rowmapC;
  nnz_lno_temp_work_view_t entriesSetIndicesC;
  nnz_lno_temp_work_view_t entriesSetsC;

  const nnz_lno_t pow2_hash_size;
  const nnz_lno_t pow2_hash_func;
  const nnz_lno_t MaxRoughNonZero;

  const size_t shared_memory_size;
  int vector_size;
  pool_memory_space m_space;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;

  /**
   * \brief Constructor.
   */

  NonzeroesC(nnz_lno_t m_, a_row_view_t row_mapA_, a_nnz_view_t entriesA_,

             b_original_row_view_t old_row_mapB_, b_compressed_row_view_t row_mapB_, b_nnz_view_t entriesSetIndicesB_,
             b_nnz_view_t entriesSetsB_,

             c_row_view_t rowmapC_, nnz_lno_temp_work_view_t entriesSetIndicesC_,

             const nnz_lno_t hash_size_, const nnz_lno_t MaxRoughNonZero_, const size_t sharedMemorySize_,
             const int vector_size_, pool_memory_space mpool_, const KokkosKernels::Impl::ExecSpaceType my_exec_space_)
      : numrows(m_),

        row_mapA(row_mapA_),
        entriesA(entriesA_),

        old_row_mapB(old_row_mapB_),
        row_mapB(row_mapB_),
        entriesSetIndicesB(entriesSetIndicesB_),
        entriesSetsB(entriesSetsB_),

        rowmapC(rowmapC_),
        entriesSetIndicesC(entriesSetIndicesC_),
        entriesSetsC(),

        pow2_hash_size(hash_size_),
        pow2_hash_func(hash_size_ - 1),
        MaxRoughNonZero(MaxRoughNonZero_),

        shared_memory_size(sharedMemorySize_),
        vector_size(vector_size_),
        m_space(mpool_),
        my_exec_space(my_exec_space_) {}

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
    nnz_lno_t row_index = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();

    using hashmapType =
        KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                        KokkosKernels::Experimental::HashOpType::bitwiseAnd>;
    if (row_index >= numrows) return;
    // get row index.

    nnz_lno_t *globally_used_hash_indices = NULL;
    nnz_lno_t globally_used_hash_count    = 0;
    nnz_lno_t used_hash_size              = 0;
    hashmapType hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    volatile nnz_lno_t *tmp = NULL;
    size_t tid              = get_thread_id(row_index);
    while (tmp == NULL) {
      tmp = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid));
    }

    globally_used_hash_indices = (nnz_lno_t *)tmp;
    tmp += pow2_hash_size;

    hm2.hash_begins = (nnz_lno_t *)(tmp);
    tmp += pow2_hash_size;

    // poins to the next elements
    hm2.hash_nexts = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;

    // holds the keys
    hm2.keys = (nnz_lno_t *)(tmp);
    tmp += MaxRoughNonZero;
    hm2.values = (nnz_lno_t *)(tmp);

    {
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t col_size  = row_mapA[row_index + 1] - col_begin;

      for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
        size_type a_col = colind + col_begin;

        nnz_lno_t rowB     = entriesA[a_col];
        size_type rowBegin = old_row_mapB(rowB);

        nnz_lno_t left_work = row_mapB(rowB) - rowBegin;

        for (nnz_lno_t i = 0; i < left_work; ++i) {
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_set_ind    = entriesSetIndicesB[adjind];
          nnz_lno_t b_set        = entriesSetsB[adjind];
          nnz_lno_t hash         = b_set_ind & pow2_hash_func;

          hm2.sequential_insert_into_hash_mergeOr_TrackHashes(b_set_ind, b_set, &used_hash_size,
                                                              &globally_used_hash_count, globally_used_hash_indices);
        }
      }

      int set_size     = sizeof(nnz_lno_t) * 8;
      nnz_lno_t num_el = rowmapC(row_index);
      for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii) {
        nnz_lno_t c_rows_setind = hm2.keys[ii];
        nnz_lno_t c_rows        = hm2.values[ii];

        int current_row = 0;
        nnz_lno_t unit  = 1;

        while (c_rows) {
          if (c_rows & unit) {
            // insert indices.
            entriesSetIndicesC(num_el++) = set_size * c_rows_setind + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit   = unit << 1;
        }
      }
      for (int i = 0; i < globally_used_hash_count; ++i) {
        nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      }
    }
    m_space.release_chunk(globally_used_hash_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag &, const team_member_t &teamMember) const {
    nnz_lno_t row_index = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank();

    using hashmapType =
        KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                        KokkosKernels::Experimental::HashOpType::bitwiseAnd>;

    if (row_index >= numrows) return;

    // printf("row:%d\n", row_index);

    int thread_memory       = ((shared_memory_size / 4 / teamMember.team_size())) * 4;
    char *all_shared_memory = (char *)(teamMember.team_shmem().get_shmem(shared_memory_size));

    // nnz_lno_t *alloc_global_memory = NULL;
    nnz_lno_t *globally_used_hash_indices = NULL;

    // shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    // used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *)(all_shared_memory);
    typedef typename std::remove_reference<decltype(*used_hash_sizes)>::type atomic_incr_type;

    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *)(all_shared_memory);

    all_shared_memory += sizeof(nnz_lno_t);
    int unit_memory                   = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
    nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory;

    nnz_lno_t *begins = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

    // poins to the next elements
    nnz_lno_t *nexts = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

    // holds the keys
    nnz_lno_t *keys = (nnz_lno_t *)(all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
    nnz_lno_t *vals = (nnz_lno_t *)(all_shared_memory);

    // printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts,
    // keys, vals); return; first level hashmap
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t, nnz_lno_t, nnz_lno_t,
                                                    KokkosKernels::Experimental::HashOpType::modulo>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

    hashmapType hm2(MaxRoughNonZero, pow2_hash_func, nullptr, nullptr, nullptr, nullptr);

    // initialize begins.
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
                         [&](int i) { begins[i] = -1; });

    // initialize hash usage sizes
    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      used_hash_sizes[0]          = 0;
      used_hash_sizes[1]          = 0;
      globally_used_hash_count[0] = 0;
    });

    bool is_global_alloced = false;

    const size_type col_end   = row_mapA[row_index + 1];
    const size_type col_begin = row_mapA[row_index];
    const nnz_lno_t col_size  = col_end - col_begin;

    for (nnz_lno_t colind = 0; colind < col_size; ++colind) {
      size_type a_col = colind + col_begin;

      nnz_lno_t rowB     = entriesA[a_col];
      size_type rowBegin = old_row_mapB(rowB);

      nnz_lno_t left_work = row_mapB(rowB) - rowBegin;

      while (left_work) {
        nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

        nnz_lno_t b_set_ind = -1, b_set = -1;
        nnz_lno_t hash = -1;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, work_to_handle), [&](nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          b_set_ind              = entriesSetIndicesB[adjind];
          b_set                  = entriesSetsB[adjind];
        });

        int num_unsuccess =
            hm.vector_atomic_insert_into_hash_mergeOr(b_set_ind, b_set, used_hash_sizes, shared_memory_hash_size);

        int overall_num_unsuccess = 0;

        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&](const int /*threadid*/, int &overall_num_unsuccess_) { overall_num_unsuccess_ += num_unsuccess; },
            overall_num_unsuccess);

        if (overall_num_unsuccess) {
          // printf("row:%d\n", row_index);
          if (!is_global_alloced) {
            volatile nnz_lno_t *tmp = NULL;
            size_t tid              = get_thread_id(row_index);
            while (tmp == NULL) {
              Kokkos::single(
                  Kokkos::PerThread(teamMember),
                  [&](volatile nnz_lno_t *&memptr) { memptr = (volatile nnz_lno_t *)(m_space.allocate_chunk(tid)); },
                  tmp);
            }
            is_global_alloced = true;

            globally_used_hash_indices = (nnz_lno_t *)tmp;
            tmp += pow2_hash_size;

            hm2.hash_begins = (nnz_lno_t *)(tmp);
            tmp += pow2_hash_size;

            // poins to the next elements
            hm2.hash_nexts = (nnz_lno_t *)(tmp);
            tmp += MaxRoughNonZero;

            // holds the keys
            hm2.keys = (nnz_lno_t *)(tmp);
            tmp += MaxRoughNonZero;
            hm2.values = (nnz_lno_t *)(tmp);
          }

          if (num_unsuccess) {
            // int insertion =
            hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                b_set_ind, b_set, used_hash_sizes + 1, globally_used_hash_count, globally_used_hash_indices);
          }
        }
        left_work -= work_to_handle;
        rowBegin += work_to_handle;
      }
    }

    Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
      if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
    });

    nnz_lno_t num_compressed_elements = used_hash_sizes[0];
    used_hash_sizes[0]                = 0;
    size_type row_begin               = rowmapC(row_index);
    int set_size                      = sizeof(nnz_lno_t) * 8;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, num_compressed_elements), [&](const nnz_lno_t ii) {
      nnz_lno_t c_rows_setind = hm.keys[ii];
      nnz_lno_t c_rows        = hm.values[ii];

      int current_row = 0;
      nnz_lno_t unit  = 1;

      while (c_rows) {
        if (c_rows & unit) {
          size_type wind                       = Kokkos::atomic_fetch_add(used_hash_sizes, atomic_incr_type(1));
          entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
        }
        current_row++;
        c_rows = c_rows & ~unit;
        unit   = unit << 1;
      }
    });

    if (is_global_alloced) {
      nnz_lno_t num_compressed_elements_ = used_hash_sizes[1];
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, num_compressed_elements_), [&](const nnz_lno_t ii) {
        nnz_lno_t c_rows_setind = hm2.keys[ii];
        nnz_lno_t c_rows        = hm2.values[ii];

        int current_row = 0;
        nnz_lno_t unit  = 1;

        while (c_rows) {
          if (c_rows & unit) {
            size_type wind                       = Kokkos::atomic_fetch_add(used_hash_sizes, atomic_incr_type(1));
            entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit   = unit << 1;
        }
      });

      // now thread leaves the memory as it finds. so there is no need to
      // initialize the hash begins
      nnz_lno_t dirty_hashes = globally_used_hash_count[0];
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, dirty_hashes), [&](nnz_lno_t i) {
        nnz_lno_t dirty_hash        = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      });

      Kokkos::single(Kokkos::PerThread(teamMember), [&]() { m_space.release_chunk(globally_used_hash_indices); });
    }
  }

  size_t team_shmem_size(int /*team_size*/) const { return shared_memory_size; }
};
}  // namespace Impl
}  // namespace KokkosSparse
