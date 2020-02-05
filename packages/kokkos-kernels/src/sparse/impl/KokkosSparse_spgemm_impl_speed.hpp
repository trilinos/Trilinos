/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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
*/

namespace KokkosSparse{

namespace Impl{

template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
          typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
          typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
          typename mpool_type>
struct KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
  NumericCMEM_CPU
{
  nnz_lno_t numrows;
  nnz_lno_t numcols;

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


  NumericCMEM_CPU(
      nnz_lno_t m_,
      nnz_lno_t k_,
      a_row_view_t row_mapA_,
      a_nnz_view_t entriesA_,
      a_scalar_view_t valuesA_,

      b_row_view_t row_mapB_,
      b_nnz_view_t entriesB_,
      b_scalar_view_t valuesB_,

      c_row_view_t rowmapC_,
      c_nnz_view_t entriesC_,
      c_scalar_view_t valuesC_,
      mpool_type memory_space_,
      const KokkosKernels::Impl::ExecSpaceType my_exec_space_,
      nnz_lno_t team_row_chunk_size):
        numrows(m_),
        numcols(k_),
        row_mapA (row_mapA_),
        entriesA(entriesA_),
        valuesA(valuesA_),

        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valuesB(valuesB_),

        rowmapC(rowmapC_),
        entriesC(entriesC_),
        valuesC(valuesC_),
        memory_space(memory_space_),
        pEntriesC(entriesC_.data()), pVals(valuesC.data()),
        my_exec_space(my_exec_space_),
        team_work_size(team_row_chunk_size){
        }


  KOKKOS_INLINE_FUNCTION
  size_t get_thread_id(const size_t row_index) const{
    switch (my_exec_space){
    default:
      return row_index;
#if defined( KOKKOS_ENABLE_SERIAL )
    case KokkosKernels::Impl::Exec_SERIAL:
      return 0;
#endif
#if defined( KOKKOS_ENABLE_OPENMP )
    case KokkosKernels::Impl::Exec_OMP:
  #ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      return Kokkos::OpenMP::hardware_thread_id();
  #else
      return Kokkos::OpenMP::impl_hardware_thread_id();
  #endif
#endif
#if defined( KOKKOS_ENABLE_THREADS )
    case KokkosKernels::Impl::Exec_PTHREADS:
  #ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      return Kokkos::Threads::hardware_thread_id();
  #else
      return Kokkos::Threads::impl_hardware_thread_id();
  #endif
#endif
#if defined( KOKKOS_ENABLE_QTHREAD)
    case KokkosKernels::Impl::Exec_QTHREADS:
      return 0; // Kokkos does not have a thread_id API for Qthreads
#endif
#if defined( KOKKOS_ENABLE_CUDA )
    case KokkosKernels::Impl::Exec_CUDA:
      return row_index;
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

    nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    scalar_t * dense_accum= NULL;
    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (dense_accum == NULL){
      dense_accum = (scalar_t * )( memory_space.allocate_chunk(tid));
    }
    char *marker = (char *) (dense_accum + numcols);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {

      const size_type c_row_begin = rowmapC[row_index];
      nnz_lno_t *myentries = pEntriesC + c_row_begin;
      scalar_t *myvals = pVals + c_row_begin;

      nnz_lno_t current_col_index = 0;
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t nnza = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

      for (nnz_lno_t colind = 0; colind < nnza; ++colind){
        size_type a_col = colind + col_begin;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;
        for (int i = 0; i < left_work; ++i){
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_col_ind = entriesB[adjind];
          scalar_t b_val = valuesB[adjind] * valA;
          if (marker[b_col_ind] == 0){
            marker[b_col_ind] = 1;
            myentries[current_col_index++] = b_col_ind;
          }
          dense_accum[b_col_ind] += b_val;
        }
      }
      for (nnz_lno_t i = 0; i < current_col_index; ++i){
        nnz_lno_t ind = myentries[i];
        myvals[i] = dense_accum[ind];
        dense_accum[ind] = 0;
        marker [ind] = 0;
      }
    });
    memory_space.release_chunk(dense_accum);
  }

};



template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__,
          typename b_row_view_t__, typename b_nnz_view_t__, typename b_scalar_view_t__,
          typename c_row_view_t__, typename c_nnz_view_t__, typename c_scalar_view_t__,
          typename c_nnz_tmp_view_t>

struct KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
  NumericCMEM
{
  nnz_lno_t numrows;

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

  const int unit_memory; //begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size;



  NumericCMEM(
      nnz_lno_t m_,
      a_row_view_t__ row_mapA_,
      a_nnz_view_t__ entriesA_,
      a_scalar_view_t__ valuesA_,

      b_row_view_t__ row_mapB_,
      b_nnz_view_t__ entriesB_,
      b_scalar_view_t__ valuesB_,

      c_row_view_t__ rowmapC_,
      c_nnz_view_t__ entriesC_,
      c_scalar_view_t__ valuesC_,

      c_nnz_tmp_view_t beginsC_,
      c_nnz_tmp_view_t nextsC_,

      const size_type sharedMemorySize_,
      const int suggested_vector_size,
      const nnz_lno_t team_row_chunk_size,
      int suggested_team_size_,
      bool KOKKOSKERNELS_VERBOSE_
      ):
        numrows(m_),
        row_mapA (row_mapA_),
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
        pbeginsC(beginsC_.data()), pnextsC(nextsC_.data()),
        pEntriesC(entriesC_.data()), pvaluesC(valuesC_.data()),
        shared_memory_size(sharedMemorySize_),

        vector_size (suggested_vector_size),
        team_work_size(team_row_chunk_size),

        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t)),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size /8 / suggested_team_size_) * 8),
        shmem_key_size(), shared_memory_hash_func(), shmem_hash_size(1)
        {
          shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 2) / unit_memory);
          if (KOKKOSKERNELS_VERBOSE_){
            std::cout << "\t\tNumericCMEM -- thread_memory:" << thread_memory  << " unit_memory:" << unit_memory <<
                " initial key size:" << shmem_key_size << std::endl;
          }
          while (shmem_hash_size * 2 <=  shmem_key_size){
            shmem_hash_size = shmem_hash_size * 2;
          }
          shared_memory_hash_func = shmem_hash_size - 1;

          shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
          shmem_key_size = (shmem_key_size >> 1) << 1;

// This guard will help ensure behavior is consistent within Trilinos
#ifdef KOKKOS_ENABLE_COMPLEX_ALIGN
          {
          // GPUTag
          // shmem allocation will be partitioned as below for hash map accumulator
          // thread_memory == 2*sizeof(nnz_lno_t) + shmem_hash_size*sizeof(nnz_lno_t) + 2*shmem_key_size*sizeof(nnz_lno_t) + rem_size*sizeof(scalar_t)

          // check that memory is partitioned into aligned chunks
          nnz_lno_t remainder_memory = thread_memory - sizeof(nnz_lno_t)*2 - shmem_hash_size*sizeof(nnz_lno_t);

          // The remainder of memory for vals must be aligned into sizeof(scalar_t) chunks, and there must be at least as many entries as keys
          nnz_lno_t val_memory = remainder_memory - 2*shmem_key_size*sizeof(nnz_lno_t);

          nnz_lno_t val_unalign_mem = val_memory % alignof(scalar_t);
          if (val_unalign_mem > 0) {
            // Redistributing between shmem_key_size and vals involves exchange of 2 "keys" (key+next) per val
            nnz_lno_t realign_chunk_mem = 2 * sizeof(nnz_lno_t);

            bool is_align_possible = (val_unalign_mem % realign_chunk_mem) == 0;
            if(!is_align_possible)
            {
              //throw std::runtime_error("NumericCMEM Ctor Error: unable to align memory for shared memory allocations. Modify your shared memory request");
              std::cout << "NumericCMEM Ctor WARNING: unable to align memory for shared memory allocations. Modify your shared memory request" << std::endl;
            }

            nnz_lno_t realign_chunks = val_unalign_mem / realign_chunk_mem; 

            shmem_key_size -= realign_chunks;
            val_memory = remainder_memory - 2*shmem_key_size*sizeof(nnz_lno_t);
            val_unalign_mem = val_memory%alignof(scalar_t);
          }

          if (val_unalign_mem > 0) {
            //throw std::runtime_error("NumericCMEM Ctor Error: shared memory realignment failed. Modify your shared memory request");
            std::cout << "NumericCMEM Ctor WARNING: shared memory realignment failed. Modify your shared memory request" << std::endl;
          }

          }
#endif

          if (KOKKOSKERNELS_VERBOSE_){
            std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << shmem_hash_size  << " shmem_key_size:" << shmem_key_size << std::endl;
          }
        }

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag&, const team_member_t & teamMember) const {


    //get the beginning and end rows of the team.
    nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);


    char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

    //shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    //used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

    //poins to the next elements
    nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;

    //holds the keys
    nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;
    scalar_t * vals = (scalar_t *) (all_shared_memory);


    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm(shmem_hash_size, shmem_key_size, begins, nexts, keys, vals);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(0, 0,
        NULL, NULL, NULL, NULL);
    /*
    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(global_memory_hash_size, global_memory_hash_size,
        pbeginsC + c_row_begin, pnextsC + c_row_begin, pEntriesC + c_row_begin, pvaluesC + c_row_begin);
        */


    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      const size_type c_row_begin = rowmapC[row_index];
      const nnz_lno_t global_memory_hash_size = nnz_lno_t(rowmapC[row_index + 1] - c_row_begin);

      hm2.hash_key_size = global_memory_hash_size;
      hm2.max_value_size = global_memory_hash_size;
      hm2.keys = pEntriesC + c_row_begin;
      hm2.values = pvaluesC + c_row_begin;
      hm2.hash_begins = pbeginsC + c_row_begin;
      hm2.hash_nexts = pnextsC + c_row_begin;

      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shmem_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);

      for (nnz_lno_t colind = 0; colind < left_work; ++colind){
        size_type a_col = colind + col_begin;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

        while (left_work_){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work_);
          nnz_lno_t b_col_ind = -1;
          scalar_t b_val = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_col_ind = entriesB[adjind];
            b_val = valuesB[adjind] * valA;
            //hash = b_col_ind % shmem_key_size;
            hash = b_col_ind & shared_memory_hash_func;
          });

          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(
              teamMember, vector_size,
              hash, b_col_ind, b_val,
              used_hash_sizes,
              shmem_key_size);

          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int /* threadid */, int &overall_num_unsuccess_) {
            overall_num_unsuccess_ += num_unsuccess;
          }, overall_num_unsuccess);

          if (overall_num_unsuccess){
            nnz_lno_t hash_ = -1;
            if (num_unsuccess) {
              hash_ = b_col_ind % global_memory_hash_size;
            }

            //int insertion =
            hm2.vector_atomic_insert_into_hash_mergeAdd(
                teamMember, vector_size,
                hash_,b_col_ind,b_val,
                used_hash_sizes + 1, hm2.max_value_size
            );

          }
          left_work_ -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
      });


      size_type num_elements = used_hash_sizes[0];


      size_type written_index = used_hash_sizes[1];
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, num_elements),
          [&] (size_type i) {
        pEntriesC[c_row_begin + written_index + i] = keys[i];
        pvaluesC[c_row_begin + written_index + i] = vals[i];
      });
    });
  }

  size_t team_shmem_size (int /* team_size */) const {
    return shared_memory_size;
  }
};


//
// * Notes on KokkosSPGEMM_numeric_speed *
//
// Prior to this routine, KokkosSPGEMM_numeric(...) was called
//
//   KokkosSPGEMM_numeric(...) :
//     if (this->spgemm_algorithm == SPGEMM_KK || SPGEMM_KK_LP == this->spgemm_algorithm) :
//       call KokkosSPGEMM_numeric_speed(...)
//     else:
//       call  KokkosSPGEMM_numeric_hash(...)
//
//
// KokkosSPGEMM_numeric_speed:
//
// Algorithm selection as follows and matching to kernel Tag:
//
//  Policy typedefs with tags found in: KokkosSparse_spgemm_impl.hpp
//
//  if Cuda enabled :
//    "KokkosSparse::NumericCMEM::KKSPEED::GPU" : gpu_team_policy_t,  i.e. GPUTag
//
//  else :
//    "KokkosSparse::NumericCMEM_CPU::DENSE::DYNAMIC" : dynamic_multicore_team_policy_t,  i.e. MultiCoreTag
//    "KokkosSparse::NumericCMEM_CPU::DENSE::STATIC" :  multicore_team_policy_t,  i.e. MultiCoreTag
//


template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void
  KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
  KokkosSPGEMM_numeric_speed(
    c_row_view_t rowmapC_,
    c_lno_nnz_view_t entriesC_,
    c_scalar_nnz_view_t valuesC_,
    KokkosKernels::Impl::ExecSpaceType my_exec_space_)
{

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\tSPEED MODE" << std::endl;
  }

  nnz_lno_t brows = row_mapB.extent(0) - 1;
  size_type bnnz =  valsB.extent(0);

  //get suggested vector size, teamsize and row chunk size.
  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

  Kokkos::Impl::Timer numeric_speed_timer_with_free;

  if (my_exec_space_ == KokkosKernels::Impl::Exec_CUDA){
    //allocate memory for begins and next to be used by the hashmap
    nnz_lno_temp_work_view_t beginsC
    (Kokkos::ViewAllocateWithoutInitializing("C keys"), valuesC_.extent(0));
    nnz_lno_temp_work_view_t nextsC
    (Kokkos::ViewAllocateWithoutInitializing("C nexts"), valuesC_.extent(0));
    Kokkos::deep_copy(beginsC, -1);

    //create the functor.
    NumericCMEM<
    const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t, nnz_lno_temp_work_view_t>
    sc(
        a_row_cnt,
        row_mapA,
        entriesA,
        valsA,

        row_mapB,
        entriesB,
        valsB,

        rowmapC_,
        entriesC_,
        valuesC_,

        beginsC, nextsC,
        shmem_size,
        suggested_vector_size,
        team_row_chunk_size,
        suggested_team_size,
        KOKKOSKERNELS_VERBOSE);

    Kokkos::Impl::Timer timer1;
    MyExecSpace().fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tGPU vector_size:" << suggested_vector_size
          <<  " team_size:" << suggested_team_size
          << " chunk_size:" << team_row_chunk_size
          << std::endl;
    }

    timer1.reset();
    //this is basically kkmem without memory pools.
    //only executed for to check the effect of memory pools.
    Kokkos::parallel_for( "KokkosSparse::NumericCMEM::KKSPEED::GPU",
        gpu_team_policy_t(
            a_row_cnt / team_row_chunk_size + 1 ,
            suggested_team_size ,
            suggested_vector_size),
            sc);
    MyExecSpace().fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
    }
  }
  else {

    Kokkos::Impl::Timer numeric_speed_timer;
    typedef KokkosKernels::Impl::UniformMemoryPool
        < MyTempMemorySpace, scalar_t> pool_memory_space;


    KokkosKernels::Impl::PoolType my_pool_type =
        KokkosKernels::Impl::OneThread2OneChunk;
    int num_chunks = concurrency;

    Kokkos::Impl::Timer timer1;
    pool_memory_space m_space
    (num_chunks, this->b_col_cnt + (this->b_col_cnt) / sizeof(scalar_t) + 1, 0,  my_pool_type);
    MyExecSpace().fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tPool Alloc Time:" << timer1.seconds() << std::endl;
      std::cout << "\tPool Size(MB):" <<
          sizeof(scalar_t) * (num_chunks *
              (this->b_col_cnt + (this->b_col_cnt) / sizeof(scalar_t) + 1))
              / 1024. / 1024.  << std::endl;
    }

    NumericCMEM_CPU<
    const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
    pool_memory_space>
    sc(
        a_row_cnt,
        b_col_cnt,
        row_mapA,
        entriesA,
        valsA,

        row_mapB,
        entriesB,
        valsB,

        rowmapC_,
        entriesC_,
        valuesC_,
        m_space,
        my_exec_space_,
        team_row_chunk_size);

    MyExecSpace().fence();
    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCPU vector_size:" << suggested_vector_size
          <<  " team_size:" << suggested_team_size
          << " chunk_size:" << team_row_chunk_size
          << std::endl;
    }
    timer1.reset();

    if (use_dynamic_schedule){
      Kokkos::parallel_for( "KokkosSparse::NumericCMEM_CPU::DENSE::DYNAMIC", dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    else {
      Kokkos::parallel_for( "KokkosSparse::NumericCMEM_CPU::DENSE::STATIC", multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }

    MyExecSpace().fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
      std::cout << "\t\tNumeric SPEED TIME:" << numeric_speed_timer.seconds() << std::endl;

    }
  }
  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tNumeric SPEED TIME WITH FREE:" << numeric_speed_timer_with_free.seconds() << std::endl;

  }
}
}
}

