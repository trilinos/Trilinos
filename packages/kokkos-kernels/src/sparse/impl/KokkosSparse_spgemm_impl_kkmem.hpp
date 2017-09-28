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
          typename pool_memory_type>
struct KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    PortableNumericCHASH{

  nnz_lno_t numrows;

  a_row_view_t row_mapA;
  a_nnz_view_t entriesA;
  a_scalar_view_t valuesA;

  b_row_view_t row_mapB;
  b_nnz_view_t entriesB;
  b_scalar_view_t valuesB;

  c_row_view_t rowmapC;
  c_nnz_view_t entriesC;
  c_scalar_view_t valuesC;


  nnz_lno_t *pEntriesC;
  scalar_t *pvaluesC;
  const size_t shared_memory_size;
  const int vector_size;
  pool_memory_type memory_space;

  //nnz_lno_t max_nnz;
  const nnz_lno_t pow2_hash_size;
  const nnz_lno_t max_nnz;
  const nnz_lno_t pow2_hash_func;
  const KokkosKernels::Impl::ExecSpaceType my_exec_space;
  const nnz_lno_t team_work_size;

  const int unit_memory; //begins, nexts, and keys. No need for vals yet.
  const int suggested_team_size;
  const int thread_memory;
  nnz_lno_t shmem_key_size;
  nnz_lno_t shared_memory_hash_func;
  nnz_lno_t shmem_hash_size;

  PortableNumericCHASH(
      nnz_lno_t m_,
      a_row_view_t row_mapA_,
      a_nnz_view_t entriesA_,
      a_scalar_view_t valuesA_,

      b_row_view_t row_mapB_,
      b_nnz_view_t entriesB_,
      b_scalar_view_t valuesB_,

      c_row_view_t rowmapC_,
      c_nnz_view_t entriesC_,
      c_scalar_view_t valuesC_,
      size_t shared_memory_size_,
      int vector_size_,
      pool_memory_type mpool_,
      nnz_lno_t min_hash_size, nnz_lno_t max_nnz_,
      int suggested_team_size_,
      const KokkosKernels::Impl::ExecSpaceType my_exec_space_,
      nnz_lno_t team_row_chunk_size,
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
        pEntriesC(entriesC_.ptr_on_device()), pvaluesC(valuesC_.ptr_on_device()),
        shared_memory_size(shared_memory_size_),
        vector_size (vector_size_),
        memory_space(mpool_),
        //max_nnz(),
        pow2_hash_size(min_hash_size), max_nnz(max_nnz_),
        pow2_hash_func(min_hash_size - 1),
        my_exec_space(my_exec_space_),
        team_work_size(team_row_chunk_size),

        unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t)),
        suggested_team_size(suggested_team_size_),
        thread_memory((shared_memory_size /8 / suggested_team_size_) * 8),
        shmem_key_size(), shared_memory_hash_func(), shmem_hash_size(1)
  {

    shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
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

    if (KOKKOSKERNELS_VERBOSE_){
      std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << shmem_hash_size  << " shmem_key_size:" << shmem_key_size << std::endl;
    }


  }
  KOKKOS_INLINE_FUNCTION
  size_t get_thread_id(const size_t row_index) const{
    switch (my_exec_space){
    default:
      return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
    case KokkosKernels::Impl::Exec_SERIAL:
      return 0;
#endif
#if defined( KOKKOS_HAVE_OPENMP )
    case KokkosKernels::Impl::Exec_OMP:
      return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
    case KokkosKernels::Impl::Exec_PTHREADS:
      return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
    case KokkosKernels::Impl::Exec_QTHREADS:
      return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
    case KokkosKernels::Impl::Exec_CUDA:
      return row_index;
#endif
    }

  }

  //assumes that the vector lane is 1, as in cpus
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(pow2_hash_size, pow2_hash_size,NULL, NULL, NULL, NULL);

    volatile nnz_lno_t * tmp = NULL;
    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL){
      tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
    }
    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
    tmp += pow2_hash_size ;

    hm2.hash_begins = (nnz_lno_t *) (tmp);
    tmp += pow2_hash_size;
    hm2.hash_nexts = (nnz_lno_t *) (tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_sizes = 0;

      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
      hm2.max_value_size = global_memory_hash_size;
      hm2.keys = pEntriesC + c_row_begin;
      hm2.values = pvaluesC + c_row_begin;

      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
      for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
        size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_workB = row_mapB(rowB + 1) - rowBegin;

        for ( nnz_lno_t i = 0; i < left_workB; ++i){
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_col_ind = entriesB[adjind];
          scalar_t b_val = valuesB[adjind] * valA;
          nnz_lno_t hash = b_col_ind & pow2_hash_func;

          //this has to be a success, we do not need to check for the success.
          //int insertion =
          hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
              hash, b_col_ind, b_val,
              &used_hash_sizes, hm2.max_value_size
              ,&globally_used_hash_count,
              globally_used_hash_indices
          );
        }
      }
      for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i){
        nnz_lno_t dirty_hash = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      }
    });
    memory_space.release_chunk(globally_used_hash_indices);
  }


  //assumes that the vector lane is 1, as in cpus
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag2&, const team_member_t & teamMember) const {

    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(pow2_hash_size, pow2_hash_size,NULL, NULL, NULL, NULL);

    volatile nnz_lno_t * tmp = NULL;
    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL){
      tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
    }
    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
    tmp += pow2_hash_size ;

    hm2.hash_begins = (nnz_lno_t *) (tmp);
    tmp += pow2_hash_size;
    hm2.hash_nexts = (nnz_lno_t *) (tmp);
    tmp += max_nnz;

    hm2.keys = (nnz_lno_t *) (tmp);
    tmp += max_nnz;
    hm2.values = (scalar_t *) (tmp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_sizes = 0;

      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
      hm2.max_value_size = global_memory_hash_size;

      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
      for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
        size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_workB = row_mapB(rowB + 1) - rowBegin;

        for ( nnz_lno_t i = 0; i < left_workB; ++i){
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_col_ind = entriesB[adjind];
          scalar_t b_val = valuesB[adjind] * valA;
          nnz_lno_t hash = b_col_ind & pow2_hash_func;

          //this has to be a success, we do not need to check for the success.
          //int insertion =
          hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
              hash, b_col_ind, b_val,
              &used_hash_sizes, hm2.max_value_size
              ,&globally_used_hash_count,
              globally_used_hash_indices
          );
        }
      }
      for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i){
        nnz_lno_t dirty_hash = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      }
      for (nnz_lno_t i = 0; i < global_memory_hash_size; ++i){
        pEntriesC [c_row_begin + i] = hm2.keys[i];
        pvaluesC [c_row_begin+i] =hm2.values[i];
      }

    });
    memory_space.release_chunk(globally_used_hash_indices);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {

    nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
    if (row_index >= numrows) return;

    nnz_lno_t globally_used_hash_count = 0;
    nnz_lno_t used_hash_sizes = 0;


    const size_type c_row_begin = rowmapC[row_index];
    const size_type c_row_end = rowmapC[row_index + 1];

    const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(pow2_hash_size, global_memory_hash_size,
          NULL, NULL, pEntriesC + c_row_begin, pvaluesC + c_row_begin);

    volatile nnz_lno_t * tmp = NULL;
    size_t tid = get_thread_id(row_index);
    while (tmp == NULL){
      tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
    }

    nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
    tmp += pow2_hash_size ;

    hm2.hash_begins = (nnz_lno_t *) (tmp);
    tmp += pow2_hash_size;
    hm2.hash_nexts = (nnz_lno_t *) (tmp);


    {
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
      for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
        size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

        for ( nnz_lno_t i = 0; i < left_work_; ++i){
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_col_ind = entriesB[adjind];
          scalar_t b_val = valuesB[adjind] * valA;
          nnz_lno_t hash = b_col_ind & pow2_hash_func;

          //this has to be a success, we do not need to check for the success.
          /*int insertion = */
          hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
              hash,b_col_ind,b_val,
              &used_hash_sizes, hm2.max_value_size
              ,&globally_used_hash_count, globally_used_hash_indices
          );
        }
      }
      for (nnz_lno_t i = 0; i < globally_used_hash_count; ++i){
        nnz_lno_t dirty_hash = globally_used_hash_indices[i];
        hm2.hash_begins[dirty_hash] = -1;
      }
    }
    memory_space.release_chunk(globally_used_hash_indices);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag&, const team_member_t & teamMember) const {

    nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);


    //int thread_memory = (shared_memory_size / 8 / teamMember.team_size()) * 8;
    char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

    //shift it to the thread private part
    all_shared_memory += thread_memory * teamMember.team_rank();

    //used_hash_sizes hold the size of 1st and 2nd level hashes
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    //int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t) ; //begins, nexts, keys and vals .
    //nnz_lno_t shmem_key_size = (thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory;
    //if (shmem_key_size & 1) shmem_key_size -= 1;

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
    hm2(pow2_hash_size, pow2_hash_size,
        NULL, NULL, NULL, NULL);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);
      hm2.max_value_size = global_memory_hash_size;
      hm2.keys = pEntriesC + c_row_begin;
      hm2.values = pvaluesC + c_row_begin;

      //initialize begins.
      Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, shmem_hash_size), [&] (nnz_lno_t i) {
          begins[i] = -1; });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });

      bool is_global_alloced = false;
      nnz_lno_t *globally_used_hash_indices = NULL;

      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
nnz_lno_t ii = left_work;
//for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
      while(ii--){
  size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];

        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;

        while (left_work_){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work_);
          nnz_lno_t b_col_ind = -1;
          scalar_t b_val = -1;
          nnz_lno_t hash = -1;

          Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, work_to_handle),
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
              [&] (const int threadid, int &overall_num_unsuccess_) {
            overall_num_unsuccess_ += num_unsuccess;
          }, overall_num_unsuccess);

          if (overall_num_unsuccess){

            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              //size_t tid = get_thread_id(row_index);
              //the code gets internal compiler error on gcc 4.7.2
              //assuming that this part only runs on GPUs for now, below fix
              //has the exact same behaviour and runs okay.
              size_t tid = row_index;

              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
                }, tmp);
              }


              is_global_alloced = true;
              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;
              hm2.hash_nexts = (nnz_lno_t *) (tmp);

            }

            nnz_lno_t hash_ = -1;
            if (num_unsuccess) {
              hash_ = b_col_ind & pow2_hash_func;
            }

            //this has to be a success, we do not need to check for the success.
            //int insertion =

        hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                teamMember, vector_size,
                hash_,b_col_ind,b_val,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
            );

          }

          left_work_ -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      if (is_global_alloced){

        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });

        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          memory_space.release_chunk(globally_used_hash_indices);
        });
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
      });

      nnz_lno_t num_elements = used_hash_sizes[0];

      nnz_lno_t written_index = used_hash_sizes[1];
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, num_elements),
          [&] (nnz_lno_t i) {
        pEntriesC[c_row_begin + written_index + i] = keys[i];
        pvaluesC[c_row_begin + written_index + i] = vals[i];
      });
    });
  }
  size_t team_shmem_size (int team_size) const {
    return shared_memory_size;
  }
};


template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void
  KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    KokkosSPGEMM_numeric_hash(
      c_row_view_t rowmapC_,
      c_lno_nnz_view_t entriesC_,
      c_scalar_nnz_view_t valuesC_,
      KokkosKernels::Impl::ExecSpaceType my_exec_space){

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\tHASH MODE" << std::endl;
  }

  nnz_lno_t brows = row_mapB.dimension_0() - 1;
  size_type bnnz =  valsB.dimension_0();

  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

  typedef KokkosKernels::Impl::UniformMemoryPool< MyTempMemorySpace, nnz_lno_t> pool_memory_space;


  nnz_lno_t max_nnz = this->handle->get_spgemm_handle()->get_max_result_nnz();
  nnz_lno_t min_hash_size = 1;
  while (max_nnz > min_hash_size){
    min_hash_size *= 4;
  }
#ifdef KOKKOSKERNELSCHANGEPARAMS
  nnz_lno_t env_hash = atoi(getenv("MINHASHSIZE"));
  nnz_lno_t env_chunksize = atoi(getenv("CHUNKSIZE"));
  nnz_lno_t env_num_chunks = atoi(getenv("NUMCHUNKS"));
#endif


  size_t chunksize = min_hash_size; //this is for used hash indices
  chunksize += min_hash_size ; //this is for the hash begins
  chunksize += max_nnz; //this is for hash nexts
  int num_chunks = concurrency / suggested_vector_size;

#ifdef KOKKOSKERNELSCHANGEPARAMS

  if (env_hash > 2) {
min_hash_size = env_hash;
  }
  if (env_chunksize > 2){
chunksize = env_chunksize;
  }
  if (env_num_chunks > 2){
num_chunks = env_num_chunks;
  }
#endif
  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\t max_nnz: " << max_nnz
              << " chunk_size:" << chunksize
              << " min_hash_size:" << min_hash_size
              << " concurrency:" << concurrency
              << " MyExecSpace::concurrency():" << MyExecSpace::concurrency()
              << " numchunks:" << num_chunks << std::endl;
  }

  KokkosKernels::Impl::PoolType my_pool_type =
      KokkosKernels::Impl::OneThread2OneChunk;
  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
    my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;
  }

  Kokkos::Impl::Timer timer1;
  pool_memory_space m_space(num_chunks, chunksize, -1,  my_pool_type);
  MyExecSpace::fence();

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tPool Alloc Time:" << timer1.seconds() << std::endl;
    std::cout << "\t\tPool Size(MB):" <<
        sizeof (nnz_lno_t) * (num_chunks * chunksize) / 1024. / 1024.  << std::endl;
    m_space.print_memory_pool();
  }

  PortableNumericCHASH<
    const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
    pool_memory_space>
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
      shmem_size,
      suggested_vector_size,
      m_space,
      min_hash_size, max_nnz,
      suggested_team_size,

      my_exec_space,
      team_row_chunk_size,KOKKOSKERNELS_VERBOSE);


  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tvector_size:" << suggested_vector_size  << " chunk_size:" << team_row_chunk_size << std::endl;
  }
  timer1.reset();

  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
    //TODO CLEAN GPU CODE
    Kokkos::parallel_for( gpu_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    MyExecSpace::fence();
  }
  else {
    if (use_dynamic_schedule){

      Kokkos::parallel_for( dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    else {

      Kokkos::parallel_for( multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    MyExecSpace::fence();
  }

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
  }

}

template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void
  KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    KokkosSPGEMM_numeric_hash2(
      c_row_view_t rowmapC_,
      c_lno_nnz_view_t entriesC_,
      c_scalar_nnz_view_t valuesC_,
      KokkosKernels::Impl::ExecSpaceType my_exec_space){
  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\tHASH MODE" << std::endl;
  }

  nnz_lno_t brows = row_mapB.dimension_0() - 1;
  size_type bnnz =  valsB.dimension_0();

  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

  typedef KokkosKernels::Impl::UniformMemoryPool< MyTempMemorySpace, nnz_lno_t> pool_memory_space;


  nnz_lno_t max_nnz = this->handle->get_spgemm_handle()->get_max_result_nnz();
  nnz_lno_t min_hash_size = 1;
  while (max_nnz > min_hash_size){
    min_hash_size *= 4;
  }

  size_t chunksize = min_hash_size; //this is for used hash indices
  chunksize += min_hash_size ; //this is for the hash begins
  chunksize += max_nnz; //this is for hash nexts
  chunksize += max_nnz; //this is for indices
  chunksize += max_nnz * (sizeof (scalar_t)/ sizeof (nnz_lno_t)); //this is for values
  int num_chunks = concurrency / suggested_vector_size;

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\t max_nnz: " << max_nnz
              << " chunk_size:" << chunksize
              << " numchunks:" << num_chunks << std::endl;
  }

  KokkosKernels::Impl::PoolType my_pool_type =
      KokkosKernels::Impl::OneThread2OneChunk;
  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
    my_pool_type = KokkosKernels::Impl::ManyThread2OneChunk;
  }

  Kokkos::Impl::Timer timer1;
  pool_memory_space m_space(num_chunks, chunksize, -1,  my_pool_type);
  MyExecSpace::fence();

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tPool Alloc Time:" << timer1.seconds() << std::endl;
    std::cout << "\t\tPool Size(MB):" <<
        sizeof (nnz_lno_t) * (num_chunks * chunksize) / 1024. / 1024.  << std::endl;
  }

  PortableNumericCHASH<
    const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
    const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
    c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t,
    pool_memory_space>
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
      shmem_size,
      suggested_vector_size,
      m_space,
      min_hash_size, max_nnz,
      suggested_team_size,

      my_exec_space,
      team_row_chunk_size,KOKKOSKERNELS_VERBOSE);


  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tvector_size:" << suggested_vector_size  << " chunk_size:" << team_row_chunk_size << std::endl;
  }
  timer1.reset();

  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
    //TODO CLEAN GPU CODE
    Kokkos::parallel_for( gpu_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    MyExecSpace::fence();
  }
  else {
    if (use_dynamic_schedule){

      Kokkos::parallel_for( dynamic_multicore_team_policy2_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    else {

      Kokkos::parallel_for( multicore_team_policy2_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    MyExecSpace::fence();
  }

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
  }

}

}
}
