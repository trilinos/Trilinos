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
#define HASHSCALAR 107
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
  nnz_lno_t thread_shmem_key_size;
  nnz_lno_t thread_shared_memory_hash_func;
  nnz_lno_t thread_shmem_hash_size;


  nnz_lno_t team_shmem_key_size;
  nnz_lno_t team_shared_memory_hash_func;
  nnz_lno_t team_shmem_hash_size;

  nnz_lno_t team_cuckoo_key_size, team_cuckoo_hash_func;

  nnz_lno_t max_first_level_hash_size;
  row_lno_persistent_work_view_t flops_per_row;
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
      nnz_lno_t team_row_chunk_size, double first_level_cut_off,
	  row_lno_persistent_work_view_t flops_per_row_,
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
        thread_shmem_key_size(),
		thread_shared_memory_hash_func(),
		thread_shmem_hash_size(1),
		team_shmem_key_size(),
		team_shared_memory_hash_func(),
		team_shmem_hash_size(1),
		team_cuckoo_key_size (1),
		team_cuckoo_hash_func(1), max_first_level_hash_size( 1), flops_per_row(flops_per_row_)


  {
	nnz_lno_t tmp_team_cuckoo_key_size = ((shared_memory_size - sizeof(nnz_lno_t) * 2) / (sizeof(nnz_lno_t) + sizeof(scalar_t )));

	while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
	team_cuckoo_hash_func = team_cuckoo_key_size - 1;
	team_shmem_key_size = ((shared_memory_size - sizeof(nnz_lno_t) * 4) / unit_memory);
    thread_shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
    if (KOKKOSKERNELS_VERBOSE_& 0 ){
      std::cout << "\t\tNumericCMEM -- thread_memory:" << thread_memory  << " unit_memory:" << unit_memory <<" initial key size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tNumericCMEM -- team_memory:" << shared_memory_size  << " unit_memory:" << unit_memory <<" initial team key size:" << team_shmem_key_size << std::endl;
    }
    while (thread_shmem_hash_size * 2 <=  thread_shmem_key_size){
      thread_shmem_hash_size = thread_shmem_hash_size * 2;
    }
    while (team_shmem_hash_size * 2 <=  team_shmem_key_size){
    	team_shmem_hash_size = team_shmem_hash_size * 2;
    }
    team_shared_memory_hash_func = team_shmem_hash_size - 1;
    thread_shared_memory_hash_func = thread_shmem_hash_size - 1;
    team_shmem_key_size = team_shmem_key_size + ((team_shmem_key_size - team_shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
    team_shmem_key_size = (team_shmem_key_size >> 1) << 1;

    thread_shmem_key_size = thread_shmem_key_size + ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
    thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;
    max_first_level_hash_size = first_level_cut_off * team_cuckoo_key_size;
    if (KOKKOSKERNELS_VERBOSE_){
      std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << thread_shmem_hash_size  << " shmem_key_size:" << thread_shmem_key_size << std::endl;
      std::cout << "\t\tNumericCMEM -- adjusted team hashsize:" << team_shmem_hash_size  << " team_shmem_key_size:" << team_shmem_key_size << std::endl;
      std::cout << "\t\tteam_cuckoo_key_size:" << team_cuckoo_key_size << " team_cuckoo_hash_func:" << team_cuckoo_hash_func << " max_first_level_hash_size:" << max_first_level_hash_size << std::endl;
      std::cout << "\t\tpow2_hash_size:" << pow2_hash_size << " pow2_hash_func:" << pow2_hash_func << std::endl;
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
#if defined( KOKKOS_ENABLE_CUDA )
    case KokkosKernels::Impl::Exec_CUDA:
      return row_index;
#endif
    }
  }

  //linear probing with tracking.
  KOKKOS_INLINE_FUNCTION
  void operator()(const MultiCoreTag4&, const team_member_t & teamMember) const {

    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);


    volatile nnz_lno_t * tmp = NULL;
    size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
    while (tmp == NULL){
      tmp = (volatile nnz_lno_t * )( memory_space.allocate_chunk(tid));
    }

    nnz_lno_t *used_indices = (nnz_lno_t *) (tmp);
    tmp += max_nnz;
    nnz_lno_t *hash_ids = (nnz_lno_t *) (tmp);
    tmp += pow2_hash_size;
    scalar_t *hash_values = (scalar_t *) (tmp);


    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      nnz_lno_t used_count = 0;

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
          nnz_lno_t hash = (b_col_ind * HASHSCALAR) & pow2_hash_func;

          while (true){
            if (hash_ids[hash] == -1){
            	used_indices[used_count++] = hash;
            	hash_ids[hash] = b_col_ind;
            	hash_values[hash] = b_val;
            	break;
            }
            else if (hash_ids[hash] == b_col_ind){
            	hash_values[hash] += b_val;
            	break;
            }
            else {
            	hash = (hash + 1) & pow2_hash_func;
            }
          }
        }
      }
      size_type c_row_begin = rowmapC[row_index];
      for (nnz_lno_t i = 0; i < used_count; ++i){
    	  nnz_lno_t used_index = used_indices[i];
    	  pEntriesC[c_row_begin] = hash_ids[used_index];
    	  pvaluesC [c_row_begin++] = hash_values[used_index];
    	  hash_ids[used_index] = -1;
      }
    });
    memory_space.release_chunk(used_indices);
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
          //nnz_lno_t hash = (b_col_ind * 107) & pow2_hash_func;
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
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_hash_size;

    //poins to the next elements
    nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_key_size;

    //holds the keys
    nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * thread_shmem_key_size;
    scalar_t * vals = (scalar_t *) (all_shared_memory);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm(thread_shmem_hash_size, thread_shmem_key_size, begins, nexts, keys, vals);

    KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
    hm2(pow2_hash_size, pow2_hash_size,
        NULL, NULL, NULL, NULL);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index) {
      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];
      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);

      bool is_global_alloced = false;
      nnz_lno_t *globally_used_hash_indices = NULL;


	  if (global_memory_hash_size > thread_shmem_key_size){
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
      hm2.max_value_size = global_memory_hash_size;
      hm2.keys = pEntriesC + c_row_begin;
      hm2.values = pvaluesC + c_row_begin;

      //initialize begins.
      Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, thread_shmem_hash_size), [&] (nnz_lno_t i) {
          begins[i] = -1; });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });


      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t left_work = row_mapA[row_index + 1] - col_begin;
      nnz_lno_t ii = left_work;
//for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
      while(ii-- > 0){
    	size_type a_col = col_begin + ii;
        nnz_lno_t rowB = entriesA[a_col];
        scalar_t valA = valuesA[a_col];


        size_type rowBegin = row_mapB(rowB);
        nnz_lno_t left_work_ = row_mapB(rowB + 1) - rowBegin;
        Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, left_work_),
            [&] (nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          nnz_lno_t b_col_ind = entriesB[adjind];
          scalar_t b_val = valuesB[adjind] * valA;
          //hash = b_col_ind % shmem_key_size;
          nnz_lno_t hash = b_col_ind & thread_shared_memory_hash_func;
          volatile int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeAdd(
              teamMember, vector_size,
              hash, b_col_ind, b_val,
              used_hash_sizes,
              thread_shmem_key_size);
          if (num_unsuccess){

        	  hash = b_col_ind & pow2_hash_func;
        	  hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
        			  teamMember, vector_size,
					  hash,b_col_ind,b_val,
					  used_hash_sizes + 1, hm2.max_value_size
					  ,globally_used_hash_count, globally_used_hash_indices
        	  );


          }

        });
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
        if (used_hash_sizes[0] > thread_shmem_key_size) used_hash_sizes[0] = thread_shmem_key_size;
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

  //one row does not fit into shmem, with thread-flat-parallel
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag6&, const team_member_t & teamMember) const {

    nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);
    char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

    const nnz_lno_t init_value = -1;
    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;
    //holds the keys
    nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * team_cuckoo_key_size;
    scalar_t * vals = (scalar_t *) (all_shared_memory);

    int thread_rank =  teamMember.team_rank();

    int vector_rank = 0;
    typedef typename std::remove_reference< decltype( *used_hash_sizes ) >::type atomic_incr_type;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
      	vector_rank = update;
      }
      update += 1;
    });
    int bs = vector_size * suggested_team_size;
    int vector_shift = thread_rank * vector_size + vector_rank;

    for (nnz_lno_t row_index = team_row_begin; row_index < team_row_end; ++row_index){

#if 1
        teamMember.team_barrier();
#endif
      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];
      const nnz_lno_t c_row_size = c_row_end -  c_row_begin;
      nnz_lno_t *c_row = entriesC.data() + c_row_begin;
      scalar_t *c_row_vals = valuesC.data() + c_row_begin;
      nnz_lno_t *global_acc_row_keys = c_row;
      scalar_t *global_acc_row_vals = c_row_vals;
	  volatile nnz_lno_t * tmp = NULL;

      if (c_row_size > max_first_level_hash_size){
    	  {
    		  while (tmp == NULL){
    			  Kokkos::single(Kokkos::PerTeam(teamMember),[=] (volatile nnz_lno_t * &memptr) {
    				  memptr = (volatile nnz_lno_t * )( memory_space.allocate_chunk(row_index));
    			  }, tmp);
    		  }
    		  global_acc_row_keys = (nnz_lno_t *) (tmp);
    		  global_acc_row_vals = (scalar_t *) (tmp + pow2_hash_size);
    	  }
          //initialize begins.
          {
        	  nnz_lno_t num_threads =  pow2_hash_size / vector_size;
        	  // not needed as team_cuckoo_key_size is always pow2. + (team_cuckoo_key_size & (vector_size - 1)) * 1;
        	  Kokkos::parallel_for( Kokkos::TeamThreadRange(teamMember, num_threads), [&] (nnz_lno_t teamind) {
        		  Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, vector_size ), [&] (nnz_lno_t i) {
        			  global_acc_row_vals[teamind * vector_size + i] = 0;
        		  });
        	  });
          }
      }

      //initialize begins.
      {
    	  nnz_lno_t num_threads =  team_cuckoo_key_size / vector_size;
    	  // not needed as team_cuckoo_key_size is always pow2. + (team_cuckoo_key_size & (vector_size - 1)) * 1;
    	  Kokkos::parallel_for( Kokkos::TeamThreadRange(teamMember, num_threads), [&] (nnz_lno_t teamind) {
    		  Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, vector_size ), [&] (nnz_lno_t i) {
    			  keys[teamind * vector_size + i] = init_value; vals[teamind * vector_size + i] = 0;
    		  });
    	  });
      }

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerTeam(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      bool insert_is_on = true;
      const size_type a_col_begin_offset = row_mapA[row_index];


      nnz_lno_t a_col_ind = entriesA[a_col_begin_offset];
      scalar_t a_col_val = valuesA[a_col_begin_offset];

      nnz_lno_t current_a_column_offset_inrow = 0;
      nnz_lno_t flops_on_the_left_of_offsett = 0;
      size_type current_b_read_offsett = row_mapB[a_col_ind];
      nnz_lno_t current_a_column_flops = row_mapB[a_col_ind + 1] - current_b_read_offsett;

      nnz_lno_t row_flops = flops_per_row(row_index);

#if 1
        teamMember.team_barrier();
#endif
        for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift< row_flops; vector_read_shift += bs){
        {
      	  nnz_lno_t my_b_col_shift = vector_read_shift - flops_on_the_left_of_offsett;
    	  nnz_lno_t my_b_col = init_value; scalar_t my_b_val = 0; nnz_lno_t hash = init_value;
    	  int fail = 0;

			  if (my_b_col_shift >= current_a_column_flops ){
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
			  //now insert it to first level hashmap accumulator.
			  hash = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
			  fail = 1;
			  bool try_to_insert = true;

			  //nnz_lno_t max_tries = team_cuckoo_key_size;
			  nnz_lno_t search_end = team_cuckoo_key_size; //KOKKOSKERNELS_MACRO_MIN(team_cuckoo_key_size, hash + max_tries);
			  for (nnz_lno_t trial = hash; trial < search_end; ){
				  if (keys[trial] == my_b_col){
					  Kokkos::atomic_add(vals + trial, my_b_val);
					  fail = 0;
					  break;
				  }
				  else if (keys[trial] == init_value){
					  if (!insert_is_on) {
						  try_to_insert = false;
						  break;
					  }
					  else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)){
						  Kokkos::atomic_add(vals + trial, my_b_val);
						  Kokkos::atomic_increment(used_hash_sizes);
						  if (used_hash_sizes[0] > max_first_level_hash_size)  insert_is_on = false;
						  fail = 0;
						  break;
					  }
				  }
				  else {
					  ++trial;
				  }
			  }
			  if (fail){
				  search_end = hash; //max_tries - (team_cuckoo_key_size -  hash);

				  for (nnz_lno_t trial = 0; try_to_insert && trial < search_end; ){
					  if (keys[trial] == my_b_col){
						  Kokkos::atomic_add(vals + trial, my_b_val);
						  fail = 0;
						  break;
					  }
					  else if (keys[trial] == init_value){
						  if (!insert_is_on) {
							  break;
						  }
						  else if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)){
							  Kokkos::atomic_add(vals + trial, my_b_val);
							  Kokkos::atomic_increment(used_hash_sizes);
							  if (used_hash_sizes[0] > max_first_level_hash_size)  insert_is_on = false;
							  fail = 0;
							  break;
						  }
					  }
					  else {
						  ++trial;
					  }
				  }


				  if (fail) {
				      nnz_lno_t new_hash = (my_b_col * HASHSCALAR) & pow2_hash_func;

					  for (nnz_lno_t trial = new_hash; trial < pow2_hash_size; ){
						  if (global_acc_row_keys[trial] == my_b_col){
							  Kokkos::atomic_add(global_acc_row_vals + trial , my_b_val);

							  //c_row_vals[trial] += my_b_val;
							  fail = 0;
							  break;
						  }
						  else if (global_acc_row_keys[trial ] == init_value){
							  if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial , init_value, my_b_col)){
								  Kokkos::atomic_add(global_acc_row_vals + trial , my_b_val);
								  //Kokkos::atomic_increment(used_hash_sizes + 1);
								  //c_row_vals[trial] = my_b_val;
								  fail = 0;
								  break;
							  }
						  }
						  else {
							  ++trial;
						  }
					  }
					  if (fail){
						  for (nnz_lno_t trial = 0; trial < new_hash; ){
							  if (global_acc_row_keys[trial ] == my_b_col){
								  //c_row_vals[trial] += my_b_val;
								  Kokkos::atomic_add(global_acc_row_vals + trial , my_b_val);

								  break;
							  }
							  else if (global_acc_row_keys[trial ] == init_value){
								  if (Kokkos::atomic_compare_exchange_strong(global_acc_row_keys + trial , init_value, my_b_col)){
									  //Kokkos::atomic_increment(used_hash_sizes + 1);
									  Kokkos::atomic_add(global_acc_row_vals + trial , my_b_val);
									  //c_row_vals[trial] = my_b_val;
									  break;
								  }
							  }
							  else {
								  ++trial;
							  }
						  }
					  }
				  }
			  }
		  }
      }

      teamMember.team_barrier();


      if (tmp != NULL){


    	  for (nnz_lno_t  my_index = vector_shift; my_index < pow2_hash_size;my_index += bs){
    		  nnz_lno_t my_b_col = global_acc_row_keys[my_index];
    		  if (my_b_col != init_value){
    			  scalar_t my_b_val = global_acc_row_vals[my_index];
				  int fail = 1;
    			  {
					  nnz_lno_t hash = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;

    				  //nnz_lno_t max_tries = team_cuckoo_key_size;
    				  nnz_lno_t search_end = team_cuckoo_key_size; //KOKKOSKERNELS_MACRO_MIN(team_cuckoo_key_size, hash + max_tries);
    				  for (nnz_lno_t trial = hash; trial < search_end;++trial){
    					  if (keys[trial] == my_b_col){
    						  vals[trial] += my_b_val;
    						  fail = 0;
    						  break;
    					  }
    					  else if (keys[trial] == init_value){
    						  break;
    					  }
    				  }
    				  search_end = hash; //max_tries - (team_cuckoo_key_size -  hash);

    				  for (nnz_lno_t trial = 0; trial < search_end; ++trial){
    					  if (keys[trial] == my_b_col){
    						  vals[trial] += my_b_val;
    						  fail = 0;
    						  break;
    					  }
    					  else if (keys[trial] == init_value){
    						  break;
    					  }
    				  }

    			  }
    			  if (fail){
    				  nnz_lno_t write_index = 0;
    				  write_index = Kokkos::atomic_fetch_add(used_hash_sizes + 1, atomic_incr_type(1));
    				  c_row[write_index] = my_b_col;
    				  c_row_vals[write_index] = my_b_val;
    			  }
    			  global_acc_row_keys[my_index] = init_value;
    		  }
    	  }

          teamMember.team_barrier();
		  Kokkos::single(Kokkos::PerTeam(teamMember),[&] () {
			  memory_space.release_chunk(global_acc_row_keys);
		  });
      }


      for (nnz_lno_t  my_index = vector_shift; my_index < team_cuckoo_key_size;my_index += bs){
    	  nnz_lno_t my_key = keys[my_index];
    	  if (my_key != init_value){
    		  scalar_t my_val = vals[my_index];
    		  nnz_lno_t write_index = 0;
    		  write_index = Kokkos::atomic_fetch_add(used_hash_sizes + 1, atomic_incr_type(1));
    		  c_row[write_index] = my_key;
    		  c_row_vals[write_index] = my_val;
    	  }
      }


    }
  }


  //In this one row fits into shmem with team-flat-parallel
  KOKKOS_INLINE_FUNCTION
  void operator()(const GPUTag4&, const team_member_t & teamMember) const {


	const nnz_lno_t init_value = -1;
    nnz_lno_t team_row_begin = teamMember.league_rank()  * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

    char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

    volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * 2;

    //holds the keys
    nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
    all_shared_memory += sizeof(nnz_lno_t) * team_cuckoo_key_size;
    scalar_t * vals = (scalar_t *) (all_shared_memory);

    int thread_rank =  teamMember.team_rank();

    int vector_rank = 0;
    typedef typename std::remove_reference< decltype( *used_hash_sizes ) >::type atomic_incr_type;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
      	vector_rank = update;
      }
      update += 1;
    });
    int bs = vector_size * suggested_team_size;
    int vector_shift = thread_rank * vector_size + vector_rank;
    for (nnz_lno_t row_index = team_row_begin; row_index < team_row_end; ++row_index){

#if 1
        teamMember.team_barrier();
#endif
      const size_type c_row_begin = rowmapC[row_index];
      //const size_type c_row_end = rowmapC[row_index + 1];
      //const nnz_lno_t c_row_size = c_row_end -  c_row_begin;
      nnz_lno_t *c_row = entriesC.data() + c_row_begin;
      scalar_t *c_row_vals = valuesC.data() + c_row_begin;

      //initialize begins.
      {
    	  nnz_lno_t num_threads =  team_cuckoo_key_size / vector_size;// not needed as team_cuckoo_key_size is always pow2. + (team_cuckoo_key_size & (vector_size - 1)) * 1;
    	  Kokkos::parallel_for( Kokkos::TeamThreadRange(teamMember, num_threads), [&] (nnz_lno_t teamind) {
    		  //nnz_lno_t team_shift = teamind * vector_size;
    		  //nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, team_shmem_hash_size - team_shift);
    		  Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, vector_size ), [&] (nnz_lno_t i) {
    			  keys[teamind * vector_size + i] = init_value; vals[teamind * vector_size + i] = 0;
    		  });
    	  });
      }


#if 0
      teamMember.team_barrier();

      Kokkos::single(Kokkos::PerTeam(teamMember),[&] () {

      for (int i = 0; i < team_shmem_hash_size; ++i){
    	  if (begins[i] != init_value){
    		  std::cout << "row_index:" << row_index << " i:" << i << " team_shmem_hash_size:" << team_shmem_hash_size << " is not init_value begins[i]:" << begins[i] << std::endl;
    	  }
      }
      });

      teamMember.team_barrier();
#endif
      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerTeam(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
#if 0
        globally_used_hash_count[0] = 0;
#endif
      });
#if 0

      teamMember.team_barrier();
#endif
#if 0
      bool is_global_alloced = false;
      nnz_lno_t *globally_used_hash_indices = NULL;
#endif
      const size_type a_col_begin_offset = row_mapA[row_index];


      nnz_lno_t a_col_ind = entriesA[a_col_begin_offset];
      scalar_t a_col_val = valuesA[a_col_begin_offset];

      nnz_lno_t current_a_column_offset_inrow = 0;
      nnz_lno_t flops_on_the_left_of_offsett = 0;
      size_type current_b_read_offsett = row_mapB[a_col_ind];
      nnz_lno_t current_a_column_flops = row_mapB[a_col_ind + 1] - current_b_read_offsett;



      //nnz_lno_t ii = left_work;
      nnz_lno_t row_flops = flops_per_row(row_index);

#if 1
        teamMember.team_barrier();
#endif

        for (nnz_lno_t vector_read_shift = vector_shift; vector_read_shift< row_flops; vector_read_shift += bs){
        {
              	  nnz_lno_t my_b_col_shift = vector_read_shift - flops_on_the_left_of_offsett;
            	  nnz_lno_t my_b_col = init_value; scalar_t my_b_val = 0; nnz_lno_t hash = init_value;
            	  int fail = 0;


			  if (my_b_col_shift >= current_a_column_flops ){

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

			  //now insert it to first level hashmap accumulator.
			  hash = (my_b_col * HASHSCALAR) & team_cuckoo_hash_func;
			  fail = 1;

			  for (nnz_lno_t trial = hash; trial < team_cuckoo_key_size; ){

				  if (keys[trial] == my_b_col){
					  Kokkos::atomic_add(vals + trial, my_b_val);
					  fail = 0;
					  break;
				  }
				  else if (keys[trial] == init_value){
					  if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)){
						  Kokkos::atomic_add(vals + trial, my_b_val);
						  fail = 0;
						  break;
					  }
				  }
				  else {
					  ++trial;
				  }
			  }
			  if (fail){
				  for (nnz_lno_t trial = 0; trial < hash; ){

					  if (keys[trial] == my_b_col){
						  Kokkos::atomic_add(vals + trial, my_b_val);
						  fail = 0;
						  break;
					  }
					  else if (keys[trial] == init_value){
						  if (Kokkos::atomic_compare_exchange_strong(keys + trial, init_value, my_b_col)){
							  Kokkos::atomic_add(vals + trial, my_b_val);
							  fail = 0;
							  break;
						  }
					  }
					  else {
						  ++trial;
					  }

				  }
			  }
		  }
      }

      teamMember.team_barrier();
      for (nnz_lno_t  my_index = vector_shift; my_index < team_cuckoo_key_size;my_index += bs){
    	  nnz_lno_t my_key = keys[my_index];
    	  if (my_key != init_value){
    		  scalar_t my_val = vals[my_index];
    		  nnz_lno_t write_index = Kokkos::atomic_fetch_add(used_hash_sizes, atomic_incr_type(1));
    		  c_row[write_index] = my_key;
    		  c_row_vals[write_index] = my_val;
    	  }
      }

    }
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
  KokkosSparse::SPGEMMAlgorithm algorithm_to_run = this->spgemm_algorithm;
  nnz_lno_t brows = row_mapB.dimension_0() - 1;
  size_type bnnz =  valsB.dimension_0();

  int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  size_t shmem_size_to_use = shmem_size;

  row_lno_persistent_work_view_t flops_per_row = this->handle->get_spgemm_handle()->row_flops;
  size_t original_overall_flops = this->handle->get_spgemm_handle()->original_overall_flops;
  nnz_lno_t max_nnz = this->handle->get_spgemm_handle()->get_max_result_nnz();
  size_type overall_nnz = this->handle->get_spgemm_handle()->get_c_nnz();

  typedef KokkosKernels::Impl::UniformMemoryPool< MyTempMemorySpace, nnz_lno_t> pool_memory_space;
  nnz_lno_t min_hash_size = 1;
  size_t chunksize  = 1;
  double first_level_cut_off  = this->handle->get_spgemm_handle()->get_first_level_hash_cut_off();
  int hash_scaler = this->handle->get_spgemm_handle()->get_min_hash_size_scale();
  nnz_lno_t tmp_max_nnz =  max_nnz;

  if (hash_scaler == 0){
	  tmp_max_nnz = KOKKOSKERNELS_MACRO_MAX(max_nnz, nnz_lno_t (this->b_col_cnt / this->concurrency + 1));
  }
  else {
	  tmp_max_nnz *= hash_scaler;
  }

  //START OF SHARED MEMORY SIZE CALCULATIONS
  nnz_lno_t unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t);
  nnz_lno_t team_shmem_key_size = ((shmem_size_to_use - sizeof(nnz_lno_t) * 4) / unit_memory);
  nnz_lno_t thread_memory = (shmem_size_to_use /8 / suggested_team_size) * 8;


  nnz_lno_t thread_shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
  if (KOKKOSKERNELS_VERBOSE){
	  std::cout << "\t\tinitial NumericCMEM -- thread_memory:" << thread_memory  << " unit_memory:" << unit_memory <<" initial key size:" << thread_shmem_key_size << std::endl;
	  std::cout << "\t\tinitial NumericCMEM -- team_memory:" << shmem_size_to_use  << " unit_memory:" << unit_memory <<" initial team key size:" << team_shmem_key_size << std::endl;
  }
  nnz_lno_t thread_shmem_hash_size = 1;
  while (thread_shmem_hash_size * 2 <=  thread_shmem_key_size){
	  thread_shmem_hash_size = thread_shmem_hash_size * 2;
  }
  nnz_lno_t team_shmem_hash_size = 1;
  while (team_shmem_hash_size * 2 <=  team_shmem_key_size){
	  team_shmem_hash_size = team_shmem_hash_size * 2;
  }
  //nnz_lno_t team_shared_memory_hash_func = team_shmem_hash_size - 1;

  team_shmem_key_size = team_shmem_key_size + ((team_shmem_key_size - team_shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
  team_shmem_key_size = (team_shmem_key_size >> 1) << 1;

  thread_shmem_key_size = thread_shmem_key_size + ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
  thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;

  //choose parameters
  if (this->spgemm_algorithm == SPGEMM_KK || SPGEMM_KK_LP == this->spgemm_algorithm){
	  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
		  //then chose the best method and parameters.
		  size_type average_row_nnz = overall_nnz / this->a_row_cnt;
		  size_t average_row_flops = original_overall_flops / this->a_row_cnt;
		  //if we have very low flops per row, or our maximum number of nnz is prett small,
		  //then we do row-base algorithm.
		  if (SPGEMM_KK_LP != this->spgemm_algorithm && (average_row_nnz < 32 || average_row_flops < 256)){
			  algorithm_to_run = SPGEMM_KK_MEMORY;
			  //if (average_row_nnz / double (thread_shmem_key_size) > 1.5){
				  while (average_row_nnz > size_type (thread_shmem_key_size) && suggested_vector_size < 32){
					  suggested_vector_size  = suggested_vector_size * 2;
					  suggested_vector_size = KOKKOSKERNELS_MACRO_MIN(32, suggested_vector_size);
					  suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
					  thread_memory = (shmem_size_to_use /8 / suggested_team_size) * 8;
					  thread_shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
					  thread_shmem_hash_size = 1;
					  while (thread_shmem_hash_size * 2 <=  thread_shmem_key_size){
						  thread_shmem_hash_size = thread_shmem_hash_size * 2;
					  }
					  thread_shmem_key_size = thread_shmem_key_size + ((thread_shmem_key_size - thread_shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
					  thread_shmem_key_size = (thread_shmem_key_size >> 1) << 1;
				  }

			  if (KOKKOSKERNELS_VERBOSE){
				  std::cout << "\t\t\tRunning KKMEM with suggested_vector_size:" << suggested_vector_size << " suggested_team_size:" << suggested_team_size << std::endl;
			  }
		  }
		  else {
			  nnz_lno_t tmp_team_cuckoo_key_size = ((shmem_size_to_use - sizeof(nnz_lno_t) * 2) / (sizeof(nnz_lno_t) + sizeof(scalar_t )));
			  int team_cuckoo_key_size = 1;
			  while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
			  suggested_vector_size = 32;
			  suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
			  algorithm_to_run = SPGEMM_KK_MEMORY_BIGSPREADTEAM;
			  while (average_row_nnz < team_cuckoo_key_size / 2 * (KOKKOSKERNELS_MACRO_MIN (first_level_cut_off + 0.05, 1))){
				  shmem_size_to_use = shmem_size_to_use / 2;
				  tmp_team_cuckoo_key_size = ((shmem_size_to_use - sizeof(nnz_lno_t) * 2) / (sizeof(nnz_lno_t) + sizeof(scalar_t )));
				  team_cuckoo_key_size = 1;
				  while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;

				  suggested_team_size = suggested_team_size / 2;
			  }
			  if (average_row_flops > size_t(2) * suggested_team_size * suggested_vector_size &&
					  average_row_nnz > size_type (team_cuckoo_key_size) * (KOKKOSKERNELS_MACRO_MIN (first_level_cut_off + 0.05, 1))){
				  shmem_size_to_use = shmem_size_to_use * 2;
				  tmp_team_cuckoo_key_size = ((shmem_size_to_use - sizeof(nnz_lno_t) * 2) / (sizeof(nnz_lno_t) + sizeof(scalar_t )));
				  team_cuckoo_key_size = 1;
				  while (team_cuckoo_key_size * 2 < tmp_team_cuckoo_key_size) team_cuckoo_key_size = team_cuckoo_key_size * 2;
				  suggested_team_size = suggested_team_size *2;
			  }
#ifdef FIRSTPARAMS
			  suggested_team_size = KOKKOSKERNELS_MACRO_MAX(4, suggested_team_size);
#else
			  suggested_team_size = KOKKOSKERNELS_MACRO_MAX(2, suggested_team_size);
#endif
			  if (max_nnz < team_cuckoo_key_size * KOKKOSKERNELS_MACRO_MIN (first_level_cut_off + 0.20, 1)){
				  algorithm_to_run = SPGEMM_KK_MEMORY_SPREADTEAM;
				  if (KOKKOSKERNELS_VERBOSE){
					  std::cout
					  	  	  << "\t\t\tRunning SPGEMM_KK_MEMORY_SPREADTEAM with suggested_vector_size:" << suggested_vector_size
							  << " suggested_team_size:" << suggested_team_size
							  << " shmem_size_to_use:" << shmem_size_to_use << std::endl;
				  }
			  }
			  else {
				  if (KOKKOSKERNELS_VERBOSE){
					  std::cout 	<< "\t\t\tRunning SPGEMM_KK_MEMORY_BIGSPREADTEAM with suggested_vector_size:" << suggested_vector_size
							  << " suggested_team_size:" << suggested_team_size
							  << " shmem_size_to_use:" << shmem_size_to_use << std::endl;
				  }
			  }
		  }
	  }
	  else {
		  bool run_dense = false;
		  nnz_lno_t max_column_cut_off = this->handle->get_spgemm_handle()->MaxColDenseAcc;
		  nnz_lno_t col_size = this->b_col_cnt;
		  if (col_size < max_column_cut_off){
			  run_dense = true;
			  if (KOKKOSKERNELS_VERBOSE){
				  std::cout << "\t\t\tRunning SPGEMM_KK_DENSE col_size:" << col_size << " max_column_cut_off:" << max_column_cut_off << std::endl;
			  }
		  }
		  else {
			  //round up maxNumRoughNonzeros to closest power of 2.
			  nnz_lno_t tmp_min_hash_size = 1;
			  while (tmp_max_nnz > tmp_min_hash_size){
			    tmp_min_hash_size *= 4;
			  }

			  size_t kkmem_chunksize = tmp_min_hash_size ; //this is for used hash indices
			  kkmem_chunksize += tmp_min_hash_size ; //this is for the hash begins
			  kkmem_chunksize += max_nnz ; //this is for hash nexts
			  kkmem_chunksize = kkmem_chunksize * sizeof (nnz_lno_t);
			  size_t dense_chunksize = (col_size + col_size / sizeof(scalar_t) + 1) * sizeof(scalar_t);


			  if (kkmem_chunksize >= dense_chunksize * 0.5){
				  run_dense = true;
				  if (KOKKOSKERNELS_VERBOSE){
					  std::cout << "\t\t\tRunning SPGEMM_KK_SPEED kkmem_chunksize:" << kkmem_chunksize << " dense_chunksize:" << dense_chunksize << std::endl;
				  }
			  }
			  else {
				  run_dense = false;
				  if (KOKKOSKERNELS_VERBOSE){
					  std::cout << "\t\t\tRunning SPGEMM_KK_MEMORY col_size:" << col_size << " max_column_cut_off:" << max_column_cut_off << std::endl;
				  }
			  }
		  }

		  if (run_dense){
			  this->KokkosSPGEMM_numeric_speed(
					  rowmapC_,
					  entriesC_,
					  valuesC_,
					  my_exec_space);
			  return;
		  }
	  }
  }
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);
  if (KOKKOSKERNELS_VERBOSE){
	  std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << thread_shmem_hash_size  << " shmem_key_size:" << thread_shmem_key_size << std::endl;
	  std::cout << "\t\tNumericCMEM -- adjusted team hashsize:" << team_shmem_hash_size  << " team_shmem_key_size:" << team_shmem_key_size << std::endl;
  }
  //END OF SHARED MEMORY SIZE CALCULATIONS


  //required memory for L2
  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){

	  if (algorithm_to_run == SPGEMM_KK_MEMORY_SPREADTEAM){
		  tmp_max_nnz = 1;
	  }
	  else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGSPREADTEAM){

	  }
	  else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGTEAM || algorithm_to_run == SPGEMM_KK_MEMORY_TEAM ){
		  //tmp_max_nnz -= team_shmem_key_size;
	  }
	  else {
		  //tmp_max_nnz -= thread_shmem_key_size;
	  }
  }

  if (algorithm_to_run == SPGEMM_KK_LP ){

	  while (tmp_max_nnz > min_hash_size){
		  min_hash_size *= 4;
	  }
	  chunksize = min_hash_size; //this is for used hash keys
	  chunksize += max_nnz; //this is for used hash keys
	  chunksize += min_hash_size * sizeof(scalar_t) / sizeof(nnz_lno_t) ; //this is for the hash values
  }
  else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGSPREADTEAM){
	  while (tmp_max_nnz > min_hash_size){
	      min_hash_size *= 2;  //try to keep it as low as possible because hashes are not tracked.
 	  }
	  chunksize = min_hash_size; //this is for used hash keys
	  chunksize += min_hash_size * sizeof(scalar_t) / sizeof(nnz_lno_t) ; //this is for the hash values
  }
  else{
	  while (tmp_max_nnz > min_hash_size){
	    min_hash_size *= 4;
	  }
	  chunksize = min_hash_size; //this is for used hash indices
	  chunksize += min_hash_size ; //this is for the hash begins
	  chunksize += max_nnz; //this is for hash nexts
  }
  int num_chunks = concurrency / suggested_vector_size;

#if defined( KOKKOS_ENABLE_CUDA )
  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA) {

    size_t free_byte ;
    size_t total_byte ;
    cudaMemGetInfo( &free_byte, &total_byte ) ;
    size_t required_size = size_t (num_chunks) * chunksize * sizeof(nnz_lno_t);
    if (KOKKOSKERNELS_VERBOSE)
      std::cout << "\tmempool required size:" << required_size << " free_byte:" << free_byte << " total_byte:" << total_byte << std::endl;
    if (required_size + num_chunks > free_byte){
      num_chunks = ((((free_byte - num_chunks)* 0.5) /8 ) * 8) / sizeof(nnz_lno_t) / chunksize;
    }
    {
      nnz_lno_t min_chunk_size = 1;
      while (min_chunk_size * 2 <= num_chunks) {
        min_chunk_size *= 2;
      }
      num_chunks = min_chunk_size;
    }
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
      shmem_size_to_use,
      suggested_vector_size,
      m_space,
      min_hash_size, max_nnz,
      suggested_team_size,

      my_exec_space,
      team_row_chunk_size,
	  first_level_cut_off, flops_per_row,
	  KOKKOSKERNELS_VERBOSE);


  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tvector_size:" << suggested_vector_size  << " chunk_size:" << team_row_chunk_size << " suggested_team_size:" << suggested_team_size<< std::endl;
  }
  timer1.reset();

  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
	  if (algorithm_to_run == SPGEMM_KK_MEMORY_SPREADTEAM){
		  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY_SPREADTEAM", gpu_team_policy4_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
		    MyExecSpace::fence();

	  }
	  else if (algorithm_to_run == SPGEMM_KK_MEMORY_BIGSPREADTEAM){
		  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY_BIGSPREADTEAM", gpu_team_policy6_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
	  }
	  else {

		  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY", gpu_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
	  }
    MyExecSpace::fence();
  }
  else {
	  if (algorithm_to_run == SPGEMM_KK_LP){
		  if (use_dynamic_schedule){
			  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::SPGEMM_KK_LP::DYNAMIC", dynamic_multicore_team_policy4_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
		  }
		  else {

			  Kokkos::parallel_for( "KOKKOSPARSE::SPGEMM::SPGEMM_KK_LP::STATIC", dynamic_multicore_team_policy4_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
		  }
	  }
	  else {
		  if (use_dynamic_schedule){

			  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::KKMEM::DYNAMIC",  dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
		  }
		  else {

			  Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::KKMEM::STATIC",  multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
		  }
	  }
	  MyExecSpace::fence();
  }

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
  }

}


//this is to isolate the memory use of accumulators and A,B,C.
//normally accumulators can use memory of C directly, but in this one we separate it for experimenting.
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
  double first_level_cut_off  = this->handle->get_spgemm_handle()->get_first_level_hash_cut_off();

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
	  team_row_chunk_size,
	  first_level_cut_off,
       this->handle->get_spgemm_handle()->row_flops, KOKKOSKERNELS_VERBOSE);


  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tvector_size:" << suggested_vector_size  << " chunk_size:" << team_row_chunk_size << std::endl;
  }
  timer1.reset();

  if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
    Kokkos::parallel_for("KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY2",  gpu_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    MyExecSpace::fence();
  }
  else {
    if (use_dynamic_schedule){

      Kokkos::parallel_for( "KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY_DYNAMIC", dynamic_multicore_team_policy2_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    else {

      Kokkos::parallel_for( "KOKKOSPARSE::SPGEMM::SPGEMM_KK_MEMORY_STATIC", multicore_team_policy2_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    MyExecSpace::fence();
  }

  if (KOKKOSKERNELS_VERBOSE){
    std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
  }

}

}
}
