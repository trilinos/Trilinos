/*
//@HEADER
// ************************************************************************
//
//          KokkosKernels: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP
//#define HASHTRACK

//#define TRACK_INSERTS
#define GPU_EXPERIMENTAL
//#define NUMERIC_USE_STATICMEM
//#define twostep
#include <KokkosKernels_Utils.hpp>
#include "../utils/KokkosKernels_SimpleUtils.hpp"
#include "../utils/KokkosKernels_SparseUtils.hpp"
#include "../utils/KokkosKernels_VectorUtils.hpp"
#include <KokkosKernels_HashmapAccumulator.hpp>
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"
#include "KokkosKernels_GraphColor.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
//#include <bitset>


//#define shmem_size 16128//12032//16128//16384

namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{


template <typename KernelHandle,
  typename alno_row_view_t_,
  typename alno_nnz_view_t_,
  typename ascalar_nnz_view_t_,
  typename blno_row_view_t_,
  typename blno_nnz_view_t_,
  typename bscalar_nnz_view_t_,
  typename clno_row_view_t_,
  typename clno_nnz_view_t_,
  typename cscalar_nnz_view_t_>
void spgemm_debug(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    alno_row_view_t_ row_mapA,
    alno_nnz_view_t_ entriesA,
    ascalar_nnz_view_t_ valuesA,

    bool transposeA,
    blno_row_view_t_ row_mapB,
    blno_nnz_view_t_ entriesB,
    bscalar_nnz_view_t_ valuesB,
    bool transposeB,
    clno_row_view_t_ row_mapC,
    clno_nnz_view_t_ &entriesC,
    cscalar_nnz_view_t_ &valuesC
    ){
  typename alno_row_view_t_::HostMirror h_rma = Kokkos::create_mirror_view (row_mapA);
  Kokkos::deep_copy (h_rma, row_mapA);
  typename alno_nnz_view_t_::HostMirror h_enta = Kokkos::create_mirror_view (entriesA);
  Kokkos::deep_copy (h_enta, entriesA);
  typename ascalar_nnz_view_t_::HostMirror h_vala = Kokkos::create_mirror_view (valuesA);
  Kokkos::deep_copy (h_vala, valuesA);

  typename blno_row_view_t_::HostMirror h_rmb = Kokkos::create_mirror_view (row_mapB);
  Kokkos::deep_copy (h_rmb, row_mapB);
  typename blno_nnz_view_t_::HostMirror h_entb = Kokkos::create_mirror_view (entriesB);
  Kokkos::deep_copy (h_entb, entriesB);
  typename bscalar_nnz_view_t_::HostMirror h_valb = Kokkos::create_mirror_view (valuesB);
  Kokkos::deep_copy (h_valb, valuesB);
  typename clno_row_view_t_::HostMirror h_rmc = Kokkos::create_mirror_view (row_mapC);


  typedef typename KernelHandle::nnz_lno_t lno_t;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_t;

  std::vector<scalar_t> accumulator(k, 0);
  std::vector<bool> acc_flag(k, false);

  std::vector<lno_t> result_c_col_indices(k);
  std::vector<scalar_t> result_c_col_values(k);

  const size_type alloc_step = 10000;
  size_type result_index = 0;
  size_type current_c_size = k;

  h_rmc(0) = 0;
  for (lno_t i = 0; i < m; ++i){
    const size_type a_row_begin = h_rma(i);
    const size_type a_row_end = h_rma(i + 1);
    lno_t a_row_size = a_row_end - a_row_begin;

    for (lno_t j = 0; j < a_row_size; ++j){
      size_type ind = a_row_begin + j;
      lno_t col = h_enta(ind);
      scalar_t val = h_vala(ind);
      //if (i == 0) std::cout << "a row:" <<  col << std::endl;
      const size_type b_row_begin = h_rmb(col);
      const size_type b_row_end = h_rmb(col + 1);
      lno_t b_row_size = b_row_end - b_row_begin;
      for (lno_t z = 0; z < b_row_size; ++z){
        size_type ind = b_row_begin + z;
        lno_t b_col = h_entb(ind);
        scalar_t b_val = h_valb(ind);
        //if (i == 0) std::cout << "\tb col:" <<  b_col << std::endl;
        if (acc_flag[b_col] == false){
          acc_flag[b_col] = true;
          result_c_col_indices[result_index++] = b_col;
          if (current_c_size == result_index){
            current_c_size += alloc_step;
            result_c_col_indices.resize(current_c_size);
            result_c_col_values.resize(current_c_size);
          }
        }
        accumulator[b_col] += b_val * val;
      }
    }
    h_rmc(i+1) = result_index;
    size_type c_row_begin = h_rmc(i);
    lno_t c_row_size =  result_index - c_row_begin;

    //if (i == 0) std::cout << "result_cols" << std::endl;

    for (lno_t j = 0; j < c_row_size; ++j){

      size_type ind = c_row_begin + j;
      lno_t result_col = result_c_col_indices[ind];
      //if (i == 0) std::cout << result_col << std::endl;
      result_c_col_values[ind] = accumulator[result_col];
      accumulator[result_col] = 0;
      acc_flag[result_col] = false;
    }

  }

  entriesC = clno_nnz_view_t_("entriesC", result_index);
  valuesC = cscalar_nnz_view_t_("entriesC", result_index);

  typename clno_nnz_view_t_::HostMirror h_entc = Kokkos::create_mirror_view (entriesC);
  typename cscalar_nnz_view_t_::HostMirror h_valc = Kokkos::create_mirror_view (valuesC);

  for (size_type i = 0; i < result_index; ++i){
    h_entc(i) = result_c_col_indices[i];
    h_valc(i) = result_c_col_values[i];
  }

  Kokkos::deep_copy (entriesC, h_entc);
  Kokkos::deep_copy (valuesC, h_valc);


}


template <typename HandleType,
  typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
  typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
class KokkosSPGEMM{
public:

  typedef a_row_view_t_ a_row_view_t;
  typedef a_lno_nnz_view_t_ a_in_lno_nnz_view_t;
  typedef a_scalar_nnz_view_t_ a_in_scalar_nnz_view_t;

  typedef b_lno_row_view_t_ b_in_lno_row_view_t;
  typedef b_lno_nnz_view_t_ b_in_lno_nnz_view_t;
  typedef b_scalar_nnz_view_t_ b_in_scalar_nnz_view_t;



  typedef typename a_row_view_t::non_const_value_type size_type;
  typedef typename a_row_view_t::const_value_type const_size_type;


  typedef typename a_in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename a_in_lno_nnz_view_t::const_value_type const_nnz_lno_t;

  typedef typename a_in_scalar_nnz_view_t::non_const_value_type scalar_t;
  typedef typename a_in_scalar_nnz_view_t::const_value_type const_scalar_t;


  typedef typename a_row_view_t::const_type const_a_lno_row_view_t;
  typedef typename a_row_view_t::non_const_type non_const_a_lno_row_view_t;

  typedef typename a_in_lno_nnz_view_t::const_type const_a_lno_nnz_view_t;
  typedef typename a_in_lno_nnz_view_t::non_const_type non_const_a_lno_nnz_view_t;

  typedef typename a_in_scalar_nnz_view_t::const_type const_a_scalar_nnz_view_t;
  typedef typename a_in_scalar_nnz_view_t::non_const_type non_const_a_scalar_nnz_view_t;


  typedef typename b_in_lno_row_view_t::const_type const_b_lno_row_view_t;
  typedef typename b_in_lno_row_view_t::non_const_type non_const_b_lno_row_view_t;

  typedef typename b_in_lno_nnz_view_t::const_type const_b_lno_nnz_view_t;
  typedef typename b_in_lno_nnz_view_t::non_const_type non_const_b_lno_nnz_view_t;

  typedef typename b_in_scalar_nnz_view_t::const_type const_b_scalar_nnz_view_t;
  typedef typename b_in_scalar_nnz_view_t::non_const_type non_const_b_scalar_nnz_view_t;

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;


  typedef typename HandleType::bool_persistent_view_t bool_persistent_view_t;
  typedef typename HandleType::bool_temp_view_t bool_temp_view_t;


  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t ;
  typedef typename team_policy_t::member_type team_member_t ;

  struct CountTag{};




  struct FillTag{};
  struct MultiCoreDenseAccumulatorTag{};
  struct MultiCoreTag{};
  struct GPUTag{};

  struct Numeric1Tag{};
  struct Numeric2Tag{};
  struct Numeric3Tag{};

  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag, MyExecSpace> multicore_dense_team_count_policy_t ;
  typedef Kokkos::TeamPolicy<MultiCoreTag, MyExecSpace> multicore_team_policy_t ;
  typedef Kokkos::TeamPolicy<GPUTag, MyExecSpace> gpu_team_policy_t ;
  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric1Tag, MyExecSpace> team_numeric1_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric2Tag, MyExecSpace> team_numeric2_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric3Tag, MyExecSpace> team_numeric3_policy_t ;


  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_multicore_dense_team_count_policy_t ;
  typedef Kokkos::TeamPolicy<MultiCoreTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_multicore_team_policy_t ;
  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_count_policy_t ;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_fill_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric1Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_numeric1_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric2Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_numeric2_policy_t ;
  typedef Kokkos::TeamPolicy<Numeric3Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_numeric3_policy_t ;


private:
  HandleType *handle;
  nnz_lno_t a_row_cnt;
  nnz_lno_t b_row_cnt;
  nnz_lno_t b_col_cnt;


  const_a_lno_row_view_t row_mapA;
  const_a_lno_nnz_view_t entriesA;
  const_a_scalar_nnz_view_t valsA;
  bool transposeA;

  const_b_lno_row_view_t row_mapB;
  const_b_lno_nnz_view_t entriesB;
  const_b_scalar_nnz_view_t valsB;
  bool transposeB;

  const size_t shmem_size;
  const size_t concurrency;
  const bool use_dynamic_schedule;
  const bool KOKKOSKERNELS_VERBOSE;
  //const int KOKKOSKERNELS_VERBOSE = 1;
public:
  KokkosSPGEMM(
      HandleType *handle_,
      nnz_lno_t m_,
      nnz_lno_t n_,
      nnz_lno_t k_,
      const_a_lno_row_view_t row_mapA_,
      const_a_lno_nnz_view_t entriesA_,
      bool transposeA_,
      const_b_lno_row_view_t row_mapB_,
      const_b_lno_nnz_view_t entriesB_,
      bool transposeB_):handle (handle_), a_row_cnt(m_), b_row_cnt(n_), b_col_cnt(k_),
          row_mapA(row_mapA_), entriesA(entriesA_), valsA(), transposeA(transposeA_),
          row_mapB(row_mapB_), entriesB(entriesB_), valsB(), transposeB(transposeB_),
          shmem_size(handle_->get_shmem_size()), concurrency(MyExecSpace::concurrency()),
          use_dynamic_schedule(handle_->is_dynamic_scheduling()), KOKKOSKERNELS_VERBOSE(handle_->get_verbose())
          //,row_mapC(), entriesC(), valsC()
          {}

  KokkosSPGEMM(
      HandleType *handle_,
      nnz_lno_t m_,
      nnz_lno_t n_,
      nnz_lno_t k_,
        const_a_lno_row_view_t row_mapA_,
        const_a_lno_nnz_view_t entriesA_,
        const_a_scalar_nnz_view_t valsA_,
        bool transposeA_,
        const_b_lno_row_view_t row_mapB_,
        const_b_lno_nnz_view_t entriesB_,
        const_b_scalar_nnz_view_t valsB_,
        bool transposeB_):handle (handle_), a_row_cnt(m_), b_row_cnt(n_), b_col_cnt(k_),
            row_mapA(row_mapA_), entriesA(entriesA_), valsA(valsA_), transposeA(transposeA_),
            row_mapB(row_mapB_), entriesB(entriesB_), valsB(valsB_), transposeB(transposeB_),
            shmem_size(handle_->get_shmem_size()), concurrency(MyExecSpace::concurrency()),
            use_dynamic_schedule(handle_->is_dynamic_scheduling()), KOKKOSKERNELS_VERBOSE(handle_->get_verbose())
            //,row_mapB(), entriesC(), valsC()
            {}


  /**
   * \brief Functor to calculate the max flops in a row of SPGEMM.
   *
   */
  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_oldrow_view_t, typename b_row_view_t>
  struct PredicMaxRowNNZ{
    nnz_lno_t m; //num rows
    a_row_view_t row_mapA;  //row pointers of a
    a_nnz_view_t entriesA;  //col
    b_oldrow_view_t row_begins_B;
    b_row_view_t row_end_indices_B;
    const size_type min_val;
    nnz_lno_t team_row_chunk_size;

    /**
     * \brief Constructor
     * \param m_: num rows in A.
     * \param row_mapA_: row pointers of A
     * \param entriesA_: col indices of A
     * \param row_begins_B_: row begin indices of B
     * \param row_end_indices_B_: row end indices of B
     * \param team_row_chunk_size_: the number of rows assigned to each team.
     */
    PredicMaxRowNNZ(
        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_oldrow_view_t row_begins_B_,
        b_row_view_t row_end_indices_B_,
        nnz_lno_t team_row_chunk_size_):
          m(m_),
          row_mapA(row_mapA_), entriesA(entriesA_),
          row_begins_B(row_begins_B_),
          row_end_indices_B(row_end_indices_B_),
          min_val(((std::numeric_limits<size_type>::lowest()))),
          team_row_chunk_size(team_row_chunk_size_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember, size_type &overal_max) const {
      //get the range of rows for team.
      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, m);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index)
      {

        const size_type col_begin = row_mapA[row_index];
        const size_type col_end = row_mapA[row_index + 1];
        const nnz_lno_t left_work = col_end - col_begin;

        size_type max_num_results_in_row = 0;

        //get the size of the rows of B, pointed by row of A
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, left_work),
            [&] (nnz_lno_t i, size_type & valueToUpdate) {
          const size_type adjind = i + col_begin;
          const nnz_lno_t colIndex = entriesA[adjind];
          valueToUpdate += row_end_indices_B (colIndex) - row_begins_B(colIndex);
        },
        max_num_results_in_row);
        //set max.
        if (overal_max < max_num_results_in_row) {
          overal_max = max_num_results_in_row;
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile size_type& dst,const volatile size_type& src) const {
      if (dst < src) { dst = src;}
    }


    KOKKOS_INLINE_FUNCTION
    void init (size_type& dst) const
    {
      dst = min_val;
    }
  };

  template <typename row_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct copyMatrix{

    const nnz_lno_t numrows;

    row_view_t row_map;
    new_nnz_view_t set_index_entries;
    new_nnz_view_t set_entries;


    new_row_view_t new_row_map;
    new_nnz_view_t out_set_index_entries;
    new_nnz_view_t out_set_entries;
    nnz_lno_t team_row_chunk_size;

    copyMatrix(
        row_view_t row_map_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_,

        new_row_view_t new_row_map_,
        new_row_view_t out_set_index_entries_,
        new_row_view_t out_set_entries_,
        nnz_lno_t team_row_chunk_size_
        ):
      row_map(row_map_),

      new_row_map(new_row_map_),
      set_index_entries(set_index_entries_),
      set_entries(set_entries_),

      out_set_index_entries(out_set_index_entries_),
      out_set_entries(out_set_entries_),
      numrows(row_map_.dimension_0() - 1),
      team_row_chunk_size(team_row_chunk_size_)
      {}


    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {



      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& i)
      {

        size_type rowBegin = new_row_map(i);
        nnz_lno_t left_work = new_row_map(i + 1) - rowBegin;
        size_type oldRowBegin= row_map(i);

        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, left_work),
            [&] (nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          const size_type oldadjind = i + oldRowBegin;

          out_set_index_entries(adjind) = set_index_entries(oldadjind);
          out_set_entries (adjind) =  set_entries(oldadjind);
        });
      });
    }

  };

  /**
   *\brief Functor to zip the B matrix.
   */
  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t>
  struct SingleStepZipMatrix{


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

    nnz_lno_t * pset_index_entries;
    nnz_lno_t * pset_entries;


    //pool_memory_space m_space;
    const int shared_memory_size;
    const nnz_lno_t team_row_chunk_size;

    SingleStepZipMatrix(
        const row_view_t row_map_,
        const nnz_view_t entries_,
        const nnz_lno_t compression_bit_mask_,
        const int compression_bit_divide_shift_,
        const int vector_size_,
        new_row_view_t new_row_map_,
        /*Global memory hash space. in the size of nnz*/
        new_nnz_view_t set_index_begins_,
        new_nnz_view_t set_index_nexts_,
        new_nnz_view_t set_index_entries_,
        new_nnz_view_t set_entries_,
        const int shared_mem,
        const nnz_lno_t team_row_chunk_size_
        ):
      numrows(row_map_.dimension_0() - 1),
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
      pset_index_begins(set_index_begins_.ptr_on_device()),
      pset_index_nexts(set_index_nexts_.ptr_on_device()),
      pset_index_entries(set_index_entries_.ptr_on_device()),
      pset_entries(set_entries_.ptr_on_device()),
      shared_memory_size(shared_mem),
      team_row_chunk_size(team_row_chunk_size_)
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_ind)
      {


      //CPU part only compares it with the previous index.
      size_type rowBegin = row_map(row_ind);
      nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

      nnz_lno_t used_size = 0;

      const nnz_lno_t n = entries(rowBegin);
      nnz_lno_t prev_nset_ind = n >> compression_bit_divide_shift;
      nnz_lno_t prev_nset = 1;
      prev_nset = prev_nset << (n & compression_bit_mask);




      for (nnz_lno_t i = 1; i < left_work; ++i){
        nnz_lno_t n_set = 1;
        const size_type adjind = i + rowBegin;
        const nnz_lno_t n = entries(adjind);
        nnz_lno_t n_set_index = n >> compression_bit_divide_shift;
        n_set = n_set << (n & compression_bit_mask);
        if (prev_nset_ind == n_set_index){
          prev_nset = prev_nset | n_set;
        } else {
          pset_index_entries[used_size + rowBegin] = prev_nset_ind ;
          pset_entries[used_size + rowBegin] = prev_nset;
          ++used_size;
          prev_nset_ind = n_set_index;
          prev_nset = n_set;
        }
      }
      pset_index_entries[used_size + rowBegin] = prev_nset_ind ;
      pset_entries[used_size + rowBegin] = prev_nset;
      ++used_size;
      new_row_map(row_ind) = rowBegin + used_size;
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const GPUTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_ind = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      //std::cout << "i:" << i << std::endl;
      if (row_ind >= numrows) return;

      //std::cout << "i:" << i << std::endl;
      //how much shared memory a thread will have in team
      int thread_memory = shared_memory_size / teamMember.team_size();
      //allocate all shared memory
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      //allocate memory in the size of vectors
      nnz_lno_t *result_keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * vector_size;
      nnz_lno_t *result_vals = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * vector_size;

      thread_memory -= vector_size * sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;

      //calculate the memory needed for 1 single hash entry.
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2; //begins, nexts, and keys. No need for vals yet.
      //how many hash entries can be held in shared memory.
      nnz_lno_t shared_memory_hash_size = thread_memory / unit_memory;


      //points to the beginning of hashes
      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);;


      //this is a hash for individual row elements. therefore, the size can be at most
      //the number of elements in a row, so the nnz_lno_t can be used instead of size_type here.

      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      size_type rowBegin = row_map(row_ind);
      size_type rowBeginP = rowBegin;


      nnz_lno_t left_work = row_map(row_ind + 1) - rowBegin;

      //same as second level hash map.
      //second level hashmap.
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm2(left_work, left_work, pset_index_begins + rowBegin, pset_index_nexts+ rowBegin, pset_index_entries+ rowBegin, pset_entries+ rowBegin);



      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });
      //initialize begins.
      //these are initialized before the loop.
      /*
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i) {
        pset_index_begins[rowBegin + i] = -1;
      });
      */

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
      });

      //lno_t neighbor_set_count = 0;

      while (left_work){

        //std::cout << "left_work:" << left_work << std::endl;
        nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

        //first get the portion of the work for the vector lane.
        nnz_lno_t n_set_index = -1;
        nnz_lno_t n_set = 1;

        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, work_to_handle),
            [&] (nnz_lno_t i) {
          const size_type adjind = i + rowBegin;
          const nnz_lno_t n = entries(adjind);
          n_set_index = n >> compression_bit_divide_shift;
          n_set = n_set << (n & compression_bit_mask);
        });



        //it is possible that multiple threads have same values.
        //first merge them, as a result of this operation we will have the n_sets merged,
        //if a thread's value merged to some other threads we have n_set = -1.
        hm.vector_mergeOr_MEM(teamMember, vector_size, n_set_index,n_set, result_keys, result_vals);



        nnz_lno_t hash = n_set_index % shared_memory_hash_size;
        if (n_set_index == -1) hash = -1;
        int overall_num_unsuccess = 0;

        int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
                                          teamMember, vector_size, hash,n_set_index, n_set, used_hash_sizes, shared_memory_hash_size);

        Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &overall_num_unsuccess) {
          overall_num_unsuccess += num_unsuccess;
        }, overall_num_unsuccess);


        //if one of the inserts was successfull, which means we run out shared memory
        if (overall_num_unsuccess){
          nnz_lno_t hash = -1;
          if (num_unsuccess) hash = n_set_index % hm2.hash_key_size;

          //int insertion =
          hm2.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size, hash,n_set_index,n_set, used_hash_sizes + 1, hm2.max_value_size);
        }

        left_work -= work_to_handle;
        rowBegin += work_to_handle;
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
        new_row_map(row_ind) = rowBeginP + used_hash_sizes[0] + used_hash_sizes[1];
      });


      size_type written_index = used_hash_sizes[1];
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, used_hash_sizes[0]),
          [&] (nnz_lno_t i) {
        pset_index_entries[rowBeginP + written_index + i] = keys[i];
        pset_entries[rowBeginP + written_index + i] = vals[i];
      });
    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };



  /**
   * \brief function return max flops for a row in the result multiplication.
   * \param m: number of rows in A
   * \param row_mapA: row pointers of A.
   * \param entriesA: column indices of A
   * \param row_pointers_begin_B: beginning of the row indices for B
   * \param row_pointers_end_B: end of the row indices for B
   */
  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_oldrow_view_t, typename b_row_view_t>
  size_t getMaxRoughRowNNZ(
      nnz_lno_t m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,

      b_oldrow_view_t row_pointers_begin_B,
      b_row_view_t row_pointers_end_B){


    //get the execution space type.
    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = this->handle->get_handle_exec_space();
    int suggested_vector_size = this->handle->get_suggested_vector_size(m, entriesA.dimension_0());
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency , m);

    PredicMaxRowNNZ<a_row_view_t, a_nnz_view_t, b_oldrow_view_t, b_row_view_t>
      pcnnnz(
        m,
        row_mapA,
        entriesA,
        row_pointers_begin_B,
        row_pointers_end_B,
        team_row_chunk_size );


    typename b_oldrow_view_t::non_const_value_type rough_size = 0;
    Kokkos::parallel_reduce( team_policy_t(m / team_row_chunk_size  + 1 , suggested_team_size, suggested_vector_size), pcnnnz, rough_size);
    MyExecSpace::fence();
    return rough_size;
  }

  template <
      typename out_row_view_t,
      typename out_nnz_view_t,
      typename in_row_map_t,
      typename in_nnz_t>
  struct unzipMatrix{
    typedef typename out_row_view_t::non_const_type non_const_c_lno_row_view_t;
    size_type m;
    in_row_map_t in_row_map;
    in_nnz_t in_set_index_entries;
    in_nnz_t in_set_entries;
    out_row_view_t out_rowmap;
    out_nnz_view_t out_entries;
    int set_size;
    size_t shared_memory_size;

    unzipMatrix(size_type m_,
                in_row_map_t row_map_c_copy_,
                in_nnz_t c_set_index_entries_,
                in_nnz_t c_set_entries_,
                out_row_view_t rowmapC_,
                out_row_view_t entriesC_, size_t shared_memory_size_): m(m_), in_row_map(row_map_c_copy_),
                    in_set_index_entries(c_set_index_entries_), in_set_entries(c_set_entries_),
                    out_rowmap(rowmapC_), out_entries(entriesC_), set_size(sizeof(typename in_nnz_t::value_type) * 8),
                    shared_memory_size(shared_memory_size_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const CountTag&, const team_member_t & teamMember) const {
      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = teamMember.league_rank()  * team_size + team_rank;

      if (row_index >= m) return;
      //printf("row_index:%d m:%d\n", row_index, m);
      //check ii is out of range. if it is, just return.
      const size_type col_begin = in_row_map[row_index];
      const size_type col_end = in_row_map[row_index + 1];



      size_type nnz = 0;

      Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
                    [&] (size_type i, size_type &nnz) {
        const size_type col_a_index = i + col_begin;
        nnz_lno_t c_rows = in_set_entries[col_a_index];
                size_type c = 0;
        while (c_rows){
          c_rows &= (c_rows-1) ;
          c++;
        }
        nnz += c;
      }, nnz);

      //out_rowmap(row_index) = nnz;
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        out_rowmap(row_index) = nnz;
      });

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const FillTag&, const team_member_t & teamMember) const {

      int team_size = teamMember.team_size();
      int team_rank = teamMember.team_rank();
      size_type row_index = teamMember.league_rank()  * team_size + team_rank;


      size_type *adj_index = (size_type *) teamMember.team_shmem().get_shmem(team_size);



      if (row_index >= m) return;

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        adj_index[team_rank] = out_rowmap(row_index);
      });
      //check ii is out of range. if it is, just return.
      const size_type col_begin = in_row_map[row_index];
      const size_type col_end = in_row_map[row_index + 1];
      size_type c = 0;
      //size_type out_row_map_index = out_rowmap(row_index);

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, col_end - col_begin),
          [&] (size_type i) {
        const size_type col_a_index = i + col_begin;

        nnz_lno_t c_rows = in_set_entries[col_a_index];
        nnz_lno_t c_rows_set_index = in_set_index_entries[col_a_index];
        int current_row = 0;
        nnz_lno_t unit = 1;


        while (c_rows){
          if (c_rows & unit){

            size_type wind = Kokkos::atomic_fetch_add(adj_index + team_rank , 1);
            out_entries(wind) = set_size * c_rows_set_index + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }

      });
    }


    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };


  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
            typename mpool_type>
  struct NumericCMEM_CPU{
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
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
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
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_,
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
          pEntriesC(entriesC_.ptr_on_device()), pVals(valuesC.ptr_on_device()),
          my_exec_space(my_exec_space_),
          team_work_size(team_row_chunk_size){
          }


    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_OPENMP )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
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

    }

  };

  template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__,
            typename b_row_view_t__, typename b_nnz_view_t__, typename b_scalar_view_t__,
            typename c_row_view_t__, typename c_nnz_view_t__, typename c_scalar_view_t__>
  struct NumericCMEM{
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

    c_nnz_view_t__ beginsC;
    c_nnz_view_t__ nextsC;

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

        c_nnz_view_t__ beginsC_,
        c_nnz_view_t__ nextsC_,

        const size_type sharedMemorySize_,
        const int suggested_vector_size,
        const nnz_lno_t team_row_chunk_size,
        int suggested_team_size_,
        bool KOKKOSKERNELS_VERBOSE):
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
          pbeginsC(beginsC_.ptr_on_device()), pnextsC(nextsC_.ptr_on_device()),
          pEntriesC(entriesC_.ptr_on_device()), pvaluesC(valuesC_.ptr_on_device()),
          shared_memory_size(sharedMemorySize_),

          vector_size (suggested_vector_size),
          team_work_size(team_row_chunk_size),

          unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t)),
          suggested_team_size(suggested_team_size_),
          thread_memory((shared_memory_size /8 / suggested_team_size_) * 8),
          shmem_key_size(), shared_memory_hash_func(), shmem_hash_size(1)
          {

            shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 2) / unit_memory);
            if (KOKKOSKERNELS_VERBOSE){
              std::cout << "\t\tNumericCMEM -- thread_memory:" << thread_memory  << " unit_memory:" << unit_memory <<
                  " initial key size:" << shmem_key_size << std::endl;
            }
            while (shmem_hash_size * 2 <=  shmem_key_size){
              shmem_hash_size = shmem_hash_size * 2;
            }
            shared_memory_hash_func = shmem_hash_size - 1;

            shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
            shmem_key_size = (shmem_key_size >> 1) << 1;

            if (KOKKOSKERNELS_VERBOSE){
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


      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm(shmem_hash_size, shmem_key_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm2(0, 0,
          NULL, NULL, NULL, NULL);
      /*
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
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
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          while (left_work){
            nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);
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
                [&] (const int threadid, int &overall_num_unsuccess) {
              overall_num_unsuccess += num_unsuccess;
            }, overall_num_unsuccess);

            if (overall_num_unsuccess){
              nnz_lno_t hash = -1;
              if (num_unsuccess) {
                hash = b_col_ind % global_memory_hash_size;
              }

              //int insertion =
              hm2.vector_atomic_insert_into_hash_mergeAdd(
                  teamMember, vector_size,
                  hash,b_col_ind,b_val,
                  used_hash_sizes + 1, hm2.max_value_size
              );

            }
            left_work -= work_to_handle;
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

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }
  };


  template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__,
            typename b_row_view_t__, typename b_nnz_view_t__, typename b_scalar_view_t__,
            typename c_row_view_t__, typename c_nnz_view_t__, typename c_scalar_view_t__>
  struct NumericCCOLOR{
    nnz_lno_t numrows;
    nnz_lno_t numcols;

    a_row_view_t__ row_mapA;
    a_nnz_view_t__ entriesA;
    a_scalar_view_t__ valuesA;

    b_row_view_t__ row_mapB;
    b_nnz_view_t__ entriesB;
    b_scalar_view_t__ valuesB;

    c_row_view_t__ rowmapC;
    c_nnz_view_t__ entriesC;
    c_scalar_view_t__ valuesC;
    nnz_lno_t *pEntriesC ;
    scalar_t *pVals;

    scalar_temp_work_view_t denseAccumulator; //initially all zeroes
    scalar_t *pdenseAccumulator; //initially all zeroes


    //bool_temp_view_t denseAccumulatorFlags; //initially all false.
    //bool * pdenseAccumulatorFlags; //initially all false.

    nnz_lno_t team_work_size;

    nnz_lno_t color_begin;
    nnz_lno_t color_end;
    nnz_lno_persistent_work_view_t color_adj;
    nnz_lno_persistent_work_view_t vertex_colors;

    nnz_lno_t consecutive_chunk_size;
    nnz_lno_t consecutive_all_color_chunk_size;

    nnz_lno_t chunk_divison;
    nnz_lno_t chunk_and;
    NumericCCOLOR(

        nnz_lno_t m_,
        nnz_lno_t k_,


        a_row_view_t__ row_mapA_,
        a_nnz_view_t__ entriesA_,
        a_scalar_view_t__ valuesA_,

        b_row_view_t__ row_mapB_,
        b_nnz_view_t__ entriesB_,
        b_scalar_view_t__ valuesB_,

        c_row_view_t__ rowmapC_,
        c_nnz_view_t__ entriesC_,
        c_scalar_view_t__ valuesC_,
        scalar_temp_work_view_t denseAccumulator_, //initially all zeroes
        //bool_temp_view_t denseAccumulatorFlags_, //initially all false.



        nnz_lno_t team_row_work_size_):
          numrows(m_), numcols(k_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          valuesA(valuesA_),

          row_mapB(row_mapB_),
          entriesB(entriesB_),
          valuesB(valuesB_),

          rowmapC(rowmapC_),
          entriesC(entriesC_),
          valuesC(valuesC_), pEntriesC(entriesC_.ptr_on_device()), pVals(valuesC.ptr_on_device()),
          denseAccumulator (denseAccumulator_), pdenseAccumulator(denseAccumulator_.ptr_on_device()),
          //denseAccumulatorFlags (denseAccumulatorFlags_), pdenseAccumulatorFlags(denseAccumulatorFlags_.ptr_on_device()),
          team_work_size(team_row_work_size_),
          color_begin(0),
          color_end(0),
          color_adj(),
          vertex_colors(), consecutive_chunk_size(numcols), consecutive_all_color_chunk_size(numcols), chunk_divison (0), chunk_and(0)
          {
          }

    //one color at a time.
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric1Tag&, const team_member_t & teamMember) const {

      const nnz_lno_t team_row_begin =
          teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end =
          KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end),
          [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        //nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator;// + numcols * color;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t acc_index = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t acc_index = my_entries[i];
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }


    //multi-color minimized writes
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {


      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator + numcols * color;

        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t acc_index = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t acc_index = my_entries[i];
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }

    //multi-color minimized reads
    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric3Tag&, const team_member_t & teamMember) const {


      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + color_begin;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, color_end);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& color_index_index) {
        nnz_lno_t row_index = color_adj(color_index_index);

        nnz_lno_t color = vertex_colors(row_index);
        scalar_t *mydenseAccumulator = pdenseAccumulator + numcols * color;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t left_work = nnz_lno_t(row_mapA[row_index + 1] - col_begin);
        const size_type c_row_begin = rowmapC[row_index];
        nnz_lno_t *my_entries = pEntriesC + c_row_begin;
        for (nnz_lno_t colind = 0; colind < left_work; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, left_work),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            nnz_lno_t col_ind = entriesB[adjind];
            nnz_lno_t acc_index = col_ind;
            //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
            scalar_t b_val = valuesB[adjind] * valA;
            mydenseAccumulator[acc_index] += b_val;
          });
        }
        scalar_t *my_vals = pVals + c_row_begin;

        nnz_lno_t row_size = rowmapC[row_index + 1] - c_row_begin;
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, row_size),
            [&] (nnz_lno_t i) {
          nnz_lno_t col_ind = my_entries[i];
          nnz_lno_t acc_index = col_ind;

          //nnz_lno_t acc_index = (col_ind  >> chunk_divison) * (consecutive_all_color_chunk_size) + (color << chunk_divison)+ (col_ind & chunk_and);
          my_vals[i] = mydenseAccumulator[acc_index];
          mydenseAccumulator[acc_index] = 0;
        });
      });
    }


    size_t team_shmem_size (int team_size) const {
      return team_size * sizeof (nnz_lno_t) * 8;
    }
  };


  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t,
            typename b_row_view_t, typename b_nnz_view_t, typename b_scalar_view_t,
            typename c_row_view_t, typename c_nnz_view_t, typename c_scalar_view_t,
            typename pool_memory_type>
  struct PortableNumericCHASH{

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
    const nnz_lno_t pow2_hash_func;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
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
        nnz_lno_t min_hash_size,
        int suggested_team_size_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_,
        nnz_lno_t team_row_chunk_size,
        bool KOKKOSKERNELS_VERBOSE):
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
          pow2_hash_size(min_hash_size),
          pow2_hash_func(min_hash_size - 1),
          my_exec_space(my_exec_space_),
          team_work_size(team_row_chunk_size),

          unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) + sizeof (scalar_t)),
          suggested_team_size(suggested_team_size_),
          thread_memory((shared_memory_size /8 / suggested_team_size_) * 8),
          shmem_key_size(), shared_memory_hash_func(), shmem_hash_size(1)
    {

      shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 4) / unit_memory);
      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tNumericCMEM -- thread_memory:" << thread_memory  << " unit_memory:" << unit_memory <<
            " initial key size:" << shmem_key_size << std::endl;
      }
      while (shmem_hash_size * 2 <=  shmem_key_size){
        shmem_hash_size = shmem_hash_size * 2;
      }
      shared_memory_hash_func = shmem_hash_size - 1;

      shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size) * sizeof(nnz_lno_t)) / (sizeof (nnz_lno_t) * 2 + sizeof(scalar_t));
      shmem_key_size = (shmem_key_size >> 1) << 1;

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tNumericCMEM -- adjusted hashsize:" << shmem_hash_size  << " shmem_key_size:" << shmem_key_size << std::endl;
      }


    }
    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_OPENMP )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }

    }

    //assumes that the vector lane is 1, as in cpus
    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, numrows);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
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

    KOKKOS_INLINE_FUNCTION
    void operator()(const Numeric2Tag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;

      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_sizes = 0;


      const size_type c_row_begin = rowmapC[row_index];
      const size_type c_row_end = rowmapC[row_index + 1];

      const nnz_lno_t global_memory_hash_size = nnz_lno_t(c_row_end - c_row_begin);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
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
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          for ( nnz_lno_t i = 0; i < left_work; ++i){
            const size_type adjind = i + rowBegin;
            nnz_lno_t b_col_ind = entriesB[adjind];
            scalar_t b_val = valuesB[adjind] * valA;
            nnz_lno_t hash = b_col_ind & pow2_hash_func;

            //this has to be a success, we do not need to check for the success.
            int insertion = hm2.sequential_insert_into_hash_mergeAdd_TrackHashes(
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

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
      hm(shmem_hash_size, shmem_key_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,scalar_t>
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

        for ( nnz_lno_t ii = 0; ii < left_work; ++ii){
          size_type a_col = col_begin + ii;
          nnz_lno_t rowB = entriesA[a_col];
          scalar_t valA = valuesA[a_col];

          size_type rowBegin = row_mapB(rowB);
          nnz_lno_t left_work = row_mapB(rowB + 1) - rowBegin;

          while (left_work){
            nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);
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
                [&] (const int threadid, int &overall_num_unsuccess) {
              overall_num_unsuccess += num_unsuccess;
            }, overall_num_unsuccess);

            if (overall_num_unsuccess){
              if (!is_global_alloced){
                volatile nnz_lno_t * tmp = NULL;
                size_t tid = get_thread_id(row_index);
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

              nnz_lno_t hash = -1;
              if (num_unsuccess) {
                hash = b_col_ind & pow2_hash_func;
              }

              //this has to be a success, we do not need to check for the success.
              //int insertion =
              hm2.vector_atomic_insert_into_hash_mergeAdd_TrackHashes(
                  teamMember, vector_size,
                  hash,b_col_ind,b_val,
                  used_hash_sizes + 1, hm2.max_value_size
                  ,globally_used_hash_count, globally_used_hash_indices
              );
            }
            left_work -= work_to_handle;
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


  /***
   * \brief Functor to calculate the row sizes of C.
   */
  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, //typename nnz_lno_temp_work_view_t,
            typename pool_memory_space>
  struct StructureC{
    const nnz_lno_t numrows; //num rows in A

    const a_row_view_t row_mapA; //A row pointers
    const a_nnz_view_t entriesA; // A column indices

    const b_original_row_view_t row_pointer_begins_B;
    const b_compressed_row_view_t row_pointer_ends_B;
    b_nnz_view_t entriesSetIndicesB;
    b_nnz_view_t entriesSetsB;

    c_row_view_t rowmapC;
    //nnz_lno_temp_work_view_t entriesSetIndicesC;
    //nnz_lno_temp_work_view_t entriesSetsC;


    const nnz_lno_t pow2_hash_size;
    const nnz_lno_t pow2_hash_func;
    const nnz_lno_t MaxRoughNonZero;

    const size_t shared_memory_size;
    const int vector_size;
    pool_memory_space m_space;
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;


    const int unit_memory; //begins, nexts, and keys. No need for vals yet.
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
    StructureC(
        const nnz_lno_t m_,
        const a_row_view_t row_mapA_,
        const a_nnz_view_t entriesA_,
        const b_original_row_view_t row_ptr_begins_B_,
        const b_compressed_row_view_t row_ptr_ends_B_,
        const b_nnz_view_t entriesSetIndicesB_,
        const b_nnz_view_t entriesSetsB_,
        c_row_view_t rowmapC_,
        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const int suggested_team_size_,
        const nnz_lno_t team_row_chunk_size_,
        const int vector_size_,
        pool_memory_space mpool_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_,
        bool KOKKOSKERNELS_VERBOSE):
          numrows(m_),
          row_mapA (row_mapA_),
          entriesA(entriesA_),
          row_pointer_begins_B(row_ptr_begins_B_),
          row_pointer_ends_B(row_ptr_ends_B_),
          entriesSetIndicesB(entriesSetIndicesB_),
          entriesSetsB(entriesSetsB_),
          rowmapC(rowmapC_),
          //entriesSetIndicesC(),
          //entriesSetsC(),
          pow2_hash_size(hash_size_),
          pow2_hash_func(hash_size_ - 1),
          MaxRoughNonZero(MaxRoughNonZero_),
          shared_memory_size(sharedMemorySize_),
          vector_size (vector_size_),
          m_space(mpool_),
          my_exec_space(my_exec_space_),

          //unit memory for a hashmap entry. assuming 1 begin, 1 next, 1 key 1 value.
          unit_memory(sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2),
          suggested_team_size(suggested_team_size_),
          thread_memory((shared_memory_size /8 / suggested_team_size_) * 8),
          shmem_key_size(),
          shared_memory_hash_func(),
          shmem_hash_size(1),
          team_row_chunk_size(team_row_chunk_size_)
    {

      //how many keys I can hold?
      //thread memory - 3 needed entry for size.
      shmem_key_size = ((thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory);

      //put the hash size closest power of 2.
      //we round down here, because we want to store more keys,
      //conflicts are cheaper.
      while (shmem_hash_size * 2 <=  shmem_key_size){
        shmem_hash_size = shmem_hash_size * 2;
      }
      //for and opeation we get -1.
      shared_memory_hash_func = shmem_hash_size - 1;

      //increase the key size wit the left over from hash size.
      shmem_key_size = shmem_key_size + ((shmem_key_size - shmem_hash_size) ) / 3;
      //round it down to 2, because of some alignment issues.
      shmem_key_size = (shmem_key_size >> 1) << 1;

      if (KOKKOSKERNELS_VERBOSE){

        std::cout << "\tStructureC "
                  << " thread_memory:" << thread_memory
                  << " unit_memory:" << unit_memory
                  << " adjusted hashsize:" << shmem_hash_size
                  << " adjusted shmem_key_size:" << shmem_key_size
                  << " using "<< (shmem_key_size * 3  + shmem_hash_size) * sizeof (nnz_lno_t) +    sizeof(nnz_lno_t) * 3
                  << " of thread_memory: " << thread_memory
                  << std::endl;
            }
    }

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_OPENMP )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreDenseAccumulatorTag&, const team_member_t & teamMember) const {
      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);

      //dense accumulators
      nnz_lno_t *indices = NULL;
      nnz_lno_t *sets = NULL;
      volatile nnz_lno_t * tmp = NULL;

      size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
      }

      //we need as much as column size for sets.
      sets = (nnz_lno_t *) tmp;
      tmp += MaxRoughNonZero; //this is set as column size before calling dense accumulators.
      //indices only needs max row size.
      indices = (nnz_lno_t *) tmp;


      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index)
      {
        nnz_lno_t index_cnt = 0;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t col_size = row_mapA[row_index + 1] - col_begin;

        //traverse columns of A
        for (nnz_lno_t colind = 0; colind < col_size; ++colind){
          size_type a_col = colind + col_begin;

          nnz_lno_t rowB = entriesA[a_col];
          size_type rowBegin = row_pointer_begins_B(rowB);
          nnz_lno_t left_work = row_pointer_ends_B(rowB ) - rowBegin;

          //traverse columns of B
          for (nnz_lno_t i = 0; i < left_work; ++i){

            const size_type adjind = i + rowBegin;
            nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
            nnz_lno_t b_set = entriesSetsB[adjind];

            //if sets are not set before, add this to indices.
            if (sets[b_set_ind] == 0){
              indices[index_cnt++] = b_set_ind;
            }
            //make a union.
            sets[b_set_ind] |= b_set;
          }
        }
        nnz_lno_t num_el = 0;
        for (nnz_lno_t ii = 0; ii < index_cnt; ++ii){
          nnz_lno_t set_ind = indices[ii];
          nnz_lno_t c_rows = sets[set_ind];
          sets[set_ind] = 0;

          //count number of set bits
          nnz_lno_t num_el2 = 0; 
          for (; c_rows; num_el2++) {
            c_rows = c_rows & (c_rows - 1); // clear the least significant bit set
          }
          num_el += num_el2;
        }
        rowmapC(row_index) = num_el;
      }
      );

      m_space.release_chunk(indices);
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {



      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_row_chunk_size;
      const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_row_chunk_size, numrows);


      //get memory from memory pool.
      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(team_row_begin + teamMember.team_rank());
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
      }

      //set first to globally used hash indices.
      nnz_lno_t *globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size;

      //create hashmap accumulator.
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

      //set memory for hash begins.
      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size ;

      hm2.hash_nexts = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;

      //holds the keys
      hm2.keys = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;
      hm2.values = (nnz_lno_t *) (tmp);

      hm2.hash_key_size = pow2_hash_size;
      hm2.max_value_size = MaxRoughNonZero;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&] (const nnz_lno_t& row_index){
        nnz_lno_t globally_used_hash_count = 0;
        nnz_lno_t used_hash_size = 0;
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t col_size = row_mapA[row_index + 1] - col_begin;
        //traverse columns of A.
        for (nnz_lno_t colind = 0; colind < col_size; ++colind){
          size_type a_col = colind + col_begin;
          nnz_lno_t rowB = entriesA[a_col];

          size_type rowBegin = row_pointer_begins_B(rowB);
          nnz_lno_t left_work = row_pointer_ends_B(rowB ) - rowBegin;
          //traverse columns of B
          for (nnz_lno_t i = 0; i < left_work; ++i){

            const size_type adjind = i + rowBegin;

            nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
            nnz_lno_t b_set = entriesSetsB[adjind];
            nnz_lno_t hash = b_set_ind & pow2_hash_func;

            //insert it to first hash.
            hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                hash,
                b_set_ind, b_set,
                &used_hash_size,
                hm2.max_value_size,&globally_used_hash_count,
                globally_used_hash_indices
            );
          }
        }

        //when done with all insertions, traverse insertions and get the size.
        nnz_lno_t num_el = 0;
        for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii){
          nnz_lno_t c_rows = hm2.values[ii];
          nnz_lno_t num_el2 = 0;

          //the number of set bits.
          for (; c_rows; num_el2++) {
            c_rows = c_rows & (c_rows - 1); // clear the least significant bit set
          }
          num_el += num_el2;
        }

        //clear the begins.
        for (int i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }
        //set the row size.
        rowmapC(row_index) = num_el;
      });

      m_space.release_chunk(globally_used_hash_indices);
    }


    KOKKOS_INLINE_FUNCTION
    void operator()(const GPUTag&, const team_member_t & teamMember) const {


      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      //int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) * 4;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //nnz_lno_t *alloc_global_memory = NULL;
      nnz_lno_t *globally_used_hash_indices = NULL;

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);

      all_shared_memory += sizeof(nnz_lno_t) ;
      //int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
      //nnz_lno_t shmem_key_size = (thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shmem_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shmem_key_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);

      //printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts, keys, vals);
      //return;
      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shmem_hash_size, shmem_key_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

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
        globally_used_hash_count[0] = 0;
      });

      bool is_global_alloced = false;

      const size_type col_end = row_mapA[row_index + 1];
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t col_size = col_end - col_begin;

      for (nnz_lno_t colind = 0; colind < col_size; ++colind){
        size_type a_col = colind + col_begin;

        nnz_lno_t rowB = entriesA[a_col];
        size_type rowBegin = row_pointer_begins_B(rowB);

        nnz_lno_t left_work = row_pointer_ends_B(rowB) - rowBegin;

        while (left_work){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

          nnz_lno_t b_set_ind = -1, b_set = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_set_ind = entriesSetIndicesB[adjind];
            b_set = entriesSetsB[adjind];
            //hash = b_set_ind % shmem_key_size;
            hash = b_set_ind & shared_memory_hash_func;
          });


          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size,
              hash, b_set_ind, b_set,
              used_hash_sizes,
              shmem_key_size);


          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int threadid, int &overall_num_unsuccess) {
            overall_num_unsuccess += num_unsuccess;
          }, overall_num_unsuccess);


          if (overall_num_unsuccess){

            //printf("row:%d\n", row_index);
            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              size_t tid = get_thread_id(row_index);
              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
                }, tmp);
              }
              is_global_alloced = true;

              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;

              //poins to the next elements
              hm2.hash_nexts = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;

              //holds the keys
              hm2.keys = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;
              hm2.values = (nnz_lno_t *) (tmp);

              hm2.hash_key_size = pow2_hash_size;
              hm2.max_value_size = MaxRoughNonZero;
            }

            nnz_lno_t hash = -1;
            if (num_unsuccess) hash = b_set_ind & pow2_hash_func;

            //int insertion =
            hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                teamMember, vector_size,
                hash,b_set_ind,b_set,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
                );

          }
          left_work -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shmem_key_size) used_hash_sizes[0] = shmem_key_size;
      });

      /*
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[1] > hm2.max_value_size) used_hash_sizes[1] = hm2.max_value_size;
      });
      */

      nnz_lno_t num_elements = 0;

      nnz_lno_t num_compressed_elements = used_hash_sizes[0];

      Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
          [&] (const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
        nnz_lno_t c_rows = hm.values[ii];
        nnz_lno_t num_el = 0;
        for (; c_rows; num_el++) {
          c_rows &= c_rows - 1; // clear the least significant bit set
        }
        num_nnz_in_row += num_el;
      }, num_elements);


      if (is_global_alloced){
        nnz_lno_t num_global_elements = 0;
        nnz_lno_t num_compressed_elements = used_hash_sizes[1];
        Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
            [&] (const nnz_lno_t ii, nnz_lno_t &num_nnz_in_row) {
          nnz_lno_t c_rows = hm2.values[ii];
          nnz_lno_t num_el = 0;
          for (; c_rows; num_el++) {
            c_rows &= c_rows - 1; // clear the least significant bit set
          }
          num_nnz_in_row += num_el;
        }, num_global_elements);


        //now thread leaves the memory as it finds. so there is no need to initialize the hash begins
        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });


        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          m_space.release_chunk(globally_used_hash_indices);
        });
        num_elements += num_global_elements;
      }

      rowmapC(row_index) = num_elements;
    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };



  template <typename a_row_view_t, typename a_nnz_view_t,
            typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t, typename nnz_lno_temp_work_view_t,
            typename pool_memory_space>
  struct NonzeroesC{
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
    const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;

    /**
     * \brief Constructor.
     */

    NonzeroesC(
        nnz_lno_t m_,
        a_row_view_t row_mapA_,
        a_nnz_view_t entriesA_,

        b_original_row_view_t old_row_mapB_,
        b_compressed_row_view_t row_mapB_,
        b_nnz_view_t entriesSetIndicesB_,
        b_nnz_view_t entriesSetsB_,

        c_row_view_t rowmapC_,
        nnz_lno_temp_work_view_t entriesSetIndicesC_,

        const nnz_lno_t hash_size_,
        const nnz_lno_t MaxRoughNonZero_,
        const size_t sharedMemorySize_,
        const int vector_size_,
        pool_memory_space mpool_,
        const KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space_):
          numrows(m_),

          row_mapA (row_mapA_),
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
          vector_size (vector_size_), m_space(mpool_), my_exec_space(my_exec_space_)
          {}

    KOKKOS_INLINE_FUNCTION
    size_t get_thread_id(const size_t row_index) const{
      switch (my_exec_space){
      default:
        return row_index;
#if defined( KOKKOS_HAVE_SERIAL )
      case KokkosKernels::Experimental::Util::Exec_SERIAL:
        return 0;
#endif
#if defined( KOKKOS_HAVE_OPENMP )
      case KokkosKernels::Experimental::Util::Exec_OMP:
        return Kokkos::OpenMP::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
      case KokkosKernels::Experimental::Util::Exec_PTHREADS:
        return Kokkos::Threads::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_QTHREAD)
      case KokkosKernels::Experimental::Util::Exec_QTHREADS:
        return Kokkos::Qthread::hardware_thread_id();
#endif
#if defined( KOKKOS_HAVE_CUDA )
      case KokkosKernels::Experimental::Util::Exec_CUDA:
        return row_index;
#endif
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MultiCoreTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;
      //get row index.

      nnz_lno_t *globally_used_hash_indices = NULL;
      nnz_lno_t globally_used_hash_count = 0;
      nnz_lno_t used_hash_size = 0;
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

      volatile nnz_lno_t * tmp = NULL;
      size_t tid = get_thread_id(row_index);
      while (tmp == NULL){
        tmp = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
      }

      globally_used_hash_indices = (nnz_lno_t *) tmp;
      tmp += pow2_hash_size ;

      hm2.hash_begins = (nnz_lno_t *) (tmp);
      tmp += pow2_hash_size ;

      //poins to the next elements
      hm2.hash_nexts = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;

      //holds the keys
      hm2.keys = (nnz_lno_t *) (tmp);
      tmp += MaxRoughNonZero;
      hm2.values = (nnz_lno_t *) (tmp);

      hm2.hash_key_size = pow2_hash_size;
      hm2.max_value_size = MaxRoughNonZero;


      {
        const size_type col_begin = row_mapA[row_index];
        const nnz_lno_t col_size = row_mapA[row_index + 1] - col_begin;

        for (nnz_lno_t colind = 0; colind < col_size; ++colind){
          size_type a_col = colind + col_begin;

          nnz_lno_t rowB = entriesA[a_col];
          size_type rowBegin = old_row_mapB(rowB);

          nnz_lno_t left_work = row_mapB(rowB ) - rowBegin;

          for (nnz_lno_t i = 0; i < left_work; ++i){

            const size_type adjind = i + rowBegin;
            nnz_lno_t b_set_ind = entriesSetIndicesB[adjind];
            nnz_lno_t b_set = entriesSetsB[adjind];
            nnz_lno_t hash = b_set_ind & pow2_hash_func;

            hm2.sequential_insert_into_hash_mergeOr_TrackHashes(
                hash,b_set_ind,b_set,
                &used_hash_size, hm2.max_value_size
                ,&globally_used_hash_count, globally_used_hash_indices
            );
          }
        }

        int set_size = sizeof(nnz_lno_t) * 8;
        nnz_lno_t num_el = rowmapC(row_index);
        for (nnz_lno_t ii = 0; ii < used_hash_size; ++ii){
          nnz_lno_t c_rows_setind = hm2.keys[ii];
          nnz_lno_t c_rows = hm2.values[ii];

          int current_row = 0;
          nnz_lno_t unit = 1;

          while (c_rows){
            if (c_rows & unit){
              //insert indices.
              entriesSetIndicesC(num_el++) = set_size * c_rows_setind + current_row;
            }
            current_row++;
            c_rows = c_rows & ~unit;
            unit = unit << 1;
          }
        }
        for (int i = 0; i < globally_used_hash_count; ++i){
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        }
      }
      m_space.release_chunk(globally_used_hash_indices);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const GPUTag&, const team_member_t & teamMember) const {

      nnz_lno_t row_index = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (row_index >= numrows) return;


      //printf("row:%d\n", row_index);

      int thread_memory = ((shared_memory_size/ 4 / teamMember.team_size())) * 4;
      char *all_shared_memory = (char *) (teamMember.team_shmem().get_shmem(shared_memory_size));

      //nnz_lno_t *alloc_global_memory = NULL;
      nnz_lno_t *globally_used_hash_indices = NULL;

      //shift it to the thread private part
      all_shared_memory += thread_memory * teamMember.team_rank();

      //used_hash_sizes hold the size of 1st and 2nd level hashes
      volatile nnz_lno_t *used_hash_sizes = (volatile nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * 2;

      nnz_lno_t *globally_used_hash_count = (nnz_lno_t *) (all_shared_memory);

      all_shared_memory += sizeof(nnz_lno_t) ;
      int unit_memory = sizeof(nnz_lno_t) * 2 + sizeof(nnz_lno_t) * 2;
      nnz_lno_t shared_memory_hash_size = (thread_memory - sizeof(nnz_lno_t) * 3) / unit_memory;

      nnz_lno_t * begins = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //poins to the next elements
      nnz_lno_t * nexts = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;

      //holds the keys
      nnz_lno_t * keys = (nnz_lno_t *) (all_shared_memory);
      all_shared_memory += sizeof(nnz_lno_t) * shared_memory_hash_size;
      nnz_lno_t * vals = (nnz_lno_t *) (all_shared_memory);

      //printf("begins:%ld, nexts:%ld, keys:%ld, vals:%ld\n", begins, nexts, keys, vals);
      //return;
      //first level hashmap
      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t>
        hm(shared_memory_hash_size, shared_memory_hash_size, begins, nexts, keys, vals);

      KokkosKernels::Experimental::UnorderedHashmap::HashmapAccumulator<nnz_lno_t,nnz_lno_t,nnz_lno_t> hm2;

      //initialize begins.
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, shared_memory_hash_size),
          [&] (int i) {
        begins[i] = -1;
      });

      //initialize hash usage sizes
      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        used_hash_sizes[0] = 0;
        used_hash_sizes[1] = 0;
        globally_used_hash_count[0] = 0;
      });

      bool is_global_alloced = false;

      const size_type col_end = row_mapA[row_index + 1];
      const size_type col_begin = row_mapA[row_index];
      const nnz_lno_t col_size = col_end - col_begin;

      for (nnz_lno_t colind = 0; colind < col_size; ++colind){
        size_type a_col = colind + col_begin;

        nnz_lno_t rowB = entriesA[a_col];
        size_type rowBegin = old_row_mapB(rowB);

        nnz_lno_t left_work = row_mapB(rowB ) - rowBegin;

        while (left_work){
          nnz_lno_t work_to_handle = KOKKOSKERNELS_MACRO_MIN(vector_size, left_work);

          nnz_lno_t b_set_ind = -1, b_set = -1;
          nnz_lno_t hash = -1;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(teamMember, work_to_handle),
              [&] (nnz_lno_t i) {
            const size_type adjind = i + rowBegin;
            b_set_ind = entriesSetIndicesB[adjind];
            b_set = entriesSetsB[adjind];
            hash = b_set_ind % shared_memory_hash_size;
          });


          int num_unsuccess = hm.vector_atomic_insert_into_hash_mergeOr(
              teamMember, vector_size,
              hash, b_set_ind, b_set,
              used_hash_sizes,
              shared_memory_hash_size);


          int overall_num_unsuccess = 0;

          Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, vector_size),
              [&] (const int threadid, int &overall_num_unsuccess) {
            overall_num_unsuccess += num_unsuccess;
          }, overall_num_unsuccess);


          if (overall_num_unsuccess){

            //printf("row:%d\n", row_index);
            if (!is_global_alloced){
              volatile nnz_lno_t * tmp = NULL;
              size_t tid = get_thread_id(row_index);
              while (tmp == NULL){
                Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile nnz_lno_t * &memptr) {
                  memptr = (volatile nnz_lno_t * )( m_space.allocate_chunk(tid));
                }, tmp);
              }
              is_global_alloced = true;

              globally_used_hash_indices = (nnz_lno_t *) tmp;
              tmp += pow2_hash_size ;

              hm2.hash_begins = (nnz_lno_t *) (tmp);
              tmp += pow2_hash_size ;

              //poins to the next elements
              hm2.hash_nexts = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;

              //holds the keys
              hm2.keys = (nnz_lno_t *) (tmp);
              tmp += MaxRoughNonZero;
              hm2.values = (nnz_lno_t *) (tmp);

              hm2.hash_key_size = pow2_hash_size;
              hm2.max_value_size = MaxRoughNonZero;
            }

            nnz_lno_t hash = -1;
            if (num_unsuccess) hash = b_set_ind & pow2_hash_func;

            //int insertion =
            hm2.vector_atomic_insert_into_hash_mergeOr_TrackHashes(
                teamMember, vector_size,
                hash,b_set_ind,b_set,
                used_hash_sizes + 1, hm2.max_value_size
                ,globally_used_hash_count, globally_used_hash_indices
                );

          }
          left_work -= work_to_handle;
          rowBegin += work_to_handle;
        }
      }

      Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
        if (used_hash_sizes[0] > shared_memory_hash_size) used_hash_sizes[0] = shared_memory_hash_size;
      });

      nnz_lno_t num_compressed_elements = used_hash_sizes[0];
      used_hash_sizes[0] = 0;
      size_type row_begin = rowmapC(row_index);
      int set_size = sizeof(nnz_lno_t) * 8;
      Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
          [&] (const nnz_lno_t ii) {
        nnz_lno_t c_rows_setind = hm.keys[ii];
        nnz_lno_t c_rows = hm.values[ii];

        int current_row = 0;
        nnz_lno_t unit = 1;

        while (c_rows){
          if (c_rows & unit){

            size_type wind = Kokkos::atomic_fetch_add(used_hash_sizes, 1);
            entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
          }
          current_row++;
          c_rows = c_rows & ~unit;
          unit = unit << 1;
        }

      });


      if (is_global_alloced){
        nnz_lno_t num_compressed_elements = used_hash_sizes[1];
        Kokkos::parallel_for( Kokkos::ThreadVectorRange(teamMember, num_compressed_elements),
            [&] (const nnz_lno_t ii) {
          nnz_lno_t c_rows_setind = hm2.keys[ii];
          nnz_lno_t c_rows = hm2.values[ii];

          int current_row = 0;
          nnz_lno_t unit = 1;

          while (c_rows){
            if (c_rows & unit){

              size_type wind = Kokkos::atomic_fetch_add(used_hash_sizes, 1);
              entriesSetIndicesC(wind + row_begin) = set_size * c_rows_setind + current_row;
            }
            current_row++;
            c_rows = c_rows & ~unit;
            unit = unit << 1;
          }
        });

        //now thread leaves the memory as it finds. so there is no need to initialize the hash begins
        nnz_lno_t dirty_hashes = globally_used_hash_count[0];
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(teamMember, dirty_hashes),
            [&] (nnz_lno_t i) {
          nnz_lno_t dirty_hash = globally_used_hash_indices[i];
          hm2.hash_begins[dirty_hash] = -1;
        });


        Kokkos::single(Kokkos::PerThread(teamMember),[&] () {
          m_space.release_chunk(globally_used_hash_indices);
        });
      }

    }

    size_t team_shmem_size (int team_size) const {
      return shared_memory_size;
    }

  };


  template <typename a_row_view_t, typename a_nnz_view_t,
              typename b_original_row_view_t,
              typename b_compressed_row_view_t, typename b_nnz_view_t,
              typename c_row_view_t>
  void symbolic_c(
      nnz_lno_t m,
      a_row_view_t row_mapA,
      a_nnz_view_t entriesA,

      b_original_row_view_t old_row_mapB,
      b_compressed_row_view_t row_mapB,
      b_nnz_view_t entriesSetIndex,
      b_nnz_view_t entriesSets,

      c_row_view_t rowmapC,
      nnz_lno_t maxNumRoughNonzeros
  ){
    typedef KokkosKernels::Experimental::Util::UniformMemoryPool< MyTempMemorySpace, nnz_lno_t> pool_memory_space;

    //get the number of rows and nonzeroes of B.
    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  entriesSetIndex.dimension_0();

    //get the SPGEMMAlgorithm to run.
    SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();

    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = this->handle->get_handle_exec_space();
    int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

    //round up maxNumRoughNonzeros to closest power of 2.
    nnz_lno_t min_hash_size = 1;
    while (maxNumRoughNonzeros > min_hash_size){
      min_hash_size *= 2;
    }

    //set the chunksize.
    size_t chunksize = min_hash_size ; //this is for used hash indices
    chunksize += min_hash_size ; //this is for the hash begins
    chunksize += maxNumRoughNonzeros ; //this is for hash nexts
    chunksize += maxNumRoughNonzeros ; //this is for hash keys
    chunksize += maxNumRoughNonzeros ; //this is for hash values

    //initizalize value for the mem pool
    int pool_init_val = -1;

    //if KKSPEED are used on CPU, or KKMEMSPEED is run with threads less than 32
    //than we use dense accumulators.
    if ((   spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMSPEED  &&
        concurrency <=  sizeof (nnz_lno_t) * 8 &&
        my_exec_space != KokkosKernels::Experimental::Util::Exec_CUDA)
        ||
        (   spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED &&
            my_exec_space != KokkosKernels::Experimental::Util::Exec_CUDA)){

      nnz_lno_t col_size = this->b_col_cnt / (sizeof (nnz_lno_t) * 8)+ 1;

      nnz_lno_t max_row_size = KOKKOSKERNELS_MACRO_MIN(col_size, maxNumRoughNonzeros);
      chunksize = col_size + max_row_size;
      //if speed is set, and exec space is cpu, then  we use dense accumulators.
      //or if memspeed is set, and concurrency is not high, we use dense accumulators.
      maxNumRoughNonzeros = col_size;
      pool_init_val = 0;
      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\tDense Acc - COLS:" << col_size << " max_row_size:" << max_row_size << std::endl;
      }
    }


    size_t num_chunks = concurrency / suggested_vector_size;

    KokkosKernels::Experimental::Util::PoolType my_pool_type = KokkosKernels::Experimental::Util::OneThread2OneChunk;
    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA) {
      my_pool_type = KokkosKernels::Experimental::Util::ManyThread2OneChunk;
    }

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tPool Size (MB):" << (num_chunks * chunksize * sizeof(nnz_lno_t)) / 1024. / 1024.  << std::endl;
    }

    Kokkos::Impl::Timer timer1;
    pool_memory_space m_space(num_chunks, chunksize, pool_init_val,  my_pool_type);
    MyExecSpace::fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tPool Alloc Time:" << timer1.seconds() << std::endl;
    }

    StructureC<a_row_view_t, a_nnz_view_t,
    b_original_row_view_t, b_compressed_row_view_t, b_nnz_view_t,
    c_row_view_t, /* nnz_lno_temp_work_view_t,*/ pool_memory_space>
    sc(
        m,
        row_mapA,
        entriesA,
        old_row_mapB,
        row_mapB,
        entriesSetIndex,
        entriesSets,
        rowmapC,
        min_hash_size,
        maxNumRoughNonzeros,
        shmem_size,
        suggested_team_size,
        team_row_chunk_size,
        suggested_vector_size,
        m_space,
        my_exec_space,KOKKOSKERNELS_VERBOSE);

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tStructureC vector_size:" << suggested_vector_size
          << " team_size:" << suggested_team_size
          << " chunk_size:" << team_row_chunk_size
          << " shmem_size:" << shmem_size << std::endl;
    }

    timer1.reset();

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA) {
      Kokkos::parallel_for( gpu_team_policy_t(m / suggested_team_size + 1 , suggested_team_size, suggested_vector_size), sc);
    }
    else {
      if (( spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMSPEED  &&
          concurrency <=  sizeof (nnz_lno_t) * 8)  ||
          spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED){

        if (use_dynamic_schedule){
          Kokkos::parallel_for( dynamic_multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
        else {
          Kokkos::parallel_for( multicore_dense_team_count_policy_t(m / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
      }
      else {
        if (use_dynamic_schedule){
          Kokkos::parallel_for( dynamic_multicore_team_policy_t(m / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
        else {
          Kokkos::parallel_for( multicore_team_policy_t(m / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
      }
    }
    MyExecSpace::fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tStructureC Kernel time:" << timer1.seconds() << std::endl<< std::endl;
    }
    //if we need to find the max nnz in a row.
    {
      Kokkos::Impl::Timer timer1;
      size_type c_max_nnz = 0;
      KokkosKernels::Experimental::Util::view_reduce_max<c_row_view_t, MyExecSpace>(m, rowmapC, c_max_nnz);
      MyExecSpace::fence();
      this->handle->get_spgemm_handle()->set_max_result_nnz(c_max_nnz);

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\tReduce Max Row Size Time:" << timer1.seconds() << std::endl;
      }
    }

    KokkosKernels::Experimental::Util::kk_exclusive_parallel_prefix_sum<row_lno_temp_work_view_t, MyExecSpace>(m+1, rowmapC);
    MyExecSpace::fence();


    auto d_c_nnz_size = Kokkos::subview(rowmapC, m);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    typename c_row_view_t::non_const_value_type c_nnz_size = h_c_nnz_size();
    this->handle->get_spgemm_handle()->set_c_nnz(c_nnz_size);


    if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
        spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
        spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\tCOLORING PHASE"<<  std::endl;
      }
      nnz_lno_temp_work_view_t entryIndicesC_ (Kokkos::ViewAllocateWithoutInitializing("entryIndicesC_"), c_nnz_size);

      //calculate the structure.
      NonzeroesC<
      a_row_view_t, a_nnz_view_t,
      b_original_row_view_t, b_compressed_row_view_t, b_nnz_view_t,
      c_row_view_t, nnz_lno_temp_work_view_t,
      pool_memory_space>
      sc( m,
          row_mapA,
          entriesA,
          old_row_mapB,
          row_mapB,
          entriesSetIndex,
          entriesSets,
          rowmapC,
          entryIndicesC_,
          min_hash_size,
          maxNumRoughNonzeros,
          shmem_size,suggested_vector_size,m_space,
          my_exec_space);

      timer1.reset();
      if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA) {
        Kokkos::parallel_for( gpu_team_policy_t(m / suggested_team_size + 1 , suggested_team_size, suggested_vector_size), sc);
      }
      else {
        if (use_dynamic_schedule){
          Kokkos::parallel_for( dynamic_multicore_team_policy_t(m / suggested_team_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
        else {
          Kokkos::parallel_for( multicore_team_policy_t(m / suggested_team_size + 1 , suggested_team_size, suggested_vector_size), sc);
        }
      }

      MyExecSpace::fence();


      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tCOLORING-NNZ-FILL-TIME:" << timer1.seconds() <<  std::endl;
      }

      nnz_lno_t original_num_colors, num_colors_in_one_step, num_multi_color_steps;
      nnz_lno_persistent_work_host_view_t h_color_xadj;
      nnz_lno_persistent_work_view_t color_adj, vertex_colors_to_store;

      //distance-2 color
      this->d2_color_c_matrix(
          rowmapC, entryIndicesC_,
          original_num_colors, h_color_xadj, color_adj , vertex_colors_to_store,
          num_colors_in_one_step, num_multi_color_steps, spgemm_algorithm);

      timer1.reset();

      //sort the color indices.
      for (nnz_lno_t i = 0; i < num_multi_color_steps; ++i){
        //sort the ones that have more than 32 rows.
        if (h_color_xadj(i+1) - h_color_xadj(i) <= 32) continue;
        auto sv = Kokkos::subview(color_adj,Kokkos::pair<nnz_lno_t, nnz_lno_t> (h_color_xadj(i), h_color_xadj(i+1)));
        //KokkosKernels::Experimental::Util::print_1Dview(sv, i ==47);
        //TODO for some reason kokkos::sort is failing on views with size 56 and 112.
        //for now we use std::sort. Delete below comment, and delete the code upto fence.
        //Kokkos::sort(sv);
        //
        auto h_sv = Kokkos::create_mirror_view (sv);
        Kokkos::deep_copy(h_sv,sv);
        auto* p_sv = h_sv.ptr_on_device();
        std::sort (p_sv, p_sv + h_color_xadj(i+1) - h_color_xadj(i));
        Kokkos::deep_copy(sv,h_sv);
        MyExecSpace::fence();
      }

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tCOLOR-SORT-TIME:" << timer1.seconds() <<  std::endl;
      }
      this->handle->get_spgemm_handle()->set_color_xadj(
          original_num_colors,
          h_color_xadj, color_adj, vertex_colors_to_store,
          num_colors_in_one_step, num_multi_color_steps);
      this->handle->get_spgemm_handle()->set_c_column_indices(entryIndicesC_);
    }

  }


  /**
   * \brief Numeric phase with speed method
   */
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_color(
      c_row_view_t rowmapC_,
      c_lno_nnz_view_t entriesC_,
      c_scalar_nnz_view_t valuesC_,
      SPGEMMAlgorithm spgemm_algorithm){

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tCOLOR MODE" << std::endl;
    }
    nnz_lno_temp_work_view_t entryIndicesC_ =
        this->handle->get_spgemm_handle()->get_c_column_indices();

    KokkosKernels::Experimental::Util::kk_copy_vector
      <nnz_lno_temp_work_view_t, c_lno_nnz_view_t, MyExecSpace>
        (entryIndicesC_.dimension_0(), entryIndicesC_, entriesC_);

    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space =
        KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();


    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();
    //get vector size, team size.
    int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);

    //get color vertices
    nnz_lno_t num_colors, num_multi_colors, num_used_colors;
    nnz_lno_persistent_work_host_view_t color_xadj;
    nnz_lno_persistent_work_view_t color_adj, vertex_colors;
    this->handle->get_spgemm_handle()->get_color_xadj(
        num_colors, color_xadj, color_adj, vertex_colors, num_multi_colors, num_used_colors);

    const nnz_lno_t block_size = 64;
    const nnz_lno_t shift_divisor = 6;
    scalar_temp_work_view_t denseAccumulator ("Scalar Accumulator", ( (this->b_col_cnt + block_size) * num_multi_colors));


    //bool_temp_view_t denseAccumulatorFlags ("Accumulator flags", ((this->k* 1.5)  * num_multi_colors));


    NumericCCOLOR<
      const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t>
      sc(
        a_row_cnt, b_col_cnt,
        row_mapA,
        entriesA,
        valsA,

        row_mapB,
        entriesB,
        valsB,

        rowmapC_,
        entriesC_,
        valuesC_,

        denseAccumulator,
        //denseAccumulatorFlags,
        -1);




    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCOLORING-num_multi_colors:" << num_multi_colors << " num_used_colors:" << num_used_colors << std::endl;
    }
    sc.color_adj = color_adj;
    sc.vertex_colors = vertex_colors;
    sc.chunk_divison = shift_divisor;
    sc.chunk_and = block_size - 1;
    sc.consecutive_chunk_size = block_size;
    sc.consecutive_all_color_chunk_size = sc.consecutive_chunk_size * num_multi_colors;

    Kokkos::Impl::Timer timer1;
    for (nnz_lno_t i = 0; i < num_used_colors; ){
      nnz_lno_t color_begin = color_xadj(i);
      nnz_lno_t lastcolor = i + 1;
      if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
        lastcolor = KOKKOSKERNELS_MACRO_MIN(i + num_multi_colors, num_used_colors );
        i += num_multi_colors;
      }
      else {
        ++i;
      }

      nnz_lno_t color_end = color_xadj(lastcolor);
      sc.color_begin = color_begin;
      sc.color_end = color_end;

      nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, color_end - color_begin);
      sc.team_work_size = team_row_chunk_size;


      if (use_dynamic_schedule){
        switch (spgemm_algorithm){
        default:
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR:
          Kokkos::parallel_for( dynamic_team_numeric1_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2:
          Kokkos::parallel_for( dynamic_team_numeric2_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR:
          Kokkos::parallel_for( dynamic_team_numeric3_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        }
      }
      else {
        switch (spgemm_algorithm){
        default:
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR:
          Kokkos::parallel_for( team_numeric1_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2:
          Kokkos::parallel_for( team_numeric2_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        case KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR:
          Kokkos::parallel_for( team_numeric3_policy_t((color_end - color_begin) / team_row_chunk_size + 1 ,
              suggested_team_size, suggested_vector_size), sc);
          break;
        }
      }
      MyExecSpace::fence();
    }

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
    }
  }

  /**
   * \brief Numeric phase with speed method
   */
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_speed(
      c_row_view_t rowmapC_,
      c_lno_nnz_view_t entriesC_,
      c_scalar_nnz_view_t valuesC_,
      KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space){

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tSPEED MODE" << std::endl;
    }

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();

    //get suggested vector size, teamsize and row chunk size.
    int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      //allocate memory for begins and next to be used by the hashmap
      c_lno_nnz_view_t beginsC
        (Kokkos::ViewAllocateWithoutInitializing("C keys"), valuesC_.dimension_0());
      c_lno_nnz_view_t nextsC
        (Kokkos::ViewAllocateWithoutInitializing("C nexts"), valuesC_.dimension_0());
      Kokkos::deep_copy(beginsC, -1);

      //create the functor.
      NumericCMEM<
      const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
      const_b_lno_row_view_t, const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t,
      c_row_view_t, c_lno_nnz_view_t, c_scalar_nnz_view_t>
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
      MyExecSpace::fence();

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tGPU vector_size:" << suggested_vector_size
                  <<  " team_size:" << suggested_team_size
                  << " chunk_size:" << team_row_chunk_size
                  << std::endl;
      }

      timer1.reset();
      Kokkos::parallel_for(
          gpu_team_policy_t(
              a_row_cnt / team_row_chunk_size + 1 ,
              suggested_team_size ,
              suggested_vector_size),
              sc);
      MyExecSpace::fence();

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
      }
    }
    else {
      typedef KokkosKernels::Experimental::Util::UniformMemoryPool
            < MyTempMemorySpace, scalar_t> pool_memory_space;


      KokkosKernels::Experimental::Util::PoolType my_pool_type =
          KokkosKernels::Experimental::Util::OneThread2OneChunk;
      int num_chunks = concurrency;

      Kokkos::Impl::Timer timer1;
      pool_memory_space m_space
        (num_chunks, this->b_col_cnt + (this->b_col_cnt) / sizeof(scalar_t) + 1, 0,  my_pool_type);
      MyExecSpace::fence();

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
          my_exec_space,
          team_row_chunk_size);

      MyExecSpace::fence();
      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tCPU vector_size:" << suggested_vector_size
                  <<  " team_size:" << suggested_team_size
                  << " chunk_size:" << team_row_chunk_size
                  << std::endl;
      }
      timer1.reset();

      if (use_dynamic_schedule){
        Kokkos::parallel_for( dynamic_multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
      }
      else {
        Kokkos::parallel_for( multicore_team_policy_t(a_row_cnt / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sc);
      }

      MyExecSpace::fence();

      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\t\tNumeric TIME:" << timer1.seconds() << std::endl;
      }
    }
  }

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_hash(
        c_row_view_t rowmapC_,
        c_lno_nnz_view_t entriesC_,
        c_scalar_nnz_view_t valuesC_,
        KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space){

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tHASH MODE" << std::endl;
    }

    nnz_lno_t brows = row_mapB.dimension_0() - 1;
    size_type bnnz =  valsB.dimension_0();

    int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size,concurrency, a_row_cnt);

    typedef KokkosKernels::Experimental::Util::UniformMemoryPool< MyTempMemorySpace, nnz_lno_t> pool_memory_space;


    nnz_lno_t max_nnz = this->handle->get_spgemm_handle()->get_max_result_nnz();
    nnz_lno_t min_hash_size = 1;
    while (max_nnz > min_hash_size){
      min_hash_size *= 4;
    }
   
    size_t chunksize = min_hash_size; //this is for used hash indices
    chunksize += min_hash_size ; //this is for the hash begins
    chunksize += max_nnz; //this is for hash nexts
    int num_chunks = concurrency / suggested_vector_size;

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\t max_nnz: " << max_nnz
                << " chunk_size:" << chunksize
                << " numchunks:" << num_chunks << std::endl;
    }

    KokkosKernels::Experimental::Util::PoolType my_pool_type =
        KokkosKernels::Experimental::Util::OneThread2OneChunk;
    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      my_pool_type = KokkosKernels::Experimental::Util::ManyThread2OneChunk;
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
        min_hash_size,
        suggested_team_size,

        my_exec_space,
        team_row_chunk_size,KOKKOSKERNELS_VERBOSE);


    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tvector_size:" << suggested_vector_size  << " chunk_size:" << team_row_chunk_size << std::endl;
    }
    timer1.reset();

    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
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

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric(
      c_row_view_t rowmapC_,
      c_lno_nnz_view_t entriesC_,
      c_scalar_nnz_view_t valuesC_){

    //get the algorithm and execution space.
    SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = KokkosKernels::Experimental::Util::get_exec_space_type<MyExecSpace>();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "Numeric PHASE" << std::endl;
    }

    if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED)
    {
      this->KokkosSPGEMM_numeric_speed(rowmapC_, entriesC_, valuesC_, my_exec_space);

    }
    else if ( spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR ||
              spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR ||
              spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
      this->KokkosSPGEMM_numeric_color(rowmapC_, entriesC_, valuesC_, spgemm_algorithm);
    }
    else {
      this->KokkosSPGEMM_numeric_hash(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }

  }


  /**
   * \brief Symbolic phase of the SPGEMM.
   * \param rowmapC_: row pointers for the result matrix. Allocated before the call with size (n+1),
   * where n is the number of rows of first matrix.
   */
  template <typename c_row_view_t>
  void KokkosSPGEMM_symbolic(c_row_view_t rowmapC_){
    //number of rows and nnzs
    nnz_lno_t n = this->row_mapB.dimension_0() - 1;
    size_type nnz = this->entriesB.dimension_0();
    //compressed b
    row_lno_temp_work_view_t new_row_mapB(Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);
    nnz_lno_temp_work_view_t set_index_entries(Kokkos::ViewAllocateWithoutInitializing("set_entries_"), nnz);
    nnz_lno_temp_work_view_t set_entries(Kokkos::ViewAllocateWithoutInitializing("set_indices_"), nnz);
    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "SYMBOLIC PHASE" << std::endl;
    }
    //First Compress B.
    Kokkos::Impl::Timer timer1;

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tCOMPRESS MATRIX-B PHASE" << std::endl;
    }
    //get the compressed matrix.
    this->compressMatrix(n, nnz, this->row_mapB, this->entriesB, new_row_mapB, set_index_entries, set_entries);

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCOMPRESS MATRIX-B overall time:" << timer1.seconds()
                    << std::endl << std::endl;
    }

    timer1.reset();

    //first get the max flops for a row, which will be used for max row size.
    nnz_lno_t maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, row_mapB, new_row_mapB);

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tMax Row Flops:" << maxNumRoughZeros  << std::endl;
      std::cout << "\tMax Row Flop Calc Time:" << timer1.seconds()  << std::endl;
    }

    //calling symbolic structure
    this->symbolic_c(a_row_cnt, row_mapA, entriesA,
                    row_mapB, new_row_mapB, set_index_entries, set_entries,
                     rowmapC_, maxNumRoughZeros);
  }

  /**
   * \brief Given a symbolic matrix (a graph), it compresses the graph using bits.
   * \param in_row_map: input row pointers.
   * \param in_entries: input column entries
   * \param out_row_map: output row pointers of the compressed matrix
   * \param out_nnz_indices: output, column set indices of the output matrix.
   * \param out_nnz_sets: output, column sets of the output matrix.
   *
   */
  template <typename in_row_view_t, typename in_nnz_view_t, typename out_rowmap_view_t, typename out_nnz_view_t>
  void compressMatrix(
      nnz_lno_t n, size_type nnz,
      in_row_view_t in_row_map,
      in_nnz_view_t in_entries,

      out_rowmap_view_t out_row_map,
      out_nnz_view_t out_nnz_indices,
      out_nnz_view_t out_nnz_sets){



    //get the execution space type.
    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = this->handle->get_handle_exec_space();

    //get the suggested vectorlane size based on the execution space, and average number of nnzs per row.
    int suggested_vector_size = this->handle->get_suggested_vector_size(n, nnz);

    //get the suggested team size.
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);

    //get the chunk size suggested by the handle.
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency , n);

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tn:" << n << " nnz:" << nnz
                << " vector_size:" << suggested_vector_size
                << " team_size:" << suggested_team_size
                << " chunk_size::" << team_row_chunk_size
                << " shmem:" << shmem_size
                << std::endl;
    }


    //size of the lno_t, how many bits it can hold.
    int lnot_size = sizeof(nnz_lno_t) * 8;
    int compression_bit_divide_shift_ = 0;
    int val = lnot_size;
    while (val > 1) {
      ++compression_bit_divide_shift_;
      val = val >> 1;
    }
    nnz_lno_t compression_bit_mask_ = lnot_size - 1;

    /*
    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tCOMPRESS MATRIX: B lnot_size:" << lnot_size
                << " compression_bit_divide_shift_:" << compression_bit_divide_shift_
                << " compression_bit_mask_:" << compression_bit_mask_
                << std::endl;
    }
    */

    Kokkos::Impl::Timer timer1;


    //Allocate memory for the linked list to be used for the hashmap
    out_nnz_view_t set_nexts_;
    out_nnz_view_t set_begins_;
    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      set_nexts_ = out_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("set_nexts_"), nnz);
      set_begins_ = out_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("set_begins_"), nnz);
      Kokkos::deep_copy (set_begins_, -1);
    }
    MyExecSpace::fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCompression Allocations:" <<  timer1.seconds() << std::endl;
    }

    //create functor to compress matrix.
    SingleStepZipMatrix <in_row_view_t, in_nnz_view_t, out_rowmap_view_t, out_nnz_view_t>
      sszm_compressMatrix(
        in_row_map, //input symbolic matrix
        in_entries, //input symbolic matrix

        compression_bit_mask_, //modular
        compression_bit_divide_shift_,  //divisor
        suggested_vector_size, //vector size.
        out_row_map, //output row map
        set_begins_, //linked list parallel array begin
        set_nexts_, //linked list parallel array next
        out_nnz_indices, //output set bits
        out_nnz_sets, //output set indices
        shmem_size, //shared memory size.
        team_row_chunk_size //chunksize.
    );

    timer1.reset();
    if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
      //call GPU kernel
      //TODO: GPU use memory pool instead.
      //move calculations out of loop, to constructor.
      Kokkos::parallel_for( gpu_team_policy_t(n / suggested_team_size + 1 , suggested_team_size, suggested_vector_size), sszm_compressMatrix);
    }
    else {
      if (use_dynamic_schedule){
        //call cpu kernel with dynamic schedule
        Kokkos::parallel_for( dynamic_multicore_team_policy_t(n / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sszm_compressMatrix);
      }
      else {
        //call cpu kernel with static schedule
        Kokkos::parallel_for( multicore_team_policy_t(n / team_row_chunk_size + 1 , suggested_team_size, suggested_vector_size), sszm_compressMatrix);
      }

    }
    MyExecSpace::fence();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCompression Kernel time:" <<  timer1.seconds() << std::endl;
    }
  }

  template <typename c_row_view_t, typename c_nnz_view_t>
  void d2_color_c_matrix(
      c_row_view_t rowmapC,
      c_nnz_view_t entryIndicesC_,

      nnz_lno_t &original_num_colors,
      nnz_lno_persistent_work_host_view_t &h_color_xadj,
      nnz_lno_persistent_work_view_t &color_adj,
      nnz_lno_persistent_work_view_t &vertex_colors_to_store,

      nnz_lno_t &num_colors_in_one_step,
      nnz_lno_t &num_multi_color_steps,
      SPGEMMAlgorithm spgemm_algorithm){

    nnz_lno_persistent_work_view_t color_xadj;


    size_type c_nnz_size = entryIndicesC_.dimension_0();

    //first we need to transpose the C graph.
    //allocate memory for that.
    row_lno_temp_work_view_t transpose_col_xadj ("transpose_col_xadj", b_col_cnt + 1);
    nnz_lno_temp_work_view_t transpose_col_adj (Kokkos::ViewAllocateWithoutInitializing("tmp_row_view"), c_nnz_size);


    KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space = this->handle->get_handle_exec_space();
    int suggested_vector_size = this->handle->get_suggested_vector_size(rowmapC.dimension_0() - 1, c_nnz_size);
    int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, concurrency,a_row_cnt);

    Kokkos::Impl::Timer timer1;
    KokkosKernels::Experimental::Util::kk_transpose_graph<
        c_row_view_t, c_nnz_view_t,
        row_lno_temp_work_view_t, nnz_lno_temp_work_view_t, row_lno_temp_work_view_t,
        MyExecSpace>
            (a_row_cnt, b_col_cnt,
                rowmapC, entryIndicesC_,
                transpose_col_xadj,transpose_col_adj,
                suggested_vector_size,
                suggested_team_size,
                team_row_chunk_size,
                use_dynamic_schedule);

    MyExecSpace::fence();
    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tTranspose Time:" << timer1.seconds() << std::endl;
    }


    {

      timer1.reset();
      this->handle->create_graph_coloring_handle();

      //for now only sequential one exists.
      handle->get_graph_coloring_handle()->set_algorithm(KokkosKernels::Experimental::Graph::COLORING_SERIAL2);
      //find distance-2 graph coloring
      KokkosKernels::Experimental::Graph::d2_graph_color
        <HandleType,
        c_row_view_t, c_nnz_view_t,
        row_lno_temp_work_view_t, nnz_lno_temp_work_view_t>
        (this->handle, a_row_cnt, b_col_cnt, rowmapC, entryIndicesC_, transpose_col_xadj, transpose_col_adj);

      original_num_colors = handle->get_graph_coloring_handle()->get_num_colors();
      num_multi_color_steps =  original_num_colors;
      num_colors_in_one_step = 1;

      if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tNum colors:" << handle->get_graph_coloring_handle()->get_num_colors() <<  " coloring time:" << timer1.seconds()  << std::endl;
      }

      typename HandleType::GraphColoringHandleType::color_view_t vertex_color_view = handle->get_graph_coloring_handle()->get_vertex_colors();

      vertex_colors_to_store = nnz_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("persistent_color_view"), a_row_cnt);


      if (KOKKOSKERNELS_VERBOSE){
        //create histogram and print it.
        nnz_lno_temp_work_view_t histogram ("histogram", handle->get_graph_coloring_handle()->get_num_colors() + 1);
        MyExecSpace::fence();
        timer1.reset();
        KokkosKernels::Experimental::Util::kk_get_histogram
        <typename HandleType::GraphColoringHandleType::color_view_t, nnz_lno_temp_work_view_t, MyExecSpace>(a_row_cnt, vertex_color_view, histogram);
        std::cout << "\t\tHistogram" << " time:" << timer1.seconds()  << std::endl << "\t\t";
        KokkosKernels::Experimental::Util::kk_print_1Dview(histogram);
      }

      {
        nnz_lno_temp_work_view_t tmp_color_view = vertex_color_view;

        //if the algorithm is spgemm, then we will have multiple colors per iteration.
        if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR){

          tmp_color_view = nnz_lno_temp_work_view_t( Kokkos::ViewAllocateWithoutInitializing("tmp_color_view"), a_row_cnt);
          //upper bound is the output size for dense acumulators.
          num_colors_in_one_step = c_nnz_size / this->b_col_cnt;

          //scale if provided.
          double scale_ = this->handle->get_spgemm_handle()->get_multi_color_scale();

          num_colors_in_one_step = scale_ * num_colors_in_one_step;

          //at the end of this tmp_color_view holds the colors that correspond the step colors.
          // that is if num_multi_colors is 32, first 32 is 1, next 32 is 2 and so on.
          if (num_colors_in_one_step > 1){
            float scale_factor = 1.0 / num_colors_in_one_step;

            //get the sets multicolors. color(i) / num_multi_colors + 1 is the new color.
            KokkosKernels::Experimental::Util::kk_a_times_x_plus_b
                <nnz_lno_temp_work_view_t, typename HandleType::GraphColoringHandleType::color_view_t, float, float, MyExecSpace>(
                    a_row_cnt, tmp_color_view, vertex_color_view, scale_factor, 0);
            num_multi_color_steps = original_num_colors / num_colors_in_one_step;
            if (original_num_colors % num_colors_in_one_step) ++num_multi_color_steps;
          }
          else {
            num_colors_in_one_step = 1;
          }
        }
        else {
          KokkosKernels::Experimental::Util::kk_a_times_x_plus_b
          <nnz_lno_temp_work_view_t, typename HandleType::GraphColoringHandleType::color_view_t, int, int, MyExecSpace>(
              a_row_cnt, tmp_color_view, tmp_color_view, 1, -1);
        }

        if (spgemm_algorithm == KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2){
          num_multi_color_steps = original_num_colors;
          num_colors_in_one_step = c_nnz_size / this->b_col_cnt;
          double scale_ = this->handle->get_spgemm_handle()->get_multi_color_scale();
          num_colors_in_one_step = scale_ * num_colors_in_one_step;
        }

        //with the modular operation, we find the colors within a step.
        // that is if num_multi_colors is 32, then color 32 will be 0, 33 will be 1, 34 will be 2.
        //it will hold their color within their multicolor step.
        KokkosKernels::Experimental::Util::kk_modular_view
          <nnz_lno_persistent_work_view_t, typename HandleType::GraphColoringHandleType::color_view_t, MyExecSpace>
          (a_row_cnt, vertex_colors_to_store, vertex_color_view, num_colors_in_one_step);
        timer1.reset();

        //allocate color xadj and adj arrays.
        color_xadj = nnz_lno_persistent_work_view_t(
            Kokkos::ViewAllocateWithoutInitializing("Reverse xadj"),
            num_multi_color_steps + 1);
        color_adj = nnz_lno_persistent_work_view_t(
            Kokkos::ViewAllocateWithoutInitializing("Reverse xadj"),
            a_row_cnt);

        //create reverse map from colors.
        KokkosKernels::Experimental::Util::kk_create_reverse_map
          <typename HandleType::GraphColoringHandleType::color_view_t,
          nnz_lno_persistent_work_view_t, MyExecSpace>
              (a_row_cnt, num_multi_color_steps, tmp_color_view, color_xadj, color_adj);
        MyExecSpace::fence();

        if (KOKKOSKERNELS_VERBOSE){
          std::cout << "\t\tReverse Map Create Time:" << timer1.seconds() << std::endl;
        }
        h_color_xadj = Kokkos::create_mirror_view (color_xadj);
        Kokkos::deep_copy (h_color_xadj, color_xadj);
        MyExecSpace::fence();
      }
      this->handle->destroy_graph_coloring_handle();
    }

  }

  template <typename c_row_view_t, typename c_nnz_view_t>
  void write_matrix_to_plot(
      nnz_lno_t &num_colors,
      nnz_lno_persistent_work_host_view_t &h_color_xadj,
      nnz_lno_persistent_work_view_t &color_adj,
      c_row_view_t &rowmapC, c_nnz_view_t &entryIndicesC_){
    std::cout << "writing to plot" << std::endl;

    nnz_lno_persistent_work_host_view_t h_color_adj = Kokkos::create_mirror_view (color_adj);
    Kokkos::deep_copy (h_color_adj, color_adj);
    auto h_rowmapC = Kokkos::create_mirror_view (rowmapC);
    Kokkos::deep_copy (h_rowmapC, rowmapC);
    auto h_entryIndicesC = Kokkos::create_mirror_view (entryIndicesC_);
    Kokkos::deep_copy (h_entryIndicesC, entryIndicesC_);

    for (nnz_lno_t i = 0; i < num_colors; ++i){
      nnz_lno_t color_begin = h_color_xadj(i);
      nnz_lno_t color_end = h_color_xadj(i + 1);

      std::string colorind = "";
      std::stringstream ss;
      ss << i;


      ss >> colorind;
      colorind += ".coords";
      std::fstream fs;
      fs.open(colorind.c_str(), std::fstream::out);

      std::cout << "COLOR:" << i << " colorbegin:" << color_begin << " colorend:" << color_end << " size:" << color_end - color_begin << std::endl;
      for (nnz_lno_t j = color_begin; j < color_end; ++j){
        nnz_lno_t row = h_color_adj(j);
        for (size_type k = h_rowmapC(row); k < h_rowmapC(row + 1); ++k){
          nnz_lno_t column = h_entryIndicesC(k);
          //std::cout << row << " " << column << std::endl;
          fs << row << " " << column << std::endl;
        }
      }
      fs.close();
    }



    std::fstream fs;
    fs.open("plot1.gnuplot", std::fstream::out);
    for (nnz_lno_t i = 0; i < num_colors; ++i){
      std::string colorind = "\"";
      std::stringstream ss;
      ss << i;

      ss >> colorind;
      colorind += ".coords\"";
      if (i > 0) fs << "re";
      fs << "plot " << colorind << std::endl;

    }
    fs << "pause -1" << std::endl;
    fs.close();
  }

};
}
}
}
}

#endif
