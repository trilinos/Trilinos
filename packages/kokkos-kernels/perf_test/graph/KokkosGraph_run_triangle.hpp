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


#include "KokkosGraph_Triangle.hpp"
#include "KokkosKernels_TestParameters.hpp"

#define TRANPOSEFIRST false
#define TRANPOSESECOND false

namespace KokkosKernels{

namespace Experiment{
template <typename crsGraph_t, typename device>
bool is_same_graph(crsGraph_t output_mat1, crsGraph_t output_mat2){

  //typedef typename crsGraph_t::StaticCrsGraphType crsGraph_t;
  typedef typename crsGraph_t::row_map_type::non_const_type lno_view_t;
  typedef typename crsGraph_t::entries_type::non_const_type   lno_nnz_view_t;
  //typedef typename crsGraph_t::values_type::non_const_type scalar_view_t;

  size_t nrows1 = output_mat1.row_map.extent(0);
  size_t nentries1 = output_mat1.entries.extent(0) ;

  size_t nrows2 = output_mat2.row_map.extent(0);
  size_t nentries2 = output_mat2.entries.extent(0) ;
  //size_t nvals2 = output_mat2.values.extent(0);


  lno_nnz_view_t h_ent1 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries1);
  lno_nnz_view_t h_vals1 (Kokkos::ViewAllocateWithoutInitializing("v1"), nentries1);


  KokkosKernels::Impl::kk_sort_graph<typename crsGraph_t::row_map_type,
    typename crsGraph_t::entries_type,
    lno_nnz_view_t,
    lno_nnz_view_t,
    lno_nnz_view_t,
    typename device::execution_space
    >(
    output_mat1.row_map, output_mat1.entries,h_vals1,
    h_ent1, h_vals1
  );

  lno_nnz_view_t h_ent2 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries2);
  lno_nnz_view_t h_vals2 (Kokkos::ViewAllocateWithoutInitializing("v1"), nentries2);

  if (nrows1 != nrows2) return false;
  if (nentries1 != nentries2) return false;

  KokkosKernels::Impl::kk_sort_graph
      <typename crsGraph_t::row_map_type,
      typename crsGraph_t::entries_type,
      lno_nnz_view_t,
      lno_nnz_view_t,
      lno_nnz_view_t,
      typename device::execution_space
      >(
      output_mat2.row_map, output_mat2.entries, h_vals2,
      h_ent2, h_vals2
    );

  bool is_identical = true;
  is_identical = KokkosKernels::Impl::kk_is_identical_view
      <typename crsGraph_t::row_map_type, typename crsGraph_t::row_map_type, typename lno_view_t::value_type,
      typename device::execution_space>(output_mat1.row_map, output_mat2.row_map, 0);
  if (!is_identical) return false;

  is_identical = KokkosKernels::Impl::kk_is_identical_view
      <lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(h_ent1, h_ent2, 0 );
  if (!is_identical) return false;

  if (!is_identical) {
    std::cout << "Incorret values" << std::endl;
  }
  return true;
}

template<size_t BufSize, typename SpaceType = Kokkos::DefaultExecutionSpace>
struct Flush {
  typedef double value_type;

  // flush a large host buffer
  Kokkos::View<value_type*,SpaceType> _buf;
  Flush(int flush_option) : _buf("Flush::buf", BufSize) {
    Kokkos::deep_copy(_buf, 1); 
    Kokkos::fence();
    if (flush_option == 2){
    for (size_t i = 0; i < BufSize; ++i){ 
      _buf(i) = rand();
    }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &update) {
    update = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type &update,
            const volatile value_type &input) {
    update += input;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &update) const {
    update += _buf[i];
  }

  void run() {
    double sum = 0;
    Kokkos::parallel_reduce("KokkosGraph::PerfTest::Flush", Kokkos::RangePolicy<SpaceType>(0,BufSize/sizeof(double)), *this, sum);
    SpaceType().fence();
    std::cout << "Flush sum:" << sum << std::endl;
    FILE *fp = fopen("/dev/null", "w");
    fprintf(fp, "%f\n", sum);
    fclose(fp);
   
/*
#pragma omp parallel 
    {
    const size_t cache_line = 64;
    const char *cp = (const char *) _buf.data();
    size_t i = 0;


    for (i = 0; i < BufSize; i += cache_line) {
            asm volatile("clflush (%0)\n\t"
                         :
                         : "r"(&cp[i])
                         : "memory");
    }

    asm volatile("sfence\n\t"
                 :
                 :
                 : "memory");
  }
*/
  }

};

template <typename ExecSpace, typename crsGraph_t, typename crsGraph_t2 , typename crsGraph_t3 , typename TempMemSpace , typename PersistentMemSpace >
void run_experiment(
    crsGraph_t crsGraph, Parameters params){
  //using namespace KokkosSparse;
  using namespace KokkosSparse;
  using namespace KokkosGraph::Experimental;
  //using namespace KokkosSparse::Experimental;

  int algorithm = params.algorithm;
  int repeat = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize = params.shmemsize;
  int team_size = params.team_size;
  int use_dynamic_scheduling = params.use_dynamic_scheduling;
  int verbose = params.verbose;

  int accumulator = params.accumulator;
  //char spgemm_step = params.spgemm_step;
  int vector_size = params.vector_size;

  //spgemm_step++;

  typedef typename crsGraph_t3::row_map_type::non_const_type lno_view_t;
  typedef typename crsGraph_t3::entries_type::non_const_type lno_nnz_view_t;

  Kokkos::View <size_t *,ExecSpace> row_mapC;
  lno_nnz_view_t entriesC;
  lno_nnz_view_t valuesC;

  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename lno_view_t::value_type size_type;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <size_type,lno_t, lno_t,
      ExecSpace, TempMemSpace,PersistentMemSpace > KernelHandle;


  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);


  if (use_dynamic_scheduling){
    kh.set_dynamic_scheduling(true);
  }
  if (verbose){
    kh.set_verbose(true);
  }
  const lno_t m = crsGraph.numRows();;



  for (int i = 0; i < repeat; ++i){
    size_type rowmap_size = crsGraph.entries.extent(0) ;
    switch (algorithm){
    case 16:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_AI);
      rowmap_size = m ;
      break;
    case 17:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA);
      std::cout << "IA" << std::endl;
      break;
    case 18:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA_UNION);
      break;
    case 19:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_LL);
      rowmap_size = m ;
      break;
    case 20:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_LU);
      rowmap_size = m ;
      break;
    default:
      kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA);
      break;
    }

    kh.get_spgemm_handle()->set_compression_steps(!params.compression2step);

    kh.get_spgemm_handle()->set_sort_lower_triangular(params.right_sort);
    kh.get_spgemm_handle()->set_create_lower_triangular(params.right_lower_triangle);
    kh.get_spgemm_handle()->set_compression(params.apply_compression);
    kh.get_spgemm_handle()->set_sort_option(params.sort_option);
    kh.get_spgemm_handle()->set_min_hash_size_scale(params.minhashscale);

    switch (accumulator){
    case 0:
    default:
      kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_DEFAULT);
      break;
    case 1:
      kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_DENSE);
      break;
    case 2:
      kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_SPARSE);
      break;
    }


    constexpr size_t LLC_CAPACITY = 256*4*1024*1024;
    if (params.cache_flush)
    {
    std::cout << "Flushing cache with option:" << params.cache_flush << std::endl;
    Flush<LLC_CAPACITY, ExecSpace> flush(params.cache_flush);
        flush.run();
    }
    if (i == 0){
      kh.get_spgemm_handle()->set_read_write_cost_calc(params.calculate_read_write_cost);
    }



    Kokkos::Impl::Timer timer1;

    row_mapC = Kokkos::View <size_t *,ExecSpace>
              ("non_const_lnow_row",
                  rowmap_size);
    entriesC = lno_nnz_view_t ("");
    valuesC  = lno_nnz_view_t ("");

    double symbolic_time = 0;
    if (params.triangle_options == 0 ){
      if (params.apply_compression){
        triangle_generic (
            &kh,
            m,
            crsGraph.row_map,
            crsGraph.entries,
            KOKKOS_LAMBDA(const lno_t& row, const lno_t &col_set_index, const lno_t &col_set,  const lno_t &thread_id) {

          //row_mapC(row) += KokkosKernels::Impl::set_bit_count<lno_t, ExecSpace>(col_set);
          row_mapC(row) += KokkosKernels::Impl::pop_count(col_set);
        }
        );
      }
      else {
        triangle_generic (
            &kh,
            m,
            crsGraph.row_map,
            crsGraph.entries,
            KOKKOS_LAMBDA(const lno_t& row, const lno_t &col_set_index, const lno_t &col_set,  const lno_t &thread_id) {

          row_mapC(row) += 1;
          //row_mapC(row) += KokkosKernels::Impl::set_bit_count<lno_t, ExecSpace>(col_set);
          //row_mapC(row) += KokkosKernels::Impl::pop_count(col_set);

        }
        );
      }

      size_t num_triangles = 0;
      KokkosKernels::Impl::kk_reduce_view< Kokkos::View <size_t *,ExecSpace>, ExecSpace>(rowmap_size, row_mapC, num_triangles);
      ExecSpace().fence();

      symbolic_time = timer1.seconds();
      std::cout << "num_triangles:" << num_triangles << std::endl;
    }
    kh.destroy_spgemm_handle();
    std::cout  << "mm_time:" << symbolic_time << std::endl;
    //only do this once
    //kh.get_spgemm_handle()->set_read_write_cost_calc(false);
  }
}


}
}
