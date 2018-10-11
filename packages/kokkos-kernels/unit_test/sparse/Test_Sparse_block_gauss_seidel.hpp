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


#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
//#include <Kokkos_Sparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <cstdlib>
#include <iostream>
#include <complex>
#include "KokkosSparse_gauss_seidel.hpp"

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
namespace Test {

template <typename crsMat_t, typename device>
int run_block_gauss_seidel_1(
    crsMat_t input_mat, int block_size,
    KokkosSparse::GSAlgorithm gs_algorithm,
    typename crsMat_t::values_type::non_const_type x_vector,
    typename crsMat_t::values_type::const_type y_vector,
    bool is_symmetric_graph,
    int apply_type = 0, // 0 for symmetric, 1 for forward, 2 for backward.
    bool skip_symbolic = false,
    bool skip_numeric = false,
    size_t shmem_size = 32128,
    typename crsMat_t::value_type omega = Kokkos::Details::ArithTraits<typename crsMat_t::value_type>::one()
    ){
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;



  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernelsHandle
      <size_type,lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;



  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_shmem_size(shmem_size);
  kh.set_dynamic_scheduling(true);
  kh.create_gs_handle(gs_algorithm);


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();
  //std::cout << "num_rows_1:" << num_rows_1 << " num_cols_1:" << num_cols_1 << std::endl;
  const int apply_count = 100;

  if (!skip_symbolic){
	  block_gauss_seidel_symbolic
      (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, is_symmetric_graph);
  }

  if (!skip_numeric){
	  block_gauss_seidel_numeric
    (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, is_symmetric_graph);
  }

  switch (apply_type){
  case 0:
    symmetric_block_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, omega, apply_count);
    break;
  case 1:
    forward_sweep_block_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, omega, apply_count);
    break;
  case 2:
    backward_sweep_block_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, omega, apply_count);
    break;
  default:
    symmetric_block_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, omega, apply_count);
    break;
  }


  kh.destroy_gs_handle();
  return 0;
}

template<typename scalar_view_t>
scalar_view_t create_x_vector(size_t nv, double max_value = 10.0){
  scalar_view_t kok_x ("X", nv);


  typename scalar_view_t::HostMirror h_x =  Kokkos::create_mirror_view (kok_x);


  for (size_t i = 0; i < nv; ++i){
    typename scalar_view_t::value_type r =
        static_cast <typename scalar_view_t::value_type> (rand()) /
        static_cast <typename scalar_view_t::value_type> (RAND_MAX / max_value);
    h_x(i) = r;
    //h_x(i) = 1;
  }
  Kokkos::deep_copy (kok_x, h_x);


  return kok_x;
}
template <typename crsMat_t, typename vector_t>
vector_t create_y_vector(crsMat_t crsMat, vector_t x_vector){
  vector_t y_vector ("Y VECTOR", crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 1, y_vector);
  return y_vector;
}

}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_block_gauss_seidel(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {

  using namespace Test;
  srand(245);
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;


  lno_t numCols = numRows;
  crsMat_t crsmat = KokkosKernels::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);


  lno_view_t pf_rm;
  lno_nnz_view_t pf_e;
  scalar_view_t pf_v;
  size_t out_r, out_c;
  int block_size = 7;

  //this makes consecutive 5 rows to have same columns.
  //it will add scalar 0's for those entries that does not exists.
  //the result is still a point crs matrix.
  KokkosKernels::Impl::kk_create_blockcrs_formated_point_crsmatrix(
		  block_size , crsmat.numRows(), crsmat.numCols(),
		  crsmat.graph.row_map, crsmat.graph.entries, crsmat.values,
		  out_r, out_c,
		  pf_rm, pf_e, pf_v);
  graph_t static_graph2 (pf_e, pf_rm);
  crsMat_t crsmat2("CrsMatrix2", out_c, pf_v, static_graph2);

  lno_view_t bf_rm;
  lno_nnz_view_t bf_e;
  scalar_view_t bf_v;
  size_t but_r, but_c;

  //this converts the previous generated matrix to block crs matrix.
  KokkosKernels::Impl::kk_create_blockcrs_from_blockcrs_formatted_point_crs(
		  block_size , out_r, out_c,
		  pf_rm, pf_e, pf_v,
		  but_r, but_c,
		  bf_rm, bf_e, bf_v);
  graph_t static_graph (bf_e, bf_rm);
  crsMat_t input_mat("CrsMatrix", but_c, bf_v, static_graph);

  lno_t nv = ((crsmat2.numRows() / block_size)+1) * block_size;


  //scalar_view_t solution_x ("sol", nv);
  //Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  //Kokkos::fill_random(solution_x,g,Kokkos::Random_XorShift64_Pool<ExecutionSpace>::generator_type::MAX_URAND);
  //values_view_t kok_x_original = create_x_vector<values_view_t>(((nv /5) + 1) * 5, MAXVAL);

  const scalar_view_t solution_x = create_x_vector<scalar_view_t>(nv);
  scalar_view_t y_vector = create_y_vector(crsmat2, solution_x);
#ifdef gauss_seidel_testmore
  GSAlgorithm gs_algorithms[] ={GS_DEFAULT, GS_TEAM, GS_PERMUTED};
  int apply_count = 3;
  for (int ii = 0; ii < 3; ++ii){
#else
  int apply_count = 1;
  GSAlgorithm gs_algorithms[] ={GS_DEFAULT};
  for (int ii = 0; ii < 1; ++ii){
#endif
    GSAlgorithm gs_algorithm = gs_algorithms[ii];
    scalar_view_t x_vector ("x vector", nv);
    const scalar_t alpha = 1.0;
    KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
    scalar_t dot_product = KokkosBlas::dot( x_vector , x_vector );
    typedef typename Kokkos::Details::ArithTraits<scalar_t>::mag_type mag_t;
    mag_t initial_norm_res = Kokkos::Details::ArithTraits<scalar_t>::abs (dot_product);
    initial_norm_res  = Kokkos::Details::ArithTraits<mag_t>::sqrt( initial_norm_res );
    Kokkos::deep_copy (x_vector , 0);

    //bool is_symmetric_graph = false;
    //int apply_type = 0;
    //bool skip_symbolic = false;
    //bool skip_numeric = false;
    scalar_t omega = 0.9;


    bool is_symmetric_graph = true;
    size_t shmem_size = 32128;
    
    for(int i = 0; i < 2; ++i)
    {
      if (i == 1) shmem_size = 2008; //make the shmem small on gpus so that it will test 2 level algorithm.
      for (int apply_type = 0; apply_type < apply_count; ++apply_type){
        for (int skip_symbolic = 0; skip_symbolic < 2; ++skip_symbolic){
          for (int skip_numeric = 0; skip_numeric < 2; ++skip_numeric){

            Kokkos::Impl::Timer timer1;
            //int res =
            run_block_gauss_seidel_1<crsMat_t, device>(input_mat, block_size, gs_algorithm, x_vector, y_vector, is_symmetric_graph, apply_type, skip_symbolic, skip_numeric, shmem_size, omega);
            //double gs = timer1.seconds();

            //KokkosKernels::Impl::print_1Dview(x_vector);
            KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
            //KokkosKernels::Impl::print_1Dview(x_vector);
            scalar_t result_dot_product = KokkosBlas::dot( x_vector , x_vector );
            mag_t result_norm_res  = Kokkos::Details::ArithTraits<scalar_t>::abs( result_dot_product );
            result_norm_res = Kokkos::Details::ArithTraits<mag_t>::sqrt(result_norm_res);
            //std::cout << "result_norm_res:" << result_norm_res << " initial_norm_res:" << initial_norm_res << std::endl;
            EXPECT_TRUE( (result_norm_res < initial_norm_res));
          }
        }
      }
    }
  }
  //device::execution_space::finalize();
}



#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ## block_gauss_seidel ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
	test_block_gauss_seidel<SCALAR,ORDINAL,OFFSET,DEVICE>(10000, 10000 * 30, 200, 10); \
}


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif




