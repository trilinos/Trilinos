
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


#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <KokkosKernels_GaussSeidel.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
//#include <Kokkos_Sparse_CrsMatrix.hpp>
#include <Kokkos_Sparse.hpp>
#include <Kokkos_Blas1_MV.hpp>

#include "TpetraKernels_ETIHelperMacros.h"
#include <cstdlib>

//const char *input_filename = "sherman1.mtx";
//const char *input_filename = "Si2.mtx";
//const char *input_filename = "wathen_30_30.mtx";
extern char * input_filename;

template <typename crsMat_t, typename device>
int run_gauss_seidel_1(
    crsMat_t input_mat,
    KokkosKernels::Experimental::Graph::GSAlgorithm gs_algorithm,
    typename crsMat_t::values_type::non_const_type x_vector,
    typename crsMat_t::values_type::const_type y_vector,
    bool is_symmetric_graph,
    int apply_type = 0, // 0 for symmetric, 1 for forward, 2 for backward.
    bool skip_symbolic = false,
    bool skip_numeric = false
    ){
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;



  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  kh.create_gs_handle(gs_algorithm);


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();
  const int apply_count = 100;

  if (!skip_symbolic){
    KokkosKernels::Experimental::Graph::gauss_seidel_symbolic
      (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, is_symmetric_graph);
  }

  if (!skip_numeric){
    KokkosKernels::Experimental::Graph::gauss_seidel_numeric
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, is_symmetric_graph);
  }

  switch (apply_type){
  case 0:
    KokkosKernels::Experimental::Graph::symmetric_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, apply_count);
    break;
  case 1:
    KokkosKernels::Experimental::Graph::forward_sweep_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, apply_count);
    break;
  case 2:
    KokkosKernels::Experimental::Graph::backward_sweep_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, apply_count);
    break;
  default:
    KokkosKernels::Experimental::Graph::symmetric_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector,false, true, apply_count);
    break;
  }


  kh.destroy_gs_handle();
  return 0;
}

template<typename scalar_view_t>
scalar_view_t create_x_vector(size_t nv, typename scalar_view_t::value_type max_value = 10.0){
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

template <typename scalar_t, typename lno_t, typename device>
void test_gauss_seidel(KokkosKernels::Experimental::Graph::GSAlgorithm gs_algorithm) {
  ASSERT_TRUE( (input_filename != NULL));
  device::execution_space::initialize();
  srand(245);
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::non_const_type   color_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename lno_view_t::non_const_value_type size_type;

  crsMat_t input_mat = KokkosKernels::Experimental::Util::read_kokkos_crst_matrix<crsMat_t>(input_filename);

  lno_t nv = input_mat.numRows();


  //scalar_view_t solution_x ("sol", nv);
  //Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  //Kokkos::fill_random(solution_x,g,Kokkos::Random_XorShift64_Pool<ExecutionSpace>::generator_type::MAX_URAND);

  const scalar_view_t solution_x = create_x_vector<scalar_view_t>(nv);
  scalar_view_t y_vector = create_y_vector(input_mat, solution_x);
  scalar_view_t x_vector ("x vector", nv);


  const scalar_t alpha = 1.0;
  KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);


  scalar_t dot_product = KokkosBlas::dot( x_vector , x_vector );
  scalar_t initial_norm_res  = sqrt( dot_product );

  Kokkos::deep_copy (x_vector , 0);

  bool is_symmetric_graph = false;
  int apply_type = 0;
  bool skip_symbolic = false;
  bool skip_numeric = false;



  for (int is_symmetric_graph = 0; is_symmetric_graph < 2; ++is_symmetric_graph){
    for (int apply_type = 0; apply_type < 3; ++apply_type){
      for (int skip_symbolic = 0; skip_symbolic < 2; ++skip_symbolic){
        for (int skip_numeric = 0; skip_numeric < 2; ++skip_numeric){

          Kokkos::Impl::Timer timer1;
          int res = run_gauss_seidel_1<crsMat_t, device>(input_mat, gs_algorithm, x_vector, y_vector, is_symmetric_graph, apply_type, skip_symbolic, skip_numeric);
          double gs = timer1.seconds();

          //KokkosKernels::Experimental::Util::print_1Dview(x_vector);
          KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
          //KokkosKernels::Experimental::Util::print_1Dview(x_vector);
          scalar_t result_dot_product = KokkosBlas::dot( x_vector , x_vector );
          scalar_t result_norm_res  = sqrt( result_dot_product );
          //std::cout << "result_norm_res:" << result_norm_res << " initial_norm_res:" << initial_norm_res << std::endl;
          EXPECT_TRUE( (result_norm_res < initial_norm_res));
        }
      }
    }
  }

  device::execution_space::finalize();
}

#define INSTMACRO( SCALAR, LO, DEVICE) \
    test_gauss_seidel<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::GS_DEFAULT); \
    test_gauss_seidel<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::GS_PERMUTED); \
    test_gauss_seidel<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::GS_TEAM);


TEST (GAUSSSEIDEL_TEST, GS) {

  TPETRAKERNELS_ETI_MANGLING_TYPEDEFS()
  TPETRAKERNELS_INSTANTIATE_SLD_NO_ORDINAL_SCALAR(INSTMACRO)
}
 

