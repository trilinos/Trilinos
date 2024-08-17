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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_TestUtils.hpp"
#include "Test_Sparse_Utils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <cstdlib>
#include <iostream>
#include <complex>
#include "KokkosSparse_gauss_seidel.hpp"

using kokkos_complex_double = Kokkos::complex<double>;
using kokkos_complex_float  = Kokkos::complex<float>;

namespace KSExp = KokkosSparse::Experimental;

namespace Test {

enum GSApplyType {
  symmetric,
  forward_sweep,
  backward_sweep,
};

template <typename lno_t, typename scalar_t, typename mag_t>
struct GSTestParams {
  // Intentionally testing block_size that's not a multiple of #rows.
  lno_t block_size = 7;
  lno_t numVecs    = 2;  // how many columns X/Y have
  scalar_t omega   = 0.9;
  mag_t tolerance  = 1e-7;  // relative error for solution x vector

  // Note: GS_DEFAULT is same as GS_TEAM and - for blocks - as GS_PERMUTED
  // Note: GS_TWOSTAGE and GS_CLUSTER are not supported for blocks
  std::vector<KokkosSparse::GSAlgorithm> gs_algorithms = {KokkosSparse::GS_DEFAULT};
  std::vector<size_t> shmem_sizes                      = {
      32128,
      2008  // make the shmem small on gpus so that it will test 2 level
            // algorithm.
  };
  std::vector<GSApplyType> apply_types = {symmetric, forward_sweep, backward_sweep};

  GSTestParams() = default;
};

template <typename mtx_t, typename vector_t, typename const_vector_t>
int run_block_gauss_seidel_1(
    mtx_t input_mat, int block_size, KokkosSparse::GSAlgorithm gs_algorithm, vector_t x_vector, const_vector_t y_vector,
    bool is_symmetric_graph, GSApplyType apply_type = Test::symmetric, bool skip_symbolic = false,
    bool skip_numeric = false, size_t shmem_size = 32128,
    typename mtx_t::value_type omega = Kokkos::ArithTraits<typename mtx_t::value_type>::one()) {
  typedef typename mtx_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename mtx_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  constexpr auto format = KokkosSparse::Impl::MatrixTraits<mtx_t>::format;

  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename mtx_t::execution_space,
                                                       typename mtx_t::memory_space, typename mtx_t::memory_space>;
  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_shmem_size(shmem_size);
  kh.set_dynamic_scheduling(true);
  kh.create_gs_handle(gs_algorithm);

  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();
  const int apply_count   = 100;

  if (!skip_symbolic) {
    KSExp::block_gauss_seidel_symbolic(&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map,
                                       input_mat.graph.entries, is_symmetric_graph);
  }

  if (!skip_numeric) {
    KSExp::block_gauss_seidel_numeric<format>(&kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map,
                                              input_mat.graph.entries, input_mat.values, is_symmetric_graph);
  }

  switch (apply_type) {
    case Test::forward_sweep:
      KSExp::forward_sweep_block_gauss_seidel_apply<format>(
          &kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values,
          x_vector, y_vector, false, true, omega, apply_count);
      break;
    case Test::backward_sweep:
      KSExp::backward_sweep_block_gauss_seidel_apply<format>(
          &kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values,
          x_vector, y_vector, false, true, omega, apply_count);
      break;
    case Test::symmetric:
    default:
      KSExp::symmetric_block_gauss_seidel_apply<format>(
          &kh, num_rows_1, num_cols_1, block_size, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values,
          x_vector, y_vector, false, true, omega, apply_count);
      break;
  }

  kh.destroy_gs_handle();
  return 0;
}

}  // namespace Test

template <KokkosSparse::SparseMatrixFormat mtx_format, typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_block_gauss_seidel_rank1(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using namespace Test;
  srand(245);
  using crsMat_t        = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using MatrixConverter = KokkosSparse::Impl::MatrixConverter<mtx_format>;
  typedef typename device::execution_space exec_space;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;

  lno_t numCols = numRows;

  const GSTestParams<lno_t, scalar_t, mag_t> params;
  lno_t block_size = params.block_size;

  crsMat_t crsmat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);

  lno_view_t pf_rm;
  lno_nnz_view_t pf_e;
  scalar_view_t pf_v;
  size_t out_r, out_c;

  // this makes consecutive 5 rows to have same columns.
  // it will add scalar 0's for those entries that does not exists.
  // the result is still a point crs matrix.
  KokkosSparse::Impl::kk_create_bsr_formated_point_crsmatrix(block_size, crsmat.numRows(), crsmat.numCols(),
                                                             crsmat.graph.row_map, crsmat.graph.entries, crsmat.values,
                                                             out_r, out_c, pf_rm, pf_e, pf_v);
  graph_t static_graph2(pf_e, pf_rm);
  crsMat_t crsmat2("CrsMatrix2", out_c, pf_v, static_graph2);

  // this converts the previous generated matrix to block matrix.
  auto input_mat = MatrixConverter::from_bsr_formated_point_crsmatrix(crsmat2, block_size);

  lno_t nv = ((crsmat2.numRows() + block_size - 1) / block_size) * block_size;

  const scalar_view_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X"), nv);
  // create_random_x_vector operates on host mirror, then copies to device. But
  // create_y does everything on device.
  create_random_x_vector(solution_x);
  exec_space().fence();
  scalar_view_t y_vector = create_random_y_vector(crsmat2, solution_x);
  mag_t initial_norm_res = KokkosBlas::nrm2(solution_x);

  for (const auto gs_algorithm : params.gs_algorithms) {
    scalar_view_t x_vector("x vector", nv);
    const scalar_t alpha = 1.0;

    bool is_symmetric_graph = true;

    for (const auto shmem_size : params.shmem_sizes) {
      for (const auto apply_type : params.apply_types) {
        for (const auto skip_symbolic : {false, true}) {
          for (const auto skip_numeric : {false, true}) {
            Kokkos::Timer timer1;
            // int res =
            run_block_gauss_seidel_1(input_mat, block_size, gs_algorithm, x_vector, y_vector, is_symmetric_graph,
                                     apply_type, skip_symbolic, skip_numeric, shmem_size, params.omega);
            // double gs = timer1.seconds();
            // KokkosKernels::Impl::print_1Dview(x_vector);
            KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
            mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
            EXPECT_LT(result_norm_res, params.tolerance * initial_norm_res);
          }
        }
      }
    }
  }
  // device::execution_space::finalize();
}

template <KokkosSparse::SparseMatrixFormat mtx_format, typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_block_gauss_seidel_rank2(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using namespace Test;
  srand(245);
  using crsMat_t        = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using MatrixConverter = KokkosSparse::Impl::MatrixConverter<mtx_format>;

  typedef typename device::execution_space exec_space;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;
  typedef Kokkos::View<scalar_t**, default_layout, device> scalar_view2d_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;

  lno_t numCols = numRows;

  const GSTestParams<lno_t, scalar_t, mag_t> params;
  lno_t block_size = params.block_size;

  crsMat_t crsmat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);

  lno_view_t pf_rm;
  lno_nnz_view_t pf_e;
  scalar_view_t pf_v;
  size_t out_r, out_c;

  // this makes consecutive 5 rows to have same columns.
  // it will add scalar 0's for those entries that does not exists.
  // the result is still a point crs matrix.
  KokkosSparse::Impl::kk_create_bsr_formated_point_crsmatrix(block_size, crsmat.numRows(), crsmat.numCols(),
                                                             crsmat.graph.row_map, crsmat.graph.entries, crsmat.values,
                                                             out_r, out_c, pf_rm, pf_e, pf_v);
  graph_t static_graph2(pf_e, pf_rm);
  crsMat_t crsmat2("CrsMatrix2", out_c, pf_v, static_graph2);

  auto input_mat = MatrixConverter::from_bsr_formated_point_crsmatrix(crsmat2, block_size);

  lno_t nv = ((crsmat2.numRows() + block_size - 1) / block_size) * block_size;

  const lno_t numVecs = params.numVecs;

  scalar_view2d_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X"), nv, params.numVecs);
  create_random_x_vector(solution_x);
  scalar_view2d_t y_vector = create_random_y_vector_mv(crsmat2, solution_x);
  exec_space().fence();
  auto solution_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solution_x);
  // Need to fence before reading from solution_host
  std::vector<mag_t> initial_norms(numVecs);
  for (lno_t i = 0; i < numVecs; i++) {
    scalar_t sum = 0;
    for (lno_t j = 0; j < nv; j++) {
      sum += solution_host(j, i) * solution_host(j, i);
    }
    initial_norms[i] = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(sum));
  }

  for (const auto gs_algorithm : params.gs_algorithms) {
    scalar_view2d_t x_vector("x vector", nv, numVecs);
    auto x_host = Kokkos::create_mirror_view(x_vector);

    bool is_symmetric_graph = true;

    scalar_view_t res_norms("Residuals", numVecs);
    auto h_res_norms = Kokkos::create_mirror_view(res_norms);

    for (const auto shmem_size : params.shmem_sizes) {
      for (const auto apply_type : params.apply_types) {
        for (const auto skip_symbolic : {false, true}) {
          for (const auto skip_numeric : {false, true}) {
            Kokkos::Timer timer1;
            // int res =
            run_block_gauss_seidel_1(input_mat, block_size, gs_algorithm, x_vector, y_vector, is_symmetric_graph,
                                     apply_type, skip_symbolic, skip_numeric, shmem_size, params.omega);
            // double gs = timer1.seconds();
            // KokkosKernels::Impl::print_1Dview(x_vector);
            Kokkos::deep_copy(x_host, x_vector);
            exec_space().fence();
            for (lno_t c = 0; c < numVecs; c++) {
              scalar_t sum = 0;
              for (lno_t r = 0; r < nv; r++) {
                scalar_t diff = x_host(r, c) - solution_host(r, c);
                sum += diff * diff;
              }
              mag_t result_res = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(sum));
              EXPECT_LT(result_res, params.tolerance * initial_norms[c]);
            }
          }
        }
      }
    }
  }
  // device::execution_space::finalize();
}

template <KokkosSparse::SparseMatrixFormat mtx_format, typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_block_gauss_seidel_empty() {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_type;
  typedef typename graph_t::entries_type::non_const_type entries_type;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                       typename device::memory_space, typename device::memory_space>;
  // The rowmap of a zero-row matrix can be length 0 or 1, so Gauss-Seidel
  // should work with both (the setup and apply are essentially no-ops but they
  // shouldn't crash or throw exceptions) For this test, create size-0 and
  // size-1 rowmaps separately. Check also 5x5 matrix with empty rows (0-nnz),
  // which can trigger different bugs.
  for (const int rowmapLen : {0, 1, 5}) {
    KernelHandle kh;
    kh.create_gs_handle(KokkosSparse::GS_DEFAULT);
    const auto num_rows    = KOKKOSKERNELS_MACRO_MAX(0, rowmapLen - 1);
    const lno_t block_size = 1;  // irrelevant (no values here)
    // initialized to 0
    row_map_type rowmap("Rowmap", rowmapLen);
    entries_type entries("Entries", 0);
    scalar_view_t values("Values", 0);
    // also, make sure graph symmetrization doesn't crash on zero rows
    KSExp::block_gauss_seidel_symbolic(&kh, num_rows, num_rows, block_size, rowmap, entries, false);
    KSExp::block_gauss_seidel_numeric<mtx_format>(&kh, num_rows, num_rows, block_size, rowmap, entries, values, false);
    scalar_view_t x("X", num_rows);
    scalar_view_t y("Y", num_rows);
    scalar_t omega(0.9);
    KSExp::symmetric_block_gauss_seidel_apply<mtx_format>(&kh, num_rows, num_rows, block_size, rowmap, entries, values,
                                                          x, y, false, true, omega, 3);
    kh.destroy_gs_handle();
  }
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                          \
  TEST_F(TestCategory, sparse_bsr_gauss_seidel_rank1_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {           \
    test_block_gauss_seidel_rank1<KokkosSparse::SparseMatrixFormat::BSR, SCALAR, ORDINAL, OFFSET, DEVICE>(   \
        500, 500 * 10, 70, 3);                                                                               \
  }                                                                                                          \
  TEST_F(TestCategory, sparse_bsr_gauss_seidel_rank2_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {           \
    test_block_gauss_seidel_rank2<KokkosSparse::SparseMatrixFormat::BSR, SCALAR, ORDINAL, OFFSET, DEVICE>(   \
        500, 500 * 10, 70, 3);                                                                               \
  }                                                                                                          \
  TEST_F(TestCategory, sparse_bsr_gauss_seidel_empty_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {           \
    test_block_gauss_seidel_empty<KokkosSparse::SparseMatrixFormat::BSR, SCALAR, ORDINAL, OFFSET, DEVICE>(); \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
