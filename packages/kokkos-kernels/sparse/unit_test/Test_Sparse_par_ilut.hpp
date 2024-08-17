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

#include <string>
#include <stdexcept>

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_par_ilut.hpp"
#include "KokkosSparse_gmres.hpp"
#include "KokkosSparse_LUPrec.hpp"
#include "KokkosSparse_SortCrs.hpp"

#include "Test_vector_fixtures.hpp"

#include <gtest/gtest.h>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

namespace Test {

namespace ParIlut {

template <class T>
struct TolMeta {
  static constexpr T value = 1e-8;
};

template <>
struct TolMeta<float> {
  static constexpr float value = 1e-5;  // Lower tolerance for floats
};

}  // namespace ParIlut

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_par_ilut() {
  using RowMapType  = Kokkos::View<size_type*, device>;
  using EntriesType = Kokkos::View<lno_t*, device>;
  using ValuesType  = Kokkos::View<scalar_t*, device>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                       typename device::memory_space, typename device::memory_space>;

  // Simple test fixture A
  std::vector<std::vector<scalar_t>> A = {
      {1., 6., 4., 7.}, {2., -5., 0., 8.}, {0.5, -3., 6., 0.}, {0.2, -0.5, -9., 0.}};

  // Allocate device CRS views for A
  RowMapType row_map("row_map", 0);
  EntriesType entries("entries", 0);
  ValuesType values("values", 0);

  compress_matrix(row_map, entries, values, A);

  const size_type nrows = A.size();

  // Make kernel handle
  KernelHandle kh;

  kh.create_par_ilut_handle();

  auto par_ilut_handle = kh.get_par_ilut_handle();
  par_ilut_handle->set_async_update(false);

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", nrows + 1);
  RowMapType U_row_map("U_row_map", nrows + 1);

  // Initial L/U approximations for A
  par_ilut_symbolic(&kh, row_map, entries, L_row_map, U_row_map);

  const size_type nnzL = par_ilut_handle->get_nnzL();
  const size_type nnzU = par_ilut_handle->get_nnzU();

  ASSERT_EQ(nnzL, 10);
  ASSERT_EQ(nnzU, 8);

  EntriesType L_entries("L_entries", nnzL);
  ValuesType L_values("L_values", nnzL);
  EntriesType U_entries("U_entries", nnzU);
  ValuesType U_values("U_values", nnzU);

  par_ilut_numeric(&kh, row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);

  // Use this to check LU
  // std::vector<std::vector<scalar_t> > expected_LU = {
  //   {1.0, 6.0, 4.0, 7.0},
  //   {2.0, 7.0, 8.0, 22.0},
  //   {0.5, 18.0, 8.0, -20.5},
  //   {0.2, 3.7, -53.2, -1.60}
  // };

  // check_matrix("LU numeric", L_row_map, L_entries, L_values, expected_LU);

  // Use these fixtures to test add_candidates
  // std::vector<std::vector<scalar_t> > expected_L_candidates = {
  //   {1., 0., 0., 0.},
  //   {2., 1., 0., 0.},
  //   {0.50, -3., 1., 0.},
  //   {0.20, -0.50, -9., 1.}
  // };

  // check_matrix("L numeric", L_row_map, L_entries, L_values,
  // expected_L_candidates);

  // std::vector<std::vector<scalar_t> > expected_U_candidates = {
  //   {1., 6., 4., 7.},
  //   {0., -5., -8., 8.},
  //   {0., 0., 6., 20.50},
  //   {0., 0., 0., 1.}
  // };

  // check_matrix("U numeric", U_row_map, U_entries, U_values,
  // expected_U_candidates);

  // Use these fixtures to test compute_l_u_factors
  // std::vector<std::vector<scalar_t> > expected_L_candidates = {
  //   {1., 0., 0., 0.},
  //   {2., 1., 0., 0.},
  //   {0.50, 0.35, 1., 0.},
  //   {0.20, 0.10, -1.32, 1.}
  // };

  // check_matrix("L numeric", L_row_map, L_entries, L_values,
  // expected_L_candidates);

  // std::vector<std::vector<scalar_t> > expected_U_candidates = {
  //   {1., 6., 4., 7.},
  //   {0., -17., -8., -6.},
  //   {0., 0., 6.82, -1.38},
  //   {0., 0., 0., -2.62}
  // };

  // check_matrix("U numeric", U_row_map, U_entries, U_values,
  // expected_U_candidates);

  // Use these fixtures to test full numeric
  std::vector<std::vector<scalar_t>> expected_L_candidates = {
      {1., 0., 0., 0.}, {2., 1., 0., 0.}, {0.50, 0.35, 1., 0.}, {0., 0., -1.32, 1.}};

  check_matrix("L numeric", L_row_map, L_entries, L_values, expected_L_candidates);

  std::vector<std::vector<scalar_t>> expected_U_candidates = {
      {1., 6., 4., 7.},
      {0., -17., -8., -6.},
      {0., 0., 6.82, 0.},
      {0., 0., 0., 0.}  // [3] = 0 for full alg, -2.62 for post-threshold only
  };

  check_matrix("U numeric", U_row_map, U_entries, U_values, expected_U_candidates);

  kh.destroy_par_ilut_handle();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_par_ilut_precond() {
  // Test using par_ilut as a preconditioner
  // Does (LU)^inv Ax = (LU)^inv b converge faster than solving Ax=b?
  using exe_space      = typename device::execution_space;
  using mem_space      = typename device::memory_space;
  using RowMapType     = Kokkos::View<size_type*, device>;
  using EntriesType    = Kokkos::View<lno_t*, device>;
  using ValuesType     = Kokkos::View<scalar_t*, device>;
  using sp_matrix_type = KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exe_space, mem_space, mem_space>;
  using float_t = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  // Create a diagonally dominant sparse matrix to test:
  //  par_ilut settings max_iters, res_delta_stop, fill_in_limit, and
  //  async_update are all left as defaults
  constexpr auto n             = 5000;
  constexpr auto m             = 15;
  constexpr auto tol           = ParIlut::TolMeta<float_t>::value;
  constexpr auto numRows       = n;
  constexpr auto numCols       = n;
  constexpr auto diagDominance = 1;
  constexpr bool verbose       = false;

  size_type nnz = 10 * numRows;
  auto A        = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<sp_matrix_type>(
      numRows, numCols, nnz, 0, lno_t(0.01 * numRows), diagDominance);

  KokkosSparse::sort_crs_matrix(A);

  // Make kernel handles
  KernelHandle kh;
  kh.create_gmres_handle(m, tol);
  auto gmres_handle = kh.get_gmres_handle();
  gmres_handle->set_verbose(verbose);
  using GMRESHandle    = typename std::remove_reference<decltype(*gmres_handle)>::type;
  using ViewVectorType = typename GMRESHandle::nnz_value_view_t;

  kh.create_par_ilut_handle();
  auto par_ilut_handle = kh.get_par_ilut_handle();
  par_ilut_handle->set_verbose(verbose);

  // Pull out views from CRS
  auto row_map = A.graph.row_map;
  auto entries = A.graph.entries;
  auto values  = A.values;

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", numRows + 1);
  RowMapType U_row_map("U_row_map", numRows + 1);

  // Initial L/U approximations for A
  par_ilut_symbolic(&kh, row_map, entries, L_row_map, U_row_map);

  const size_type nnzL = par_ilut_handle->get_nnzL();
  const size_type nnzU = par_ilut_handle->get_nnzU();

  EntriesType L_entries("L_entries", nnzL);
  ValuesType L_values("L_values", nnzL);
  EntriesType U_entries("U_entries", nnzU);
  ValuesType U_values("U_values", nnzU);

  par_ilut_numeric(&kh, row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);

  // Create CRSs
  sp_matrix_type L("L", numRows, numCols, L_values.extent(0), L_values, L_row_map, L_entries),
      U("U", numRows, numCols, U_values.extent(0), U_values, U_row_map, U_entries);

  // Set initial vectors:
  ViewVectorType X("X", n);    // Solution and initial guess
  ViewVectorType Wj("Wj", n);  // For checking residuals at end.
  ViewVectorType B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"),
                   n);  // right-hand side vec
  // Make rhs ones so that results are repeatable:
  Kokkos::deep_copy(B, 1.0);

  int num_iters_plain(0), num_iters_precond(0);

  // Solve Ax = b
  {
    gmres(&kh, A, B, X);

    // Double check residuals at end of solve:
    float_t nrmB = KokkosBlas::nrm2(B);
    KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
    KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
    float_t endRes = KokkosBlas::nrm2(B) / nrmB;

    const auto conv_flag = gmres_handle->get_conv_flag_val();
    num_iters_plain      = gmres_handle->get_num_iters();

    EXPECT_GT(num_iters_plain, 0);
    EXPECT_LT(endRes, gmres_handle->get_tol());
    EXPECT_EQ(conv_flag, GMRESHandle::Flag::Conv);
  }

  // Solve Ax = b with LU preconditioner.
  {
    gmres_handle->reset_handle(m, tol);
    gmres_handle->set_verbose(verbose);

    // Make precond
    KokkosSparse::Experimental::LUPrec<sp_matrix_type, KernelHandle> myPrec(L, U);

    // reset X for next gmres call
    Kokkos::deep_copy(X, 0.0);

    gmres(&kh, A, B, X, &myPrec);

    // Double check residuals at end of solve:
    float_t nrmB = KokkosBlas::nrm2(B);
    KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
    KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
    float_t endRes = KokkosBlas::nrm2(B) / nrmB;

    const auto conv_flag = gmres_handle->get_conv_flag_val();
    num_iters_precond    = gmres_handle->get_num_iters();

    EXPECT_LT(endRes, gmres_handle->get_tol());
    EXPECT_EQ(conv_flag, GMRESHandle::Flag::Conv);
    EXPECT_LT(num_iters_precond, num_iters_plain);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_par_ilut_zerorow_A() {
  using RowMapType  = Kokkos::View<size_type*, device>;
  using EntriesType = Kokkos::View<lno_t*, device>;
  using ValuesType  = Kokkos::View<scalar_t*, device>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                       typename device::memory_space, typename device::memory_space>;

  const size_type nrows = 0;

  // Allocate device CRS views for A
  RowMapType row_map("row_map", 0);
  EntriesType entries("entries", 0);
  ValuesType values("values", 0);

  // Create host mirror views for CRS A
  auto hrow_map = Kokkos::create_mirror_view(row_map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);

  // Make kernel handle
  KernelHandle kh;

  kh.create_par_ilut_handle();

  auto par_ilut_handle = kh.get_par_ilut_handle();

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", nrows + 1);
  RowMapType U_row_map("U_row_map", nrows + 1);

  // Initial L/U approximations for A
  par_ilut_symbolic(&kh, row_map, entries, L_row_map, U_row_map);

  const size_type nnzL = par_ilut_handle->get_nnzL();
  const size_type nnzU = par_ilut_handle->get_nnzU();

  ASSERT_EQ(nnzL, 0);
  ASSERT_EQ(nnzU, 0);

  EntriesType L_entries("L_entries", nnzL);
  ValuesType L_values("L_values", nnzL);
  EntriesType U_entries("U_entries", nnzU);
  ValuesType U_values("U_values", nnzU);

  par_ilut_numeric(&kh, row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);

  const auto itrs        = par_ilut_handle->get_num_iters();
  const auto end_rel_res = par_ilut_handle->get_end_rel_res();

  EXPECT_EQ(itrs, 0);
  EXPECT_EQ(end_rel_res, scalar_t(0.));

  kh.destroy_par_ilut_handle();
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_par_ilut() {
  Test::run_test_par_ilut<scalar_t, lno_t, size_type, device>();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_par_ilut_precond() {
  Test::run_test_par_ilut_precond<scalar_t, lno_t, size_type, device>();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_par_ilut_zerorow_A() {
  Test::run_test_par_ilut_zerorow_A<scalar_t, lno_t, size_type, device>();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                  \
  TEST_F(TestCategory, sparse##_##par_ilut##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {           \
    test_par_ilut<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                \
  }                                                                                                  \
  TEST_F(TestCategory, sparse##_##par_ilut_zerorow_A##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_par_ilut_zerorow_A<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
  }                                                                                                  \
  TEST_F(TestCategory, sparse##_##par_ilut_precond##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {   \
    test_par_ilut_precond<SCALAR, ORDINAL, OFFSET, DEVICE>();                                        \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX
