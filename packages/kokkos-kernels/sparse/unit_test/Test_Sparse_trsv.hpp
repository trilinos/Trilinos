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
#if !defined(TEST_HIP_SPARSE_CPP) && !defined(TEST_SYCL_SPARSE_CPP) && \
    !defined(TEST_OPENMPTARGET_BATCHED_DENSE_CPP) &&                   \
    (!defined(TEST_CUDA_SPARSE_CPP) ||                                 \
     (defined(TEST_CUDA_SPARSE_CPP) && defined(KOKKOS_ENABLE_CUDA_UVM)))

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosSparse_trsv.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>

#include <KokkosKernels_Utils.hpp>

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

namespace Test {
// TODO: remove this once MD develop branch is merge.
// The below functionolity exists in SparseUtils.

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_trsv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type b,
                   y_vector_type expected_x, int numMV, const char uplo[],
                   const char trans[]) {
  // typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename scalar_view_t::value_type ScalarA;
  double eps = (std::is_same<ScalarA, float>::value
                    ? 2 * 1e-2
                    : (std::is_same<ScalarA, std::complex<float>>::value ||
                       std::is_same<ScalarA, Kokkos::complex<float>>::value)
                          ? 2 * 1e-1
                          : 1e-7);

  Kokkos::fence();
  KokkosSparse::trsv(uplo, trans, "N", input_mat, b, x);

  for (int i = 0; i < numMV; ++i) {
    auto x_i = Kokkos::subview(x, Kokkos::ALL(), i);

    auto expected_x_i = Kokkos::subview(expected_x, Kokkos::ALL(), i);

    EXPECT_NEAR_KK_1DVIEW(expected_x_i, x_i, eps);
  }
}
}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type,
          typename layout, class Device>
void test_trsv_mv(lno_t numRows, size_type nnz, lno_t bandwidth,
                  lno_t row_size_variance, int numMV) {
  lno_t numCols = numRows;

  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;
  // typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeY;

  ViewTypeX b_x("A", numRows, numMV);
  ViewTypeY b_y("B", numCols, numMV);
  ViewTypeX b_x_copy("B", numCols, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);
  Kokkos::fill_random(b_x_copy, rand_pool, scalar_t(10));

  typename ViewTypeY::non_const_value_type alpha = 1;
  typename ViewTypeY::non_const_value_type beta  = 0;

  // this function creates a dense lower and upper triangular matrix.
  // TODO: SHOULD CHANGE IT TO SPARSE
  crsMat_t lower_part =
      KokkosSparse::Impl::kk_generate_triangular_sparse_matrix<crsMat_t>(
          'L', numRows, numCols, nnz, row_size_variance, bandwidth);

  Test::shuffleMatrixEntries(lower_part.graph.row_map, lower_part.graph.entries,
                             lower_part.values);

  KokkosSparse::spmv("N", alpha, lower_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(lower_part, b_x, b_y, b_x_copy, numMV, "L", "N");

  KokkosSparse::spmv("T", alpha, lower_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(lower_part, b_x, b_y, b_x_copy, numMV, "L", "T");
  // typedef typename Kokkos::View<lno_t*, layout, Device> indexview;

  crsMat_t upper_part =
      KokkosSparse::Impl::kk_generate_triangular_sparse_matrix<crsMat_t>(
          'U', numRows, numCols, nnz, row_size_variance, bandwidth);

  Test::shuffleMatrixEntries(upper_part.graph.row_map, upper_part.graph.entries,
                             upper_part.values);

  KokkosSparse::spmv("N", alpha, upper_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(upper_part, b_x, b_y, b_x_copy, numMV, "U", "N");

  KokkosSparse::spmv("T", alpha, upper_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(upper_part, b_x, b_y, b_x_copy, numMV, "U", "T");
}

// Note BMK 7-22: the matrix generator used by this test always
// generates a dense triangle. It ignores bandwidth, nnz and row size variance.

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                    \
  TEST_F(                                                                           \
      TestCategory,                                                                 \
      sparse##_##trsv_mv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_trsv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        1000, 1000 * 30, 200, 10, 1);                                               \
    test_trsv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        800, 800 * 30, 100, 10, 5);                                                 \
    test_trsv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        400, 400 * 20, 100, 5, 10);                                                 \
  }

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestExecSpace)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutRight, TestExecSpace)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT

#undef EXECUTE_TEST_MV

#endif  // check for CUDA and UVM
