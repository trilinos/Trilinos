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
    (!defined(TEST_CUDA_SPARSE_CPP) || (defined(TEST_CUDA_SPARSE_CPP) && defined(KOKKOS_ENABLE_CUDA_UVM)))

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

template <typename Crs, typename LUType, typename size_type,
          typename std::enable_if<is_crs_matrix<LUType>::value>::type* = nullptr>
LUType get_LU(char l_or_u, int n, size_type& nnz, int row_size_variance, int bandwidth, int) {
  auto LU =
      KokkosSparse::Impl::kk_generate_triangular_sparse_matrix<Crs>(l_or_u, n, n, nnz, row_size_variance, bandwidth);

  return LU;
}

template <typename Crs, typename LUType, typename size_type,
          typename std::enable_if<is_bsr_matrix<LUType>::value>::type* = nullptr>
LUType get_LU(char l_or_u, int n, size_type& nnz, int row_size_variance, int bandwidth, int block_size) {
  auto LU_unblocked =
      KokkosSparse::Impl::kk_generate_triangular_sparse_matrix<Crs>(l_or_u, n, n, nnz, row_size_variance, bandwidth);

  // Convert to BSR
  LUType LU(LU_unblocked, block_size);

  return LU;
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, typename device>
struct TrsvTest {
  using View2D          = Kokkos::View<scalar_t**, layout, device>;
  using execution_space = typename device::execution_space;

  using Crs = CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using Bsr = BsrMatrix<scalar_t, lno_t, device, void, size_type>;

  // TODO: remove this once MD develop branch is merge.
  // The below functionolity exists in SparseUtils.
  template <bool UseBlocks, typename sp_matrix_type>
  static void check_trsv_mv(sp_matrix_type input_mat, View2D x, View2D b, View2D expected_x, int numMV,
                            const char uplo[], const char trans[]) {
    double eps =
        (std::is_same<scalar_t, float>::value ? 2 * 1e-2
         : (std::is_same<scalar_t, std::complex<float>>::value || std::is_same<scalar_t, Kokkos::complex<float>>::value)
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

  template <bool UseBlocks>
  static void test_trsv_mv(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, int numMV) {
    using sp_matrix_type = std::conditional_t<UseBlocks, Bsr, Crs>;

    constexpr auto block_size = UseBlocks ? 10 : 1;

    lno_t numCols = numRows;

    View2D b_x("A", numRows, numMV);
    View2D b_y("B", numCols, numMV);
    View2D b_x_copy("B", numCols, numMV);

    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
    Kokkos::fill_random(b_x_copy, rand_pool, scalar_t(10));

    scalar_t alpha = 1;
    scalar_t beta  = 0;

    // this function creates a dense lower and upper triangular matrix.
    auto lower_part =
        get_LU<Crs, sp_matrix_type, size_type>('L', numRows, nnz, row_size_variance, bandwidth, block_size);

    Test::shuffleMatrixEntries(lower_part.graph.row_map, lower_part.graph.entries, lower_part.values, block_size);

    KokkosSparse::spmv("N", alpha, lower_part, b_x_copy, beta, b_y);
    check_trsv_mv<UseBlocks>(lower_part, b_x, b_y, b_x_copy, numMV, "L", "N");

    if (!UseBlocks) {
      KokkosSparse::spmv("T", alpha, lower_part, b_x_copy, beta, b_y);
      check_trsv_mv<UseBlocks>(lower_part, b_x, b_y, b_x_copy, numMV, "L", "T");
    }

    auto upper_part =
        get_LU<Crs, sp_matrix_type, size_type>('U', numRows, nnz, row_size_variance, bandwidth, block_size);

    Test::shuffleMatrixEntries(upper_part.graph.row_map, upper_part.graph.entries, upper_part.values, block_size);

    KokkosSparse::spmv("N", alpha, upper_part, b_x_copy, beta, b_y);
    check_trsv_mv<UseBlocks>(upper_part, b_x, b_y, b_x_copy, numMV, "U", "N");

    if (!UseBlocks) {
      KokkosSparse::spmv("T", alpha, upper_part, b_x_copy, beta, b_y);
      check_trsv_mv<UseBlocks>(upper_part, b_x, b_y, b_x_copy, numMV, "U", "T");
    }
  }
};

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename layout, typename device>
void test_trsv_mv() {
  using TestStruct = Test::TrsvTest<scalar_t, lno_t, size_type, layout, device>;
  TestStruct::template test_trsv_mv<false>(1000, 1000 * 30, 200, 10, 1);
  TestStruct::template test_trsv_mv<false>(800, 800 * 30, 100, 10, 5);
  TestStruct::template test_trsv_mv<false>(400, 400 * 20, 100, 5, 10);
  TestStruct::template test_trsv_mv<true>(1000, 1000 * 30, 200, 10, 1);
  TestStruct::template test_trsv_mv<true>(800, 800 * 30, 100, 10, 5);
  TestStruct::template test_trsv_mv<true>(400, 400 * 20, 100, 5, 10);
}

// Note BMK 7-22: the matrix generator used by this test always
// generates a dense triangle. It ignores bandwidth, nnz and row size variance.

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                                     \
  TEST_F(TestCategory, sparse##_##trsv_mv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_trsv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>();                                 \
  }

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutRight, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT

#undef EXECUTE_TEST_MV

#endif  // check for CUDA and UVM
