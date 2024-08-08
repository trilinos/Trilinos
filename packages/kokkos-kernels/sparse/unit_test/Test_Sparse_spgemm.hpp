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

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_SortCrs.hpp"
// For Test::is_same_matrix
#include "Test_Sparse_Utils.hpp"
#include <string>
#include <stdexcept>

#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>

// This file contains the matrix for test_issue402
#include "matrixIssue402.hpp"

// const char *input_filename = "sherman1.mtx";
// const char *input_filename = "Si2.mtx";
// const char *input_filename = "wathen_30_30.mtx";
// const size_t expected_num_cols = 9906;

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

namespace Test {

// 3 ways to call SpGEMM:
// - symbolic/numeric with Views
// - symbolic/numeric with CrsMatrices
// - non-reuse with CrsMatrices
enum spgemm_call_mode { spgemm_reuse_view, spgemm_reuse_matrix, spgemm_noreuse };

// Randomize matrix values again from the same uniform distribution as
// kk_generate_sparse_matrix uses.
template <typename Values>
void randomize_matrix_values(const Values &v) {
  using ScalarType = typename Values::value_type;
  ScalarType randStart, randEnd;
  KokkosKernels::Impl::getRandomBounds(50.0, randStart, randEnd);
  Kokkos::Random_XorShift64_Pool<typename Values::execution_space> pool(13718);
  // Instead of sampling from [-50, 50] or [-50-50i, 50+50i],
  // sample from [1, 50] or [1+i, 50+50i]. That way relative
  // error between values can't become large if values happen to sum close to 0.
  Kokkos::fill_random(v, pool, randEnd / 50.0, randEnd);
}

template <typename crsMat_t>
void run_spgemm_noreuse(crsMat_t A, crsMat_t B, crsMat_t &C) {
  C = KokkosSparse::spgemm<crsMat_t>(A, false, B, false);
}

template <typename crsMat_t, typename device>
int run_spgemm(crsMat_t &A, crsMat_t &B, KokkosSparse::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &C, bool testReuse) {
  typedef typename crsMat_t::size_type size_type;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::value_type scalar_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                           typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  kh.create_spgemm_handle(spgemm_algorithm);
  {
    auto sh = kh.get_spgemm_handle();

    EXPECT_FALSE(sh->is_symbolic_called());
    EXPECT_FALSE(sh->is_numeric_called());
    EXPECT_FALSE(sh->are_rowptrs_computed());
    EXPECT_FALSE(sh->are_entries_computed());

    KokkosSparse::spgemm_symbolic(kh, A, false, B, false, C);

    EXPECT_TRUE(sh->is_symbolic_called());

    KokkosSparse::spgemm_numeric(kh, A, false, B, false, C);

    EXPECT_TRUE(sh->are_entries_computed());
    EXPECT_TRUE(sh->is_numeric_called());

    if (testReuse) {
      // Give A and B completely new random values (changing both the pointer
      // and contents), and re-run just numeric.
      A.values = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new A values"), A.nnz());
      B.values = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new B values"), B.nnz());
      randomize_matrix_values(A.values);
      randomize_matrix_values(B.values);
      KokkosSparse::spgemm_numeric(kh, A, false, B, false, C);
      EXPECT_TRUE(sh->are_entries_computed());
      EXPECT_TRUE(sh->is_numeric_called());
    }
  }
  kh.destroy_spgemm_handle();

  return 0;
}

template <typename crsMat_t, typename device>
int run_spgemm_old_interface(crsMat_t &A, crsMat_t &B, KokkosSparse::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &result,
                             bool testReuse) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                           typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  // kh.set_verbose(true);

  kh.create_spgemm_handle(spgemm_algorithm);
  {
    auto sh = kh.get_spgemm_handle();

    const size_t num_rows_A = A.numRows();
    const size_t num_rows_B = B.numRows();
    const size_t num_cols_B = B.numCols();

    lno_view_t row_mapC("non_const_lnow_row", num_rows_A + 1);
    lno_nnz_view_t entriesC;
    scalar_view_t valuesC;

    EXPECT_FALSE(sh->is_symbolic_called());
    EXPECT_FALSE(sh->is_numeric_called());
    EXPECT_FALSE(sh->are_rowptrs_computed());
    EXPECT_FALSE(sh->are_entries_computed());

    KokkosSparse::Experimental::spgemm_symbolic(&kh, num_rows_A, num_rows_B, num_cols_B, A.graph.row_map,
                                                A.graph.entries, false, B.graph.row_map, B.graph.entries, false,
                                                row_mapC);

    EXPECT_TRUE(sh->is_symbolic_called());

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    entriesC          = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC           = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
    KokkosSparse::Experimental::spgemm_numeric(&kh, num_rows_A, num_rows_B, num_cols_B, A.graph.row_map,
                                               A.graph.entries, A.values, false, B.graph.row_map, B.graph.entries,
                                               B.values, false, row_mapC, entriesC, valuesC);

    EXPECT_TRUE(sh->are_entries_computed());
    EXPECT_TRUE(sh->is_numeric_called());

    if (testReuse) {
      // Give A and B completely new random values (changing both the pointer
      // and contents), and re-run just numeric.
      A.values = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new A values"), A.nnz());
      B.values = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new B values"), B.nnz());
      randomize_matrix_values(A.values);
      randomize_matrix_values(B.values);
      KokkosSparse::Experimental::spgemm_numeric(&kh, num_rows_A, num_rows_B, num_cols_B, A.graph.row_map,
                                                 A.graph.entries, A.values, false, B.graph.row_map, B.graph.entries,
                                                 B.values, false, row_mapC, entriesC, valuesC);
      EXPECT_TRUE(sh->are_entries_computed());
      EXPECT_TRUE(sh->is_numeric_called());
    }

    graph_t static_graph(entriesC, row_mapC);
    result = crsMat_t("CrsMatrix", num_cols_B, valuesC, static_graph);
  }
  kh.destroy_spgemm_handle();
  return 0;
}
}  // namespace Test

// Generate matrices and test all supported spgemm algorithms.
// C := AB, where A is m*k, B is k*n, and C is m*n.
template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spgemm(lno_t m, lno_t k, lno_t n, size_type nnz, lno_t bandwidth, lno_t row_size_variance,
                 Test::spgemm_call_mode callMode, bool testReuse = false) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  {
    std::cerr << "TEST SKIPPED: See "
                 "https://github.com/kokkos/kokkos-kernels/issues/1542 for details."
              << std::endl;
    return;
  }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
  using namespace Test;
  // device::execution_space::initialize();
  // device::execution_space::print_configuration(std::cout);

  typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  // typedef typename crsMat_t::StaticCrsGraphType graph_t;
  // typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  // typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  // typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  // Generate random compressed sparse row matrix. Randomly generated (non-zero)
  // values are stored in a 1-D (1 rank) array.
  crsMat_t A = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, k, nnz, row_size_variance, bandwidth);
  crsMat_t B = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(k, n, nnz, row_size_variance, bandwidth);
  randomize_matrix_values(A.values);
  randomize_matrix_values(B.values);

  KokkosSparse::sort_crs_matrix(A);
  KokkosSparse::sort_crs_matrix(B);

  crsMat_t output_mat2;
  // If this test won't reuse symbolic, we can compute the reference matrix once
  // here
  if (!testReuse) {
    run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, output_mat2, false);
  }

  std::vector<SPGEMMAlgorithm> algorithms;
  if (callMode == spgemm_noreuse) {
    // No-reuse interface always uses the default algorithm
    algorithms = {SPGEMM_KK};
  } else {
    algorithms = {
        SPGEMM_KK, SPGEMM_KK_LP, SPGEMM_KK_MEMORY /* alias SPGEMM_KK_MEMSPEED */,
        SPGEMM_KK_SPEED /* alias SPGEMM_KK_DENSE */
    };
  }

  for (auto spgemm_algorithm : algorithms) {
    std::string algo         = "UNKNOWN";
    bool is_expected_to_fail = false;

    switch (spgemm_algorithm) {
      case SPGEMM_KK: algo = "SPGEMM_KK"; break;
      case SPGEMM_KK_LP: algo = "SPGEMM_KK_LP"; break;
      case SPGEMM_KK_MEMSPEED: algo = "SPGEMM_KK_MEMSPEED"; break;
      case SPGEMM_KK_SPEED: algo = "SPGEMM_KK_SPEED"; break;
      case SPGEMM_KK_MEMORY: algo = "SPGEMM_KK_MEMORY"; break;
      default: algo = "!!! UNKNOWN ALGO !!!";
    }

    Kokkos::Timer timer1;
    crsMat_t output_mat;

    bool failed = false;
    int res     = 0;
    try {
      switch (callMode) {
        case spgemm_reuse_view:
          res = run_spgemm_old_interface<crsMat_t, device>(A, B, spgemm_algorithm, output_mat, testReuse);
          break;
        case spgemm_reuse_matrix:
          res = run_spgemm<crsMat_t, device>(A, B, spgemm_algorithm, output_mat, testReuse);
          break;
        case spgemm_noreuse: run_spgemm_noreuse(A, B, output_mat); break;
      }
    } catch (const char *message) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << message;
      failed = true;
    } catch (std::string message) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << message;
      failed = true;
    } catch (std::exception &e) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << e.what();
      failed = true;
    }
    EXPECT_EQ(is_expected_to_fail, failed);

    // If this is testing reuse, the values of A and B changed so
    // the reference matrix must be recomputed
    if (testReuse) {
      run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, output_mat2, false);
    }

    // double spgemm_time = timer1.seconds();

    timer1.reset();
    if (!is_expected_to_fail) {
      EXPECT_TRUE((res == 0)) << algo;
      bool is_identical = is_same_matrix<crsMat_t, device>(output_mat, output_mat2);
      EXPECT_TRUE(is_identical) << algo;
      // EXPECT_TRUE( equal) << algo;
    }
    // std::cout << "algo:" << algo << " spgemm_time:" << spgemm_time << "
    // output_check_time:" << timer1.seconds() << std::endl;
  }
  // device::execution_space::finalize();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spgemm_symbolic(bool callSymbolicFirst, bool testEmpty) {
  using crsMat_t       = CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using graph_t        = typename crsMat_t::StaticCrsGraphType;
  using values_t       = typename crsMat_t::values_type;
  using entries_t      = typename graph_t::entries_type;
  using rowmap_t       = typename graph_t::row_map_type::non_const_type;
  using const_rowmap_t = typename graph_t::row_map_type;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                       typename device::memory_space, typename device::memory_space>;
  // A is m*n, B is n*k, C is m*k
  int m = 100;
  int n = 300;
  int k = 200;
  crsMat_t A, B;
  // Target 1000 total nonzeros in both A and B.
  if (testEmpty) {
    // Create A,B with the same dimensions, but zero entries
    values_t emptyValues;
    entries_t emptyEntries;
    // Initialize these to 0
    rowmap_t A_rowmap("A rowmap", m + 1);
    rowmap_t B_rowmap("B rowmap", n + 1);
    A = crsMat_t("A", m, n, 0, emptyValues, A_rowmap, emptyEntries);
    B = crsMat_t("B", n, k, 0, emptyValues, B_rowmap, emptyEntries);
  } else {
    size_type nnz = 1000;
    A             = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, n, nnz, 10, 50);
    nnz           = 1000;
    B             = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(n, k, nnz, 10, 50);
    KokkosSparse::sort_crs_matrix(A);
    KokkosSparse::sort_crs_matrix(B);
  }
  // Call reference impl to get complete product
  crsMat_t C_reference;
  Test::run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, C_reference, false);
  // Now call just symbolic, and specifically request that rowptrs be populated
  // Make sure this never depends on C_rowmap being initialized
  rowmap_t C_rowmap(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmapC"), m + 1);
  Kokkos::deep_copy(C_rowmap, size_type(123));
  KernelHandle kh;
  kh.create_spgemm_handle();
  if (callSymbolicFirst) {
    KokkosSparse::Experimental::spgemm_symbolic(&kh, m, n, k, A.graph.row_map, A.graph.entries, false, B.graph.row_map,
                                                B.graph.entries, false, C_rowmap);
  }
  KokkosSparse::Experimental::spgemm_symbolic(&kh, m, n, k, A.graph.row_map, A.graph.entries, false, B.graph.row_map,
                                              B.graph.entries, false, C_rowmap, true);
  kh.destroy_spgemm_handle();
  bool isCorrect = KokkosKernels::Impl::kk_is_identical_view<const_rowmap_t, const_rowmap_t, size_type,
                                                             typename device::execution_space>(
      C_rowmap, C_reference.graph.row_map, 0);
  EXPECT_TRUE(isCorrect) << " spgemm_symbolic produced incorrect rowptrs - callSymbolicFirst = " << callSymbolicFirst
                         << ", empty A/B = " << testEmpty;
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_issue402() {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  {
    std::cerr << "TEST SKIPPED: See "
                 "https://github.com/kokkos/kokkos-kernels/issues/1542 for details."
              << std::endl;
    return;
  }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
  using namespace Test;
  typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;

  // this specific matrix (from a circuit simulation) reliably replicated issue
  // #402 (incorrect/crashing SPGEMM KKMEM)
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  const lno_t numRows = 1813;
  const size_type nnz = 11156;
  lno_view_t Arowmap("A rowmap", numRows + 1);
  lno_nnz_view_t Aentries("A entries", nnz);
  scalar_view_t Avalues("A values", nnz);
  // Read out the matrix from the header file "matrixIssue402.hpp"
  {
    auto rowmapHost  = Kokkos::create_mirror_view(Arowmap);
    auto entriesHost = Kokkos::create_mirror_view(Aentries);
    auto valuesHost  = Kokkos::create_mirror_view(Avalues);
    for (lno_t i = 0; i < numRows + 1; i++) rowmapHost(i) = MatrixIssue402::rowmap[i];
    for (size_type i = 0; i < nnz; i++) {
      entriesHost(i) = MatrixIssue402::entries[i];
      valuesHost(i)  = MatrixIssue402::values[i];
    }
    Kokkos::deep_copy(Arowmap, rowmapHost);
    Kokkos::deep_copy(Aentries, entriesHost);
    Kokkos::deep_copy(Avalues, valuesHost);
  }
  crsMat_t A("A", numRows, numRows, nnz, Avalues, Arowmap, Aentries);
  // compute explicit transpose: the bug was replicated by computing AA'
  lno_view_t Browmap("B = A^T rowmap", numRows + 1);
  lno_nnz_view_t Bentries("B = A^T entries", nnz);
  scalar_view_t Bvalues("B = A^T values", nnz);
  KokkosSparse::Impl::transpose_matrix<lno_view_t, lno_nnz_view_t, scalar_view_t, lno_view_t, lno_nnz_view_t,
                                       scalar_view_t, lno_view_t, typename device::execution_space>(
      numRows, numRows, Arowmap, Aentries, Avalues, Browmap, Bentries, Bvalues);
  crsMat_t B("B=A^T", numRows, numRows, nnz, Bvalues, Browmap, Bentries);
  KokkosSparse::sort_crs_matrix(A);
  KokkosSparse::sort_crs_matrix(B);
  crsMat_t Cgold;
  run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, Cgold, false);
  crsMat_t C;
  bool success = true;
  std::string errMsg;
  try {
    int res = run_spgemm<crsMat_t, device>(A, B, SPGEMM_KK_MEMORY, C, false);
    if (res) throw "run_spgemm returned error code";
  } catch (const char *message) {
    errMsg  = message;
    success = false;
  } catch (std::string message) {
    errMsg  = message;
    success = false;
  } catch (std::exception &e) {
    errMsg  = e.what();
    success = false;
  }
  EXPECT_TRUE(success) << "SpGEMM still has issue 402 bug! Error message:\n" << errMsg << '\n';
  bool correctResult = is_same_matrix<crsMat_t, device>(C, Cgold);
  EXPECT_TRUE(correctResult) << "SpGEMM still has issue 402 bug; C=AA' is incorrect!\n";
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_issue1738() {
  // Make sure that std::invalid_argument is thrown if you:
  //  - call numeric where an input matrix's entries have changed.
  //  - try to reuse an spgemm handle by calling symbolic with new input
  //  matrices
  // This check is only enabled in debug builds.
#ifndef NDEBUG
  using crsMat_t = CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                                       typename device::memory_space, typename device::memory_space>;
  crsMat_t A1 = KokkosSparse::Impl::kk_generate_diag_matrix<crsMat_t>(100);
  crsMat_t B1 = KokkosSparse::Impl::kk_generate_diag_matrix<crsMat_t>(100);
  crsMat_t A2 = KokkosSparse::Impl::kk_generate_diag_matrix<crsMat_t>(50);
  crsMat_t B2 = KokkosSparse::Impl::kk_generate_diag_matrix<crsMat_t>(50);
  {
    KernelHandle kh;
    kh.create_spgemm_handle();
    crsMat_t C1;
    KokkosSparse::spgemm_symbolic(kh, A1, false, B1, false, C1);
    KokkosSparse::spgemm_numeric(kh, A1, false, B1, false, C1);
    crsMat_t C2;
    EXPECT_THROW(KokkosSparse::spgemm_symbolic(kh, A2, false, B2, false, C2), std::invalid_argument);
  }
  {
    KernelHandle kh;
    kh.create_spgemm_handle();
    crsMat_t C1;
    KokkosSparse::spgemm_symbolic(kh, A1, false, B1, false, C1);
    // Note: A1 is a 100x100 diagonal matrix, so the first entry in the first
    // row is 0. Change it to a 1 and make sure spgemm_numeric notices that it
    // changed.
    Kokkos::deep_copy(Kokkos::subview(A1.graph.entries, 0), 1);
    EXPECT_THROW(KokkosSparse::spgemm_numeric(kh, A1, false, B1, false, C1), std::invalid_argument);
  }
#endif
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                   \
  TEST_F(TestCategory, sparse##_##spgemm##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {                              \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 8000, 6000, 8000 * 20, 500, 10, ::Test::spgemm_reuse_matrix); \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 8000, 6000, 8000 * 20, 500, 10, ::Test::spgemm_reuse_view);   \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 500, 1600, 1000 * 20, 500, 10, ::Test::spgemm_reuse_matrix,    \
                                                 true);                                                               \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 500, 1600, 1000 * 20, 500, 10, ::Test::spgemm_reuse_view,      \
                                                 true);                                                               \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 0, 0, 0, 10, 10, ::Test::spgemm_reuse_matrix);                    \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 0, 0, 0, 10, 10, ::Test::spgemm_reuse_view);                      \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 12, 5, 0, 10, 0, ::Test::spgemm_reuse_matrix);                    \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 12, 5, 0, 10, 0, ::Test::spgemm_reuse_view);                      \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 0, 0, 10, 10, ::Test::spgemm_reuse_matrix);                  \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 0, 0, 10, 10, ::Test::spgemm_reuse_view);                    \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 10, 0, 0, 0, ::Test::spgemm_reuse_matrix);                   \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 10, 0, 0, 0, ::Test::spgemm_reuse_view);                     \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 8000, 6000, 8000 * 20, 500, 10, ::Test::spgemm_noreuse);      \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 500, 1600, 1000 * 20, 500, 10, ::Test::spgemm_noreuse);        \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 0, 0, 0, 10, 10, ::Test::spgemm_noreuse);                         \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 12, 5, 0, 10, 0, ::Test::spgemm_noreuse);                         \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 0, 0, 10, 10, ::Test::spgemm_noreuse);                       \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 10, 0, 0, 0, ::Test::spgemm_noreuse);                        \
    test_spgemm_symbolic<SCALAR, ORDINAL, OFFSET, DEVICE>(true, true);                                                \
    test_spgemm_symbolic<SCALAR, ORDINAL, OFFSET, DEVICE>(false, true);                                               \
    test_spgemm_symbolic<SCALAR, ORDINAL, OFFSET, DEVICE>(true, false);                                               \
    test_spgemm_symbolic<SCALAR, ORDINAL, OFFSET, DEVICE>(false, false);                                              \
    test_issue402<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                                 \
    test_issue1738<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                                \
  }

// test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 100, 10);
// test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 200, 10);

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
