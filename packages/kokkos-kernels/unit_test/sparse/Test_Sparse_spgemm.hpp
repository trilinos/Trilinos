/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#include "KokkosKernels_SparseUtils.hpp"
#include "KokkosKernels_Sorting.hpp"
#include <Kokkos_Concepts.hpp>
#include <string>
#include <stdexcept>

#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <KokkosKernels_IOUtils.hpp>

// This file contains the matrix for test_issue402
#include "matrixIssue402.hpp"

// const char *input_filename = "sherman1.mtx";
// const char *input_filename = "Si2.mtx";
// const char *input_filename = "wathen_30_30.mtx";
// const size_t expected_num_cols = 9906;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

namespace Test {

template <typename crsMat_t, typename device>
int run_spgemm(crsMat_t A, crsMat_t B,
               KokkosSparse::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &C) {
  typedef typename crsMat_t::size_type size_type;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::value_type scalar_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, typename device::execution_space,
      typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  kh.create_spgemm_handle(spgemm_algorithm);

  KokkosSparse::spgemm_symbolic(kh, A, false, B, false, C);
  KokkosSparse::spgemm_numeric(kh, A, false, B, false, C);
  kh.destroy_spgemm_handle();

  return 0;
}

template <typename crsMat_t, typename device>
int run_spgemm_old_interface(crsMat_t input_mat, crsMat_t input_mat2,
                             KokkosSparse::SPGEMMAlgorithm spgemm_algorithm,
                             crsMat_t &result) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, typename device::execution_space,
      typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  // kh.set_verbose(true);

  kh.create_spgemm_handle(spgemm_algorithm);

  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_rows_2 = input_mat2.numRows();
  const size_t num_cols_2 = input_mat2.numCols();

  const size_t num_cols_1 = input_mat.numCols();
  bool equal              = num_rows_2 == num_cols_1;
  if (!equal) return 1;

  lno_view_t row_mapC("non_const_lnow_row", num_rows_1 + 1);
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  spgemm_symbolic(&kh, num_rows_1, num_rows_2, num_cols_2,
                  input_mat.graph.row_map, input_mat.graph.entries, false,
                  input_mat2.graph.row_map, input_mat2.graph.entries, false,
                  row_mapC);

  size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  entriesC          = lno_nnz_view_t(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
  valuesC = scalar_view_t(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
  spgemm_numeric(&kh, num_rows_1, num_rows_2, num_cols_2,
                 input_mat.graph.row_map, input_mat.graph.entries,
                 input_mat.values, false,

                 input_mat2.graph.row_map, input_mat2.graph.entries,
                 input_mat2.values, false, row_mapC, entriesC, valuesC);

  graph_t static_graph(entriesC, row_mapC);
  result = crsMat_t("CrsMatrix", num_cols_2, valuesC, static_graph);
  kh.destroy_spgemm_handle();

  return 0;
}
template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat_actual, crsMat_t output_mat_reference) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows_actual    = output_mat_actual.numRows();
  size_t nentries_actual = output_mat_actual.graph.entries.extent(0);
  size_t nvals_actual    = output_mat_actual.values.extent(0);

  size_t nrows_reference    = output_mat_reference.numRows();
  size_t nentries_reference = output_mat_reference.graph.entries.extent(0);
  size_t nvals_reference    = output_mat_reference.values.extent(0);

  if (nrows_actual != nrows_reference) {
    std::cout << "nrows_actual:" << nrows_actual
              << " nrows_reference:" << nrows_reference << std::endl;
    return false;
  }
  if (nentries_actual != nentries_reference) {
    std::cout << "nentries_actual:" << nentries_actual
              << " nentries_reference:" << nentries_reference << std::endl;
    return false;
  }
  if (nvals_actual != nvals_reference) {
    std::cout << "nvals_actual:" << nvals_actual
              << " nvals_reference:" << nvals_reference << std::endl;
    return false;
  }

  KokkosKernels::sort_crs_matrix(output_mat_actual);
  KokkosKernels::sort_crs_matrix(output_mat_reference);

  bool is_identical = true;
  is_identical      = KokkosKernels::Impl::kk_is_identical_view<
      typename graph_t::row_map_type, typename graph_t::row_map_type,
      typename lno_view_t::value_type, typename device::execution_space>(
      output_mat_actual.graph.row_map, output_mat_reference.graph.row_map, 0);

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    std::cout << "Actual rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.row_map);
    std::cout << "Correct rowmap (SPGEMM_DEBUG):\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.row_map);
    return false;
  }

  is_identical = KokkosKernels::Impl::kk_is_identical_view<
      lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(output_mat_actual.graph.entries,
                                        output_mat_reference.graph.entries, 0);

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.entries);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.entries);
    return false;
  }

  typedef typename Kokkos::Details::ArithTraits<
      typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type, float>::value ? 2 * 1e-3 : 1e-7;

  is_identical = KokkosKernels::Impl::kk_is_relatively_identical_view<
      scalar_view_t, scalar_view_t, eps_type, typename device::execution_space>(
      output_mat_actual.values, output_mat_reference.values, eps);

  if (!is_identical) {
    std::cout << "values are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.values);

    return false;
  }
  return true;
}
}  // namespace Test

// Generate matrices and test all supported spgemm algorithms.
// C := AB, where A is m*k, B is k*n, and C is m*n.
template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_spgemm(lno_t m, lno_t k, lno_t n, size_type nnz, lno_t bandwidth,
                 lno_t row_size_variance, bool oldInterface = false) {
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
  crsMat_t A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      m, k, nnz, row_size_variance, bandwidth);
  crsMat_t B = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      k, n, nnz, row_size_variance, bandwidth);

  crsMat_t output_mat2;
  if (oldInterface)
    run_spgemm_old_interface<crsMat_t, device>(A, B, SPGEMM_DEBUG, output_mat2);
  else
    run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, output_mat2);

  std::vector<SPGEMMAlgorithm> algorithms = {
      SPGEMM_KK, SPGEMM_KK_LP, SPGEMM_KK_MEMORY /* alias SPGEMM_KK_MEMSPEED */,
      SPGEMM_KK_SPEED /* alias SPGEMM_KK_DENSE */
  };

#ifdef HAVE_KOKKOSKERNELS_MKL
  algorithms.push_back(SPGEMM_MKL);
#endif

  for (auto spgemm_algorithm : algorithms) {
    const uint64_t max_integer = 2147483647;
    std::string algo           = "UNKNOWN";
    bool is_expected_to_fail   = false;

    switch (spgemm_algorithm) {
      case SPGEMM_CUSPARSE:
        // TODO: add these test failure cases for cusparse too.
        algo = "SPGEMM_CUSPARSE";
#if !defined(KERNELS_HAVE_CUSPARSE) && \
    !defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
        is_expected_to_fail = true;
#endif
        break;

      case SPGEMM_MKL:
        algo = "SPGEMM_MKL";
        // MKL requires scalar to be either float or double
        if (!(std::is_same<float, scalar_t>::value ||
              std::is_same<double, scalar_t>::value)) {
          is_expected_to_fail = true;
        }
        // mkl requires local ordinals to be int.
        if (!(std::is_same<int, lno_t>::value)) {
          is_expected_to_fail = true;
        }
        // if size_type is larger than int, mkl casts it to int.
        // it will fail if casting cause overflow.
        if (A.values.extent(0) > max_integer) {
          is_expected_to_fail = true;
        }

        if (!(Kokkos::SpaceAccessibility<
                typename Kokkos::HostSpace::execution_space,
                typename device::memory_space>::accessible)) {
          is_expected_to_fail = true;
        }
        break;

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
      if (oldInterface)
        res = run_spgemm_old_interface<crsMat_t, device>(A, B, spgemm_algorithm,
                                                         output_mat);
      else
        res = run_spgemm<crsMat_t, device>(A, B, spgemm_algorithm, output_mat);
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
    EXPECT_TRUE((failed == is_expected_to_fail));

    // double spgemm_time = timer1.seconds();

    timer1.reset();
    if (!is_expected_to_fail) {
      EXPECT_TRUE((res == 0)) << algo;
      bool is_identical =
          is_same_matrix<crsMat_t, device>(output_mat, output_mat2);
      EXPECT_TRUE(is_identical) << algo;
      // EXPECT_TRUE( equal) << algo;
    }
    // std::cout << "algo:" << algo << " spgemm_time:" << spgemm_time << "
    // output_check_time:" << timer1.seconds() << std::endl;
  }
  // device::execution_space::finalize();
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_issue402() {
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
    for (lno_t i = 0; i < numRows + 1; i++)
      rowmapHost(i) = MatrixIssue402::rowmap[i];
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
  KokkosKernels::Impl::transpose_matrix<
      lno_view_t, lno_nnz_view_t, scalar_view_t, lno_view_t, lno_nnz_view_t,
      scalar_view_t, lno_view_t, typename device::execution_space>(
      numRows, numRows, Arowmap, Aentries, Avalues, Browmap, Bentries, Bvalues);
  crsMat_t B("B=A^T", numRows, numRows, nnz, Bvalues, Browmap, Bentries);
  crsMat_t Cgold;
  run_spgemm<crsMat_t, device>(A, B, SPGEMM_DEBUG, Cgold);
  crsMat_t C;
  bool success = true;
  std::string errMsg;
  try {
    int res = run_spgemm<crsMat_t, device>(A, B, SPGEMM_KK_MEMORY, C);
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
  EXPECT_TRUE(success) << "KKMEM still has issue 402 bug! Error message:\n"
                       << errMsg << '\n';
  bool correctResult = is_same_matrix<crsMat_t, device>(C, Cgold);
  EXPECT_TRUE(correctResult)
      << "KKMEM still has issue 402 bug; C=AA' is incorrect!\n";
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                          \
  TEST_F(TestCategory,                                                         \
         sparse##_##spgemm##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {     \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 10000, 10000,          \
                                                 10000 * 20, 500, 10, false);  \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 10000, 10000,          \
                                                 10000 * 20, 500, 10, true);   \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 0, 0, 0, 10, 10, false);   \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 0, 0, 0, 10, 10, true);    \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 12, 5, 0, 10, 0, false);   \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(0, 12, 5, 0, 10, 0, true);    \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 0, 0, 10, 10, false); \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 0, 0, 10, 10, true);  \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 10, 0, 0, 0, false);  \
    test_spgemm<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 10, 10, 0, 0, 0, true);   \
    test_issue402<SCALAR, ORDINAL, OFFSET, DEVICE>();                          \
  }

// test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 100, 10);
// test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 200, 10);

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&         \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#undef EXECUTE_TEST
