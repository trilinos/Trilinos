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
#include <string>
#include <stdexcept>

#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>

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

template <typename crsMat_t, typename device, typename scalar_type, typename dinv_view_t>
int run_spgemm_jacobi(crsMat_t input_mat, crsMat_t input_mat2, scalar_type omega, dinv_view_t dinv,
                      KokkosSparse::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &result) {
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

  spgemm_symbolic(&kh, num_rows_1, num_rows_2, num_cols_2, input_mat.graph.row_map, input_mat.graph.entries, false,
                  input_mat2.graph.row_map, input_mat2.graph.entries, false, row_mapC);

  size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size) {
    entriesC = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC  = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
  }
  spgemm_jacobi(&kh, num_rows_1, num_rows_2, num_cols_2, input_mat.graph.row_map, input_mat.graph.entries,
                input_mat.values, false, input_mat2.graph.row_map, input_mat2.graph.entries, input_mat2.values, false,
                row_mapC, entriesC, valuesC, omega, dinv);

  graph_t static_graph(entriesC, row_mapC);
  crsMat_t crsmat("CrsMatrix", num_cols_2, valuesC, static_graph);
  result = crsmat;
  kh.destroy_spgemm_handle();

  return 0;
}

template <typename crsMat_t, typename device>
bool is_same_mat(crsMat_t output_mat1, crsMat_t output_mat2) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows1    = output_mat1.graph.row_map.extent(0);
  size_t nentries1 = output_mat1.graph.entries.extent(0);
  size_t nvals1    = output_mat1.values.extent(0);

  size_t nrows2    = output_mat2.graph.row_map.extent(0);
  size_t nentries2 = output_mat2.graph.entries.extent(0);
  size_t nvals2    = output_mat2.values.extent(0);

  if (nrows1 != nrows2) {
    std::cout << "nrows1:" << nrows1 << " nrows2:" << nrows2 << std::endl;
    return false;
  }
  if (nentries1 != nentries2) {
    std::cout << "nentries1:" << nentries1 << " nentries2:" << nentries2 << std::endl;
    return false;
  }
  if (nvals1 != nvals2) {
    std::cout << "nvals1:" << nvals1 << " nvals2:" << nvals2 << std::endl;
    return false;
  }

  bool is_identical = true;
  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<typename graph_t::row_map_type, typename graph_t::row_map_type,
                                                typename lno_view_t::value_type, typename device::execution_space>(
          output_mat1.graph.row_map, output_mat2.graph.row_map, 0);

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat1.graph.row_map);
    KokkosKernels::Impl::kk_print_1Dview(output_mat2.graph.row_map);
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
                                                typename device::execution_space>(output_mat1.graph.entries,
                                                                                  output_mat2.graph.entries, 0);

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat1.graph.entries);
    KokkosKernels::Impl::kk_print_1Dview(output_mat2.graph.entries);
    return false;
  }

  typedef typename Kokkos::ArithTraits<typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type, float>::value ? 2 * 1e-3 : 1e-7;

  is_identical = KokkosKernels::Impl::kk_is_relatively_identical_view<scalar_view_t, scalar_view_t, eps_type,
                                                                      typename device::execution_space>(
      output_mat1.values, output_mat2.values, eps);

  if (!is_identical) {
    std::cout << "values are different for eps: " << eps << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat1.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat2.values);

    return false;
  }
  return true;
}
}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spgemm_jacobi(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
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

  lno_t numCols      = numRows;
  crsMat_t input_mat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);

  KokkosSparse::sort_crs_matrix(input_mat);

  crsMat_t output_mat2;
  scalar_t omega = 3.0;

  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename device::execution_space c_exec_t;
  typedef typename device::memory_space c_temp_t;
  typedef typename Kokkos::Device<c_exec_t, c_temp_t> UniformDevice_t;
  typedef
      typename Kokkos::View<scalar_t **, typename KokkosKernels::Impl::GetUnifiedLayout<scalar_view_t>::array_layout,
                            UniformDevice_t>
          view_t;

  view_t dinv("Dinv", numRows, 1);
  Kokkos::deep_copy(dinv, 2.0);

  run_spgemm_jacobi<crsMat_t, device, scalar_t, view_t>(input_mat, input_mat, omega, dinv, SPGEMM_SERIAL, output_mat2);

  SPGEMMAlgorithm spgemm_algorithm = SPGEMM_KK_MEMORY;  // should we test other SpGEMM algorithms as well?

  crsMat_t output_mat;

  run_spgemm_jacobi<crsMat_t, device>(input_mat, input_mat, omega, dinv, spgemm_algorithm, output_mat);
  // Sort the reference output_mat2, but not output_mat. It should already be
  // soted.
  KokkosSparse::sort_crs_matrix(output_mat2);
  bool is_identical = is_same_mat<crsMat_t, device>(output_mat, output_mat2);
  EXPECT_TRUE(is_identical);
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                             \
  TEST_F(TestCategory, sparse##_##spgemm_jacobi##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_spgemm_jacobi<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 10, 50, 10);               \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
