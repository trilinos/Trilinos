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

#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_default_types.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

namespace Test {
template <typename crsMat_t, typename device>
int run_graphcolor(crsMat_t input_mat, ColoringAlgorithm coloring_algorithm, size_t &num_colors,
                   typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type &vertex_colors) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                              typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  kh.create_graph_coloring_handle(coloring_algorithm);

  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();

  graph_color<KernelHandle, lno_view_t, lno_nnz_view_t>(&kh, num_rows_1, num_cols_1, input_mat.graph.row_map,
                                                        input_mat.graph.entries);

  num_colors    = kh.get_graph_coloring_handle()->get_num_colors();
  vertex_colors = kh.get_graph_coloring_handle()->get_vertex_colors();
  kh.destroy_graph_coloring_handle();
  return 0;
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_coloring(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::non_const_type color_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  // typedef typename lno_view_t::non_const_value_type size_type;

  lno_t numCols = numRows;
  crsMat_t input_mat =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numCols, nnz, row_size_variance, bandwidth);

  typename lno_view_t::non_const_type sym_xadj;
  typename lno_nnz_view_t::non_const_type sym_adj;

  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<
      lno_view_t, lno_nnz_view_t, typename lno_view_t::non_const_type, typename lno_nnz_view_t::non_const_type,
      typename device::execution_space>(numRows, input_mat.graph.row_map, input_mat.graph.entries, sym_xadj, sym_adj);
  size_type numentries = sym_adj.extent(0);
  scalar_view_t newValues("vals", numentries);

  graph_t static_graph(sym_adj, sym_xadj);
  input_mat = crsMat_t("CrsMatrix", numCols, newValues, static_graph);

  std::vector<ColoringAlgorithm> coloring_algorithms = {COLORING_DEFAULT, COLORING_SERIAL, COLORING_VB, COLORING_VBBIT,
                                                        COLORING_VBCS};

  // FIXME: VBD sometimes fails on CUDA and HIP
#if defined(KOKKOS_ENABLE_CUDA)
  if (!std::is_same<typename device::execution_space, Kokkos::Cuda>::value) {
    coloring_algorithms.push_back(COLORING_VBD);
  }
#elif defined(KOKKOS_ENABLE_HIP)
  if (!std::is_same<typename device::execution_space, Kokkos::HIP>::value) {
    coloring_algorithms.push_back(COLORING_VBD);
  }
#else
  coloring_algorithms.push_back(COLORING_VBD);
#endif

  // FIXME SYCL: re-enable this when EB is working
#ifdef KOKKOS_ENABLE_SYCL
  if (!std::is_same<typename device::execution_space, Kokkos::Experimental::SYCL>::value) {
    coloring_algorithms.push_back(COLORING_EB);
  }
#else
  coloring_algorithms.push_back(COLORING_EB);
#endif

  for (size_t ii = 0; ii < coloring_algorithms.size(); ++ii) {
    ColoringAlgorithm coloring_algorithm = coloring_algorithms[ii];
    color_view_t vector_colors;
    size_t num_colors;

    Kokkos::Timer timer1;
    crsMat_t output_mat;
    int res = run_graphcolor<crsMat_t, device>(input_mat, coloring_algorithm, num_colors, vector_colors);
    // double coloring_time = timer1.seconds();
    EXPECT_TRUE((res == 0));

    const lno_t num_rows_1 = input_mat.numRows();
    const lno_t num_cols_1 = input_mat.numCols();
    lno_t num_conflict     = KokkosSparse::Impl::kk_is_d1_coloring_valid<lno_view_t, lno_nnz_view_t, color_view_t,
                                                                     typename device::execution_space>(
        num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, vector_colors);

    lno_t conf = 0;
    {
      // also check the correctness of the validation code :)
      typename lno_view_t::HostMirror hrm          = Kokkos::create_mirror_view(input_mat.graph.row_map);
      typename lno_nnz_view_t::HostMirror hentries = Kokkos::create_mirror_view(input_mat.graph.entries);
      typename color_view_t::HostMirror hcolor     = Kokkos::create_mirror_view(vector_colors);
      Kokkos::deep_copy(hrm, input_mat.graph.row_map);
      Kokkos::deep_copy(hentries, input_mat.graph.entries);
      Kokkos::deep_copy(hcolor, vector_colors);

      for (lno_t i = 0; i < num_rows_1; ++i) {
        const size_type b = hrm(i);
        const size_type e = hrm(i + 1);
        for (size_type j = b; j < e; ++j) {
          lno_t d = hentries(j);
          if (i != d) {
            if (hcolor(d) == hcolor(i)) {
              conf++;
            }
          }
        }
      }
    }
    EXPECT_TRUE((num_conflict == conf)) << "Coloring algo " << (int)coloring_algorithm
                                        << ": kk_is_d1_coloring_valid returned incorrect number of conflicts ("
                                        << num_conflict << ", should be " << conf << ")";

    EXPECT_TRUE((num_conflict == 0)) << "Coloring algo " << (int)coloring_algorithm
                                     << ": D1 coloring produced invalid coloring (" << num_conflict << " conflicts)";
  }
  // device::execution_space::finalize();
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                        \
  TEST_F(TestCategory, graph##_##graph_color##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_coloring<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 30, 200, 10);              \
    test_coloring<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 30, 100, 10);              \
  }

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(default_scalar, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(default_scalar, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(default_scalar, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(default_scalar, int64_t, size_t, TestDevice)
#endif

#undef EXECUTE_TEST
