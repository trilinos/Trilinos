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
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_default_types.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

namespace Test {
template <typename crsMat_t, typename device>
int run_graphcolor_deter(crsMat_t input_mat, ColoringAlgorithm coloring_algorithm, size_t &num_colors,
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
void test_coloring_deterministic(lno_t numRows, size_type nnz) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::non_const_type color_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  // typedef typename lno_view_t::non_const_value_type size_type;

  lno_t numCols = numRows;

  typename lno_view_t::non_const_type xadj("xadj", numRows + 1);
  typename lno_view_t::non_const_type::HostMirror h_xadj = Kokkos::create_mirror_view(xadj);
  typename lno_nnz_view_t::non_const_type adj("adj", nnz);
  typename lno_nnz_view_t::non_const_type::HostMirror h_adj = Kokkos::create_mirror_view(adj);

  // Fill up the rowPtr array
  h_xadj(0)  = 0;
  h_xadj(1)  = 3;
  h_xadj(2)  = 7;
  h_xadj(3)  = 11;
  h_xadj(4)  = 14;
  h_xadj(5)  = 18;
  h_xadj(6)  = 23;
  h_xadj(7)  = 29;
  h_xadj(8)  = 33;
  h_xadj(9)  = 37;
  h_xadj(10) = 42;
  h_xadj(11) = 47;
  h_xadj(12) = 51;
  h_xadj(13) = 55;
  h_xadj(14) = 58;
  h_xadj(15) = 62;
  h_xadj(16) = 66;
  h_xadj(17) = 70;
  h_xadj(18) = 74;
  Kokkos::deep_copy(xadj, h_xadj);

  // Fill up the column indices array
  h_adj(0)  = 0;
  h_adj(1)  = 1;
  h_adj(2)  = 4;
  h_adj(3)  = 0;
  h_adj(4)  = 1;
  h_adj(5)  = 2;
  h_adj(6)  = 5;
  h_adj(7)  = 1;
  h_adj(8)  = 2;
  h_adj(9)  = 3;
  h_adj(10) = 6;
  h_adj(11) = 2;
  h_adj(12) = 3;
  h_adj(13) = 7;
  h_adj(14) = 0;
  h_adj(15) = 4;
  h_adj(16) = 5;
  h_adj(17) = 8;
  h_adj(18) = 1;
  h_adj(19) = 4;
  h_adj(20) = 5;
  h_adj(21) = 6;
  h_adj(22) = 9;
  h_adj(23) = 2;
  h_adj(24) = 5;
  h_adj(25) = 6;
  h_adj(26) = 7;
  h_adj(27) = 10;
  h_adj(28) = 12;
  h_adj(29) = 3;
  h_adj(30) = 6;
  h_adj(31) = 7;
  h_adj(32) = 17;
  h_adj(33) = 4;
  h_adj(34) = 8;
  h_adj(35) = 9;
  h_adj(36) = 13;
  h_adj(37) = 5;
  h_adj(38) = 8;
  h_adj(39) = 9;
  h_adj(40) = 10;
  h_adj(41) = 14;
  h_adj(42) = 6;
  h_adj(43) = 9;
  h_adj(44) = 10;
  h_adj(45) = 11;
  h_adj(46) = 15;
  h_adj(47) = 10;
  h_adj(48) = 11;
  h_adj(49) = 12;
  h_adj(50) = 16;
  h_adj(51) = 6;
  h_adj(52) = 11;
  h_adj(53) = 12;
  h_adj(54) = 17;
  h_adj(55) = 8;
  h_adj(56) = 13;
  h_adj(57) = 14;
  h_adj(58) = 9;
  h_adj(59) = 13;
  h_adj(60) = 14;
  h_adj(61) = 15;
  h_adj(62) = 10;
  h_adj(63) = 14;
  h_adj(64) = 15;
  h_adj(65) = 16;
  h_adj(66) = 11;
  h_adj(67) = 15;
  h_adj(68) = 16;
  h_adj(69) = 17;
  h_adj(70) = 7;
  h_adj(71) = 12;
  h_adj(72) = 16;
  h_adj(73) = 17;
  Kokkos::deep_copy(adj, h_adj);

  size_type numentries = adj.extent(0);
  scalar_view_t newValues("vals", numentries);

  graph_t static_graph(adj, xadj);
  crsMat_t input_mat("CrsMatrix", numCols, newValues, static_graph);

  std::vector<ColoringAlgorithm> coloring_algorithms;

  coloring_algorithms.push_back(COLORING_VBD);
  coloring_algorithms.push_back(COLORING_VBDBIT);

  for (size_t ii = 0; ii < coloring_algorithms.size(); ++ii) {
    ColoringAlgorithm coloring_algorithm = coloring_algorithms[ii];
    color_view_t vector_colors;
    size_t num_colors;

    Kokkos::Timer timer1;
    int res = run_graphcolor_deter<crsMat_t, device>(input_mat, coloring_algorithm, num_colors, vector_colors);
    EXPECT_TRUE((res == 0));

    EXPECT_TRUE((num_colors == 2));

    size_type num_conflict                            = 0;
    typename color_view_t::HostMirror h_vector_colors = Kokkos::create_mirror_view(vector_colors);
    Kokkos::deep_copy(h_vector_colors, vector_colors);
    int exact_colors[18] = {2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1};

    for (lno_t vertexIdx = 0; vertexIdx < numRows; ++vertexIdx) {
      if (h_vector_colors(vertexIdx) != exact_colors[vertexIdx]) {
        ++num_conflict;
      }
    }

    EXPECT_TRUE((num_conflict == 0));
    // device::execution_space::finalize();
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                      \
  TEST_F(TestCategory, graph##_##graph_color_deterministic##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_coloring_deterministic<SCALAR, ORDINAL, OFFSET, DEVICE>(18, 74);                                  \
    test_coloring_deterministic<SCALAR, ORDINAL, OFFSET, DEVICE>(18, 74);                                  \
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
