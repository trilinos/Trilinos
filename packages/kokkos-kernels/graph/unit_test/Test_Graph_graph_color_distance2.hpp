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
#include <random>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_Distance2Color.hpp"
#include "KokkosGraph_MIS2.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

namespace Test {

// Verify that a distance-2 coloring is correct (all views must be hostspace)
template <typename lno_t, typename size_type, typename rowmap_t, typename entries_t, typename colors_t>
bool verifyD2Coloring(lno_t numVerts, const rowmap_t& rowmap, const entries_t& entries, const colors_t& colors) {
  // Just do the simplest possible neighbors-of-neighbors loop to find conflicts
  for (lno_t v = 0; v < numVerts; v++) {
    if (colors(v) == 0) {
      std::cout << "Vertex " << v << " is uncolored.\n";
      return false;
    }
    size_type rowBegin = rowmap(v);
    size_type rowEnd   = rowmap(v + 1);
    for (size_type i = rowBegin; i < rowEnd; i++) {
      lno_t nei1 = entries(i);
      if (nei1 < numVerts && nei1 != v) {
        // check for dist-1 conflict
        if (colors(v) == colors(nei1)) {
          std::cout << "Dist-1 conflict between " << v << " and " << nei1 << '\n';
          return false;
        }
        // iterate over dist-2 neighbors
        size_type colBegin = rowmap(nei1);
        size_type colEnd   = rowmap(nei1 + 1);
        for (size_type j = colBegin; j < colEnd; j++) {
          lno_t nei2 = entries(j);
          if (nei2 < numVerts && nei2 != v) {
            if (colors(v) == colors(nei2)) {
              std::cout << "Dist-2 conflict between " << v << " and " << nei2 << '\n';
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

template <typename lno_t, typename size_type, typename rowmap_t, typename entries_t, typename colors_t>
bool verifyBipartitePartialColoring(lno_t numRows, lno_t numCols, const rowmap_t& rowmap, const entries_t& entries,
                                    const rowmap_t& t_rowmap, const entries_t& t_entries, const colors_t& colors) {
  // Just do the simplest possible neighbors-of-neighbors loop to find conflicts
  for (lno_t v = 0; v < numRows; v++) {
    if (colors(v) == 0) {
      std::cout << "Vertex " << v << " is uncolored.\n";
      return false;
    }
    size_type rowBegin = rowmap(v);
    size_type rowEnd   = rowmap(v + 1);
    for (size_type i = rowBegin; i < rowEnd; i++) {
      lno_t nei1 = entries(i);
      if (nei1 < numCols) {
        // iterate over dist-2 neighbors
        size_type colBegin = t_rowmap(nei1);
        size_type colEnd   = t_rowmap(nei1 + 1);
        for (size_type j = colBegin; j < colEnd; j++) {
          lno_t nei2 = t_entries(j);
          if (nei2 < numRows && nei2 != v) {
            if (colors(v) == colors(nei2)) {
              std::cout << "Hyperedge conflict between " << v << " and " << nei2 << '\n';
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}
}  // namespace Test

template <typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_dist2_coloring(lno_t numVerts, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using execution_space = typename device::execution_space;
  using memory_space    = typename device::memory_space;
  using crsMat          = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type      = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t      = typename graph_type::row_map_type;
  using c_entries_t     = typename graph_type::entries_type;
  using rowmap_t        = typename c_rowmap_t::non_const_type;
  using entries_t       = typename c_entries_t::non_const_type;
  using KernelHandle    = KokkosKernelsHandle<size_type, lno_t, double, execution_space, memory_space, memory_space>;
  // Generate graph, and add some out-of-bounds columns
  crsMat A =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
  auto G = A.graph;
  // Symmetrize the graph
  rowmap_t symRowmap;
  entries_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<c_rowmap_t, c_entries_t, rowmap_t, entries_t, execution_space>(
      numVerts, G.row_map, G.entries, symRowmap, symEntries);
  auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symRowmap);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symEntries);
  std::vector<GraphColoringAlgorithmDistance2> algos = {COLORING_D2_DEFAULT, COLORING_D2_SERIAL,    COLORING_D2_VB,
                                                        COLORING_D2_VB_BIT,  COLORING_D2_VB_BIT_EF, COLORING_D2_NB_BIT};
  for (auto algo : algos) {
    KernelHandle kh;
    kh.create_distance2_graph_coloring_handle(algo);
    // Compute the Distance-2 graph coloring.
    graph_color_distance2<KernelHandle, c_rowmap_t, c_entries_t>(&kh, numVerts, symRowmap, symEntries);
    execution_space().fence();
    auto coloring_handle = kh.get_distance2_graph_coloring_handle();
    auto colors          = coloring_handle->get_vertex_colors();
    auto colorsHost      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);
    auto numColors       = coloring_handle->get_num_colors();
    EXPECT_LE(numColors, numVerts);
    bool success =
        Test::verifyD2Coloring<lno_t, size_type, decltype(rowmapHost), decltype(entriesHost), decltype(colorsHost)>(
            numVerts, rowmapHost, entriesHost, colorsHost);
    EXPECT_TRUE(success) << "Dist-2: algorithm " << coloring_handle->getD2AlgorithmName()
                         << " produced invalid coloring";
    kh.destroy_distance2_graph_coloring_handle();
  }
}

template <typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_bipartite_symmetric(lno_t numVerts, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using execution_space = typename device::execution_space;
  using memory_space    = typename device::memory_space;
  using crsMat          = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type      = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t      = typename graph_type::row_map_type;
  using c_entries_t     = typename graph_type::entries_type;
  using rowmap_t        = typename c_rowmap_t::non_const_type;
  using entries_t       = typename c_entries_t::non_const_type;
  using KernelHandle    = KokkosKernelsHandle<size_type, lno_t, double, execution_space, memory_space, memory_space>;
  // Generate graph, and add some out-of-bounds columns
  crsMat A =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
  auto G = A.graph;
  // Symmetrize the graph
  rowmap_t symRowmap;
  entries_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<c_rowmap_t, c_entries_t, rowmap_t, entries_t, execution_space>(
      numVerts, G.row_map, G.entries, symRowmap, symEntries);
  auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symRowmap);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symEntries);
  std::vector<GraphColoringAlgorithmDistance2> algos = {COLORING_D2_DEFAULT, COLORING_D2_SERIAL,    COLORING_D2_VB,
                                                        COLORING_D2_VB_BIT,  COLORING_D2_VB_BIT_EF, COLORING_D2_NB_BIT};
  for (auto algo : algos) {
    KernelHandle kh;
    kh.create_distance2_graph_coloring_handle(algo);
    // Compute the Distance-2 graph coloring.
    bipartite_color_rows<KernelHandle, c_rowmap_t, c_entries_t>(&kh, numVerts, numVerts, symRowmap, symEntries, true);
    execution_space().fence();
    auto coloring_handle = kh.get_distance2_graph_coloring_handle();
    auto colors          = coloring_handle->get_vertex_colors();
    auto colorsHost      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);
    auto numColors       = coloring_handle->get_num_colors();
    EXPECT_LE(numColors, numVerts);
    bool success = Test::verifyBipartitePartialColoring<lno_t, size_type, decltype(rowmapHost), decltype(entriesHost),
                                                        decltype(colorsHost)>(
        numVerts, numVerts, rowmapHost, entriesHost, rowmapHost, entriesHost, colorsHost);
    EXPECT_TRUE(success) << "Dist-2: algorithm " << coloring_handle->getD2AlgorithmName()
                         << " produced invalid coloring";
    kh.destroy_distance2_graph_coloring_handle();
  }
}

template <typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_bipartite(lno_t numRows, lno_t numCols, size_type nnz, lno_t bandwidth, lno_t row_size_variance,
                    bool colorRows) {
  using execution_space = typename device::execution_space;
  using memory_space    = typename device::memory_space;
  using crsMat          = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type      = typename crsMat::StaticCrsGraphType;
  using rowmap_t        = typename graph_type::row_map_type::non_const_type;
  using entries_t       = typename graph_type::entries_type::non_const_type;
  using c_rowmap_t      = typename graph_type::row_map_type;
  using c_entries_t     = typename graph_type::entries_type;
  using KernelHandle    = KokkosKernelsHandle<size_type, lno_t, double, execution_space, memory_space, memory_space>;
  // Generate graph
  crsMat A = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat>(numRows, numCols, nnz, row_size_variance, bandwidth);
  auto G   = A.graph;
  rowmap_t t_rowmap("rowmap^T", numCols + 1);
  entries_t t_entries("entries^T", G.entries.extent(0));
  KokkosSparse::Impl::transpose_graph<c_rowmap_t, c_entries_t, rowmap_t, entries_t, rowmap_t, execution_space>(
      numRows, numCols, G.row_map, G.entries, t_rowmap, t_entries);
  // TODO: remove me, shouldn't be needed even with UVM
  execution_space().fence();
  auto rowmapHost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), G.row_map);
  auto entriesHost   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), G.entries);
  auto t_rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), t_rowmap);
  auto t_entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), t_entries);
  std::vector<GraphColoringAlgorithmDistance2> algos = {COLORING_D2_DEFAULT, COLORING_D2_SERIAL,    COLORING_D2_VB,
                                                        COLORING_D2_VB_BIT,  COLORING_D2_VB_BIT_EF, COLORING_D2_NB_BIT};
  for (auto algo : algos) {
    KernelHandle kh;
    kh.create_distance2_graph_coloring_handle(algo);
    // Compute the one-sided bipartite coloring.
    if (colorRows) {
      bipartite_color_rows<KernelHandle, c_rowmap_t, c_entries_t>(&kh, numRows, numCols, G.row_map, G.entries);
    } else {
      bipartite_color_columns<KernelHandle, c_rowmap_t, c_entries_t>(&kh, numRows, numCols, G.row_map, G.entries);
    }
    execution_space().fence();
    auto coloring_handle = kh.get_distance2_graph_coloring_handle();
    auto colors          = coloring_handle->get_vertex_colors();
    auto colorsHost      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);
    auto numColors       = coloring_handle->get_num_colors();
    bool success;
    if (colorRows) {
      EXPECT_LE(numColors, numRows);
      success = Test::verifyBipartitePartialColoring<lno_t, size_type, decltype(rowmapHost), decltype(entriesHost),
                                                     decltype(colorsHost)>(numRows, numCols, rowmapHost, entriesHost,
                                                                           t_rowmapHost, t_entriesHost, colorsHost);
    } else {
      EXPECT_LE(numColors, numCols);
      success = Test::verifyBipartitePartialColoring<lno_t, size_type, decltype(rowmapHost), decltype(entriesHost),
                                                     decltype(colorsHost)>(
          numCols, numRows, t_rowmapHost, t_entriesHost, rowmapHost, entriesHost, colorsHost);
    }
    EXPECT_TRUE(success) << "Bipartite " << (colorRows ? "row" : "column") << " coloring: algorithm "
                         << coloring_handle->getD2AlgorithmName() << " produced invalid coloring";
    kh.destroy_distance2_graph_coloring_handle();
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                      \
  TEST_F(TestCategory, graph##_##graph_color_distance2##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {     \
    test_dist2_coloring<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 5000 * 20, 1000, 10);                       \
    test_dist2_coloring<SCALAR, ORDINAL, OFFSET, DEVICE>(50, 50 * 10, 40, 10);                             \
  }                                                                                                        \
  TEST_F(TestCategory, graph##_##graph_color_bipartite_sym##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_bipartite_symmetric<SCALAR, ORDINAL, OFFSET, DEVICE>(50, 50 * 5, 30, 1);                          \
    test_bipartite_symmetric<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 800, 10);                   \
  }                                                                                                        \
  TEST_F(TestCategory, graph##_##graph_color_bipartite_row##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_bipartite<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 4000, 3000 * 20, 800, 10, true);                 \
    test_bipartite<SCALAR, ORDINAL, OFFSET, DEVICE>(4000, 2000, 3000 * 20, 800, 10, true);                 \
  }                                                                                                        \
  TEST_F(TestCategory, graph##_##graph_color_bipartite_col##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_bipartite<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 4000, 3000 * 20, 800, 10, false);                \
    test_bipartite<SCALAR, ORDINAL, OFFSET, DEVICE>(4000, 2000, 3000 * 20, 800, 10, false);                \
  }

#if defined(KOKKOSKERNELS_INST_DOUBLE)
#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestDevice)
#endif
#endif

#undef EXECUTE_TEST
