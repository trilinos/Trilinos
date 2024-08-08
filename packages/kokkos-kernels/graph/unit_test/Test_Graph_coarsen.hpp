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
#include <set>
#include <list>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_CoarsenConstruct.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

// namespace Test {

template <class coarsener_t>
bool verify_coarsening(typename coarsener_t::coarse_level_triple fine_l,
                       typename coarsener_t::coarse_level_triple coarse_l) {
  using crsMat      = typename coarsener_t::matrix_t;
  using graph_type  = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t  = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t    = typename c_rowmap_t::non_const_type;
  using entries_t   = typename c_entries_t::non_const_type;
  using svt         = typename coarsener_t::wgt_view_t;
  using ordinal_t   = typename entries_t::value_type;
  using edge_t      = typename rowmap_t::value_type;

  crsMat A         = fine_l.mtx;
  crsMat coarse_A  = coarse_l.mtx;
  auto f_rowmap    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto c_rowmap    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarse_A.graph.row_map);
  auto f_entries   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto vcmap       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarse_l.interp_mtx.graph.entries);
  auto few         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  auto cew         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarse_A.values);
  auto fvw         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), fine_l.vtx_wgts);
  auto cvw         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarse_l.vtx_wgts);
  ordinal_t f_size = 0;
  ordinal_t c_size = 0;
  for (ordinal_t i = 0; i < static_cast<ordinal_t>(fvw.extent(0)); i++) {
    f_size += fvw(i);
  }
  for (ordinal_t i = 0; i < static_cast<ordinal_t>(cvw.extent(0)); i++) {
    c_size += cvw(i);
  }
  // number of columns in interpolation matrix should give number of rows in
  // coarse matrix
  if (coarse_l.interp_mtx.numCols() != coarse_l.mtx.numRows()) {
    return false;
  }
  // sum of vertex weights in each graph should be equal
  if (f_size != c_size) {
    return false;
  }
  typename svt::value_type f_edges = 0, c_edges = 0;
  for (ordinal_t i = 0; i < A.numRows(); i++) {
    for (edge_t j = f_rowmap(i); j < f_rowmap(i + 1); j++) {
      ordinal_t v = f_entries(j);
      if (vcmap(i) != vcmap(v)) {
        f_edges += few(j);
      }
    }
  }
  for (ordinal_t i = 0; i < coarse_A.numRows(); i++) {
    for (edge_t j = c_rowmap(i); j < c_rowmap(i + 1); j++) {
      c_edges += cew(j);
    }
  }
  // sum of inter-aggregate edges in fine graph should be sum of all edges in
  // coarse graph
  if (f_edges != c_edges) {
    return false;
  }
  return true;
}

template <class crsMat>
bool verify_is_graph(crsMat A) {
  using graph_type  = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t  = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t    = typename c_rowmap_t::non_const_type;
  using entries_t   = typename c_entries_t::non_const_type;
  using ordinal_t   = typename entries_t::value_type;
  using edge_t      = typename rowmap_t::value_type;
  auto rowmap       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto entries      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

  for (ordinal_t i = 0; i < A.numRows(); i++) {
    std::set<ordinal_t> adjset;
    for (edge_t j = rowmap(i); j < rowmap(i + 1); j++) {
      ordinal_t v = entries(j);
      // A should not contain out-of-bounds columns
      if (v >= A.numRows()) {
        return false;
      }
      // Each row should not contain duplicate columns
      if (adjset.find(v) != adjset.end()) {
        return false;
      }
      adjset.insert(v);
    }
  }
  return true;
}

template <class crsMat>
bool verify_aggregator(crsMat A, crsMat agg) {
  using graph_type  = typename crsMat::StaticCrsGraphType;
  using c_entries_t = typename graph_type::entries_type;
  using entries_t   = typename c_entries_t::non_const_type;
  using ordinal_t   = typename entries_t::value_type;

  // aggregator should have as many rows as A
  if (A.numRows() != agg.numRows()) {
    return false;
  }
  // aggregator should have as many entries as A has rows
  if (A.numRows() != static_cast<ordinal_t>(agg.nnz())) {
    return false;
  }
  // column count gives number of vertices in coarse graph
  // aggregator should not have more columns than A has rows
  // it is valid for coarse graph to be same size as input graph
  // some graphs (for example a graph with only self-loops) may produce
  // coarsenings that don't shrink the graph
  if (A.numRows() < agg.numCols()) {
    return false;
  }
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), agg.graph.entries);

  std::vector<int> aggregateSizes(agg.numCols(), 0);
  for (ordinal_t i = 0; i < static_cast<ordinal_t>(agg.nnz()); i++) {
    ordinal_t v = entries(i);
    // aggregator should not have out-of-bounds columns
    if (v >= agg.numCols()) {
      return false;
    }
    aggregateSizes[v]++;
  }
  for (ordinal_t i = 0; i < agg.numCols(); i++) {
    // each aggregate label should contain at least one fine vertex
    if (aggregateSizes[i] == 0) {
      return false;
    }
  }
  return true;
}

template <class crsMat>
crsMat gen_grid() {
  using graph_type  = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t  = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t    = typename c_rowmap_t::non_const_type;
  using entries_t   = typename c_entries_t::non_const_type;
  using svt         = typename crsMat::values_type;
  using lno_t       = typename crsMat::ordinal_type;
  using scalar      = typename crsMat::value_type;
  lno_t rows        = 200;
  lno_t cols        = 300;
  lno_t total_vtx   = rows * cols;

  // create a 2D-mesh
  rowmap_t row_map("rowmap", total_vtx + 1);
  auto rm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_map);
  rm(0)   = 0;
  for (lno_t i = 0; i < rows; i++) {
    for (lno_t j = 0; j < cols; j++) {
      lno_t vtx_id     = i * cols + j;
      lno_t edge_count = 0;
      if (i > 0) edge_count++;
      if (i + 1 < rows) edge_count++;
      if (j > 0) edge_count++;
      if (j + 1 < cols) edge_count++;
      rm(vtx_id + 1) = rm(vtx_id) + edge_count;
    }
  }
  lno_t total_edges = rm(total_vtx);
  Kokkos::deep_copy(row_map, rm);
  entries_t entries("entries", total_edges);
  auto e = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries);

  for (lno_t i = 0; i < rows; i++) {
    for (lno_t j = 0; j < cols; j++) {
      lno_t vtx_id       = i * cols + j;
      lno_t write_offset = rm(vtx_id);
      if (i > 0) {
        e(write_offset) = vtx_id - cols;
        write_offset++;
      }
      if (i + 1 < rows) {
        e(write_offset) = vtx_id + cols;
        write_offset++;
      }
      if (j > 0) {
        e(write_offset) = vtx_id - 1;
        write_offset++;
      }
      if (j + 1 < cols) {
        e(write_offset) = vtx_id + 1;
      }
    }
  }
  Kokkos::deep_copy(entries, e);
  graph_type g(entries, row_map);
  svt values("values", total_edges);
  Kokkos::deep_copy(values, static_cast<scalar>(1));
  crsMat A("A", total_vtx, values, g);
  return A;
}

template <typename scalar, typename lno_t, typename size_type, typename device>
void test_multilevel_coarsen_grid() {
  using crsMat      = KokkosSparse::CrsMatrix<scalar, lno_t, device, void, size_type>;
  crsMat A          = gen_grid<crsMat>();
  using coarsener_t = coarse_builder<crsMat>;
  typename coarsener_t::coarsen_handle handle;
  using clt = typename coarsener_t::coarse_level_triple;
  handle.h  = coarsener_t::HECv1;
  handle.b  = coarsener_t::Hybrid;
  coarsener_t::generate_coarse_graphs(handle, A, true);
  std::list<clt> levels = handle.results;
  auto fine             = levels.begin();
  auto coarse           = fine;
  coarse++;
  while (coarse != levels.end()) {
    bool correct_aggregator = verify_aggregator(fine->mtx, coarse->interp_mtx);
    EXPECT_TRUE(correct_aggregator) << "Multilevel coarsening produced invalid aggregator on level "
                                    << coarse->level - 1;
    bool correct_graph      = verify_is_graph<crsMat>(coarse->mtx);
    bool correct_coarsening = verify_coarsening<coarsener_t>(*fine, *coarse);
    EXPECT_TRUE(correct_graph) << "Multilevel coarsening produced invalid graph on level " << coarse->level;
    EXPECT_TRUE(correct_coarsening) << "Multilevel coarsening produced invalid coarsening on level " << coarse->level;
    fine++;
    coarse++;
  }
}

template <typename scalar, typename lno_t, typename size_type, typename device>
void test_coarsen_grid() {
  using crsMat      = KokkosSparse::CrsMatrix<scalar, lno_t, device, void, size_type>;
  using graph_type  = typename crsMat::StaticCrsGraphType;
  using c_entries_t = typename graph_type::entries_type;
  using entries_t   = typename c_entries_t::non_const_type;
  crsMat A          = gen_grid<crsMat>();
  entries_t vWgts("vertex weights", A.numRows());
  Kokkos::deep_copy(vWgts, static_cast<typename entries_t::value_type>(1));
  using coarsener_t = coarse_builder<crsMat>;
  typename coarsener_t::coarsen_handle handle;
  using clt = typename coarsener_t::coarse_level_triple;
  clt fine_A;
  fine_A.mtx                                              = A;
  fine_A.vtx_wgts                                         = vWgts;
  fine_A.level                                            = 0;
  fine_A.uniform_weights                                  = true;
  std::vector<typename coarsener_t::Heuristic> heuristics = {coarsener_t::HECv1,   coarsener_t::Match,
                                                             coarsener_t::MtMetis, coarsener_t::MIS2,
                                                             coarsener_t::GOSHv1,  coarsener_t::GOSHv2};
  std::vector<typename coarsener_t::Builder> builders = {coarsener_t::Sort, coarsener_t::Hashmap, coarsener_t::Hybrid,
                                                         coarsener_t::Spgemm, coarsener_t::Spgemm_transpose_first};
  for (auto h : heuristics) {
    handle.h                = h;
    crsMat aggregator       = coarsener_t::generate_coarse_mapping(handle, fine_A.mtx, true);
    bool correct_aggregator = verify_aggregator(fine_A.mtx, aggregator);
    EXPECT_TRUE(correct_aggregator) << "Aggregation heuristic " << static_cast<int>(h)
                                    << " produced invalid aggregator.";
    for (auto b : builders) {
      handle.b                = b;
      clt coarse_A            = coarsener_t::build_coarse_graph(handle, fine_A, aggregator);
      bool correct_graph      = verify_is_graph<crsMat>(coarse_A.mtx);
      bool correct_coarsening = verify_coarsening<coarsener_t>(fine_A, coarse_A);
      EXPECT_TRUE(correct_graph) << "Coarsening with dedupe method " << static_cast<int>(b)
                                 << " produced invalid graph with aggregation heuristic " << static_cast<int>(h) << ".";
      EXPECT_TRUE(correct_coarsening) << "Coarsening with dedupe method " << static_cast<int>(b)
                                      << " produced invalid coarsening with aggregation heuristic "
                                      << static_cast<int>(h) << ".";
    }
  }
}

template <typename scalar, typename lno_t, typename size_type, typename device>
void test_coarsen_random(lno_t numVerts, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using execution_space = typename device::execution_space;
  using crsMat          = KokkosSparse::CrsMatrix<scalar, lno_t, device, void, size_type>;
  using graph_type      = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t      = typename graph_type::row_map_type;
  using c_entries_t     = typename graph_type::entries_type;
  using rowmap_t        = typename c_rowmap_t::non_const_type;
  using entries_t       = typename c_entries_t::non_const_type;
  using svt             = typename crsMat::values_type;
  // Generate graph
  crsMat A =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
  auto G = A.graph;
  // Symmetrize the graph
  rowmap_t symRowmap;
  entries_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<c_rowmap_t, c_entries_t, rowmap_t, entries_t, execution_space>(
      numVerts, G.row_map, G.entries, symRowmap, symEntries);
  graph_type GS(symEntries, symRowmap);
  svt symValues("sym values", symEntries.extent(0));
  entries_t vWgts("vertex weights", numVerts);
  Kokkos::deep_copy(symValues, static_cast<scalar>(2.5));
  Kokkos::deep_copy(vWgts, static_cast<typename entries_t::value_type>(1));
  crsMat AS("A symmetric", numVerts, symValues, GS);
  using coarsener_t = coarse_builder<crsMat>;
  typename coarsener_t::coarsen_handle handle;
  using clt = typename coarsener_t::coarse_level_triple;
  clt fine_A;
  fine_A.mtx                                              = AS;
  fine_A.vtx_wgts                                         = vWgts;
  fine_A.level                                            = 0;
  fine_A.uniform_weights                                  = true;
  std::vector<typename coarsener_t::Heuristic> heuristics = {coarsener_t::HECv1,   coarsener_t::Match,
                                                             coarsener_t::MtMetis, coarsener_t::MIS2,
                                                             coarsener_t::GOSHv1,  coarsener_t::GOSHv2};
  std::vector<typename coarsener_t::Builder> builders = {coarsener_t::Sort, coarsener_t::Hashmap, coarsener_t::Hybrid,
                                                         coarsener_t::Spgemm, coarsener_t::Spgemm_transpose_first};
  for (auto h : heuristics) {
    handle.h                = h;
    crsMat aggregator       = coarsener_t::generate_coarse_mapping(handle, fine_A.mtx, true);
    bool correct_aggregator = verify_aggregator(fine_A.mtx, aggregator);
    EXPECT_TRUE(correct_aggregator) << "Aggregation heuristic " << static_cast<int>(h)
                                    << " produced invalid aggregator.";
    for (auto b : builders) {
      handle.b                = b;
      clt coarse_A            = coarsener_t::build_coarse_graph(handle, fine_A, aggregator);
      bool correct_graph      = verify_is_graph<crsMat>(coarse_A.mtx);
      bool correct_coarsening = verify_coarsening<coarsener_t>(fine_A, coarse_A);
      EXPECT_TRUE(correct_graph) << "Coarsening with dedupe method " << static_cast<int>(b)
                                 << " produced invalid graph with aggregation heuristic " << static_cast<int>(h) << ".";
      EXPECT_TRUE(correct_coarsening) << "Coarsening with dedupe method " << static_cast<int>(b)
                                      << " produced invalid coarsening with aggregation heuristic "
                                      << static_cast<int>(h) << ".";
    }
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                          \
  TEST_F(TestCategory, graph##_##random_graph_coarsen##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {          \
    test_coarsen_random<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 5000 * 20, 1000, 10);                           \
    test_coarsen_random<SCALAR, ORDINAL, OFFSET, DEVICE>(50, 50 * 10, 40, 10);                                 \
    test_coarsen_random<SCALAR, ORDINAL, OFFSET, DEVICE>(5, 5 * 3, 5, 0);                                      \
  }                                                                                                            \
  TEST_F(TestCategory, graph##_##grid_graph_coarsen##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {            \
    test_coarsen_grid<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                      \
  }                                                                                                            \
  TEST_F(TestCategory, graph##_##grid_graph_multilevel_coarsen##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_multilevel_coarsen_grid<SCALAR, ORDINAL, OFFSET, DEVICE>();                                           \
  }

// FIXME_SYCL
#ifndef KOKKOS_ENABLE_SYCL
#if defined(KOKKOSKERNELS_INST_DOUBLE)
#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestDevice)
#endif
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
