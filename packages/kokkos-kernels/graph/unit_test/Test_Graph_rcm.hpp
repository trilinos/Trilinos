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

#include "KokkosGraph_RCM.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

#include <vector>

// Generates a graph from 3D 7-pt stencil. Slices grid into 2 connected
// components near the middle of X dimension.
template <typename rowmap_t, typename entries_t>
void generate7pt(rowmap_t& rowmapView, entries_t& entriesView, int gridX, int gridY, int gridZ) {
  using size_type   = typename rowmap_t::non_const_value_type;
  using lno_t       = typename entries_t::non_const_value_type;
  auto getVertexID  = [=](lno_t x, lno_t y, lno_t z) -> lno_t { return x + y * gridX + z * gridX * gridY; };
  lno_t numVertices = gridX * gridY * gridZ;
  // Generate the graph on host (use std::vector to not need to know
  // how many entries ahead of time)
  std::vector<size_type> rowmap(numVertices + 1);
  std::vector<lno_t> entries;
  rowmap[0]    = 0;
  lno_t xslice = gridX / 2;
  for (lno_t k = 0; k < gridZ; k++) {
    for (lno_t j = 0; j < gridY; j++) {
      for (lno_t i = 0; i < gridX; i++) {
        lno_t v = getVertexID(i, j, k);
        if (i != 0 && i != xslice + 1) entries.push_back(getVertexID(i - 1, j, k));
        if (i != gridX - 1 && i != xslice) entries.push_back(getVertexID(i + 1, j, k));
        if (j != 0) entries.push_back(getVertexID(i, j - 1, k));
        if (j != gridY - 1) entries.push_back(getVertexID(i, j + 1, k));
        if (k != 0) entries.push_back(getVertexID(i, j, k - 1));
        if (k != gridZ - 1) entries.push_back(getVertexID(i, j, k + 1));
        rowmap[v + 1] = entries.size();
      }
    }
  }
  size_type numEdges = entries.size();
  // Now that the graph is formed, copy rowmap and entries to Kokkos::Views in
  // device memory The nonowning host views just alias the std::vectors.
  Kokkos::View<size_type*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> rowmapHost(rowmap.data(),
                                                                                                  numVertices + 1);
  Kokkos::View<lno_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> entriesHost(entries.data(),
                                                                                               numEdges);
  // Allocate owning views on device with the correct size.
  rowmapView  = rowmap_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Rowmap"), numVertices + 1);
  entriesView = entries_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Colinds"), numEdges);
  // Copy the graph from host to device
  Kokkos::deep_copy(rowmapView, rowmapHost);
  Kokkos::deep_copy(entriesView, entriesHost);
}

template <typename rowmap_t, typename entries_t, typename labels_t>
int maxBandwidth(const rowmap_t& rowmap, const entries_t& entries, const labels_t& invPerm, const labels_t& perm) {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  lno_t numVerts  = std::max(1, rowmap.extent_int(0)) - 1;
  int bw          = 0;
  for (lno_t i = 0; i < numVerts; i++) {
    lno_t origRow = perm(i);
    for (size_type j = rowmap(origRow); j < rowmap(origRow + 1); j++) {
      lno_t origNei = entries(j);
      lno_t nei     = invPerm(origNei);
      if (nei > i) {
        lno_t thisBW = nei - i;
        if (thisBW > bw) bw = thisBW;
      }
    }
  }
  return bw;
}

template <typename device, typename rowmap_t, typename entries_t>
void test_rcm(const rowmap_t& rowmap, const entries_t& entries, bool expectBandwidthReduced) {
  using lno_t      = typename entries_t::non_const_value_type;
  auto rcm         = KokkosGraph::Experimental::graph_rcm<device, rowmap_t, entries_t>(rowmap, entries);
  auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmap);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries);
  auto rcmHost     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rcm);
  lno_t numVerts   = std::max(rowmap.extent_int(0), 1) - 1;
  decltype(rcmHost) rcmPermHost(Kokkos::view_alloc(Kokkos::WithoutInitializing, "RCMPerm"), numVerts);
  for (lno_t i = 0; i < numVerts; i++) rcmPermHost(rcmHost(i)) = i;
  // make sure each row index shows up exactly once
  {
    std::vector<int> counts(numVerts);
    for (lno_t i = 0; i < numVerts; i++) {
      lno_t orig = rcmHost(i);
      ASSERT_GE(orig, 0);
      ASSERT_LT(orig, numVerts);
      counts[orig]++;
    }
    for (lno_t i = 0; i < numVerts; i++) ASSERT_EQ(counts[i], 1);
  }
  if (expectBandwidthReduced) {
    Kokkos::View<lno_t*, Kokkos::HostSpace> identityOrder(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Identity"),
                                                          numVerts);
    for (lno_t i = 0; i < numVerts; i++) identityOrder(i) = i;
    size_t origBW = maxBandwidth(rowmapHost, entriesHost, identityOrder, identityOrder);
    size_t rcmBW  = maxBandwidth(rowmapHost, entriesHost, rcmHost, rcmPermHost);
    EXPECT_LE(rcmBW, origBW);
  }
}

template <typename lno_t, typename size_type, typename device>
void test_rcm_zerorows() {
  using graph_t   = Kokkos::StaticCrsGraph<lno_t, default_layout, device, void, size_type>;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  rowmap_t rowmap;
  entries_t entries;
  test_rcm<device>(rowmap, entries, false);
}

template <typename lno_t, typename size_type, typename device>
void test_rcm_7pt(lno_t gridX, lno_t gridY, lno_t gridZ, bool expectBandwidthReduced) {
  using graph_t   = Kokkos::StaticCrsGraph<lno_t, default_layout, device, void, size_type>;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  rowmap_t rowmap;
  entries_t entries;
  generate7pt(rowmap, entries, gridX, gridY, gridZ);
  test_rcm<device>(rowmap, entries, expectBandwidthReduced);
}

template <typename lno_t, typename size_type, typename device>
void test_rcm_4clique() {
  using graph_t   = Kokkos::StaticCrsGraph<lno_t, default_layout, device, void, size_type>;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  rowmap_t rowmap("rowmap", 5);
  entries_t entries("entries", 16);
  auto rowmap_host  = Kokkos::create_mirror_view(rowmap);
  auto entries_host = Kokkos::create_mirror_view(entries);
  for (lno_t i = 0; i < 5; i++) rowmap_host(i) = i * 4;
  for (lno_t i = 0; i < 16; i++) entries_host(i) = i % 4;
  Kokkos::deep_copy(rowmap, rowmap_host);
  Kokkos::deep_copy(entries, entries_host);
  test_rcm<device>(rowmap, entries, false);
}

template <typename lno_t, typename size_type, typename device>
void test_rcm_multiple_components() {
  using graph_t   = Kokkos::StaticCrsGraph<lno_t, default_layout, device, void, size_type>;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  // Generate a single 3D grid first
  rowmap_t rowmap_cube;
  entries_t entries_cube;
  generate7pt(rowmap_cube, entries_cube, 7, 7, 7);
  auto rowmap_cube_host  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmap_cube);
  auto entries_cube_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries_cube);
  lno_t nv_cube          = 7 * 7 * 7;
  lno_t ne_cube          = entries_cube.extent(0);
  // Now replicate the graph twice, so there are 2 disconnected copies of the
  // cube
  rowmap_t rowmap("rowmap", nv_cube * 2 + 1);
  entries_t entries("entries", ne_cube * 2);
  auto rowmap_host  = Kokkos::create_mirror_view(rowmap);
  auto entries_host = Kokkos::create_mirror_view(entries);
  for (lno_t i = 0; i <= nv_cube * 2; i++) {
    if (i < nv_cube)
      rowmap_host(i) = rowmap_cube_host(i);
    else
      rowmap_host(i) = ne_cube + rowmap_cube_host(i - nv_cube);
  }
  for (lno_t i = 0; i < ne_cube * 2; i++) {
    if (i < ne_cube)
      entries_host(i) = entries_cube_host(i);
    else
      entries_host(i) = nv_cube + entries_cube_host(i - ne_cube);
  }
  Kokkos::deep_copy(rowmap, rowmap_host);
  Kokkos::deep_copy(entries, entries_host);
  test_rcm<device>(rowmap, entries, true);
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                    \
  TEST_F(TestCategory, graph##_##rcm_zerorows##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {            \
    test_rcm_zerorows<ORDINAL, OFFSET, DEVICE>();                                                        \
  }                                                                                                      \
  TEST_F(TestCategory, graph##_##rcm_7pt##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {                 \
    test_rcm_7pt<ORDINAL, OFFSET, DEVICE>(1, 1, 1, false);                                               \
    test_rcm_7pt<ORDINAL, OFFSET, DEVICE>(2, 1, 1, false);                                               \
    test_rcm_7pt<ORDINAL, OFFSET, DEVICE>(6, 3, 3, true);                                                \
    test_rcm_7pt<ORDINAL, OFFSET, DEVICE>(20, 20, 20, true);                                             \
    test_rcm_7pt<ORDINAL, OFFSET, DEVICE>(100, 100, 1, true);                                            \
  }                                                                                                      \
  TEST_F(TestCategory, graph##_##rcm_4clique##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {             \
    test_rcm_4clique<ORDINAL, OFFSET, DEVICE>();                                                         \
  }                                                                                                      \
  TEST_F(TestCategory, graph##_##rcm_multiple_components##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_rcm_multiple_components<ORDINAL, OFFSET, DEVICE>();                                             \
  }

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

#undef EXECUTE_TEST
