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

#include "KokkosGraph_RCM.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include <vector>

// Generates a graph from 3D 7-pt stencil. Slices grid into 2 connected
// components near the middle of X dimension.
template <typename rowmap_t, typename entries_t>
void generate7pt(rowmap_t& rowmapView, entries_t& entriesView, int gridX,
                 int gridY, int gridZ) {
  using size_type  = typename rowmap_t::non_const_value_type;
  using lno_t      = typename entries_t::non_const_value_type;
  auto getVertexID = [=](lno_t x, lno_t y, lno_t z) -> lno_t {
    return x + y * gridX + z * gridX * gridY;
  };
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
        if (i != 0 && i != xslice + 1)
          entries.push_back(getVertexID(i - 1, j, k));
        if (i != gridX - 1 && i != xslice)
          entries.push_back(getVertexID(i + 1, j, k));
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
  Kokkos::View<size_type*, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      rowmapHost(rowmap.data(), numVertices + 1);
  Kokkos::View<lno_t*, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      entriesHost(entries.data(), numEdges);
  // Allocate owning views on device with the correct size.
  rowmapView =
      rowmap_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Rowmap"),
               numVertices + 1);
  entriesView = entries_t(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Colinds"), numEdges);
  // Copy the graph from host to device
  Kokkos::deep_copy(rowmapView, rowmapHost);
  Kokkos::deep_copy(entriesView, entriesHost);
}

template <typename rowmap_t, typename entries_t, typename labels_t>
int maxBandwidth(const rowmap_t& rowmap, const entries_t& entries,
                 const labels_t& invPerm, const labels_t& perm) {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  lno_t numVerts  = rowmap.extent(0) - 1;
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

template <typename lno_t, typename size_type, typename device>
void test_rcm(lno_t gridX, lno_t gridY, lno_t gridZ) {
  typedef
      typename KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>
          crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type rowmap_t;
  typedef typename graph_t::entries_type entries_t;
  lno_t numVerts = gridX * gridY * gridZ;
  typename rowmap_t::non_const_type rowmap;
  typename entries_t::non_const_type entries;
  generate7pt(rowmap, entries, gridX, gridY, gridZ);
  auto rcm = KokkosGraph::Experimental::graph_rcm<device, rowmap_t, entries_t>(
      rowmap, entries);
  auto rowmapHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmap);
  auto entriesHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries);
  auto rcmHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rcm);
  decltype(rcmHost) rcmPermHost(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "RCMPerm"), numVerts);
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
  Kokkos::View<lno_t*, Kokkos::HostSpace> identityOrder(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Identity"), numVerts);
  for (lno_t i = 0; i < numVerts; i++) identityOrder(i) = i;
  size_t origBW =
      maxBandwidth(rowmapHost, entriesHost, identityOrder, identityOrder);
  size_t rcmBW = maxBandwidth(rowmapHost, entriesHost, rcmHost, rcmPermHost);
  EXPECT_LE(rcmBW, origBW);
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                  \
  TEST_F(TestCategory,                                                 \
         graph##_##rcm##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_rcm<ORDINAL, OFFSET, DEVICE>(6, 3, 3);                        \
    test_rcm<ORDINAL, OFFSET, DEVICE>(20, 20, 20);                     \
    test_rcm<ORDINAL, OFFSET, DEVICE>(100, 100, 1);                    \
  }

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#undef EXECUTE_TEST
