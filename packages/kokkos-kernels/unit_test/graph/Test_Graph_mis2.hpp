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
// Questions? Contact Brian Kelley (bmkelle@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>
#include <random>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_MIS2.hpp"
#include "KokkosGraph_ExplicitCoarsening.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

namespace Test {

template<typename lno_t, typename size_type, typename rowmap_t, typename entries_t, typename mis_t>
bool verifyD2MIS(
    lno_t numVerts,
    const rowmap_t& rowmap, const entries_t& entries,
    const mis_t& misArray)
{
  //set a std::set of the mis, for fast membership test
  std::set<lno_t> mis;
  for(size_t i = 0; i < misArray.extent(0); i++)
    mis.insert(misArray(i));
  for(lno_t i = 0; i < numVerts; i++)
  {
    //determine whether another vertex in the set is
    //within 2 hops of i.
    bool misIn2Hops = false;
    for(size_type j = rowmap(i); j < rowmap(i + 1); j++)
    {
      lno_t nei1 = entries(j);
      if(nei1 == i || nei1 >= numVerts)
        continue;
      if(mis.find(nei1) != mis.end())
      {
        misIn2Hops = true;
        break;
      }
      for(size_type k = rowmap(nei1); k < rowmap(nei1 + 1); k++)
      {
        lno_t nei2 = entries(k);
        if(nei2 == i || nei2 >= numVerts)
          continue;
        if(mis.find(nei2) != mis.end())
        {
          misIn2Hops = true;
          break;
        }
      }
    }
    if(mis.find(i) == mis.end())
    {
      //i is not in the set
      if(!misIn2Hops)
      {
        std::cout << "INVALID D2 MIS: vertex " << i << " is not in the set,\n";
        std::cout << "but there are no vertices in the set within 2 hops.\n";
        return false;
      }
    }
    else
    {
      //i is in the set
      if(misIn2Hops)
      {
        std::cout << "INVALID D2 MIS: vertex " << i << " is in the set,\n";
        std::cout << "but there is another vertex within 2 hops which is also in the set.\n";
        return false;
      }
    }
  }
  return true;
}
}

template<typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_mis2(lno_t numVerts, size_type nnz, lno_t bandwidth, lno_t row_size_variance)
{
  using execution_space = typename device::execution_space;
  using crsMat = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t = typename c_rowmap_t::non_const_type;
  using entries_t = typename c_entries_t::non_const_type;
  //Generate graph, and add some out-of-bounds columns
  crsMat A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
  auto G = A.graph;
  //Symmetrize the graph
  rowmap_t symRowmap;
  entries_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap
    <c_rowmap_t, c_entries_t,
    rowmap_t, entries_t, execution_space>
      (numVerts, G.row_map, G.entries, symRowmap, symEntries);
  auto rowmapHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symRowmap);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symEntries);
  //For each algorithm, compute and verify the MIS
  std::vector<MIS2_Algorithm> algos
    = {MIS2_FAST, MIS2_QUALITY};
  for(auto algo : algos)
  {
    auto mis = graph_d2_mis<device, rowmap_t, entries_t>(symRowmap, symEntries, algo);
    auto misHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mis);
    bool success = Test::verifyD2MIS
      <lno_t, size_type, decltype(rowmapHost), decltype(entriesHost), decltype(misHost)>
      (numVerts, rowmapHost, entriesHost, misHost);
    EXPECT_TRUE(success) << "Dist-2 MIS (algo " << (int) algo << ") produced invalid set.";
  }
}

template<typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_mis2_coarsening(lno_t numVerts, size_type nnz, lno_t bandwidth, lno_t row_size_variance)
{
  using execution_space = typename device::execution_space;
  using crsMat = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t = typename c_rowmap_t::non_const_type;
  using entries_t = typename c_entries_t::non_const_type;
  //Generate graph, and add some out-of-bounds columns
  crsMat A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
  auto G = A.graph;
  //Symmetrize the graph
  rowmap_t symRowmap;
  entries_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap
    <c_rowmap_t, c_entries_t,
    rowmap_t, entries_t, execution_space>
      (numVerts, G.row_map, G.entries, symRowmap, symEntries);
  auto rowmapHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symRowmap);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), symEntries);
  //For each algorithm, compute and verify the MIS
  std::vector<MIS2_Algorithm> algos
    = {MIS2_FAST, MIS2_QUALITY};
  for(auto algo : algos)
  {
    lno_t numClusters = 0;
    auto labels = graph_mis2_coarsen<device, rowmap_t, entries_t>(symRowmap, symEntries, numClusters, algo);
    auto labelsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), labels);
    //Not a strong test, but sanity check the number of clusters returned
    EXPECT_TRUE(numClusters >= 1 && numClusters <= numVerts);
    //Check that every label is in the range [0, numClusters)
    for(lno_t i = 0; i < numVerts; i++)
      EXPECT_TRUE(0 <= labelsHost(i) && labelsHost(i) < numClusters);
    //Test explicit coarsening given the labels, with and without compressing the result
    rowmap_t coarseRowmapNC, coarseRowmapC;
    entries_t coarseEntriesNC, coarseEntriesC;
    KokkosGraph::Experimental::graph_explicit_coarsen<device, rowmap_t, entries_t, entries_t, rowmap_t, entries_t>
      (symRowmap, symEntries, labels, numClusters, coarseRowmapNC, coarseEntriesNC, false);
    KokkosGraph::Experimental::graph_explicit_coarsen<device, rowmap_t, entries_t, entries_t, rowmap_t, entries_t>
      (symRowmap, symEntries, labels, numClusters, coarseRowmapC, coarseEntriesC, true);
    EXPECT_EQ(coarseRowmapC.extent(0), numClusters + 1);
    EXPECT_EQ(coarseRowmapNC.extent(0), numClusters + 1);
    //Check that coarse graph doesn't have more edges than fine graph
    EXPECT_LE(coarseEntriesC.extent(0), symEntries.extent(0));
    EXPECT_LE(coarseEntriesNC.extent(0), symEntries.extent(0));
    //Verify compression is working.
    auto hostRowmapNC = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarseRowmapNC);
    auto hostEntriesNC = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarseEntriesNC);
    auto hostRowmapC = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarseRowmapC);
    auto hostEntriesC = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarseEntriesC);
    for(lno_t i = 0; i < numClusters; i++)
    {
      //std::set maintains uniqueness as well as ascending order of elements.
      //So it should exactly match the entries in the compressed version.
      std::set<lno_t> uniqueEntries;
      for(size_type j = hostRowmapNC(i); j < hostRowmapNC(i + 1); j++)
      {
        uniqueEntries.insert(hostEntriesNC(j));
      }
      size_type compressedRowLen = hostRowmapC(i + 1) - hostRowmapC(i);
      ASSERT_EQ(uniqueEntries.size(), compressedRowLen);
      auto it = uniqueEntries.begin();
      for(size_type j = hostRowmapC(i); j < hostRowmapC(i + 1); j++)
      {
        EXPECT_EQ(*it, hostEntriesC(j));
        it++;
      }
    }
  }
}

template<typename scalar_unused, typename lno_t, typename size_type, typename device>
void test_mis2_coarsening_zero_rows()
{
  using crsMat = KokkosSparse::CrsMatrix<double, lno_t, device, void, size_type>;
  using graph_type = typename crsMat::StaticCrsGraphType;
  using c_rowmap_t = typename graph_type::row_map_type;
  using c_entries_t = typename graph_type::entries_type;
  using rowmap_t = typename c_rowmap_t::non_const_type;
  using entries_t = typename c_entries_t::non_const_type;
  rowmap_t fineRowmap;
  entries_t fineEntries;
  //note: MIS2 coarsening first calls MIS2 on the fine graph, so this covers the zero-row case for MIS2 alone.
  lno_t numClusters;
  auto labels = graph_mis2_coarsen<device, rowmap_t, entries_t>(fineRowmap, fineEntries, numClusters, KokkosGraph::MIS2_FAST);
  EXPECT_EQ(numClusters, 0);
  EXPECT_EQ(labels.extent(0), 0);
  //coarsen, should also produce a graph with 0 rows/entries
  rowmap_t coarseRowmap;
  entries_t coarseEntries;
  KokkosGraph::Experimental::graph_explicit_coarsen<device, rowmap_t, entries_t, entries_t, rowmap_t, entries_t>
    (fineRowmap, fineEntries, labels, 0, coarseRowmap, coarseEntries, false);
  EXPECT_LE(coarseRowmap.extent(0), 1);
  EXPECT_EQ(coarseEntries.extent(0), 0);
  KokkosGraph::Experimental::graph_explicit_coarsen<device, rowmap_t, entries_t, entries_t, rowmap_t, entries_t>
    (fineRowmap, fineEntries, labels, 0, coarseRowmap, coarseEntries, true);
  EXPECT_LE(coarseRowmap.extent(0), 1);
  EXPECT_EQ(coarseEntries.extent(0), 0);
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
    TEST_F(TestCategory, graph##_##graph_mis2##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) \
    { \
      test_mis2<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 5000 * 20, 1000, 10); \
      test_mis2<SCALAR, ORDINAL, OFFSET, DEVICE>(50, 50 * 10, 40, 10); \
      test_mis2<SCALAR, ORDINAL, OFFSET, DEVICE>(5, 5 * 3, 5, 0); \
    } \
    TEST_F(TestCategory, graph##_##graph_mis2_coarsening##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) \
    { \
      test_mis2_coarsening<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 5000 * 200, 2000, 10); \
      test_mis2_coarsening<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 5000 * 20, 1000, 10); \
      test_mis2_coarsening<SCALAR, ORDINAL, OFFSET, DEVICE>(50, 50 * 10, 40, 10); \
      test_mis2_coarsening<SCALAR, ORDINAL, OFFSET, DEVICE>(5, 5 * 3, 5, 0); \
      test_mis2_coarsening_zero_rows<SCALAR, ORDINAL, OFFSET, DEVICE>(); \
    }

#if defined(KOKKOSKERNELS_INST_DOUBLE)
#if(defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) \
  || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if(defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) \
  || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if(defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) \
  || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if(defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) \
  || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif
#endif
