// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <Kokkos_Sort.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <tuple>
#include <vector>
#include "stk_search/SearchMethod.hpp"

namespace std {
template <typename Ident, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::IdentProc<Ident,Proc>,stk::search::IdentProc<Ident,Proc> > const& ip)
{
  return out << "[" << ip.first << ":" << ip.second << "]";
}
} // namespace std

namespace {

template <typename SearchResultsView>
void expect_search_results_with_views(int num_procs, int proc_id, SearchResultsView const& searchResults)
{
  if (num_procs == 1) {
    ASSERT_EQ(searchResults.extent(0), 2u);
    EXPECT_EQ(searchResults[0], (IdentProcIntersection{IdentProc(0, 0), IdentProc(2, 0)}) );
    EXPECT_EQ(searchResults[1], (IdentProcIntersection{IdentProc(1, 0), IdentProc(3, 0)}) );
  } else {
    if (proc_id == 0) {
      ASSERT_EQ(searchResults.size(), 4u);
      EXPECT_EQ(searchResults[0], (IdentProcIntersection{IdentProc(0, 0), IdentProc(2, 0)}) );
      EXPECT_EQ(searchResults[1], (IdentProcIntersection{IdentProc(1, 0), IdentProc(3, 0)}) );
      EXPECT_EQ(searchResults[2], (IdentProcIntersection{IdentProc(4, 1), IdentProc(2, 0)}) );
      EXPECT_EQ(searchResults[3], (IdentProcIntersection{IdentProc(5, 1), IdentProc(3, 0)}) );
    } else if (proc_id == num_procs - 1) {
      ASSERT_EQ(searchResults.size(), 4u);
      int prev = proc_id - 1;
      EXPECT_EQ(searchResults[0], (IdentProcIntersection{IdentProc(proc_id * 4, proc_id), IdentProc(prev * 4 + 2, prev)}) );
      EXPECT_EQ(searchResults[1], (IdentProcIntersection{IdentProc(proc_id * 4, proc_id), IdentProc(proc_id * 4 + 2, proc_id)}) );
      EXPECT_EQ(searchResults[2], (IdentProcIntersection{IdentProc(proc_id * 4 + 1, proc_id), IdentProc(prev * 4 + 3, prev)}) );
      EXPECT_EQ(searchResults[3], (IdentProcIntersection{IdentProc(proc_id * 4 + 1, proc_id), IdentProc(proc_id * 4 + 3, proc_id)}) );
    } else {
      ASSERT_EQ(searchResults.size(), 6u);
      int prev = proc_id - 1;
      int next = proc_id + 1;
      EXPECT_EQ(searchResults[0], (IdentProcIntersection{IdentProc(proc_id * 4, proc_id), IdentProc(prev * 4 + 2, prev)}) );
      EXPECT_EQ(searchResults[1], (IdentProcIntersection{IdentProc(proc_id * 4, proc_id), IdentProc(proc_id * 4 + 2, proc_id)}) );
      EXPECT_EQ(searchResults[2], (IdentProcIntersection{IdentProc(proc_id * 4 + 1, proc_id), IdentProc(prev * 4 + 3, prev)}) );
      EXPECT_EQ(searchResults[3], (IdentProcIntersection{IdentProc(proc_id * 4 + 1, proc_id), IdentProc(proc_id * 4 + 3, proc_id)}) );
      EXPECT_EQ(searchResults[4], (IdentProcIntersection{IdentProc(next * 4, next), IdentProc(proc_id * 4 + 2, proc_id)}) );
      EXPECT_EQ(searchResults[5], (IdentProcIntersection{IdentProc(next * 4 + 1, next), IdentProc(proc_id * 4 + 3, proc_id)}) );
    }
  }
}

void expect_search_results(int num_procs, int proc_id, const SearchResults&  searchResults)
{
  if (num_procs == 1) {
    ASSERT_EQ(searchResults.size(), 2u);
    EXPECT_EQ(searchResults[0], std::make_pair(IdentProc(0, 0), IdentProc(2, 0)));
    EXPECT_EQ(searchResults[1], std::make_pair(IdentProc(1, 0), IdentProc(3, 0)));
  } else {
    if (proc_id == 0) {
      ASSERT_EQ(searchResults.size(), 4u);
      EXPECT_EQ(searchResults[0], std::make_pair(IdentProc(0, 0), IdentProc(2, 0)));
      EXPECT_EQ(searchResults[1], std::make_pair(IdentProc(1, 0), IdentProc(3, 0)));
      EXPECT_EQ(searchResults[2], std::make_pair(IdentProc(4, 1), IdentProc(2, 0)));
      EXPECT_EQ(searchResults[3], std::make_pair(IdentProc(5, 1), IdentProc(3, 0)));
    } else if (proc_id == num_procs - 1) {
      ASSERT_EQ(searchResults.size(), 4u);
      int prev = proc_id - 1;
      EXPECT_EQ(searchResults[0], std::make_pair(IdentProc(proc_id * 4, proc_id), IdentProc(prev * 4 + 2, prev)));
      EXPECT_EQ(searchResults[1], std::make_pair(IdentProc(proc_id * 4, proc_id), IdentProc(proc_id * 4 + 2, proc_id)));
      EXPECT_EQ(searchResults[2], std::make_pair(IdentProc(proc_id * 4 + 1, proc_id), IdentProc(prev * 4 + 3, prev)));
      EXPECT_EQ(searchResults[3], std::make_pair(IdentProc(proc_id * 4 + 1, proc_id), IdentProc(proc_id * 4 + 3, proc_id)));
    } else {
      ASSERT_EQ(searchResults.size(), 6u);
      int prev = proc_id - 1;
      int next = proc_id + 1;
      EXPECT_EQ(searchResults[0], std::make_pair(IdentProc(proc_id * 4, proc_id), IdentProc(prev * 4 + 2, prev)));
      EXPECT_EQ(searchResults[1], std::make_pair(IdentProc(proc_id * 4, proc_id), IdentProc(proc_id * 4 + 2, proc_id)));
      EXPECT_EQ(searchResults[2], std::make_pair(IdentProc(proc_id * 4 + 1, proc_id), IdentProc(prev * 4 + 3, prev)));
      EXPECT_EQ(searchResults[3], std::make_pair(IdentProc(proc_id * 4 + 1, proc_id), IdentProc(proc_id * 4 + 3, proc_id)));
      EXPECT_EQ(searchResults[4], std::make_pair(IdentProc(next * 4, next), IdentProc(proc_id * 4 + 2, proc_id)));
      EXPECT_EQ(searchResults[5], std::make_pair(IdentProc(next * 4 + 1, next), IdentProc(proc_id * 4 + 3, proc_id)));
    }
  }
}

template <typename FloatType, typename ExecSpace=Kokkos::DefaultExecutionSpace>
void test_coarse_search_for_algorithm_with_views(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  using HostSpace = Kokkos::DefaultHostExecutionSpace;
  using BoxType = stk::search::Box<FloatType>;
  using PointType = stk::search::Point<FloatType>;
  using BoxIdentProcType = stk::search::BoxIdentProc<BoxType, IdentProc>;
  using BoxIdentProcViewType = Kokkos::View<BoxIdentProcType*, ExecSpace>;
  using SearchResultsViewType = Kokkos::View<stk::search::IdentProcIntersection<IdentProc, IdentProc>*, ExecSpace>;

  BoxType box;
  IdentProc identProc;
  auto domain = BoxIdentProcViewType("domain test view", 2);
  auto range = BoxIdentProcViewType("range test view", 2);

  auto domainHost = Kokkos::create_mirror_view_and_copy(HostSpace{}, domain);
  auto rangeHost = Kokkos::create_mirror_view_and_copy(HostSpace{}, range);

  box = BoxType( PointType(proc_id + 0.1, 0.0, 0.0), PointType(proc_id + 0.9, 1.0, 1.0));
  identProc = IdentProc(proc_id * 4, proc_id);
  domainHost(0) = BoxIdentProcType{box, identProc};

  box = BoxType( PointType(proc_id + 0.1, 2.0, 0.0), PointType(proc_id + 0.9, 3.0, 1.0));
  identProc = IdentProc(proc_id * 4+1, proc_id);
  domainHost(1) = BoxIdentProcType{box, identProc};

  box = BoxType( PointType(proc_id + 0.6, 0.5, 0.0), PointType(proc_id + 1.4, 1.5, 1.0));
  identProc = IdentProc(proc_id * 4+2, proc_id);
  rangeHost(0) = BoxIdentProcType{box, identProc};

  box = BoxType( PointType(proc_id + 0.6, 2.5, 0.0), PointType(proc_id + 1.4, 3.5, 1.0));
  identProc = IdentProc(proc_id * 4+3, proc_id);
  rangeHost(1) = BoxIdentProcType{box, identProc};

  Kokkos::deep_copy(domain, domainHost);
  Kokkos::deep_copy(range, rangeHost);

  SearchResultsViewType searchResults;

  auto execSpace = ExecSpace{};
  bool enforceSearchResultSymmetry = true;
  bool autoSwapDomainAndRange = true;
  bool sortSearchResults = true;

  stk::search::coarse_search(domain, range, algorithm, comm, searchResults, execSpace,
                                     enforceSearchResultSymmetry, autoSwapDomainAndRange, sortSearchResults);

  auto searchResultsHost = Kokkos::create_mirror_view_and_copy(HostSpace{}, searchResults);

  expect_search_results_with_views(num_procs, proc_id, searchResultsHost);
}

template <typename FloatType>
void test_coarse_search_for_algorithm(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  using BoxType = stk::search::Box<FloatType>;
  using PointType = stk::search::Point<FloatType>;
  using BoxIdentProcVectorType = std::vector<std::pair<BoxType, IdentProc>>;

  BoxIdentProcVectorType domain, range;
  BoxType box;
  IdentProc identProc;

  box = BoxType( PointType(proc_id + 0.1, 0.0, 0.0), PointType(proc_id + 0.9, 1.0, 1.0));
  identProc = IdentProc(proc_id * 4, proc_id);
  domain.push_back(std::make_pair(box, identProc));

  box = BoxType( PointType(proc_id + 0.1, 2.0, 0.0), PointType(proc_id + 0.9, 3.0, 1.0));
  identProc = IdentProc(proc_id * 4+1, proc_id);
  domain.push_back(std::make_pair(box, identProc));

  box = BoxType( PointType(proc_id + 0.6, 0.5, 0.0), PointType(proc_id + 1.4, 1.5, 1.0));
  identProc = IdentProc(proc_id * 4+2, proc_id);
  range.push_back(std::make_pair(box, identProc));

  box = BoxType( PointType(proc_id + 0.6, 2.5, 0.0), PointType(proc_id + 1.4, 3.5, 1.0));
  identProc = IdentProc(proc_id * 4+3, proc_id);
  range.push_back(std::make_pair(box, identProc));

  SearchResults searchResults;

  bool enforceSearchResultSymmetry = true;
  bool autoSwapDomainAndRange = true;
  bool sortSearchResults = true;

  stk::search::coarse_search(domain, range, algorithm, comm, searchResults, enforceSearchResultSymmetry, autoSwapDomainAndRange, sortSearchResults);

  expect_search_results(num_procs, proc_id, searchResults);
}

#ifdef KOKKOS_ENABLE_CUDA
void test_coarse_search_with_non_default_view(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  using Box = stk::search::Box<double>;
  using IdentProc = stk::search::IdentProc<int, int>;
  using BoxIdentProc = stk::search::BoxIdentProc<Box, IdentProc>;
  using ViewType1 = Kokkos::View<BoxIdentProc*, Kokkos::CudaSpace>;
  using ViewType2 = Kokkos::View<BoxIdentProc*, Kokkos::CudaUVMSpace>;
  Kokkos::Cuda execSpace{};

  int myrank = stk::parallel_machine_rank(comm);
  ViewType1 domainView("domainView", 1);
  ViewType2 rangeView("rangeView", 1);
  auto createBoxes = KOKKOS_LAMBDA(int idx)
  {
    Box box(2*myrank, 0, 0, 2*myrank+1, 1, 1);
    domainView(idx) = BoxIdentProc{box, IdentProc(0, myrank)};
    rangeView(idx)  = BoxIdentProc{box, IdentProc(0, myrank)};
  };

  Kokkos::parallel_for("create_boxes", 1, createBoxes);

  Kokkos::View<stk::search::IdentProcIntersection<IdentProc, IdentProc>*, Kokkos::CudaSpace> results;
  stk::search::coarse_search(domainView, rangeView, algorithm, comm, results, execSpace);

  auto resultsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);

  EXPECT_EQ(resultsHost.size(), 1U);
  EXPECT_EQ(resultsHost(0).domainIdentProc, IdentProc(0, myrank));
  EXPECT_EQ(resultsHost(0).rangeIdentProc, IdentProc(0, myrank));
}

void test_local_coarse_search_with_non_default_view(stk::search::SearchMethod algorithm)
{
  using Box = stk::search::Box<double>;
  using Ident = int;
  using BoxIdent = stk::search::BoxIdent<Box, Ident>;
  using ViewType1 = Kokkos::View<BoxIdent*, Kokkos::CudaSpace>;
  using ViewType2 = Kokkos::View<BoxIdent*, Kokkos::CudaUVMSpace>;
  Kokkos::Cuda execSpace{};

  ViewType1 domainView("domainView", 1);
  ViewType2 rangeView("rangeView", 1);
  auto createBoxes = KOKKOS_LAMBDA(int idx)
  {
    Box box(0, 0, 0, 1, 1, 1);
    domainView(idx) = BoxIdent{box, Ident(0) };
    rangeView(idx)  = BoxIdent{box, Ident(0)};
  };

  Kokkos::parallel_for("create_boxes", 1, createBoxes);

  Kokkos::View<stk::search::IdentIntersection<Ident, Ident>*, Kokkos::CudaSpace> results;
  stk::search::local_coarse_search(domainView, rangeView, algorithm, results, execSpace);

  auto resultsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);

  EXPECT_EQ(resultsHost.size(), 1U);
  EXPECT_EQ(resultsHost(0).domainIdent, Ident(0));
  EXPECT_EQ(resultsHost(0).rangeIdent,  Ident(0));
}


#endif

TEST(stk_search, coarseSearchDoubleBoxes_KDTREE)
{
  test_coarse_search_for_algorithm<double>(stk::search::KDTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarseSearchFloatBoxes_KDTREE)
{
  test_coarse_search_for_algorithm<float>(stk::search::KDTREE, MPI_COMM_WORLD);
}

TEST(CoarseSearchCorrectness, coarseSearchDoubleBoxes_MORTON_LBVH)
{
  test_coarse_search_for_algorithm<double>(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<double, Kokkos::DefaultExecutionSpace>(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<double, Kokkos::DefaultHostExecutionSpace>(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
}

TEST(CoarseSearchCorrectness, coarseSearchFloatBoxes_MORTON_LBVH)
{
  test_coarse_search_for_algorithm<float>(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<float>(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
}

TEST(CoarseSearchCorrectness, coarseSearchDoubleBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_for_algorithm<double>(stk::search::ARBORX, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<double, Kokkos::DefaultExecutionSpace>(stk::search::ARBORX, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<double, Kokkos::DefaultHostExecutionSpace>(stk::search::ARBORX, MPI_COMM_WORLD);

}

TEST(CoarseSearchCorrectness, coarseSearchFloatBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_for_algorithm<float>(stk::search::ARBORX, MPI_COMM_WORLD);
  test_coarse_search_for_algorithm_with_views<float>(stk::search::ARBORX, MPI_COMM_WORLD);
}

template <typename HostIntersectionType>
void local_expect_search_results(const HostIntersectionType & searchResults)
{
  if constexpr (Kokkos::is_view_v<HostIntersectionType>) {
    ASSERT_EQ(searchResults.extent(0), 2u);

    EXPECT_EQ(searchResults(0).domainIdent, 0);
    EXPECT_EQ(searchResults(0).rangeIdent, 2);

    EXPECT_EQ(searchResults(1).domainIdent, 1);
    EXPECT_EQ(searchResults(1).rangeIdent, 3);
  } else {
    ASSERT_EQ(searchResults.size(), 2u);

    EXPECT_EQ(searchResults[0].first, 0);
    EXPECT_EQ(searchResults[0].second, 2);

    EXPECT_EQ(searchResults[1].first, 1);
    EXPECT_EQ(searchResults[1].second, 3);
  }
}

template <typename FloatType>
void host_local_test_coarse_search_for_algorithm(stk::search::SearchMethod algorithm)
{
  using BoxType = stk::search::Box<FloatType>;
  using PointType = stk::search::Point<FloatType>;
  using BoxIdentType = std::pair<BoxType, Ident>;

  std::vector<BoxIdentType> domain;
  std::vector<BoxIdentType> range;

  domain.emplace_back(BoxType(PointType(0.1, 0.0, 0.0), PointType(0.9, 1.0, 1.0)), 0);
  domain.emplace_back(BoxType(PointType(0.1, 2.0, 0.0), PointType(0.9, 3.0, 1.0)), 1);
  range.emplace_back(BoxType(PointType(0.6, 0.5, 0.0), PointType(1.4, 1.5, 1.0)), 2);
  range.emplace_back(BoxType(PointType(0.6, 2.5, 0.0), PointType(1.4, 3.5, 1.0)), 3);

  LocalSearchResults intersections;

  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, algorithm, intersections, sortSearchResults);

  local_expect_search_results(intersections);
}

template <typename FloatType>
void device_local_test_coarse_search_for_algorithm(stk::search::SearchMethod algorithm)
{
  using BoxType = stk::search::Box<FloatType>;
  using PointType = stk::search::Point<FloatType>;
  using BoxIdentType = stk::search::BoxIdent<BoxType, Ident>;

  auto domain = Kokkos::View<BoxIdentType*, stk::ngp::ExecSpace>("domain box-ident", 2);
  auto range = Kokkos::View<BoxIdentType*, stk::ngp::ExecSpace>("range box-ident", 2);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const unsigned & i) {
      domain[0] = {BoxType(PointType(0.1, 0.0, 0.0), PointType(0.9, 1.0, 1.0)), 0};
      domain[1] = {BoxType(PointType(0.1, 2.0, 0.0), PointType(0.9, 3.0, 1.0)), 1};
      range[0]  = {BoxType(PointType(0.6, 0.5, 0.0), PointType(1.4, 1.5, 1.0)), 2};
      range[1]  = {BoxType(PointType(0.6, 2.5, 0.0), PointType(1.4, 3.5, 1.0)), 3};
    });

  auto intersections = Kokkos::View<IdentIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  auto execSpace = stk::ngp::ExecSpace{};
  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, algorithm, intersections, execSpace, sortSearchResults);

  Kokkos::View<IdentIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  local_expect_search_results(hostIntersections);
}

TEST(CoarseSearchCorrectness, Ngp_Local_CoarseSearchDoubleBoxes_MORTON_LBVH)
{
  host_local_test_coarse_search_for_algorithm<double>(stk::search::MORTON_LBVH);
std::cout<<"finished local test"<<std::endl;
  device_local_test_coarse_search_for_algorithm<double>(stk::search::MORTON_LBVH);
std::cout<<"finished device test"<<std::endl;
}

TEST(CoarseSearchCorrectness, Ngp_Local_CoarseSearchFloatBoxes_MORTON_LBVH)
{
  host_local_test_coarse_search_for_algorithm<float>(stk::search::MORTON_LBVH);
std::cout<<"finished local test"<<std::endl;
  device_local_test_coarse_search_for_algorithm<float>(stk::search::MORTON_LBVH);
std::cout<<"finished device test"<<std::endl;
}

TEST(CoarseSearchCorrectness, Ngp_Local_CoarseSearchDoubleBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  host_local_test_coarse_search_for_algorithm<double>(stk::search::ARBORX);
  device_local_test_coarse_search_for_algorithm<double>(stk::search::ARBORX);
}

TEST(CoarseSearchCorrectness, Ngp_Local_CoarseSearchFloatBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  host_local_test_coarse_search_for_algorithm<float>(stk::search::ARBORX);
  device_local_test_coarse_search_for_algorithm<float>(stk::search::ARBORX);
}


TEST(stk_search, Local_CoarseSearchDoubleBoxes_KDTREE)
{
  host_local_test_coarse_search_for_algorithm<double>(stk::search::KDTREE);
}

TEST(stk_search, Local_CoarseSearchFloatBoxes_KDTREE)
{
  host_local_test_coarse_search_for_algorithm<float>(stk::search::KDTREE);
}

template <typename FloatType, typename ExecSpace>
void local_test_coarse_search_for_algorithm_with_views(stk::search::SearchMethod algorithm)
{
  using BoxType = stk::search::Box<FloatType>;
  using PointType = stk::search::Point<FloatType>;
  using BoxIdentType = stk::search::BoxIdent<BoxType, Ident>;
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  auto domain = Kokkos::View<BoxIdentType*, ExecSpace>("domain box-ident", 2);
  auto range = Kokkos::View<BoxIdentType*, ExecSpace>("range box-ident", 2);

  auto domainHost = Kokkos::create_mirror_view(HostSpace{}, domain);
  auto rangeHost  = Kokkos::create_mirror_view(HostSpace{}, range);


  domainHost(0) = {BoxType(PointType(0.1, 0.0, 0.0), PointType(0.9, 1.0, 1.0)), 0};
  domainHost(1) = {BoxType(PointType(0.1, 2.0, 0.0), PointType(0.9, 3.0, 1.0)), 1};
  rangeHost(0)  = {BoxType(PointType(0.6, 0.5, 0.0), PointType(1.4, 1.5, 1.0)), 2};
  rangeHost(1)  = {BoxType(PointType(0.6, 2.5, 0.0), PointType(1.4, 3.5, 1.0)), 3};

  Kokkos::deep_copy(domain, domainHost);
  Kokkos::deep_copy(range, rangeHost);
  auto intersections = Kokkos::View<IdentIntersection*, ExecSpace>("intersections", 0);

  auto execSpace = ExecSpace{};
  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, algorithm, intersections, execSpace, sortSearchResults);

  auto hostIntersections = Kokkos::create_mirror_view(HostSpace{}, intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  local_expect_search_results(hostIntersections);
}

TEST(stk_search, Local_CoarseSearchWithViews_MORTON_LBVH)
{
  local_test_coarse_search_for_algorithm_with_views<float, Kokkos::DefaultExecutionSpace>(stk::search::MORTON_LBVH);
  local_test_coarse_search_for_algorithm_with_views<float, Kokkos::DefaultHostExecutionSpace>(stk::search::MORTON_LBVH);
}

TEST(stk_search, Local_CoarseSearchWithViews_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif

  local_test_coarse_search_for_algorithm_with_views<float, Kokkos::DefaultExecutionSpace>(stk::search::ARBORX);
  local_test_coarse_search_for_algorithm_with_views<float, Kokkos::DefaultHostExecutionSpace>(stk::search::ARBORX);
}


std::pair<StkBoxIdentProcVector, StkBoxIdentProcVector> build_range_boxes_and_nested_domain_boxes(int num_procs, int proc_id, int sizeParam=1)
{
  StkBoxIdentProcVector local_domain, local_range;

  int startID = 0;
  if (proc_id == 0) {
    for (int i = 0; i < sizeParam; ++i) {
      for (int j = 0; j < sizeParam; ++j, ++startID) {
        StkBox box( Point(num_procs*i, num_procs*j, 0.0), Point(num_procs*(i+1), num_procs*(j+1), 1.0));
        IdentProc id(startID, proc_id);
        local_range.push_back(std::make_pair(box,id));
      }
    }
  }

  for (int procShift = 0; procShift < num_procs; ++procShift) {
    for (int i = 0; i < sizeParam; ++i) {
      for (int j = 0; j < sizeParam; ++j, ++startID) {
        StkBox box( Point(procShift*sizeParam+i, sizeParam*proc_id+j, 0.0), Point(procShift*sizeParam+i+1.0, sizeParam*proc_id+j+1, 1.0));
        IdentProc id(startID, proc_id);
        local_domain.push_back(std::make_pair(box,id));
      }
    }
  }

  return std::make_pair(local_domain, local_range);

}

void expect_coarse_search_range_box_communication(int num_procs,
    int proc_id,
    const SearchResults& searchResultsCommunicateOn,
    const SearchResults& searchResultsCommunicateOff)
{
  if (num_procs == 1) {
    ASSERT_EQ(searchResultsCommunicateOn.size(), 1u);
    ASSERT_EQ(searchResultsCommunicateOff.size(), 1u);
    EXPECT_EQ(searchResultsCommunicateOn[0], std::make_pair(IdentProc(1, 0), IdentProc(0, 0)));
    EXPECT_EQ(searchResultsCommunicateOff[0], std::make_pair(IdentProc(1, 0), IdentProc(0, 0)));
  }

  else {
    bool onRangeBoxProc = (proc_id == 0);
    if (onRangeBoxProc) {
      ASSERT_GT(searchResultsCommunicateOn.size(), searchResultsCommunicateOff.size());
      ASSERT_EQ(searchResultsCommunicateOn.size(), static_cast<size_t>(num_procs * num_procs));
    } else {
      ASSERT_EQ(searchResultsCommunicateOn.size(), static_cast<size_t>(num_procs));
    }
    ASSERT_EQ(searchResultsCommunicateOff.size(), static_cast<size_t>(num_procs));

    for (auto searchResult : searchResultsCommunicateOff) {
      EXPECT_EQ(searchResult.first.proc(), proc_id);
    }
  }
}

void test_coarse_search_range_box_communication(stk::search::SearchMethod algorithm, MPI_Comm comm)
{

  StkBoxIdentProcVector local_domain, local_range;
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  std::tie(local_domain, local_range) = build_range_boxes_and_nested_domain_boxes(num_procs, proc_id);

  SearchResults searchResultsCommunicateOn;
  SearchResults searchResultsCommunicateOff;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsCommunicateOn, true, false);
  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsCommunicateOff, false, false);

  expect_coarse_search_range_box_communication(num_procs, proc_id, searchResultsCommunicateOn, searchResultsCommunicateOff);
}

TEST(stk_search, coarse_search_range_box_communication_KDTREE)
{
  test_coarse_search_range_box_communication(stk::search::KDTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_range_box_communication_MORTON_LBVH)
{
  test_coarse_search_range_box_communication(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
}

void test_coarse_search_determine_domain_and_range_communicate_on(stk::search::SearchMethod algorithm, MPI_Comm comm)
{

  StkBoxIdentProcVector local_domain, local_range;
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  std::tie(local_domain, local_range) = build_range_boxes_and_nested_domain_boxes(num_procs, proc_id);

  SearchResults searchResultsDetermineOn;
  SearchResults searchResultsDetermineOff;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsDetermineOn, true, true, true);
  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsDetermineOff, true, false, true);

  EXPECT_EQ(searchResultsDetermineOn, searchResultsDetermineOff);
}

TEST(stk_search, coarse_search_determine_domain_and_range_communicate_on_KDTREE)
{
  test_coarse_search_determine_domain_and_range_communicate_on(stk::search::KDTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_determine_domain_and_range_communicate_on_MORTON_LBVH)
{
  test_coarse_search_determine_domain_and_range_communicate_on(stk::search::MORTON_LBVH, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_determine_domain_and_range_communicate_on_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_determine_domain_and_range_communicate_on(stk::search::ARBORX, MPI_COMM_WORLD);
}

void test_coarse_search_two_pass(stk::search::SearchMethod algorithm, MPI_Comm comm, int sizeParam)
{
  StkBoxIdentProcVector local_domain, local_range, additional_domain;
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  std::tie(local_domain, local_range) = build_range_boxes_and_nested_domain_boxes(num_procs, proc_id, sizeParam);
  if (proc_id == num_procs - 1) {
    StkBox box( Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0));
    IdentProc id(local_domain.size(), proc_id);
    additional_domain.push_back(std::make_pair(box,id));
  }

  SearchResults searchResultsPassOne;
  SearchResults searchResultsPassTwo;
  SearchResults searchResultsAllBoxes;

  bool communicate = false;
  bool determine = false;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsPassOne, communicate, determine);
  stk::search::coarse_search(additional_domain, local_range, algorithm, comm, searchResultsPassTwo, communicate, determine);

  local_domain.insert(local_domain.end(), additional_domain.begin(), additional_domain.end());

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResultsAllBoxes, communicate, determine);

  searchResultsPassOne.insert(searchResultsPassOne.end(), searchResultsPassTwo.begin(), searchResultsPassTwo.end());

  stk::util::sort_and_unique(searchResultsPassOne);
  stk::util::sort_and_unique(searchResultsAllBoxes);

  EXPECT_EQ(searchResultsPassOne, searchResultsAllBoxes);
}

TEST(stk_search, coarse_search_two_pass_KDTREE)
{
  test_coarse_search_two_pass(stk::search::KDTREE, MPI_COMM_WORLD, 2);
}

TEST(stk_search, coarse_search_two_pass_MORTON_LBVH)
{
  test_coarse_search_two_pass(stk::search::MORTON_LBVH, MPI_COMM_WORLD, 2);
}

TEST(stk_search, coarse_search_two_pass_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_two_pass(stk::search::ARBORX, MPI_COMM_WORLD, 2);
}

void test_ident_proc_with_search_with_views(stk::search::SearchMethod searchMethod)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using BoxIdentProcViewType = Kokkos::View<FloatBoxIdentProc*, ExecSpace>;
  using SearchResultsViewType = Kokkos::View<IdentProcIntersection*, ExecSpace>;

  MPI_Comm comm = MPI_COMM_WORLD;
  int procId = -1;
  MPI_Comm_rank(comm, &procId);
  int numProcs = -1;
  MPI_Comm_size(comm, &numProcs);

  if (numProcs != 1) {
    FloatBox box1(0, 0, 0, 1, 1, 1);
    FloatBox box2(0.5, 0.5, 0.5, 1.5, 1.5, 1.5);
    IdentProc id1(1, 0);
    IdentProc id2(1, 1);

    BoxIdentProcViewType boxes("", 1);
    if (procId == 0) {
      boxes(0) = {box1, id1};
    } else if (procId == 1) {
      boxes(0) = {box2, id2};
    }

    SearchResultsViewType searchResults("", 3);

    auto execSpace = ExecSpace{};
    bool enforceSearchResultSymmetry = true;
    bool autoSwapDomainAndRange = true;
    bool sortSearchResults = true;
    coarse_search(boxes, boxes, searchMethod, comm, searchResults, execSpace,
                         enforceSearchResultSymmetry, autoSwapDomainAndRange, sortSearchResults);

    SearchResultsViewType goldResults("", 3);

    IdentProc goldId1(1, 0);
    IdentProc goldId2(1, 1);

    if (procId == 0) {
      goldResults(0) = {goldId1, goldId1};
      goldResults(1) = {goldId1, goldId2};
      goldResults(2) = {goldId2, goldId1};
      ASSERT_EQ(goldResults.extent(0), searchResults.extent(0));
    } else if (procId == 1) {
      goldResults(0) = {goldId1, goldId2};
      goldResults(1) = {goldId2, goldId1};
      goldResults(2) = {goldId2, goldId2};
      ASSERT_EQ(3u, searchResults.extent(0));
    }

    Kokkos::sort(goldResults, stk::search::Comparator<typename SearchResultsViewType::value_type>());

    for (size_t i = 0; i < goldResults.extent(0); i++) {
      EXPECT_EQ(goldResults[i], searchResults[i])
          << "Test comparison for proc " << procId << " failed for comparsion #" << i << std::endl;
    }
  }
}

void test_ident_proc_with_search(stk::search::SearchMethod searchMethod)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int procId = -1;
  MPI_Comm_rank(comm, &procId);
  int numProcs = -1;
  MPI_Comm_size(comm, &numProcs);

  if (numProcs != 1) {
    FloatBox box1(0, 0, 0, 1, 1, 1);
    FloatBox box2(0.5, 0.5, 0.5, 1.5, 1.5, 1.5);
    IdentProc id1(1, 0);
    IdentProc id2(1, 1);

    FloatBoxIdentProcVector boxes;
    if (procId == 0) {
      boxes.push_back(std::make_pair(box1, id1));
    } else if (procId == 1) {
      boxes.push_back(std::make_pair(box2, id2));
    }

    SearchResults searchResults;

    bool enforceSearchResultSymmetry = true;
    bool autoSwapDomainAndRange = true;
    bool sortSearchResults = true;
    coarse_search(boxes, boxes, searchMethod, comm, searchResults,
                         enforceSearchResultSymmetry, autoSwapDomainAndRange, sortSearchResults);

    SearchResults goldResults;

    IdentProc goldId1(1, 0);
    IdentProc goldId2(1, 1);

    if (procId == 0) {
      goldResults.push_back(std::make_pair(goldId1, goldId1));
      goldResults.push_back(std::make_pair(goldId1, goldId2));
      goldResults.push_back(std::make_pair(goldId2, goldId1));
      ASSERT_EQ(goldResults.size(), searchResults.size());
    } else if (procId == 1) {
      goldResults.push_back(std::make_pair(goldId1, goldId2));
      goldResults.push_back(std::make_pair(goldId2, goldId1));
      goldResults.push_back(std::make_pair(goldId2, goldId2));
      ASSERT_EQ(3u, searchResults.size());
    }

    for (size_t i = 0; i < goldResults.size(); i++) {
      EXPECT_EQ(goldResults[i], searchResults[i])
          << "Test comparison for proc " << procId << " failed for comparsion #" << i << std::endl;
    }
  }
}

TEST(stk_search, coarse_search_ident_proc_switch_KDTREE)
{
  test_ident_proc_with_search(stk::search::KDTREE);
}

TEST(stk_search, coarse_search_ident_proc_switch_MORTON_LBVH)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { GTEST_SKIP(); }
  test_ident_proc_with_search(stk::search::MORTON_LBVH);
  test_ident_proc_with_search_with_views(stk::search::MORTON_LBVH);
}

TEST(stk_search, coarse_search_ident_proc_switch_ARBORX)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { GTEST_SKIP(); }
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_ident_proc_with_search(stk::search::ARBORX);
  test_ident_proc_with_search_with_views(stk::search::ARBORX);
}

void test_coarse_search_one_point(stk::search::SearchMethod searchMethod)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    //int num_procs = stk::parallel_machine_size(comm);
    int proc_id   = stk::parallel_machine_rank(comm);

    Point min_corner, max_corner;

    StkBoxIdentProcVector local_domain, local_range;
    // what if identifier is NOT unique
    // x_min <= x_max
    // y_min <= y_max
    // z_min <= z_max

    min_corner[0] = 0.0; min_corner[1] = 0.0; min_corner[2] = 0.0;
    max_corner[0] = 1.0; max_corner[1] = 1.0; max_corner[2] = 1.0;

    // One bounding box on processor 0 with the label:  0
    // All other processors have empty domain.
    IdentProc domainBox1(0, 0);
    if (proc_id == 0) {
      local_domain.push_back(std::make_pair(StkBox(min_corner, max_corner), domainBox1));
    }

    min_corner[0] = 0.5; min_corner[1] = 0.5; min_corner[2] = 0.5;
    max_corner[0] = 0.5; max_corner[1] = 0.5; max_corner[2] = 0.5;

    // One range target on processor 0 with the label:  1
    // All other processors have empty range.
    IdentProc rangeBox1(1, 0);
    if (proc_id == 0) {
      local_range.push_back(std::make_pair(StkBox(min_corner, max_corner), rangeBox1));
    }

    SearchResults searchResults;

    stk::search::coarse_search(local_domain, local_range, searchMethod, comm, searchResults);

    if (proc_id == 0) {
      ASSERT_EQ(searchResults.size(), 1u);
      EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
    } else {
      ASSERT_EQ(searchResults.size(), 0u);
    }
}

TEST(stk_search, coarse_search_one_point_KDTREE)
{
  test_coarse_search_one_point(stk::search::KDTREE);
}

TEST(stk_search, coarse_search_one_point_MORTON_LBVH)
{
  test_coarse_search_one_point(stk::search::MORTON_LBVH);
}

TEST(stk_search, coarse_search_one_point_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_one_point(stk::search::ARBORX);
}

void test_coarse_search_for_determining_sharing_all_all_case(stk::search::SearchMethod searchMethod)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    typedef std::vector< std::pair<Sphere,IdentProc> > SphereVector;

    SphereVector source_bbox_vector;

    Point coords(1.0, 1.0, 1.0);
    double radius = 1.0e-6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    IdentProc id = IdentProc(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    const bool communicateRangeBoxInfo = true;
    const bool dontCommunicateRangeBoxInfo = false;

    SearchResults searchResults;
    stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                               communicateRangeBoxInfo);
    EXPECT_EQ(2*p_size - 1, static_cast<int>(searchResults.size()));

    searchResults.clear();
    stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                               dontCommunicateRangeBoxInfo);
    EXPECT_EQ(p_size, static_cast<int>(searchResults.size()));

    std::set<int> procs;

    for(size_t i=0;i<searchResults.size();++i)
    {
        procs.insert(searchResults[i].second.proc());
        procs.insert(searchResults[i].first.proc());
    }

    std::set<int>::iterator iter = procs.begin();

    int procCounter = 0;
    for (;iter!=procs.end();++iter)
    {
        EXPECT_EQ(procCounter, *iter);
        procCounter++;
    }
}

TEST(CoarseSearch, forDeterminingSharingAllAllCase_KDTREE)
{
  test_coarse_search_for_determining_sharing_all_all_case(stk::search::KDTREE);
}

TEST(CoarseSearch, forDeterminingSharingAllAllCase_MORTON_LBVH)
{
  test_coarse_search_for_determining_sharing_all_all_case(stk::search::MORTON_LBVH);
}

TEST(CoarseSearch, forDeterminingSharingAllAllCase_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_for_determining_sharing_all_all_case(stk::search::ARBORX);
}

void test_coarse_search_for_determining_sharing_linear_adjacent_case(
    stk::search::SearchMethod searchMethod, int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    typedef std::vector< std::pair<Sphere,IdentProc> > SphereVector;

    SphereVector source_bbox_vector;

    Point coords(1.0 + p_rank, 1.0, 1.0);
    double radius = 0.6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    IdentProc id = IdentProc(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    const bool communicateRangeBoxInfo = true;
    const bool dontCommunicateRangeBoxInfo = false;

    SearchResults searchResults;

    double markTime = 0.0;
    double totalTime = 0.0;

    for (int count = 0; count < numLoops; ++count) {

      const bool lastLoop = (count + 1 == numLoops);

      searchResults.clear();
      markTime = stk::wall_time();
      stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                                 communicateRangeBoxInfo);
      totalTime += stk::wall_time() - markTime;

      if (lastLoop) {
        if (p_size == 1) {
          EXPECT_EQ(1, static_cast<int>(searchResults.size()));
        }
        else {
          const int expectedSize = ((p_rank == 0) || (p_rank == p_size - 1) ? 3 : 5);
          EXPECT_EQ(expectedSize, static_cast<int>(searchResults.size()));
        }
      }

      searchResults.clear();
      markTime = stk::wall_time();
      stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                                 dontCommunicateRangeBoxInfo);
      totalTime += stk::wall_time() - markTime;

      if (lastLoop) {
        if (p_size == 1) {
          EXPECT_EQ(1, static_cast<int>(searchResults.size()));
        }
        else {
          const int expectedSize = ((p_rank == 0) || (p_rank == p_size - 1) ? 2 : 3);
          EXPECT_EQ(expectedSize, static_cast<int>(searchResults.size()));
        }
      }
    }

    if (p_rank == 0) {
      double avgTime =  totalTime / numLoops;
      std::cout << "Average search time measured on rank 0 of " << p_size << " ranks is " << avgTime << " through " << numLoops << " loops" << std::endl;
    }

    // Do more correctness checking after the last loop.
    std::set<int> procs;
    for(size_t i=0;i<searchResults.size();++i)
    {
      procs.insert(searchResults[i].second.proc());
      procs.insert(searchResults[i].first.proc());
    }
    std::set<int>::iterator iter = procs.begin();

    int minNeighbor = p_rank - 1;
    int maxNeighbor = p_rank + 1;

    int procCounter = 0;
    for (;iter!=procs.end();++iter)
    {
      EXPECT_LE(minNeighbor, *iter);
      EXPECT_GE(maxNeighbor, *iter);
      ++procCounter;
    }
    if (p_size == 1) {
      EXPECT_EQ(1, procCounter);
    }
    else {
      const int expectedCount = ((p_rank == 0) || (p_rank == p_size - 1) ? 2 : 3);
      EXPECT_EQ(expectedCount, static_cast<int>(procCounter));
    }
}

TEST(CoarseSearch, forDeterminingSharingLinearAdjacentCase_KDTREE)
{
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::KDTREE);
}

TEST(CoarseSearchScaling, forDeterminingSharingLinearAdjacentCase_KDTREE)
{
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::KDTREE, 1000);
}

TEST(CoarseSearch, forDeterminingSharingLinearAdjacentCase_MORTON_LBVH)
{
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::MORTON_LBVH);
}

TEST(CoarseSearchScaling, forDeterminingSharingLinearAdjacentCase_MORTON_LBVH)
{
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::MORTON_LBVH, 1000);
}

TEST(CoarseSearch, forDeterminingSharingLinearAdjacentCase_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::ARBORX);
}

TEST(CoarseSearchScaling, forDeterminingSharingLinearAdjacentCase_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  test_coarse_search_for_determining_sharing_linear_adjacent_case(stk::search::ARBORX, 1000);
}

#ifdef KOKKOS_ENABLE_CUDA
TEST(CoarseSearch, nonDefaultView_MORTON_LBVH)
{
  test_coarse_search_with_non_default_view(stk::search::MORTON_LBVH, stk::parallel_machine_world());
}

TEST(CoarseSearch, nonDefaultView_ARBORX)
{
  test_coarse_search_with_non_default_view(stk::search::ARBORX, stk::parallel_machine_world());
}

TEST(LocalCoarseSearch, nonDefaultView_MORTON_LBVH)
{
  test_local_coarse_search_with_non_default_view(stk::search::MORTON_LBVH);
}

TEST(LocalCoarseSearch, nonDefaultView_ARBORX)
{
  test_local_coarse_search_with_non_default_view(stk::search::ARBORX);
}
#endif

} //namespace
