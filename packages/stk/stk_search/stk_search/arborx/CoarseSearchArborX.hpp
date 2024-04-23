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

#ifndef COARSE_SEARCH_ARBORX_HPP
#define COARSE_SEARCH_ARBORX_HPP

#include <tuple>
#include <vector>

#include "ArborX.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/arborx/StkToArborX.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk::search
{
namespace impl
{
template <typename DomainBoxPairType, typename RangeBoxPairType, typename SearchResultsListType>
void fill_search_result(DomainBoxPairType const& domainBox,
    RangeBoxPairType const& rangeBox,
    SearchResultsListType& searchResults,
    [[maybe_unused]] int& searchResultIdx)
{
  constexpr bool usesSpheres =
      is_stk_sphere<typename DomainBoxPairType::first_type> || is_stk_sphere<typename RangeBoxPairType::first_type>;
  if constexpr (!Kokkos::is_view_v<SearchResultsListType>) {
    if constexpr (usesSpheres) {
      if (intersects(domainBox.first, rangeBox.first)) searchResults.emplace_back(domainBox.second, rangeBox.second);
    } else {
      searchResults.emplace_back(domainBox.second, rangeBox.second);
    }
  } else {
    if constexpr (usesSpheres) {
      if (intersects(domainBox.box, rangeBox.box))
        searchResults[searchResultIdx++] = {domainBox.identProc, rangeBox.identProc};
    } else {
      searchResults[searchResultIdx++] = {domainBox.identProc, rangeBox.identProc};
    }
  }
}

template <typename LocalDomainListType,
    typename LocalRangeListType,
    typename ValuesHostViewType,
    typename OffsetsHostViewType,
    typename SearchResultsListType>
void fill_search_results(LocalDomainListType const& localDomain,
    LocalRangeListType const& localRange,
    ValuesHostViewType valuesHost,
    OffsetsHostViewType offsetsHost,
    MPI_Comm comm,
    SearchResultsListType& searchResults)
{
  int searchResultIdx = 0;
  int myRank = 0;
  int commSize = 0;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &commSize);

  if (stk::parallel_machine_size(comm) > 1) {
    stk::CommSparse commSparse(comm);
    stk::CommSparse commSparse2(comm);

    for (std::size_t i = 0; i < offsetsHost.extent(0) - 1; ++i) {
      int rangeIndexBegin = offsetsHost(i);
      int rangeIndexEnd = offsetsHost(i + 1);

      for (auto j = rangeIndexBegin; j < rangeIndexEnd; ++j) {
        auto domainIdx = i;
        auto rangeIdx = valuesHost(j).index;
        auto rangeProc = valuesHost(j).rank;
        if (rangeProc == myRank) {
          fill_search_result(localDomain[domainIdx], localRange[rangeIdx], searchResults, searchResultIdx);
        }
      }
    }

    stk::pack_and_communicate(commSparse, [&]() {
      for (std::size_t i = 0; i < offsetsHost.extent(0) - 1; ++i) {
        int rangeIndexBegin = offsetsHost(i);
        int rangeIndexEnd = offsetsHost(i + 1);

        for (auto j = rangeIndexBegin; j < rangeIndexEnd; ++j) {
          auto domainIdx = i;
          auto rangeIdx = valuesHost(j).index;
          auto rangeProc = valuesHost(j).rank;

          if (rangeProc == myRank) {
            continue;
          } else {
            commSparse.send_buffer(rangeProc).pack((int) domainIdx);
            commSparse.send_buffer(rangeProc).pack(rangeIdx);
          }
        }
      }
    });

    // domainIdx, domainProc, rangeIdx
    std::vector<std::tuple<int, int, int>> pairingInfo;

    for (int proc = 0; proc < commSize; ++proc) {
      if (proc == myRank) {
        continue;
      }
      while (commSparse.recv_buffer(proc).remaining()) {
        int domainIdx;
        int rangeIdx;
        commSparse.recv_buffer(proc).unpack(domainIdx);
        commSparse.recv_buffer(proc).unpack(rangeIdx);
        pairingInfo.emplace_back(domainIdx, proc, rangeIdx);
      }
    }

    stk::pack_and_communicate(commSparse2, [&]() {
      for (auto info : pairingInfo) {
        commSparse2.send_buffer(std::get<1>(info))
            .pack(std::make_pair(std::get<0>(info), localRange[std::get<2>(info)]));
      }
    });

    for (int proc = 0; proc < commSize; ++proc) {
      if (proc == myRank) {
        continue;
      }
      while (commSparse2.recv_buffer(proc).remaining()) {
        std::pair<int, typename LocalRangeListType::value_type> recv;
        commSparse2.recv_buffer(proc).unpack(recv);

        auto domainIdx = recv.first;
        auto recvRangeBox = recv.second;

        fill_search_result(localDomain[domainIdx], recvRangeBox, searchResults, searchResultIdx);
      }
    }

  } else {
    for (std::size_t i = 0; i < offsetsHost.extent(0) - 1; ++i) {
      int rangeIndexBegin = offsetsHost(i);
      int rangeIndexEnd = offsetsHost(i + 1);

      for (auto j = rangeIndexBegin; j < rangeIndexEnd; ++j) {
        auto domainIdx = i;
        auto rangeIdx = valuesHost(j).index;

        fill_search_result(localDomain[domainIdx], localRange[rangeIdx], searchResults, searchResultIdx);
      }
    }
  }

  if constexpr (Kokkos::is_view_v<SearchResultsListType>) {
    Kokkos::resize(searchResults, searchResultIdx);
  }

}
}  // namespace impl

template <typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
inline void coarse_search_arborx(std::vector<std::pair<DomainBoxType, DomainIdentProcType>> const& localDomain,
    std::vector<std::pair<RangeBoxType, RangeIdentProcType>> const& localRange,
    MPI_Comm comm,
    std::vector<std::pair<DomainIdentProcType, RangeIdentProcType>>& searchResults,
    bool enforceSearchResultSymmetry = true)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using MemSpace = typename ExecSpace::memory_space;
  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType = typename RangeBoxType::value_type;
  using ArborXDomainType = typename impl::StkToArborX<DomainBoxType>::ArborXType;
  using ArborXRangeType = typename impl::StkToArborX<RangeBoxType>::ArborXType;

  STK_ThrowRequireMsg((std::is_same_v<DomainValueType, RangeValueType>),
      "The domain and range boxes must have the same floating-point precision");

  Kokkos::View<ArborXDomainType*, MemSpace> domainBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Domain_BBs"), localDomain.size());
  Kokkos::View<ArborXRangeType*, MemSpace> rangeBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Range_BBs"), localRange.size());
  Kokkos::View<ArborX::PairIndexRank*, MemSpace> values(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Indices_And_Ranks"), 0);
  Kokkos::View<ArborX::Intersects<ArborXDomainType>*, MemSpace> queries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Queries"), localDomain.size());
  Kokkos::View<int*, MemSpace> offsets(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Offsets"), 0);

  impl::convert_and_init_bounding_boxes(localDomain, domainBoundingBoxes);
  impl::convert_and_init_bounding_boxes(localRange, rangeBoundingBoxes);

  ArborX::DistributedTree<MemSpace> tree(comm, ExecSpace{}, rangeBoundingBoxes);

  Kokkos::parallel_for(
      "setup_queries", Kokkos::RangePolicy<ExecSpace>(0, localDomain.size()),
      KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(domainBoundingBoxes(i)); });
  tree.query(ExecSpace{}, queries, values, offsets);

  Kokkos::fence();

  auto valuesHost = Kokkos::create_mirror_view_and_copy(ExecSpace{}, values);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(ExecSpace{}, offsets);

  const int numCollisions = values.extent(0);
  searchResults.clear();
  searchResults.reserve(numCollisions);

  impl::fill_search_results(localDomain, localRange, valuesHost, offsetsHost, comm, searchResults);

  if (enforceSearchResultSymmetry) {
    communicate_vector(comm, searchResults, enforceSearchResultSymmetry);
    std::sort(searchResults.begin(), searchResults.end());
  }
}

template <typename DomainBoxType,
    typename DomainIdentProcType,
    typename RangeBoxType,
    typename RangeIdentProcType,
    typename ExecutionSpace>
inline void coarse_search_arborx(
    Kokkos::View<BoxIdentProc<DomainBoxType, DomainIdentProcType>*, ExecutionSpace> const& localDomain,
    Kokkos::View<BoxIdentProc<RangeBoxType, RangeIdentProcType>*, ExecutionSpace> const& localRange,
    MPI_Comm comm,
    Kokkos::View<IdentProcIntersection<DomainIdentProcType, RangeIdentProcType>*, ExecutionSpace>& searchResults,
    bool enforceSearchResultSymmetry = true)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using MemSpace = typename ExecSpace::memory_space;
  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType = typename RangeBoxType::value_type;
  using ArborXDomainType = typename impl::StkToArborX<DomainBoxType>::ArborXType;
  using ArborXRangeType = typename impl::StkToArborX<RangeBoxType>::ArborXType;

  STK_ThrowRequireMsg((std::is_same_v<DomainValueType, RangeValueType>),
      "The domain and range boxes must have the same floating-point precision");

  auto execSpace = ExecSpace{};

  Kokkos::View<ArborXDomainType*, MemSpace> domainBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Domain_BBs"), localDomain.extent(0));
  Kokkos::View<ArborXRangeType*, MemSpace> rangeBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Range_BBs"), localRange.extent(0));
  Kokkos::View<ArborX::PairIndexRank*, MemSpace> values(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Indices_And_Ranks"), 0);
  Kokkos::View<ArborX::Intersects<ArborXDomainType>*, MemSpace> queries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Queries"), localDomain.extent(0));
  Kokkos::View<int*, MemSpace> offsets(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Offsets"), 0);

  impl::convert_and_init_bounding_boxes(localDomain, domainBoundingBoxes);
  impl::convert_and_init_bounding_boxes(localRange, rangeBoundingBoxes);

  ArborX::DistributedTree<MemSpace> tree(comm, execSpace, rangeBoundingBoxes);

  Kokkos::parallel_for(
      "setup_queries", Kokkos::RangePolicy<ExecSpace>(0, localDomain.extent(0)),
      KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(domainBoundingBoxes(i)); });
  tree.query(execSpace, queries, values, offsets);

  Kokkos::fence();

  auto valuesHost = Kokkos::create_mirror_view_and_copy(execSpace, values);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(execSpace, offsets);

  const int numCollisions = values.extent(0);
  searchResults = Kokkos::View<IdentProcIntersection<DomainIdentProcType, RangeIdentProcType>*, ExecSpace>(
      Kokkos::ViewAllocateWithoutInitializing(searchResults.label()), numCollisions);

  auto localDomainHost = Kokkos::create_mirror_view_and_copy(execSpace, localDomain);
  auto localRangeHost = Kokkos::create_mirror_view_and_copy(execSpace, localRange);
  auto searchResultsHost = Kokkos::create_mirror_view_and_copy(execSpace, searchResults);

  impl::fill_search_results(localDomainHost, localRangeHost, valuesHost, offsetsHost, comm, searchResultsHost);
  Kokkos::sort(searchResultsHost);

  if (enforceSearchResultSymmetry) {
    communicate_views(comm, searchResultsHost, enforceSearchResultSymmetry);
    Kokkos::resize(Kokkos::WithoutInitializing, searchResults, searchResultsHost.extent(0));
  }

  
  Kokkos::deep_copy(searchResults, searchResultsHost);
}

}  // namespace stk::search

#endif
