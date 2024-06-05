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

#ifndef COARSESEARCHMORTONLBVH_HPP
#define COARSESEARCHMORTONLBVH_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include "Kokkos_Sort.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Search.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Tree.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_ParallelConsistencyUtils.hpp"
#include "stk_search/IdentProc.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/BoundingBox.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <vector>
#include <utility>
#include <type_traits>

namespace stk::search {

template <typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
inline void coarse_search_morton_lbvh(std::vector<std::pair<DomainBoxType, DomainIdentProcType>> const & localDomain,
                                      std::vector<std::pair<RangeBoxType, RangeIdentProcType>> const & localRange,
                                      MPI_Comm comm,
                                      std::vector<std::pair<DomainIdentProcType, RangeIdentProcType>> & searchResults,
                                      bool enforceSearchResultSymmetry = true)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");

  using ValueType = typename DomainBoxType::value_type;

  Kokkos::Profiling::pushRegion("Parallel consistency: extend range box list");
  const auto [extendedRangeBoxes, remoteRangeIdentProcs] =
      morton_extend_local_range_with_remote_boxes_that_might_intersect<HostSpace>(localDomain, localRange, comm, HostSpace{});
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Fill domain and range trees");
  stk::search::MortonAabbTree<ValueType, HostSpace> domainTree("Domain Tree", localDomain.size());
  stk::search::MortonAabbTree<ValueType, HostSpace> rangeTree("Range Tree", extendedRangeBoxes.size());

  stk::search::export_from_box_ident_proc_vec_to_morton_tree(localDomain, domainTree);
  stk::search::export_from_box_vec_to_morton_tree(extendedRangeBoxes, rangeTree);
  domainTree.sync_to_device();
  rangeTree.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Perform Morton query");
  stk::search::CollisionList<HostSpace> collisionList("Collision List");
  stk::search::morton_lbvh_search<ValueType, HostSpace>(domainTree, rangeTree, collisionList, HostSpace{});
  collisionList.sync_from_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Aggregate search results");
  searchResults.clear();
  const unsigned numCollisions = collisionList.hm_idx();
  searchResults.reserve(numCollisions);

  const unsigned numLocalRange = localRange.size();

  auto insert_into_results = [&, &remoteRangeIdentProcs=remoteRangeIdentProcs](unsigned domainIdx, unsigned rangeIdx) {
    if (rangeIdx < numLocalRange) {
      searchResults.emplace_back(localDomain[domainIdx].second, localRange[rangeIdx].second);
    }
    else {
      searchResults.emplace_back(localDomain[domainIdx].second, remoteRangeIdentProcs[rangeIdx - numLocalRange]);
    }
  };

  for (unsigned i = 0; i < numCollisions; ++i) {
    const unsigned domainIdx = collisionList.hm_data(i, 0);
    const unsigned rangeIdx = collisionList.hm_data(i, 1);

    if constexpr ((std::is_same_v<DomainBoxType, Box<ValueType>> || std::is_same_v<DomainBoxType, Point<ValueType>>) &&
                  (std::is_same_v<RangeBoxType, Box<ValueType>> || std::is_same_v<RangeBoxType, Point<ValueType>>)) {
      insert_into_results(domainIdx, rangeIdx);
    }
    else {
      if (intersects(localDomain[domainIdx].first, extendedRangeBoxes[rangeIdx])) {
        insert_into_results(domainIdx, rangeIdx);
      }
    }
  }
  Kokkos::Profiling::popRegion();

  if (enforceSearchResultSymmetry) {
    Kokkos::Profiling::pushRegion("Enforce results symmetry");
    stk::search::communicate_vector(comm, searchResults, enforceSearchResultSymmetry);
    std::sort(searchResults.begin(), searchResults.end());
    Kokkos::Profiling::popRegion();
  }
}

template <typename DomainBoxType,
    typename DomainIdentProcType,
    typename RangeBoxType,
    typename RangeIdentProcType,
    typename ExecutionSpace>
inline void coarse_search_morton_lbvh(
    Kokkos::View<BoxIdentProc<DomainBoxType, DomainIdentProcType>*, ExecutionSpace> const& localDomain,
    Kokkos::View<BoxIdentProc<RangeBoxType, RangeIdentProcType>*, ExecutionSpace> const& localRange,
    MPI_Comm comm,
    Kokkos::View<IdentProcIntersection<DomainIdentProcType, RangeIdentProcType>*, ExecutionSpace>& searchResults,
    ExecutionSpace const& execSpace = ExecutionSpace{},
    bool enforceSearchResultSymmetry = true)
{

  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");

  using ValueType = typename DomainBoxType::value_type;

  Kokkos::Profiling::pushRegion("Move device results to host and convert into compatible data type.");
  auto localDomainHost = Kokkos::create_mirror_view_and_copy(HostSpace{}, localDomain);
  auto localRangeHost  = Kokkos::create_mirror_view_and_copy(HostSpace{}, localRange);

  std::vector<std::pair<DomainBoxType, DomainIdentProcType>> localDomainVec(localDomainHost.size());
  std::vector<std::pair<RangeBoxType, RangeIdentProcType>> localRangeVec(localRangeHost.size());

  for (unsigned i = 0; i < localDomainHost.size(); i++) {
    auto hostPair = localDomainHost(i);
    std::pair<DomainBoxType, DomainIdentProcType> domainPair{hostPair.box, hostPair.identProc};
    localDomainVec[i] = domainPair;
  }

  for (unsigned i = 0; i < localRangeHost.size(); i++) {
    auto hostPair = localRangeHost(i);
    std::pair<RangeBoxType, RangeIdentProcType> rangePair{hostPair.box, hostPair.identProc};
    localRangeVec[i] = rangePair;
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Parallel consistency: extend range box list");
  const auto [extendedRangeBoxes, remoteRangeIdentProcs] =
      morton_extend_local_range_with_remote_boxes_that_might_intersect<ExecutionSpace>(localDomainVec, localRangeVec, comm, execSpace);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Fill domain and range trees");
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> domainTree("Domain Tree", localDomainHost.size());
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> rangeTree("Range Tree", extendedRangeBoxes.size());

  stk::search::export_from_box_ident_proc_vec_to_morton_tree(localDomainVec, domainTree);
  stk::search::export_from_box_vec_to_morton_tree<ValueType, ExecutionSpace>(extendedRangeBoxes, rangeTree);

  domainTree.sync_to_device();
  rangeTree.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Perform Morton query");
  stk::search::CollisionList<ExecutionSpace> collisionList("Collision List");
  stk::search::morton_lbvh_search<ValueType, ExecutionSpace>(domainTree, rangeTree, collisionList, execSpace);
  collisionList.sync_from_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Aggregate search results");
  const unsigned numCollisions = collisionList.hm_idx();
  searchResults = Kokkos::View<IdentProcIntersection<DomainIdentProcType, RangeIdentProcType>*, ExecutionSpace>(
      Kokkos::ViewAllocateWithoutInitializing(searchResults.label()), numCollisions);

  auto searchResultsHost = Kokkos::create_mirror_view_and_copy(HostSpace{}, searchResults);

  const unsigned numLocalRange = localRangeHost.size();
  unsigned searchResultIndex = 0;

  auto insert_into_results = [&, &remoteRangeIdentProcs=remoteRangeIdentProcs](unsigned domainIdx, unsigned rangeIdx, unsigned& searchResultIdx) {
    if (rangeIdx < numLocalRange) {
      searchResultsHost(searchResultIdx) = {localDomainHost(domainIdx).identProc, localRangeHost(rangeIdx).identProc};
      searchResultIdx++;
    }
    else {
      searchResultsHost(searchResultIdx) = {localDomainHost(domainIdx).identProc, remoteRangeIdentProcs[rangeIdx - numLocalRange]};
      searchResultIdx++;
    }
  };

  for (unsigned i = 0; i < numCollisions; ++i) {
    const unsigned domainIdx = collisionList.hm_data(i, 0);
    const unsigned rangeIdx = collisionList.hm_data(i, 1);

    if constexpr ((std::is_same_v<DomainBoxType, Box<ValueType>> || std::is_same_v<DomainBoxType, Point<ValueType>>) &&
                  (std::is_same_v<RangeBoxType, Box<ValueType>> || std::is_same_v<RangeBoxType, Point<ValueType>>)) {
      insert_into_results(domainIdx, rangeIdx, searchResultIndex);
    }
    else {
      if (intersects(localDomainHost(domainIdx).box, extendedRangeBoxes[rangeIdx])) {
        insert_into_results(domainIdx, rangeIdx, searchResultIndex);
      }
    }
  }

  Kokkos::resize(searchResultsHost, searchResultIndex);

  Kokkos::Profiling::popRegion();

  if (enforceSearchResultSymmetry) {
    Kokkos::Profiling::pushRegion("Enforce results symmetry");
    communicate_views(comm, searchResultsHost, enforceSearchResultSymmetry);
    Kokkos::Profiling::popRegion();
  }

  Kokkos::resize(Kokkos::WithoutInitializing, searchResults, searchResultsHost.extent(0));
  Kokkos::deep_copy(searchResults, searchResultsHost);
  Kokkos::sort(searchResults);
}

}

#endif // COARSESEARCHMORTONLBVH_HPP
