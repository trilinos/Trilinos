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

#ifndef LOCALCOARSESEARCHMORTONLBVH_HPP
#define LOCALCOARSESEARCHMORTONLBVH_HPP

#include "stk_search/Box.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Search.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Tree.hpp"

#include <limits>
#include <type_traits>

#include "Kokkos_Core.hpp"

namespace stk::search {

template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ExecutionSpace>
void insert_intersections_into_results(const std::vector<std::pair<DomainBoxType, DomainIdentType>>& domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>>& range,
    const stk::search::CollisionList<ExecutionSpace> rawIntersections,
    std::vector<std::pair<DomainIdentType, RangeIdentType>>& intersections)
{
  const int numCollisions = rawIntersections.get_num_collisions();

  for (int i = 0; i < numCollisions; ++i) {
    const unsigned domainIdx = rawIntersections.m_data(i, 0);
    const unsigned rangeIdx = rawIntersections.m_data(i, 1);
    intersections.emplace_back(domain[domainIdx].second, range[rangeIdx].second);
  };
}

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType,
          typename ExecutionSpace>
void insert_intersections_into_results(
    const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace> & domain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace> & range,
    const stk::search::CollisionList<ExecutionSpace> rawIntersections,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace> & intersections,
    ExecutionSpace const& execSpace)
{
  const int numCollisions = rawIntersections.get_num_collisions();

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, numCollisions),
      KOKKOS_LAMBDA(int index) {
        const unsigned domainIdx = rawIntersections.m_data(index, 0);
        const unsigned rangeIdx = rawIntersections.m_data(index, 1);
        intersections[index] = {domain[domainIdx].ident, range[rangeIdx].ident};
      });
}

template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ExecutionSpace>
void insert_only_confirmed_intersections_into_results(
    const std::vector<std::pair<DomainBoxType, DomainIdentType>>& domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>>& range,
    const stk::search::CollisionList<ExecutionSpace> rawIntersections,
    std::vector<std::pair<DomainIdentType, RangeIdentType>>& intersections)
{
  constexpr DomainIdentType INVALID_DOMAIN_IDENT = std::numeric_limits<DomainIdentType>::max();
  constexpr RangeIdentType INVALID_RANGE_IDENT = std::numeric_limits<RangeIdentType>::max();
  const int numCollisions = rawIntersections.get_num_collisions();

  for (int index = 0; index < numCollisions; ++index) {
    const unsigned domainIdx = rawIntersections.m_data(index, 0);
    const unsigned rangeIdx = rawIntersections.m_data(index, 1);
    const auto& domainBoxIdent = domain[domainIdx];
    const auto& rangeBoxIdent = range[rangeIdx];
    if (intersects(domainBoxIdent.first, rangeBoxIdent.first)) {
      intersections.emplace_back(domainBoxIdent.second, rangeBoxIdent.second);
    } else {
      intersections.emplace_back(INVALID_DOMAIN_IDENT, INVALID_RANGE_IDENT);
    }
  };

  int numActualIntersections = 0;
  int destIndex = 0;
  for (int sourceIndex = 0; sourceIndex < numCollisions; ++sourceIndex) {
    if (intersections[sourceIndex].first != INVALID_DOMAIN_IDENT) {
      intersections[destIndex++] = intersections[sourceIndex];
      numActualIntersections++;
    }
  }

  intersections.resize(numActualIntersections);
}

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType,
          typename ExecutionSpace>
void insert_only_confirmed_intersections_into_results(
    const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace> & domain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace> & range,
    const stk::search::CollisionList<ExecutionSpace> rawIntersections,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace> & intersections,
    ExecutionSpace const& execSpace)
{
  static bool constexpr isSphere = impl::is_stk_sphere<DomainBoxType> || impl::is_stk_sphere<RangeBoxType>;

  const int numCollisions = rawIntersections.get_num_collisions();

  Kokkos::View<int, ExecutionSpace> counter("counter");
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, numCollisions),
      KOKKOS_LAMBDA(int index) {
        const unsigned domainIdx = rawIntersections.m_data(index, 0);
        const unsigned rangeIdx = rawIntersections.m_data(index, 1);
        const auto domainBoxIdent = domain[domainIdx];
        const auto rangeBoxIdent = range[rangeIdx];

        if (isSphere)
        {
          if (intersects(domainBoxIdent.box, rangeBoxIdent.box))
          {
            int outputIdx = Kokkos::atomic_fetch_add(&(counter()), 1);
            intersections[outputIdx] = {domainBoxIdent.ident, rangeBoxIdent.ident};
          }
        } else
        {
          intersections[index] = {domainBoxIdent.ident, rangeBoxIdent.ident};
        }
      });

  if constexpr (isSphere)
  {
    int numActualIntersections;
    Kokkos::deep_copy(numActualIntersections, counter);
    Kokkos::resize(intersections, numActualIntersections);
  }
}

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
void local_coarse_search_morton_lbvh(
    const std::vector<std::pair<DomainBoxType, DomainIdentType>> & domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>> & range,
    std::vector<std::pair<DomainIdentType, RangeIdentType>> & searchResults)
{
  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");

  using ValueType = typename DomainBoxType::value_type;
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  Kokkos::Profiling::pushRegion("Fill domain and range trees");
  const bool supportHostBoxes = false;
  stk::search::MortonAabbTree<ValueType, HostSpace> domainTree("Domain Tree", domain.size(), supportHostBoxes);
  stk::search::MortonAabbTree<ValueType, HostSpace> rangeTree("Range Tree", range.size(), supportHostBoxes);

  stk::search::export_from_box_ident_vector_to_morton_tree(domain, domainTree);
  stk::search::export_from_box_ident_vector_to_morton_tree(range, rangeTree);
  Kokkos::Profiling::popRegion();

  stk::search::CollisionList<HostSpace> collisionList("Collision List");
  stk::search::morton_lbvh_search<ValueType, HostSpace>(domainTree, rangeTree, collisionList, HostSpace{});

  Kokkos::Profiling::pushRegion("Aggregate search results");
  const int numCollisions = collisionList.get_num_collisions();
  searchResults.reserve(numCollisions);

  if constexpr ((std::is_same_v<DomainBoxType, Box<ValueType>> || std::is_same_v<DomainBoxType, Point<ValueType>>) &&
                (std::is_same_v<RangeBoxType, Box<ValueType>> || std::is_same_v<RangeBoxType, Point<ValueType>>))
  {
    insert_intersections_into_results(domain, range, collisionList, searchResults);
  }
  else {
    insert_only_confirmed_intersections_into_results(domain, range, collisionList, searchResults);
  }

  Kokkos::Profiling::popRegion();
}

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType,
          typename ExecutionSpace>
void local_coarse_search_morton_lbvh(
    const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace> & domain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace> & range,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace> & searchResults,
    ExecutionSpace const& execSpace = ExecutionSpace{})
{
  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");

  using ValueType = typename DomainBoxType::value_type;

  Kokkos::Profiling::pushRegion("STK Fill domain and range trees");
  const bool supportHostBoxes = false;
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> domainTree("Domain Tree", domain.extent(0), supportHostBoxes);
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> rangeTree("Range Tree", range.extent(0), supportHostBoxes);

  stk::search::export_from_box_ident_view_to_morton_tree(domain, domainTree);
  stk::search::export_from_box_ident_view_to_morton_tree(range, rangeTree);
  domainTree.sync_to_device();
  rangeTree.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("STK Morton Search");
  stk::search::CollisionList<ExecutionSpace> collisionList("Collision List");
  stk::search::morton_lbvh_search<ValueType, ExecutionSpace>(domainTree, rangeTree, collisionList, execSpace);

  Kokkos::Profiling::pushRegion("STK Aggregate search results");
  const int numCollisions = collisionList.get_num_collisions();
  searchResults = Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace>(
        Kokkos::ViewAllocateWithoutInitializing(searchResults.label()), numCollisions);

  if constexpr ((std::is_same_v<DomainBoxType, Box<ValueType>> || std::is_same_v<DomainBoxType, Point<ValueType>>) &&
                (std::is_same_v<RangeBoxType, Box<ValueType>> || std::is_same_v<RangeBoxType, Point<ValueType>>))
  {
    insert_intersections_into_results(domain, range, collisionList, searchResults, execSpace);
  }
  else {
    insert_only_confirmed_intersections_into_results(domain, range, collisionList, searchResults, execSpace);
  }

  Kokkos::Profiling::popRegion();
}

}

#endif // LOCALCOARSESEARCHMORTONLBVH_HPP
