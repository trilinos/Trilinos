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

#ifndef LOCAL_COARSE_SEARCH_ARBORX_HPP
#define LOCAL_COARSE_SEARCH_ARBORX_HPP

#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

#include "ArborX.hpp"
#include "Kokkos_Core.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/arborx/StkToArborX.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk::search
{
namespace impl
{
  
template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ValuesViewType,
    typename OffsetsViewType>
void fill_local_search_results(const std::vector<std::pair<DomainBoxType, DomainIdentType>>& localDomain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>>& localRange,
    ValuesViewType values,
    OffsetsViewType offsets,
    std::vector<std::pair<DomainIdentType, RangeIdentType>>& searchResults)
{
  bool constexpr isSphere = is_stk_sphere<DomainBoxType> || is_stk_sphere<RangeBoxType>;
  searchResults.clear();

  for (std::size_t i = 0; i < offsets.extent(0) - 1; ++i) {
    int rangeIndexBegin = offsets(i);
    int rangeIndexEnd = offsets(i + 1);

    for (auto j = rangeIndexBegin; j < rangeIndexEnd; ++j) {
      auto domainIdx = i;
      auto rangeIdx = values(j);

      if (!(isSphere) || intersects(localDomain[domainIdx].first, localRange[rangeIdx].first)) {
        searchResults.emplace_back(localDomain[domainIdx].second, localRange[rangeIdx].second);
      }
    }
  }
}

template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ValuesViewType,
    typename OffsetsViewType,
    typename ExecutionSpace>
void fill_local_search_results(const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace> &  localDomain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace> & localRange,
    ValuesViewType values,
    OffsetsViewType offsets,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace> & searchResults)
{ 

  bool constexpr isSphere = is_stk_sphere<DomainBoxType> || is_stk_sphere<RangeBoxType>;

  int numActualSearchResults = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0, 1), KOKKOS_LAMBDA(int /*index*/, int & searchResultSum) {
      for (std::size_t i = 0; i < offsets.extent(0) - 1; ++i) {
        int rangeIndexBegin = offsets(i);
        int rangeIndexEnd = offsets(i + 1);

      for (auto j = rangeIndexBegin; j < rangeIndexEnd; ++j) {
        auto domainIdx = i;
        auto rangeIdx = values(j);

        if (!(isSphere) || intersects(localDomain[domainIdx].box, localRange[rangeIdx].box)) {
            searchResults(searchResultSum++) = {localDomain[domainIdx].ident, localRange[rangeIdx].ident};
        }
      }
    }
    },
      numActualSearchResults);

    Kokkos::resize(searchResults, numActualSearchResults);

}

}  // namespace impl

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
inline void local_coarse_search_arborx(const std::vector<std::pair<DomainBoxType, DomainIdentType>>& localDomain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>>& localRange,
    std::vector<std::pair<DomainIdentType, RangeIdentType>>& searchResults)
{
  using ExecSpace = Kokkos::DefaultHostExecutionSpace;
  using MemSpace = typename ExecSpace::memory_space;
  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType = typename RangeBoxType::value_type;
  using ArborXDomainType = typename impl::StkToArborX<DomainBoxType>::ArborXType;
  using ArborXRangeType = typename impl::StkToArborX<RangeBoxType>::ArborXType;

  STK_ThrowRequireMsg((std::is_same_v<DomainValueType, RangeValueType>),
      "The domain and range boxes must have the same floating-point precision");

  auto execSpace = ExecSpace{};
  Kokkos::View<ArborXDomainType*, MemSpace> domainBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Domain_BBs"), localDomain.size());
  Kokkos::View<ArborXRangeType*, MemSpace> rangeBoundingBoxes(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "ArborX_Range_BBs"), localRange.size());

  Kokkos::View<int*, MemSpace> values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values"), 0);
  Kokkos::View<ArborX::Intersects<ArborXDomainType>*, MemSpace> queries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Queries"), localDomain.size());
  Kokkos::View<int*, MemSpace> offsets(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Offsets"), 0);

  impl::convert_and_init_bounding_boxes(localDomain, domainBoundingBoxes);
  impl::convert_and_init_bounding_boxes(localRange, rangeBoundingBoxes);

  ArborX::BVH<MemSpace> bvh(execSpace, rangeBoundingBoxes);

  Kokkos::parallel_for(
      "setup_queries", Kokkos::RangePolicy<ExecSpace>(0, localDomain.size()),
      KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(domainBoundingBoxes(i)); });
  bvh.query(execSpace, queries, values, offsets);

  execSpace.fence();

  auto valuesHost = Kokkos::create_mirror_view_and_copy(execSpace, values);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(execSpace, offsets);

  const int numCollisions = values.extent(0);
  searchResults.reserve(numCollisions);

  impl::fill_local_search_results(localDomain, localRange, valuesHost, offsetsHost, searchResults);
}

template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ExecutionSpace>
inline void local_coarse_search_arborx(
    const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace>& localDomain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace>& localRange,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace>& searchResults)
{
  using ExecSpace = ExecutionSpace;
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

  Kokkos::View<int*, MemSpace> values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values"), 0);
  Kokkos::View<ArborX::Intersects<ArborXDomainType>*, MemSpace> queries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Queries"), localDomain.extent(0));
  Kokkos::View<int*, MemSpace> offsets(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Offsets"), 0);

  impl::convert_and_init_bounding_boxes(localDomain, domainBoundingBoxes);
  impl::convert_and_init_bounding_boxes(localRange, rangeBoundingBoxes);

  ArborX::BVH<MemSpace> bvh(execSpace, rangeBoundingBoxes);

  Kokkos::parallel_for(
      "setup_queries", Kokkos::RangePolicy<ExecSpace>(0, localDomain.extent(0)),
      KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(domainBoundingBoxes(i)); });
  bvh.query(execSpace, queries, values, offsets);

  Kokkos::fence();

  const int numCollisions = values.extent(0);
  searchResults = Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace>(
        Kokkos::ViewAllocateWithoutInitializing(searchResults.label()), numCollisions);

  impl::fill_local_search_results(localDomain, localRange, values, offsets, searchResults);

}

}  // namespace stk::search

#endif
