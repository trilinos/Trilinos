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
#include "stk_search/arborx/AccessTraits.hpp"

namespace stk::search
{
namespace impl
{

template <typename DomainViewType, typename RangeViewType>
class SearchResultsInserter
{
  public:
    static_assert(Kokkos::is_view_v<DomainViewType>);
    static_assert(Kokkos::is_view_v<RangeViewType>);

    using DomainBoxIdent = typename DomainViewType::value_type;
    using RangeBoxIdent  = typename RangeViewType::value_type;
    using DomainBox      = typename DomainBoxIdent::box_type;
    using DomainIdent    = typename DomainBoxIdent::ident_type;
    using RangeBox       = typename RangeBoxIdent::box_type;
    using RangeIdent     = typename RangeBoxIdent::ident_type;
    using ExecutionSpace = typename DomainViewType::execution_space;
    using SearchResultsView = Kokkos::View<IdentIntersection<DomainIdent, RangeIdent>*, ExecutionSpace>;
    using ArborXPredicateWithIndex = typename ArborX::AccessTraits<impl::ViewWrapperForArborXTraits<DomainViewType>, ArborX::PredicatesTag>::ArborXPredicateWithIndex;

    static bool constexpr isSphere = impl::is_stk_sphere<DomainBox> || impl::is_stk_sphere<RangeBox>;

    SearchResultsInserter(DomainViewType localDomain, RangeViewType localRange, SearchResultsView searchResults) :
      m_localDomain(localDomain),
      m_localRange(localRange),
      m_searchResults(searchResults),
      m_counter("counter"),
      m_searchResultsSize(0)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const ArborXPredicateWithIndex& predicate, int rangeBoxIdx) const
    {
      int domainBoxIdx = ArborX::getData(predicate);
      const DomainBoxIdent domainBoxIdent = m_localDomain(domainBoxIdx);
      const RangeBoxIdent rangeBoxIdent   = m_localRange(rangeBoxIdx);
      if (!(isSphere) || intersects(domainBoxIdent.box, rangeBoxIdent.box))
      {
        int idx = Kokkos::atomic_fetch_add(&(m_counter()), 1);
        if (idx < m_searchResultsSize)
        {
          m_searchResults(idx) = {domainBoxIdent.ident, rangeBoxIdent.ident};
        }
      }
    }

    bool resizeSearchResults()
    {
      int initialResultsSize = m_searchResultsSize;
      Kokkos::deep_copy(m_searchResultsSize, m_counter);
      bool needToRunAgain = m_searchResultsSize > initialResultsSize;

      Kokkos::resize(m_searchResults, m_searchResultsSize);
      if (needToRunAgain)
      {
        Kokkos::deep_copy(m_counter, 0);
      }

      return needToRunAgain;
    }

    SearchResultsView getSearchResults() { return m_searchResults; }

  private:
    DomainViewType m_localDomain;
    RangeViewType m_localRange;
    SearchResultsView m_searchResults;
    Kokkos::View<int, ExecutionSpace> m_counter;
    int m_searchResultsSize;
};

template <typename DomainViewType, typename RangeViewType>
SearchResultsInserter<DomainViewType, RangeViewType> create_results_inserter(DomainViewType localDomain, RangeViewType localRange,
                                                       typename SearchResultsInserter<DomainViewType, RangeViewType>::SearchResultsView searchResults)
{
  return SearchResultsInserter<DomainViewType, RangeViewType>(localDomain, localRange, searchResults);
}

}  // namespace impl

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
inline void local_coarse_search_arborx(const std::vector<std::pair<DomainBoxType, DomainIdentType>>& localDomain,
    const std::vector<std::pair<RangeBoxType, RangeIdentType>>& localRange,
    std::vector<std::pair<DomainIdentType, RangeIdentType>>& searchResults)
{
  using ExecSpace  = Kokkos::DefaultHostExecutionSpace;
  using MemSpace   = typename ExecSpace::memory_space;
  using DomainView = Kokkos::View<const std::pair<DomainBoxType, DomainIdentType>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using RangeView  = Kokkos::View<const std::pair<RangeBoxType, RangeIdentType>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using Access     = typename ArborX::AccessTraits<impl::ViewWrapperForArborXTraits<DomainView>, ArborX::PredicatesTag>;
  using ArborXPredicateWithIndex = typename Access::ArborXPredicateWithIndex;

  static bool constexpr isSphere = impl::is_stk_sphere<DomainBoxType> || impl::is_stk_sphere<RangeBoxType>;

  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
      "The domain and range boxes must have the same floating-point precision");

  auto execSpace = ExecSpace{};

  DomainView localDomainView(Kokkos::view_wrap(localDomain.data()), localDomain.size());
  RangeView localRangeView(Kokkos::view_wrap(localRange.data()), localRange.size());
  auto localDomainViewWrapped = impl::wrap_view_for_arborx(localDomainView);
  auto localRangeViewWrapped  = impl::wrap_view_for_arborx(localRangeView);

  searchResults.clear();
  auto callback = [&](const ArborXPredicateWithIndex& predicate, int rangeBoxIdx)
  {
    int domainBoxIdx = ArborX::getData(predicate);
    const auto& domainBoxIdent = localDomain[domainBoxIdx];
    const auto& rangeBoxIdent   = localRange[rangeBoxIdx];

    if (!(isSphere) || intersects(domainBoxIdent.first, rangeBoxIdent.first))
    {
      searchResults.emplace_back(domainBoxIdent.second, rangeBoxIdent.second);
    }
  };

  ArborX::BVH<MemSpace> bvh(execSpace, localRangeViewWrapped);
  bvh.query(execSpace, localDomainViewWrapped, callback);
  execSpace.fence();
}

template <typename DomainBoxType,
    typename DomainIdentType,
    typename RangeBoxType,
    typename RangeIdentType,
    typename ExecutionSpace>
inline void local_coarse_search_arborx(
    const Kokkos::View<BoxIdent<DomainBoxType, DomainIdentType>*, ExecutionSpace>& localDomain,
    const Kokkos::View<BoxIdent<RangeBoxType, RangeIdentType>*, ExecutionSpace>& localRange,
    Kokkos::View<IdentIntersection<DomainIdentType, RangeIdentType>*, ExecutionSpace>& searchResults,
    ExecutionSpace const& execSpace = ExecutionSpace{})
{
  using ExecSpace       = ExecutionSpace;
  using MemSpace        = typename ExecSpace::memory_space;
  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType  = typename RangeBoxType::value_type;

  STK_ThrowRequireMsg((std::is_same_v<DomainValueType, RangeValueType>),
      "The domain and range boxes must have the same floating-point precision");

  Kokkos::Profiling::pushRegion("STK call arborx");

  auto localRangeWrapped = impl::wrap_view_for_arborx(localRange);
  auto localDomainWrapped = impl::wrap_view_for_arborx(localDomain);
  auto callback = impl::create_results_inserter(localDomain, localRange, searchResults);

  ArborX::BVH<MemSpace> bvh(execSpace, localRangeWrapped);
  bvh.query(execSpace, localDomainWrapped, callback);

  bool runSecondPass = callback.resizeSearchResults();
  if (runSecondPass)
  {
    bvh.query(execSpace, localDomainWrapped, callback);
  }

  searchResults = callback.getSearchResults();

  Kokkos::fence();
  Kokkos::Profiling::popRegion();
}

}  // namespace stk::search

#endif
