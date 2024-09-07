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
#include "stk_search/HelperTraits.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Search.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Tree.hpp"

#include <limits>
#include <type_traits>

#include "Kokkos_Core.hpp"

namespace stk::search {

namespace impl {
template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
class LocalMortonCoarseSearchVectorCallback
{
  using DomainVec = std::vector<std::pair<DomainBoxType, DomainIdentType>>;
  using RangeVec  = std::vector<std::pair<RangeBoxType, RangeIdentType>>;
  using ResultVec = std::vector<std::pair<DomainIdentType, RangeIdentType>>;

  static bool constexpr isSearchExact = !(impl::is_stk_sphere<DomainBoxType> || impl::is_stk_sphere<RangeBoxType>);

  public:
    LocalMortonCoarseSearchVectorCallback(const DomainVec& domain, const RangeVec& range, ResultVec& searchResults) :
      m_domain(domain),
      m_range(range),
      m_searchResults(searchResults)
    {
      m_searchResults.resize(0);
    }

    void operator()(int domainIdx, int rangeIdx) const
    {
#ifdef _OPENMP      
      #pragma omp critical
      {
#endif
      if constexpr (isSearchExact)
      {
        m_searchResults.push_back({m_domain[domainIdx].second, m_range[rangeIdx].second});
      } else
      {
        if (intersects(m_domain[domainIdx].first, m_range[rangeIdx].first))
        {
          m_searchResults.push_back({m_domain[domainIdx].second, m_range[rangeIdx].second});
        }
      }
#ifdef _OPENMP
      }
#endif
    }

    bool resize_for_second_pass()
    {
      return false;
    }

  private:
    const DomainVec& m_domain;
    const RangeVec& m_range;
    ResultVec& m_searchResults;
};

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
  using Callback = impl::LocalMortonCoarseSearchVectorCallback<DomainBoxType, DomainIdentType, RangeBoxType, RangeIdentType>;

  Kokkos::Profiling::pushRegion("Fill domain and range trees");
  const bool supportHostBoxes = false;
  stk::search::MortonAabbTree<ValueType, HostSpace> domainTree("Domain Tree", domain.size(), supportHostBoxes);
  stk::search::MortonAabbTree<ValueType, HostSpace> rangeTree("Range Tree", range.size(), supportHostBoxes);

  stk::search::export_from_box_ident_vector_to_morton_tree(domain, domainTree);
  stk::search::export_from_box_ident_vector_to_morton_tree(range, rangeTree);
  Kokkos::Profiling::popRegion();

  stk::search::CollisionList<HostSpace> collisionList("Collision List");
  Callback callback(domain, range, searchResults);
  stk::search::morton_lbvh_search<ValueType, HostSpace>(domainTree, rangeTree, callback, HostSpace{});
}

namespace impl {
template <typename DomainView, typename RangeView, typename ResultView, typename ExecutionSpace>
class LocalMortonCoarseSearchViewCallback
{

  using DomainBoxIdent = typename DomainView::value_type;
  using RangeBoxIdent  = typename RangeView::value_type;
  using DomainBox      = typename DomainBoxIdent::box_type;
  using RangeBox       = typename RangeBoxIdent::box_type;
  static bool constexpr isSearchExact = !(impl::is_stk_sphere<DomainBox> || impl::is_stk_sphere<RangeBox>);

  public:
    LocalMortonCoarseSearchViewCallback(const DomainView& domain, const RangeView& range, ResultView& searchResults) :
      m_domain(domain),
      m_range(range),
      m_searchResults(searchResults),
      m_idx(Kokkos::ViewAllocateWithoutInitializing("result_idx"))
    {
      check_coarse_search_types_local<DomainView, RangeView, ResultView, ExecutionSpace>();
      Kokkos::deep_copy(m_idx, 0);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int domainIdx, int rangeIdx) const
    {
      if constexpr (isSearchExact)
      {
        unsigned int idx = Kokkos::atomic_fetch_inc(&m_idx());
        if (idx < m_searchResults.size())
        {
          m_searchResults(idx) = {m_domain(domainIdx).ident, m_range(rangeIdx).ident};
        }
      } else
      {
        DomainBoxIdent domainBoxIdent = m_domain(domainIdx);
        RangeBoxIdent rangeBoxIdent   = m_range(rangeIdx);
        if (intersects(domainBoxIdent.box, rangeBoxIdent.box))
        {
          unsigned int idx = Kokkos::atomic_fetch_inc(&m_idx());
          if (idx < m_searchResults.size())
          {
            m_searchResults(idx) = {domainBoxIdent.ident, rangeBoxIdent.ident};
          }
        }
      }
    }

    bool resize_for_second_pass()
    {
      unsigned int numResults = 0;
      Kokkos::deep_copy(numResults, m_idx);
      bool needSecondPass = numResults > m_searchResults.size();
      Kokkos::resize(Kokkos::WithoutInitializing, m_searchResults, numResults);
      Kokkos::deep_copy(m_idx, 0);

      return needSecondPass;
    }

    ResultView get_search_results() const { return m_searchResults; }

  private:
    DomainView m_domain;
    RangeView m_range;
    ResultView m_searchResults;
    Kokkos::View<unsigned int, ExecutionSpace> m_idx;
};

}

template <typename DomainView, typename RangeView, typename ResultView,
          typename ExecutionSpace>
void local_coarse_search_morton_lbvh(
    const DomainView & domain,
    const RangeView & range,
    ResultView & searchResults,
    ExecutionSpace const& execSpace = ExecutionSpace{})
{
  Kokkos::Profiling::pushRegion("local_coarse_search_morton_lbvh");
  check_coarse_search_types_local<DomainView, RangeView, ResultView, ExecutionSpace>();
  using DomainBoxType = typename DomainView::value_type::box_type;
  using RangeBoxType  = typename RangeView::value_type::box_type;
  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");
  using ValueType = typename DomainBoxType::value_type;
  using Callback = impl::LocalMortonCoarseSearchViewCallback<DomainView, RangeView, ResultView, ExecutionSpace>;

  Kokkos::Profiling::pushRegion("STK Fill domain and range trees");
  const bool supportHostBoxes = false;
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> domainTree("Domain Tree", domain.extent(0), supportHostBoxes);
  stk::search::MortonAabbTree<ValueType, ExecutionSpace> rangeTree("Range Tree", range.extent(0), supportHostBoxes);

  Kokkos::Profiling::pushRegion("STK Export box ident views to trees");
  stk::search::export_box_ident_view_to_morton_tree(domain, domainTree, execSpace);
  stk::search::export_box_ident_view_to_morton_tree(range, rangeTree, execSpace);
  execSpace.fence();
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("STK Morton Search");
  if (searchResults.size() == 0)
  {
    size_t sizeGuess = std::max(domain.size(), range.size()) * COLLISION_SCALE_FACTOR;
    Kokkos::resize(Kokkos::WithoutInitializing, searchResults, sizeGuess);
  }

  Callback callback(domain, range, searchResults);
  stk::search::morton_lbvh_search<ValueType, ExecutionSpace>(domainTree, rangeTree, callback, execSpace);
  searchResults = callback.get_search_results();
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();
}

}

#endif // LOCALCOARSESEARCHMORTONLBVH_HPP
