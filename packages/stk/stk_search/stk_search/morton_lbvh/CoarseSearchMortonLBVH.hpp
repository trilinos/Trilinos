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
#include "stk_search/HelperTraits.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <vector>
#include <utility>
#include <type_traits>

namespace stk::search {

namespace impl {
template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
class MortonCoarseSearchVectorCallback
{
  using DomainVec = std::vector<std::pair<DomainBoxType, DomainIdentType>>;
  using RangeVec  = std::vector<std::pair<RangeBoxType, RangeIdentType>>;
  using ExtendedRangeVec = std::vector<RangeBoxType>;
  using RemoteRangeIdentProcVec = std::vector<RangeIdentType>;
  using ResultVec = std::vector<std::pair<DomainIdentType, RangeIdentType>>;

  static bool constexpr isSearchExact = !(impl::is_stk_sphere<DomainBoxType> || impl::is_stk_sphere<RangeBoxType>);

  public:
    MortonCoarseSearchVectorCallback(const DomainVec& localDomain, const RangeVec& localRange,
                                     const ExtendedRangeVec& extendedRangeBoxes,
                                     const RemoteRangeIdentProcVec& remoteRangeIdentProcs,
                                     ResultVec& searchResults) :
      m_localDomain(localDomain),
      m_localRange(localRange),
      m_extendedRangeBoxes(extendedRangeBoxes),
      m_remoteRangeIdentProcs(remoteRangeIdentProcs),
      m_numLocalRange(localRange.size()),
      m_searchResults(searchResults)
    {
      m_searchResults.resize(0);
    }

    void operator()(int domainIdx, int rangeIdx) const
    {
      if constexpr (isSearchExact)
      {
        insert_result(domainIdx, rangeIdx);
      } else
      {
        if (intersects(m_localDomain[domainIdx].first, m_extendedRangeBoxes[rangeIdx]))
        {
          insert_result(domainIdx, rangeIdx);
        }
      }
    }

    bool resize_for_second_pass()
    {
      return false;
    }

  private:
    void insert_result(int domainIdx, int rangeIdx) const
    {
#ifdef _OPENMP      
      #pragma omp critical
      {
#endif        
      if (rangeIdx < m_numLocalRange) {
        m_searchResults.emplace_back(m_localDomain[domainIdx].second, m_localRange[rangeIdx].second);
      }
      else {
        m_searchResults.emplace_back(m_localDomain[domainIdx].second, m_remoteRangeIdentProcs[rangeIdx - m_numLocalRange]);
      }
#ifdef _OPENMP
      }
#endif
    }

    const DomainVec& m_localDomain;
    const RangeVec& m_localRange;
    const ExtendedRangeVec& m_extendedRangeBoxes;
    const RemoteRangeIdentProcVec& m_remoteRangeIdentProcs;
    int m_numLocalRange;
    ResultVec& m_searchResults;
};
}

template <typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
inline void coarse_search_morton_lbvh(std::vector<std::pair<DomainBoxType, DomainIdentProcType>> const & localDomain,
                                      std::vector<std::pair<RangeBoxType, RangeIdentProcType>> const & localRange,
                                      MPI_Comm comm,
                                      std::vector<std::pair<DomainIdentProcType, RangeIdentProcType>> & searchResults,
                                      bool enforceSearchResultSymmetry = true,
                                      bool sortSearchResults = false)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;
  using Callback =  impl::MortonCoarseSearchVectorCallback<DomainBoxType, DomainIdentProcType, RangeBoxType, RangeIdentProcType>;

  STK_ThrowRequireMsg((std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>),
                      "The domain and range boxes must have the same floating-point precision");

  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType = typename RangeBoxType::value_type;

  Kokkos::Profiling::pushRegion("Parallel consistency: extend range box list");
  const auto [extendedRangeBoxes, remoteRangeIdentProcs] =
      morton_extend_local_range_with_remote_boxes_that_might_intersect<HostSpace>(localDomain, localRange, comm, HostSpace{});
  Kokkos::Profiling::popRegion();

  using StkDomainBoxType = stk::search::Box<DomainValueType>;
  using StkRangeBoxType = stk::search::Box<RangeValueType>;

  Kokkos::Profiling::pushRegion("Fill domain and range trees");
  using DomainViewType = Kokkos::View<BoxIdentProc<StkDomainBoxType,DomainIdentProcType>*,HostSpace>;
  using RangeViewType = Kokkos::View<BoxIdentProc<StkRangeBoxType,RangeIdentProcType>*,HostSpace>;
  using DomainTreeType = stk::search::MortonAabbTree<DomainViewType, HostSpace>;
  using RangeTreeType = stk::search::MortonAabbTree<RangeViewType, HostSpace>;
  DomainTreeType domainTree("Domain Tree", localDomain.size());
  RangeTreeType rangeTree("Range Tree", extendedRangeBoxes.size());

  stk::search::export_from_box_ident_proc_vec_to_morton_tree<DomainTreeType,DomainBoxType,DomainIdentProcType>(localDomain, domainTree);
  stk::search::export_from_box_vec_to_morton_tree<RangeTreeType,HostSpace,RangeBoxType>(extendedRangeBoxes, rangeTree);
  domainTree.sync_to_device();
  rangeTree.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Perform Morton query");
  Callback callback(localDomain, localRange, extendedRangeBoxes, remoteRangeIdentProcs, searchResults);
  stk::search::morton_lbvh_search<DomainViewType,RangeViewType,HostSpace,Callback>(domainTree, rangeTree, callback, HostSpace{});
  Kokkos::Profiling::popRegion();

  if (enforceSearchResultSymmetry) {
    Kokkos::Profiling::pushRegion("Enforce results symmetry");
    stk::search::communicate_vector(comm, searchResults, enforceSearchResultSymmetry);
    Kokkos::Profiling::popRegion();
  }

  if (sortSearchResults) {
    Kokkos::Profiling::pushRegion("Sort searchResults");
    std::sort(searchResults.begin(), searchResults.end());
    Kokkos::Profiling::popRegion();
  }
}

namespace impl {
template <typename DomainView,
          typename RangeView,
          typename ResultView,
          typename ExecutionSpace>
class BoundingShapeIntersectionCheckFunctor
{
  public:
    using DomainBoxType       = typename DomainView::value_type::box_type;
    using RangeBoxType        = typename RangeView::value_type::box_type;
    using DomainIdentProcType = typename DomainView::value_type::ident_proc_type;
    using RangeIdentProcType  = typename RangeView::value_type::ident_proc_type;

    using ValueType = typename DomainBoxType::value_type;

    static constexpr bool isBoundingBoxSearchExact = (std::is_base_of_v<DomainBoxType, Box<ValueType>> || std::is_base_of_v<DomainBoxType, Point<ValueType>>) &&
                                                     (std::is_base_of_v<RangeBoxType, Box<ValueType>>  || std::is_base_of_v<RangeBoxType, Point<ValueType>>);

    BoundingShapeIntersectionCheckFunctor(
      DomainView localDomain,
      RangeView localRange,
      const Kokkos::View<RangeBoxType*, ExecutionSpace> extendedRangeBoxes,
      const Kokkos::View<RangeIdentProcType*, ExecutionSpace> remoteIdentProc,
      ResultView searchResults) :
      m_localDomain(localDomain),
      m_localRange(localRange),
      m_extendedRangeBoxes(extendedRangeBoxes),
      m_remoteIdentProc(remoteIdentProc),
      m_searchResults(searchResults),
      m_searchResultIdx("index_in_search_results")
    {
      check_coarse_search_types_parallel<DomainView, RangeView, ResultView, ExecutionSpace>();
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int domainIdx, int rangeIdx) const
    {
      if constexpr (isBoundingBoxSearchExact)
      {
        assign_to_results(domainIdx, rangeIdx);
      } else
      {
        DomainBoxType domainBox = m_localDomain(domainIdx).box;
        RangeBoxType rangeBox   = m_extendedRangeBoxes(rangeIdx);
        if (intersects(domainBox, rangeBox))
        {
          assign_to_results(domainIdx, rangeIdx);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void assign_to_results(int domainIdx, int rangeIdx) const
    {
      unsigned int idx = Kokkos::atomic_fetch_inc(&m_searchResultIdx());
      if (idx < m_searchResults.size())
      {
        if (size_t(rangeIdx) < m_localRange.extent(0))
        {
          m_searchResults(idx) = {m_localDomain(domainIdx).identProc, m_localRange(rangeIdx).identProc};
        } else
        {
          m_searchResults(idx) = {m_localDomain(domainIdx).identProc, m_remoteIdentProc(rangeIdx - m_localRange.extent(0))};
        }
      }
    }

    bool resize_for_second_pass()
    {
      unsigned int numResults = 0;
      Kokkos::deep_copy(numResults, m_searchResultIdx);
      bool needSecondPass = numResults > m_searchResults.size();
      Kokkos::resize(m_searchResults, numResults);
      Kokkos::deep_copy(m_searchResultIdx, 0);

      return needSecondPass;
    }

    ResultView get_search_results() const { return m_searchResults; }

  private:
    DomainView m_localDomain;
    RangeView m_localRange;
    Kokkos::View<RangeBoxType*, ExecutionSpace> m_extendedRangeBoxes;
    Kokkos::View<RangeIdentProcType*, ExecutionSpace> m_remoteIdentProc;
    ResultView m_searchResults;

    Kokkos::View<unsigned int, ExecutionSpace> m_searchResultIdx;
};

}

template <typename DomainView,
          typename RangeView,
          typename ResultView,
          typename ExecutionSpace = typename DomainView::execution_space>
inline void coarse_search_morton_lbvh(
    DomainView const& localDomain,
    RangeView const& localRange,
    MPI_Comm comm,
    ResultView& searchResults,
    ExecutionSpace const& execSpace = ExecutionSpace{},
    bool enforceSearchResultSymmetry = true,
    bool sortSearchResults = false,
    bool doParallelConsistencyOnHost = false)
{
  check_coarse_search_types_parallel<DomainView, RangeView, ResultView, ExecutionSpace>();

  using HostSpace = Kokkos::DefaultHostExecutionSpace;
  using DomainBoxType       = typename DomainView::value_type::box_type;
  using RangeBoxType        = typename RangeView::value_type::box_type;
  using RangeIdentProcType  = typename RangeView::value_type::ident_proc_type;
  using BoundingShapeIntersectionChecker = impl::BoundingShapeIntersectionCheckFunctor<DomainView, RangeView, ResultView, ExecutionSpace>;
  using ExtendedRangeBoxView           = Kokkos::View<RangeBoxType*, ExecutionSpace>;
  using ExtendedRangeIdentProcView     = Kokkos::View<RangeIdentProcType*, ExecutionSpace>;

  static_assert(std::is_same_v<typename DomainBoxType::value_type, typename RangeBoxType::value_type>,
                      "The domain and range boxes must have the same floating-point precision");

  using MDomainViewType = typename BoxIdentProcViewTrait<DomainView>::ViewType;
  using MRangeViewType = typename BoxIdentProcViewTrait<RangeView>::ViewType;

  Kokkos::Profiling::pushRegion("Parallel consistency: extend range box list");
  ExtendedRangeBoxView extendedRangeBoxes;
  ExtendedRangeIdentProcView remoteRangeIdentProcs;

  bool expectedParallelConsistencyFasterOnHost = (localDomain.size() + localRange.size()) < 30000U;
  if (expectedParallelConsistencyFasterOnHost || doParallelConsistencyOnHost)
  {
    HostSpace hostSpace{};
    auto localDomainHost = Kokkos::create_mirror_view_and_copy(hostSpace, localDomain);
    auto localRangeHost = Kokkos::create_mirror_view_and_copy(hostSpace, localRange);

    const auto [extendedRangeBoxesHost, remoteRangeIdentProcsHost] =
        morton_extend_local_range_with_remote_boxes_that_might_intersect(localDomainHost, localRangeHost, hostSpace, comm);

    Kokkos::resize(extendedRangeBoxes, extendedRangeBoxesHost.size());
    Kokkos::resize(remoteRangeIdentProcs, remoteRangeIdentProcsHost.size());

    Kokkos::deep_copy(execSpace, extendedRangeBoxes, extendedRangeBoxesHost);
    Kokkos::deep_copy(execSpace, remoteRangeIdentProcs, remoteRangeIdentProcsHost);
    execSpace.fence();
  } else
  {
    std::tie(extendedRangeBoxes, remoteRangeIdentProcs) =
        morton_extend_local_range_with_remote_boxes_that_might_intersect(localDomain, localRange, execSpace, comm);
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("STK Fill domain and range trees");

  const bool setBoxesOnHost = false;
  using DomainTreeType = stk::search::MortonAabbTree<MDomainViewType, ExecutionSpace>;
  using RangeTreeType = stk::search::MortonAabbTree<MRangeViewType, ExecutionSpace>;
  DomainTreeType domainTree("Domain Tree", localDomain.extent(0), setBoxesOnHost);
  RangeTreeType rangeTree("Range Tree", extendedRangeBoxes.extent(0), setBoxesOnHost);

  if constexpr (std::is_same_v<DomainView, MDomainViewType>) {
    domainTree.m_minMaxs = localDomain;
  }
  else {
    stk::search::export_box_ident_view_to_morton_tree(localDomain, domainTree, execSpace);
  }

  stk::search::export_box_view_to_morton_tree(extendedRangeBoxes, rangeTree, execSpace);
  execSpace.fence();

  domainTree.sync_to_device();
  rangeTree.sync_to_device();
  Kokkos::Profiling::popRegion();

  if (searchResults.size() == 0)
  {
    size_t sizeGuess = std::max(localDomain.size(), extendedRangeBoxes.size()) * COLLISION_SCALE_FACTOR;
    Kokkos::resize(Kokkos::WithoutInitializing, searchResults, sizeGuess);
  }

  Kokkos::Profiling::pushRegion("Inner morton search");
  BoundingShapeIntersectionChecker intersectionChecker(localDomain, localRange, extendedRangeBoxes,
                                                       remoteRangeIdentProcs, searchResults);
  stk::search::morton_lbvh_search<MDomainViewType, MRangeViewType, ExecutionSpace, BoundingShapeIntersectionChecker>(domainTree, rangeTree, intersectionChecker, execSpace);
  searchResults = intersectionChecker.get_search_results();
  Kokkos::Profiling::popRegion();


  if (enforceSearchResultSymmetry) {
    Kokkos::Profiling::pushRegion("Enforce results symmetry");
    SearchResultCommunication<ResultView, ExecutionSpace> resultComm(comm, searchResults, execSpace, enforceSearchResultSymmetry);
    searchResults = resultComm.run();
    Kokkos::Profiling::popRegion();
  }

  if (sortSearchResults) {
    Kokkos::Profiling::pushRegion("Sort searchResults");
    Kokkos::sort(searchResults, Comparator<typename ResultView::value_type>());
    Kokkos::Profiling::popRegion();
  }

}

}

#endif // COARSESEARCHMORTONLBVH_HPP
