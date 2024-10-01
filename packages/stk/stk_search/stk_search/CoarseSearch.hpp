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

#ifndef stk_search_CoarseSearch_hpp
#define stk_search_CoarseSearch_hpp

#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/kdtree/CoarseSearchKdTree.hpp>
#include <stk_search/morton_lbvh/CoarseSearchMortonLBVH.hpp>
#ifdef STK_HAS_ARBORX
#include <stk_search/arborx/CoarseSearchArborX.hpp>
#endif

#include <stk_search/SearchMethod.hpp>
#include <stk_search/CommonSearchUtil.hpp>
#include <stk_search/HelperTraits.hpp>
#include <Kokkos_Core.hpp>
#include <vector>
#include <utility>

namespace stk::search {

// After providing a vector of local domain and range search box information
// (that includes both the box itself and the target entity ID and owning
// processor), this function will compute a parallel-consistent list of
// all box intersections between local entities (in the domain vector) and
// range entities from all processors.
//
// The "enforceSearchResultSymmetry" argument indicates whether you want to
// perform an additional round of parallel communication to make sure that
// no intersections with anything in the local domain or range vectors have
// been missed by only being detected on another processor.  Intersections
// can be missed if, say, some local entities exist only in the range vector
// and not the domain vector.  This flag will add another round of parallel
// communication that makes sure all local interactions that involve an
// off-processor entity are shared with that remote processor.  Results are
// sorted if this flag is set.
//
// Search performance will generally be better if the global domain vector is
// larger than the range vector.  The "autoSwapDomainAndRange" argument
// controls if a global input length comparison is made, followed by
// potentially swapping the domain and range vectors.  The final results should
// be independent of this flag.
//
template <typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
void coarse_search(std::vector<std::pair<DomainBoxType, DomainIdentProcType>> const & domain,
                   std::vector<std::pair<RangeBoxType, RangeIdentProcType>> const & range,
                   SearchMethod method,
                   stk::ParallelMachine comm,
                   std::vector<std::pair<DomainIdentProcType, RangeIdentProcType>>& intersections,
                   bool enforceSearchResultSymmetry = true,
                   bool autoSwapDomainAndRange = true,
                   bool sortSearchResults = false)
{
  switch (method) {
    case ARBORX: {
#ifdef STK_HAS_ARBORX
      coarse_search_arborx(domain, range, comm, intersections, enforceSearchResultSymmetry, sortSearchResults);
#else
      STK_ThrowErrorMsg("STK(stk_search) was not configured with ARBORX enabled. Please use KDTREE or MORTON_LBVH.");
#endif
      break;
    }
    case KDTREE: {
      if (autoSwapDomainAndRange) {
        coarse_search_kdtree_driver(domain, range, comm, intersections, enforceSearchResultSymmetry, sortSearchResults);
      }
      else {
        coarse_search_kdtree(domain, range, comm, intersections, enforceSearchResultSymmetry, sortSearchResults);
      }
      break;
    }
    case MORTON_LBVH: {
      coarse_search_morton_lbvh(domain, range, comm, intersections, enforceSearchResultSymmetry, sortSearchResults);
      break;
    }
    default: {
      STK_ThrowErrorMsg("Unsupported coarse_search method supplied. Choices are: KDTREE, MORTON_LBVH, or ARBORX.");
    }
  }
}

template <typename DomainView, typename RangeView, typename ResultView, typename ExecutionSpace = typename DomainView::execution_space>
void coarse_search(DomainView const & domain,
                   RangeView const & range,
                   SearchMethod method,
                   stk::ParallelMachine comm,
                   ResultView& intersections,
                   ExecutionSpace const& execSpace = ExecutionSpace{},
                   bool enforceSearchResultSymmetry = true,
                   bool autoSwapDomainAndRange = true,
                   bool sortSearchResults = false)
{
  check_coarse_search_types_parallel<DomainView, RangeView, ResultView, ExecutionSpace>();
  Kokkos::Profiling::pushRegion("STK coarse search with Views");

  switch (method) {
    case ARBORX: {
#ifdef STK_HAS_ARBORX
      coarse_search_arborx(domain, range, comm, intersections, execSpace, enforceSearchResultSymmetry, sortSearchResults);
#else
      STK_ThrowErrorMsg("STK(stk_search) was not configured with ARBORX enabled. Please use KDTREE or MORTON_LBVH.");
#endif
      break;
    }
    case KDTREE: {
      STK_ThrowErrorMsg("The KDTREE search method is not supported on GPUs.  Please use MORTON_LBVH or ARBORX instead.");
      break;
    }
    case MORTON_LBVH: {
      coarse_search_morton_lbvh(domain, range, comm, intersections, execSpace, enforceSearchResultSymmetry, sortSearchResults);
      break;
    }
    default: {
      STK_ThrowErrorMsg("Unsupported coarse_search method supplied. Choices are: KDTREE, MORTON_LBVH, or ARBORX.");
    }
  }

  Kokkos::Profiling::popRegion();
}

} // namespace stk::search

#endif // stk_search_CoarseSearch_hpp
