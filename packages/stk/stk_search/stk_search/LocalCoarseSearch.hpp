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

#ifndef LOCALCOARSESEARCH_HPP
#define LOCALCOARSESEARCH_HPP

#include "stk_util/stk_config.h"
#include "stk_search/SearchMethod.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/kdtree/LocalCoarseSearchKdTree.hpp"
#include "stk_search/morton_lbvh/LocalCoarseSearchMortonLBVH.hpp"
#ifdef STK_HAS_ARBORX
#include "stk_search/arborx/LocalCoarseSearchArborX.hpp"
#endif
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"

namespace stk::search {

template <typename DomainBoxType, typename DomainIdentType, typename RangeBoxType, typename RangeIdentType>
void local_coarse_search(
    std::vector<std::pair<DomainBoxType, DomainIdentType>> const & domain,
    std::vector<std::pair<RangeBoxType, RangeIdentType>> const & range,
    SearchMethod method,
    std::vector<std::pair<DomainIdentType, RangeIdentType>> & intersections,
    bool sortSearchResults = false)
{
  switch (method) {
    case ARBORX: {
#ifdef STK_HAS_ARBORX
      local_coarse_search_arborx(domain, range, intersections, sortSearchResults);
#else
      STK_ThrowErrorMsg("STK(stk_search) was not configured with ARBORX. Please use KDTREE or MORTON_LBVH.");
#endif
      break;
    }
    case KDTREE: {
      local_coarse_search_kdtree_driver(domain, range, intersections, sortSearchResults);
      break;
    }
    case MORTON_LBVH: {
      local_coarse_search_morton_lbvh(domain, range, intersections, sortSearchResults);
      break;
    }
    default: {
      STK_ThrowErrorMsg("Unsupported local_coarse_search method supplied. Choices are: KDTREE, MORTON_LVBH, or ARBORX.");
    }
  }

}


template <typename DomainView, typename RangeView, typename ResultView, typename ExecutionSpace = typename DomainView::execution_space>
std::shared_ptr<SearchData>
 local_coarse_search(
    DomainView const & domain,
    RangeView const & range,
    SearchMethod method,
    ResultView & intersections,
    ExecutionSpace const& execSpace = ExecutionSpace{},
    bool sortSearchResults = false,
    std::shared_ptr<SearchData> searchData = nullptr)
{
  check_coarse_search_types_local<DomainView, RangeView, ResultView, ExecutionSpace>();

  switch (method) {
    case ARBORX: {
#ifdef STK_HAS_ARBORX
      local_coarse_search_arborx(domain, range, intersections, execSpace, sortSearchResults);
#else
      STK_ThrowErrorMsg("STK(stk_search) was not configured with ARBORX. Please use KDTREE or MORTON_LBVH.");
#endif
      break;
    }
    case KDTREE: {
      STK_ThrowErrorMsg("The KDTREE search method is not supported on GPUs.  Please use MORTON_LBVH or ARBORX instead.");
      break;
    }
    case MORTON_LBVH: {
      searchData = local_coarse_search_morton_lbvh(domain, range, intersections, execSpace, sortSearchResults, searchData);
      break;
    }
    default: {
      STK_ThrowErrorMsg("Unsupported local_coarse_search method supplied. Choices are: KDTREE, MORTON_LVBH, or ARBORX.");
    }
  }

  return searchData;
}

}

#endif // LOCALCOARSESEARCH_HPP
