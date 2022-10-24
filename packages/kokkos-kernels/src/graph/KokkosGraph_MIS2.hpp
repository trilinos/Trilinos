/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Brian Kelley (bmkelle@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSGRAPH_DISTANCE2_MIS_HPP
#define _KOKKOSGRAPH_DISTANCE2_MIS_HPP

#include "KokkosGraph_Distance2MIS_impl.hpp"

namespace KokkosGraph {

enum MIS2_Algorithm { MIS2_QUALITY, MIS2_FAST };

// Compute a distance-2 maximal independent set, given a symmetric CRS graph.
// Returns a list of the vertices in the set.
//
// Column indices >= num_verts are ignored.

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename lno_view_t = typename colinds_t::non_const_type>
lno_view_t graph_d2_mis(const rowmap_t& rowmap, const colinds_t& colinds,
                        MIS2_Algorithm algo = MIS2_FAST) {
  if (rowmap.extent(0) <= 1) {
    // zero vertices means the MIS is empty.
    return lno_view_t();
  }
  switch (algo) {
    case MIS2_QUALITY: {
      Impl::D2_MIS_FixedPriority<device_t, rowmap_t, colinds_t, lno_view_t> mis(
          rowmap, colinds);
      return mis.compute();
    }
    case MIS2_FAST: {
      Impl::D2_MIS_RandomPriority<device_t, rowmap_t, colinds_t, lno_view_t>
          mis(rowmap, colinds);
      return mis.compute();
    }
  }
  throw std::invalid_argument("graph_d2_mis: invalid algorithm");
}

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename labels_t = typename colinds_t::non_const_type>
labels_t graph_mis2_coarsen(
    const rowmap_t& rowmap, const colinds_t& colinds,
    typename colinds_t::non_const_value_type& numClusters) {
  if (rowmap.extent(0) <= 1) {
    // there are no vertices to label
    numClusters = 0;
    return labels_t();
  }
  Impl::D2_MIS_Aggregation<device_t, rowmap_t, colinds_t, labels_t> aggregation(
      rowmap, colinds);
  aggregation.compute(false);
  numClusters = aggregation.numAggs;
  return aggregation.labels;
}

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename labels_t = typename colinds_t::non_const_type>
labels_t graph_mis2_aggregate(
    const rowmap_t& rowmap, const colinds_t& colinds,
    typename colinds_t::non_const_value_type& numAggregates) {
  if (rowmap.extent(0) <= 1) {
    // there are no vertices to label
    numAggregates = 0;
    return labels_t();
  }
  Impl::D2_MIS_Aggregation<device_t, rowmap_t, colinds_t, labels_t> aggregation(
      rowmap, colinds);
  aggregation.compute(true);
  numAggregates = aggregation.numAggs;
  return aggregation.labels;
}

inline const char* mis2_algorithm_name(MIS2_Algorithm algo) {
  switch (algo) {
    case MIS2_QUALITY: return "MIS2_QUALITY";
    case MIS2_FAST: return "MIS2_FAST";
  }
  return "*** Invalid MIS2 algo enum value.\n";
}

}  // end namespace KokkosGraph

// For backward compatibility
namespace KokkosGraph {
namespace Experimental {

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename lno_view_t = typename colinds_t::non_const_type>
[[deprecated]] lno_view_t graph_d2_mis(const rowmap_t& rowmap,
                                       const colinds_t& colinds,
                                       MIS2_Algorithm algo = MIS2_FAST) {
  return KokkosGraph::graph_d2_mis<device_t, rowmap_t, colinds_t, lno_view_t>(
      rowmap, colinds, algo);
}

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename labels_t = typename colinds_t::non_const_type>
[[deprecated]] labels_t graph_mis2_coarsen(
    const rowmap_t& rowmap, const colinds_t& colinds,
    typename colinds_t::non_const_value_type& numClusters) {
  return KokkosGraph::graph_mis2_coarsen<device_t, rowmap_t, colinds_t,
                                         labels_t>(rowmap, colinds,
                                                   numClusters);
}

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename labels_t = typename colinds_t::non_const_type>
[[deprecated]] labels_t graph_mis2_aggregate(
    const rowmap_t& rowmap, const colinds_t& colinds,
    typename colinds_t::non_const_value_type& numAggregates) {
  return KokkosGraph::graph_mis2_aggregate<device_t, rowmap_t, colinds_t,
                                           labels_t>(rowmap, colinds,
                                                     numAggregates);
}

[[deprecated]] inline const char* mis2_algorithm_name(MIS2_Algorithm algo) {
  return KokkosGraph::mis2_algorithm_name(algo);
}

}  // namespace Experimental
}  // namespace KokkosGraph

#endif
