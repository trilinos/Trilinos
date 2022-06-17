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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSGRAPH_EXPLICIT_COARSEN_HPP
#define KOKKOSGRAPH_EXPLICIT_COARSEN_HPP

#include "KokkosGraph_ExplicitCoarsening_impl.hpp"
#include "KokkosKernels_Sorting.hpp"

namespace KokkosGraph {
namespace Experimental {

// Given a CRS graph and coarse labels, produce a new CRS graph representing the
// coarsened graph. If A is nonsquare, entries in columns >= numVerts are
// discarded. The labels should be in the range [0, numCoarseVerts), and the
// output graph wil have numCoarseVerts.
//
// If compress, sort and merge entries in each row.
// An uncompressed graph will still work as input to some things like D1 graph
// coloring.

template <typename device_t, typename fine_rowmap_t, typename fine_entries_t,
          typename labels_t, typename coarse_rowmap_t,
          typename coarse_entries_t>
void graph_explicit_coarsen(
    const fine_rowmap_t& fineRowmap, const fine_entries_t& fineEntries,
    const labels_t& labels,
    typename fine_entries_t::non_const_value_type numCoarseVerts,
    coarse_rowmap_t& coarseRowmap, coarse_entries_t& coarseEntries,
    bool compress = true) {
  using size_type  = typename fine_rowmap_t::non_const_value_type;
  using lno_t      = typename fine_entries_t::non_const_value_type;
  using exec_space = typename device_t::execution_space;
  static_assert(
      std::is_same<lno_t,
                   typename coarse_entries_t::non_const_value_type>::value,
      "graph_explicit_coarsen: The coarse and fine entry Views have different "
      "value types.");
  KokkosGraph::Impl::ExplicitGraphCoarsening<
      lno_t, size_type, device_t, fine_rowmap_t, fine_entries_t, labels_t,
      coarse_rowmap_t, coarse_entries_t, coarse_entries_t>
      egc(fineRowmap, fineEntries, labels, numCoarseVerts);
  coarseRowmap  = egc.coarseRowmap;
  coarseEntries = egc.coarseEntries;
  if (compress) {
    coarse_rowmap_t mergedRowmap;
    coarse_entries_t mergedEntries;
    KokkosKernels::sort_and_merge_graph<exec_space, coarse_rowmap_t,
                                        coarse_entries_t>(
        coarseRowmap, coarseEntries, mergedRowmap, mergedEntries);
    coarseRowmap  = mergedRowmap;
    coarseEntries = mergedEntries;
  }
}

// Same as above, but also produce the map from coarse vertices to fine vertices
// (inverse map of labels)
template <typename device_t, typename fine_rowmap_t, typename fine_entries_t,
          typename labels_t, typename coarse_rowmap_t,
          typename coarse_entries_t, typename ordinal_view_t>
void graph_explicit_coarsen_with_inverse_map(
    const fine_rowmap_t& fineRowmap, const fine_entries_t& fineEntries,
    const labels_t& labels,
    typename fine_entries_t::non_const_value_type numCoarseVerts,
    coarse_rowmap_t& coarseRowmap, coarse_entries_t& coarseEntries,
    ordinal_view_t& inverseOffsets, ordinal_view_t& inverseLabels,
    bool compress = true) {
  using size_type  = typename fine_rowmap_t::non_const_value_type;
  using lno_t      = typename fine_entries_t::non_const_value_type;
  using exec_space = typename device_t::execution_space;
  static_assert(
      std::is_same<lno_t,
                   typename coarse_entries_t::non_const_value_type>::value,
      "graph_explicit_coarsen: The coarse and fine entry Views have different "
      "value types.");
  KokkosGraph::Impl::ExplicitGraphCoarsening<
      lno_t, size_type, device_t, fine_rowmap_t, fine_entries_t, labels_t,
      coarse_rowmap_t, coarse_entries_t, ordinal_view_t>
      egc(fineRowmap, fineEntries, labels, numCoarseVerts);
  coarseRowmap   = egc.coarseRowmap;
  coarseEntries  = egc.coarseEntries;
  inverseOffsets = egc.clusterOffsets;
  inverseLabels  = egc.clusterVerts;
  if (compress) {
    coarse_rowmap_t mergedRowmap;
    coarse_entries_t mergedEntries;
    KokkosKernels::sort_and_merge_graph<exec_space, coarse_rowmap_t,
                                        coarse_entries_t>(
        coarseRowmap, coarseEntries, mergedRowmap, mergedEntries);
    coarseRowmap  = mergedRowmap;
    coarseEntries = mergedEntries;
  }
}

}  // namespace Experimental
}  // namespace KokkosGraph

#endif
