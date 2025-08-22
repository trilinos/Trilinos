//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSGRAPH_RCM_HPP
#define KOKKOSGRAPH_RCM_HPP

#include "KokkosGraph_BFS_impl.hpp"

namespace KokkosGraph {
namespace Experimental {

// Compute the reverse Cuthill-McKee ordering of a graph.
// The graph must be symmetric, but it may have any number of connected
// components. This function returns a list of vertices in RCM order.

template <typename device_t, typename rowmap_t, typename colinds_t,
          typename labels_t = typename colinds_t::non_const_type>
labels_t graph_rcm(const rowmap_t& rowmap, const colinds_t& colinds) {
  using lno_t = typename colinds_t::non_const_value_type;
  if (rowmap.extent(0) <= 2) {
    // there are 0 or 1 vertices - return trivial ordering
    lno_t numVerts = rowmap.extent(0);
    if (numVerts) numVerts--;
    return labels_t("RCM Labels", numVerts);
  }
  Impl::SerialRCM<rowmap_t, colinds_t, labels_t> algo(rowmap, colinds);
  return algo.rcm();
}

}  // namespace Experimental
}  // namespace KokkosGraph

#endif
