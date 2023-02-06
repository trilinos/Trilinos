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

#ifndef _KOKKOSGRAPH_BFS_IMPL_HPP
#define _KOKKOSGRAPH_BFS_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosKernels_Utils.hpp"
#include <vector>
#include <algorithm>

namespace KokkosGraph {
namespace Experimental {
namespace Impl {

template <typename rowmap_t, typename entries_t, typename lno_view_t>
struct SerialRCM {
  using size_type       = typename rowmap_t::non_const_value_type;
  using lno_t           = typename entries_t::non_const_value_type;
  using host_rowmap_t   = Kokkos::View<size_type*, Kokkos::HostSpace>;
  using host_lno_view_t = Kokkos::View<lno_t*, Kokkos::HostSpace>;

  lno_t numVerts;
  host_rowmap_t rowmap;
  host_lno_view_t entries;

  SerialRCM(const rowmap_t& rowmap_, const entries_t& entries_)
      : numVerts(rowmap_.extent(0) - 1),
        rowmap(Kokkos::view_alloc(Kokkos::WithoutInitializing, "HostRowmap"),
               rowmap_.extent(0)),
        entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "HostEntries"),
                entries_.extent(0)) {
    Kokkos::deep_copy(rowmap, rowmap_);
    Kokkos::deep_copy(entries, entries_);
  }

  lno_t findPseudoPeripheral() {
    // Choose vertex with smallest degree
    lno_t periph    = -1;
    lno_t periphDeg = numVerts;
    for (lno_t i = 0; i < numVerts; i++) {
      lno_t deg = rowmap(i + 1) - rowmap(i);
      if (deg < periphDeg) {
        periph    = i;
        periphDeg = deg;
        if (deg == 0) break;
      }
    }
    return periph;
  }

  lno_view_t rcm() {
    lno_t start = findPseudoPeripheral();
    host_lno_view_t q(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Queue"),
                      numVerts);
    host_lno_view_t label(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "Permutation"),
        numVerts);
    for (lno_t i = 0; i < numVerts; i++) label(i) = -1;
    lno_t qhead  = 0;
    lno_t qtail  = 0;
    label(start) = qtail;
    q(qtail++)   = start;
    std::vector<lno_t> neighbors;
    lno_t outerQueue = 0;
    while (true) {
      lno_t v = q(qhead++);
      neighbors.clear();
      for (size_type j = rowmap(v); j < rowmap(v + 1); j++) {
        lno_t nei = entries(j);
        if (nei == v || nei >= numVerts) continue;
        if (label(nei) == -1) {
          neighbors.push_back(nei);
        }
      }
      std::sort(neighbors.begin(), neighbors.end(),
                [&](lno_t n1, lno_t n2) -> bool {
                  // return true if n1 has a lower degree than n2
                  return (rowmap(n1 + 1) - rowmap(n1)) <
                         (rowmap(n2 + 1) - rowmap(n2));
                });
      // label and enqueue all unlabeled neighbors
      for (lno_t nei : neighbors) {
        label(nei) = qtail;
        q(qtail++) = nei;
      }
      if (qtail == numVerts) {
        // have labeled all vertices
        break;
      } else if (qhead == qtail) {
        // have exhausted this connected component, but others remain unlabeled
        while (label(outerQueue) != -1) outerQueue++;
        label(outerQueue) = qtail;
        q(qtail++)        = outerQueue;
      }
    }
    lno_view_t labelOut(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "RCM Permutation"),
        numVerts);
    // reverse the labels
    for (lno_t i = 0; i < numVerts; i++) label(i) = numVerts - label(i) - 1;
    Kokkos::deep_copy(labelOut, label);
    return labelOut;
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosGraph
#endif
