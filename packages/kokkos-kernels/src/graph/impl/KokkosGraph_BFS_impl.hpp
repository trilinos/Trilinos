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
