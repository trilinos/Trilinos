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
#include "KokkosGraph_wiki_9pt_stencil.hpp"
#include "KokkosGraph_RCM.hpp"

template <typename rowmap_t, typename entries_t, typename labels_t>
void printReorderedMatrix(const rowmap_t& rowmapIn, const entries_t& entriesIn, const labels_t& invPermIn) {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  auto rowmap     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmapIn);
  auto entries    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entriesIn);
  auto invPerm    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), invPermIn);
  lno_t numVerts  = rowmap.extent(0) - 1;
  decltype(invPerm) perm(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Perm"), numVerts);
  for (lno_t i = 0; i < numVerts; i++) perm(invPerm(i)) = i;
  std::vector<lno_t> neighbors;
  for (lno_t i = 0; i < numVerts; i++) {
    lno_t origRow = perm(i);
    neighbors.clear();
    for (size_type j = rowmap(origRow); j < rowmap(origRow + 1); j++) {
      lno_t origNei = entries(j);
      lno_t nei     = invPerm(origNei);
      neighbors.push_back(nei);
    }
    std::sort(neighbors.begin(), neighbors.end());
    size_t it = 0;
    for (lno_t j = 0; j < numVerts; j++) {
      if (it < neighbors.size() && j == neighbors[it]) {
        std::cout << '*';
        it++;
      } else
        std::cout << ' ';
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

int main() {
  Kokkos::initialize();
  {
    using GraphDemo::numVertices;
    GraphDemo::setGridDimensions(6, 6);
    RowmapType rowmapDevice;
    ColindsType colindsDevice;
    // Make the graph smaller so the matrix can be printed easily
    // Step 1: Generate the graph on host, allocate space on device, and copy.
    // See function "generate9pt" below.
    GraphDemo::generate9pt(rowmapDevice, colindsDevice);
    // Step 2: Run RCM and print the reordered matrix
    {
      auto rcmDevice =
          KokkosGraph::Experimental::graph_rcm<ExecSpace, RowmapType, ColindsType>(rowmapDevice, colindsDevice);
      std::cout << "Graph reordered by reverse Cuthill-McKee:\n";
      printReorderedMatrix(rowmapDevice, colindsDevice, rcmDevice);
    }
  }
  Kokkos::finalize();
  return 0;
}
