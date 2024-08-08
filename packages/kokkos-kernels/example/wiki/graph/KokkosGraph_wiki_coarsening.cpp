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
#include "KokkosGraph_MIS2.hpp"

int main() {
  Kokkos::initialize();
  {
    using GraphDemo::numVertices;
    RowmapType rowmapDevice;
    ColindsType colindsDevice;
    // Step 1: Generate the graph on host, allocate space on device, and copy.
    // See function "generate9pt" below.
    GraphDemo::generate9pt(rowmapDevice, colindsDevice);
    // Step 2: Run MIS-2 based coarsening and print the result
    {
      std::cout << "Coarsened vertex labels:\n";
      Ordinal numClusters = 0;
      auto labels = KokkosGraph::graph_mis2_aggregate<ExecSpace, RowmapType, ColindsType>(rowmapDevice, colindsDevice,
                                                                                          numClusters);
      // coarsening labels can be printed in the same way as colors
      GraphDemo::printColoring(labels, numClusters);
      putchar('\n');
    }
  }
  Kokkos::finalize();
  return 0;
}
