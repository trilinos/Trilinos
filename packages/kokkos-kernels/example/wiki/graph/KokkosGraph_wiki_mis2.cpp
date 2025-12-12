// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
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
    // Step 2: Run distance-2 MIS and print the results, with three different
    // algorithms
    {
      // Run coloring
      auto misDevice = KokkosGraph::graph_d2_mis<ExecSpace, RowmapType, ColindsType>(rowmapDevice, colindsDevice,
                                                                                     KokkosGraph::MIS2_FAST);
      std::cout << "Distance-2 MIS, FAST algorithm: contains " << misDevice.extent(0) << " out of "
                << GraphDemo::numVertices << " vertices.\n";
      GraphDemo::printMIS(misDevice);
      putchar('\n');
      misDevice = KokkosGraph::graph_d2_mis<ExecSpace, RowmapType, ColindsType>(rowmapDevice, colindsDevice,
                                                                                KokkosGraph::MIS2_QUALITY);
      std::cout << "Distance-2 MIS, QUALITY algorithm: contains " << misDevice.extent(0) << " out of "
                << GraphDemo::numVertices << " vertices.\n";
      GraphDemo::printMIS(misDevice);
      putchar('\n');
    }
  }
  Kokkos::finalize();
  return 0;
}
