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
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance2Color.hpp"

// Greedy Graph Coloring
//  -Generate the graph for a rectangular grid, with a 9-point stencil
//   (each vertex is adjacent to the 8 vertices around it within 1 grid square)
//  -Run Distance-1 coloring (usual coloring: adjacent vertices must have
//  different colors)
//    -Print out the colors of each vertex in a grid
//  -Run Distance-2 coloring, and print out the colors
//    -Different constraint: two vertices separated by a path of length 1 OR 2
//     must have different colors)

int main() {
  Kokkos::initialize();
  {
    using GraphDemo::numVertices;
    RowmapType rowmapDevice;
    ColindsType colindsDevice;
    // Step 1: Generate the graph on host, allocate space on device, and copy.
    // See function "generate9pt" below.
    GraphDemo::generate9pt(rowmapDevice, colindsDevice);
    // Step 2: Create handle and run distance-1 coloring.
    {
      Handle handle;
      // Use the default algorithm (chosen based on ExecSpace)
      handle.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
      // Run coloring (graph is square and symmetric)
      KokkosGraph::Experimental::graph_color(&handle, numVertices, numVertices, rowmapDevice, colindsDevice);
      // Get the colors array, and the number of colors used from the handle.
      auto colors       = handle.get_graph_coloring_handle()->get_vertex_colors();
      Ordinal numColors = handle.get_graph_coloring_handle()->get_num_colors();
      printf("9-pt stencil: Distance-1 Colors (used %d):\n", (int)numColors);
      GraphDemo::printColoring(colors, numColors);
      putchar('\n');
      // Clean up
      handle.destroy_graph_coloring_handle();
    }
    // Step 3: Create handle and run distance-2 coloring.
    {
      Handle handle;
      // Use the default algorithm (chosen based on ExecSpace)
      handle.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_DEFAULT);
      // Run coloring
      KokkosGraph::Experimental::graph_color_distance2(&handle, numVertices, rowmapDevice, colindsDevice);
      // Get the colors array, and the number of colors used from the handle.
      auto colors       = handle.get_distance2_graph_coloring_handle()->get_vertex_colors();
      Ordinal numColors = handle.get_distance2_graph_coloring_handle()->get_num_colors();
      printf("9-pt stencil: Distance-2 Colors (used %d):\n", (int)numColors);
      GraphDemo::printColoring(colors, numColors);
      putchar('\n');
      // Clean up
      handle.destroy_distance2_graph_coloring_handle();
    }
  }
  Kokkos::finalize();
  return 0;
}
