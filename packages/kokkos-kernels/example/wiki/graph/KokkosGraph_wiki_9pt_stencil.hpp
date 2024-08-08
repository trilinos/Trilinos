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
#ifndef WIKI_9PT_STENCIL_H
#define WIKI_9PT_STENCIL_H

#include "Kokkos_Core.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Handle.hpp"
#include <vector>
#include <set>
#include <cstdio>
#include <cmath>
#include <sstream>

using Ordinal     = default_lno_t;
using Offset      = default_size_type;
using Layout      = default_layout;
using ExecSpace   = Kokkos::DefaultExecutionSpace;
using DeviceSpace = typename ExecSpace::memory_space;
using Kokkos::HostSpace;
using RowmapType  = Kokkos::View<Offset*, DeviceSpace>;
using ColindsType = Kokkos::View<Ordinal*, DeviceSpace>;
using Handle = KokkosKernels::Experimental::KokkosKernelsHandle<Offset, Ordinal, default_scalar, ExecSpace, DeviceSpace,
                                                                DeviceSpace>;

namespace GraphDemo {
Ordinal gridX       = 15;
Ordinal gridY       = 25;
Ordinal numVertices = gridX * gridY;

void setGridDimensions(Ordinal newX, Ordinal newY) {
  gridX       = newX;
  gridY       = newY;
  numVertices = gridX * gridY;
}

// Helper to get the vertex ID given grid coordinates
Ordinal getVertexID(Ordinal x, Ordinal y) { return y * gridX + x; }

// Inverse of getVertexID
void getVertexPos(Ordinal vert, Ordinal& x, Ordinal& y) {
  x = vert % gridX;
  y = vert / gridX;
}

// Helper to print out colors in the shape of the grid
template <typename ColorView>
void printColoring(ColorView colors, Ordinal numColors) {
  // Read colors on host
  auto colorsHost = Kokkos::create_mirror_view_and_copy(HostSpace(), colors);
  int numDigits   = ceil(log10(numColors + 1));
  // Print out the grid, with columns aligned and at least one space between
  // numbers
  std::ostringstream numFmtStream;
  numFmtStream << '%' << numDigits + 1 << 'd';
  std::string numFmt = numFmtStream.str();
  for (Ordinal y = 0; y < gridY; y++) {
    for (Ordinal x = 0; x < gridX; x++) {
      Ordinal vertex = getVertexID(x, y);
      int color      = colorsHost(vertex);
      printf(numFmt.c_str(), color);
    }
    putchar('\n');
  }
}

template <typename MISView>
void printMIS(MISView misList) {
  // Read colors on host
  auto misHost = Kokkos::create_mirror_view_and_copy(HostSpace(), misList);
  std::set<Ordinal> mis;
  for (Offset i = 0; i < (Offset)misList.extent(0); i++) mis.insert(misHost(i));
  for (Ordinal y = 0; y < gridY; y++) {
    for (Ordinal x = 0; x < gridX; x++) {
      Ordinal vertex = getVertexID(x, y);
      if (mis.find(vertex) == mis.end())
        printf(". ");
      else
        printf("# ");
    }
    putchar('\n');
  }
}

// Build the graph on host, allocate these views on device and copy the graph to
// them. Both rowmapDevice and colindsDevice are output parameters and should
// default-initialized (empty) on input.
void generate9pt(RowmapType& rowmapDevice, ColindsType& colindsDevice) {
  // Generate the graph on host (use std::vector to not need to know
  // how many entries ahead of time)
  std::vector<Offset> rowmap(numVertices + 1);
  std::vector<Ordinal> colinds;
  rowmap[0] = 0;
  for (Ordinal vert = 0; vert < numVertices; vert++) {
    Ordinal x, y;
    getVertexPos(vert, x, y);
    // Loop over the neighbors in a 3x3 region
    for (Ordinal ny = y - 1; ny <= y + 1; ny++) {
      for (Ordinal nx = x - 1; nx <= x + 1; nx++) {
        // exclude the edge to self
        if (nx == x && ny == y) continue;
        // exclude vertices that would be outside the grid
        if (nx < 0 || nx >= gridX || ny < 0 || ny >= gridY) continue;
        // add the neighbor to colinds, forming an edge
        colinds.push_back(getVertexID(nx, ny));
      }
    }
    // mark where the current row ends
    rowmap[vert + 1] = colinds.size();
  }
  Offset numEdges = colinds.size();
  // Now that the graph is formed, copy rowmap and colinds to Kokkos::Views in
  // device memory The nonowning host views just alias the std::vectors.
  Kokkos::View<Offset*, HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> rowmapHost(rowmap.data(), numVertices + 1);
  Kokkos::View<Ordinal*, HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> colindsHost(colinds.data(), numEdges);
  // Allocate owning views on device with the correct size.
  rowmapDevice  = RowmapType("Rowmap", numVertices + 1);
  colindsDevice = ColindsType("Colinds", numEdges);
  // Copy the graph from host to device
  Kokkos::deep_copy(rowmapDevice, rowmapHost);
  Kokkos::deep_copy(colindsDevice, colindsHost);
}
}  // namespace GraphDemo

#endif
