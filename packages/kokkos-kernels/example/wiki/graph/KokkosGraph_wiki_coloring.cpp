#include <vector>
#include <cstdio>
#include <cmath>
#include <sstream>
#include "Kokkos_Core.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance2Color.hpp"

//Greedy Graph Coloring
//  -Generate the graph for a rectangular grid, with a 9-point stencil
//   (each vertex is adjacent to the 8 vertices around it within 1 grid square)
//  -Run Distance-1 coloring (usual coloring: adjacent vertices must have different colors)
//    -Print out the colors of each vertex in a grid
//  -Run Distance-2 coloring, and print out the colors
//    -Different constraint: two vertices separated by a path of length 1 OR 2
//     must have different colors)

using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;
using ExecSpace = Kokkos::DefaultExecutionSpace;
using DeviceSpace = typename ExecSpace::memory_space;
using Kokkos::HostSpace;
using RowmapType = Kokkos::View<Offset*, DeviceSpace>;
using ColindsType = Kokkos::View<Ordinal*, DeviceSpace>;
using Handle  = KokkosKernels::Experimental::
  KokkosKernelsHandle<Offset, Ordinal, default_scalar, ExecSpace, DeviceSpace, DeviceSpace>;

namespace ColoringDemo
{
  constexpr Ordinal gridX = 15;
  constexpr Ordinal gridY = 25;
  constexpr Ordinal numVertices = gridX * gridY;

  //Helper to get the vertex ID given grid coordinates
  Ordinal getVertexID(Ordinal x, Ordinal y)
  {
    return y * gridX + x;
  }

  //Inverse of getVertexID
  void getVertexPos(Ordinal vert, Ordinal& x, Ordinal& y)
  {
    x = vert % gridX;
    y = vert / gridX;
  }

  //Helper to print out colors in the shape of the grid
  template<typename ColorView>
  void printColoring(ColorView colors, Ordinal numColors)
  {
    //Read colors on host
    auto colorsHost = Kokkos::create_mirror_view_and_copy(HostSpace(), colors);
    int numDigits = ceil(log10(numColors + 1));
    //Print out the grid, with columns aligned and at least one space between numbers
    std::ostringstream numFmtStream;
    numFmtStream << '%' << numDigits + 1 << 'd';
    std::string numFmt = numFmtStream.str();
    for(Ordinal y = 0; y < gridY; y++)
    {
      for(Ordinal x = 0; x < gridX; x++)
      {
        Ordinal vertex = getVertexID(x, y);
        int color = colorsHost(vertex);
        printf(numFmt.c_str(), color);
      }
      putchar('\n');
    }
  }

  //Build the graph on host, allocate these views on device and copy the graph to them.
  //Both rowmapDevice and colindsDevice are output parameters and should default-initialized (empty) on input.
  void generate9pt(RowmapType& rowmapDevice, ColindsType& colindsDevice)
  {
    //Generate the graph on host (use std::vector to not need to know
    //how many entries ahead of time)
    std::vector<Offset> rowmap(numVertices + 1);
    std::vector<Ordinal> colinds;
    rowmap[0] = 0;
    for(Ordinal vert = 0; vert < numVertices; vert++)
    {
      Ordinal x, y;
      getVertexPos(vert, x, y);
      //Loop over the neighbors in a 3x3 region
      for(Ordinal ny = y - 1; ny <= y + 1; ny++)
      {
        for(Ordinal nx = x - 1; nx <= x + 1; nx++)
        {
          //exclude the edge to self
          if(nx == x && ny == y)
            continue;
          //exclude vertices that would be outside the grid
          if(nx < 0 || nx >= gridX || ny < 0 || ny >= gridY)
            continue;
          //add the neighbor to colinds, forming an edge
          colinds.push_back(getVertexID(nx, ny));
        }
      }
      //mark where the current row ends
      rowmap[vert + 1] = colinds.size();
    }
    Offset numEdges = colinds.size();
    //Now that the graph is formed, copy rowmap and colinds to Kokkos::Views in device memory
    //The nonowning host views just alias the std::vectors.
    Kokkos::View<Offset*, HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> rowmapHost(rowmap.data(), numVertices + 1);
    Kokkos::View<Ordinal*, HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> colindsHost(colinds.data(), numEdges);
    //Allocate owning views on device with the correct size.
    rowmapDevice = RowmapType("Rowmap", numVertices + 1);
    colindsDevice = ColindsType("Colinds", numEdges);
    //Copy the graph from host to device
    Kokkos::deep_copy(rowmapDevice, rowmapHost);
    Kokkos::deep_copy(colindsDevice, colindsHost);
  }
}

int main(int argc, char* argv[])
{
  Kokkos::initialize();
  {
    using ColoringDemo::numVertices;
    RowmapType rowmapDevice;
    ColindsType colindsDevice;
    //Step 1: Generate the graph on host, allocate space on device, and copy.
    //See function "generate9pt" below.
    ColoringDemo::generate9pt(rowmapDevice, colindsDevice);
    //Step 2: Create handle and run distance-1 coloring.
    {
      Handle handle;
      //Use the default algorithm (chosen based on ExecSpace)
      handle.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
      //Run coloring (graph is square and symmetric)
      KokkosGraph::Experimental::graph_color(&handle, numVertices, numVertices, rowmapDevice, colindsDevice);
      //Get the colors array, and the number of colors used from the handle.
      auto colors = handle.get_graph_coloring_handle()->get_vertex_colors();
      Ordinal numColors = handle.get_graph_coloring_handle()->get_num_colors();
      printf("9-pt stencil: Distance-1 Colors (used %d):\n", (int) numColors);
      ColoringDemo::printColoring(colors, numColors);
      putchar('\n');
      //Clean up
      handle.destroy_graph_coloring_handle();
    }
    //Step 3: Create handle and run distance-2 coloring.
    {
      Handle handle;
      //Use the default algorithm (chosen based on ExecSpace)
      handle.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_DEFAULT);
      //Run coloring
      KokkosGraph::Experimental::graph_color_distance2(&handle, numVertices, rowmapDevice, colindsDevice);
      //Get the colors array, and the number of colors used from the handle.
      auto colors = handle.get_distance2_graph_coloring_handle()->get_vertex_colors();
      Ordinal numColors = handle.get_distance2_graph_coloring_handle()->get_num_colors();
      printf("9-pt stencil: Distance-2 Colors (used %d):\n", (int) numColors);
      ColoringDemo::printColoring(colors, numColors);
      putchar('\n');
      //Clean up
      handle.destroy_distance2_graph_coloring_handle();
    }
  }
  Kokkos::finalize();
  return 0;
}

