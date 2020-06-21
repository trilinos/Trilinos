Graph Examples
==============
This directory contains Kokkos-Kernels graph routine examples.


Distance 2 Graph Coloring (D2GC)
--------------------------------
The distance-2 graph coloring example is contained in the file
`KokkosKernels_Example_Distance2GraphColor.cpp`.

This example demonstrates the use of the Distance 2 Graph Coloring code and
some of the associated helper functions for it.  This example is a simplified
version of the test that can be found in the `perf_test` directories.

The following source files contain the D2GC implementation:
- `src/graph/KokkosGraph_Distance2Color.hpp`
- `src/graph/KokkosGraph_Distance2ColorHandle.hpp`
- `src/graph/impl/KokkosGraph_Distance2Color_impl.hpp`
- `src/graph/impl/KokkosGraph_Distance2Color_MatrixSquared_impl.hpp`

The following methods from D2GC are demonstrated by this example:
- `graph_compute_distance2_color()` - Executes the D2GC routine.
- `graph_verify_distance2_color()`  - Verify the correctness of the coloring produced by D2GC.
- `graph_print_distance2_color_histogram()` - Prints out a histogram of the colors produced.
- `dump_graphviz()` - Save the graph with coloring labels to a GraphViz `.dot` file.

This example allows the user to select from the following variants of the D2GC algorithm:
- `COLORING_D2_MATRIX_SQUARED`: Matrix-squared + distance-1 method.
- `COLORING_D2_SERIAL`:         Serial Distance-2 coloring (must be used with 'serial' mode).
- `COLORING_D2_VB`:             Vertex-based method using an array of bools for the forbidden array (default).
- `COLORING_D2_VB_BIT`:         VB method using a bit vector for the forbidden array.
- `COLORING_D2_VB_BIT_EF`:      `VB_BIT` method with some additional edge filtering.

### Building
Examples are enabled by the makefile generated from `scripts/generate_makefile.bash` by passing
the option `build-all` or `build-example` to the make command.

### Usage
Once the examples are built, the binaries will be located in the `example` subdirectory of your 
build directory. The D2GC example application will be in `KokkosKernels_Example_Distance2GraphColor.exe`.
Executing the application with no parameters will print out the usage options.

```
Usage:
  ./KokkosKernels_Example_Distance2GraphColor.exe [parameters]

Parameters:
  Parallelism (select one of the following):
      --serial <N>        Execute serially.
      --threads <N>       Use N posix threads.
      --openmp <N>        Use OpenMP with N threads.
      --cuda              Use CUDA

  Required Parameters:
      --amtx <filename>   Input file in Matrix Market format (.mtx).

      --algorithm <algorithm_name>   Set the algorithm to use.  Allowable values are:
                 COLORING_D2_MATRIX_SQUARED  - Matrix-squared + Distance-1 method.
                 COLORING_D2_SERIAL          - Serial algorithm (must use with 'serial' mode)
                 COLORING_D2_VB              - Vertex Based method using boolean forbidden array (Default).
                 COLORING_D2_VB_BIT          - VB with Bitvector Forbidden Array
                 COLORING_D2_VB_BIT_EF       - VB_BIT with Edge Filtering

  Optional Parameters:
      --output-histogram              Print out a histogram of the colors.
      --output-graphviz               Write the output to a graphviz file (G.dot).
                                      Note: Vertices with color 0 will be filled in and colored
      --output-graphviz-vert-max <N>  Upper limit of vertices in G to allow graphviz output. Default=1500.
                                      Requires --output-graphviz to also be enabled.
      --validate                      Check that the coloring is a valid distance-2 graph coloring
      --verbose-level <N>             Set verbosity level [0..5] where N > 0 means print verbose messags.
                                      Default: 0
      --help                          Print out command line help.
```


### Example Use

Running the example with the following parameters:
```
./KokkosKernels_Example_Distance2GraphColor.exe \
        --openmp 8 \
        --algorithm COLORING_D2_VB_BIT \
        --output-histogram \
        --validate \
        --amtx ${HOME}/Data/MatrixMarket/vanHeukelum/cage11.mtx
```
runs the Vertex-Based + BitVector based forbidden-array variation using 8 OpenMP threads
on the [vanHeukelum][1] [Cage11][2] graph. The coloring will be checked for correctness and a histogram 
summary will be printed out. This results in the following output:
```
Run Graph Color D2 (COLORING_D2_VB_BIT)

Distance-2 Graph Coloring is VALID

Distance-2 Color Histogram:
1767 1737 1774 1748 1619 1610 1596 1494 1467 1395 1358 1288 1272 1188 1161 1082 1050 1016 962 ... ... ... 30 26 22 21 16 18 13 13 12 8 8 7 6 5 5 4 3 2 2 2

Summary
-------
    KExecSName     : OpenMP
    Filename       : cage11.mtx
    Num Verts      : 39082
    Num Edges      : 559722
    Concurrency    : 8
    Algorithm      : COLORING_D2_VB_BIT
Coloring Stats
    Num colors     : 76
    Num Phases     : 2
    Validation     : VALID
```




[1]: https://sparse.tamu.edu/vanHeukelum
[2]: https://sparse.tamu.edu/vanHeukelum/cage11


