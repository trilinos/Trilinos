Finite Element Mesh Assembly
----------------------------


Files
------

### Common Files
- `fem_assembly_typedefs.hpp`
  - Type definitions for the FEM Assembly examples.
- `fem_assembly_Element.hpp`
  - Provides functions that provide a rough approximation of a 2D finite-difference stencil.
- `fem_assembly_MeshDatabase.hpp`
  - Provides the core 2D mesh database for this set of examples.
- `fem_assembly_commandLineOpts.hpp`
  - Provides helpers for processing command line options.
- `fem_assembly_main`
  - Contains the main() function for the application.

### Finite Element Mesh Assembly Example Sources

- `fem_assembly_TotalElementLoop.hpp`
  - Each process contains information from both its owned elements and its ghost elements so mesh construction
    is possible without requiring communication.
- `fem_assembly_InsertGlobalIndices_FE.hpp`
  - Constructs an FECrsGraph by looping over the _owned_ elements and inserting the representations of the
    connectivity of each element into the graph using their global ids. Uses a Kokkos kernel locally
    and lets Tpetra handle the communication.


Running The Example
-------------------

### Command Line Options
- `--help` Print help and exit.
- `--num-elements-x` Number of elements to generate in the x-axis of the 2D grid.
- `--num-elements-y` Number of elements to generate in the y-axis of the 2D grid.
- `--verbose` Print out extra information to STDOUT.
- `--timing`  Print out timing information at the end of execution.
- `--with-insert-global-indices-fe` Execute the Insert FECrsMatrix Global Indices FEM Assembly kernel.
- `--with-total-element-loop`  Execute the Total Element Loop FEM Assembly kernel.

