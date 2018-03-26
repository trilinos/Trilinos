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

### Finite Element Mesh Assembly Examples
- `fem_assembly_InsertGlobalIndices_DP.hpp`
  - Constructs the crsGraph by looping over the _owned_ elements and inserting the representations of the
    connectivity of each element into the graph using their global ids and lets Tpetra handle the exchanges.
  - Performance:  :star:
  - Ease of use:  :star: :star: :star:
  - DynamicProfile
- `fem_assembly_LocalElementLoop_DP.hpp`
  - This method distinguishes owned from overlapping nodes.
    - Owned nodes are those that _only_ touch elements owned by the same process.
    - Overlapping nodes are nodes that touch elements owned by a different process.
  - Performance:  :star: :star:
  - Ease of use:  :star:
  - DynamicProfile
- `fem_assembly_TotalElementLoop_DP.hpp`
  - Each process contains information from both its owned elements and its ghost elements so mesh construction
    is possible without requiring communication.
  - Performance:  :star: :star: :star:
  - Ease of use:  :star: :star:
  - DynamicProfile
- `fem_assembly_TotalElementLoop_SP.hpp`
  - StaticProfile version of the Total Element Loop example.


Running The Example
-------------------

### Command Line Options
- `--help` Print help and exit.
- `--num-elements-x` Number of elements to generate in the x-axis of the 2D grid.
- `--num-elements-y` Number of elements to generate in the y-axis of the 2D grid.
- `--verbose` Print out extra information to STDOUT.
- `--timing`  Print out timing information at the end of execution.
- `--with-insert-global-indices-dp` Execute the "Insert Global Indices" example using DynamicProfile.
- `--with-local-element-loop-dp` Execute the "Local Element Loop" example using DynamicProfile.
- `--with-total-element-loop-dp` Execute the "Total Element Loop" example using DynamicProfile.
- `--with-total-element-loop-sp` Execute the "Total Element Loop" example using StaticProfile.



