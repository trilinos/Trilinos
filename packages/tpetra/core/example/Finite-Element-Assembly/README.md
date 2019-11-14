Finite Element Mesh Assembly
----------------------------


Files
------

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
- `fem_assembly_TotalElementLoop.hpp`
  - StaticProfile version of the Total Element Loop example.

Running The Example
-------------------

### Command Line Options
- `--help` Print help and exit.
- `--num-elements-x` Number of elements to generate in the x-axis of the 2D grid.
- `--num-elements-y` Number of elements to generate in the y-axis of the 2D grid.
- `--verbose` Print out extra information to STDOUT.
- `--timing`  Print out timing information at the end of execution.

