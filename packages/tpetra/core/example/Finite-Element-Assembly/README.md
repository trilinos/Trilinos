Finite Element Mesh Assembly
----------------------------


Files
------
- `typedefs.hpp`
  - Type definitions for the FEM Assembly examples.
- `Element.hpp`
  - Provides functions that provide a rough approximation of a 2D finite-difference stencil.
- `MeshDatabase.hpp`
  - Provides the core 2D mesh database for this set of examples.
- `fem_assembly_InsertGlobalIndices_DP.cpp`
  - Constructs the crsGraph by looping over the _owned_ elements and inserting the representations of the
    connectivity of each element into the graph using their global ids and lets Tpetra handle the exchanges.
  - Performance:  :star:
  - Ease of use:  :star: :star: :star:
  - DynamicProfile
- `fem_assembly_LocalElementLoop_DP.cpp`
  - This method distinguishes owned from overlapping nodes.
    - Owned nodes are those that _only_ touch elements owned by the same process.
    - Overlapping nodes are nodes that touch elements owned by a different process.
  - Performance:  :star: :star:
  - Ease of use:  :star:
  - DynamicProfile
- `fem_assembly_TotalElementLoop_DP.cpp`
  - Each process contains information from both its owned elements and its ghost elements so mesh construction
    is possible without requiring communication.
  - Performance:  :star: :star: :star:
  - Ease of use:  :star: :star:
  - DynamicProfile
- `fem_assembly_TotalElementLoop_SP.cpp`
  - StaticProfile version of the Total Element Loop example.
