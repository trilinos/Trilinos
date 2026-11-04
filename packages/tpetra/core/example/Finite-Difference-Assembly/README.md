# Finite Difference Assembly

This example constructs the graph of a finite-difference operator directly
from a stencil.  It uses a two-dimensional Cartesian grid with rank-blocked,
locally lexicographic global numbering and the five-point stencil

```text
                 (0, 1)
                    |
(-1, 0) -------- (0, 0) -------- (1, 0)
                    |
                 (0,-1)
```

The example passes the row map, local Cartesian grid dimensions, and stencil
offsets to `Tpetra::Experimental::createCrsGraphFromStencil`.  The routine
constructs the fill-complete `Tpetra::CrsGraph` directly in device-resident CSR
form.  Physical-boundary neighbors are omitted.

For MPI runs, ranks are arranged into the same two-dimensional Cartesian
process grid inferred by the graph-construction routine.  `nx` and `ny` are
global dimensions and must be divisible by the corresponding process-grid
dimensions.

The example complements `Finite-Element-Assembly`: finite-element graph
assembly inserts the connectivity contributed by each element, whereas
finite-difference graph assembly inserts the neighbors prescribed for each
grid point by a stencil.

Options:

- `--nx`: number of grid points in the x direction;
- `--ny`: number of grid points in the y direction;
- `--verbose`: describe the completed graph.
