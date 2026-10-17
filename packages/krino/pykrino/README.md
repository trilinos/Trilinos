# pykrino

Python bindings and utilities for Krino: bounding box and file meshes, refinement (uniform, or to a size or metric
given as a Python function), mesh smoothing objectives and optimizers, analytic and STL surfaces, and log capture.

## Building

pykrino needs Python 3, pybind11, numpy and [mpi4py](https://mpi4py.readthedocs.io), with mpi4py built against the
same MPI as Krino.

- **In a Trilinos build with Krino:** configure with `-DKrino_ENABLE_PyKrino:BOOL=ON` and
  `-DCMAKE_POSITION_INDEPENDENT_CODE=ON` (and `-Dpybind11_DIR=...` if CMake cannot find pybind11), as
  `krino/cmake_install_test/run_cmake_krino` does with `ENABLE_PYKRINO=ON`.  The module is built in
  `<build>/packages/krino/pykrino/src` and installed to `<prefix>/lib/pythonX.Y/site-packages/pykrino`; the tests run
  with `ctest -R PYTHON`.
- **Against an installed Krino** (e.g. the `py-krino` Spack package): `pykrino/SierraCmakeCode.cmake` builds the module
  with `find_package(krino)`.

Put the directory containing the `pykrino` package on `PYTHONPATH`, e.g.
`export PYTHONPATH=<build>/packages/krino/pykrino/src`.

## Setup, MPI and parallel runs

Importing pykrino imports `mpi4py.MPI`, which initializes MPI.  The first mesh then sets up Kokkos and pykrino's
communicator (a duplicate of `MPI.COMM_WORLD`), and `pykrino.finalize()` runs at exit.  So no setup is needed:

```python
import pykrino

mesh = pykrino.BoundingBoxMesh("TET_4", [0, 0, 0], [1, 1, 1], 0.25, 0)
mesh.populate_mesh()
print(pykrino.parallel_rank(), pykrino.parallel_size())
```

Run in parallel with `mpirun -np 4 python3 script.py`.  To use another communicator, call
`pykrino.initialize(comm)` (an mpi4py communicator) before creating the first mesh.  `pykrino.EnvData(comm)` from
older scripts still works and is equivalent to `pykrino.initialize(comm)`.

## Documentation

- `examples/QUICK_START.md`: the main classes and how to use them.
- `examples/README.md`: the example scripts.
- `test/`: unit tests, which also show usage.
