# PyKrino Examples

Example scripts for the pykrino Python bindings.  See `../README.md` for building pykrino and putting it on
`PYTHONPATH`, and `QUICK_START.md` for the classes used here.  All run serially or in parallel
(`mpirun -np N python3 <script>`); no MPI or Kokkos setup is needed in the scripts.

| script | shows | needs |
|---|---|---|
| `mesh_smoothing_example.py` | smoothing a mesh from a file with `MeanRatioElementObjective` and the LBFGS, SteepestDescent and ROL optimizers, with and without a gradient filter | a mesh file `example_mesh.e`; prints the API pattern otherwise |
| `kinematic_mesh_smoothing_example.py` | the Seth-Hill and Modified Neo-Hookean kinematic element objectives and `SethHillParams` | nothing (constructs objectives; shows the mesh usage pattern) |
| `gradient_filter_example.py` | creating gradient filters and attaching them to mesh objectives | nothing |
| `log_redirecter_example.py` | capturing krino log output with `LogRedirecter` | nothing |

The tests in `../test/` also show how each class is used.
