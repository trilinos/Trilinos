# PyKrino Quick Start Guide

See `../README.md` for building pykrino and running in parallel.  No MPI or Kokkos setup is needed: importing pykrino
initializes MPI (through mpi4py), and the first mesh sets up Kokkos.

```python
import pykrino
```

## Meshes

```python
mesh = pykrino.MeshFromFile("your_mesh.e")          # optional: decomp_method, use_all_sides_for_shells, read_fields
# or a generated mesh of a box: element type TET_4, HEX_8, TRI_3 or QUAD_4
mesh = pykrino.BoundingBoxMesh("TET_4", [0, 0, 0], [1, 1, 1], 0.25)   # optional: pad_cells
mesh.set_mesh_structure_type("BCC")                 # optional, before populate_mesh: CUBIC, BCC, ...
mesh.populate_mesh()

mesh.write_mesh("out.e")                            # compose=False to write without composing child elements
```

## Refinement

`Refinement` adds fields to the mesh, so create it before `populate_mesh()`:

```python
mesh = pykrino.BoundingBoxMesh("TET_4", [0, 0, 0], [1, 1, 1], 0.25)
refinement = pykrino.Refinement(mesh)
mesh.populate_mesh()

refinement.do_uniform_refinement(1)
# or refine every edge longer than sqrt(2) in a size or metric given as a function of coordinates;
# each returns the number of edges still too long (0 if converged)
refinement.refine_edges_too_long_in_size(lambda x: 0.1 + 0.2 * x[0])
refinement.refine_edges_too_long_in_metric(lambda x: [1/0.05**2, 0, 0, 0, 1/0.2**2, 0, 0, 0, 1/0.2**2])

# smoothing uses every element of the mesh, so delete the refinement parents first
refinement.delete_parent_elements()
```

## Smoothing

A mesh objective combines an element objective with the mesh; `SmoothMesh` minimizes it with an optimizer.

```python
objective = pykrino.MeshObjective(mesh, pykrino.MeanRatioElementObjective())

# kinematic objectives, with an element objective SethHillElementObjective([SethHillParams(mu, kappa, m, n)])
# or ModifiedNeoHookeanElementObjective(), and a reference element size:
kinematic = pykrino.ModifiedNeoHookeanElementObjective()
objective = pykrino.KinematicConstantIsotropicMeshObjective(mesh, kinematic, 0.1)        # one size
objective = pykrino.KinematicIsotropicMeshObjective(mesh, kinematic, pykrino.SmoothReference.EqualVolume)
objective = pykrino.KinematicAnisotropicMeshObjective(mesh, kinematic)  # then set a metric, see below

# optional: restrict node motion, e.g. keep nodes on the faces of the cube [-1,1]^3 on those faces
objective.set_gradient_filter(pykrino.GradientFilter.create_cube_boundary_filter(-1.0, 1.0))

optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS)  # or SteepestDescent, ROL; optional stepTol, gradTol, maxIter
pykrino.SmoothMesh().smooth_mesh(mesh, optimizer, objective, print_diagnostics=True)
```

For `KinematicAnisotropicMeshObjective`, give the metric M at each node (spatial_dimension x spatial_dimension
values, row-major) in the order of `get_node_ids_and_coordinates`, e.g. from a metric function:

```python
_, coords = objective.get_node_ids_and_coordinates()
metrics = [m for i in range(0, len(coords), 3) for m in metric_function(coords[i:i + 3])]
objective.set_element_metric_map_data_from_nodal_metrics(metrics)
```

or from a nodal field of the mesh with `set_element_metric_map_data_from_nodal_metric_field(field_name)`.

## Surfaces

```python
sphere = pykrino.Sphere([0, 0, 0], 0.5)
surface = pykrino.STLSurface("part.stl")     # also Cylinder(point, axis, radius)
distances = sphere.signed_distance(points)   # points: numpy array of shape (N, 3)
bbox = surface.bounding_box(scale=1.2)       # numpy array of shape (2, 3)
```

## Capturing krino log output

```python
log = pykrino.LogRedirecter()   # captures krinolog output while it exists
...
text = log.get_log()
log.clear()
```

Only output written to `krinolog` is captured, and only one `LogRedirecter` should be active at a time.

## Examples

| script | shows |
|---|---|
| `mesh_smoothing_example.py` | smoothing a mesh from a file with different optimizers |
| `kinematic_mesh_smoothing_example.py` | the Seth-Hill and Modified Neo-Hookean kinematic objectives |
| `gradient_filter_example.py` | gradient filters |
| `log_redirecter_example.py` | capturing krino log output |
