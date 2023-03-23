# Creating nonconformal interface mesh for Exodus files

This document describes how to create the nonconformal interface mesh
for meshes stored in Exodus files.

The overall procedure is
  1. Create nonconforming volume mesh in mesh generator of your choice
  2. Create sidesets for both sides of the nonconformal interface
    * Remember the names of the sideset for future step
  3. Export mesh in Exodus format
    * Ensure there are no mesh vertices in common between the two sides
      of the interface.
  4. Pass the name of the Exodus file, the names of corresponding sidesets
     (from Step 2), and the filenames the output will be written to to
     the executable `nc_generator`

Details for each step are given below.

The usage of `nc_generator` is:

```
  nc_generator input_fname sideset_names output_fname output_fname2
```

The `sideset_names` is a list of pairs of sideset names.  The
pairs are separated by semicolons, and the two names within each
pair are separated by a colon.  Leading and trailing whitespace is
removed from the sideset names (to allow for more readable formatting).

The process of computing the nonconformal interface may modify the
coordinates of the nodes in the sideset (to ensures the boundaries of
the two sidesets are coincident everywhere).
`output_fname` is the name of the file the input mesh will be written
to, with updated coordinates.
`output_fname2` is the name of a new Exodus file that contains only the
nonconformal interface mesh and the information relating it to the
sidesets in `output_fname`.

Examples:
```
  nc_generator input.exo foo:bar output.exo output_surface.exo
```
  will create a nonconformal interface between sidesets `foo` and `bar`.
  `output.exo` contains the same mesh as `input.exo` with updated
  coordinates`.  The nonconformal interface mesh is in `output_surface.exo`

```
  nc_generator input.exo "foo : bar" output.exo output_surface.exo
```

is equivalent (the whitespace will be removed).

To specify more than one sideset

```
  nc_generator input.exo foo:bar;abc:def output.exo output_surface.exo
```

will generate nonconformal interfaces between sidesets `foo` and `bar`,
as well as between `abc` and `def`.
An equivalent (but more readable) invocation is

```
  nc_generator input.exo foo:bar ; abc:def output.exo output_surface.exo
```

`nc_generator` will throw an error if the sidesets do not exist.
It does *not* verify that the nodes on sideset pairs lie the same surface
(because this is very difficult to do for arbitrary geometries).
It is recommended that the user visualize the `output.exo` in a
visualization program (ex. Paraview) and verify the newly created
interface mesh is located where they expect it do be (see below for
the name of the interface mesh).


## Steps 1 and 2

  At the end of this step, you should have a mesh with at least two
  blocks and one pair of coincident faces
  (one face of the first block should be
  coincident with one face of the second block).  We will assume the
  meshes on the coincident faces are nonconformal for generality,
  (although this is not required, this package has been tested and
  works on conformal meshes).

  There are (at least) two ways to get CAD-based mesh generators to
  generate nonconformal meshes:
    1. Create nonconformal geometry
    2. Create two geometries, mesh them independently, and combine the
       results (for example, using the `ejoin` program from SEACAS)

  Option 1 is CAD program specific.  In Solidworks, this can be
  done by creating two Extrusion features and unchecking the "Merge Result"
  option for the second one.  This will create two separate volume object
  in the CAD file, which the mesher can then mesh independently.

  Option 2 can be accomplished by creating two separate CAD models with
  coincident faces, or by creating a single CAD model and exporting
  some of the entities to one file and the remaining entities to another,
  if the CAD software supports this.

  Note that Option 1 is more general than Option 2 for the following
  reason: Option 2 requires the nonconformal interface to partition
  the geometry (ie. each geometric entity must exist on exactly one
  of the CAD models), otherwise the combined mesh will have overlapping
  volume meshes.  This can be worked around by disabling meshing of
  certain volumes in either of the CAD models, but this may create
  additional nonconformal faces.

  In either case, you must create two sidesets for each nonconformal
  faces, one for each side of the face.  You must also remember the
  names given to the sidesets.  While it might be more convenient to
  refer to sidesets by number rather than name, some mesh generators
  (Pointwise) do not output sidesets any recognizable order, so it
  is impossible to know what the sidesets numbers in the exported
  Exodus file will be.

## Step 3

  When exporting the mesh to an Exodus file, there is one important
  criteria that must be satisifed: there can be no mesh vertices
  in common between the faces on the nonconformal interface.
  This can occur with Option 1 (above) if the CAD geometry uses the
  same geometric vertices as part of both volumes (ie. the CAD model
  is non-manifold and two geometric faces are created with the same
  geometric vertices).

  One solution to this problem is to export
  the volume meshes to two separate files and then use `ejoin`
  to combine them into one file.  In Pointwise, you can select
  Blocks in List View and then do the usual `File -> Export CAE` and
  it will export only the selected blocks.

# Step 4

  This is the step where the nonconformal interface mesh is created.
  The information required is:

    * The name of the input Exodus file
    * For each nonconformal interface, a pair of names that identify the two sideset created in Step 2
    * A filename to write the updated version of the input Exodus file to
    * A filename to write a new Exodus file to (containing the nonconformal interface mesh)


  The second output mesh will contain one new volume block for each nonconformal
  interface, as well as global attributes giving the names of the sidesets
  on the first output mesh that form the interface.  The attribute names are

  ```
    ci7Hl3_interface$i_sideset_names
  ```

  where `$i` goes from 0 to the number of sideset pairs minus 1.  The values
  are

  ```
  $name1\t$name2
  ```

  where `$name1` is the name of the first sideset in the pair and `$name2`
  is the name of the second sideset in the pair.  Note that the number of
  sideset pairs is equal to the number of blocks in the file.


  The block name will be computed as:

  ```
    ci7Hl3_nonconformal_interface_$name1_and_$name2
  ```

  The part topology with `SHELL_TRI3`,
  and will contain a unique set of nodes (ie. the nodes are not shared
  with the volumes block associated with the sidesets).
  Currently, the mesh is not fully connected, in that there is one set of
  nodes and elements for each face on the first sideset.  The adjacent
  face on the first sideset has a different set of nodes and elements
  (some of the nodes are coincident with the ones on adjacent faces).

  The newly created blocks have an `Attribute` field defined on the elements, which
  specifies the faces on the left and right sidesets that each shell
  element is contained in.  The name of the field is

  ```
  ci7Hl3_nonconformal_interface_gid_field
  ```

  There are 4 values per element, and they are:

  `sideset 1 element GID, sideset 1 face ordinal, sideset 2 element GID,
  sideset 2 face ordinal`.

  These values can be used to match with the Exodus sideset information
  in the first output mesh,
  (which IOSS stores in the `element_side`  field on each `SideBlock`.
  Currently the values are stored as `double`s, because IOSS does not support
  integer-valued `Attribute` fields.  This may change when IOSS support
  improves.

  Note that all names associated with the nonconformal interface have the
  prefix `ci7Hl3`, which was generated randomly.  The intention is that
  this will prevent name collisions with other blocks and fields, and also
  provide an easy means of filtering blocks when they are loaded by the
  solver.
