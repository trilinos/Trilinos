# Zellij

Zellij is a "mesh concatenation" application for generating a mesh
consisting of a "lattice" containing one or more "unit cell" template
meshes.  The lattice is a two-dimensional arrangement of the unit cell
template meshes.

The unit cell template meshes are placed by zellij into the specified
locations in the lattice and the nodes on the boundaries of the unit
cell meshes are united or coincident.  Each unit cell mesh must have
the same exterior boundary meshes and coordinate extents on the X and
Y coordinate faces, but the Z faces are only required to have the same
coordinate extent; the Z face meshes are not required to be the same
among the different unit cells.

The lattice can be represented as a IxJ regular grid with each "cell"
in the grid or lattice containing one of the unit cell template
meshes.

[TOC]

## Execution

Executing zellij with the `-help` option will result in output similar to the following:

```
Zellij
        (A code for tiling 1 or more template databases into a single output database.)
        (Version: 1.4.1) Modified: 2021/03/16
        Parallel Capability Not Enabled.

usage: zellij [options] -lattice <lattice_definition_file>
        -lattice <$val> (Name of file to read lattice definition from. [required])
        -output <$val> (Name of output file to create. Default is `zellij-out.e`)

        -rcb (Use recursive coordinate bisection method to decompose the input mesh in a parallel run.)
        -rib (Use recursive inertial bisection method to decompose the input mesh in a parallel run.)
        -hsfc (Use hilbert space-filling curve method to decompose the input mesh in a parallel run. [default])
        -linear (Use the linear method to decompose the input mesh in a parallel run.
                Elements in order first n/p to proc 0, next to proc 1.)
        -cyclic (Use the cyclic method to decompose the input mesh in a parallel run.
                Elements handed out to id % proc_count)
        -random (Use the random method to decompose the input mesh in a parallel run.
                Elements assigned randomly to processors in a way that preserves balance
                (do _not_ use for a real run))

        -ranks <$val> (Number of ranks to decompose mesh across)
        -start_rank <$val> (In partial output mode, start outputting decomposed files at this rank)
        -rank_count <$val> (In partial output or subcycle modes, output this number of ranks)
        -subcycle (Process cells in groups of '-rank_count'.  Helps minimize open files,
                but is faster than only having a single file open.)
        -scale <$val> (Scale the output mesh coordinates by the specified value)
        -minimize_open_files [$val] (Close files after accessing them to avoid issues with too many open files.
                If argument is 'output' then close output, if 'unit' then close unit cells;
                if 'all' or no argument close all.
                Should not need to use this option unless you get an error message indicating this issue.)

        -ignore_sidesets (Do not copy any sidesets in the unit cells to the output file.)
        -generate_sidesets <$val> (Which surfaces on the output mesh should have sidesets generated,
                 Valid options are:
                 'x' or 'i' for surface on minimum X coordinate, default name = `min_i`
                 'y' or 'j' for surface on minimum Y coordinate, default name = `min_j`
                 'z' or 'k' for surface on minimum Z coordinate, default name = `min_k`
                 'X' or 'I' for surface on maximum X coordinate, default name = `max_i`
                 'Y' or 'J' for surface on maximum Y coordinate, default name = `max_j`
                 'Z' or 'K' for surface on maximum Z coordinate, default name = `max_k`
                 For example `xyXY` would generate sidesets on min/max X and Y surfaces.)
        -sideset_names <$val> (Specify names for one or more of the generated sidesets.
                 Form is `axis:name,axis:name,...`
                 where 'axis' is one of 'ijkIJKxyzXYZ', and 'name' is the name of the sideset.
                 The default names are 'min_i', 'max_i', 'min_j', 'max_j', 'min_k', 'max_k'.
                 For example `x:left,X:right` would name the sideset on the min x face 'left' and the max X face 'right'.)

        -netcdf3 (Output database will be a netcdf3 native classical netcdf file format (32-bit only))
        -netcdf4 (Output database will be a netcdf4 hdf5-based file instead of the classical netcdf file format (default))
        -netcdf5 (Output database will be a netcdf5 (CDF5) file instead of the classical netcdf file format)

        -32-bit (True if forcing the use of 32-bit integers for the output file)
        -64-bit (True if forcing the use of 64-bit integers for the output file (default))

        -zlib (Use the Zlib / libz compression method if compression is enabled (default) [exodus only].)
        -szip (Use SZip compression. [exodus only, enables netcdf-4])
        -compress <$val> (Specify the hdf5 zlib compression level [0..9] or szip [even, 4..32] to be used on the output file.)

        -separate_cells (Do not equivalence the nodes between adjacent unit cells.)
        -repeat <$val> (Each lattice entry will be used the specified number of times as will
                each row in the lattice (for debugging). `-repeat 2` would double the lattice.)
        -skip <$val> (Skip the specified number of lattice entries and rows. For example, -skip 1
                would read every other entry on the row and every other row. (for debugging))
        -help (Print this summary and exit)
        -version (Print version and exit)
        -debug <$val> (debug level (values are or'd)
                   1 = Time stamp information.
                   2 = Memory information.
                   4 = Verbose Unit Cell information.
                   8 = Verbose output of Grid finalization calculations.
                  16 = Put exodus library into verbose mode.
                  32 = Verbose decomposition information.
                  64 = Verbose output database summary information.
                 128 = Verbose sideset generation information.)
        -copyright (Show copyright and license data.)

        Can also set options via ZELLIJ_OPTIONS environment variable.

        ->->-> Send email to sierra-help@sandia.gov for zellij support.<-<-<-
```
The only required option is `-lattice` followed by the name of the file containing the lattice description.  The other options are used to specify compression of the output file; the format of the output file; or to request additional debug output.

If the `-output <filename>` option is not specified, then the output mesh will be named `zellij-out.e`.

## Lattice Description File Format
The format of the lattice description file is fairly simple, but is also very rigid.  There are two sections of the file -- the _unit cell_ dictionary and the _lattice_ definition.

### Unit Cell Dictionary
The unit cell dictionary defines the unit cell template meshes that will be placed in the lattice.  The dictionary begins with a line containing `BEGIN_DICTIONARY` followed by one or more lines defining the unit cells and is then ended with a line containing `END_DICTIONARY`

The syntax of the lines defining the unit cells consists of two fields -- an arbitrary _key_ and the filename containing the Exodus file defining the mesh for this unit cell.  The only restriction on the _key_ is that it must be unique in the dictionary.  The filenames must specify the path (either absolute or relative to the current execution directory) to the Exodus file; it can optionally be delimited by double quotes (`"`).  The filenames do not need to be unique, but it is more efficient in both memory and time if each unit cell template mesh is unique.

As an example, here is a valid dictionary definition:

```
BEGIN_DICTIONARY
  0001 "../zellij-example/xatom-1b.e"
  0002 "../zellij-example/xatom-Y.e"
  0003 "../zellij-example/xatom-X.e"
  0004 "../zellij-example/xatom-2b.e"
END_DICTIONARY
```
The unit cell dictionary must appear before the lattice definition in the lattice description file.

If an error is detected during the parsing of the unit cell dictionary, the code will output an error message and terminate.  Errors can be incorrect syntax, missing unit cell template meshes, duplicate keys, or problems reading the mesh description from a unit cell template mesh.  The unit cell template mesh file is accessed and partially read at the time that zellij parses the corresponding unit cell dictionary line.

### Lattice Definition

The lattice definition specifies the size of the lattice and the distribution of the unit cell(s) within that lattice.  The lattice definition must follow the unit cell dictionary in the lattice description file.

The first line of the lattice definition begins with the line `BEGIN_LATTICE {i} {j} 1` where `{i}` and `{j}` specify the size of the `IxJ` arrangement of unit cells.  For example, the line `BEGIN_LATTICE 5 5 1` would define a lattice containing 25 unit cell instances arranged in a 5 by 5 regular grid.

The last line of the lattice definition is the line `END_LATTICE`.  When that line is encountered, zellij will begin outputting the mesh.

Between the `BEGIN_LATTICE` and `END_LATTICE` are `{j}` lines with `{i}` entries per line.  The entries are any of the _key_s that were specified in the unit cell dictionary.

As an example, here is a valid lattice definition using the keys of the example dictionary from the previous section:

```
BEGIN_LATTICE 5  5  1
  0001 0002 0003 0002 0001
  0002 0003 0003 0003 0002
  0003 0003 0004 0003 0003
  0002 0003 0003 0003 0002
  0001 0002 0003 0002 0001
END_LATTICE
```

Although the lattice is typically symmetric and square, this is not a requirement and is not checked.

If an error is detected during the parsing of the _lattice_, the code will output an error message and terminate.  Errors can include invalid keys, incorrect number of lattice definition lines, or incorrect number of keys on a definition line.

Note that zellij does not require that the unit cell keys be numeric; the following example shows a different method for specifying the same lattice definition file as the previous example:

```
BEGIN_DICTIONARY
  - "../zellij-example/xatom-1b.e"
  | "../zellij-example/xatom-Y.e"
  + "../zellij-example/xatom-X.e"
  * "../zellij-example/xatom-2b.e"
END_DICTIONARY

BEGIN_LATTICE 5  5  1
  - | + | -
  | + + + |
  + + * + +
  | + + + |
  - | + | -
END_LATTICE
```

## Unit Cell Template Mesh Requirements
Zellij requires that the boundary mesh (`X` and `Y` faces) of each of the unit cell templates be a _regular_ "structured" mesh.  Basically this means that the faces of the mesh elements on the boundary are in a regular rectangular grid such that each mesh face is rectangular (90 degree corners) and that the boundary mesh on the minimum `X` face is the same as that on the maximum `X` face and similarly for the minimum `Y` face and the maximum `Y` face.

Additionally, the X faces on _all_ unit cells must match and the Y faces on _all_
unit cells must match both in structure and in coordinate extent. This requirement is verified during execution. The `Z` faces are less constrained with the only requirement being that the coordinate extents of all `Z` faces must be the same (which follows from the `X` and `Y` face requirement); the structure of the mesh on the `Z` faces is arbitrary.

The unit cell meshes can contain any number of element blocks; however, each element block _must_ contain hexahedral elements with 8-nodes per element.  The element blocks do not need to be the same in each unit cell mesh, but if they do share the same element block `id`, then those elements will be combined into the same element block in the output mesh with the same `id`.

The output mesh will contain the union of all element blocks existing on the input mesh unit cells.  For example, if:

*  unit cell `0001` has element blocks `1 10 100`
*  unit cell `0002` has element blocks `2 20 200`
*  unit cell `0003` has element blocks `1 2 10 20`
*  unit cell `0004` has element blocks `10 20 100 200`

The output mesh will have element blocks `1 2 10 20 100 200`

## Sideset Handling
By default, zellij will replicate any sidesets that are defined on the
input unit cell meshes to the output mesh file.  The sidesets will
have the same names and ids as the sidesets on the input unit cell
meshes.  If you do not want the sidesets replicated, you can add the
command line option `-ignore_sidesets` and any sidesets on the input
unit cell meshes will be ignored.

Zellij can also generate new sidesets on the boundaries of the output
mesh via the command line option `-generate_sidesets <axes>` where
`axes` is one or more letters specifying the face of the output mesh
on which to generate a sideset.  Valid letters are `xyzXYZ` or
`ijkIJK` which correspond to:

*  `x` or `i` for surface on minimum X coordinate (default name = `min_i`)
*  `y` or `j` for surface on minimum Y coordinate (default name = `min_j`)
*  `z` or `k` for surface on minimum Z coordinate (default name = `min_k`)
*  `X` or `I` for surface on maximum X coordinate (default name = `max_i`)
*  `Y` or `J` for surface on maximum Y coordinate (default name = `max_j`)
*  `Z` or `K` for surface on maximum Z coordinate (default name = `max_k`)

For example `-generate_sidesets xyXY` would generate sideset on the
surfaces corresponding to the minimum and maximum X and Y coordinates
on the output mesh.

By default, the generated sidesets will be named as shown in the list
above.  The names can be changed with the `-sideset_names <arg>`
command line option.  The syntax of `<arg>` is
`axis:name,axis:name,...` where `axis` is one of `ijkIJK` or `xyzXYZ`
and `name` is the name of the specified sideset.  For example,
`-sideset_names x:left,X:right`  would name the sidesets on the
minimum x and maximum X faces `left` and `right` respectively.  There
will be an error if two or more sidesets have the same name.

## Parallel Execution
Zellij can produce a mesh decomposed into a file-per-rank for use in a
parallel analysis application.  Note that Zellij itself is run
serially.  The main option that tells Zellij to produce the decomposed
files is `-ranks <number_of_ranks>`.  If this is specified, then
Zellij will create `number_of_ranks` individual files each containing
a portion of the complete model.  The files will have the information
needed by a parallel application to read the data and set up the
correct communication paths and identify the nodes that are shared
across processor boundaries.

The decomposition method can also be specified.  This determines the
algorithm that is used to break the lattice into `number_of_ranks`
pieces each with approximately the same computational complexity.  The
decomposition methods are:

*   `-rcb` Use recursive coordinate bisection method to decompose the input mesh in a parallel run.
*   `-rib` Use recursive inertial bisection method to decompose the input mesh in a parallel run.
*   `-hsfc` Use hilbert space-filling curve method to decompose the input mesh in a parallel run.
*   `-linear` Use the linear method to decompose the input mesh in a parallel run. Elements in order first `n/p` to proc 0, next to proc 1.
*   `-cyclic` Use the cyclic method to decompose the input mesh in a parallel run.      Elements handed out to `id % proc_count`.
*   `-random` Use the random method to decompose the input mesh in a parallel run.      Elements are assigned randomly to processors in a way that preserves balance (do _not_ use for a real run))

The `-hsfc` method is the default if no other decomposition method is
specified. Note that the decomposition occurs at the _grid_ level so
the elements of each grid cell will not be split across multiple ranks.  The grid
cells are weighted by the number of elements in the cell which should
produce a balanced decomposition if there are unit cells of varying
element counts.

The `-linear`, `-cyclic`, and `-random` methods are typically used for
debugging and testing Zellij and should not be used in a production
run, especially the `-random` method.

### Partial Parallel Output Mode
There is a _partial parallel output_ mode in which you can tell Zellij to only output a portion of the parallel decomposed files.  This is selected with the `-start_rank <rank>` and `-rank_count <count>` options.  In this case, Zellij will only output the ranks from `rank` up to `rank+count-1`.  For example, if you run `zellij -ranks 10 -start_rank 5 -rank_count 3`, then zellij would output files for ranks 5, 6, and 7.  This is somewhat inefficient since zellij will do many of the calculations for all ranks and only output the specified ranks; however, it does allow you to run multiple copies of zellij simultaneously.  For example, you could run:

```
 zellij -ranks 16 --start_rank  0 --rank_count 4
 zellij -ranks 16 --start_rank  4 --rank_count 4
 zellij -ranks 16 --start_rank  8 --rank_count 4
 zellij -ranks 16 --start_rank 12 --rank_count 4
```

simultaneously and all 16 files should be output faster than running a single execution that wrote all of the files.

### Parallel Capable Parallel Execution

If Zellij is compiled with parallel capability enabled (This is shown at the beginning of the `-help` output or the version
information output when zellij begins executing as `Parallel Capability Enabled`), then you can run Zellij in parallel using the
normal `mpiexec -np <#> zellij <normal zellij options>` command.  In this case, there will be `#` copies of zellij running
simultaneously and each copy will divide up the output files and work among each process/copy.

For example, if you run `mpiexec -np 8 zellij -ranks 1024 -latice lattice.txt`, then there will be 8 copies of zellij running
and each will output `1024/8 = 128` output files.

### Maximum Open File Complications

Most compute systems have a limit on the number of files that a program can have open simultaneously.  For many systems, this
limit is 1024.  The files that zellij deals with are (1) the unit cell meshes and (2) the per-rank output files, and (3) the
standard input, output, and error files. Because of this, it is somewhat easy for a zellij execution to exceed the open file
limit.  Zellij attempts to handle this automatically using logic similar to:

*  If the unit cell count exceeds the open file limit, then close each unit cell after each access before opening the next unit
cell mesh.

*  If the number of `-ranks` that zellij is creating exceeds the open file count, then determine how many output files can be
open at one time (max_open = open file limit - 3 - number of unit cells open simultaneously) and run zellij in a `subcycle` mode
where it is only writing to `max_open` files at one time.

*  If the `max_open` calculated in the above bullet is too small, then set the mode to only open a single unit cell mesh at a
time and redo the calculation.

*  If all else fails, run with only a single unit cell file open and only a single output mesh rank file open.

If the above logic fails and Zellij is unable to run without exceeding the open file count, you can specify the behavior
manually using a combination of the `-minimize_open_files=<UNIT|OUTPUT|ALL>` option and the `-subcycle` and `-rank_count <#>`
options.

The options to `-minimize_open_files` are:

*  `UNIT` - only have a single unit cell mesh open at one time; close before accessing another unit cell mesh.
*  `OUTPUT` - only have a single output rank mesh file open at one time.
*  `ALL` - both of the above options.

The `-subcycle` and `-rank_count <#>` options cause zellij to output `#` output files at a time and then cycle to the next `#`
output files until all files have been output.  For example, `zellij -ranks 1024 -subcycle -rank_count 256` would do the
following:

*  First cycle would output ranks 0 to 255,
*  Second cycle would output ranks 256 to 511,
*  Third cycle would output ranks 512 to 767,
*  Fourth cycle would output ranks 768 to 1023.

In this mode, there will the `#` output files open simultaneously (unless
`-minimize_open_files=OUTPUT|ALL` was specified also).  So the total number of open files will be `unit cell count + 3 + #` or
`1 + 3 + #` if `-minimize_open_files=UNIT` was specified.

## Execution Complexity

Zellij is intended to produce extremely large meshes and is therefore very concerned with both memory efficiency and execution
time efficiency.

### Memory Complexity

Zellij stores the following data:

*  For each unit cell template mesh:
  *  metadata
  *  64-bit Ids of nodes on each min_I, max_I, min_J, max_J face
*  For each entry in the lattice definition:
  *  metadata (approximately 1KiByte)
  *  temporarily it will hold 64-bit Ids of nodes on the max_I and max_J faces. This will be deleted once the upper `I` and upper `J` "neighbor" entry has been processed (see below)
*  For the lattice:
  *  vector containing the lattice definition.

The main memory use once the output file is being processed is the temporary storage containing the nodes on the `max_I` and `max_J` faces.  The lattice is processed cell by cell.  For an `II by JJ` sized grid, the cells are processed in the order `(1,1), (2,1), ... , (II, 1), (1,2), (2,2), ..., (II, JJ)`.  The temporary storage on the `max_I` face is only needed until the next cell is processed.  That is, for cell `(i,j)`, its `max_I` nodes will be used during the processing of cell `(i+1, j)` and then deleted.

The temporary storage on the `max_J` face is retained for a longer time.  For cell `(i,j)`, the `max_J` storage is needed for cell `(i, j+1)` and then deleted.

For a grid of size `(II, JJ)`, there will at most be:

*  1 temporary vector of size `max_I` nodes
*  `II` temporary vectors of size `max_J` nodes.

If you have a lattice that is rectangular (`II != JJ`), then it is more efficient for memory usage to make the `I` direction the smallest value if possible.

In addition to the above memory usage, zellij must also transfer the mesh coordinate data and element block connectivity data for each lattice entry to the output file.  Zellij outputs the model using the following pseudo-code:

```
for each j : J
  for each i : I
     read cell(i,j)  x, y, and z local coordinates
     map coordinates to offset in output mesh
     eliminate nodes that join to an already output neighbor cell
     write cell(i,j) x, y, and z global coordinates

for each j : J
   for each i : I
      for each element block in cell(i,j) mesh
         read block connectivity
         map local node ids to global node ids
         write block connectivity
```

The maximum memory use will be the size of storage needed for the `x` `y` and `z` coordinates of a unit cell mesh or the storage needed to hold the connectivity for a single unit cell element block.

Note that the memory requirements are proportional to the size of an individual unit cell mesh and not a function of the size of the output mesh.  It is possible to create meshes which are much larger than the amount of memory present on the compute system running zellij.

The memory being used by zellij during execution will be output if the `--debug 2` argument is specified at execution time.

### Execution Time Complexity

For a large model, the majority of the execution time is related to:

*  Read/process/write element block connectivity
*  Read/process/write nodal coordinates
*  Categorize boundary nodes on each unit cell mesh

### Efficiency at the NetCDF level
The Exodus format which is used for the unit cell template meshes and the output mesh uses the NetCDF library for on-disk storage.  There are several variants of the NetCDF on-disk storage including the format: `netcdf3`, `netcdf4`, and `netcdf5` and the integer size (32-bit integers or 64-bit integers).  Although these details are usually transparent to the user, they can affect the execution time especially when very large meshes are being processed.

#### Format
The `netcdf3` format is the original native NetCDF format.  At the time the library was being developed, the `byte endianness` of data stored on disk was not standard among the computes in use at that time and the NetCDF developers had to pick an `endianness` for the data.  They picked the XDR standard which stood for _eXternal Data Representation_ which was used for communicating between different computer systems.  Regretfully, the representation used by XDR turned out to be opposite of the representation used by (almost?) all systems in use today, so each read and write of data in the `netcdf3` format results in a translation of the endianness.  This translation is very fast, but is overhead that would not be needed if the on-disk format was the opposite representation.  This representation is also used by the `netcdf5` format.

However, the NetCDF `netcdf4` format is based on using the HDF5 library to manage the underlying data format on disk and it can read and write data using the native endianness of the system on which the data is being read and written and therefore does not incur the cost of transforming the data's endianness.

#### Integer Size
By default, most current mesh generators will output a mesh using 32-bit integer data.  This is sufficient to represent a mesh with up to approximately 2.1 billion nodes and elements.

If the input mesh and the output mesh have the same integer size, then there is no data conversion needed.  The data will be read as `N`-bit integers, processed as `N`-bit integers, and written as `N`-bit integers.  However, if the input mesh is `N`-bit integers and the output mesh is `M`-bit integers, then the NetCDF library will convert all integer data (element block connectivity typically) from `N` bits to `M` bits which for large meshes can incur an execution time overhead.

#### Compression
The NetCDF library supports compression of the output file.  Typically, the `zlib` compression algorithm is used, but recently NetCDF begain supporting the `szip` compression and a few more algorithms are starting to be supported.

The benefit of the compression is that it can result in much smaller output (and input) mesh files; the disadvantage is that the default `zlib` compression algorithm is not very fast and can increase the execution time of zellij.  The `szip` compression algorithm is faster with typically (but not always) slightly less compression, but it still will incur an overhead in execution time.

#### Recommendations
For minimal overhead, it is recommended that:

*  Use the `netcdf4` format for all input and output meshes
*  Use the same integer size for all input and output meshes
  *  The integer size of the output mesh can be specified using the `-32` or `-64` options.
  *  The `-64` option is the default.

It is most efficient if the format and integer size of the input mesh matches the output mesh.  The format of the input meshes can be converted using the `io_shell` application with the `-netcdf4` and `-64` or `-32` options. The format and integer size of a mesh can be queried using the `exo_format` application.

For illustration, here is the execution time for several runs with different format and integer size.  In all cases, the input and output mesh sizes are the same:

| input   | output  | integer input | integer output | execution time |
|:-------:|:-------:|:-------------:|:--------------:|:--------------:|
| netcdf3 | netcdf3 |       32      |       32       |  7.0           |
| netcdf3 | netcdf4 |       32      |       32       |  2.6           |
| netcdf3 | netcdf4 |       32      |       64       |  3.8           |
| netcdf4 | netcdf3 |       32      |       32       |  6.5           |
| netcdf4 | netcdf3 |       64      |       32       |  7.4           |
| netcdf4 | netcdf5 |       64      |       64       |  9.4           |
| netcdf4 | netcdf4 |       32      |       32       |  2.4           |
| netcdf4 | netcdf4 |       32      |       64       |  3.6           |
| netcdf4 | netcdf4 |       64      |       32       |  3.2           |
| netcdf4 | netcdf4 |       64      |       64       |  3.3           |

The fastest option is both input and output using 32-bit integers and the `netcdf4` format.  Almost as fast is the case where the input format is `netcdf3` and the output `netcdf4`.  The `64-bit` integer options with both input and output using `netcdf4` are slightly slower, but this is probably due to the doubling of the size of the integer data being read and written.

The output mesh in this case consisted of 37.3 million elements and 38.5 million nodes in a grid of 46 x 46 unit cells.  There were 56 unit cell template meshes.
