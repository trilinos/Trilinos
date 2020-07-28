## Copyright(C) 1999-2020 National Technology & Engineering Solutions
## of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
## NTESS, the U.S. Government retains certain rights in this software.
##
## See packages/seacas/LICENSE for details
## General Properties

  Property | Value    | Description
 ----------|----------|------------
 LOGGING   | on/\[off] | enable/disable logging of field input/output
 LOWER\_CASE\_VARIABLE\_NAMES | \[on]/off | Convert all variable names on database to lowercase; replace ' ' with '_'
 VARIABLE\_NAME\_CASE | upper/lower | Convert all variable names on output database to upper or lower case
 USE\_GENERIC\_CANONICAL\_NAMES | on/\[off]  | use `block_{id}` as canonical name of an element block instead of the name (if any) stored on the database. The database name will be an alias.
 ENABLE\_FIELD\_RECOGNITION | \[on]/off | try to combine scalar fields with common basename and recognized suffix into vector, tensor, ...
 FIELD\_SUFFIX\_SEPARATOR  | character \['_'] | use this suffix as separtor between field basename and suffices when recognizing fields
 MINIMIZE\_OPEN\_FILES | on/\[off] | If on, then close file after each timestep and then reopen on next output
 CYCLE\_COUNT | {int}\[infinite] | See `OVERLAY_COUNT` description below
 OVERLAY\_COUNT | {int}\[0] | Output `OVERLAY_COUNT` steps to database on top of each other.  `DB_STEP = (((IOSS_STEP-1) / (OVERLAY_COUNT+1)) % CYCLE_COUNT) +1`
## Auto-Decomposition-Related Properties

 Property        | Value  | Description
-----------------|--------|-----------------------------------------------------------
MODEL\_DECOMPOSITION\_METHOD | {method} | Decompose a DB with type `MODEL` using `method`
RESTART\_DECOMPOSITION\_METHOD | {method} | Decompose a DB with type `RESTART_IN` using `method`
DECOMPOSITION\_METHOD | {method} | Decompose all input DB using `method`
PARALLEL\_CONSISTENCY | \[on]/off | On if the client will call Ioss functions consistently on all processors. If off, then the auto-decomp and auto-join cannot be used.
RETAIN\_FREE\_NODES | \[on]/off | In auto-decomp, will nodes not connected to any elements be retained.
RETAIN\_EMPTY\_BLOCKS | on/\[off] | Empty blocks will / won't be retained in model. If retained, will have topology type "unknown".
LOAD\_BALANCE\_THRESHOLD | {real} \[1.4] | CGNS-Structured only -- Load imbalance permitted Load on Proc / Avg Load
LINE\_DECOMPOSITION | string | a list of comma-separated BC names. Zone with this bc will not be decomposed perpendicular to this surface. If name is `__ordinal_{ijk}` then use {ijk} as ordinal not to decompose.

### Valid values for Decomposition Method

Method     | Description
------------|-------------------
rcb        | recursive coordinate bisection
rib        | recursive inertial bisection
hsfc       | hilbert space-filling curve
metis\_sfc  | metis space-filling-curve
kway       | metis kway graph-based
kway\_geom  | metis kway graph-based method with geometry speedup
linear     | elements in order first n/p to proc 0, next to proc 1.
cyclic     | elements handed out to id % proc\_count
random     | elements assigned randomly to processors in a way that preserves balance (do not use for a real run)
external   | Files are decomposed externally into a file-per-processor in a parallel run.

## Output File Composition -- Single File output from parallel run instead of file-per-processor

 Property        | Value  | Description
-----------------|--------|-----------------------------------------------------------
COMPOSE\_RESTART  | on/\[off] |
COMPOSE\_RESULTS  | on/\[off] |
PARALLEL\_IO\_MODE | netcdf4, hdf5, pnetcdf | mpiio and mpiposix are deprecated hdf5=netcdf4

## Properties Related to byte size of reals and integers

 Property              | Value  | Description
-----------------------|--------|-----------------------------------------------------------
 INTEGER\_SIZE\_DB       | \[4] / 8 | byte size of integers stored on the database.
 INTEGER\_SIZE\_API      | \[4] / 8 | byte size of integers used in api functions.
 REAL\_SIZE\_DB          | 4 / \[8] | byte size of floating point stored on the database.
 REAL\_SIZE\_API         | 4 / \[8] | byte size of floating point used in api functions.

## Properties related to underlying file type (exodus only)

 Property              | Value  | Description
-----------------------|--------|-----------------------------------------------------------
  FILE\_TYPE            | \[netcdf], netcdf4, netcdf-4, hdf5 |
 COMPRESSION\_LEVEL     | \[0]-9    | In the range \[0..9]. A value of 0 indicates no compression, will automatically set `file_type=netcdf4`, recommend <=4
 COMPRESSION\_SHUFFLE   | on/\[off] |to enable/disable hdf5's shuffle compression algorithm.
 MAXIMUM\_NAME\_LENGTH   | \[32]     | Maximum length of names that will be returned/passed via api call.
 APPEND\_OUTPUT         | on/\[off] | Append output to end of existing output database
 APPEND\_OUTPUT\_AFTER\_STEP | {step}| Max step to read from an input db or a db being appended to (typically used with APPEND\_OUTPUT)
 APPEND\_OUTPUT\_AFTER\_TIME | {time}| Max time to read from an input db or a db being appended to (typically used with APPEND\_OUTPUT)

## Properties for the heartbeat output
 Property              | Value  | Description
-----------------------|--------|-----------------------------------------------------------
  FLUSH\_INTERVAL       | int   | For heartbeat, the minimum time interval in seconds between flushing heartbeat data to disk.  Default is 10 seconds
  FLUSH\_INTERVAL       | int   | For non-heartbeat, the number of output steps between flushing data to disk; if 0, then no flush
  TIME\_STAMP\_FORMAT    | \[%H:%M:%S] | Format used to format time stamp.  See strftime man page
  SHOW\_TIME\_STAMP      | on/off | Should the output lines be preceded by the timestamp
  PRECISION            | 0..16 \[5] | Precision used for floating point output.
  FIELD\_WIDTH          | 0.. |  Width of an output field. If 0, then use natural width.
  SHOW\_LABELS          | on/\[off]  | Should each field be preceded by its name (ke=1.3e9, ie=2.0e9)
  SHOW\_LEGEND          | \[on]/off  | Should a legend be printed at the beginning of the output showing the field names for each column of data.
  SHOW\_TIME\_FIELD      | on/\[off]  | Should the current analysis time be output as the first field.

## Experimental

 Property              | Value  | Description
-----------------------|--------|-----------------------------------------------------------
MEMORY\_READ        | on/\[off]   | experimental
MEMORY\_WRITE       | on/\[off]   | experimental
ENABLE\_FILE\_GROUPS | on/\[off]   | experimental

## Debugging / Profiling

  Property | Value    | Description
 ----------|----------|------------
 SHOW_CONFIG | ignored | output the build configuration of IOSS and TPL. Show supported database types.
 LOGGING   | on/\[off] | enable/disable logging of field input/output
 DECOMP\_SHOW\_PROGRESS | on/\[off] | show memory and elapsed time during autodecomp.
 DECOMP\_SHOW\_HWM      | on/\[off] | show high-water memory during autodecomp
 IOSS\_TIME\_FILE\_OPEN\_CLOSE | on/\[off] | show elapsed time during parallel-io file open/close/create
 CHECK\_PARALLEL\_CONSISTENCY | ignored | check Ioss::GroupingEntity parallel consistency

## Setting properties via an environment variable

Although the properties are usually accessed internally in the
application calling the IOSS library, it is possible to set the
properties externally prior to running the application via the setting
of the environment variable `IOSS_PROPERTIES`.  The value of the
variable is one or more colon-separated property/property-value pairs.
For example, to set the `DECOMPOSITION_METHOD` and the `FILE_TYPE`
externally, the following would be used:
```
    export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rib:FILE_TYPE=netcdf4"
```
If the environment variable is set correctly, there should be an
informational message output during running of the application similar
to:
```
	IOSS: Adding property 'DECOMPOSITION_METHOD' with value 'rib'
	IOSS: Adding property 'FILE_TYPE' with value 'netcdf4'
```
