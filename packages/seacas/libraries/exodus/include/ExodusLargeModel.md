\page large_model Large model modifications to the Exodus library
\section large_model_ Large model modifications to the Exodus library

\note The changes documented below are now the default method for
storing an Exodus database. There is no reason to write a "non-large"
model database.

The changes are to support the storage of larger models. There are two
pieces of this. The first is the setting of the type of netCDF file
that will be created; either one with 32-bit offsets or one with
64-bit offsets. (Note that this is *not* related to storing 64-bit
integers). This can be specified in a couple ways:

1. Pass the `EX_LARGE_MODEL` flag in the mode argument to ex_create.
2. Set the environment variable `EXODUS_LARGE_MODEL`.

If either of these are set, then the library will pass the
`NC_64BIT_OFFSET` flag to the netCDF library. See the commit log for
netCDF library for more details.

The other change is to reduce the size of some of the datasets in an
Exodus library. Even with the netCDF changes, the maximum dataset
size is still 2 GiB. To reduce the size of the datasets, the nodal
coordinates and the nodal variables have been split to store by
component. The old behavior stored all x, y, and z coordinates in a
single dataset; in the new behavior, each component is stored
separately -- there is a coordx, coordy, and coordz dataset.

The nodal variables used to be stored in a single blob of dimension
`(#times,#nodes,#variables)`. This has now been split into `#variable`
datasets of size `(#times,#nodes)`.

These two changes should increase the maximum model sizes
significantly. Prior to the change, the maximum number of nodes that
could be stored in the coordinate dataset was about 90 Million nodes;
the new storage permits 270 Million nodes in double precision. The old
model was more restrictive if there were multiple nodal variables, but
the new storage should not depend on the number of nodal variables.

These changes were made such that the new library would create
old-style files by default and would read either old or new style
files. The version has been changed to 3.01 for the file version and
4.01 for the API version.

An additional attribute is now written to the file. It is called
"file_size" or `ATT_FILESIZE`. If it is 0 or not present, then the old
format is assumed; if it is 1, then the new format is assumed.

There is also a new internal function called `ex_large_model(int exoid)`
which will return 1 if new version; 0 if old version.

If the function is passed a negative exoid, then it will check the
environment variable `EXODUS_LARGE_MODEL` and return 1 if it is
defined. It also currently prints a warning message saying that the
large model size was selected via the environment variable.

If you are using the Exodus api, then the only change to the client
application is the passing of the `EX_LARGE_MODEL` flag to `ex_create` or
the setting of the `EXODUS_LARGE_MODEL` environment variable. If your
client application is reading the database, no changes are needed.

========================================================================
If your client application bypasses some or all of the Exodus API
and makes direct netCDF calls, you will need to modify the calls.  The
changes that were made are shown below along with the name of the
Exodus API function in which the changes were made.

Alternatively, you can look at the changes that were made to the API at
http://github.com/sandialabs/seacas/packages/seacas/libraries/exodus/src.
The files that were changed are:

\note The filenames have been changed since this was written.
*   exgnvt.c
*   expcor.c
*   expini.c
*   expnv.c
*   expvp.c
*   expvpc.c
*   ex_utils.c
*   excopy.c
*   excre.c
*   exgcor.c
*   exgnv.c

\section ex_create ex_create():
-- Check whether the `EX_LARGE_MODEL` mode was set.  If so, then the
mode passed to nccreate must have the `NC_64BIT_OFFSET` bit set.  For
example, `mode |= NC_64BIT_OFFSET;`

  NOTE: `NC_64BIT_OFFSET` is defined in the Sandia's netCDF version
        "3.4-snl10X".  It should also be in netCDF-3.6.0 once it is released.

-- Write the exodus file size `ATT_FILESIZE` attribute (1=large, 0=normal):

```
filesiz = (nclong)(((cmode & EX_LARGE_MODEL) != 0) || (ex_large_model(-1) == 1));
   if (ncattput (exoid, NC_GLOBAL, ATT_FILESIZE, NC_LONG, 1, &filesiz) == -1)
    ... handle errors...
```

\section ex_put_init ex_put_init():
-- If writing a "large model" capable database, then the coordinates
are defined as components instead of an array.  The variables are
`VAR_COORD_X`, `VAR_COORD_Y` (if 2D or 3D), `VAR_COORD_Z` (if 3D). If not,
define the `VAR_COORD` variable as is currently done.

```
     if (ex_large_model(exoid) == 1) {
       /* node coordinate arrays -- separate storage... */
       dim[0] = numnoddim;
       if (ncvardef (exoid, VAR_COORD_X, nc_flt_code(exoid), 1, dim) == -1)
         { ... handle error }

       if (num_dim > 1) {
         if (ncvardef (exoid, VAR_COORD_Y, nc_flt_code(exoid), 1, dim) == -1)
           { ... handle error }
       }

       if (num_dim > 2) {
         if (ncvardef (exoid, VAR_COORD_Z, nc_flt_code(exoid), 1, dim) == -1)
           { ... handle error }
       }
     } else {
       /* node coordinate arrays: -- all stored together (old method) */
       .... define the old way...
     }
```

\section ex_put_coord ex_put_coord():
-- If writing a "large model" capable database, then the coordinates
are written a component at a time, otherwise write the old way as a single blob.

```
  if (ex_large_model(exoid) == 0) {
    ... write coordinates old way...
  } else {
    if ((coordidx = ncvarid (exoid, VAR_COORD_X)) == -1)
      { ... handle error }

    if (num_dim > 1) {
      if ((coordidy = ncvarid (exoid, VAR_COORD_Y)) == -1)
        { ... handle error }
    } else {
      coordidy = 0;
    }
    if (num_dim > 2) {
      if ((coordidz = ncvarid (exoid, VAR_COORD_Z)) == -1)
        { ... handle error }
    } else {
      coordidz = 0;
    }
    /* write out the coordinates  */
    for (i=0; i<= num_vars; i++) {
      dims[0] = time_dim;
      dims[1] = num_nod_dim;
      if ((ncvardef (exoid, VAR_NOD_VAR_NEW(i),
                     nc_flt_code(exoid), 2, dims)) == -1)
        {  ... handle error ... }
    }
  }
```

\section ex_put_nodal_var ex_put_nodal_var():
 -- If the large model method, write the nodal variable data to the correct variable;
    if the old method, determine the location within the blob

```
   if (ex_large_model(exoid) == 0) {
     /* write values of the nodal variable */
     if ((varid = ncvarid (exoid, VAR_NOD_VAR)) == -1) {
       ... handle error...
     }
     start[0] = --time_step;
     start[1] = --nodal_var_index;
     start[2] = 0;

     count[0] = 1;
     count[1] = 1;
     count[2] = num_nodes;
   } else {
     /* nodal variables stored separately, find variable for this variable
        index */
     if ((varid = ncvarid (exoid, VAR_NOD_VAR_NEW(nodal_var_index))) == -1) {
       ... handle error ...
     }

     start[0] = --time_step;
     start[1] = 0;

     count[0] = 1;
     count[1] = num_nodes;
   }

   if (ncvarput (exoid, varid, start, count,
                 ex_conv_array(exoid,WRITE_CONVERT,nodal_var_vals,num_nodes)) == -1) {
     ...handle error ...
   }
```

There are similar modifications to the reading of the nodal coordinates
and the reading of nodal variables. If interested in those functions, see the viewcvs URL
listed above.

If there are any questions, add an issue at
  https://github.com/sandialabs/seacas/issues
