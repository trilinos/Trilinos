/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*****************************************************************************
 *
 * exopen - ex_open
 *
 * entry conditions -
 *   input parameters:
 *       char*   path                    exodus filename path
 *       int     mode                    access mode w/r
 *
 * exit conditions -
 *       int     exoid                   exodus file id
 *       int*    comp_ws                 computer word size
 *       int*    io_ws                   storage word size
 *       float*  version                 EXODUS interface version number
 *
 * revision history -
 *
 *
 *****************************************************************************/
/* Determine whether compiling against a parallel netcdf... */
#include "exodusII.h"
#if defined(PARALLEL_AWARE_EXODUS)

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, etc
#include <mpi.h>          // for MPI_Comm, MPI_Info, etc
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <string.h>
/*!

\note The ex_open_par_int() is an internal function called by
ex_open_par(). The user should call ex_open_par() and not ex_open_par_int().

The function ex_open_par() opens an existing exodus file and returns
an ID that can subsequently be used to refer to the file, the word
size of the floating point values stored in the file, and the version
of the exodus database (returned as a ``float'', regardless of the
compute or I/O word size). Multiple files may be ``open'' simultaneously.

\return In case of an error, ex_open_par() returns a negative
number. Possible causes of errors include:
  -  The specified file does not exist.
  -  The mode specified is something other than the predefined constant
\fparam{EX_READ} or \fparam{EX_WRITE}.
  -  Database version is earlier than 2.0.

\param path The file name of the exodus file. This can be given as either an
            absolute path name (from the root of the file system) or a relative
            path name (from the current directory).

\param mode Access mode. Use one of the following predefined constants:
        -  \fparam{EX_READ} To open the file just for reading.
        -  \fparam{EX_WRITE} To open the file for writing and reading.

\param[in,out] comp_ws The word size in bytes (0, 4 or 8) of the floating point
variables
               used in the application program. If 0 (zero) is passed, the
default
               size of floating point values for the machine will be used and
               returned in this variable. WARNING: all exodus functions
requiring
               reals must be passed reals declared with this passed in or
returned
               compute word size (4 or 8).


\param[in,out] io_ws The word size in bytes (0, 4 or 8) of the floating
                    point data as they are stored in the exodus file. If the
word
                    size does not match the word size of data stored in the
file,
                    a fatal error is returned. If this argument is 0, the word
size
                    of the floating point data already stored in the file is
returned.

\param[out] version  Returned exodus database version number.

The following opens an exodus file named \file{test.exo} for parallel
read only, using default settings for compute and I/O word sizes:

~~~{.c}
int CPU_word_size,IO_word_size, exoid;
float version;

CPU_word_size = sizeof(float);   \co{float or double}
IO_word_size = 0;                \co{use what is stored in file}

\comment{open exodus files}
exoid = ex_open_par ("test.exo",     \co{filename path}
                     EX_READ,        \co{access mode = READ}
                     &CPU_word_size, \co{CPU word size}
                     &IO_word_size,  \co{IO word size}
                     &version,       \co{ExodusII library version
                     MPI_COMM_WORLD,
                     MPI_INFO_NULL);}
~~~
 */

static int warning_output = 0;

struct ncvar /* variable */
{
  char    name[MAX_VAR_NAME_LENGTH];
  nc_type type;
  int     ndims;
  int     dims[NC_MAX_VAR_DIMS];
  int     natts;
};

/* NOTE: Do *not* call `ex_open_par_int()` directly.  The public API
 *       function name is `ex_open_par()` which is a wrapper that
 *       calls `ex_open_par_int` with an additional argument to make
 *       sure library and include file are consistent
 */
int ex_open_par_int(const char *path, int mode, int *comp_ws, int *io_ws, float *version,
                    MPI_Comm comm, MPI_Info info, int run_version)
{
  int     exoid;
  int     status, stat_att, stat_dim;
  nc_type att_type = NC_NAT;
  size_t  att_len  = 0;
  int     nc_mode  = 0;
  int     old_fill;
  int     file_wordsize;
  int     dim_str_name;
  int     int64_status = 0;
  int     is_hdf5      = 0;
  int     is_pnetcdf   = 0;
  int     in_redef     = 0;

  char errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();

  /* set error handling mode to no messages, non-fatal errors */
  ex_opts(exoptval); /* call required to set ncopts first time through */

  if (run_version != EX_API_VERS_NODOT && warning_output == 0) {
    int run_version_major = run_version / 100;
    int run_version_minor = run_version % 100;
    int lib_version_major = EX_API_VERS_NODOT / 100;
    int lib_version_minor = EX_API_VERS_NODOT % 100;
    fprintf(stderr,
            "EXODUS: Warning: This code was compiled with exodus "
            "version %d.%02d,\n          but was linked with exodus "
            "library version %d.%02d\n          This is probably an "
            "error in the build process of this code.\n",
            run_version_major, run_version_minor, lib_version_major, lib_version_minor);
    warning_output = 1;
  }

  if ((mode & EX_READ) && (mode & EX_WRITE)) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Cannot specify both EX_READ and EX_WRITE");
    ex_err(__func__, errmsg, EX_BADFILEMODE);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (mode & EX_WRITE) {
    nc_mode = (NC_WRITE | NC_MPIIO);
  }
  else {
    nc_mode = (NC_NOWRITE | NC_SHARE | NC_MPIIO);
  }
  if ((status = nc_open_par(path, nc_mode, comm, info, &exoid)) != NC_NOERR) {
    /* It is possible that the user is trying to open a netcdf4
       file, but the netcdf4 capabilities aren't available in the
       netcdf linked to this library. Note that we can't just use a
       compile-time define since we could be using a shareable
       netcdf library, so the netcdf4 capabilities aren't known
       until runtime...

       Later versions of netcdf-4.X have a function that can be
       queried to determine whether the library being used was
       compiled with --enable-netcdf4, but not everyone is using that
       version yet, so we may have to do some guessing...

       At this time, query the beginning of the file and see if it
       is an HDF-5 file and if it is assume that the open failure
       is due to the netcdf library not enabling netcdf4 features unless
       we have the define that shows it is enabled, then assume other error...
    */
    int type = 0;
    ex_check_file_type(path, &type);

    if (type == 5) {
#if NC_HAS_HDF5
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the netcdf-4 "
              "file:\n\t'%s'\n\t failed. The netcdf library supports "
              "netcdf-4 so there must be a filesystem or some other "
              "issue \n",
              path);
#else
      /* This is an hdf5 (netcdf4) file. If NC_HAS_HDF5 is not defined,
         then we either don't have hdf5 support in this netcdf version,
         OR this is an older netcdf version that doesn't provide that define.

         In either case, we don't have enough information, so we
         assume that the netcdf doesn't have netcdf4 capabilities
         enabled.  Tell the user...
      */
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the netcdf-4 "
              "file:\n\t'%s'\n\tEither the netcdf library does not "
              "support netcdf-4 or there is a filesystem or some "
              "other issue \n",
              path);
#endif
    }
    else if (type == 4) {
#if NC_HAS_CDF5
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the CDF5 "
              "file:\n\t'%s'\n\t failed. The netcdf library supports "
              "CDF5-type files so there must be a filesystem or some other "
              "issue \n",
              path);
#else
      /* This is an cdf5 (64BIT_DATA) file. If NC_64BIT_DATA is not defined,
         then we either don't have cdf5 support in this netcdf version,
         OR this is an older netcdf version that doesn't provide that define.

         In either case, we don't have enough information, so we
         assume that the netcdf doesn't have cdf5 capabilities
         enabled.  Tell the user...
      */
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the CDF5 "
              "file:\n\t'%s'\n\tEither the netcdf library does not "
              "support CDF5 or there is a filesystem or some "
              "other issue \n",
              path);

#endif
    }
    else if (type == 1 || type == 2) {
#if NC_HAS_PNETCDF
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the classic NetCDF "
              "file:\n\t'%s'\n\t failed. The netcdf library supports "
              "PNetCDF files as required for parallel readinf of this "
              "file type, so there must be a filesystem or some other "
              "issue \n",
              path);
#else
      /* This is an normal NetCDF format file, for parallel reading, the PNetCDF
         library is required but that is not compiled into this version.
      */
      fprintf(stderr,
              "EXODUS: ERROR: Attempting to open the NetCDF "
              "file:\n\t'%s'\n\tThe NetCDF library was not "
              "built with PNetCDF support as required for parallel access to this file.\n",
              path);

#endif
    }

    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to open %s of type %d for reading. Either the "
             "file does not exist, or there is a permission or file format "
             "issue.",
             path, type);
    ex_err(__func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* File opened correctly */
  int type = 0;
  ex_check_file_type(path, &type);
  if (type == 5) {
    is_hdf5 = 1;
  }
  else if (type == 1 || type == 2 || type == 4) {
    is_pnetcdf = 1;
  }

  if (mode & EX_WRITE) { /* Appending */
    /* turn off automatic filling of netCDF variables */
    if (is_pnetcdf) {
      if ((status = nc_redef(exoid)) != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to put file id %d into define mode", exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
      in_redef = 1;
    }

    if ((status = nc_set_fill(exoid, NC_NOFILL, &old_fill)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to set nofill mode in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    stat_att = nc_inq_att(exoid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, &att_type, &att_len);
    stat_dim = nc_inq_dimid(exoid, DIM_STR_NAME, &dim_str_name);
    if (stat_att != NC_NOERR || stat_dim != NC_NOERR) {
      if (!in_redef) {
        if ((status = nc_redef(exoid)) != NC_NOERR) {
          snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to put file id %d into define mode",
                   exoid);
          ex_err_fn(exoid, __func__, errmsg, status);
          EX_FUNC_LEAVE(EX_FATAL);
        }
        in_redef = 1;
      }
      if (stat_att != NC_NOERR) {
        int max_so_far = 32;
        nc_put_att_int(exoid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &max_so_far);
      }

      /* If the DIM_STR_NAME variable does not exist on the database, we need to
       * add it now. */
      if (stat_dim != NC_NOERR) {
        /* Not found; set to default value of 32+1. */
        int max_name = ex_default_max_name_length < 32 ? 32 : ex_default_max_name_length;
        nc_def_dim(exoid, DIM_STR_NAME, max_name + 1, &dim_str_name);
      }
    }

    if (in_redef) {
      if ((status = ex_leavedef(exoid, __func__)) != NC_NOERR) {
        EX_FUNC_LEAVE(EX_FATAL);
      }
      in_redef = 0;
    }

    /* If this is a parallel execution and we are appending, then we
     * need to set the parallel access method for all transient variables to NC_COLLECTIVE since
     * they will be being extended.
     */
    int ndims;    /* number of dimensions */
    int nvars;    /* number of variables */
    int ngatts;   /* number of global attributes */
    int recdimid; /* id of unlimited dimension */

    int varid;

    /* Determine number of variables on the database... */
    nc_inq(exoid, &ndims, &nvars, &ngatts, &recdimid);

    for (varid = 0; varid < nvars; varid++) {
      struct ncvar var;
      nc_inq_var(exoid, varid, var.name, &var.type, &var.ndims, var.dims, &var.natts);

      if ((strcmp(var.name, VAR_GLO_VAR) == 0) || (strncmp(var.name, "vals_elset_var", 14) == 0) ||
          (strncmp(var.name, "vals_sset_var", 13) == 0) ||
          (strncmp(var.name, "vals_fset_var", 13) == 0) ||
          (strncmp(var.name, "vals_eset_var", 13) == 0) ||
          (strncmp(var.name, "vals_nset_var", 13) == 0) ||
          (strncmp(var.name, "vals_nod_var", 12) == 0) ||
          (strncmp(var.name, "vals_edge_var", 13) == 0) ||
          (strncmp(var.name, "vals_face_var", 13) == 0) ||
          (strncmp(var.name, "vals_elem_var", 13) == 0) ||
          (strcmp(var.name, VAR_WHOLE_TIME) == 0)) {
        nc_var_par_access(exoid, varid, NC_COLLECTIVE);
      }
    }
  } /* End of (mode & EX_WRITE) */

  /* determine version of EXODUS file, and the word size of
   * floating point and integer values stored in the file
   */

  if ((status = nc_get_att_float(exoid, NC_GLOBAL, ATT_VERSION, version)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get database version for file id: %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* check ExodusII file version - old version 1.x files are not supported */
  if (*version < 2.0) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Unsupported file version %.2f in file id: %d",
             *version, exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (nc_get_att_int(exoid, NC_GLOBAL, ATT_FLT_WORDSIZE, &file_wordsize) !=
      NC_NOERR) { /* try old (prior to db version 2.02) attribute name */
    if ((status = nc_get_att_int(exoid, NC_GLOBAL, ATT_FLT_WORDSIZE_BLANK, &file_wordsize)) !=
        NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get file wordsize from file id: %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }

  /* See if int64 status attribute exists and if so, what data is stored as
   * int64
   * Older files don't have the attribute, so it is not an error if it is
   * missing
   */
  if (nc_get_att_int(exoid, NC_GLOBAL, ATT_INT64_STATUS, &int64_status) != NC_NOERR) {
    int64_status = 0; /* Just in case it gets munged by a failed get_att_int call */
  }

  /* Merge in API int64 status flags as specified by caller of function... */
  int64_status |= (mode & EX_ALL_INT64_API);

  /* Verify that there is not an existing file_item struct for this
     exoid This could happen (and has) when application calls
     ex_open(), but then closes file using nc_close() and then reopens
     file.  NetCDF will possibly reuse the exoid which results in
     internal corruption in exodus data structures since exodus does
     not know that file was closed and possibly new file opened for
     this exoid
  */
  if (ex_find_file_item(exoid) != NULL) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: There is an existing file already using the file "
             "id %d which was also assigned to file %s.\n\tWas "
             "nc_close() called instead of ex_close() on an open Exodus "
             "file?\n",
             exoid, path);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    nc_close(exoid);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* initialize floating point and integer size conversion. */
  if (ex_conv_ini(exoid, comp_ws, io_ws, file_wordsize, int64_status, 1, is_hdf5, is_pnetcdf) !=
      EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to initialize conversion routines in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_LASTERR);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  EX_FUNC_LEAVE(exoid);
}
#else
/*
 * Prevent warning in some versions of ranlib(1) because the object
 * file has no symbols.
 */
const char exodus_unused_symbol_dummy_1;
#endif
