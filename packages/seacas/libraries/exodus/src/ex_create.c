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

/*!
\ingroup Utilities

\note The ex_create_int() is an internal function called by
ex_create(). The user should call ex_create() and not ex_create_int().

The function ex_create() creates a new exodus file and returns an ID
that can subsequently be used to refer to the file.

All floating point values in an exodus file are stored as either
4-byte (float) or 8-byte (double) numbers; no mixing of 4- and
8-byte numbers in a single file is allowed. An application code can
compute either 4- or 8-byte values and can designate that the values
be stored in the exodus file as either 4- or 8-byte numbers;
conversion between the 4- and 8-byte values is performed automatically
by the API routines. Thus, there are four possible combinations of
compute word size and storage (or I/O) word size.

\return In case of an error, ex_create() returns a negative number. Possible
causes of errors include:
  -  Passing a file name that includes a directory that does not
 exist.
  -  Specifying a file name of a file that exists and also
 specifying a no clobber option.
  -  Attempting to create a file in a directory without permission
 to create files there.
  -  Passing an invalid file clobber mode.


\param path The file name of the new exodus file. This can be given as either an
            absolute path name (from the root of the file system) or a relative
            path name (from the current directory).

\param cmode Mode. Use one of the following predefined constants:
\arg #EX_NOCLOBBER  To create the new file only if the given file name does
not refer to a
                      file that already exists.

\arg #EX_CLOBBER    To create the new file, regardless of whether a file with
the same
                      name already exists. If a file with the same name does
exist, its
                      contents will be erased.

\arg #EX_64BIT_OFFSET To create a model that can store individual datasets
larger than
                        2 gigabytes. This modifies the internal storage used by
exodusII and
                        also puts the underlying NetCDF file into the \e 64-bit
offset'
                        mode. See largemodel for more details on this
                        mode. A large model file will also be created if the
                        environment variable EXODUS_LARGE_MODEL is defined
                        in the users environment. A message will be printed to
standard output
                        if this environment variable is found. #EX_LARGE_MODEL is
alias.

\arg #EX_NORMAL_MODEL Create a standard model.

\arg #EX_64BIT_DATA      To create a model using the CDF5 format which uses the
                        classic model but has 64-bit dimensions and sizes.
                        This type will also be created if the
                        environment variable EXODUS_NETCDF5 is defined in the
                        users environment. A message will be printed to standard
                        output if
                        this environment variable is found.

\arg #EX_NETCDF4 To create a model using the HDF5-based NetCDF-4
                        output. An HDF5-based NetCDF-4 file will also be created
if the
                        environment variable EXODUS_NETCDF4 is defined in the
                        users environment. A message will be printed to standard
output if
                        this environment variable is found.

\arg #EX_NOSHARE Do not open the underlying NetCDF file in \e share
mode. See the
                        NetCDF documentation for more details.

\arg #EX_SHARE   Do open the underlying NetCDF file in \e share mode. See
the NetCDF
                        documentation for more details.

\param[in,out] comp_ws  The word size in bytes (0, 4 or 8) of the floating point
variables
                        used in the application program. If 0 (zero) is passed,
the default
                        sizeof(float) will be used and returned in this
variable. WARNING: all
                        exodus functions requiring floats must be passed floats
declared with
                        this passed in or returned compute word size (4 or 8).}

\param io_ws            The word size in bytes (4 or 8) of the floating point
                        data as they are to be stored in the exodus file.

\param run_version (internally generated) used to verify compatibility of library
and include files.

The following code segment creates an exodus file called \file{test.exo}:

~~~{.c}
#include "exodusII.h"
int CPU_word_size, IO_word_size, exoid;
CPU_word_size = sizeof(float);      \comment{use float or double}
IO_word_size = 8;                   \comment{store variables as doubles}

\comment{create exodus file}
exoid = ex_create ("test.exo"       \comment{filename path}
                    EX_CLOBBER,     \comment{create mode}
                    &CPU_word_size, \comment{CPU float word size in bytes}
                    &IO_word_size); \comment{I/O float word size in bytes}
~~~

*/
#include "exodusII.h"
#include "exodusII_int.h"

/* NOTE: Do *not* call `ex_create_int()` directly.  The public API
 *       function name is `ex_create()` which is a wrapper that calls
 *       `ex_create_int` with an additional argument to make sure
 *       library and include file are consistent
 */
int ex_create_int(const char *path, int cmode, int *comp_ws, int *io_ws, int run_version)
{
  int  exoid;
  int  status;
  char errmsg[MAX_ERR_LENGTH];
  int  nc_mode = 0;

  unsigned int my_mode     = cmode;
  int          is_parallel = 0;

  EX_FUNC_ENTER();

  nc_mode = ex__handle_mode(my_mode, is_parallel, run_version);

  if ((status = nc_create(path, nc_mode, &exoid)) != NC_NOERR) {
#if NC_HAS_HDF5
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: file create failed for %s", path);
#else
    if (my_mode & EX_NETCDF4) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: file create failed for %s in NETCDF4 "
               "mode.\n\tThis library does not support netcdf-4 files.",
               path);
    }
    else {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: file create failed for %s", path);
    }
#endif
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  status = ex__populate_header(exoid, path, my_mode, is_parallel, comp_ws, io_ws);
  if (status != EX_NOERR) {
    EX_FUNC_LEAVE(status);
  }

  EX_FUNC_LEAVE(exoid);
}
