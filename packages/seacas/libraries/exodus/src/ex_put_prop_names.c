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

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, etc
#include "netcdf.h"       // for NC_NOERR, nc_enddef, etc
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <string.h> // for strlen

/*!

The function ex_put_prop_names() writes object property names and
allocates space for object property arrays used to assign integer
properties to element blocks, node sets, or side sets. The property
arrays are initialized to zero (0). Although this function is
optional, since ex_put_prop() will allocate space within the data file
if it hasn't been previously allocated, it is more efficient to use
ex_put_prop_names() if there is more than one property to store.

\return In case of an error, ex_put_prop_names() returns a negative number; a
warning will return a positive number.  Possible causes of errors
include:
  -  data file not properly opened with call to ex_create() or ex_open()
  -  data file opened for read only.
  -  data file not initialized properly with call to ex_put_init().
  -  invalid object type specified.
  -  no object of the specified type is stored in the file.

\param[in] exoid       exodus file ID returned from a previous call to
ex_create() or ex_open().
\param[in] obj_type    Type of object; use one of the options in the table
below.
\param[in] num_props   The number of integer properties to be assigned to all of
the objects
                       of the type specified (element blocks, node sets, or side
sets).
\param[in] prop_names  Array containing  num_props names (of maximum length
                       of \p MAX_STR_LENGTH ) of properties to be stored.

| ex_entity_type | description               |
| -------------- | ------------------------- |
|  EX_NODE_SET   |  Node Set entity type     |
|  EX_EDGE_BLOCK |  Edge Block entity type   |
|  EX_EDGE_SET   |  Edge Set entity type     |
|  EX_FACE_BLOCK |  Face Block entity type   |
|  EX_FACE_SET   |  Face Set entity type     |
|  EX_ELEM_BLOCK |  Element Block entity type|
|  EX_ELEM_SET   |  Element Set entity type  |
|  EX_SIDE_SET   |  Side Set entity type     |
|  EX_ELEM_MAP   |  Element Map entity type  |
|  EX_NODE_MAP   |  Node Map entity type     |
|  EX_EDGE_MAP   |  Edge Map entity type     |
|  EX_FACE_MAP   |  Face Map entity type     |

For instance, suppose a user wanted to assign the 1st, 3rd, and 5th
element blocks (those element blocks stored 1st, 3rd, and 5th,
regardless of their ID) to a group (property) called \b TOP, and the
2nd, 3rd, and 4th element blocks to a group called \b LSIDE. This
could be accomplished with the following code:

~~~{.c}

char* prop_names[2];
int top_part[]   = {1,0,1,0,1};
int lside_part[] = {0,1,1,1,0};

int id[] = {10, 20, 30, 40, 50};

prop_names[0] = "TOP";
prop_names[1] = "LSIDE";

\comment{This call to ex_put_prop_names is optional, but more efficient}
ex_put_prop_names (exoid, EX_ELEM_BLOCK, 2, prop_names);

\comment{The property values can be output individually thus}
for (i=0; i < 5; i++) {
   ex_put_prop (exoid, EX_ELEM_BLOCK, id[i], prop_names[0],
                top_part[i]);

   ex_put_prop (exoid, EX_ELEM_BLOCK, id[i], prop_names[1],
                lside_part[i]);
}

\comment{Alternatively, the values can be output as an array}
ex_put_prop_array (exoid, EX_ELEM_BLOCK, prop_names[0],
                   top_part);

ex_put_prop_array (exoid, EX_ELEM_BLOCK, prop_names[1],
                   lside_part);

~~~

*/

int ex_put_prop_names(int exoid, ex_entity_type obj_type, int num_props, char **prop_names)
{
  int       status;
  int       oldfill, temp;
  int       i, propid, dimid, dims[1];
  size_t    name_length, prop_name_len;
  char *    name;
  long long vals[1];
  int       max_name_len = 0;
  int       int_type     = NC_INT;

  char errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  ex_check_valid_file_id(exoid, __func__);

  if (ex_int64_status(exoid) & EX_IDS_INT64_DB) {
    int_type = NC_INT64;
  }

  /* Get the name string length */
  name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH) + 1;

  /* inquire id of previously defined dimension (number of objects) */
  if ((status = nc_inq_dimid(exoid, ex_dim_num_objects(obj_type), &dimid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate number of %s in file id %d",
             ex_name_of_object(obj_type), exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  nc_set_fill(exoid, NC_FILL, &oldfill); /* fill with zeros per routine spec */

  /* put netcdf file into define mode  */
  if ((status = nc_redef(exoid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to place file id %d into define mode", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* define num_props variables; we postpend the netcdf variable name with  */
  /* a counter starting at 2 because "xx_prop1" is reserved for the id array*/
  dims[0] = dimid;

  for (i = 0; i < num_props; i++) {
    switch (obj_type) {
    case EX_ELEM_BLOCK: name = VAR_EB_PROP(i + 2); break;
    case EX_FACE_BLOCK: name = VAR_FA_PROP(i + 2); break;
    case EX_EDGE_BLOCK: name = VAR_ED_PROP(i + 2); break;
    case EX_NODE_SET: name = VAR_NS_PROP(i + 2); break;
    case EX_SIDE_SET: name = VAR_SS_PROP(i + 2); break;
    case EX_EDGE_SET: name = VAR_ES_PROP(i + 2); break;
    case EX_FACE_SET: name = VAR_FS_PROP(i + 2); break;
    case EX_ELEM_SET: name = VAR_ELS_PROP(i + 2); break;
    case EX_ELEM_MAP: name = VAR_EM_PROP(i + 2); break;
    case EX_FACE_MAP: name = VAR_FAM_PROP(i + 2); break;
    case EX_EDGE_MAP: name = VAR_EDM_PROP(i + 2); break;
    case EX_NODE_MAP: name = VAR_NM_PROP(i + 2); break;
    default:
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported; file id %d", obj_type,
               exoid);
      ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
      goto error_ret; /* Exit define mode and return */
    }

    if ((status = nc_def_var(exoid, name, int_type, 1, dims, &propid)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to create property array variable in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      goto error_ret; /* Exit define mode and return */
    }

    vals[0] = 0; /* fill value */

    /*   create attribute to cause variable to fill with zeros per routine spec
     */
    if ((status = nc_put_att_longlong(exoid, propid, _FillValue, int_type, 1, vals)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to create property name fill attribute in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      goto error_ret; /* Exit define mode and return */
    }

    /*   Check that the property name length is less than MAX_NAME_LENGTH */
    prop_name_len = strlen(prop_names[i]) + 1;
    if (prop_name_len > name_length) {
      fprintf(stderr,
              "Warning: The property name '%s' is too long.\n\tIt will "
              "be truncated from %d to %d characters\n",
              prop_names[i], (int)prop_name_len - 1, (int)name_length - 1);
      prop_name_len = name_length;
    }

    if (prop_name_len > max_name_len) {
      max_name_len = prop_name_len;
    }

    /*   store property name as attribute of property array variable */
    if ((status = nc_put_att_text(exoid, propid, ATT_PROP_NAME, prop_name_len, prop_names[i])) !=
        NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to store property name %s in file id %d",
               prop_names[i], exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      goto error_ret; /* Exit define mode and return */
    }
  }

  /* leave define mode  */
  if ((status = nc_enddef(exoid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to leave define mode in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* Update the maximum_name_length attribute on the file. */
  ex_update_max_name_length(exoid, max_name_len - 1);

  nc_set_fill(exoid, oldfill, &temp); /* default: turn off fill */
  EX_FUNC_LEAVE(EX_NOERR);

/* Fatal error: exit definition mode and return */
error_ret:
  if ((status = nc_enddef(exoid)) != NC_NOERR) { /* exit define mode */
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to complete definition for file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
  }
  EX_FUNC_LEAVE(EX_FATAL);
}
