/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 *
 *      ex_get_elem_type()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      elem_blk_id        - The element block ID to find the element type
 *                           for.
 *      elem_type          - The returned name of the element.
 *
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include "exodusII.h"     // for ex_err, exerrval, etc
#include "exodusII_int.h" // for EX_FATAL, ATT_NAME_ELB, etc
#include "netcdf.h"       // for NC_NOERR, nc_get_att_text, etc
#include <inttypes.h>     // for PRId64
#include <stddef.h>       // for size_t
#include <stdio.h>

int ex_get_elem_type(int exoid, ex_entity_id elem_blk_id, char *elem_type)
/*
 *      Reads the element type for a specific element block
 *           elem_type is assumed to have a length of MAX_STR_LENGTH+1
 */
{
  const char *func_name = "ex_get_elem_type";

  int    connid, el_blk_id_ndx, status;
  size_t len;
  char   errmsg[MAX_ERR_LENGTH];

  /*****************************************************************************/

  ex_check_valid_file_id(exoid);

  /* inquire id's of previously defined dimensions */
  if ((el_blk_id_ndx = ex_id_lkup(exoid, EX_ELEM_BLOCK, elem_blk_id)) == -1) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to find element block ID %" PRId64 " in file %d", elem_blk_id, exoid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(exoid, VAR_CONN(el_blk_id_ndx), &connid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find connectivity variable in file ID %d",
             exoid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* get the element type name */
  if ((status = nc_inq_attlen(exoid, connid, ATT_NAME_ELB, &len)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find attribute in file ID %d", exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (len > (MAX_STR_LENGTH + 1)) {
    exerrval = EX_MSG;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Element type must be of length %d in file ID %d",
             (int)len, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Make sure the end of the string is terminated with a null character */
  elem_type[MAX_STR_LENGTH] = '\0';

  if ((status = nc_get_att_text(exoid, connid, ATT_NAME_ELB, elem_type)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get attribute \"%s\" in file ID %d",
             ATT_NAME_ELB, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }
  return (EX_NOERR);
}
