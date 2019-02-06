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
#include "exodusII_int.h" // for EX_FATAL, ex_id_lkup, etc
#include <inttypes.h>     // for PRId64
#include <stddef.h>       // for size_t
#include <stdio.h>

/*! write out the connectivity array */
#define EX_WRITE_CONN(TNAME, VARCONN, VARCONNVAL)                                                  \
  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {                                                \
    status = nc_put_var_longlong(exoid, VARCONN, VARCONNVAL);                                      \
  }                                                                                                \
  else {                                                                                           \
    status = nc_put_var_int(exoid, VARCONN, VARCONNVAL);                                           \
  }                                                                                                \
  if (status != NC_NOERR) {                                                                        \
    snprintf(errmsg, MAX_ERR_LENGTH,                                                               \
             "ERROR: failed to write connectivity array for %s block %" PRId64 " in file id %d",   \
             TNAME, blk_id, exoid);                                                                \
    ex_err_fn(exoid, __func__, errmsg, status);                                                    \
    EX_FUNC_LEAVE(EX_FATAL);                                                                       \
  }

/*!
 * writes the connectivity array for a block
 * \param exoid           exodus file id
 * \param blk_type        type of block
 * \param blk_id          id of block
 * \param node_conn       node-element connectivity
 * \param elem_edge_conn  element-edge connectivity (NULL if none)
 * \param elem_face_conn  element-face connectivity (NULL if none)
 */

int ex_put_conn(int exoid, ex_entity_type blk_type, ex_entity_id blk_id, const void_int *node_conn,
                const void_int *elem_edge_conn, const void_int *elem_face_conn)
{
  int  connid = -1, blk_id_ndx, status;
  char errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  ex_check_valid_file_id(exoid, __func__);

  blk_id_ndx = ex_id_lkup(exoid, blk_type, blk_id);
  if (blk_id_ndx <= 0) {
    ex_get_err(NULL, NULL, &status);

    if (status != 0) {
      if (status == EX_NULLENTITY) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "Warning: connectivity array not allowed for NULL %s %" PRId64 " in file id %d",
                 ex_name_of_object(blk_type), blk_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, EX_NULLENTITY);
        EX_FUNC_LEAVE(EX_WARN);
      }
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate %s id %" PRId64 " in id array in file id %d",
               ex_name_of_object(blk_type), blk_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }

  /* inquire id's of previously defined dimensions  */
  if (node_conn) {
    switch (blk_type) {
    case EX_ELEM_BLOCK: status = nc_inq_varid(exoid, VAR_CONN(blk_id_ndx), &connid); break;
    case EX_FACE_BLOCK: status = nc_inq_varid(exoid, VAR_FBCONN(blk_id_ndx), &connid); break;
    case EX_EDGE_BLOCK: status = nc_inq_varid(exoid, VAR_EBCONN(blk_id_ndx), &connid); break;
    default:
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Internal ERROR: unrecognized block type in switch: %d in file id %d", blk_type,
               exoid);
      ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
    }
    if (status != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate connectivity array for %s %" PRId64 " in file id %d",
               ex_name_of_object(blk_type), blk_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    EX_WRITE_CONN(ex_name_of_object(blk_type), connid, node_conn);
  }

  /* If there are edge and face connectivity arrays that belong with the element
   * block, write them now. Warn if they are required but not specified or
   * specified but not required.
   */
  if (blk_type == EX_ELEM_BLOCK) {
    int    nedpereldim, nfapereldim;
    size_t num_ed_per_elem, num_fa_per_elem;

    if (elem_edge_conn != 0) {
      status = nc_inq_dimid(exoid, DIM_NUM_EDG_PER_EL(blk_id_ndx), &nedpereldim);
      if (status != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: edge connectivity specified but failed to "
                 "locate number of edges/element in block %" PRId64 " in file id %d",
                 blk_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }

    if (elem_face_conn != 0) {
      status = nc_inq_dimid(exoid, DIM_NUM_FAC_PER_EL(blk_id_ndx), &nfapereldim);
      if (status != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: face connectivity specified but failed to "
                 "locate number of faces/element in block %" PRId64 " in file id %d",
                 blk_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }

    num_ed_per_elem = 0;
    if ((elem_edge_conn != 0) &&
        ((status = nc_inq_dimlen(exoid, nedpereldim, &num_ed_per_elem)) != NC_NOERR)) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to get number of edges/elem in block %" PRId64 " in file id %d",
               blk_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    num_fa_per_elem = 0;
    if ((elem_face_conn != 0) &&
        ((status = nc_inq_dimlen(exoid, nfapereldim, &num_fa_per_elem)) != NC_NOERR)) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to get number of faces/elem in block %" PRId64 " in file id %d",
               blk_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if ((num_ed_per_elem == 0 && elem_edge_conn != 0) ||
        (num_ed_per_elem != 0 && elem_edge_conn == 0)) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: number of edges per element (%ld) doesn't "
               "agree with elem_edge_conn (0x%p)",
               (long)num_ed_per_elem, (void *)elem_edge_conn);
      ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if ((num_fa_per_elem == 0 && elem_face_conn != 0) ||
        (num_fa_per_elem != 0 && elem_face_conn == 0)) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: number of faces per element (%ld) doesn't "
               "agree with elem_face_conn (0x%p)",
               (long)num_fa_per_elem, (void *)elem_face_conn);
      ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if (num_ed_per_elem != 0) {
      status = nc_inq_varid(exoid, VAR_ECONN(blk_id_ndx), &connid);
      if (status != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to locate connectivity array for "
                 "element edge block %" PRId64 " in file id %d",
                 blk_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
      EX_WRITE_CONN("element edge", connid, elem_edge_conn);
    }

    if (num_fa_per_elem != 0) {
      status = nc_inq_varid(exoid, VAR_FCONN(blk_id_ndx), &connid);
      if (status != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to locate connectivity array for "
                 "element face block %" PRId64 " in file id %d",
                 blk_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
      EX_WRITE_CONN("element face", connid, elem_face_conn);
    }
  }

  EX_FUNC_LEAVE(EX_NOERR);
}
