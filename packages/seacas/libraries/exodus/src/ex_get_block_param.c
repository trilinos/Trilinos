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
 *
 * \undoc exgblk - read block parameters
 *
 * entry conditions -
 *   input parameters:
 *       int     idexo                   exodus file id
 *       int     blk_type                block type (edge,face,element)
 *       int     blk_id                  block id
 *
 * exit conditions -
 *       char*   elem_type               element type name
 *       int*    num_entries_this_blk    number of elements in this element block
 *       int*    num_nodes_per_entry     number of nodes per element block
 *       int*    num_attr_per_entry      number of attributes
 *
 * revision history -
 *
 *
 */

#include "exodusII.h"     // for ex_block, ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, ATT_NAME_ELB, etc
#include <inttypes.h>     // for PRId64
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <string.h> // for strcpy

/*
 * reads the parameters used to describe an edge, face, or element block
 */

int ex_get_block_param(int exoid, ex_block *block)
{
  int         dimid, connid, blk_id_ndx;
  size_t      len, i;
  char        errmsg[MAX_ERR_LENGTH];
  int         status;
  const char *dnument = NULL;
  const char *dnumnod = NULL;
  const char *dnumedg = NULL;
  const char *dnumfac = NULL;
  const char *dnumatt = NULL;
  const char *ablknam = NULL;
  const char *vblkcon = NULL;

  struct ex_file_item *file = NULL;

  EX_FUNC_ENTER();

  file = ex_find_file_item(exoid);
  if (!file) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d in ex_get_block_param().", exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* First, locate index of element block id in VAR_ID_EL_BLK array */
  blk_id_ndx = ex_id_lkup(exoid, block->type, block->id);
  if (blk_id_ndx <= 0) {
    ex_get_err(NULL, NULL, &status);
    if (status != 0) {
      strcpy(block->topology, "NULL"); /* NULL element type name */
      block->num_entry           = 0;  /* no elements            */
      block->num_nodes_per_entry = 0;  /* no nodes               */
      block->num_edges_per_entry = 0;
      block->num_faces_per_entry = 0;
      block->num_attribute       = 0; /* no attributes          */
      if (status == EX_NULLENTITY) {  /* NULL element block?    */
        EX_FUNC_LEAVE(EX_NOERR);
      }
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate %s id  %" PRId64 " in id array in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }

  switch (block->type) {
  case EX_EDGE_BLOCK:
    dnument = DIM_NUM_ED_IN_EBLK(blk_id_ndx);
    dnumnod = DIM_NUM_NOD_PER_ED(blk_id_ndx);
    dnumedg = 0;
    dnumfac = 0;
    dnumatt = DIM_NUM_ATT_IN_EBLK(blk_id_ndx);
    vblkcon = VAR_EBCONN(blk_id_ndx);
    ablknam = ATT_NAME_ELB;
    break;
  case EX_FACE_BLOCK:
    dnument = DIM_NUM_FA_IN_FBLK(blk_id_ndx);
    dnumnod = DIM_NUM_NOD_PER_FA(blk_id_ndx);
    dnumedg = 0; /* it is possible this might be non-NULL some day */
    dnumfac = 0;
    dnumatt = DIM_NUM_ATT_IN_FBLK(blk_id_ndx);
    vblkcon = VAR_FBCONN(blk_id_ndx);
    ablknam = ATT_NAME_ELB;
    break;
  case EX_ELEM_BLOCK:
    dnument = DIM_NUM_EL_IN_BLK(blk_id_ndx);
    dnumnod = DIM_NUM_NOD_PER_EL(blk_id_ndx);
    dnumedg = DIM_NUM_EDG_PER_EL(blk_id_ndx);
    dnumfac = DIM_NUM_FAC_PER_EL(blk_id_ndx);
    dnumatt = DIM_NUM_ATT_IN_BLK(blk_id_ndx);
    vblkcon = VAR_CONN(blk_id_ndx);
    ablknam = ATT_NAME_ELB;
    break;
  default:
    snprintf(errmsg, MAX_ERR_LENGTH, "Bad block type parameter (%d) specified for file id %d.",
             block->type, exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire values of some dimensions */
  if ((status = nc_inq_dimid(exoid, dnument, &dimid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate number of entities in %s  %" PRId64 " in file id %d",
             ex_name_of_object(block->type), block->id, exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &len)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to get number of %ss in block  %" PRId64 " in file id %d",
             ex_name_of_object(block->type), block->id, exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }
  block->num_entry = len;

  if ((status = nc_inq_dimid(exoid, dnumnod, &dimid)) != NC_NOERR) {
    /* undefined => no node entries per element */
    len = 0;
  }
  else {
    if ((status = nc_inq_dimlen(exoid, dimid, &len)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to get number of nodes/entity in %s  %" PRId64 " in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }
  block->num_nodes_per_entry = len;

  if (!file->has_edges || block->type != EX_ELEM_BLOCK) {
    block->num_edges_per_entry = 0;
  }
  else {
    if ((status = nc_inq_dimid(exoid, dnumedg, &dimid)) != NC_NOERR) {
      /* undefined => no edge entries per element */
      len = 0;
    }
    else {
      if ((status = nc_inq_dimlen(exoid, dimid, &len)) != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to get number of edges/entry in %s  %" PRId64 " in file id %d",
                 ex_name_of_object(block->type), block->id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }
    block->num_edges_per_entry = len;
  }

  if (!file->has_faces || block->type != EX_ELEM_BLOCK) {
    block->num_faces_per_entry = 0;
  }
  else {
    if ((status = nc_inq_dimid(exoid, dnumfac, &dimid)) != NC_NOERR) {
      /* undefined => no face entries per element */
      len = 0;
    }
    else {
      if ((status = nc_inq_dimlen(exoid, dimid, &len)) != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to get number of faces/entity in %s  %" PRId64 " in file id %d",
                 ex_name_of_object(block->type), block->id, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }
    block->num_faces_per_entry = len;
  }

  if ((status = nc_inq_dimid(exoid, dnumatt, &dimid)) != NC_NOERR) {
    /* dimension is undefined */
    block->num_attribute = 0;
  }
  else {
    if ((status = nc_inq_dimlen(exoid, dimid, &len)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to get number of attributes in %s  %" PRId64 " in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
    block->num_attribute = len;
  }

  if (block->num_nodes_per_entry == 0 && block->num_edges_per_entry == 0 &&
      block->num_faces_per_entry == 0) {
    vblkcon = NULL;
  }

  if (block->num_nodes_per_entry > 0) {
    ; /* Do nothing, vblkcon should be correctly set already */
  }
  else if (block->num_edges_per_entry > 0) {
    vblkcon = VAR_EBCONN(blk_id_ndx);
  }
  else if (block->num_faces_per_entry > 0) {
    vblkcon = VAR_FCONN(blk_id_ndx);
  }

  if (vblkcon) {
    /* look up connectivity array for this element block id */
    if ((status = nc_inq_varid(exoid, vblkcon, &connid)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate connectivity array for %s  %" PRId64 " in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if ((status = nc_inq_attlen(exoid, connid, ablknam, &len)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get %s  %" PRId64 " type in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if (len > (MAX_STR_LENGTH + 1)) {
      len = MAX_STR_LENGTH;
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Warning: %s  %" PRId64 " type will be truncated to %ld chars",
               ex_name_of_object(block->type), block->id, (long)len);
      ex_err_fn(exoid, __func__, errmsg, EX_MSG);
    }

    for (i = 0; i < MAX_STR_LENGTH + 1; i++) {
      block->topology[i] = '\0';

      /* get the element type name */
    }
    if ((status = nc_get_att_text(exoid, connid, ablknam, block->topology)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get %s  %" PRId64 " type in file id %d",
               ex_name_of_object(block->type), block->id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    /* get rid of trailing blanks */
    ex_trim_internal(block->topology);
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
