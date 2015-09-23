/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
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
/*****************************************************************************
*
* exgsnl - ex_get_side_set_node_list_len
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*       int     side_set_id             side set id
*
* exit conditions - 
*       int     *side_set_node_list_len length of node list
*
* revision history - 
*
*
*****************************************************************************/

#include <ctype.h>                      // for toupper
#include <inttypes.h>                   // for PRId64
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf
#include <stdlib.h>                     // for malloc, NULL, free
#include <string.h>                     // for strncmp, strlen
#include <sys/types.h>                  // for int64_t
#include "exodusII.h"                   // for ex_err, exerrval, etc
#include "exodusII_int.h"               // for elem_blk_parm, EX_FATAL, etc
#include "netcdf.h"                     // for NC_NOERR

/*!
 * This routine is designed to read the Exodus II V 2.0 side set side 
 * definition  and return the length of a ExodusI style side set node list.
 * \param           exoid                   exodus file id
 * \param           side_set_id             side set id
 * \param[out]     *side_set_node_list_len length of node list
 */

int ex_get_side_set_node_list_len(int exoid,
				  ex_entity_id side_set_id,
				  void_int *side_set_node_list_len)
{
  size_t i, j;
  size_t m;
  int64_t num_side_sets, num_elem_blks, num_df, ndim;
  size_t list_len = 0;
  int64_t tot_num_elem = 0, tot_num_ss_elem = 0; 
  void_int *elem_blk_ids;
  int *ss_elem_ndx = NULL;
  int64_t *ss_elem_ndx_64 = NULL;
  
  void_int *side_set_elem_list;
  void_int *side_set_side_list;
  int elem_ctr; 

  int err_stat = EX_NOERR;
  int status;
  
  struct elem_blk_parm  *elem_blk_parms;

  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

    if (ex_int64_status(exoid) & EX_BULK_INT64_API)
      *(int64_t*)side_set_node_list_len = 0; /* default value */
    else
      *(int*)side_set_node_list_len = 0; /* default value */
      
  /* first check if any side sets are specified */
  /* inquire how many side sets have been stored */

  /* get the dimensionality of the coordinates;  this is necessary to
     distinguish between 2d TRIs and 3d TRIs */

  ndim = ex_inquire_int(exoid, EX_INQ_DIM);
  if (ndim < 0)  {
    sprintf(errmsg,
           "Error: failed to get dimensionality in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    return(EX_FATAL);
  }

  tot_num_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
  if (tot_num_elem < 0) {
    sprintf(errmsg,
           "Error: failed to get total number of elements in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    return(EX_FATAL);
  }

  num_elem_blks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
  if (num_elem_blks < 0) {
    sprintf(errmsg,
           "Error: failed to get number of element blocks in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    return(EX_FATAL);
  }

  num_side_sets = ex_inquire_int(exoid, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    sprintf(errmsg,
           "Error: failed to get number of side sets in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    return(EX_FATAL);
  }

  if (num_side_sets == 0) {
    sprintf(errmsg,
           "Warning: no side sets defined in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,EX_WARN);
    return(EX_WARN);
  }

  /* First determine the  # of elements in the side set*/
  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    status = ex_get_side_set_param(exoid,side_set_id,&tot_num_ss_elem,&num_df);
  } else {
    int tot;
    int df;
    status = ex_get_side_set_param(exoid,side_set_id,&tot,&df);
    tot_num_ss_elem = tot;
    num_df = df;
  }

  if (status != NC_NOERR) {
    sprintf(errmsg,
         "Error: failed to get number of elements in side set %"PRId64" in file id %d",
            side_set_id, exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    return(EX_FATAL);
  }

  if (tot_num_ss_elem == 0) /* NULL side set? */
    return (EX_NOERR); /* return zero */

  /* Minor optimization/kluge -- If num_df is nonzero, or 1 per face
     then assume that it matches the number of nodes in the sideset... */
  if (num_df > 0 && num_df != tot_num_ss_elem) {
    if (ex_int64_status(exoid) & EX_BULK_INT64_API)
      *(int64_t*)side_set_node_list_len = num_df;
    else
      *(int*)side_set_node_list_len = num_df;
    return(EX_NOERR);
  }

  /* Allocate space for the side set element list */
  {
    int int_size = sizeof(int);
    if (ex_int64_status(exoid) & EX_BULK_INT64_API)
      int_size = sizeof(int64_t);
    if (!(side_set_elem_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set element list for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
      return (EX_FATAL);
    }

    /* Allocate space for the side set side list */
    if (!(side_set_side_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set side list for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
      err_stat = EX_FATAL;
      goto cleanup;
    }

    if (ex_get_side_set(exoid, side_set_id, 
			side_set_elem_list, side_set_side_list) != NC_NOERR) {
      sprintf(errmsg,
	      "Error: failed to get side set %"PRId64" in file id %d",
	      side_set_id, exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
      err_stat = EX_FATAL;
      goto cleanup;
    }
    
    /* Allocate space for the ss element index array */
    if (int_size == sizeof(int64_t)) {
      ss_elem_ndx_64=malloc(tot_num_ss_elem*int_size);
    } else {
      ss_elem_ndx   =malloc(tot_num_ss_elem*int_size);
    }

    if (ss_elem_ndx_64==NULL && ss_elem_ndx == NULL) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set elem sort array for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
      err_stat = EX_FATAL;
      goto cleanup;
    }
  }

  /* Sort side set element list into index array  - non-destructive */
  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    for (i=0;i<tot_num_ss_elem;i++)
      ss_elem_ndx_64[i] = i; /* init index array to current position */
    ex_iqsort64(side_set_elem_list, ss_elem_ndx_64,tot_num_ss_elem);
  } else {
    for (i=0;i<tot_num_ss_elem;i++)
      ss_elem_ndx[i] = i; /* init index array to current position */
    ex_iqsort(side_set_elem_list, ss_elem_ndx,tot_num_ss_elem);
  }


  /* Allocate space for the element block ids */
  {
    int int_size = sizeof(int);
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      int_size = sizeof(int64_t);
    }

    if (!(elem_blk_ids=malloc(num_elem_blks*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for element block ids for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
      err_stat = EX_FATAL;
      goto cleanup;
    }
  }
  
  if (ex_get_elem_blk_ids(exoid, elem_blk_ids)) {
    sprintf(errmsg,
	    "Error: failed to get element block ids in file id %d",
            exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,EX_MSG);
    err_stat = EX_FATAL;
    goto cleanup;
  } 

  /* Allocate space for the element block params */
  if (!(elem_blk_parms=malloc(num_elem_blks*sizeof(struct elem_blk_parm)))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
      "Error: failed to allocate space for element block params for file id %d",
            exoid);
    ex_err("ex_get_side_set_node_list_len",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }

  elem_ctr = 0;
  for (i=0; i<num_elem_blks; i++) {
    ex_entity_id id;
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      id = ((int64_t*)elem_blk_ids)[i];
    } else {
      id = ((int*)elem_blk_ids)[i];
    }

    err_stat = ex_int_get_block_param(exoid, id, ndim, &elem_blk_parms[i]);
    if (err_stat != EX_NOERR) {
      goto cleanup;
    }

    elem_ctr += elem_blk_parms[i].num_elem_in_blk;
    elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */
  }

/* Walk through element list and keep a running count of the node length */

  list_len = 0;
  for (i=0;i<tot_num_ss_elem;i++)
  {
    size_t elem;
    size_t side;
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      elem = ((int64_t*)side_set_elem_list)[i];
      side = ((int64_t*)side_set_side_list)[i];
    } else {
      elem = ((int*)side_set_elem_list)[i];
      side = ((int*)side_set_side_list)[i];
    }

    for (j=0; j<num_elem_blks; j++)
    {
      if (elem_blk_parms[j].elem_type_val != EX_EL_NULL_ELEMENT)
        if (elem <= elem_blk_parms[j].elem_ctr)
          break; /* stop because we found the parameters for this element */
    }
    if (j >= num_elem_blks)
    {
      exerrval = EX_BADPARAM;
      sprintf(errmsg,
             "Error: Invalid element number %"ST_ZU" found in side set %"PRId64" in file %d",
              elem, side_set_id, exoid);
      ex_err("ex_get_side_set_node_list_len",errmsg,EX_MSG);
      err_stat = EX_FATAL;
      goto cleanup;
    }
    list_len += elem_blk_parms[j].num_nodes_per_side[side-1];
  }

  if (ex_int64_status(exoid) & EX_BULK_INT64_API)
    *(int64_t*)side_set_node_list_len = list_len;
  else
    *(int*)side_set_node_list_len = list_len;

  /* All done: release element block ids array,
     element block parameters array, and side set element index array */
 cleanup:
  ex_safe_free(elem_blk_ids);
  ex_safe_free(elem_blk_parms);
  ex_safe_free(ss_elem_ndx);
  ex_safe_free(ss_elem_ndx_64);
  ex_safe_free(side_set_side_list);
  ex_safe_free(side_set_elem_list);

  return(err_stat);
}

