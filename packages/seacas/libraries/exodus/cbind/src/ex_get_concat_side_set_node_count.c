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
* exgcssc - ex_get_concat_side_set_node_count
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*
* exit conditions - 
*       int     *side_set_node_cnt_list returned array of number of nodes for
*                                       side or face for all sidesets
* revision history - 
*
*****************************************************************************/

#include <assert.h>                     // for assert
#include <ctype.h>                      // for toupper
#include <inttypes.h>                   // for PRId64
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf
#include <stdlib.h>                     // for malloc, NULL, free
#include <string.h>                     // for strncmp, strlen
#include <sys/types.h>                  // for int64_t
#include "exodusII.h"                   // for ex_err, exerrval, EX_MSG, etc
#include "exodusII_int.h"               // for elem_blk_parm, EX_FATAL, etc

/*! \endcond */

/*! \undoc */
int ex_get_concat_side_set_node_count(int exoid,
				      int *side_set_node_cnt_list)
{
  int ii, i, j, iss, ioff; 
  ex_entity_id side_set_id;
  int  num_side_sets, num_elem_blks, num_df, ndim;
  int64_t tot_num_elem = 0, tot_num_ss_elem = 0, side, elem;
  void_int *elem_blk_ids       = NULL;
  void_int *side_set_ids       = NULL;
  void_int *ss_elem_ndx        = NULL;
  void_int *side_set_elem_list = NULL;
  void_int *side_set_side_list = NULL;
  size_t elem_ctr;
  int int_size, ids_size;
  int status;
  
  struct elem_blk_parm  *elem_blk_parms = NULL;

  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any side sets are specified */
  /* inquire how many side sets have been stored */
  num_side_sets = ex_inquire_int(exoid, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    sprintf(errmsg,
           "Error: failed to get number of side sets in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  if (num_side_sets == 0) {
    sprintf(errmsg,
           "Warning: no side sets defined in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_WARN);
    return(EX_WARN);
  }

  num_elem_blks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
  if (num_elem_blks < 0) {
    sprintf(errmsg,
           "Error: failed to get number of element blocks in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  tot_num_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
  if (tot_num_elem < 0) {
    sprintf(errmsg,
           "Error: failed to get total number of elements in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  /* get the dimensionality of the coordinates;  this is necessary to
     distinguish between 2d TRIs and 3d TRIs */
  ndim = ex_inquire_int(exoid, EX_INQ_DIM);
  if (ndim < 0) {
    sprintf(errmsg,
           "Error: failed to get dimensionality in file id %d",exoid);
    ex_err("ex_cvt_nodes_to_sides",errmsg,exerrval);
    return(EX_FATAL);
  }

  int_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_BULK_INT64_API)
    int_size = sizeof(int64_t);
  
  /* Allocate space for the element block ids */
  ids_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
    ids_size = sizeof(int64_t);
  }

  if (!(elem_blk_ids=malloc(num_elem_blks*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
            "Error: failed to allocate space for element block ids for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  if (ex_get_ids(exoid, EX_ELEM_BLOCK, elem_blk_ids) != EX_NOERR) {
    sprintf(errmsg,
            "Error: failed to get element block ids in file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
    goto error_ret;
  } 

  /* Allocate space for the element block params */
  if (!(elem_blk_parms=malloc(num_elem_blks*sizeof(struct elem_blk_parm)))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
      "Error: failed to allocate space for element block params for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  elem_ctr = 0;
  for (i=0; i<num_elem_blks; i++)
  {
    ex_entity_id id;
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      id = ((int64_t*)elem_blk_ids)[i];
    } else {
      id = ((int*)elem_blk_ids)[i];
    }

    if (ex_int_get_block_param(exoid, id, ndim, &elem_blk_parms[i]) != EX_NOERR) {
      goto error_ret;
    }

    elem_ctr += elem_blk_parms[i].num_elem_in_blk;
    elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */
  }


  /* Finally... Create the list of node counts for each face in the
   * side set.
   */
  /* Allocate space for the sideset ids */
  if (!(side_set_ids=malloc(num_side_sets*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
            "Error: failed to allocate space for side set ids for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  if (ex_get_ids(exoid, EX_SIDE_SET, side_set_ids) != EX_NOERR) {
    sprintf(errmsg,
            "Error: failed to get side set ids in file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
    goto error_ret;
  } 

  /* Lookup index of side set id in VAR_SS_IDS array */
  ioff = 0;
  for (iss=0; iss<num_side_sets; iss++) {
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      side_set_id = ((int64_t*)side_set_ids)[iss];
    } else {
      side_set_id = ((int*)side_set_ids)[iss];
    }
    
    /* First determine the  # of elements in the side set*/
    if (int_size == sizeof(int64_t)) {
      status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot_num_ss_elem,&num_df);
    } else {
      int tot, df;
      status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot,&df);
      tot_num_ss_elem = tot;
      num_df = df;
    }

    if (status != EX_NOERR) {
      sprintf(errmsg,
	      "Error: failed to get number of elements in side set  %"PRId64" in file id %d",
	      side_set_id, exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    if (tot_num_ss_elem == 0)
      continue;
   
    /* Allocate space for the side set element list */
    if (!(side_set_elem_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set element list for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Allocate space for the side set side list */
    if (!(side_set_side_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set side list for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    if (ex_get_set(exoid, EX_SIDE_SET, side_set_id, 
		   side_set_elem_list, side_set_side_list) != EX_NOERR) {
      sprintf(errmsg,
	      "Error: failed to get side set  %"PRId64" in file id %d",
	      side_set_id, exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Allocate space for the ss element index array */
    if (!(ss_elem_ndx=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set elem sort array for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Sort side set element list into index array  - non-destructive */
    if (int_size == sizeof(int64_t)) {
      /* Sort side set element list into index array  - non-destructive */
      int64_t *elems = (int64_t*)ss_elem_ndx;
      for (i=0;i<tot_num_ss_elem;i++) {
	elems[i] = i; /* init index array to current position */
      }
      ex_iqsort64(side_set_elem_list, ss_elem_ndx,tot_num_ss_elem);
    } else {
      /* Sort side set element list into index array  - non-destructive */
      int *elems = (int*)ss_elem_ndx;
      for (i=0;i<tot_num_ss_elem;i++) {
	elems[i] = i; /* init index array to current position */
      }
      ex_iqsort(side_set_elem_list, ss_elem_ndx,tot_num_ss_elem);
    }
    
    j = 0; /* The current element block... */
    for (ii=0;ii<tot_num_ss_elem;ii++) {
      
      int64_t elem_ndx;
      if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
	elem_ndx = ((int64_t*)ss_elem_ndx)[ii];
	elem     = ((int64_t*)side_set_elem_list)[elem_ndx];
	side     = ((int64_t*)side_set_side_list)[elem_ndx]-1;
      } else {
	elem_ndx = ((int*)ss_elem_ndx)[ii];
	elem     = ((int*)side_set_elem_list)[elem_ndx];
	side     = ((int*)side_set_side_list)[elem_ndx]-1;
      }
      
      /*
       * Since the elements are being accessed in sorted, order, the
       * block that contains the elements must progress sequentially
       * from block 0 to block[num_elem_blks-1]. Once we find an element
       * not in this block, find a following block that contains it...
       */       
      for ( ; j<num_elem_blks; j++) {
	if (elem <= elem_blk_parms[j].elem_ctr) {
	  break;
	}
      }
      
      if (j < num_elem_blks) {
	assert(side < elem_blk_parms[j].num_sides);  
	side_set_node_cnt_list[elem_ndx+ioff] = elem_blk_parms[j].num_nodes_per_side[side];
      } else {
	exerrval = EX_BADPARAM;
	sprintf(errmsg,
		"Error: Invalid element number  %"PRId64" found in side set  %"PRId64" in file %d",
		elem, side_set_id, exoid);
	ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
	goto error_ret;
      }
    }
    ss_elem_ndx        = ex_safe_free(ss_elem_ndx);
    side_set_elem_list = ex_safe_free(side_set_elem_list);
    side_set_side_list = ex_safe_free(side_set_side_list);
    ioff += tot_num_ss_elem;
  }
    
  /* All done: release allocated memory */
  ex_safe_free(elem_blk_ids);
  ex_safe_free(side_set_ids);
  ex_safe_free(ss_elem_ndx);
  ex_safe_free(side_set_elem_list);
  ex_safe_free(side_set_side_list);
  ex_safe_free(elem_blk_parms);
  return(EX_NOERR);

 error_ret:
  ex_safe_free(elem_blk_ids);
  ex_safe_free(side_set_ids);
  ex_safe_free(ss_elem_ndx);
  ex_safe_free(side_set_elem_list);
  ex_safe_free(side_set_side_list);
  ex_safe_free(elem_blk_parms);
  return (EX_FATAL);
}
