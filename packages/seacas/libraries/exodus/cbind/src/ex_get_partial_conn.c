/*
 * Copyright (c) 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
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
* exgnconn - exodusII read partial edge/face/element block connectivity
*
* entry conditions - 
*   expelb must be called first to establish element block parameters.
*   input parameters:
*       int     exoid           exodus file id
*       int     blk_type        block type (edge, face, element)
*       int     blk_id          block id
*       int     start_num       The starting index (1-based) of entities to read
*       int     num_ent         The number of entities to read connectivity info for.
* 
* exit conditions - 
*       int*    nodeconn        nodal connectivity array
*       int*    edgeconn        edge connectivity array (where applicable)
*       int*    faceconn        face connectivity array (where applicable)
*
* revision history - 
*
*/
#include <stdlib.h>
#include "exodusII.h"
#include "exodusII_int.h"

/*
 * reads the connectivity array for an element block
 */

int ex_get_partial_conn( int   exoid,
		   ex_entity_type blk_type,
		   ex_entity_id   blk_id,
		   int64_t   start_num,
		   int64_t   num_ent,
		   void_int*  nodeconn,
		   void_int*  edgeconn,
		   void_int*  faceconn )
{
  int connid = -1;
  int econnid = -1;
  int fconnid = -1;

  int blk_id_ndx, status;

  int numnodperentdim = -1;
  int numedgperentdim = -1;
  int numfacperentdim = -1;

  int iexit = (EX_NOERR); /* exit status */

  size_t num_nodes_per_entry = 0;
  size_t num_edges_per_entry = 0;
  size_t num_faces_per_entry = 0;
  char errmsg[MAX_ERR_LENGTH];

  const char* dnumnodent = NULL;
  const char* dnumedgent = NULL;
  const char* dnumfacent = NULL;
  const char* vnodeconn = NULL;
  const char* vedgeconn = NULL;
  const char* vfaceconn = NULL;

  /* The partial connectivity input function can currently only handle nodal
   * connectivity.  Print a warning if edgeconn or faceconn is non-NULL
   */
  if (edgeconn != NULL || faceconn != NULL) {
    exerrval = 1005;
    sprintf(errmsg,
	    "Warning: ex_get_partial_conn only supports nodal connectivity at this time. %s %"PRId64" in file id %d",
	    ex_name_of_object(blk_type),blk_id,exoid);
    ex_err("ex_get_partial_conn",errmsg,EX_MSG);
  }

  exerrval = 0; /* clear error code */

  /* Locate index of element block id in VAR_ID_EL_BLK array */

  blk_id_ndx = ex_id_lkup(exoid,blk_type,blk_id);
  if (exerrval != 0) 
    {
      if (exerrval == EX_NULLENTITY)
	{
	  sprintf(errmsg,
		  "Warning: no connectivity array for NULL %s %"PRId64" in file id %d",
		  ex_name_of_object(blk_type),blk_id,exoid);
	  ex_err("ex_get_partial_conn",errmsg,EX_MSG);
	  return (EX_WARN); /* no connectivity array for this element block */
	}
      else
	{
	  sprintf(errmsg,
		  "Error: failed to locate %s id %"PRId64" in id array in file id %d",
		  ex_name_of_object(blk_type),blk_id,exoid);
	  ex_err("ex_get_partial_conn",errmsg,exerrval);
	  return (EX_FATAL);
	}
    }

  switch (blk_type) {
  case EX_EDGE_BLOCK:
    dnumnodent = DIM_NUM_NOD_PER_ED(blk_id_ndx);
    dnumedgent = 0;
    dnumfacent = 0;
    vnodeconn = VAR_EBCONN(blk_id_ndx);
    vedgeconn = 0;
    vfaceconn = 0;
    break;
  case EX_FACE_BLOCK:
    dnumnodent = DIM_NUM_NOD_PER_FA(blk_id_ndx);
    dnumedgent = 0;
    dnumfacent = 0;
    vnodeconn = VAR_FBCONN(blk_id_ndx);
    vedgeconn = 0;
    vfaceconn = 0;
    break;
  case EX_ELEM_BLOCK:
    dnumnodent = DIM_NUM_NOD_PER_EL(blk_id_ndx);
    dnumedgent = DIM_NUM_EDG_PER_EL(blk_id_ndx);
    dnumfacent = DIM_NUM_FAC_PER_EL(blk_id_ndx);
    vnodeconn = VAR_CONN(blk_id_ndx);
    vedgeconn = VAR_ECONN(blk_id_ndx);
    vfaceconn = VAR_FCONN(blk_id_ndx);
    break;
  default:
    exerrval = 1005;
    sprintf(errmsg,
	    "Internal Error: unrecognized block type in switch: %d in file id %d",
	    blk_type,exoid);
    ex_err("ex_get_partial_conn",errmsg,EX_MSG);
    return (EX_FATAL);              /* number of attributes not defined */
  }
  /* inquire id's of previously defined dimensions  */

  if ((status = nc_inq_dimid (exoid, dnumnodent, &numnodperentdim)) != NC_NOERR)
    {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate number of nodes/entity for %s %"PRId64" in file id %d",
	      ex_name_of_object(blk_type),blk_id,exoid);
      ex_err("ex_get_partial_conn",errmsg,exerrval);
      return(EX_FATAL);
    }

  if (nc_inq_dimlen(exoid, numnodperentdim, &num_nodes_per_entry) != NC_NOERR)
    {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get number of nodes/entity for %s %"PRId64" in file id %d",
	      ex_name_of_object(blk_type),blk_id,exoid);
      ex_err("ex_get_partial_conn",errmsg, exerrval);
      return(EX_FATAL);
    }

  if ( edgeconn && dnumedgent ) {
    num_edges_per_entry = 0;
    if ((status = nc_inq_dimid(exoid, dnumedgent, &numedgperentdim)) != NC_NOERR) {
      numedgperentdim = -1;
    } else {
      if ((status = nc_inq_dimlen(exoid, numedgperentdim, &num_edges_per_entry)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get number of edges/entry for %s %"PRId64" in file id %d",
		ex_name_of_object(blk_type),blk_id,exoid);
	ex_err("ex_get_partial_conn",errmsg, exerrval);
	return(EX_FATAL);
      }
    }
  }

  if ( faceconn && dnumfacent ) {
    num_faces_per_entry = 0;
    if ((status = nc_inq_dimid(exoid, dnumfacent, &numfacperentdim)) != NC_NOERR) {
      numfacperentdim = -1;
    } else {
      if ((status = nc_inq_dimlen(exoid, numfacperentdim, &num_faces_per_entry)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get number of faces/entry for %s %"PRId64" in file id %d",
		ex_name_of_object(blk_type),blk_id,exoid);
	ex_err("ex_get_partial_conn",errmsg, exerrval);
	return(EX_FATAL);
      }
    }
  }


  if ((status = nc_inq_varid(exoid, vnodeconn, &connid)) != NC_NOERR)
    {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate connectivity array for %s %"PRId64" in file id %d",
	      ex_name_of_object(blk_type),blk_id,exoid);
      ex_err("ex_get_partial_conn",errmsg, exerrval);
      return(EX_FATAL);
    }

  status = 0;
  if (edgeconn && (numedgperentdim > 0) &&
      ((status = nc_inq_varid (exoid, vedgeconn, &econnid)) != NC_NOERR))
    {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate edge connectivity array for %s %"PRId64" in file id %d",
	      ex_name_of_object(blk_type),blk_id,exoid);
      ex_err("ex_get_partial_conn",errmsg, exerrval);
      return(EX_FATAL);
    }

  if (faceconn && (numfacperentdim > 0) &&
      ((status = nc_inq_varid (exoid, vfaceconn, &fconnid)) != NC_NOERR))
    {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate face connectivity array for %s %"PRId64" in file id %d",
	      ex_name_of_object(blk_type),blk_id,exoid);
      ex_err("ex_get_partial_conn",errmsg, exerrval);
      return(EX_FATAL);
    }


  /* read in the connectivity array */
  if ( edgeconn && num_edges_per_entry > 0) {
    size_t start[2], count[2];
    
    start[0] = (start_num-1);
    start[1] = 0;
    count[0] = num_ent;
    count[1] = num_edges_per_entry;

    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      status = nc_get_vara_longlong(exoid, econnid, start, count, edgeconn);
    } else {
      status = nc_get_vara_int(exoid, econnid, start, count, edgeconn);
    }

    if (status != NC_NOERR)
      {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get edge connectivity array for %s %"PRId64" in file id %d",
		ex_name_of_object(blk_type),blk_id,exoid);
	ex_err("ex_get_partial_conn",errmsg, exerrval);
	return(EX_FATAL);
      }
  }

  if ( faceconn && num_faces_per_entry > 0) {
    size_t start[2], count[2];
    
    start[0] = start_num-1;
    start[1] = 0;
    count[0] = num_ent;
    count[1] = num_faces_per_entry;
     
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      status = nc_get_vara_longlong(exoid, fconnid, start, count, faceconn);
    } else {
      status = nc_get_vara_int(exoid, fconnid, start, count, faceconn);
    }
     
    if (status != NC_NOERR)
      {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get face connectivity array for %s %"PRId64" in file id %d",
		ex_name_of_object(blk_type),blk_id,exoid);
	ex_err("ex_get_partial_conn",errmsg, exerrval);
	return(EX_FATAL);
      }
  }

  if ( nodeconn && num_nodes_per_entry > 0) {
    size_t start[2], count[2];
     
    start[0] = start_num-1;
    start[1] = 0;
    count[0] = num_ent;
    count[1] = num_nodes_per_entry;
     
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      status = nc_get_vara_longlong(exoid, connid, start, count, nodeconn);
    } else {
      status = nc_get_vara_int(exoid, connid, start, count, nodeconn);
    }

    if (status != NC_NOERR)
      {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get connectivity array for %s %"PRId64" in file id %d",
		ex_name_of_object(blk_type),blk_id,exoid);
	ex_err("ex_get_partial_conn",errmsg, exerrval);
	return(EX_FATAL);
      }
  }

  return iexit;
}
