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
* expblk - ex_put_block: write edge, face, or element block parameters
*
* entry conditions - 
*   input parameters:
*       int     idexo                   exodus file id
*       int     blk_type                type of block (edge, face, or element)
*       int     blk_id                  block identifer
*       char*   entry_descrip           string describing shape of entries in the block
*       int     num_entries_this_blk    number of entries(records) in the block
*       int     num_nodes_per_entry     number of nodes per block entry
*       int     block.num_edges_per_entry     number of edges per block entry
*       int     block.num_faces_per_entry     number of faces per block entry
*       int     num_attribute      number of attributes per block entry
*
* exit conditions - 
*
*
*****************************************************************************/

#include "exodusII.h"
#include "exodusII_int.h"
#include <string.h>

/*!
 * writes the parameters used to describe an element/face/edge block
 * \param   exoid                   exodus file id
 * \param   block                   ex_block structure describing block counts
 */

int ex_put_block_param( int         exoid,
			const ex_block block)
{
  int conn_int_type;
  int status;
  int arbitrary_polyhedra = 0; /* 1 if block is arbitrary 2d polyhedra type; 2 if 3d polyhedra */
  int att_name_varid = -1;
  int varid, dimid, dims[2], blk_id_ndx, blk_stat, strdim;
  size_t start[2];
  int num_blk;
  size_t temp;
  int cur_num_blk, numblkdim, numattrdim;
  int nnodperentdim = -1;
  int nedgperentdim = -1;
  int nfacperentdim = -1;
  int connid;
  int npeid;
  char errmsg[MAX_ERR_LENGTH];
  char entity_type1[5];
  char entity_type2[5];
  const char* dnumblk = NULL;
  const char* vblkids = NULL;
  const char* vblksta = NULL;
  const char* vnodcon = NULL;
  const char* vnpecnt = NULL;
  const char* vedgcon = NULL;
  const char* vfaccon = NULL;
  const char* vconn   = NULL;
  const char* vattnam = NULL;
  const char* vblkatt = NULL;
  const char* dneblk  = NULL;
  const char* dnape   = NULL;
  const char* dnnpe   = NULL;
  const char* dnepe   = NULL;
  const char* dnfpe   = NULL;

  exerrval  = 0; /* clear error code */

  switch (block.type) {
  case EX_EDGE_BLOCK:
    dnumblk = DIM_NUM_ED_BLK;
    vblkids = VAR_ID_ED_BLK;
    vblksta = VAR_STAT_ED_BLK;
    break;
  case EX_FACE_BLOCK:
    dnumblk = DIM_NUM_FA_BLK;
    vblkids = VAR_ID_FA_BLK;
    vblksta = VAR_STAT_FA_BLK;
    break;
  case EX_ELEM_BLOCK:
    dnumblk = DIM_NUM_EL_BLK;
    vblkids = VAR_ID_EL_BLK;
    vblksta = VAR_STAT_EL_BLK;
    break;
  default:
    exerrval = EX_BADPARAM;
    sprintf( errmsg, "Error: Bad block type (%d) specified for file id %d",
	     block.type, exoid );
    ex_err( "ex_put_block", errmsg, exerrval );
    return (EX_FATAL);
  }

  /* first check if any element blocks are specified */

  if ((status = ex_get_dimension(exoid, dnumblk, ex_name_of_object(block.type),
				 &temp, &dimid, "ex_put_block")) != NC_NOERR) {
    return EX_FATAL;
  }
  num_blk = temp;
  
  /* Next: Make sure that this is not a duplicate element block id by
     searching the vblkids array.
     WARNING: This must be done outside of define mode because id_lkup accesses
     the database to determine the position
  */

  if ((status = nc_inq_varid(exoid, vblkids, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate %s ids in file id %d",
	    ex_name_of_object(block.type), exoid);
    ex_err("ex_put_block",errmsg,exerrval);
  }

  ex_id_lkup(exoid,block.type,block.id); /* Error value used, but don't need return value */
  if (exerrval != EX_LOOKUPFAIL) {   /* found the element block id */
    exerrval = EX_FATAL;
    sprintf(errmsg,
            "Error: %s id %"PRId64" already exists in file id %d",
	    ex_name_of_object(block.type), block.id,exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Keep track of the total number of element blocks defined using a counter 
     stored in a linked list keyed by exoid.
     NOTE: ex_get_file_item  is a function that finds the number of element 
     blocks for a specific file and returns that value incremented.
  */
  cur_num_blk=ex_get_file_item(exoid, ex_get_counter_list(block.type));
  if (cur_num_blk >= num_blk) {
    exerrval = EX_FATAL;
    sprintf(errmsg,
	    "Error: exceeded number of %ss (%d) defined in file id %d",
	    ex_name_of_object(block.type), num_blk,exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }


  /*   NOTE: ex_get_file_item  is a function that finds the number of element
       blocks for a specific file and returns that value incremented. */
  cur_num_blk=ex_inc_file_item(exoid, ex_get_counter_list(block.type));
  start[0] = cur_num_blk;

  /* write out block id to previously defined id array variable*/
  status = nc_put_var1_longlong(exoid, varid, start, (long long*)&block.id);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to store %s id to file id %d",
	    ex_name_of_object(block.type), exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  blk_id_ndx = start[0]+1; /* element id index into vblkids array*/

  if (block.num_entry == 0) /* Is this a NULL element block? */
    blk_stat = 0; /* change element block status to NULL */
  else
    blk_stat = 1; /* change element block status to TRUE */

  if ((status = nc_inq_varid (exoid, vblksta, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate %s status in file id %d",
	    ex_name_of_object(block.type), exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_put_var1_int(exoid, varid, start, &blk_stat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to store %s id %"PRId64" status to file id %d",
	    ex_name_of_object(block.type), block.id, exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (block.num_entry == 0) {/* Is this a NULL element block? */
    return(EX_NOERR);
  }


  /* put netcdf file into define mode  */
  if ((status=nc_redef (exoid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,"Error: failed to place file id %d into define mode",exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  switch (block.type) {
  case EX_EDGE_BLOCK:
    dneblk = DIM_NUM_ED_IN_EBLK(blk_id_ndx);
    dnnpe = DIM_NUM_NOD_PER_ED(blk_id_ndx);
    dnepe = 0;
    dnfpe = 0;
    dnape = DIM_NUM_ATT_IN_EBLK(blk_id_ndx);
    vblkatt = VAR_EATTRIB(blk_id_ndx);
    vattnam = VAR_NAME_EATTRIB(blk_id_ndx);
    vnodcon = VAR_EBCONN(blk_id_ndx);
    vedgcon = 0;
    vfaccon = 0;
    break;
  case EX_FACE_BLOCK:
    dneblk = DIM_NUM_FA_IN_FBLK(blk_id_ndx);
    dnnpe = DIM_NUM_NOD_PER_FA(blk_id_ndx);
    dnepe = 0;
    dnfpe = 0;
    dnape = DIM_NUM_ATT_IN_FBLK(blk_id_ndx);
    vblkatt = VAR_FATTRIB(blk_id_ndx);
    vattnam = VAR_NAME_FATTRIB(blk_id_ndx);
    vnodcon = VAR_FBCONN(blk_id_ndx);
    vnpecnt = VAR_FBEPEC(blk_id_ndx);
    vedgcon = 0;
    vfaccon = 0;
    break;
  case EX_ELEM_BLOCK:
    dneblk = DIM_NUM_EL_IN_BLK(blk_id_ndx);
    dnnpe = DIM_NUM_NOD_PER_EL(blk_id_ndx);
    dnepe = DIM_NUM_EDG_PER_EL(blk_id_ndx);
    dnfpe = DIM_NUM_FAC_PER_EL(blk_id_ndx);
    dnape = DIM_NUM_ATT_IN_BLK(blk_id_ndx);
    vblkatt = VAR_ATTRIB(blk_id_ndx);
    vattnam = VAR_NAME_ATTRIB(blk_id_ndx);
    vnodcon = VAR_CONN(blk_id_ndx);
    vnpecnt = VAR_EBEPEC(blk_id_ndx);
    vedgcon = VAR_ECONN(blk_id_ndx);
    vfaccon = VAR_FCONN(blk_id_ndx);
    break;
  default:
    exerrval = 1005;
    sprintf(errmsg,
	    "Internal Error: unrecognized block type in switch: %d in file id %d",
	    block.type,exoid);
    ex_err("ex_put_block",errmsg,EX_MSG);
    return (EX_FATAL);              /* number of attributes not defined */
  }
  /* define some dimensions and variables*/

  if ((status = nc_def_dim(exoid,dneblk,block.num_entry, &numblkdim )) != NC_NOERR) {
    if (status == NC_ENAMEINUSE) {        /* duplicate entry */
      exerrval = status;
      sprintf(errmsg,
	      "Error: %s %"PRId64" already defined in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
    } else {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define number of entities/block for %s %"PRId64" file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
    }
    goto error_ret;         /* exit define mode and return */
  }

  if ( dnnpe && block.num_nodes_per_entry > 0) {
    /* A nfaced block would not have any nodes defined... */
    if ((status = nc_def_dim(exoid,dnnpe,block.num_nodes_per_entry, &nnodperentdim)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define number of nodes/entity for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
  }

  if (dnepe && block.num_edges_per_entry > 0 ) {
    if ((status = nc_def_dim (exoid,dnepe,block.num_edges_per_entry, &nedgperentdim)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define number of edges/entity for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
  }

  if ( dnfpe && block.num_faces_per_entry > 0 ) {
    if ((status = nc_def_dim(exoid,dnfpe,block.num_faces_per_entry, &nfacperentdim)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define number of faces/entity for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
  }

  /* element attribute array */
  if (block.num_attribute > 0) {

    if ((status = nc_def_dim(exoid, dnape, block.num_attribute, &numattrdim)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define number of attributes in %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }

    dims[0] = numblkdim;
    dims[1] = numattrdim;

    if ((status = nc_def_var(exoid, vblkatt, nc_flt_code(exoid), 2, dims, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error:  failed to define attributes for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
    ex_compress_variable(exoid, varid, 2);

    /* inquire previously defined dimensions  */
    if ((status = nc_inq_dimid(exoid, DIM_STR_NAME, &strdim)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get string length in file id %d",exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      return (EX_FATAL);
    }
     
    /* Attribute names... */
    dims[0] = numattrdim;
    dims[1] = strdim;
	    
    if ((status = nc_def_var(exoid, vattnam, NC_CHAR, 2, dims, &att_name_varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to define %s attribute name array in file id %d",
	      ex_name_of_object(block.type), exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
  }

  conn_int_type = NC_INT;
  if (ex_int64_status(exoid) & EX_BULK_INT64_DB) {
    conn_int_type = NC_INT64;
  }

  /* See if storing an 'nsided' element block (arbitrary 2d polyhedra or super element) */
  if (strlen(block.topology) >= 3) {
    if ((block.topology[0] == 'n' || block.topology[0] == 'N') &&
	(block.topology[1] == 's' || block.topology[1] == 'S') &&
	(block.topology[2] == 'i' || block.topology[2] == 'I'))
      arbitrary_polyhedra = 1;
    else if ((block.topology[0] == 'n' || block.topology[0] == 'N') &&
	     (block.topology[1] == 'f' || block.topology[1] == 'F') &&
	     (block.topology[2] == 'a' || block.topology[2] == 'A'))
      /* If a FACE_BLOCK, then we are dealing with the faces of the nfaced block. */
      arbitrary_polyhedra = block.type == EX_FACE_BLOCK ? 1 : 2;
  }

  /* element connectivity array */
  if (arbitrary_polyhedra > 0) {
    if (block.type != EX_FACE_BLOCK && block.type != EX_ELEM_BLOCK) {
      exerrval = EX_BADPARAM;
      sprintf( errmsg, "Error: Bad block type (%d) for nsided/nfaced block in file id %d",
	       block.type, exoid );
      ex_err( "ex_put_block", errmsg, exerrval );
      return (EX_FATAL);
    }

    if (arbitrary_polyhedra == 1) {
      dims[0] = nnodperentdim;
      vconn = vnodcon;

      /* store entity types as attribute of npeid variable -- node/elem, node/face, face/elem*/
      strcpy(entity_type1, "NODE");
      if (block.type == EX_ELEM_BLOCK)
	strcpy(entity_type2, "ELEM");
      else
	strcpy(entity_type2, "FACE");
    } else if (arbitrary_polyhedra == 2) {
      dims[0] = nfacperentdim;
      vconn = vfaccon;

      /* store entity types as attribute of npeid variable -- node/elem, node/face, face/elem*/
      strcpy(entity_type1, "FACE");
      strcpy(entity_type2, "ELEM");
    }

    if ((status = nc_def_var(exoid, vconn, conn_int_type, 1, dims, &connid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to create connectivity array for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
    
    /* element face-per-element or node-per-element count array */
    dims[0] = numblkdim;
    
    if ((status = nc_def_var(exoid, vnpecnt, conn_int_type, 1, dims, &npeid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to create face- or node- per-entity count array for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id, exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }

    if ((status = nc_put_att_text(exoid, npeid, "entity_type1", strlen(entity_type1)+1,
				  entity_type1)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to store entity type attribute text for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id, exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
    if ((status = nc_put_att_text(exoid, npeid, "entity_type2", strlen(entity_type2)+1,
				  entity_type2)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to store entity type attribute text for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id, exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
  } else {
    /* "Normal" (non-polyhedra) element block type */
    dims[0] = numblkdim;
    dims[1] = nnodperentdim;
    
    if ((status = nc_def_var(exoid, vnodcon, conn_int_type, 2, dims, &connid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to create connectivity array for %s %"PRId64" in file id %d",
	      ex_name_of_object(block.type), block.id,exoid);
      ex_err("ex_put_block",errmsg,exerrval);
      goto error_ret;         /* exit define mode and return */
    }
    ex_compress_variable(exoid, connid, 1);
  }
  /* store element type as attribute of connectivity variable */
  if ((status = nc_put_att_text(exoid, connid, ATT_NAME_ELB, strlen(block.topology)+1, 
				block.topology)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to store %s type name %s in file id %d",
	    ex_name_of_object(block.type), block.topology,exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    goto error_ret;         /* exit define mode and return */
  }

  if (arbitrary_polyhedra == 0) {
    if (vedgcon && block.num_edges_per_entry ) {
      dims[0] = numblkdim;
      dims[1] = nedgperentdim;
      
      if ((status = nc_def_var(exoid, vedgcon, conn_int_type, 2, dims, &varid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to create edge connectivity array for %s %"PRId64" in file id %d",
		ex_name_of_object(block.type), block.id,exoid);
	ex_err("ex_put_block",errmsg,exerrval);
	goto error_ret;         /* exit define mode and return */
      }
    }
    
    if ( vfaccon && block.num_faces_per_entry ) {
      dims[0] = numblkdim;
      dims[1] = nfacperentdim;
      
      if ((status = nc_def_var(exoid, vfaccon, conn_int_type, 2, dims, &varid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to create face connectivity array for %s %"PRId64" in file id %d",
		ex_name_of_object(block.type), block.id,exoid);
	ex_err("ex_put_block",errmsg,exerrval);
	goto error_ret;         /* exit define mode and return */
      }
    }
  }
  /* leave define mode  */

  if ((exerrval=nc_enddef (exoid)) != NC_NOERR) {
    sprintf(errmsg,
	    "Error: failed to complete %s definition in file id %d", 
	    ex_name_of_object(block.type), exoid);
    ex_err("ex_put_block",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Output a dummy empty attribute name in case client code doesn't
     write anything; avoids corruption in some cases.
  */
  if (block.num_attribute > 0 && att_name_varid >= 0) {
    size_t  count[2];
    char *text = "";
    size_t i;

    count[0] = 1;
    start[1] = 0;
    count[1] = strlen(text)+1;
  
    for (i = 0; i < block.num_attribute; i++) {
      start[0] = i;
      nc_put_vara_text(exoid, att_name_varid, start, count, text);
    }
  }

  return (EX_NOERR);

  /* Fatal error: exit definition mode and return */
 error_ret:
  if (nc_enddef (exoid) != NC_NOERR) {    /* exit define mode */
    sprintf(errmsg,
	    "Error: failed to complete definition for file id %d",
	    exoid);
    ex_err("ex_put_block",errmsg,exerrval);
  }
  return (EX_FATAL);
}

