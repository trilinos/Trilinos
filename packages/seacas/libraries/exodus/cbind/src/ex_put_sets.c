/*
 * Copyright (c) 2012 Sandia Corporation. Under the terms of Contract
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

#include "exodusII.h"
#include "exodusII_int.h"
#include <stdlib.h> /* for free() */

/*!
 * writes the set parameters and optionally set data for 1 or more sets
 * \param   exoid                   exodus file id
 * \param   set_count               number of sets to write
 * \param  *sets                    array of ex_set structures
 */

int ex_put_sets (int   exoid,
		 size_t set_count,
		 const struct ex_set *sets)
{
  size_t i;
  int needs_define = 0;
  int set_stat;
  int dimid, varid, status, dims[1];
  int set_id_ndx;
  size_t start[1]; 
  int cur_num_sets;
  char errmsg[MAX_ERR_LENGTH];
  int* sets_to_define = NULL;
  char* numentryptr 	= NULL;
  char* entryptr = NULL;
  char* extraptr = NULL;
  char* idsptr = NULL;
  char* statptr = NULL;
  char* numdfptr = NULL;
  char* factptr = NULL;

  size_t int_size;
  
  exerrval = 0; /* clear error code */

  sets_to_define = malloc(set_count*sizeof(int));
  
  /* Note that this routine can be called:
     1) just define the sets
     2) just output the set data (after a previous call to define)
     3) define and output the set data in one call.
  */
  for (i=0; i < set_count; i++) {
    /* first check if any sets are specified */
    if ((status = nc_inq_dimid(exoid, ex_dim_num_objects(sets[i].type), &dimid)) != NC_NOERR) {
      if (status == NC_EBADDIM) {
	exerrval = status;
	sprintf(errmsg,
		"Error: no %ss defined for file id %d", ex_name_of_object(sets[i].type), exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
      } else {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to locate %ss defined in file id %d",
		ex_name_of_object(sets[i].type), exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
      }
      return (EX_FATAL);
    }

    set_id_ndx = ex_id_lkup(exoid, sets[i].type, sets[i].id);
    if (exerrval != EX_LOOKUPFAIL) {  /* found the side set id, so set is already defined... */
      sets_to_define[i] = 0;
      continue;
    } else {
      needs_define++;
      sets_to_define[i] = 1;
    }
  }
    
  if (needs_define > 0) {
    /* put netcdf file into define mode  */
    if ((status = nc_redef (exoid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to put file id %d into define mode",
	      exoid);
      ex_err("ex_put_sets",errmsg,exerrval);
      return (EX_FATAL);
    }
    
    for (i=0; i < set_count; i++) {
      if (sets_to_define[i] == 0)
	continue;
      
      /*   NOTE: ex_inc_file_item finds the current number of sets defined
	   for a specific file and returns that value incremented. */
      cur_num_sets=ex_inc_file_item(exoid, ex_get_counter_list(sets[i].type));
      set_id_ndx = cur_num_sets + 1;
      sets_to_define[i] = set_id_ndx;
      
      if (sets[i].num_entry == 0)
	continue;
      
      /* setup pointers based on set_type */
      if (sets[i].type == EX_NODE_SET) {
	numentryptr = DIM_NUM_NOD_NS(set_id_ndx);
	entryptr = VAR_NODE_NS(set_id_ndx);
	extraptr = NULL;
	/* note we are using DIM_NUM_NODE_NS instead of DIM_NUM_DF_NS */
	numdfptr = DIM_NUM_NOD_NS(set_id_ndx);
	factptr = VAR_FACT_NS(set_id_ndx);
      }
      else if (sets[i].type == EX_EDGE_SET) {
	numentryptr = DIM_NUM_EDGE_ES(set_id_ndx);
	entryptr = VAR_EDGE_ES(set_id_ndx);
	extraptr = VAR_ORNT_ES(set_id_ndx);
	numdfptr = DIM_NUM_DF_ES(set_id_ndx);
	factptr = VAR_FACT_ES(set_id_ndx);
      }
      else if (sets[i].type == EX_FACE_SET) {
	numentryptr = DIM_NUM_FACE_FS(set_id_ndx);
	entryptr = VAR_FACE_FS(set_id_ndx);
	extraptr = VAR_ORNT_FS(set_id_ndx);
	numdfptr = DIM_NUM_DF_FS(set_id_ndx);
	factptr = VAR_FACT_FS(set_id_ndx);
      }
      else if (sets[i].type == EX_SIDE_SET) {
	numentryptr = DIM_NUM_SIDE_SS(set_id_ndx);
	entryptr = VAR_ELEM_SS(set_id_ndx);
	extraptr = VAR_SIDE_SS(set_id_ndx);
	numdfptr = DIM_NUM_DF_SS(set_id_ndx);
	factptr = VAR_FACT_SS(set_id_ndx);
      }
      else if (sets[i].type == EX_ELEM_SET) {
	numentryptr = DIM_NUM_ELE_ELS(set_id_ndx);
	entryptr = VAR_ELEM_ELS(set_id_ndx);
	extraptr = NULL;
	numdfptr = DIM_NUM_DF_ELS(set_id_ndx);
	factptr = VAR_FACT_ELS(set_id_ndx);
      }

      /* define dimensions and variables */
      if ((status = nc_def_dim(exoid, numentryptr,
			       sets[i].num_entry, &dimid)) != NC_NOERR) {
	exerrval = status;
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: %s %"PRId64" -- size already defined in file id %d",
		  ex_name_of_object(sets[i].type), sets[i].id,exoid);
	  ex_err("ex_put_sets",errmsg,exerrval);
	}
	else {
	  sprintf(errmsg,
		  "Error: failed to define number of entries in %s %"PRId64" in file id %d",
		  ex_name_of_object(sets[i].type), sets[i].id,exoid);
	  ex_err("ex_put_sets",errmsg,exerrval);
	}
	goto error_ret;
      }
      
      int_size = sizeof(int);
      if (ex_int64_status(exoid) & EX_BULK_INT64_DB) {
	int_size = sizeof(int64_t);
      }
      
      /* create variable array in which to store the entry lists */
      dims[0] = dimid;
      if ((status = nc_def_var(exoid, entryptr, int_size, 1, dims, &varid)) != NC_NOERR) {
	exerrval = status;
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: entry list already exists for %s %"PRId64" in file id %d",
		  ex_name_of_object(sets[i].type), sets[i].id,exoid);
	  ex_err("ex_put_sets",errmsg,exerrval);
	} else {
	  sprintf(errmsg,
		  "Error: failed to create entry list for %s %"PRId64" in file id %d",
		  ex_name_of_object(sets[i].type), sets[i].id,exoid);
	  ex_err("ex_put_sets",errmsg,exerrval);
	}
	goto error_ret;            /* exit define mode and return */
      }
      ex_compress_variable(exoid, varid, 1);
      
      if (extraptr) {
	if ((status = nc_def_var(exoid, extraptr, int_size, 1, dims, &varid)) != NC_NOERR) {
	  exerrval = status;
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: extra list already exists for %s %"PRId64" in file id %d",
		    ex_name_of_object(sets[i].type), sets[i].id, exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create extra list for %s %"PRId64" in file id %d",
		    ex_name_of_object(sets[i].type), sets[i].id,exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	  }
	  goto error_ret;         /* exit define mode and return */
	}
	ex_compress_variable(exoid, varid, 1);
      }

      /* Create distribution factors variable if required */
      if (sets[i].num_distribution_factor > 0) {
	if (sets[i].type != EX_SIDE_SET) {
	  /* but sets[i].num_distribution_factor must equal number of nodes */
	  if (sets[i].num_distribution_factor != sets[i].num_entry) {
	    exerrval = EX_FATAL;
	    sprintf(errmsg,
		    "Error: # dist fact (%"PRId64") not equal to # nodes (%"PRId64") in node  set %"PRId64" file id %d",
		    sets[i].num_distribution_factor, sets[i].num_entry, sets[i].id, exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	    goto error_ret;    /* exit define mode and return */
	  }
	} else {
	  /* resuse dimid from entry lists */
	  if ((status = nc_def_dim(exoid, numdfptr, 
				   sets[i].num_distribution_factor, &dimid)) != NC_NOERR) {
	    exerrval = status;
	    sprintf(errmsg,
		    "Error: failed to define number of dist factors in %s %"PRId64" in file id %d",
		    ex_name_of_object(sets[i].type), sets[i].id,exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	    goto error_ret;          /* exit define mode and return */
	  }
	}
	
	/* create variable array in which to store the set distribution factors
	 */
	dims[0] = dimid;
	if ((status = nc_def_var(exoid, factptr, nc_flt_code(exoid), 1, dims, &varid)) != NC_NOERR) {
	  exerrval = status;
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: dist factors list already exists for %s %"PRId64" in file id %d",
		    ex_name_of_object(sets[i].type), sets[i].id,exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create dist factors list for %s %"PRId64" in file id %d",
		    ex_name_of_object(sets[i].type), sets[i].id,exoid);
	    ex_err("ex_put_sets",errmsg,exerrval);
	  }
	  goto error_ret;            /* exit define mode and return */
	}
	ex_compress_variable(exoid, varid, 2);
      }
    }

    /* leave define mode  */
    if ((status = nc_enddef (exoid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to complete definition in file id %d", exoid);
      ex_err("ex_put_sets",errmsg,exerrval);
      return (EX_FATAL);
    }

    /* Output the set ids and status... */
    for (i=0; i < set_count; i++) {
    /* setup pointers based on sets[i].type */
      if (sets[i].type == EX_NODE_SET) {
	idsptr = VAR_NS_IDS;
	statptr = VAR_NS_STAT;
      }
      else if (sets[i].type == EX_EDGE_SET) {
	idsptr = VAR_ES_IDS;
	statptr = VAR_ES_STAT;
      }
      else if (sets[i].type == EX_FACE_SET) {
	idsptr = VAR_FS_IDS;
	statptr = VAR_FS_STAT;
      }
      else if (sets[i].type == EX_SIDE_SET) {
	idsptr = VAR_SS_IDS;
	statptr = VAR_SS_STAT;
      }
      else if (sets[i].type == EX_ELEM_SET) {
	idsptr = VAR_ELS_IDS;
	statptr = VAR_ELS_STAT;
      }
      
      /* first: get id of set id variable */
      if ((status = nc_inq_varid(exoid, idsptr, &varid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to locate %s %"PRId64" in file id %d", ex_name_of_object(sets[i].type),
		sets[i].id, exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
	return (EX_FATAL);
      }
      
      /* write out set id */
      start[0] = sets_to_define[i]-1;
      status = nc_put_var1_longlong(exoid, varid, start, (long long*)&sets[i].id);
    
      if (status != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to store %s id %"PRId64" in file id %d", ex_name_of_object(sets[i].type),
		sets[i].id, exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
	return (EX_FATAL);
      }
      
      set_stat = (sets[i].num_entry == 0) ? 0 : 1;
      
      if ((status = nc_inq_varid(exoid, statptr, &varid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to locate %s status in file id %d", ex_name_of_object(sets[i].type),
		exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
	return (EX_FATAL);
      }
      
      if ((status = nc_put_var1_int(exoid, varid, start, &set_stat)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to store %s %"PRId64" status to file id %d", ex_name_of_object(sets[i].type),
		sets[i].id, exoid);
	ex_err("ex_put_sets",errmsg,exerrval);
	return (EX_FATAL);
      }
    }
    free(sets_to_define);
  }
  
  /* Sets are now all defined; see if any set data needs to be output... */
  status = EX_NOERR;
  for (i=0; i < set_count; i++) {
    int stat;
    if (sets[i].entry_list != NULL || sets[i].extra_list != NULL) {
      /* NOTE: ex_put_set will write the warning/error message... */
      stat = ex_put_set(exoid, sets[i].type, sets[i].id, sets[i].entry_list, sets[i].extra_list);
      if (stat != EX_NOERR) status = EX_FATAL;
    }
    if (sets[i].distribution_factor_list != NULL) {
      /* NOTE: ex_put_set_dist_fact will write the warning/error message... */
      stat = ex_put_set_dist_fact(exoid, sets[i].type, sets[i].id, sets[i].distribution_factor_list);
      if (stat != EX_NOERR) status = EX_FATAL;
    }
  }  
  return (status);

  /* Fatal error: exit definition mode and return */
 error_ret:
  free(sets_to_define);
  
  if (nc_enddef (exoid) != NC_NOERR) {    /* exit define mode */
    sprintf(errmsg,
	    "Error: failed to complete definition for file id %d",
	    exoid);
    ex_err("ex_put_sets",errmsg,exerrval);
  }
  return (EX_FATAL);
}
