/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 *
 *     ne_leavedef()
 *     ne_catstr2()
 *     ne_id_lkup()
 *     ne_get_file_type()
 *     ne_put_version()
 *     ne_check_file_version()
 *     ne_get_idx()
 *
 *****************************************************************************
 * Much of this code is a modified version of what is found in NemesisI.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <string.h>
#include <math.h>

#include <netcdf.h>
#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

/* Global variables */
char *ne_ret_string;

int ne_leavedef(int neid,
                char *call_rout
                )
{
  char errmsg[MAX_ERR_LENGTH];
  int status;
  
  if ((status = nc_enddef(neid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to end define mode for file id %d",
            neid);
    ex_err(call_rout, errmsg, exerrval);

    return (EX_FATAL);
  }
  return (EX_NOERR);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
char *ne_catstr2(char *name,
                 int   num1,
                 int   num2
                 )
{
  char *func_name="ne_catstr2";

  char  errmsg[MAX_ERR_LENGTH];

  exerrval = 0;  /* clear error code */

  if (ne_ret_string == NULL) {
    ne_ret_string = (char *)malloc((NC_MAX_NAME+1)*sizeof(char));
    if (ne_ret_string == NULL) {
      exerrval = EX_MSG;
      sprintf(errmsg, "Error: Insufficient memory!\n");
      ex_err(func_name, errmsg, exerrval);
      return NULL;
    }
  }

  if (strlen(name) > NC_MAX_NAME) {
    exerrval = EX_MSG;
    sprintf(errmsg, "Error: name too long!");
    ex_err(func_name, errmsg, exerrval);

    return (NULL);
  }

  sprintf(ne_ret_string, "%s%d-%d", name, num1, num2);

  return ne_ret_string;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Note: This function assumes a 1-d vector of data for "ne_var_name".
 */
/*****************************************************************************/
int ne_id_lkup(int neid, char *ne_var_name, int64_t *idx, int ne_var_id)
{
  char    *func_name="ne_id_lkup";

  int      status;
  int      varid, ndims, dimid[1], ret=-1;
  nc_type  var_type;
  size_t   index, length, start[1];
  size_t    begin, end;
  int   id_val;

  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  if ((status = nc_inq_varid(neid, ne_var_name, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            ne_var_name, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* check if I need the length for this varible */
  if (idx[1] == -1) {
    /* Get the dimension IDs for this variable */
    if ((status = nc_inq_var(neid, varid, (char *) 0, &var_type, &ndims,
			     dimid, (int *) 0)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
          "Error: failed to find dimension ID for variable \"%s\" in file ID %d",
              ne_var_name, neid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    /* Get the length of this variable */
    if ((status = nc_inq_dimlen(neid, dimid[0], &length)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
          "Error: failed to find dimension for variable \"%s\" in file ID %d",
              ne_var_name, neid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    idx[1] = length;
  } /* End "if (idx[1] == -1)" */

  begin = idx[0];
  end = idx[1];

  /* Find the index by looping over each entry */
  for(index=begin; index < end; index++) {
    start[0] = index;
    if ((status = nc_get_var1_int(neid, varid, start, &id_val)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable \"%s\" in file ID %d",
              ne_var_name, neid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    if (id_val == ne_var_id) {
      ret = (int) index;
      break;
    }
  }

  return ret;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function retrieves the file type from a Nemesis file.
 */
/*****************************************************************************/
int ne_get_file_type(int neid,
                     char *ftype
                     )
{
  char   *func_name="ne_get_file_type";

  int  status;
  int     varid;
  int  lftype;

  char    errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  if ((status = nc_inq_varid(neid, VAR_FILE_TYPE, &varid)) != NC_NOERR) {

    /* If no file type is found, assume parallel */
    ftype[0] = 'p';
    ftype[1] = '\0';

    return (EX_NOERR);
  }

  if ((status = nc_get_var1_int(neid, varid, NULL, &lftype)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get variable \"%s\" from file ID %d",
            VAR_FILE_TYPE, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Set the appropriate character */
  if (lftype == 0)       strcpy(ftype, "p");
  else if (lftype == 1)  strcpy(ftype, "s");

  return (EX_NOERR);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs the Nemesis version information to the file.
 */
/*****************************************************************************/
int ne_put_version(int neid)
{
  char  *func_name="ne_put_version";
  int    status;
  float  file_ver, api_ver;

  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear the error code */

  file_ver = NEMESIS_FILE_VERSION;
  api_ver  = NEMESIS_API_VERSION;

  /* Check to see if the nemesis file version is already in the file */
  if (nc_get_att_float(neid, NC_GLOBAL, "nemesis_file_version", &file_ver) != NC_NOERR) {

    /* Output the Nemesis file version */
    if ((status = nc_put_att_float(neid, NC_GLOBAL, "nemesis_file_version", NC_FLOAT,
				   1, &file_ver)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output nemesis file version in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the Nemesis API version */
    if ((status = nc_put_att_float(neid, NC_GLOBAL, "nemesis_api_version", NC_FLOAT,
				   1, &api_ver)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output nemesis api version in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  return (EX_NOERR);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function checks that the version info is correct.
 */
/*****************************************************************************/
int ne_check_file_version(int neid)
{
  char  *func_name="ne_check_file_version";

  float  file_ver;

  int    status;
  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0;  /* clear error code */

  /* Get the file version */
  if ((status = nc_get_att_float(neid, NC_GLOBAL, "nemesis_file_version", &file_ver)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get the nemesis file version from file ID %d",
            neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (fabs(NEMESIS_FILE_VERSION-file_ver) > 0.001) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: Nemesis version mismatch in file ID %d!\n", neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function gets the index for the given variable at the
 * position given.
 */
/*****************************************************************************/
int ne_get_idx(int neid, char *ne_var_name, int64_t *index, int pos)
{
  char  *func_name="ne_get_idx";

  int      status;
  int      varid;
  size_t   start[1], count[1];
  int64_t  varidx[2];

  char   errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* set default values for idx */
  index[0] = 0;
  index[1] = -1;

  /*
   * assume that if there is an error returned, that this
   * means that this is a parallel file, and the index does
   * not exists. This is not an error
   */
  if ((status = nc_inq_varid(neid, ne_var_name, &varid)) == NC_NOERR) {
    /* check if we are at the beginning of the index vector */
    if (pos == 0) {
      start[0] = pos;
      count[0] = 1;
    } else {
      start[0] = pos - 1;
      count[0] = 2;
    }

    if ((status = nc_get_vara_longlong(neid, varid, start, count, varidx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable \"%s\" in file ID %d",
              ne_var_name, neid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    if (pos == 0) {
      index[0] = 0;
      index[1] = varidx[0];
    } else {
      index[0] = varidx[0];
      index[1] = varidx[1];
    }
  }

  return 1;
}

int ne_get_proc_count(int neid)
{
  int status;
  int dimid;
  size_t ltemp;
  int proc_count = -1;

  if ((status = nc_inq_dimid(neid, DIM_NUM_PROCS, &dimid)) != NC_NOERR) {
    return (-1);
  }

  /* Get the value of the number of processors */
  if ((status = nc_inq_dimlen(neid, dimid, &ltemp)) != NC_NOERR) {
    return (-1);
  }
  proc_count = ltemp;
  return proc_count;
}

struct ne_stat_struct {
  int neid;
  int *stats;
  size_t num;
} ne_stat_struct;

static struct ne_stat_struct *nem_int_e_maps = 0;
static struct ne_stat_struct *nem_bor_e_maps = 0;
static struct ne_stat_struct *nem_int_n_maps = 0;
static struct ne_stat_struct *nem_bor_n_maps = 0;
static struct ne_stat_struct *nem_ext_n_maps = 0;

int ne_fill_map_status(int neid,
			char *stat_var,
			int *stat_array)
{
  int status;
  int varid;
  char   errmsg[MAX_ERR_LENGTH];
  char *func_name = "ne_fill_map_status(internal function)";
  
  /* Now query the map status and fill the array. */
  if ((status = nc_inq_varid(neid, stat_var, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" from file ID %d",
	    VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_get_var_int(neid, varid, stat_array)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file %d",
	    VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}

int ne_get_map_status(int neid,
		      char *stat_var,
		      int proc_id,
		      int *stat)
{
  int  num_procs = 0;

  struct ne_stat_struct *maps = 0;
  if (strcmp(stat_var, VAR_INT_N_STAT ) == 0) {
    if (nem_int_n_maps == 0) {
      nem_int_n_maps = malloc(sizeof(ne_stat_struct));
      num_procs = ne_get_proc_count(neid);
      nem_int_n_maps->stats = malloc(sizeof(int) * num_procs);

      ne_fill_map_status(neid, stat_var, nem_int_n_maps->stats);
    }
    maps = nem_int_n_maps;
    maps->neid = neid;
    maps->num  = num_procs;
  } 
  else if (strcmp(stat_var, VAR_BOR_N_STAT ) == 0) {
    if (nem_bor_n_maps == 0) {
      nem_bor_n_maps = malloc(sizeof(ne_stat_struct));
      num_procs = ne_get_proc_count(neid);
      nem_bor_n_maps->stats = malloc(sizeof(int) * num_procs);

      ne_fill_map_status(neid, stat_var, nem_bor_n_maps->stats);
    }
    maps = nem_bor_n_maps;
    maps->neid = neid;
    maps->num  = num_procs;
  } 
  else if (strcmp(stat_var, VAR_EXT_N_STAT ) == 0) {
    if (nem_ext_n_maps == 0) {
      nem_ext_n_maps = malloc(sizeof(ne_stat_struct));
      num_procs = ne_get_proc_count(neid);
      nem_ext_n_maps->stats = malloc(sizeof(int) * num_procs);

      ne_fill_map_status(neid, stat_var, nem_ext_n_maps->stats);
    }
    maps = nem_ext_n_maps;
    maps->neid = neid;
    maps->num  = num_procs;
  } 
  else if (strcmp(stat_var, VAR_INT_E_STAT ) == 0) {
    if (nem_int_e_maps == 0) {
      nem_int_e_maps = malloc(sizeof(ne_stat_struct));
      num_procs = ne_get_proc_count(neid);
      nem_int_e_maps->stats = malloc(sizeof(int) * num_procs);

      ne_fill_map_status(neid, stat_var, nem_int_e_maps->stats);
    }
    maps = nem_int_e_maps;
    maps->neid = neid;
    maps->num  = num_procs;
  } 
  else if (strcmp(stat_var, VAR_BOR_E_STAT ) == 0) {
    if (nem_bor_e_maps == 0) {
      nem_bor_e_maps = malloc(sizeof(ne_stat_struct));
      num_procs = ne_get_proc_count(neid);
      nem_bor_e_maps->stats = malloc(sizeof(int) * num_procs);

      ne_fill_map_status(neid, stat_var, nem_bor_e_maps->stats);
    }
    maps = nem_bor_e_maps;
    maps->neid = neid;
    maps->num  = num_procs;
  } 

  exerrval = 0; /* clear error code */

  if (maps != 0) {
    *stat = maps->stats[proc_id];
    return(EX_NOERR);
  } else {
    return(EX_FATAL);
  }
}

