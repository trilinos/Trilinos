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
 *     ex_leavedef()
 *     ex_catstrn12()
 *     ne_id_lkup()
 *     ex_get_file_type()
 *     ex_put_nemesis_version()
 *     ne_check_file_version()
 *     ex_get_idx()
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

/* Global variables */
char *ne_ret_string;

int ex_leavedef(int exoid,
                const char *call_rout
                )
{
  char errmsg[MAX_ERR_LENGTH];
  int status;
  
  if ((status = nc_enddef(exoid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to end define mode for file id %d",
            exoid);
    ex_err(call_rout, errmsg, exerrval);

    return (EX_FATAL);
  }
  return (EX_NOERR);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
char *ex_catstrn12(char *name,
                 int   num1,
                 int   num2
                 )
{
  const char *func_name="ex_catstrn12";

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
int ne_id_lkup(int exoid, const char *ne_var_name, int64_t *idx, ex_entity_id ne_var_id)
{
  const char    *func_name="ne_id_lkup";

  int      status;
  int      varid, ndims, dimid[1], ret=-1;
  nc_type  var_type;
  size_t   length, start[1];
  int64_t  my_index, begin, end;
  long long id_val;

  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  if ((status = nc_inq_varid(exoid, ne_var_name, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            ne_var_name, exoid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* check if I need the length for this varible */
  if (idx[1] == -1) {
    /* Get the dimension IDs for this variable */
    if ((status = nc_inq_var(exoid, varid, (char *) 0, &var_type, &ndims,
			     dimid, (int *) 0)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
          "Error: failed to find dimension ID for variable \"%s\" in file ID %d",
              ne_var_name, exoid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    /* Get the length of this variable */
    if ((status = nc_inq_dimlen(exoid, dimid[0], &length)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
          "Error: failed to find dimension for variable \"%s\" in file ID %d",
              ne_var_name, exoid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    idx[1] = length;
  } /* End "if (idx[1] == -1)" */

  begin = idx[0];
  end = idx[1];

  /* Find the index by looping over each entry */
  for(my_index=begin; my_index < end; my_index++) {
    start[0] = my_index;
    status = nc_get_var1_longlong(exoid, varid, start, &id_val);

    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable \"%s\" in file ID %d",
              ne_var_name, exoid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    if (id_val == ne_var_id) {
      ret = (int) my_index;
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
int ex_get_file_type(int exoid,
                     char *ftype
                     )
{
  const char   *func_name="ex_get_file_type";

  int  status;
  int     varid;
  int  lftype;

  char    errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  if ((status = nc_inq_varid(exoid, VAR_FILE_TYPE, &varid)) != NC_NOERR) {

    /* If no file type is found, assume parallel */
    ftype[0] = 'p';
    ftype[1] = '\0';

    return (EX_NOERR);
  }

  if ((status = nc_get_var1_int(exoid, varid, NULL, &lftype)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get variable \"%s\" from file ID %d",
            VAR_FILE_TYPE, exoid);
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
int ex_put_nemesis_version(int exoid)
{
  const char  *func_name="ex_put_nemesis_version";
  int    status;
  float  file_ver, api_ver;

  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear the error code */

  file_ver = NEMESIS_FILE_VERSION;
  api_ver  = NEMESIS_API_VERSION;

  /* Check to see if the nemesis file version is already in the file */
  if (nc_get_att_float(exoid, NC_GLOBAL, "nemesis_file_version", &file_ver) != NC_NOERR) {

    /* Output the Nemesis file version */
    if ((status = nc_put_att_float(exoid, NC_GLOBAL, "nemesis_file_version", NC_FLOAT,
				   1, &file_ver)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output nemesis file version in file ID %d",
              exoid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the Nemesis API version */
    if ((status = nc_put_att_float(exoid, NC_GLOBAL, "nemesis_api_version", NC_FLOAT,
				   1, &api_ver)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output nemesis api version in file ID %d",
              exoid);
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
int ne_check_file_version(int exoid)
{
#if 0
  const char  *func_name="ne_check_file_version";

  float  file_ver;

  int    status;
  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0;  /* clear error code */

  /* Get the file version */
  if ((status = nc_get_att_float(exoid, NC_GLOBAL, "nemesis_file_version", &file_ver)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get the nemesis file version from file ID %d",
            exoid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (fabs(NEMESIS_FILE_VERSION-file_ver) > 0.001) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: Nemesis version mismatch in file ID %d!\n", exoid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }
#endif
  return (EX_NOERR);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function gets the index for the given variable at the
 * position given.
 */
/*****************************************************************************/
int ex_get_idx(int exoid, const char *ne_var_name, int64_t *my_index, int pos)
{
  const char  *func_name="ex_get_idx";

  int      status;
  int      varid;
  size_t   start[1], count[1];
#if defined(NC_NETCDF4)
  long long varidx[2];
#else
  int varidx[2];
#endif
  char   errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* set default values for idx */
  my_index[0] = 0;
  my_index[1] = -1;

  /*
   * assume that if there is an error returned, that this
   * means that this is a parallel file, and the index does
   * not exists. This is not an error
   */
  if ((status = nc_inq_varid(exoid, ne_var_name, &varid)) == NC_NOERR) {
    /* check if we are at the beginning of the index vector */
    if (pos == 0) {
      start[0] = pos;
      count[0] = 1;
    } else {
      start[0] = pos - 1;
      count[0] = 2;
    }

#if defined(NC_NETCDF4)
    status = nc_get_vara_longlong(exoid, varid, start, count, varidx);
#else
    status = nc_get_vara_int(exoid, varid, start, count, varidx);
#endif
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable \"%s\" in file ID %d",
              ne_var_name, exoid);
      ex_err(func_name, errmsg, exerrval);
      return -1;
    }

    if (pos == 0) {
      my_index[0] = 0;
      my_index[1] = varidx[0];
    } else {
      my_index[0] = varidx[0];
      my_index[1] = varidx[1];
    }
  }

  return 1;
}
