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
/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: ne_jack.src,v $
 *
 * $Author: gdsjaar $
 *
 * $Date: 2007/10/31 21:39:17 $
 *
 * $Revision: 1.14 $
 *
 *====================================================================*/
/*
 * OVERVIEW
 *
 * This file contains jacket routines written in C for interfacing Fortran
 * NemesisI function calls to the actual C binding for NemsisI.  This code
 * is written explicitly for DARWIN.  In general, these functions handle
 * character-string parameter conventions, convert between
 * column-major-order arrays and row-major-order arrays, and map between
 * array indices beginning at one and array indices beginning at zero.
 *
 */

/* LINTLIBRARY */
#include        <ctype.h>
#include        <string.h>
#include        <stdio.h>
#include        <stdlib.h>
#include        "netcdf.h"
#include        "exodusII.h"
#include        "exodusII_int.h"
#include        "ne_nemesisI.h"
#include        "ne_nemesisI_int.h"


/*
 * The Build64 is for the "normal" SEACAS build which uses compiler
 * options to change reals and integers into 8-byte quantities.  The
 * routines in addrwrap.F are used to down-convert the 8-byte integers
 * into 4-byte integers which then call through to the routines in
 * this file which have a '4' or '4_' appended to the routine name.
 * These routines then call through to the C API routines.
 *
 * If DEFAULT_REAL_INT is defined, then the build is to build a
 * fortran library interface that takes 4-byte ints and either 4-byte
 * or 8-byte floating point (real/double) variables. In this case, the
 * addrwrap routines are not built and a fortran client will call the
 * routines in this file directly.
 *
 */

#if defined(Build64) && !defined(DEFAULT_REAL_INT)
/* 64-bit */
#define real double
#ifdef ADDC_
#define F2C(name) name##4_
#else
#define F2C(name) name##4
#endif

#else
/* 32-bit */
#define real float
#ifdef ADDC_
#define F2C(name) name##_
#else
#define F2C(name) name
#endif
#endif


extern int ncopts;   /* default is (NC_FATAL | NC_VERBOSE) */
extern int exerrval; /* global int that contains a Exodus-specific error code */

/* blank fill C string to make FORTRAN string */
static void
ne_fcdcpy (fstring, fslen, sstring)
    char *fstring;              /* output string to be blank-filled */
    int fslen;                  /* length of output string */
    char *sstring;              /* input string, null-terminated */
{
    int i, len = strlen(sstring);

    for (i = 0; i < len; i++)
        *(fstring + i) = *(sstring + i);
    for (i = len; i < fslen; i++)
        *(fstring + i) = ' ';
}

/* ne_lenstr - string length (w/o trailing blanks) */
static int ne_lenstr (char *string)
{
  char *ptr;

  ptr=string+strlen(string);    /* start at end of string including blanks */
  while (*(--ptr) == ' ');      /* skip blanks */
  return(ptr-string+1);         /* return trimmed length of string */
}

/* copy function used to copy strings and strip trailing blanks */
static void
ne_fstrncpy (target, source, maxlen)
    char *target;               /* space to be copied into */
    char *source;               /* string to be copied */
    int maxlen;                 /* maximum length of *source */
{
    while (maxlen-- && *source != '\0')
        *target++ = *source++;
    while (*(--target) == ' '); /* strip blanks */
    *(++target) = '\0';         /* insert new EOS marker */
}

/* Above are utility functions used below                                   */
/* ======================================================================== */
/* Below are the nemesis API functions                                      */
/*
 * Adding a new function:
 * +  Protect the name with the f2c (uppercase) macro which will add/not add '4' and or '_'
 *    depending on the compilation mode.
 *
 * +  float/double arguments are declared as 'real' which will be replaced with float or double.
 *
 * +  If there are any character arguments 'X', then add an int* argument 'Xlen' at end of argument list
 *    This will contain the length of the passed in character argument.
 *
 * +  Look at existing functions for guidance...
 */

/*
 *  Get initial information from nemesis file
 */
void
F2C(negii)(idne, nproc, nproc_in_f, ftype, ierr, ftypelen)
    int		*idne;	
    int		*nproc;	
    int		*nproc_in_f;	
    char	*ftype;	
    int		ftypelen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  int slen;
  char *file_type;

  /* WARNING: ftypelen SHOULD be 1, but may not be depending on how
              the Fortran programmer passed it. It is best at
              this time to hard code it per NEPII spec. */
  slen = 1;
  if (ftypelen != 1)
  {
    slen = ftypelen;
#if defined(EXODUS_STRING_LENGTH_WARNING)
    sprintf(errmsg,"Warning: file type string length is %d in file id %d\n",
            ftypelen, *idne);
    ex_err("negii",errmsg,EX_MSG);
#endif
  }

  file_type = (char *) malloc((slen+1)*sizeof(char));

  if ((*ierr = ne_get_init_info(*idne, nproc, nproc_in_f, file_type)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to get initial information from file id %d",
	    *idne);
    ex_err("negii",errmsg,EX_MSG);
  }

  if (*ierr == 0)
    ne_fcdcpy (ftype, slen, file_type);

  free(file_type);
}

/*
 *  Write initial information from nemesis file
 */
void
F2C(nepii)(idne, nproc, nproc_in_f, ftype, ierr, ftypelen)
    int		*idne;	
    int		*nproc;	
    int		*nproc_in_f;	
    char	*ftype;	
    int		ftypelen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  int slen;
  char *file_type;

  /* WARNING: ftypelen SHOULD be 1, but may not be depending on how
              the Fortran programmer passed it. It is best at
              this time to hard code it per NEPII spec. */
  slen = 1;
  if (ftypelen != 1)
  {
    slen = ftypelen;
#if defined(EXODUS_STRING_LENGTH_WARNING)
    sprintf(errmsg,"Warning: file type string length is %d in file id %d\n",
            ftypelen, *idne);
    ex_err("nepii",errmsg,EX_MSG);
#endif
  }

  file_type = (char *) malloc((slen+1)*sizeof(char));

  (void) ne_fstrncpy (file_type, ftype, slen);

  if ((*ierr = ne_put_init_info(*idne, *nproc, *nproc_in_f, file_type)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to put initial information in file id %d",
	    *idne);
    ex_err("nepii",errmsg,EX_MSG);
  }

  free(file_type);
}

/*
 * Read initial global information
 */
void
F2C(negig)(idne, nnodes_g, nelems_g, nelem_blks_g, nnode_sets_g, nside_sets_g, ierr)
    int		*idne;	
    int		*nnodes_g;	
    int		*nelems_g;	
    int		*nelem_blks_g;	
    int		*nnode_sets_g;	
    int		*nside_sets_g;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_init_global(*idne, nnodes_g, nelems_g, nelem_blks_g,
                                  nnode_sets_g, nside_sets_g)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read initial global information from file id %d",
	    *idne);
    ex_err("negig",errmsg,EX_MSG);
  }
}

/*
 * Write initial global information
 */
void
F2C(nepig)(idne, nnodes_g, nelems_g, nelem_blks_g, nnode_sets_g, nside_sets_g, ierr)
    int		*idne;	
    int		*nnodes_g;	
    int		*nelems_g;	
    int		*nelem_blks_g;	
    int		*nnode_sets_g;	
    int		*nside_sets_g;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr =  ne_put_init_global(*idne, *nnodes_g, *nelems_g, *nelem_blks_g,
                                   *nnode_sets_g, *nside_sets_g)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store initial global information in file id %d",
	    *idne);
    ex_err("nepig",errmsg,EX_MSG);
  }
}

/*
 * Read load balance parameters
 */
void
F2C(neglbp)(idne, nint_nodes, nbor_nodes, next_nodes, nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, processor, ierr)
    int		*idne;	
    int		*nint_nodes;	
    int		*nbor_nodes;	
    int		*next_nodes;	
    int		*nint_elems;	
    int		*nbor_elems;	
    int		*nnode_cmaps;	
    int		*nelem_cmaps;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_loadbal_param(*idne, nint_nodes, nbor_nodes,
                                    next_nodes, nint_elems, nbor_elems,
                                    nnode_cmaps, nelem_cmaps, *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read load balance parameters from file id %d",
	    *idne);
    ex_err("neglbp",errmsg,EX_MSG);
  }
}

/*
 * Write load balance parameters
 */
void
F2C(neplbp)(idne, nint_nodes, nbor_nodes, next_nodes, nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, processor, ierr)
    int		*idne;	
    int		*nint_nodes;	
    int		*nbor_nodes;	
    int		*next_nodes;	
    int		*nint_elems;	
    int		*nbor_elems;	
    int		*nnode_cmaps;	
    int		*nelem_cmaps;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_loadbal_param(*idne, *nint_nodes, *nbor_nodes,
                                    *next_nodes, *nint_elems, *nbor_elems,
                                    *nnode_cmaps, *nelem_cmaps,
                                    *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store load balance parameters in file id %d",
	    *idne);
    ex_err("neplbp",errmsg,EX_MSG);
  }
}

/*
 * Write concatenated load balance parameters
 */
void
F2C(neplbpc)(idne, nint_nodes, nbor_nodes, next_nodes, nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, ierr)
    int		*idne;	
    int		*nint_nodes;	
    int		*nbor_nodes;	
    int		*next_nodes;	
    int		*nint_elems;	
    int		*nbor_elems;	
    int		*nnode_cmaps;	
    int		*nelem_cmaps;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_loadbal_param_cc(*idne, nint_nodes, nbor_nodes,
                                       next_nodes, nint_elems, nbor_elems,
                                       nnode_cmaps, nelem_cmaps)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store load balance parameters in file id %d",
	    *idne);
    ex_err("neplbpc",errmsg,EX_MSG);
  }
}

/*
 * Read global node set parameters
 */
void
F2C(negnspg)(idne, ns_ids_glob, ns_n_cnt_glob, ns_df_cnt_glob, ierr)
    int		*idne;	
    int		*ns_ids_glob;	
    int		*ns_n_cnt_glob;	
    int		*ns_df_cnt_glob;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_ns_param_global(*idne, ns_ids_glob, ns_n_cnt_glob,
                                      ns_df_cnt_glob)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read global node set parameters from file id %d",
	    *idne);
    ex_err("negnspg",errmsg,EX_MSG);
  }
}

/*
 * Write global node set parameters
 */
void
F2C(nepnspg)(idne, global_ids, global_n_cnts, global_df_cnts, ierr)
    int		*idne;	
    int		*global_ids;	
    int		*global_n_cnts;	
    int		*global_df_cnts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_ns_param_global(*idne, global_ids, global_n_cnts,
                                      global_df_cnts)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store global node set parameters in file id %d",
	    *idne);
    ex_err("nepnspg",errmsg,EX_MSG);
  }
}

/*
 * Read global side set parameters
 */
void
F2C(negsspg)(idne, ss_ids_glob, ss_n_cnt_glob, ss_df_cnt_glob, ierr)
    int		*idne;	
    int		*ss_ids_glob;	
    int		*ss_n_cnt_glob;	
    int		*ss_df_cnt_glob;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_ss_param_global(*idne, ss_ids_glob, ss_n_cnt_glob,
                                      ss_df_cnt_glob)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read global side set parameters from file id %d",
	    *idne);
    ex_err("negsspg",errmsg,EX_MSG);
  }
}

/*
 * Write global side set parameters
 */
void
F2C(nepsspg)(idne, global_ids, global_el_cnts, global_df_cnts, ierr)
    int		*idne;	
    int		*global_ids;	
    int		*global_el_cnts;	
    int		*global_df_cnts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_ss_param_global(*idne, global_ids, global_el_cnts,
                                      global_df_cnts)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store global side set parameters in file id %d",
	    *idne);
    ex_err("nepsspg",errmsg,EX_MSG);
  }
}

/*
 * Read global element block information
 */
void
F2C(negebig)(idne, el_blk_ids, el_blk_cnts, ierr)
    int		*idne;	
    int		*el_blk_ids;	
    int		*el_blk_cnts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read global element block info from file id %d",
	    *idne);
    ex_err("negebig",errmsg,EX_MSG);
  }
}

/*
 * Write global element block information
 */
void
F2C(nepebig)(idne, el_blk_ids, el_blk_cnts, ierr)
    int		*idne;	
    int		*el_blk_ids;	
    int		*el_blk_cnts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to store global element block info in file id %d",
	    *idne);
    ex_err("nepebig",errmsg,EX_MSG);
  }
}

/*
 * Read side set element list and side set side list
 */
void
F2C(negnss)(idne, ss_id, start_side_num, num_sides, ss_elem_list, ss_side_list, ierr)
    int		*idne;	
    int		*ss_id;	
    int		*start_side_num;	
    int		*num_sides;	
    int		*ss_elem_list;	
    int		*ss_side_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_side_set(*idne, *ss_id, *start_side_num, *num_sides,
                                 ss_elem_list, ss_side_list)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to read side set element list from file id %d",
	    *idne);
    ex_err("negnss",errmsg,EX_MSG);
  }
}

/*
 * Write side set element list and side set side list
 */
void
F2C(nepnss)(idne, ss_id, start_side_num, num_sides, ss_elem_list, ss_side_list, ierr)
    int		*idne;	
    int		*ss_id;	
    int		*start_side_num;	
    int		*num_sides;	
    int		*ss_elem_list;	
    int		*ss_side_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_side_set(*idne, *ss_id, *start_side_num, *num_sides,
                                 ss_elem_list, ss_side_list)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to write side set element list to file id %d",
	    *idne);
    ex_err("nepnss",errmsg,EX_MSG);
  }
}

/*
 * Read side set distribution factor
 */
void
F2C(negnssd)(idne, ss_id, start_num, num_df_to_get, ss_df, ierr)
    int		*idne;	
    int		*ss_id;	
    int		*start_num;	
    int		*num_df_to_get;	
    real	*ss_df;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_side_set_df(*idne, *ss_id, *start_num,
                                    *num_df_to_get, ss_df)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to read side set dist factor from file id %d",
	    *idne);
    ex_err("negnssd",errmsg,EX_MSG);
  }
}

/*
 * Write side set distribution factor
 */
void
F2C(nepnssd)(idne, ss_id, start_num, num_df_to_get, ss_df, ierr)
    int		*idne;	
    int		*ss_id;	
    int		*start_num;	
    int		*num_df_to_get;	
    real	*ss_df;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_side_set_df(*idne, *ss_id, *start_num,
                                    *num_df_to_get, ss_df)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to write side set dist factor to file id %d",
	    *idne);
    ex_err("nepnssd",errmsg,EX_MSG);
  }
}

/*
 * Read node set list for a single node set
 */
void
F2C(negnns)(idne, ns_id, start_node_num, num_node, ns_node_list, ierr)
    int		*idne;	
    int		*ns_id;	
    int		*start_node_num;	
    int		*num_node;	
    int		*ns_node_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_node_set(*idne, *ns_id, *start_node_num,
                                 *num_node, ns_node_list)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to read node set node list from file id %d",
	    *idne);
    ex_err("negnns",errmsg,EX_MSG);
  }
}

/*
 * Write node set list for a single node set
 */
void
F2C(nepnns)(idne, ns_id, start_node_num, num_node, ns_node_list, ierr)
    int		*idne;	
    int		*ns_id;	
    int		*start_node_num;	
    int		*num_node;	
    int		*ns_node_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_node_set(*idne, *ns_id, *start_node_num,
                                 *num_node, ns_node_list)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to write node set node list to file id %d",
	    *idne);
    ex_err("nepnns",errmsg,EX_MSG);
  }
}

/*
 * Read node set distribution factor
 */
void
F2C(negnnsd)(idne, ns_id, start_num, num_df_to_get, ns_df, ierr)
    int		*idne;	
    int		*ns_id;	
    int		*start_num;	
    int		*num_df_to_get;	
    real	*ns_df;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_node_set_df(*idne, *ns_id, *start_num,
                                    *num_df_to_get, ns_df)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to read node set dist factor from file id %d",
	    *idne);
    ex_err("negnnsd",errmsg,EX_MSG);
  }
}

/*
 * Write node set distribution factor
 */
void
F2C(nepnnsd)(idne, ns_id, start_num, num_df_to_get, ns_df, ierr)
    int		*idne;	
    int		*ns_id;	
    int		*start_num;	
    int		*num_df_to_get;	
    real	*ns_df;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_node_set_df(*idne, *ns_id, *start_num,
                                    *num_df_to_get, ns_df)) != 0)

  {
    sprintf(errmsg,
	    "Error: failed to write node set dist factor to file id %d",
	    *idne);
    ex_err("nepnnsd",errmsg,EX_MSG);
  }
}

/*
 * Read coordinates of the nodes
 */
void
F2C(negcor)(idne, start_node_num, num_nodes, x_coor, y_coor, z_coor, ierr)
    int		*idne;	
    int		*start_node_num;	
    int		*num_nodes;	
    real	*x_coor;	
    real	*y_coor;	
    real	*z_coor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_coord(*idne, *start_node_num, *num_nodes,
                              x_coor, y_coor, z_coor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read node coordinates from file id %d",
	    *idne);
    ex_err("negcor",errmsg,EX_MSG);
  }
}

/*
 * Write coordinates of the nodes
 */
void
F2C(nepcor)(idne, start_node_num, num_nodes, x_coor, y_coor, z_coor, ierr)
    int		*idne;	
    int		*start_node_num;	
    int		*num_nodes;	
    real	*x_coor;	
    real	*y_coor;	
    real	*z_coor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_coord(*idne, *start_node_num, *num_nodes,
                              x_coor, y_coor, z_coor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write node coordinates to file id %d",
	    *idne);
    ex_err("nepcor",errmsg,EX_MSG);
  }
}

/*
 * Read an element block's connectivity list
 */
void
F2C(negnec)(idne, elem_blk_id, start_elem_num, num_elems, connect, ierr)
    int		*idne;	
    int		*elem_blk_id;	
    int		*start_elem_num;	
    int		*num_elems;	
    int		*connect;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_elem_conn(*idne, *elem_blk_id, *start_elem_num,
                                  *num_elems, connect)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read element block connectivity from file id %d",
	    *idne);
    ex_err("negnec",errmsg,EX_MSG);
  }
}

/*
 * Write an element block's connectivity list
 */
void
F2C(nepnec)(idne, elem_blk_id, start_elem_num, num_elems, connect, ierr)
    int		*idne;	
    int		*elem_blk_id;	
    int		*start_elem_num;	
    int		*num_elems;	
    int		*connect;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_elem_conn(*idne, *elem_blk_id, *start_elem_num,
                                  *num_elems, connect)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write element block connectivity to file id %d",
	    *idne);
    ex_err("negnec",errmsg,EX_MSG);
  }
}

/*
 * Read an element block's attributes
 */
void
F2C(negneat)(idne, elem_blk_id, start_elem_num, num_elems, attrib, ierr)
    int		*idne;	
    int		*elem_blk_id;	
    int		*start_elem_num;	
    int		*num_elems;	
    real	*attrib;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_elem_attr(*idne, *elem_blk_id, *start_elem_num,
                                  *num_elems, attrib)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read element block attribute from file id %d",
	    *idne);
    ex_err("negneat",errmsg,EX_MSG);
  }
}

/*
 * Write an element block's attributes
 */
void
F2C(nepneat)(idne, elem_blk_id, start_elem_num, num_elems, attrib, ierr)
    int		*idne;	
    int		*elem_blk_id;	
    int		*start_elem_num;	
    int		*num_elems;	
    real	*attrib;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_elem_attr(*idne, *elem_blk_id, *start_elem_num,
                                  *num_elems, attrib)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write element block attribute to file id %d",
	    *idne);
    ex_err("nepneat",errmsg,EX_MSG);
  }
}

/*
 * Read the element type for a specific element block
 */
void
F2C(negelt)(idne, elem_blk_id, elem_type, ierr, elem_typelen)
    int		*idne;	
    int		*elem_blk_id;	
    char	*elem_type;	
    int		elem_typelen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  int slen;
  char *etype;

  /* WARNING: ftypelen SHOULD be MAX_STR_LENGTH, but may not be depending
              on how the Fortran programmer passed it. It is best at
              this time to hard code it per NEMESIS spec. */
  slen = MAX_STR_LENGTH;
  if (elem_typelen != MAX_STR_LENGTH)
  {
    slen = elem_typelen;
#if defined(EXODUS_STRING_LENGTH_WARNING)
    sprintf(errmsg,"Warning: element type string length is %d in file id %d\n",
            elem_typelen, *idne);
    ex_err("negelt",errmsg,EX_MSG);
#endif
  }

  etype = (char *) malloc((slen+1)*sizeof(char));

  if ((*ierr = ne_get_elem_type(*idne, *elem_blk_id, etype)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read element block type from file id %d",
	    *idne);
    ex_err("negelt",errmsg,EX_MSG);
  }

  if (*ierr == 0)
    ne_fcdcpy (elem_type, slen, etype);

  free(etype);
}

/*
 * Read a variable for an element block
 */
void
F2C(negnev)(idne, time_step, elem_var_index, elem_blk_id, num_elem_this_blk, start_elem_num, num_elem, elem_var_vals, ierr)
    int		*idne;	
    int		*time_step;	
    int		*elem_var_index;	
    int		*elem_blk_id;	
    int		*num_elem_this_blk;	
    int		*start_elem_num;	
    int		*num_elem;	
    real	*elem_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_elem_var(*idne, *time_step, *elem_var_index,
                                  *elem_blk_id, *num_elem_this_blk,
                                  *start_elem_num, *num_elem,
                                  elem_var_vals)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read element block variable from file id %d",
	    *idne);
    ex_err("negnec",errmsg,EX_MSG);
  }
}

/*
 * Write a variable slab for an element block
 */
void
F2C(nepevs)(idne, time_step, elem_var_index, elem_blk_id, start_pos, num_vals, elem_var_vals, ierr)
    int		*idne;	
    int		*time_step;	
    int		*elem_var_index;	
    int		*elem_blk_id;	
    int		*start_pos;	
    int		*num_vals;	
    real	*elem_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_elem_var_slab(*idne, *time_step, *elem_var_index,
                                    *elem_blk_id, *start_pos, *num_vals,
                                    elem_var_vals)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write elem block variable slab to file id %d",
	    *idne);
    ex_err("negnec",errmsg,EX_MSG);
  }
}

/*
 * Read the values of a single nodal variable for a single time step
 */
void
F2C(negnnv)(idne, time_step, nodal_var_index, start_node_num, num_nodes, nodal_vars, ierr)
    int		*idne;	
    int		*time_step;	
    int		*nodal_var_index;	
    int		*start_node_num;	
    int		*num_nodes;	
    real	*nodal_vars;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_nodal_var(*idne, *time_step, *nodal_var_index,
                                  *start_node_num, *num_nodes,
                                  nodal_vars)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read nodal variable from file id %d",
	    *idne);
    ex_err("negnnv",errmsg,EX_MSG);
  }
}

/*
 * Write nodal variable slab
 */
void
F2C(nepnvs)(idne, time_step, nodal_var_index, start_pos, num_vals, nodal_var_vals, ierr)
    int		*idne;	
    int		*time_step;	
    int		*nodal_var_index;	
    int		*start_pos;	
    int		*num_vals;	
    real	*nodal_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_nodal_var_slab(*idne, *time_step, *nodal_var_index,
                                     *start_pos, *num_vals,
                                     nodal_var_vals)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write nodal variable slab to file id %d",
	    *idne);
    ex_err("nepnvs",errmsg,EX_MSG);
  }
}

/*
 * Read the element numbering map
 */
void
F2C(negnenm)(idne, starte, num_ent, elem_map, ierr)
    int		*idne;	
    int		*starte;	
    int		*num_ent;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_elem_num_map(*idne, *starte, *num_ent, elem_map)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read element numbering map from file id %d",
	    *idne);
    ex_err("negnenm",errmsg,EX_MSG);
  }
}

/*
 * Write the element numbering map
 */
void
F2C(nepnenm)(idne, starte, num_ent, elem_map, ierr)
    int		*idne;	
    int		*starte;	
    int		*num_ent;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_elem_num_map(*idne, *starte, *num_ent, elem_map)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write element numbering map to file id %d",
	    *idne);
    ex_err("nepnenm",errmsg,EX_MSG);
  }
}

/*
 * Read the node numbering map
 */
void
F2C(negnnnm)(idne, startn, num_ent, node_map, ierr)
    int		*idne;	
    int		*startn;	
    int		*num_ent;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_n_node_num_map(*idne, *startn, *num_ent, node_map)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read node numbering map from file id %d",
	    *idne);
    ex_err("negnnnm",errmsg,EX_MSG);
  }
}

/*
 * Write the node numbering map
 */
void
F2C(nepnnnm)(idne, startn, num_ent, node_map, ierr)
    int		*idne;	
    int		*startn;	
    int		*num_ent;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_n_node_num_map(*idne, *startn, *num_ent, node_map)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write node numbering map to file id %d",
	    *idne);
    ex_err("nepnnnm",errmsg,EX_MSG);
  }
}

/*
 * Read the node map for a processor
 */
void
F2C(negnm)(idne, node_mapi, node_mapb, node_mape, processor, ierr)
    int		*idne;	
    int		*node_mapi;	
    int		*node_mapb;	
    int		*node_mape;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_node_map(*idne, node_mapi, node_mapb, node_mape,
                               *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read processor node map from file id %d",
	    *idne);
    ex_err("negnm",errmsg,EX_MSG);
  }
}

/*
 * Write a node map for a processor
 */
void
F2C(nepnm)(idne, node_mapi, node_mapb, node_mape, processor, ierr)
    int		*idne;	
    int		*node_mapi;	
    int		*node_mapb;	
    int		*node_mape;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_node_map(*idne, node_mapi, node_mapb, node_mape,
                               *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write processor node map to file id %d",
	    *idne);
    ex_err("nepnm",errmsg,EX_MSG);
  }
}

/*
 * Read the element map for a processor
 */
void
F2C(negem)(idne, elem_mapi, elem_mapb, processor, ierr)
    int		*idne;	
    int		*elem_mapi;	
    int		*elem_mapb;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_elem_map(*idne, elem_mapi, elem_mapb, *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read processor element map from file id %d",
	    *idne);
    ex_err("negem",errmsg,EX_MSG);
  }
}

/*
 * Write the element map for a processor
 */
void
F2C(nepem)(idne, elem_mapi, elem_mapb, processor, ierr)
    int		*idne;	
    int		*elem_mapi;	
    int		*elem_mapb;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_elem_map(*idne, elem_mapi, elem_mapb, *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write processor element map to file id %d",
	    *idne);
    ex_err("nepem",errmsg,EX_MSG);
  }
}

/*
 * Read the communications map parameters for a single processor
 */
void
F2C(negcmp)(idne, ncmap_ids, ncmap_node_cnts, ecmap_ids, ecmap_elem_cnts, processor, ierr)
    int		*idne;	
    int		*ncmap_ids;	
    int		*ncmap_node_cnts;	
    int		*ecmap_ids;	
    int		*ecmap_elem_cnts;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_cmap_params(*idne, ncmap_ids, ncmap_node_cnts,
                                  ecmap_ids, ecmap_elem_cnts, *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read comm map parameters from file id %d",
	    *idne);
    ex_err("negcmp",errmsg,EX_MSG);
  }
}

/*
 * Write the communications map parameters for a single processor
 */
void
F2C(nepcmp)(idne, nmap_ids, nmap_node_cnts, emap_ids, emap_elem_cnts, processor, ierr)
    int		*idne;	
    int		*nmap_ids;	
    int		*nmap_node_cnts;	
    int		*emap_ids;	
    int		*emap_elem_cnts;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_cmap_params(*idne, nmap_ids, nmap_node_cnts,
                                  emap_ids, emap_elem_cnts, *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write comm map parameters to file id %d",
	    *idne);
    ex_err("nepcmp",errmsg,EX_MSG);
  }
}

/*
 * Write the communications map parameters for all processors
 */
void
F2C(nepcmpc)(idne, nmap_ids, nmap_node_cnts, nproc_ptrs, emap_ids, emap_elem_cnts, eproc_ptrs, ierr)
    int		*idne;	
    int		*nmap_ids;	
    int		*nmap_node_cnts;	
    int		*nproc_ptrs;	
    int		*emap_ids;	
    int		*emap_elem_cnts;	
    int		*eproc_ptrs;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_cmap_params_cc(*idne, nmap_ids, nmap_node_cnts,
                                     nproc_ptrs, emap_ids, emap_elem_cnts,
                                     eproc_ptrs)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write comm map parameters to file id %d",
	    *idne);
    ex_err("nepcmpc",errmsg,EX_MSG);
  }
}

/*
 * Read the nodal communications map for a single processor
 */
void
F2C(negncm)(idne, map_id, node_ids, proc_ids, processor, ierr)
    int		*idne;	
    int		*map_id;	
    int		*node_ids;	
    int		*proc_ids;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_node_cmap(*idne, *map_id, node_ids, proc_ids,
                                *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read nodal communications map from file id %d",
	    *idne);
    ex_err("negncm",errmsg,EX_MSG);
  }
}

/*
 * Write the nodal communications map for a single processor
 */
void
F2C(nepncm)(idne, map_id, node_ids, proc_ids, processor, ierr)
    int		*idne;	
    int		*map_id;	
    int		*node_ids;	
    int		*proc_ids;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_node_cmap(*idne, *map_id, node_ids, proc_ids,
                                *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write nodal communications map to file id %d",
	    *idne);
    ex_err("nepncm",errmsg,EX_MSG);
  }
}

/*
 * Read the elemental communications map for a single processor
 */
void
F2C(negecm)(idne, map_id, elem_ids, side_ids, proc_ids, processor, ierr)
    int		*idne;	
    int		*map_id;	
    int		*elem_ids;	
    int		*side_ids;	
    int		*proc_ids;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_get_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids,
                                *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to read elemental comm map from file id %d",
	    *idne);
    ex_err("negecm",errmsg,EX_MSG);
  }
}

/*
 * Write the elemental communications map for a single processor
 */
void
F2C(nepecm)(idne, map_id, elem_ids, side_ids, proc_ids, processor, ierr)
    int		*idne;	
    int		*map_id;	
    int		*elem_ids;	
    int		*side_ids;	
    int		*proc_ids;	
    int		*processor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  if ((*ierr = ne_put_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids,
                                *processor)) != 0)
  {
    sprintf(errmsg,
	    "Error: failed to write elemental comm map to file id %d",
	    *idne);
    ex_err("nepecm",errmsg,EX_MSG);
  }
}

