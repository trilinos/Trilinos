/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
#include "exodusII.h"
#include "exodusII_int.h"
#include "netcdf.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
typedef double       real;
typedef ex_entity_id entity_id;
#ifdef ADDC_
#define F2C(name) name##4_
#else
#define F2C(name) name##4
#endif

#else
/* 32-bit */
typedef float real;
typedef int   entity_id;
#ifdef ADDC_
#define F2C(name) name##_
#else
#define F2C(name) name
#endif
#endif

/* blank fill C string to make FORTRAN string */
static void ex_fcdcpy(char *fstring, size_t fslen, char *sstring)
/* output string to be blank-filled */
/* length of output string */
/* input string, null-terminated */
{
  size_t i, len = strlen(sstring);

  for (i = 0; i < len; i++)
    *(fstring + i) = *(sstring + i);
  for (i = len; i < fslen; i++)
    *(fstring + i) = ' ';
}

/* copy function used to copy strings and strip trailing blanks */
static void ex_fstrncpy(char *target, char *source, size_t maxlen)
/* space to be copied into */
/* string to be copied */
/* maximum length of *source */
{
  while (maxlen-- && *source != '\0')
    *target++ = *source++;
  while (*(--target) == ' ')
    ;                 /* strip blanks */
  *(++target) = '\0'; /* insert new EOS marker */
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
 * +  If there are any character arguments 'X', then add an int* argument 'Xlen' at end of argument
 * list
 *    This will contain the length of the passed in character argument.
 *
 * +  Look at existing functions for guidance...
 */

/*
 *  Get initial information from nemesis file
 */
void F2C(negii)(int *idne, int *nproc, int *nproc_in_f, char *ftype, int *ierr, size_t ftypelen)
{
  size_t slen = 1;
  char * file_type;

  /* WARNING: ftypelen SHOULD be 1, but may not be depending on how
              the Fortran programmer passed it. It is best at
              this time to hard code it per NEPII spec. */
  if (ftypelen != 1) {
#if defined(EXODUS_STRING_LENGTH_WARNING)
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Warning: file type string length is %lu in file id %d\n", ftypelen, *idne);
    ex_err(__func__, errmsg, EX_MSG);
#endif
    slen = ftypelen;
  }

  file_type = (char *)malloc((slen + 1) * sizeof(char));

  if ((*ierr = ex_get_init_info(*idne, nproc, nproc_in_f, file_type)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to get initial information from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }

  if (*ierr == 0)
    ex_fcdcpy(ftype, slen, file_type);

  free(file_type);
}

/*
 *  Write initial information from nemesis file
 */
void F2C(nepii)(int *idne, int *nproc, int *nproc_in_f, char *ftype, int *ierr, size_t ftypelen)
{

  char errmsg[MAX_ERR_LENGTH];

  size_t slen = 1;
  char * file_type;

  /* WARNING: ftypelen SHOULD be 1, but may not be depending on how
              the Fortran programmer passed it. It is best at
              this time to hard code it per NEPII spec. */
  if (ftypelen != 1) {
    slen = ftypelen;
#if defined(EXODUS_STRING_LENGTH_WARNING)
    sprintf(errmsg, "Warning: file type string length is %lu in file id %d\n", ftypelen, *idne);
    ex_err(__func__, errmsg, EX_MSG);
#endif
  }

  file_type = (char *)malloc((slen + 1) * sizeof(char));

  ex_fstrncpy(file_type, ftype, slen);

  if ((*ierr = ex_put_init_info(*idne, *nproc, *nproc_in_f, file_type)) != 0) {
    sprintf(errmsg, "Error: failed to put initial information in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }

  free(file_type);
}

/*
 * Read initial global information
 */
void F2C(negig)(int *idne, void_int *nnodes_g, void_int *nelems_g, void_int *nelem_blks_g,
                void_int *nnode_sets_g, void_int *nside_sets_g, int *ierr)
{
  if ((*ierr = ex_get_init_global(*idne, nnodes_g, nelems_g, nelem_blks_g, nnode_sets_g,
                                  nside_sets_g)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read initial global information from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write initial global information
 */
void F2C(nepig)(int *idne, void_int *nnodes_g, void_int *nelems_g, void_int *nelem_blks_g,
                void_int *nnode_sets_g, void_int *nside_sets_g, int *ierr)
{
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    int64_t *n_nnodes_g     = (int64_t *)nnodes_g;
    int64_t *n_nelems_g     = (int64_t *)nelems_g;
    int64_t *n_nelem_blks_g = (int64_t *)nelem_blks_g;
    int64_t *n_nnode_sets_g = (int64_t *)nnode_sets_g;
    int64_t *n_nside_sets_g = (int64_t *)nside_sets_g;
    *ierr = ex_put_init_global(*idne, *n_nnodes_g, *n_nelems_g, *n_nelem_blks_g, *n_nnode_sets_g,
                               *n_nside_sets_g);
  }
  else {
    int *n_nnodes_g     = (int *)nnodes_g;
    int *n_nelems_g     = (int *)nelems_g;
    int *n_nelem_blks_g = (int *)nelem_blks_g;
    int *n_nnode_sets_g = (int *)nnode_sets_g;
    int *n_nside_sets_g = (int *)nside_sets_g;
    *ierr = ex_put_init_global(*idne, *n_nnodes_g, *n_nelems_g, *n_nelem_blks_g, *n_nnode_sets_g,
                               *n_nside_sets_g);
  }

  if (*ierr != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store initial global information in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read load balance parameters
 */
void F2C(neglbp)(int *idne, void_int *nint_nodes, void_int *nbor_nodes, void_int *next_nodes,
                 void_int *nint_elems, void_int *nbor_elems, void_int *nnode_cmaps,
                 void_int *nelem_cmaps, int *processor, int *ierr)
{
  if ((*ierr = ex_get_loadbal_param(*idne, nint_nodes, nbor_nodes, next_nodes, nint_elems,
                                    nbor_elems, nnode_cmaps, nelem_cmaps, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read load balance parameters from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write load balance parameters
 */
void F2C(neplbp)(int *idne, void_int *nint_nodes, void_int *nbor_nodes, void_int *next_nodes,
                 void_int *nint_elems, void_int *nbor_elems, void_int *nnode_cmaps,
                 void_int *nelem_cmaps, int *processor, int *ierr)
{
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    int64_t *n_nint_nodes  = (int64_t *)nint_nodes;
    int64_t *n_nbor_nodes  = (int64_t *)nbor_nodes;
    int64_t *n_next_nodes  = (int64_t *)next_nodes;
    int64_t *n_nint_elems  = (int64_t *)nint_elems;
    int64_t *n_nbor_elems  = (int64_t *)nbor_elems;
    int64_t *n_nnode_cmaps = (int64_t *)nnode_cmaps;
    int64_t *n_nelem_cmaps = (int64_t *)nelem_cmaps;
    *ierr = ex_put_loadbal_param(*idne, *n_nint_nodes, *n_nbor_nodes, *n_next_nodes, *n_nint_elems,
                                 *n_nbor_elems, *n_nnode_cmaps, *n_nelem_cmaps, *processor);
  }
  else {
    int *n_nint_nodes  = (int *)nint_nodes;
    int *n_nbor_nodes  = (int *)nbor_nodes;
    int *n_next_nodes  = (int *)next_nodes;
    int *n_nint_elems  = (int *)nint_elems;
    int *n_nbor_elems  = (int *)nbor_elems;
    int *n_nnode_cmaps = (int *)nnode_cmaps;
    int *n_nelem_cmaps = (int *)nelem_cmaps;
    *ierr = ex_put_loadbal_param(*idne, *n_nint_nodes, *n_nbor_nodes, *n_next_nodes, *n_nint_elems,
                                 *n_nbor_elems, *n_nnode_cmaps, *n_nelem_cmaps, *processor);
  }
  if (*ierr != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store load balance parameters in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write concatenated load balance parameters
 */
void F2C(neplbpc)(int *idne, void_int *nint_nodes, void_int *nbor_nodes, void_int *next_nodes,
                  void_int *nint_elems, void_int *nbor_elems, void_int *nnode_cmaps,
                  void_int *nelem_cmaps, int *ierr)
{
  if ((*ierr = ex_put_loadbal_param_cc(*idne, nint_nodes, nbor_nodes, next_nodes, nint_elems,
                                       nbor_elems, nnode_cmaps, nelem_cmaps)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store load balance parameters in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read global node set parameters
 */
void F2C(negnspg)(int *idne, void_int *ns_ids_glob, void_int *ns_n_cnt_glob,
                  void_int *ns_df_cnt_glob, int *ierr)
{
  if ((*ierr = ex_get_ns_param_global(*idne, ns_ids_glob, ns_n_cnt_glob, ns_df_cnt_glob)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read global node set parameters from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write global node set parameters
 */
void F2C(nepnspg)(int *idne, void_int *global_ids, void_int *global_n_cnts,
                  void_int *global_df_cnts, int *ierr)
{
  if ((*ierr = ex_put_ns_param_global(*idne, global_ids, global_n_cnts, global_df_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store global node set parameters in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read global side set parameters
 */
void F2C(negsspg)(int *idne, void_int *ss_ids_glob, void_int *ss_n_cnt_glob,
                  void_int *ss_df_cnt_glob, int *ierr)
{

  if ((*ierr = ex_get_ss_param_global(*idne, ss_ids_glob, ss_n_cnt_glob, ss_df_cnt_glob)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read global side set parameters from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write global side set parameters
 */
void F2C(nepsspg)(int *idne, void_int *global_ids, void_int *global_el_cnts,
                  void_int *global_df_cnts, int *ierr)
{
  if ((*ierr = ex_put_ss_param_global(*idne, global_ids, global_el_cnts, global_df_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store global side set parameters in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read global element block information
 */
void F2C(negebig)(int *idne, void_int *el_blk_ids, void_int *el_blk_cnts, int *ierr)
{
  if ((*ierr = ex_get_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read global element block info from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write global element block information
 */
void F2C(nepebig)(int *idne, void_int *el_blk_ids, void_int *el_blk_cnts, int *ierr)
{
  if ((*ierr = ex_put_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to store global element block info in file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read side set element list and side set side list
 */
void F2C(negnss)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
                 void_int *ss_elem_list, void_int *ss_side_list, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_side_set(*idne, *ss_id, st, cnt, ss_elem_list, ss_side_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read side set element list from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write side set element list and side set side list
 */
void F2C(nepnss)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
                 void_int *ss_elem_list, void_int *ss_side_list, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_side_set(*idne, *ss_id, st, cnt, ss_elem_list, ss_side_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write side set element list to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read side set distribution factor
 */
void F2C(negnssd)(int *idne, entity_id *ss_id, void_int *start, void_int *count, real *ss_df,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_side_set_df(*idne, *ss_id, st, cnt, ss_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read side set dist factor from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write side set distribution factor
 */
void F2C(nepnssd)(int *idne, entity_id *ss_id, void_int *start, void_int *count, real *ss_df,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_side_set_df(*idne, *ss_id, st, cnt, ss_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write side set dist factor to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read node set list for a single node set
 */
void F2C(negnns)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
                 void_int *ns_node_list, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_node_set(*idne, *ns_id, st, cnt, ns_node_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read node set node list from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write node set list for a single node set
 */
void F2C(nepnns)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
                 void_int *ns_node_list, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_node_set(*idne, *ns_id, st, cnt, ns_node_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write node set node list to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read node set distribution factor
 */
void F2C(negnnsd)(int *idne, entity_id *ns_id, void_int *start, void_int *count, real *ns_df,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_node_set_df(*idne, *ns_id, st, cnt, ns_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read node set dist factor from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write node set distribution factor
 */
void F2C(nepnnsd)(int *idne, entity_id *ns_id, void_int *start, void_int *count, real *ns_df,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_node_set_df(*idne, *ns_id, st, cnt, ns_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write node set dist factor to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read coordinates of the nodes
 */
void F2C(negcor)(int *idne, void_int *start, void_int *count, real *x_coor, real *y_coor,
                 real *z_coor, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_coord(*idne, st, cnt, x_coor, y_coor, z_coor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read node coordinates from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write coordinates of the nodes
 */
void F2C(nepcor)(int *idne, void_int *start, void_int *count, real *x_coor, real *y_coor,
                 real *z_coor, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_coord(*idne, st, cnt, x_coor, y_coor, z_coor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write node coordinates to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read an element block's connectivity list
 */
void F2C(negnec)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
                 void_int *connect, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_elem_conn(*idne, *elem_blk_id, st, cnt, connect)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read element block connectivity from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write an element block's connectivity list
 */
void F2C(nepnec)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
                 void_int *connect, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_elem_conn(*idne, *elem_blk_id, st, cnt, connect)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write element block connectivity to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read an element block's attributes
 */
void F2C(negneat)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count, real *attrib,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_elem_attr(*idne, *elem_blk_id, st, cnt, attrib)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read element block attribute from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write an element block's attributes
 */
void F2C(nepneat)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count, real *attrib,
                  int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_n_elem_attr(*idne, *elem_blk_id, st, cnt, attrib)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write element block attribute to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element type for a specific element block
 */
void F2C(negelt)(int *idne, entity_id *elem_blk_id, char *elem_type, int *ierr, size_t elem_typelen)
{
  size_t slen = MAX_STR_LENGTH;
  char * etype;

  /* WARNING: ftypelen SHOULD be MAX_STR_LENGTH, but may not be depending
              on how the Fortran programmer passed it. It is best at
              this time to hard code it per NEMESIS spec. */
  if (elem_typelen != MAX_STR_LENGTH) {
#if defined(EXODUS_STRING_LENGTH_WARNING)
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Warning: element type string length is %lu in file id %d\n", elem_typelen,
            *idne);
    ex_err(__func__, errmsg, EX_MSG);
#endif
    slen = elem_typelen;
  }

  etype = (char *)malloc((slen + 1) * sizeof(char));

  if ((*ierr = ex_get_elem_type(*idne, *elem_blk_id, etype)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read element block type from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }

  if (*ierr == 0)
    ex_fcdcpy(elem_type, slen, etype);

  free(etype);
}

/*
 * Read a variable for an element block
 */
void F2C(negnev)(int *idne, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
                 void_int *num_elem_this_blk, void_int *start, void_int *count, real *elem_var_vals,
                 int *ierr)
{
  int64_t st, cnt, ne;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
    ne  = *(int64_t *)num_elem_this_blk;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
    ne  = *(int *)num_elem_this_blk;
  }

  if ((*ierr = ex_get_n_elem_var(*idne, *time_step, *elem_var_index, *elem_blk_id, ne, st, cnt,
                                 elem_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read element block variable from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write a variable slab for an element block
 */
void F2C(nepevs)(int *idne, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
                 void_int *start, void_int *count, real *elem_var_vals, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_elem_var_slab(*idne, *time_step, *elem_var_index, *elem_blk_id, st, cnt,
                                    elem_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write elem block variable slab to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the values of a single nodal variable for a single time step
 */
void F2C(negnnv)(int *idne, int *time_step, int *nodal_var_index, void_int *start, void_int *count,
                 real *nodal_vars, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_get_n_nodal_var(*idne, *time_step, *nodal_var_index, st, cnt, nodal_vars)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read nodal variable from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write nodal variable slab
 */
void F2C(nepnvs)(int *idne, int *time_step, int *nodal_var_index, void_int *start, void_int *count,
                 real *nodal_var_vals, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  if ((*ierr = ex_put_nodal_var_slab(*idne, *time_step, *nodal_var_index, st, cnt,
                                     nodal_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write nodal variable slab to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element numbering map
 */
void F2C(negnenm)(int *idne, void_int *starte, void_int *num_ent, void_int *elem_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)starte;
    cnt = *(int64_t *)num_ent;
  }
  else {
    st  = *(int *)starte;
    cnt = *(int *)num_ent;
  }

  if ((*ierr = ex_get_n_elem_num_map(*idne, st, cnt, elem_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read element numbering map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the element numbering map
 */
void F2C(nepnenm)(int *idne, void_int *starte, void_int *num_ent, void_int *elem_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)starte;
    cnt = *(int64_t *)num_ent;
  }
  else {
    st  = *(int *)starte;
    cnt = *(int *)num_ent;
  }

  if ((*ierr = ex_put_n_elem_num_map(*idne, st, cnt, elem_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write element numbering map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the node numbering map
 */
void F2C(negnnnm)(int *idne, void_int *startn, void_int *num_ent, void_int *node_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)startn;
    cnt = *(int64_t *)num_ent;
  }
  else {
    st  = *(int *)startn;
    cnt = *(int *)num_ent;
  }

  if ((*ierr = ex_get_n_node_num_map(*idne, st, cnt, node_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read node numbering map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the node numbering map
 */
void F2C(nepnnnm)(int *idne, void_int *startn, void_int *num_ent, void_int *node_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idne) & EX_BULK_INT64_API) {
    st  = *(int64_t *)startn;
    cnt = *(int64_t *)num_ent;
  }
  else {
    st  = *(int *)startn;
    cnt = *(int *)num_ent;
  }

  if ((*ierr = ex_put_n_node_num_map(*idne, st, cnt, node_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write node numbering map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the node map for a processor
 */
void F2C(negnm)(int *idne, void_int *node_mapi, void_int *node_mapb, void_int *node_mape,
                int *processor, int *ierr)
{
  if ((*ierr = ex_get_processor_node_maps(*idne, node_mapi, node_mapb, node_mape, *processor)) !=
      0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read processor node map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write a node map for a processor
 */
void F2C(nepnm)(int *idne, void_int *node_mapi, void_int *node_mapb, void_int *node_mape,
                int *processor, int *ierr)
{
  if ((*ierr = ex_put_processor_node_maps(*idne, node_mapi, node_mapb, node_mape, *processor)) !=
      0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write processor node map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element map for a processor
 */
void F2C(negem)(int *idne, void_int *elem_mapi, void_int *elem_mapb, int *processor, int *ierr)
{
  if ((*ierr = ex_get_processor_elem_maps(*idne, elem_mapi, elem_mapb, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read processor element map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the element map for a processor
 */
void F2C(nepem)(int *idne, void_int *elem_mapi, void_int *elem_mapb, int *processor, int *ierr)
{
  if ((*ierr = ex_put_processor_elem_maps(*idne, elem_mapi, elem_mapb, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write processor element map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the communications map parameters for a single processor
 */
void F2C(negcmp)(int *idne, void_int *ncmap_ids, void_int *ncmap_node_cnts, void_int *ecmap_ids,
                 void_int *ecmap_elem_cnts, int *processor, int *ierr)
{
  if ((*ierr = ex_get_cmap_params(*idne, ncmap_ids, ncmap_node_cnts, ecmap_ids, ecmap_elem_cnts,
                                  *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read comm map parameters from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the communications map parameters for a single processor
 */
void F2C(nepcmp)(int *idne, void_int *nmap_ids, void_int *nmap_node_cnts, void_int *emap_ids,
                 void_int *emap_elem_cnts, int *processor, int *ierr)
{
  if ((*ierr = ex_put_cmap_params(*idne, nmap_ids, nmap_node_cnts, emap_ids, emap_elem_cnts,
                                  *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write comm map parameters to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the communications map parameters for all processors
 */
void F2C(nepcmpc)(int *idne, void_int *nmap_ids, void_int *nmap_node_cnts, void_int *nproc_ptrs,
                  void_int *emap_ids, void_int *emap_elem_cnts, void_int *eproc_ptrs, int *ierr)
{
  if ((*ierr = ex_put_cmap_params_cc(*idne, nmap_ids, nmap_node_cnts, nproc_ptrs, emap_ids,
                                     emap_elem_cnts, eproc_ptrs)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write comm map parameters to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the nodal communications map for a single processor
 */
void F2C(negncm)(int *idne, entity_id *map_id, void_int *node_ids, void_int *proc_ids,
                 int *processor, int *ierr)
{
  if ((*ierr = ex_get_node_cmap(*idne, *map_id, node_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read nodal communications map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the nodal communications map for a single processor
 */
void F2C(nepncm)(int *idne, entity_id *map_id, void_int *node_ids, void_int *proc_ids,
                 int *processor, int *ierr)
{
  if ((*ierr = ex_put_node_cmap(*idne, *map_id, node_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write nodal communications map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Read the elemental communications map for a single processor
 */
void F2C(negecm)(int *idne, entity_id *map_id, void_int *elem_ids, void_int *side_ids,
                 void_int *proc_ids, int *processor, int *ierr)
{
  if ((*ierr = ex_get_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to read elemental comm map from file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}

/*
 * Write the elemental communications map for a single processor
 */
void F2C(nepecm)(int *idne, entity_id *map_id, void_int *elem_ids, void_int *side_ids,
                 void_int *proc_ids, int *processor, int *ierr)
{
  if ((*ierr = ex_put_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg, "Error: failed to write elemental comm map to file id %d", *idne);
    ex_err(__func__, errmsg, EX_MSG);
  }
}
