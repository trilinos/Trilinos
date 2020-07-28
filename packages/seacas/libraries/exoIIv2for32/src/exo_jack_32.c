/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*
 * OVERVIEW
 *
 * This file contains jacket routines written in C for interfacing Fortran
 * Exodus function calls to the actual C binding for Exodus.

 * In general, these functions handle
 * character-string parameter conventions, convert between
 * column-major-order arrays and row-major-order arrays, and map between
 * array indices beginning at one and array indices beginning at zero.
 *
 */

/* LINTLIBRARY */
#include "exodusII.h" // for void_int, ex_err, etc
#include "exodusII_int.h"
#include <stddef.h> // for size_t
#include <stdint.h> // for int64_t
#include <stdio.h>  // for snprintf, NULL
#include <stdlib.h> // for free, malloc, calloc
#include <string.h> // for memset, strlen

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
#define entity_id ex_entity_id

#ifdef ADDC_
#define F2C(name, NAME) name##4_
#else
#ifdef _MSC_VER
#define F2C(name, NAME) NAME##4
#else
#define F2C(name, NAME) name##4
#endif
#endif

#else
/* 32-bit */
#define real float
#define entity_id int
#ifdef ADDC_
#define F2C(name, NAME) name##_
#else
#ifdef _MSC_VER
#define F2C(name, NAME) NAME
#else
#define F2C(name, NAME) name
#endif
#endif

#endif /* 64 vs 32 bit build */

/* blank fill C string to make FORTRAN string */
static void ex_fcdcpy(char *fstring, /* output string to be blank-filled */
                      int   fslen,   /* length of output string */
                      char *sstring)
{ /* input string, null-terminated */
  int i, len;

  if (sstring != NULL) {
    len = strlen(sstring);
    if (len > fslen) {
      len = fslen;
    }

    for (i = 0; i < len; i++) {
      *(fstring + i) = *(sstring + i);
    }
    for (i = len; i < fslen; i++) {
      *(fstring + i) = ' ';
    }
  }
  else {
    for (i = 0; i < fslen; i++) {
      *(fstring + i) = ' ';
    }
  }
}

/* copy function used to copy strings and strip trailing blanks */
static void ex_fstrncpy(char *target, /* space to be copied into */
                        char *source, /* string to be copied */
                        int   maxlen)
{ /* maximum length of *source */
  int len = maxlen;

  while (len-- && *source != '\0') {
    *target++ = *source++;
  }

  len = maxlen;
  while (len-- && *(--target) == ' ') {
    ; /* strip blanks */
  }
  *(++target) = '\0'; /* insert new EOS marker */
}

/* copy function used to copy strings terminated with blanks */
static void ex_nstrncpy(char *target, /* space to be copied into */
                        char *source, /* string to be copied */
                        int   maxlen)
{ /* maximum length of *source */
  while (maxlen-- && *source != ' ') {
    *target++ = *source++;
  }
  *target = '\0';
}

/* Above are utility functions used below                                   */
/* ======================================================================== */
/* Below are the exodus API functions                                       */
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
 * create an EXODUS II file
 */
int F2C(excre, EXCRE)(char *path, int *clobmode, int *cpu_word_size, int *io_word_size, int *ierr,
                      int pathlen)
{
  char *name;
  int   idexo;

  if (!(name = malloc((pathlen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return (EX_FATAL);
  }
  (void)ex_nstrncpy(name, path, pathlen);

  if ((idexo = ex_create(name, *clobmode, cpu_word_size, io_word_size)) != EX_FATAL) {
    free(name);
    *ierr = 0;
    return (idexo);
  }
  free(name);
  *ierr = exerrval;
  return (EX_FATAL);
}

/*
 * open an EXODUS II file
 */
int F2C(exopen, EXOPEN)(char *path, int *mode, int *cpu_word_size, int *io_word_size,
                        float *version, /* This is float always; not real */
                        int *ierr, int pathlen)
{
  char *name;
  int   idexo;

  if (!(name = malloc((pathlen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return (EX_FATAL);
  }
  (void)ex_nstrncpy(name, path, pathlen);
  if ((idexo = ex_open(name, *mode, cpu_word_size, io_word_size, version)) != EX_FATAL) {
    free(name);
    *ierr = 0;
    return (idexo);
  }
  free(name);
  *ierr = EX_FATAL;
  return (EX_FATAL);
}

/*
 * close an EXODUS II file
 */
void F2C(exclos, EXCLOS)(int *idexo, int *ierr) { *ierr = ex_close(*idexo); }

/*
 * update an EXODUS II file
 */
void F2C(exupda, EXUPDA)(int *idexo, int *ierr) { *ierr = ex_update(*idexo); }

/*
 * write initialization parameters
 */
void F2C(expini, EXPINI)(int *idexo, char *title, void_int *num_dim, void_int *num_nodes,
                         void_int *num_elem, void_int *num_elem_blk, void_int *num_node_sets,
                         void_int *num_side_sets, int *ierr, int titlelen)
{
  int   slen;
  char *name;

  slen = MAX_LINE_LENGTH; /* max line size */
  if (titlelen != MAX_LINE_LENGTH) {
    slen = titlelen;
  }
  name = malloc((slen + 1) * sizeof(char));
  (void)ex_fstrncpy(name, title, slen);

  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int64_t *n_dim       = num_dim;
    int64_t *n_nodes     = num_nodes;
    int64_t *n_elem      = num_elem;
    int64_t *n_elem_blk  = num_elem_blk;
    int64_t *n_node_sets = num_node_sets;
    int64_t *n_side_sets = num_side_sets;

    *ierr = ex_put_init(*idexo, name, *n_dim, *n_nodes, *n_elem, *n_elem_blk, *n_node_sets,
                        *n_side_sets);
  }
  else {
    int *n_dim       = num_dim;
    int *n_nodes     = num_nodes;
    int *n_elem      = num_elem;
    int *n_elem_blk  = num_elem_blk;
    int *n_node_sets = num_node_sets;
    int *n_side_sets = num_side_sets;

    *ierr = ex_put_init(*idexo, name, *n_dim, *n_nodes, *n_elem, *n_elem_blk, *n_node_sets,
                        *n_side_sets);
  }
  free(name);
}

/*
 * read initialization parameters
 */
void F2C(exgini, EXGINI)(int *idexo, char *title, void_int *num_dim, void_int *num_nodes,
                         void_int *num_elem, void_int *num_elem_blk, void_int *num_node_sets,
                         void_int *num_side_sets, int *ierr, int titlelen)
{
  int   slen;
  char *name;

  *ierr = 0;
  slen  = MAX_LINE_LENGTH; /* max line size */
  if (titlelen != MAX_LINE_LENGTH) {
    slen = titlelen;
  }
  if (!(name = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  memset(name, 0, slen + 1);

  *ierr = ex_get_init(*idexo, name, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
                      num_side_sets);

  ex_fcdcpy(title, slen, name);
  free(name);
}

/*
 * write QA records
 */
void F2C(expqa, EXPQA)(int *idexo, int *num_qa_records, char *qa_record, int *ierr,
                       int qa_recordlen)
{
  char   errmsg[MAX_ERR_LENGTH];
  char **sptr; /* internal string pointer array for malloc
                * use */
  int i, ii, iii, slen, alen;

  *ierr = 0; /* default no error */

  slen = MAX_STR_LENGTH; /* max str size */
  if (qa_recordlen != MAX_STR_LENGTH) {
    slen = qa_recordlen;
  }
  alen = 4; /* qa records are 4 strings deep */

  /* Allocate space for the name ptr array */
  if (!(sptr = malloc(((*num_qa_records) * alen + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /*
   * Allocate space for each of the strings, where size = slen, place
   * ptr into str ptr array,  and Copy Fortran qa records to staging
   * space
   */
  iii = 0; /* offset counter */
  for (i = 0; i < *num_qa_records; i++) {
    for (ii = 0; ii < alen; ii++) {
      *(sptr + iii) = malloc((slen + 1) * sizeof(char));
      if (*(sptr + iii) == 0) {
        free(sptr); /* free up array ptr space */
        *ierr = EX_MEMFAIL;
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "Error: failed to allocate space for qa record %d for file id %d", i, *idexo);
        ex_err_fn(*idexo, __func__, errmsg, EX_MEMFAIL);
        return;
      }
      /* copy fortran string into allocated space */
      ex_fstrncpy(*(sptr + iii), qa_record + iii * qa_recordlen, slen);
      iii++; /* bump char array pointer */
    }
  }
  *(sptr + iii) = 0; /* set last pointer to null */

  if (ex_put_qa(*idexo, *num_qa_records, (void *)sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }

  /* Free up the space we used */
  iii = 0;
  for (i = 0; i < *num_qa_records; i++) {
    for (ii = 0; ii < alen; ii++) {
      free(*(sptr + iii)); /* First free up string space */
      iii++;
    }
  }
  free(sptr); /* Then free up array ptr space */
}

/*
 * read QA records
 */
void F2C(exgqa, EXGQA)(int *idexo, char *qa_record, int *ierr, int qa_recordlen)
{
  int    num_qa_records;
  char **sptr; /* internal string pointer array for malloc
                * use */
  int i, ii, iii, slen, alen;

  *ierr = 0; /* default no error */

  slen = MAX_STR_LENGTH; /* max str size */
  if (qa_recordlen != MAX_STR_LENGTH) {
    slen = qa_recordlen;
  }
  alen = 4; /* qa records are 4 strings deep */

  /* do Exodus C call to find out how many qa records are avail */
  num_qa_records = ex_inquire_int(*idexo, EX_INQ_QA);
  if (num_qa_records < 0) {
    *ierr = EX_FATAL;
    return;
  }

  /* Allocate space for the QA string ptr array */
  if (!(sptr = malloc((num_qa_records * alen + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /*
   * Step 1: Allocate space for each of the strings, where size = slen,
   * place string ptr into str ptr array. Step 2: Call routine to get
   * qa records Step 3: Copy C qa records to passed Fortran array space
   */

  iii = 0;                               /* offset counter */
  for (i = 0; i < num_qa_records; i++) { /* pointer allocation loop */
    for (ii = 0; ii < alen; ii++) {
      *(sptr + iii) = malloc((slen + 1) * sizeof(char));
      if (*(sptr + iii) == 0) {
        free(sptr);
        *ierr = EX_MEMFAIL;
        return;
      }
      iii++; /* bump char array pointer */
    }
  }
  *(sptr + iii) = 0; /* null out last pointer */

  /* do Exodus C call to get qa records */
  if (ex_get_qa(*idexo, (void *)sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    goto error_ret;
  }
  iii = 0;                               /* offset counter */
  for (i = 0; i < num_qa_records; i++) { /* string copy loop */
    for (ii = 0; ii < alen; ii++) {
      /* copy fortran string into allocated space */
      ex_fcdcpy(qa_record + iii * qa_recordlen, slen, *(sptr + iii));
      iii++; /* bump char array pointer */
    }
  }

error_ret:
  /* Free up the space we used */
  iii = 0;
  for (i = 0; i < num_qa_records; i++) {
    for (ii = 0; ii < alen; ii++) {
      free(*(sptr + iii)); /* First free up string space */
      iii++;
    }
  }
  free(sptr); /* Then free up array ptr space */
}

/*
 * write information records
 */
void F2C(expinf, EXPINF)(int *idexo, int *num_info, char *info, int *ierr, int infolen)
{
  char **aptr; /* internal string array pointer for malloc
                * use */
  char *sptr;  /* internal string pointer for malloc use */
  int   i, slen;

  *ierr = 0;               /* default no error */
  slen  = MAX_LINE_LENGTH; /* max str size */
  if (infolen != MAX_LINE_LENGTH) {
    slen = infolen;
  }
  /* Allocate space for the string ptr array */
  if (!(aptr = malloc(((*num_info) + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the info records */
  if (!(sptr = malloc(*num_info * (slen + 1) * sizeof(char)))) {
    free(aptr); /* Free up string ptr array */
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy Fortran info records to staging space */
  for (i = 0; i < *num_info; i++) {
    *(aptr + i) = sptr + i * (slen + 1);                /* put address into ptr array */
    ex_fstrncpy(*(aptr + i), info + i * infolen, slen); /* copy string into
                                                         * buffer */
  }
  *(aptr + i) = 0; /* null out last ptr */
  if (ex_put_info(*idexo, *num_info, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* Free up string staging area */
    free(aptr); /* Free up string ptr array */
    return;
  }
  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * read information records
 */
void F2C(exginf, EXGINF)(int *idexo, char *info, int *ierr, int infolen)
{
  char **aptr; /* internal string array pointer for malloc
                * use */
  char *sptr;  /* internal string pointer for malloc use */
  int   i, slen, num_info;

  *ierr = 0; /* default no error */

  /* do exodus C call to find out how many info records are avail */
  num_info = ex_inquire_int(*idexo, EX_INQ_INFO);
  if (num_info < 0) {
    *ierr = EX_FATAL;
    return;
  }

  slen = MAX_LINE_LENGTH; /* max str size */
  if (infolen != MAX_LINE_LENGTH) {
    slen = infolen;
  }
  /*
   * Step 1: Allocate space for string ptr array Step 2: Allocate space
   * for info record strings, and put pointers into str ptr array Step
   * 3: Do Exodus call to get records Step 4: Copy strings into
   * passed Fortran buffer space
   */

  /* Allocate space for the string ptr array */
  if (!(aptr = malloc((num_info + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate block of space for info strings */
  if (!(sptr = malloc(num_info * (slen + 1) * sizeof(char)))) {
    free(aptr); /* Free up string ptr array */
    *ierr = EX_MEMFAIL;
    return;
  }
  for (i = 0; i < num_info; i++) { /* Put pointers to the info records in ptr
                                    * array */
    *(aptr + i) = sptr + i * (slen + 1);
  }                /* put ptr in string ptr
                    * array */
  *(aptr + i) = 0; /* null out last pointer */

  /* Do exodus call to get info records */
  if (ex_get_info(*idexo, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr);
    free(aptr);
    return;
  }
  for (i = 0; i < num_info; i++) {                    /* Copy Fortran info records to
                                                       * staging space */
    ex_fcdcpy(info + i * infolen, slen, *(aptr + i)); /* copy string into
                                                       * buffer */
    /** printf("[exginf] rec: %d , %s\n",i,*(aptr+i)); **/
  }

  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * write nodal coordinates
 */
void F2C(expcor, EXPCOR)(int *idexo, real *x_coor, real *y_coor, real *z_coor, int *ierr)
{
  *ierr = ex_put_coord(*idexo, x_coor, y_coor, z_coor);
}

/*
 * read nodal coordinates
 */
void F2C(exgcor, EXGCOR)(int *idexo, real *x_coor, real *y_coor, real *z_coor, int *ierr)
{
  *ierr = ex_get_coord(*idexo, x_coor, y_coor, z_coor);
}

/*
 * write coordinate names
 */
void F2C(expcon, EXPCON)(int *idexo, char *coord_names, int *ierr, int coord_nameslen)
{
  char **aptr; /* internal array of string pointers for
                * malloc use */
  char *sptr;  /* internal string pointer for malloc use */
  int   i, ndim, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (coord_nameslen < slen) {
    slen = coord_nameslen;
  }
  /* do Exodus C call to find out how many dimensions  */
  ndim = ex_inquire_int(*idexo, EX_INQ_DIM);
  if (ndim < 0) {
    *ierr = EX_FATAL;
    return;
  }

  /* Allocate space for the name ptr array */
  if (!(aptr = malloc((ndim + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /*
   * Allocate a block of space for the strings, where size = slen,
   * place ptrs into str ptr array,  and Copy Fortran coordinate names
   * to staging space
   */

  if (!(sptr = malloc(ndim * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr);
    return;
  }
  for (i = 0; i < ndim; i++) {
    *(aptr + i) = sptr + i * (slen + 1);
    /* copy fortran string into allocated space */
    ex_fstrncpy(*(aptr + i), coord_names + i * coord_nameslen, slen);
  }
  *(aptr + i) = 0; /* set last pointer to null */

  if (ex_put_coord_names(*idexo, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  /* Free up the space we used */
  free(sptr); /* First free up string space */
  free(aptr); /* Then free up array ptr space */
}
/*
 * read coordinate names
 */
void F2C(exgcon, EXGCON)(int *idexo, char *coord_names, int *ierr, int coord_nameslen)
{
  char **aptr; /* internal string array pointer for malloc
                * use */
  char *sptr;  /* internal string pointer for malloc use */
  int   ndim;
  int   i, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (coord_nameslen < slen) {
    slen = coord_nameslen;
  }
  /* do Exodus C call to find out how many dimensions */
  ndim = ex_inquire_int(*idexo, EX_INQ_DIM);
  if (ndim < 0) {
    *ierr = EX_FATAL;
    return;
  }

  /* allocate memory to stage the coordinate name ptrs into */
  if (!(aptr = malloc((ndim + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* allocate a block of memory to stage the coordinate names into */
  if (!(sptr = malloc(ndim * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr); /* free up array ptr space */
    return;
  }
  for (i = 0; i < ndim; i++) { /* put pointers to staging space into ptr
                                * array */
    *(aptr + i) = sptr + i * (slen + 1);
  }

  /* do Exodus C call to get coord name records */
  if (ex_get_coord_names(*idexo, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up string space */
    free(aptr); /* free up array ptr space */
    return;
  }
  /* copy C strings to Fortran arrays */
  memset(coord_names, 0, ndim * coord_nameslen);
  for (i = 0; i < ndim; i++) {
    ex_fcdcpy(coord_names + i * coord_nameslen, slen, *(aptr + i)); /* copy and blank fill */
  }

  free(sptr); /* Free up string buffer space */
  free(aptr);
}

/*
 * write element order map
 */
void F2C(expmap, EXPMAP)(int *idexo, void_int *elem_map, int *ierr)
{
  *ierr = ex_put_map(*idexo, elem_map);
}

/*
 * read element order map
 */
void F2C(exgmap, EXGMAP)(int *idexo, void_int *elem_map, int *ierr)
{
  *ierr = ex_get_map(*idexo, elem_map);
}

/*
 * write concatenated element block parameters
 */
void F2C(expclb, EXPCLB)(int *idexo, void_int *elem_blk_id, char *elem_type,
                         void_int *num_elem_this_blk, void_int *num_nodes_per_elem,
                         void_int *num_attr, int *create_maps, int *ierr, int elem_typelen)
{
  int num_elem_blk;

  char **aptr; /* ptr to temp staging space for string array
                * ptrs */
  char * sptr; /* ptr to temp staging space for strings */
  int    i;
  size_t slen;

  *ierr = 0; /* default no error */

  num_elem_blk = ex_inquire_int(*idexo, EX_INQ_ELEM_BLK);
  if (num_elem_blk < 0) {
    *ierr = EX_FATAL;
    return;
  }

  slen = MAX_STR_LENGTH; /* max str size */
  if (elem_typelen != MAX_STR_LENGTH) {
    slen = elem_typelen;
  }
  /* allocate memory for pointer array */
  if (!(aptr = malloc((num_elem_blk + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* allocate memory to stage the element type name into */
  if (!(sptr = malloc(num_elem_blk * (slen + 1) * sizeof(char)))) {
    free(aptr);
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy element type names from Fortran array to staging area */
  for (i = 0; i < num_elem_blk; i++) {
    *(aptr + i) = sptr + i * (slen + 1);                          /* put address into ptr array */
    ex_fstrncpy(*(aptr + i), elem_type + i * elem_typelen, slen); /* copy string into
                                                                   * buffer */
  }
  *(aptr + i) = 0; /* null out last ptr */

  if (ex_put_concat_elem_block(*idexo, elem_blk_id, aptr, num_elem_this_blk, num_nodes_per_elem,
                               num_attr, *create_maps) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr);
  free(aptr);
}

/*
 * write element block parameters
 */
void F2C(expelb, EXPELB)(int *idexo, entity_id *elem_blk_id, char *elem_type,
                         void_int *num_elem_this_blk, void_int *num_nodes_per_elem,
                         void_int *num_attr, int *ierr, int elem_typelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0; /* default no error */

  slen = MAX_STR_LENGTH; /* max str size */
  if (elem_typelen != MAX_STR_LENGTH) {
    slen = elem_typelen;
  }
  /* allocate memory to stage the element type name into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy element type names from Fortran array to staging area */
  ex_fstrncpy(sptr, elem_type, slen);

  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int64_t *n_elem_this_blk  = num_elem_this_blk;
    int64_t *n_nodes_per_elem = num_nodes_per_elem;
    int64_t *n_attr           = num_attr;

    *ierr = ex_put_block(*idexo, EX_ELEM_BLOCK, *elem_blk_id, sptr, *n_elem_this_blk,
                         *n_nodes_per_elem, 0, 0, *n_attr);
  }
  else {
    int *n_elem_this_blk  = num_elem_this_blk;
    int *n_nodes_per_elem = num_nodes_per_elem;
    int *n_attr           = num_attr;

    *ierr = ex_put_block(*idexo, EX_ELEM_BLOCK, *elem_blk_id, sptr, *n_elem_this_blk,
                         *n_nodes_per_elem, 0, 0, *n_attr);
  }
  free(sptr);
}

/*
 * read element block parameters
 */
void F2C(exgelb, EXGELB)(int *idexo, entity_id *elem_blk_id, char *elem_type,
                         void_int *num_elem_this_blk, void_int *num_nodes_per_elem,
                         void_int *num_attr, int *ierr, int elem_typelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0;

  slen = MAX_STR_LENGTH; /* max str size */
  if (elem_typelen != MAX_STR_LENGTH) {
    slen = elem_typelen;
  }
  /* allocate memory to stage the element type names into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  if (ex_get_block(*idexo, EX_ELEM_BLOCK, *elem_blk_id, sptr, num_elem_this_blk, num_nodes_per_elem,
                   NULL, NULL, num_attr) == EX_FATAL) {
    *ierr = EX_FATAL;
    return;
  }
  /* Copy element type name from staging area to Fortran array */
  memset(elem_type, 0, elem_typelen);
  ex_fcdcpy(elem_type, slen, sptr);
  free(sptr);
}

/*
 * read element blocks IDs
 */
void F2C(exgebi, EXGEBI)(int *idexo, void_int *elem_blk_ids, int *ierr)
{
  *ierr = ex_get_ids(*idexo, EX_ELEM_BLOCK, elem_blk_ids);
}

/*
 * write element block connectivity
 */
void F2C(expelc, EXPELC)(int *idexo, entity_id *elem_blk_id, void_int *connect, int *ierr)
{
  *ierr = ex_put_conn(*idexo, EX_ELEM_BLOCK, *elem_blk_id, connect, NULL, NULL);
}

/*
 * read element block connectivity
 */
void F2C(exgelc, EXGELC)(int *idexo, entity_id *elem_blk_id, void_int *connect, int *ierr)
{
  *ierr = ex_get_conn(*idexo, EX_ELEM_BLOCK, *elem_blk_id, connect, NULL, NULL);
}

/*
 * write entity count-per-polyhedra information for nsided block
 */
void F2C(expecpp, EXPECPP)(int *idexo, int *obj_type, entity_id *elem_blk_id, int *counts,
                           int *ierr)
{
  *ierr =
      ex_put_entity_count_per_polyhedra(*idexo, (ex_entity_type)*obj_type, *elem_blk_id, counts);
}

/*
 * read entity count-per-polyhedra information for nsided block
 */
void F2C(exgecpp, EXGECPP)(int *idexo, int *obj_type, entity_id *elem_blk_id, int *counts,
                           int *ierr)
{
  *ierr =
      ex_get_entity_count_per_polyhedra(*idexo, (ex_entity_type)*obj_type, *elem_blk_id, counts);
}

/*
 * write element block attributes
 */
void F2C(expeat, EXPEAT)(int *idexo, entity_id *elem_blk_id, real *attrib, int *ierr)
{
  *ierr = ex_put_attr(*idexo, EX_ELEM_BLOCK, *elem_blk_id, attrib);
}

/*
 * read element block attributes
 */
void F2C(exgeat, EXGEAT)(int *idexo, entity_id *elem_blk_id, real *attrib, int *ierr)
{
  *ierr = ex_get_attr(*idexo, EX_ELEM_BLOCK, *elem_blk_id, attrib);
}

/*
 * read element block attribute names
 */
void F2C(exgean, EXGEAN)(int *idexo, entity_id *elem_blk_id, int *num_attr, char *names, int *ierr,
                         int nameslen)
{
  char **aptr; /* ptr to temp staging space for string array
                * ptrs */
  char *sptr;  /* ptr to temp staging space for strings */
  int   i, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (nameslen < slen) {
    slen = nameslen;
  }
  /* allocate memory to for pointer array */
  if (!(aptr = malloc((*num_attr + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_attr * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr); /* Free up string ptr array */
    return;
  }
  for (i = 0; i < *num_attr; i++) {
    *(aptr + i) = sptr + i * (slen + 1); /* put address into ptr array */
  }
  *(aptr + i) = 0; /* null out last ptr */

  *ierr = 0;
  if (ex_get_attr_names(*idexo, EX_ELEM_BLOCK, *elem_blk_id, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    free(aptr);
    return;
  }
  /* Copy Fortran names from staging space */
  memset(names, 0, *num_attr * nameslen);
  for (i = 0; i < *num_attr; i++) {
    ex_fcdcpy(names + i * nameslen, slen, *(aptr + i)); /* copy str into Fortran
                                                         * buffer */
  }

  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * write element block attribute names
 */
void F2C(expean, EXPEAN)(int *idexo, entity_id *elem_blk_id, int *num_attr, char *names, int *ierr,
                         int nameslen)
{
  char **aptr; /* ptr to temp staging space for string array
                * ptrs */
  char *sptr;  /* ptr to temp staging space for strings */
  int   i, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (nameslen < slen) {
    slen = nameslen;
  }
  /* allocate memory to for pointer array */
  if (!(aptr = malloc((*num_attr + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_attr * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr); /* Free up string ptr array */
    return;
  }
  /* Copy Fortran names to staging space */
  for (i = 0; i < *num_attr; i++) {
    *(aptr + i) = sptr + i * (slen + 1);                  /* put address into ptr array */
    ex_fstrncpy(*(aptr + i), names + i * nameslen, slen); /* copy string into
                                                           * buffer */
  }
  *(aptr + i) = 0; /* null out last ptr */

  *ierr = 0;
  if (ex_put_attr_names(*idexo, EX_ELEM_BLOCK, *elem_blk_id, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * write object names
 */
void F2C(expnams, EXPNAMS)(int *idexo, int *type, int *num_obj, char *names, int *ierr,
                           int nameslen)
{
  char **aptr; /* ptr to temp staging space for string array
                * ptrs */
  char *sptr;  /* ptr to temp staging space for strings */
  int   i, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (nameslen < slen) {
    slen = nameslen;
  }
  /* allocate memory for pointer array */
  if (!(aptr = malloc((*num_obj + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_obj * (slen + 1) * sizeof(char)))) {
    free(aptr); /* Free up string ptr array */
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy Fortran names to staging space */
  for (i = 0; i < *num_obj; i++) {
    *(aptr + i) = sptr + i * (slen + 1);                  /* put address into ptr array */
    ex_fstrncpy(*(aptr + i), names + i * nameslen, slen); /* copy string into
                                                           * buffer */
  }
  *(aptr + i) = 0; /* null out last ptr */
  /* do Exodus C call to write results variables names */
  if (ex_put_names(*idexo, (ex_entity_type)*type, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * write single object names
 */
void F2C(expnam, EXPNAM)(int *idexo, int *type, int *id, char *name, int *ierr, int namelen)
{
  char *sptr; /* ptr to temp staging space for string */
  int   slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (namelen < slen) {
    slen = namelen;
  }

  /* Allocate staging space for the variable name */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy Fortran names to staging space */
  ex_fstrncpy(sptr, name, slen); /* copy string into buffer */

  /* do Exodus C call to write results variables names */
  if (ex_put_name(*idexo, (ex_entity_type)*type, *id, sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr); /* Free up string staging area */
}

/*
 * read object names
 */
void F2C(exgnams, EXGNAMS)(int *idexo, int *type, int *num_obj, char *names, int *ierr,
                           int nameslen)
{
  char **aptr; /* ptr to temp staging space for string array
                * ptrs */
  char *sptr;  /* ptr to temp staging space for strings */
  int   i, slen;

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (nameslen < slen) {
    slen = nameslen;
  }
  /* allocate memory to for pointer array */
  if (!(aptr = malloc((*num_obj + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_obj * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr); /* Free up string ptr array */
    return;
  }
  for (i = 0; i < *num_obj; i++) {
    *(aptr + i) = sptr + i * (slen + 1); /* put address into ptr array */
  }
  *(aptr + i) = 0; /* null out last ptr */

  /* do Exodus C call to read results variables names */
  if (ex_get_names(*idexo, (ex_entity_type)*type, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    free(aptr);
    return;
  }
  /* Copy Fortran names from staging space */
  memset(names, 0, *num_obj * nameslen);
  for (i = 0; i < *num_obj; i++) {
    ex_fcdcpy(names + i * nameslen, slen, *(aptr + i)); /* copy str into Fortran
                                                         * buffer */
  }

  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * read single object name
 */
void F2C(exgnam, EXGNAM)(int *idexo, int *type, int *id, char *name, int *ierr, int namelen)
{
  char *sptr; /* ptr to temp staging space for string */
  int   slen;
  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (namelen < slen) {
    slen = namelen;
  }
  /* Allocate staging space for the object name */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* do Exodus C call to read object  name */
  if (ex_get_name(*idexo, *type, *id, sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    return;
  }
  /* Copy Fortran variable names to staging space */
  memset(name, 0, namelen);
  ex_fcdcpy(name, slen, sptr); /* copy string into Fortran buffer */

  free(sptr); /* Free up string staging area */
}

/*
 * write property array names
 */
void F2C(exppn, EXPPN)(int *idexo, int *obj_type, int *num_props, char *prop_names, int *ierr,
                       int prop_nameslen)
{
  char **aptr; /* internal string array pointer for malloc
                * use */
  char *sptr;  /* internal string pointer for malloc use */
  int   i, slen;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_nameslen < slen) {
    slen = prop_nameslen;
  }
  /* Allocate space for the name ptr array */
  if (!(aptr = malloc((*num_props + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /*
   * Allocate a block of space for the strings, where size = slen,
   * place ptrs into str ptr array,  and Copy Fortran coordinate names
   * to staging space
   */

  if (!(sptr = malloc((*num_props) * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr);
    return;
  }
  for (i = 0; i < *num_props; i++) {
    *(aptr + i) = sptr + i * (slen + 1);
    /* copy fortran string into allocated space */
    ex_fstrncpy(*(aptr + i), prop_names + i * prop_nameslen, slen);
  }
  *(aptr + i) = 0; /* set last pointer to null */

  if (ex_put_prop_names(*idexo, (ex_entity_type)*obj_type, *num_props, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  /* Free up the space we used */
  free(sptr); /* First free up string space */
  free(aptr); /* Then free up array ptr space */
}

/*
 * read property array names
 */
void F2C(exgpn, EXGPN)(int *idexo, int *obj_type, char *prop_names, int *ierr, int prop_nameslen)
{
  char   errmsg[MAX_ERR_LENGTH];
  char **aptr;     /* internal string array pointer for malloc
                    * use */
  char *     sptr; /* internal string pointer for malloc use */
  int        i, slen;
  ex_inquiry inq_code;
  int        num_props;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_nameslen < slen) {
    slen = prop_nameslen;
  }
  switch ((ex_entity_type)*obj_type) {
  case EX_ELEM_BLOCK: inq_code = EX_INQ_EB_PROP; break;
  case EX_NODE_SET: inq_code = EX_INQ_NS_PROP; break;
  case EX_SIDE_SET: inq_code = EX_INQ_SS_PROP; break;
  case EX_ELEM_MAP: inq_code = EX_INQ_EM_PROP; break;
  case EX_NODE_MAP: inq_code = EX_INQ_NM_PROP; break;
  default:
    exerrval = EX_BADPARAM;
    *ierr    = EX_BADPARAM;
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: object type %d not supported; file id %d", *obj_type,
             *idexo);
    ex_err_fn(*idexo, __func__, errmsg, exerrval);
    return;
  }

  /* do Exodus C call to find out how many properties */
  num_props = ex_inquire_int(*idexo, inq_code);
  if (num_props < 0) {
    *ierr = EX_FATAL;
    return;
  }

  /* Allocate space for the name ptr array */
  if (!(aptr = malloc((num_props + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /*
   * Allocate a block of space for the strings, where size = slen,
   * place ptrs into str ptr array,  and Copy Fortran coordinate names
   * to staging space
   */

  if (!(sptr = malloc(num_props * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr);
    return;
  }
  memset(sptr, 0, num_props * (slen + 1));

  for (i = 0; i < num_props; i++) {
    *(aptr + i) = sptr + i * (slen + 1);
  }                /* put ptrs to staging space
                    * into ptr array */
  *(aptr + i) = 0; /* set last pointer to null */

  /* do Exodus C call to get property name records */
  if (ex_get_prop_names(*idexo, (ex_entity_type)*obj_type, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up string space */
    free(aptr); /* free up array ptr space */
    return;
  }
  /* copy C strings to Fortran arrays */
  memset(prop_names, 0, num_props * prop_nameslen);
  for (i = 0; i < num_props; i++) {
    ex_fcdcpy(prop_names + i * prop_nameslen, slen, *(aptr + i)); /* copy and blank fill */
  }

  /* Free up the space we used */
  free(sptr); /* First free up string space */
  free(aptr); /* Then free up array ptr space */
}

/*
 * write object property
 */
void F2C(expp, EXPP)(int *idexo, int *obj_type, entity_id *obj_id, char *prop_name,
                     entity_id *value, int *ierr, int prop_namelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_namelen < slen) {
    slen = prop_namelen;
  }
  /* allocate memory to stage the property name into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy property name from Fortran string to staging area */
  ex_fstrncpy(sptr, prop_name, slen);

  *ierr = ex_put_prop(*idexo, (ex_entity_type)*obj_type, *obj_id, sptr, *value);

  free(sptr);
}

/*
 * read object property
 */
void F2C(exgp, EXGP)(int *idexo, int *obj_type, entity_id *obj_id, char *prop_name, void_int *value,
                     int *ierr, int prop_namelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_namelen < slen) {
    slen = prop_namelen;
  }
  /* allocate memory to stage the property name into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy property name from Fortran string to staging area */
  ex_fstrncpy(sptr, prop_name, slen);

  /* use exodus C routine to get the property value */
  if (ex_get_prop(*idexo, (ex_entity_type)*obj_type, *obj_id, sptr, value) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr);
}

/*
 * read object property array
 */
void F2C(exgpa, EXGPA)(int *idexo, int *obj_type, char *prop_name, void_int *values, int *ierr,
                       int prop_namelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_namelen < slen) {
    slen = prop_namelen;
  }
  /* allocate memory to stage the property name into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  memset(sptr, 0, slen + 1);

  /* Copy property name from Fortran string to staging area */
  ex_fstrncpy(sptr, prop_name, slen);

  /* use exodus C routine to get the values array */
  if (ex_get_prop_array(*idexo, (ex_entity_type)*obj_type, sptr, values) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr);
}

/*
 * write object property array
 */
void F2C(exppa, EXPPA)(int *idexo, int *obj_type, char *prop_name, void_int *values, int *ierr,
                       int prop_namelen)
{
  char *sptr; /* internal string pointer for malloc use */
  int   slen;

  *ierr = 0;

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (prop_namelen < slen) {
    slen = prop_namelen;
  }
  /* allocate memory to stage the property name into */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy property name from Fortran string to staging area */
  ex_fstrncpy(sptr, prop_name, slen);

  /* Use exodus C routine to store the property values */
  if (ex_put_prop_array(*idexo, (ex_entity_type)*obj_type, sptr, values) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr);
}

/*
 * write node set parameters
 */
void F2C(expnp, EXPNP)(int *idexo, entity_id *node_set_id, void_int *num_nodes_in_set,
                       void_int *num_dist_in_set, int *ierr)
{
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int64_t nnis = *(int64_t *)num_nodes_in_set;
    int64_t ndis = *(int64_t *)num_dist_in_set;
    *ierr        = ex_put_set_param(*idexo, EX_NODE_SET, *node_set_id, nnis, ndis);
  }
  else {
    int nnis = *(int *)num_nodes_in_set;
    int ndis = *(int *)num_dist_in_set;
    *ierr    = ex_put_set_param(*idexo, EX_NODE_SET, *node_set_id, nnis, ndis);
  }
}

/*
 * read node set parameters
 */
void F2C(exgnp, EXGNP)(int *idexo, entity_id *node_set_id, void_int *num_nodes_in_set,
                       void_int *num_dist_in_set, int *ierr)
{
  *ierr = ex_get_set_param(*idexo, EX_NODE_SET, *node_set_id, num_nodes_in_set, num_dist_in_set);
}

/*
 * write node set
 */
void F2C(expns, EXPNS)(int *idexo, entity_id *node_set_id, void_int *node_set_node_list, int *ierr)
{
  *ierr = ex_put_set(*idexo, EX_NODE_SET, *node_set_id, node_set_node_list, NULL);
}

/*
 * write node set dist factors
 */
void F2C(expnsd, EXPNSD)(int *idexo, entity_id *node_set_id, real *node_set_dist_fact, int *ierr)
{
  *ierr = ex_put_set_dist_fact(*idexo, EX_NODE_SET, *node_set_id, node_set_dist_fact);
}

/*
 * read node set
 */
void F2C(exgns, EXGNS)(int *idexo, entity_id *node_set_id, void_int *node_set_node_list, int *ierr)
{
  *ierr = ex_get_set(*idexo, EX_NODE_SET, *node_set_id, node_set_node_list, NULL);
}

/*
 * read node set dist factors
 */
void F2C(exgnsd, EXGNSD)(int *idexo, entity_id *node_set_id, real *node_set_dist_fact, int *ierr)
{
  *ierr = ex_get_set_dist_fact(*idexo, EX_NODE_SET, *node_set_id, node_set_dist_fact);
}

/*
 * read node sets IDs
 */
void F2C(exgnsi, EXGNSI)(int *idexo, void_int *node_set_ids, int *ierr)
{
  *ierr = ex_get_ids(*idexo, EX_NODE_SET, node_set_ids);
}

/*
 * write concatenated node sets
 */
void F2C(expcns, EXPCNS)(int *idexo, void_int *node_set_ids, void_int *num_nodes_per_set,
                         void_int *num_dist_per_set, void_int *node_sets_node_index,
                         void_int *node_sets_dist_index, void_int *node_sets_node_list,
                         real *node_sets_dist_fact, int *ierr)
{
  struct ex_set_specs set_specs;
  int                 num_node_sets, i;
  int                 int_size;

  void_int *node_index_ptr, *dist_index_ptr;

  *ierr = 0;

  num_node_sets = ex_inquire_int(*idexo, EX_INQ_NODE_SETS);
  if (num_node_sets < 0) {
    *ierr = EX_FATAL;
    return;
  }

  int_size = sizeof(int);
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int_size = sizeof(int64_t);
  }

  /* allocate memory for C node index array */
  if (!(node_index_ptr = malloc(num_node_sets * int_size))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* allocate memory for C dist factor index array */
  if (!(dist_index_ptr = malloc(num_node_sets * int_size))) {
    free(node_index_ptr);
    *ierr = EX_MEMFAIL;
    return;
  }

  if (int_size == sizeof(int64_t)) {
    for (i = 0; i < num_node_sets; i++) { /* change from 1-based to 0 index */
      ((int64_t *)node_index_ptr)[i] = ((int64_t *)node_sets_node_index)[i] - 1;
      ((int64_t *)dist_index_ptr)[i] = ((int64_t *)node_sets_dist_index)[i] - 1;
    }
  }
  else {
    for (i = 0; i < num_node_sets; i++) { /* change from 1-based to 0 index */
      ((int *)node_index_ptr)[i] = ((int *)node_sets_node_index)[i] - 1;
      ((int *)dist_index_ptr)[i] = ((int *)node_sets_dist_index)[i] - 1;
    }
  }

  set_specs.sets_ids            = node_set_ids;
  set_specs.num_entries_per_set = num_nodes_per_set;
  set_specs.num_dist_per_set    = num_dist_per_set;
  set_specs.sets_entry_index    = node_index_ptr;
  set_specs.sets_dist_index     = dist_index_ptr;
  set_specs.sets_entry_list     = node_sets_node_list;
  set_specs.sets_extra_list     = NULL;
  set_specs.sets_dist_fact      = node_sets_dist_fact;

  *ierr = ex_put_concat_sets(*idexo, EX_NODE_SET, &set_specs);

  free(node_index_ptr);
  free(dist_index_ptr);
}

/*
 * read concatenated node sets
 */
void F2C(exgcns, EXGCNS)(int *idexo, void_int *node_set_ids, void_int *num_nodes_per_set,
                         void_int *num_dist_per_set, void_int *node_sets_node_index,
                         void_int *node_sets_dist_index, void_int *node_sets_node_list,
                         real *node_sets_dist_fact, int *ierr)
{
  struct ex_set_specs set_specs;
  int                 num_node_sets, i;

  set_specs.sets_ids            = node_set_ids;
  set_specs.num_entries_per_set = num_nodes_per_set;
  set_specs.num_dist_per_set    = num_dist_per_set;
  set_specs.sets_entry_index    = node_sets_node_index;
  set_specs.sets_dist_index     = node_sets_dist_index;
  set_specs.sets_entry_list     = node_sets_node_list;
  set_specs.sets_extra_list     = NULL;
  set_specs.sets_dist_fact      = node_sets_dist_fact;

  *ierr = ex_get_concat_sets(*idexo, EX_NODE_SET, &set_specs);

  num_node_sets = ex_inquire_int(*idexo, EX_INQ_NODE_SETS);
  if (num_node_sets < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    for (i = 0; i < num_node_sets; i++) { /* change from 0-based to 1 index */
      ((int64_t *)node_sets_node_index)[i] += 1;
      ((int64_t *)node_sets_dist_index)[i] += 1;
    }
  }
  else {
    for (i = 0; i < num_node_sets; i++) { /* change from 0-based to 1 index */
      ((int *)node_sets_node_index)[i] += 1;
      ((int *)node_sets_dist_index)[i] += 1;
    }
  }
}

/*
 * write side set parameters
 */
void F2C(expsp, EXPSP)(int *idexo, entity_id *side_set_id, void_int *num_sides_in_set,
                       void_int *num_df_in_set, int *ierr)
{
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int64_t nsis = *(int64_t *)num_sides_in_set;
    int64_t ndis = *(int64_t *)num_df_in_set;
    *ierr        = ex_put_set_param(*idexo, EX_SIDE_SET, *side_set_id, nsis, ndis);
  }
  else {
    int nsis = *(int *)num_sides_in_set;
    int ndis = *(int *)num_df_in_set;
    *ierr    = ex_put_set_param(*idexo, EX_SIDE_SET, *side_set_id, nsis, ndis);
  }
}

/*
 * read side set parameters
 */
void F2C(exgsp, EXGSP)(int *idexo, entity_id *side_set_id, void_int *num_sides_in_set,
                       void_int *num_df_in_set, int *ierr)
{
  *ierr = ex_get_set_param(*idexo, EX_SIDE_SET, *side_set_id, num_sides_in_set, num_df_in_set);
}

/*
 * get side set node list length
 */
void F2C(exgsnl, EXGSNL)(int *idexo, entity_id *side_set_id, void_int *num_nodes_in_set, int *ierr)
{
  *ierr = ex_get_side_set_node_list_len(*idexo, *side_set_id, num_nodes_in_set);
}

/*
 * write side set
 */
void F2C(expss, EXPSS)(int *idexo, entity_id *side_set_id, void_int *side_set_elem_list,
                       void_int *side_set_side_list, int *ierr)
{
  *ierr = ex_put_set(*idexo, EX_SIDE_SET, *side_set_id, side_set_elem_list, side_set_side_list);
}

/*
 * read side set
 */
void F2C(exgss, EXGSS)(int *idexo, entity_id *side_set_id, void_int *side_set_elem_list,
                       void_int *side_set_side_list, int *ierr)
{
  *ierr = ex_get_set(*idexo, EX_SIDE_SET, *side_set_id, side_set_elem_list, side_set_side_list);
}

/*
 * write side set distribution factors
 */
void F2C(expssd, EXPSSD)(int *idexo, entity_id *side_set_id, real *side_set_dist_fact, int *ierr)
{
  *ierr = ex_put_set_dist_fact(*idexo, EX_SIDE_SET, *side_set_id, side_set_dist_fact);
}

/*
 * read side set distribution factors
 */
void F2C(exgssd, EXGSSD)(int *idexo, entity_id *side_set_id, real *side_set_dist_fact, int *ierr)
{
  *ierr = ex_get_set_dist_fact(*idexo, EX_SIDE_SET, *side_set_id, side_set_dist_fact);
}

/*
 * read side sets IDs
 */
void F2C(exgssi, EXGSSI)(int *idexo, void_int *side_set_ids, int *ierr)
{
  *ierr = ex_get_ids(*idexo, EX_SIDE_SET, side_set_ids);
}

/*
 * write concatenated side sets
 */
void F2C(expcss, EXPCSS)(int *idexo, void_int *side_set_ids, void_int *num_elem_per_set,
                         void_int *num_dist_per_set, void_int *side_sets_elem_index,
                         void_int *side_sets_dist_index, void_int *side_sets_elem_list,
                         void_int *side_sets_side_list, real *side_sets_dist_fact, int *ierr)
{
  int                 num_side_sets, i;
  void_int *          elem_index_ptr, *dist_index_ptr;
  int                 int_size;
  struct ex_set_specs set_specs;

  *ierr = 0;

  num_side_sets = ex_inquire_int(*idexo, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    *ierr = EX_FATAL;
    return;
  }

  int_size = sizeof(int);
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    int_size = sizeof(int64_t);
  }

  /* allocate memory for C element index array */
  if (!(elem_index_ptr = malloc(num_side_sets * int_size))) {
    *ierr = EX_MEMFAIL;
    return;
  }

  /* allocate memory for C dist factor index array */
  if (!(dist_index_ptr = malloc(num_side_sets * int_size))) {
    free(elem_index_ptr);
    *ierr = EX_MEMFAIL;
    return;
  }

  if (int_size == sizeof(int64_t)) {
    for (i = 0; i < num_side_sets; i++) { /* change from 1-based to 0 index */
      ((int64_t *)elem_index_ptr)[i] = ((int64_t *)side_sets_elem_index)[i] - 1;
      ((int64_t *)dist_index_ptr)[i] = ((int64_t *)side_sets_dist_index)[i] - 1;
    }
  }
  else {
    for (i = 0; i < num_side_sets; i++) { /* change from 1-based to 0 index */
      ((int *)elem_index_ptr)[i] = ((int *)side_sets_elem_index)[i] - 1;
      ((int *)dist_index_ptr)[i] = ((int *)side_sets_dist_index)[i] - 1;
    }
  }

  set_specs.sets_ids            = side_set_ids;
  set_specs.num_entries_per_set = num_elem_per_set;
  set_specs.num_dist_per_set    = num_dist_per_set;
  set_specs.sets_entry_index    = elem_index_ptr;
  set_specs.sets_dist_index     = dist_index_ptr;
  set_specs.sets_entry_list     = side_sets_elem_list;
  set_specs.sets_extra_list     = side_sets_side_list;
  set_specs.sets_dist_fact      = side_sets_dist_fact;

  *ierr = ex_put_concat_sets(*idexo, EX_SIDE_SET, &set_specs);
  free(elem_index_ptr);
  free(dist_index_ptr);
}

/*
 * read concatenated side sets
 */
void F2C(exgcss, EXGCSS)(int *idexo, void_int *side_set_ids, void_int *num_elem_per_set,
                         void_int *num_dist_per_set, void_int *side_sets_elem_index,
                         void_int *side_sets_dist_index, void_int *side_sets_elem_list,
                         void_int *side_sets_side_list, real *side_sets_dist_fact, int *ierr)
{
  int                 i, num_side_sets;
  struct ex_set_specs set_specs;

  *ierr = 0;

  num_side_sets = ex_inquire_int(*idexo, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    *ierr = EX_FATAL;
    return;
  }

  set_specs.sets_ids            = side_set_ids;
  set_specs.num_entries_per_set = num_elem_per_set;
  set_specs.num_dist_per_set    = num_dist_per_set;
  set_specs.sets_entry_index    = side_sets_elem_index;
  set_specs.sets_dist_index     = side_sets_dist_index;
  set_specs.sets_entry_list     = side_sets_elem_list;
  set_specs.sets_extra_list     = side_sets_side_list;
  set_specs.sets_dist_fact      = side_sets_dist_fact;

  *ierr = ex_get_concat_sets(*idexo, EX_SIDE_SET, &set_specs);

  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    for (i = 0; i < num_side_sets; i++) { /* change from 0-based to 1 index */
      ((int64_t *)side_sets_elem_index)[i] += 1;
      ((int64_t *)side_sets_dist_index)[i] += 1;
    }
  }
  else {
    for (i = 0; i < num_side_sets; i++) { /* change from 0-based to 1 index */
      ((int *)side_sets_elem_index)[i] += 1;
      ((int *)side_sets_dist_index)[i] += 1;
    }
  }
}

/*
 * read concatenated side sets (no dist factors)
 */
void F2C(exgcssf, EXGCSSF)(int *idexo, void_int *side_set_ids, void_int *num_elem_per_set,
                           void_int *num_dist_per_set, void_int *side_sets_elem_index,
                           void_int *side_sets_dist_index, void_int *side_sets_elem_list,
                           void_int *side_sets_side_list, int *ierr)
{
  int                 i, num_side_sets;
  struct ex_set_specs set_specs;

  num_side_sets = ex_inquire_int(*idexo, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    *ierr = EX_FATAL;
    return;
  }

  set_specs.sets_ids            = side_set_ids;
  set_specs.num_entries_per_set = num_elem_per_set;
  set_specs.num_dist_per_set    = num_dist_per_set;
  set_specs.sets_entry_index    = side_sets_elem_index;
  set_specs.sets_dist_index     = side_sets_dist_index;
  set_specs.sets_entry_list     = side_sets_elem_list;
  set_specs.sets_extra_list     = side_sets_side_list;
  set_specs.sets_dist_fact      = NULL;

  *ierr = ex_get_concat_sets(*idexo, EX_SIDE_SET, &set_specs);

  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    for (i = 0; i < num_side_sets; i++) { /* change from 0-based to 1 index */
      ((int64_t *)side_sets_elem_index)[i] += 1;
      ((int64_t *)side_sets_dist_index)[i] += 1;
    }
  }
  else {
    for (i = 0; i < num_side_sets; i++) { /* change from 0-based to 1 index */
      ((int *)side_sets_elem_index)[i] += 1;
      ((int *)side_sets_dist_index)[i] += 1;
    }
  }
}

/*
 * write results variables parameters
 */
void F2C(expvp, EXPVP)(int *idexo, char *var_type, int *num_vars, int *ierr, int var_typelen)
{
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);
  *ierr                   = ex_put_variable_param(*idexo, obj_type, *num_vars);
  EX_UNUSED(var_typelen);
}

/*
 * read results variables parameters
 */
void F2C(exgvp, EXGVP)(int *idexo, char *var_type, int *num_vars, int *ierr, int var_typelen)
{
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);
  *ierr                   = ex_get_variable_param(*idexo, obj_type, num_vars);
  EX_UNUSED(var_typelen);
}

/*
 * write results variables names
 */
void F2C(expvan, EXPVAN)(int *idexo, char *var_type, int *num_vars, char *var_names, int *ierr,
                         int var_typelen, int var_nameslen)
{
  char **aptr;         /* ptr to temp staging space for string array
                        * ptrs */
  char *         sptr; /* ptr to temp staging space for strings */
  int            i, slen;
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);

  *ierr = 0; /* default no error */
  EX_UNUSED(var_typelen);

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (var_nameslen < slen) {
    slen = var_nameslen;
  }
  /* allocate memory for pointer array */
  if (!(aptr = malloc((*num_vars + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_vars * (slen + 1) * sizeof(char)))) {
    free(aptr); /* Free up string ptr array */
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Copy Fortran variable names to staging space */
  for (i = 0; i < *num_vars; i++) {
    *(aptr + i) = sptr + i * (slen + 1);                          /* put address into ptr array */
    ex_fstrncpy(*(aptr + i), var_names + i * var_nameslen, slen); /* copy string into
                                                                   * buffer */
  }
  *(aptr + i) = 0; /* null out last ptr */
  /* do Exodus C call to write results variables names */
  if (ex_put_variable_names(*idexo, obj_type, *num_vars, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
  }
  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}
/*
 * read results variables names
 */
void F2C(exgvan, EXGVAN)(int *idexo, char *var_type, int *num_vars, char *var_names, int *ierr,
                         int var_typelen, int var_nameslen)
{
  char **aptr;         /* ptr to temp staging space for string array
                        * ptrs */
  char *         sptr; /* ptr to temp staging space for strings */
  int            i, slen;
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);
  EX_UNUSED(var_typelen);

  *ierr = 0; /* default no error */

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (var_nameslen < slen) {
    slen = var_nameslen;
  }
  /* allocate memory to for pointer array */
  if (!(aptr = malloc((*num_vars + 1) * sizeof(char *)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* Allocate staging space for the variable names */
  if (!(sptr = malloc(*num_vars * (slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    free(aptr); /* Free up string ptr array */
    return;
  }
  for (i = 0; i < *num_vars; i++) {
    *(aptr + i) = sptr + i * (slen + 1); /* put address into ptr array */
  }
  *(aptr + i) = 0; /* null out last ptr */

  /* do Exodus C call to read results variables names */
  if (ex_get_variable_names(*idexo, obj_type, *num_vars, aptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    free(aptr);
    return;
  }
  /* Copy Fortran variable names to staging space */
  memset(var_names, 0, *num_vars * var_nameslen);
  for (i = 0; i < *num_vars; i++) {
    ex_fcdcpy(var_names + i * var_nameslen, slen, *(aptr + i)); /* copy str into Fortran
                                                                 * buffer */
  }

  free(sptr); /* Free up string staging area */
  free(aptr); /* Free up string ptr array */
}

/*
 * write element variable truth table
 */
void F2C(expvtt, EXPVTT)(int *idexo, int *num_elem_blk, int *num_elem_var, int *elem_var_tab,
                         int *ierr)
{
  *ierr = ex_put_truth_table(*idexo, EX_ELEM_BLOCK, *num_elem_blk, *num_elem_var, elem_var_tab);
}

/*
 * write nodeset variable truth table
 */
void F2C(expnstt, EXPNSTT)(int *idexo, int *num_entity, int *num_var, int *var_tab, int *ierr)
{
  *ierr = ex_put_truth_table(*idexo, EX_NODE_SET, *num_entity, *num_var, var_tab);
}

/*
 * write sideset variable truth table
 */
void F2C(expsstt, EXPSSTT)(int *idexo, int *num_entity, int *num_var, int *var_tab, int *ierr)
{
  *ierr = ex_put_truth_table(*idexo, EX_SIDE_SET, *num_entity, *num_var, var_tab);
}

/*
 * read element variable truth table
 */
void F2C(exgvtt, EXGVTT)(int *idexo, int *num_elem_blk, int *num_elem_var, int *elem_var_tab,
                         int *ierr)
{
  *ierr = ex_get_truth_table(*idexo, EX_ELEM_BLOCK, *num_elem_blk, *num_elem_var, elem_var_tab);
}

/*
 * read nodeset variable truth table
 */
void F2C(exgnstt, EXGNSTT)(int *idexo, int *num_entity, int *num_var, int *var_tab, int *ierr)
{
  *ierr = ex_get_truth_table(*idexo, EX_NODE_SET, *num_entity, *num_var, var_tab);
}

/*
 * read sideset variable truth table
 */
void F2C(exgsstt, EXGSSTT)(int *idexo, int *num_entity, int *num_var, int *var_tab, int *ierr)
{
  *ierr = ex_get_truth_table(*idexo, EX_SIDE_SET, *num_entity, *num_var, var_tab);
}

/*
 * write global variable values at time step
 */
void F2C(expgv, EXPGV)(int *idexo, int *time_step, int *num_glob_vars, real *glob_var_vals,
                       int *ierr)
{
  *ierr = ex_put_var(*idexo, *time_step, EX_GLOBAL, 1, 0, *num_glob_vars, glob_var_vals);
}

/*
 * read global variable values at a time step
 */
void F2C(exggv, EXGGV)(int *idexo, int *time_step, int *num_glob_vars, real *glob_var_vals,
                       int *ierr)
{
  *ierr = ex_get_var(*idexo, *time_step, EX_GLOBAL, 1, 0, *num_glob_vars, glob_var_vals);
}

/*
 * read global variable values through time
 */
void F2C(exggvt, EXGGVT)(int *idexo, int *glob_var_index, int *beg_time_step, int *end_time_step,
                         real *glob_var_vals, int *ierr)
{
  *ierr = ex_get_var_time(*idexo, EX_GLOBAL, *glob_var_index, 1, *beg_time_step, *end_time_step,
                          glob_var_vals);
}

/*
 * write nodal variable values at a time step
 */
void F2C(expnv, EXPNV)(int *idexo, int *time_step, int *nodal_var_index, void_int *num_nodes,
                       real *nodal_var_vals, int *ierr)
{
  int64_t nnodes;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    nnodes = *(int64_t *)num_nodes;
  }
  else {
    nnodes = *(int *)num_nodes;
  }

  *ierr = ex_put_var(*idexo, *time_step, EX_NODAL, *nodal_var_index, 1, nnodes, nodal_var_vals);
}

/*
 * read nodal variable values at a time step
 */
void F2C(exgnv, EXGNV)(int *idexo, int *time_step, int *nodal_var_index, void_int *num_nodes,
                       real *nodal_var_vals, int *ierr)
{
  int64_t nnodes;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    nnodes = *(int64_t *)num_nodes;
  }
  else {
    nnodes = *(int *)num_nodes;
  }

  *ierr = ex_get_var(*idexo, *time_step, EX_NODAL, *nodal_var_index, 1, nnodes, nodal_var_vals);
}

/*
 * read nodal variable values through time
 */
void F2C(exgnvt, EXGNVT)(int *idexo, int *nodal_var_index, void_int *node_number,
                         int *beg_time_step, int *end_time_step, real *nodal_var_vals, int *ierr)
{
  int64_t nnode;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    nnode = *(int64_t *)node_number;
  }
  else {
    nnode = *(int *)node_number;
  }

  *ierr = ex_get_var_time(*idexo, EX_NODAL, *nodal_var_index, nnode, *beg_time_step, *end_time_step,
                          nodal_var_vals);
}

/*
 * write element variable values at a time step
 */
void F2C(expev, EXPEV)(int *idexo, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
                       void_int *num_elem_this_blk, real *elem_var_vals, int *ierr)
{
  int64_t neblk;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    neblk = *(int64_t *)num_elem_this_blk;
  }
  else {
    neblk = *(int *)num_elem_this_blk;
  }

  *ierr = ex_put_var(*idexo, *time_step, EX_ELEM_BLOCK, *elem_var_index, *elem_blk_id, neblk,
                     elem_var_vals);
}

/*
 * read element variable values at a time step
 */
void F2C(exgev, EXGEV)(int *idexo, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
                       void_int *num_elem_this_blk, real *elem_var_vals, int *ierr)
{
  int64_t neblk;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    neblk = *(int64_t *)num_elem_this_blk;
  }
  else {
    neblk = *(int *)num_elem_this_blk;
  }

  *ierr = ex_get_var(*idexo, *time_step, EX_ELEM_BLOCK, *elem_var_index, *elem_blk_id, neblk,
                     elem_var_vals);
}

/*
 * read element variable values through time
 */
void F2C(exgevt, EXGEVT)(int *idexo, int *elem_var_index, void_int *elem_number, int *beg_time_step,
                         int *end_time_step, real *elem_var_vals, int *ierr)
{
  int64_t el_num;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    el_num = *(int64_t *)elem_number;
  }
  else {
    el_num = *(int *)elem_number;
  }

  *ierr = ex_get_var_time(*idexo, EX_ELEM_BLOCK, *elem_var_index, el_num, *beg_time_step,
                          *end_time_step, elem_var_vals);
}

/*
 * write nodeset variable values at a time step
 */
void F2C(expnsv, EXPNSV)(int *idexo, int *time_step, int *var_index, entity_id *id,
                         void_int *num_entity, real *var_vals, int *ierr)
{
  int64_t n_entity;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    n_entity = *(int64_t *)num_entity;
  }
  else {
    n_entity = *(int *)num_entity;
  }

  *ierr = ex_put_var(*idexo, *time_step, EX_NODE_SET, *var_index, *id, n_entity, var_vals);
}

/*
 * read nodeset variable values at a time step
 */
void F2C(exgnsv, EXGNSV)(int *idexo, int *time_step, int *var_index, entity_id *id,
                         void_int *num_entity, real *var_vals, int *ierr)
{
  int64_t n_entity;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    n_entity = *(int64_t *)num_entity;
  }
  else {
    n_entity = *(int *)num_entity;
  }

  *ierr = ex_get_var(*idexo, *time_step, EX_NODE_SET, *var_index, *id, n_entity, var_vals);
}

/*
 * write sideset variable values at a time step
 */
void F2C(expssv, EXPSSV)(int *idexo, int *time_step, int *var_index, entity_id *id,
                         void_int *num_entity, real *var_vals, int *ierr)
{
  int64_t n_entity;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    n_entity = *(int64_t *)num_entity;
  }
  else {
    n_entity = *(int *)num_entity;
  }

  *ierr = ex_put_var(*idexo, *time_step, EX_SIDE_SET, *var_index, *id, n_entity, var_vals);
}

/*
 * read sideset variable values at a time step
 */
void F2C(exgssv, EXGSSV)(int *idexo, int *time_step, int *var_index, entity_id *id,
                         void_int *num_entity, real *var_vals, int *ierr)
{
  int64_t n_entity;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    n_entity = *(int64_t *)num_entity;
  }
  else {
    n_entity = *(int *)num_entity;
  }

  *ierr = ex_get_var(*idexo, *time_step, EX_SIDE_SET, *var_index, *id, n_entity, var_vals);
}

/*
 * write time value for a time step
 */
void F2C(exptim, EXPTIM)(int *idexo, int *time_step, real *time_value, int *ierr)
{
  *ierr = ex_put_time(*idexo, *time_step, time_value);
}

/*
 * read time value for a time step
 */
void F2C(exgtim, EXGTIM)(int *idexo, int *time_step, real *time_value, int *ierr)
{
  *ierr = ex_get_time(*idexo, *time_step, time_value);
}

/*
 * read all time values
 */
void F2C(exgatm, EXGATM)(int *idexo, real *time_values, int *ierr)
{
  *ierr = ex_get_all_times(*idexo, time_values);
}

/*
 * inquire EXODUS parameters
 */
void F2C(exinq, EXINQ)(int *idexo, int *req_info, void_int *ret_int, float *ret_float,
                       char *ret_char, int *ierr, int ret_charlen)
{
  *ierr = ex_inquire(*idexo, (ex_inquiry)*req_info, ret_int, ret_float, ret_char);
  EX_UNUSED(ret_charlen);
}

/*
 * inquire integer EXODUS parameters
 */
int64_t F2C(exinqi, EXINQI)(int *idexo, int *req_info)
{
  return ex_inquire_int(*idexo, (ex_inquiry)*req_info);
}

/*
 * convert side set node lists to side set side lists
 */
void F2C(excn2s, EXCN2S)(int *idexo, void_int *num_elem_per_set, void_int *num_nodes_per_set,
                         void_int *side_sets_elem_index, void_int *side_sets_node_index,
                         void_int *side_sets_elem_list, void_int *side_sets_node_list,
                         void_int *side_sets_side_list, int *ierr)
{
  EX_UNUSED(side_sets_elem_index);
  EX_UNUSED(side_sets_node_index);

  *ierr = ex_cvt_nodes_to_sides(*idexo, num_elem_per_set, num_nodes_per_set, NULL, /* unused */
                                NULL,                                              /* unused */
                                side_sets_elem_list, side_sets_node_list, side_sets_side_list);
}

/*
 * read side set node list
 */
void F2C(exgssn, EXGSSN)(int *idexo, entity_id *side_set_id, int *side_set_node_cnt_list,
                         void_int *side_set_node_list, int *ierr)
{
  *ierr =
      ex_get_side_set_node_list(*idexo, *side_set_id, side_set_node_cnt_list, side_set_node_list);
}

/*
 * read side set node count
 */
void F2C(exgssc, EXGSSC)(int *idexo, entity_id *side_set_id, int *side_set_node_cnt_list, int *ierr)
{
  *ierr = ex_get_side_set_node_count(*idexo, *side_set_id, side_set_node_cnt_list);
}

/*
 * read concatenated side set node count
 */
void F2C(exgcssc, EXGCSSC)(int *idexo, int *side_set_node_cnt_list, int *ierr)
{
  *ierr = ex_get_concat_side_set_node_count(*idexo, side_set_node_cnt_list);
}

/* ex_get_coordinate_frames -- read coordinate frames */
void F2C(exgfrm, EXGFRM)(int *idexo, int *nframeo, void_int *cfids, real *coord, int *tags,
                         int *ierr)
{
  int   i;
  char *ctags = NULL;

  /* Determine number of coordinate frames stored in file */
  int nframe = ex_inquire_int(*idexo, EX_INQ_COORD_FRAMES);
  if (nframe < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (nframe != *nframeo) {
    *ierr = EX_FATAL;
    return;
  }
  /* Create array of characters to store tags... */
  if (nframe > 0) {
    if (!(ctags = calloc(nframe, sizeof(char)))) {
      *ierr = EX_MEMFAIL;
      return;
    }
    *ierr = 0;

    if (ex_get_coordinate_frames(*idexo, &nframe, cfids, coord, ctags) == EX_FATAL) {
      *ierr = EX_FATAL;
      return;
    }
    /* Convert character tags back to integer tags for fortran */
    for (i = 0; i < nframe; i++) {
      if (ctags[i] == 'R' || ctags[i] == 'r') {
        tags[i] = EX_CF_RECTANGULAR;
      }
      else if (ctags[i] == 'C' || ctags[i] == 'c') {
        tags[i] = EX_CF_CYLINDRICAL;
      }
      else if (ctags[i] == 'S' || ctags[i] == 's') {
        tags[i] = EX_CF_SPHERICAL;
      }
    }
    free(ctags);
  }
}

/* ex_put_coordinate_frames -- define/write coordinate frames */
void F2C(expfrm, EXPFRM)(int *idexo, int *nframe, void_int *cfids, real *coord, int *tags,
                         int *ierr)
{
  int   i;
  char *ctags = NULL;

  /* Create array of characters to store tags... */
  if (*nframe > 0) {
    if (!(ctags = calloc(*nframe, sizeof(char)))) {
      *ierr = EX_MEMFAIL;
      return;
    }
    /* Convert fortran integer tags to C API character tags */
    for (i = 0; i < *nframe; i++) {
      if (tags[i] == EX_CF_RECTANGULAR) {
        ctags[i] = 'R';
      }
      else if (tags[i] == EX_CF_CYLINDRICAL) {
        ctags[i] = 'C';
      }
      else if (tags[i] == EX_CF_SPHERICAL) {
        ctags[i] = 'S';
      }
    }

    *ierr = 0;

    if (ex_put_coordinate_frames(*idexo, *nframe, cfids, coord, ctags) == EX_FATAL) {
      *ierr = EX_FATAL;
    }
    free(ctags);
  }
}

/* Routine to return floating point word size */
int F2C(excpws, EXCPWS)() { return (ex__get_cpu_ws()); }

/* Routine to return large model setting */
int F2C(exlgmd, EXLGMD)(int *idexo) { return (ex_large_model(*idexo)); }

/* Generalized error handling function */
void F2C(exerr, EXERR)(char *pname, char *err_string, int *errcode, int pnamelen, int err_stringlen)
{

  char *proc_name, *error_string;
  if (!(proc_name = malloc((pnamelen + 1) * sizeof(char)))) {
    ex_err(__func__, "Error: failed to allocate space for process name buffer", EX_MEMFAIL);
    *errcode = EX_MEMFAIL;
    return;
  }
  if (!(error_string = malloc((err_stringlen + 1) * sizeof(char)))) {
    free(proc_name);
    ex_err(__func__, "Error: failed to allocate space for error msg buffer", EX_MEMFAIL);
    *errcode = EX_MEMFAIL;
    return;
  }
  ex_fstrncpy(proc_name, pname, pnamelen);
  ex_fstrncpy(error_string, err_string, err_stringlen);
  ex_err(proc_name, error_string, *errcode);
  free(proc_name);
  free(error_string);
}

/* Error message reporting options setting function */
void F2C(exopts, EXOPTS)(int *option_val, int *ierr)
{
  *ierr = 0;
  ex_opts((ex_options)*option_val);
  if (exerrval != 0) {
    *ierr = EX_FATAL;
  }
}

void F2C(exsi64, EXSI64)(int *idexo, int *mode, int *ierr)
{
  *ierr = ex_set_int64_status(*idexo, *mode);
}

int64_t F2C(exgi64, EXGI64)(int *idexo) { return ex_int64_status(*idexo); }

void F2C(exmxnm, EXMXNM)(int *idexo, int *length, int *ierr)
{
  *ierr = ex_set_max_name_length(*idexo, *length);
}

void F2C(exsetopt, EXSETOPT)(int *idexo, int *option, int *value, int *ierr)
{
  *ierr = ex_set_option(*idexo, *option, *value);
}

/*
 * copy EXODUS file
 */
void F2C(excopy, EXCOPY)(int *idexo_in, int *idexo_out, int *ierr)
{
  *ierr = ex_copy(*idexo_in, *idexo_out);
}

/*
 * get element map
 */
void

    F2C(exgem, EXGEM)(int *idexo, entity_id *map_id, void_int *elem_map, int *ierr)
{
  *ierr = ex_get_num_map(*idexo, EX_ELEM_MAP, *map_id, elem_map);
}

/*
 * get partial_element map
 */
void F2C(exgpem, EXGPEM)(int *idexo, entity_id *map_id, void_int *start, void_int *count,
                         void_int *elem_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }

  *ierr = ex_get_partial_num_map(*idexo, EX_ELEM_MAP, *map_id, st, cnt, elem_map);
}

/*
 * get element number map
 */
void F2C(exgenm, EXGENM)(int *idexo, void_int *elem_map, int *ierr)
{
  *ierr = ex_get_id_map(*idexo, EX_ELEM_MAP, elem_map);
}

/*
 * get map parameters
 */
void F2C(exgmp, EXGMP)(int *idexo, int *num_node_maps, int *num_elem_maps, int *ierr)
{
  *ierr = ex_get_map_param(*idexo, num_node_maps, num_elem_maps);
}

/*
 * get node map
 */
void F2C(exgnm, EXGNM)(int *idexo, entity_id *map_id, void_int *node_map, int *ierr)
{
  *ierr = ex_get_num_map(*idexo, EX_NODE_MAP, *map_id, node_map);
}

/*
 * get node number map
 */
void F2C(exgnnm, EXGNNM)(int *idexo, void_int *node_map, int *ierr)
{
  *ierr = ex_get_id_map(*idexo, EX_NODE_MAP, node_map);
}

/*
 * read single results variables name
 */
void F2C(exgvnm, EXGVNM)(int *idexo, char *var_type, int *var_index, char *var_name, int *ierr,
                         int var_typelen, int var_namelen)
{
  char *sptr; /* ptr to temp staging space for string */
  int   slen;
  *ierr                   = 0; /* default no error */
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);

  slen = ex_inquire_int(*idexo, EX_INQ_MAX_READ_NAME_LENGTH); /* max string size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (var_namelen < slen) {
    slen = var_namelen;
  }
  /* Allocate staging space for the variable name */
  if (!(sptr = malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  /* do Exodus C call to read results variables names */
  if (ex_get_variable_name(*idexo, obj_type, *var_index, sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    return;
  }
  /* Copy Fortran variable names to staging space */
  /** printf("[exgvnm] var_name(%d): %s\n",*var_index,sptr)); **/
  memset(var_name, 0, var_namelen);
  ex_fcdcpy(var_name, slen, sptr); /* copy string into Fortran buffer */

  free(sptr); /* Free up string staging area */
  EX_UNUSED(var_typelen);
}

/*
 * put element map
 */
void F2C(expem, EXPEM)(int *idexo, entity_id *map_id, void_int *elem_map, int *ierr)
{
  *ierr = ex_put_num_map(*idexo, EX_ELEM_MAP, *map_id, elem_map);
}

/*
 * put partial element map
 */
void F2C(exppem, EXPPEM)(int *idexo, entity_id *map_id, void_int *start, void_int *count,
                         void_int *elem_map, int *ierr)
{
  int64_t st, cnt;
  if (ex_int64_status(*idexo) & EX_BULK_INT64_API) {
    st  = *(int64_t *)start;
    cnt = *(int64_t *)count;
  }
  else {
    st  = *(int *)start;
    cnt = *(int *)count;
  }
  *ierr = ex_put_partial_num_map(*idexo, EX_ELEM_MAP, *map_id, st, cnt, elem_map);
}

/*
 * put element number map
 */
void F2C(expenm, EXPENM)(int *idexo, void_int *elem_map, int *ierr)
{
  *ierr = ex_put_id_map(*idexo, EX_ELEM_MAP, elem_map);
}

/*
 * put map parameters
 */
void F2C(expmp, EXPMP)(int *idexo, int *num_node_maps, int *num_elem_maps, int *ierr)
{
  *ierr = ex_put_map_param(*idexo, *num_node_maps, *num_elem_maps);
}

/*
 * put node map
 */
void F2C(expnm, EXPNM)(int *idexo, entity_id *map_id, void_int *node_map, int *ierr)
{
  *ierr = ex_put_num_map(*idexo, EX_NODE_MAP, *map_id, node_map);
}

/*
 * put node number map
 */
void F2C(expnnm, EXPNNM)(int *idexo, void_int *node_map, int *ierr)
{
  *ierr = ex_put_id_map(*idexo, EX_NODE_MAP, node_map);
}

/*
 * write results variable name
 */
void F2C(expvnm, EXPVNM)(int *idexo, char *var_type, int *var_index, char *var_name, int *ierr,
                         int var_typelen, int var_namelen)
{
  char *sptr; /* ptr to temp staging space for string */
  int   slen;
  *ierr                   = 0; /* default no error */
  ex_entity_type obj_type = ex_var_type_to_ex_entity_type(*var_type);

  slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH); /* max str size */
  if (slen < 0) {
    *ierr = EX_FATAL;
    return;
  }

  if (var_namelen < slen) {
    slen = var_namelen;
  }
  /* Allocate staging space for the variable name */
  if (!(sptr = (char *)malloc((slen + 1) * sizeof(char)))) {
    *ierr = EX_MEMFAIL;
    return;
  }
  ex_fstrncpy(sptr, var_name, slen); /* copy string into buffer */

  /* do Exodus C call to write results variable name */
  if (ex_put_variable_name(*idexo, obj_type, *var_index, sptr) == EX_FATAL) {
    *ierr = EX_FATAL;
    free(sptr); /* free up allocated space */
    return;
  }
  free(sptr); /* Free up string staging area */
  EX_UNUSED(var_typelen);
}

/*
 *  Get initial information from nemesis file
 */
void F2C(exgii, EXGII)(int *idne, int *nproc, int *nproc_in_f, char *ftype, int *ierr,
                       size_t ftypelen)
{
  size_t slen = 1;
  char * file_type;

  /* WARNING: ftypelen SHOULD be 1, but may not be depending on how
              the Fortran programmer passed it. It is best at
              this time to hard code it per NEPII spec. */
  if (ftypelen != 1) {
#if defined(EXODUS_STRING_LENGTH_WARNING)
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: file type string length is %lu in file id %d\n",
             ftypelen, *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
#endif
    slen = ftypelen;
  }

  file_type = (char *)malloc((slen + 1) * sizeof(char));

  if ((*ierr = ex_get_init_info(*idne, nproc, nproc_in_f, file_type)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to get initial information from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }

  if (*ierr == 0) {
    ex_fcdcpy(ftype, slen, file_type);
  }

  free(file_type);
}

/*
 *  Write initial information from nemesis file
 */
void F2C(expii, EXPII)(int *idne, int *nproc, int *nproc_in_f, char *ftype, int *ierr,
                       size_t ftypelen)
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
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: file type string length is %lu in file id %d\n",
             ftypelen, *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
#endif
  }

  file_type = (char *)malloc((slen + 1) * sizeof(char));

  ex_fstrncpy(file_type, ftype, slen);

  if ((*ierr = ex_put_init_info(*idne, *nproc, *nproc_in_f, file_type)) != 0) {
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to put initial information in file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }

  free(file_type);
}

/*
 * Read initial global information
 */
void F2C(exgig, EXGIG)(int *idne, void_int *nnodes_g, void_int *nelems_g, void_int *nelem_blks_g,
                       void_int *nnode_sets_g, void_int *nside_sets_g, int *ierr)
{
  if ((*ierr = ex_get_init_global(*idne, nnodes_g, nelems_g, nelem_blks_g, nnode_sets_g,
                                  nside_sets_g)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read initial global information from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write initial global information
 */
void F2C(expig, EXPIG)(int *idne, void_int *nnodes_g, void_int *nelems_g, void_int *nelem_blks_g,
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
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to store initial global information in file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read load balance parameters
 */
void F2C(exglbp, EXGLBP)(int *idne, void_int *nint_nodes, void_int *nbor_nodes,
                         void_int *next_nodes, void_int *nint_elems, void_int *nbor_elems,
                         void_int *nnode_cmaps, void_int *nelem_cmaps, int *processor, int *ierr)
{
  if ((*ierr = ex_get_loadbal_param(*idne, nint_nodes, nbor_nodes, next_nodes, nint_elems,
                                    nbor_elems, nnode_cmaps, nelem_cmaps, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read load balance parameters from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write load balance parameters
 */
void F2C(explbp, EXPLBP)(int *idne, void_int *nint_nodes, void_int *nbor_nodes,
                         void_int *next_nodes, void_int *nint_elems, void_int *nbor_elems,
                         void_int *nnode_cmaps, void_int *nelem_cmaps, int *processor, int *ierr)
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
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to store load balance parameters in file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write concatenated load balance parameters
 */
void F2C(explbpc, EXPLBPC)(int *idne, void_int *nint_nodes, void_int *nbor_nodes,
                           void_int *next_nodes, void_int *nint_elems, void_int *nbor_elems,
                           void_int *nnode_cmaps, void_int *nelem_cmaps, int *ierr)
{
  if ((*ierr = ex_put_loadbal_param_cc(*idne, nint_nodes, nbor_nodes, next_nodes, nint_elems,
                                       nbor_elems, nnode_cmaps, nelem_cmaps)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to store load balance parameters in file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read global node set parameters
 */
void F2C(exgnspg, EXGNSPG)(int *idne, void_int *ns_ids_glob, void_int *ns_n_cnt_glob,
                           void_int *ns_df_cnt_glob, int *ierr)
{
  if ((*ierr = ex_get_ns_param_global(*idne, ns_ids_glob, ns_n_cnt_glob, ns_df_cnt_glob)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read global node set parameters from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write global node set parameters
 */
void F2C(expnspg, EXPNSPG)(int *idne, void_int *global_ids, void_int *global_n_cnts,
                           void_int *global_df_cnts, int *ierr)
{
  if ((*ierr = ex_put_ns_param_global(*idne, global_ids, global_n_cnts, global_df_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to store global node set parameters in file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read global side set parameters
 */
void F2C(exgsspg, EXGSSPG)(int *idne, void_int *ss_ids_glob, void_int *ss_n_cnt_glob,
                           void_int *ss_df_cnt_glob, int *ierr)
{

  if ((*ierr = ex_get_ss_param_global(*idne, ss_ids_glob, ss_n_cnt_glob, ss_df_cnt_glob)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read global side set parameters from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write global side set parameters
 */
void F2C(expsspg, EXPSSPG)(int *idne, void_int *global_ids, void_int *global_el_cnts,
                           void_int *global_df_cnts, int *ierr)
{
  if ((*ierr = ex_put_ss_param_global(*idne, global_ids, global_el_cnts, global_df_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to store global side set parameters in file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read global element block information
 */
void F2C(exgebig, EXGEBIG)(int *idne, void_int *el_blk_ids, void_int *el_blk_cnts, int *ierr)
{
  if ((*ierr = ex_get_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read global element block info from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write global element block information
 */
void F2C(expebig, EXPEBIG)(int *idne, void_int *el_blk_ids, void_int *el_blk_cnts, int *ierr)
{
  if ((*ierr = ex_put_eb_info_global(*idne, el_blk_ids, el_blk_cnts)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to store global element block info in file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read side set element list and side set side list
 */
void F2C(exgnss, EXGNSS)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_get_partial_set(*idne, EX_SIDE_SET, *ss_id, st, cnt, ss_elem_list,
                                  ss_side_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read side set element list from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write side set element list and side set side list
 */
void F2C(expnss, EXPNSS)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_put_partial_set(*idne, EX_SIDE_SET, *ss_id, st, cnt, ss_elem_list,
                                  ss_side_list)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write side set element list to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read side set distribution factor
 */
void F2C(exgnssd, EXGNSSD)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
                           real *ss_df, int *ierr)
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

  if ((*ierr = ex_get_partial_set_dist_fact(*idne, EX_SIDE_SET, *ss_id, st, cnt, ss_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read side set dist factor from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write side set distribution factor
 */
void F2C(expnssd, EXPNSSD)(int *idne, entity_id *ss_id, void_int *start, void_int *count,
                           real *ss_df, int *ierr)
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

  if ((*ierr = ex_put_partial_set_dist_fact(*idne, EX_SIDE_SET, *ss_id, st, cnt, ss_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write side set dist factor to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read node set list for a single node set
 */
void F2C(exgnns, EXGNNS)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_get_partial_set(*idne, EX_NODE_SET, *ns_id, st, cnt, ns_node_list, NULL)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read node set node list from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write node set list for a single node set
 */
void F2C(expnns, EXPNNS)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_put_partial_set(*idne, EX_NODE_SET, *ns_id, st, cnt, ns_node_list, NULL)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write node set node list to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read node set distribution factor
 */
void F2C(exgnnsd, EXGNNSD)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
                           real *ns_df, int *ierr)
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

  if ((*ierr = ex_get_partial_set_dist_fact(*idne, EX_NODE_SET, *ns_id, st, cnt, ns_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read node set dist factor from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write node set distribution factor
 */
void F2C(expnnsd, EXPNNSD)(int *idne, entity_id *ns_id, void_int *start, void_int *count,
                           real *ns_df, int *ierr)
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

  if ((*ierr = ex_put_partial_set_dist_fact(*idne, EX_NODE_SET, *ns_id, st, cnt, ns_df)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write node set dist factor to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read coordinates of the nodes
 */
void F2C(exgncor, EXGNCOR)(int *idne, void_int *start, void_int *count, real *x_coor, real *y_coor,
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

  if ((*ierr = ex_get_partial_coord(*idne, st, cnt, x_coor, y_coor, z_coor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read node coordinates from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write coordinates of the nodes
 */
void F2C(expncor, EXPNCOR)(int *idne, void_int *start, void_int *count, real *x_coor, real *y_coor,
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

  if ((*ierr = ex_put_partial_coord(*idne, st, cnt, x_coor, y_coor, z_coor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write node coordinates to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read an element block's connectivity list
 */
void F2C(exgnec, EXGNEC)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_get_partial_conn(*idne, EX_ELEM_BLOCK, *elem_blk_id, st, cnt, connect, NULL,
                                   NULL)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read element block connectivity from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write an element block's connectivity list
 */
void F2C(expnec, EXPNEC)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
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

  if ((*ierr = ex_put_partial_conn(*idne, EX_ELEM_BLOCK, *elem_blk_id, st, cnt, connect, NULL,
                                   NULL)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to write element block connectivity to file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read an element block's attributes
 */
void F2C(exgneat, EXGNEAT)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
                           real *attrib, int *ierr)
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

  if ((*ierr = ex_get_partial_attr(*idne, EX_ELEM_BLOCK, *elem_blk_id, st, cnt, attrib)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read element block attribute from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write an element block's attributes
 */
void F2C(expneat, EXPNEAT)(int *idne, entity_id *elem_blk_id, void_int *start, void_int *count,
                           real *attrib, int *ierr)
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

  if ((*ierr = ex_put_partial_attr(*idne, EX_ELEM_BLOCK, *elem_blk_id, st, cnt, attrib)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write element block attribute to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element type for a specific element block
 */
void F2C(exgelt, EXGELT)(int *idne, entity_id *elem_blk_id, char *elem_type, int *ierr,
                         size_t elem_typelen)
{
  size_t slen = MAX_STR_LENGTH;
  char * etype;

  /* WARNING: ftypelen SHOULD be MAX_STR_LENGTH, but may not be depending
              on how the Fortran programmer passed it. It is best at
              this time to hard code it per NEMESIS spec. */
  if (elem_typelen != MAX_STR_LENGTH) {
#if defined(EXODUS_STRING_LENGTH_WARNING)
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: element type string length is %lu in file id %d\n",
             elem_typelen, *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
#endif
    slen = elem_typelen;
  }

  etype = (char *)malloc((slen + 1) * sizeof(char));

  if ((*ierr = ex_get_elem_type(*idne, *elem_blk_id, etype)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read element block type from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }

  if (*ierr == 0) {
    ex_fcdcpy(elem_type, slen, etype);
  }

  free(etype);
}

/*
 * Read a variable for an element block
 */
void F2C(exgnev, EXGNEV)(int *idne, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
                         void_int *num_elem_this_blk, void_int *start, void_int *count,
                         real *elem_var_vals, int *ierr)
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

  if ((*ierr = ex_get_partial_var(*idne, *time_step, EX_ELEM_BLOCK, *elem_var_index, *elem_blk_id,
                                  st, cnt, elem_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read element block variable from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
  EX_UNUSED(num_elem_this_blk);
}

/*
 * Write a variable slab for an element block
 */
void F2C(expevs, EXPEVS)(int *idne, int *time_step, int *elem_var_index, entity_id *elem_blk_id,
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

  if ((*ierr = ex_put_partial_var(*idne, *time_step, EX_ELEM_BLOCK, *elem_var_index, *elem_blk_id,
                                  st, cnt, elem_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to write elem block variable slab to file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the values of a single nodal variable for a single time step
 */
void F2C(exgnnv, EXGNNV)(int *idne, int *time_step, int *nodal_var_index, void_int *start,
                         void_int *count, real *nodal_vars, int *ierr)
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

  if ((*ierr = ex_get_partial_var(*idne, *time_step, EX_NODAL, *nodal_var_index, 1, st, cnt,
                                  nodal_vars)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read nodal variable from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write nodal variable slab
 */
void F2C(expnvs, EXPNVS)(int *idne, int *time_step, int *nodal_var_index, void_int *start,
                         void_int *count, real *nodal_var_vals, int *ierr)
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

  if ((*ierr = ex_put_partial_var(*idne, *time_step, EX_NODAL, *nodal_var_index, 1, st, cnt,
                                  nodal_var_vals)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write nodal variable slab to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element numbering map
 */
void F2C(exgnenm, EXGNENM)(int *idne, void_int *starte, void_int *num_ent, void_int *elem_map,
                           int *ierr)
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

  if ((*ierr = ex_get_partial_id_map(*idne, EX_ELEM_MAP, st, cnt, elem_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read element numbering map from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the element numbering map
 */
void F2C(expnenm, EXPNENM)(int *idne, void_int *starte, void_int *num_ent, void_int *elem_map,
                           int *ierr)
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

  if ((*ierr = ex_put_partial_id_map(*idne, EX_ELEM_MAP, st, cnt, elem_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write element numbering map to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the node numbering map
 */
void F2C(exgnnnm, EXGNNNM)(int *idne, void_int *startn, void_int *num_ent, void_int *node_map,
                           int *ierr)
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

  if ((*ierr = ex_get_partial_id_map(*idne, EX_NODE_MAP, st, cnt, node_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read node numbering map from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the node numbering map
 */
void F2C(expnnnm, EXPNNNM)(int *idne, void_int *startn, void_int *num_ent, void_int *node_map,
                           int *ierr)
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

  if ((*ierr = ex_put_partial_id_map(*idne, EX_NODE_MAP, st, cnt, node_map)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write node numbering map to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the node map for a processor
 */
void F2C(exgnmp, EXGNMP)(int *idne, void_int *node_mapi, void_int *node_mapb, void_int *node_mape,
                         int *processor, int *ierr)
{
  if ((*ierr = ex_get_processor_node_maps(*idne, node_mapi, node_mapb, node_mape, *processor)) !=
      0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read processor node map from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write a node map for a processor
 */
void F2C(expnmp, EXPNMP)(int *idne, void_int *node_mapi, void_int *node_mapb, void_int *node_mape,
                         int *processor, int *ierr)
{
  if ((*ierr = ex_put_processor_node_maps(*idne, node_mapi, node_mapb, node_mape, *processor)) !=
      0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write processor node map to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the element map for a processor
 */
void F2C(exgemp, EXGEMP)(int *idne, void_int *elem_mapi, void_int *elem_mapb, int *processor,
                         int *ierr)
{
  if ((*ierr = ex_get_processor_elem_maps(*idne, elem_mapi, elem_mapb, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read processor element map from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the element map for a processor
 */
void F2C(expemp, EXPEMP)(int *idne, void_int *elem_mapi, void_int *elem_mapb, int *processor,
                         int *ierr)
{
  if ((*ierr = ex_put_processor_elem_maps(*idne, elem_mapi, elem_mapb, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write processor element map to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the communications map parameters for a single processor
 */
void F2C(exgcmp, EXGCMP)(int *idne, void_int *ncmap_ids, void_int *ncmap_node_cnts,
                         void_int *ecmap_ids, void_int *ecmap_elem_cnts, int *processor, int *ierr)
{
  if ((*ierr = ex_get_cmap_params(*idne, ncmap_ids, ncmap_node_cnts, ecmap_ids, ecmap_elem_cnts,
                                  *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read comm map parameters from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the communications map parameters for a single processor
 */
void F2C(expcmp, EXPCMP)(int *idne, void_int *nmap_ids, void_int *nmap_node_cnts,
                         void_int *emap_ids, void_int *emap_elem_cnts, int *processor, int *ierr)
{
  if ((*ierr = ex_put_cmap_params(*idne, nmap_ids, nmap_node_cnts, emap_ids, emap_elem_cnts,
                                  *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write comm map parameters to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the communications map parameters for all processors
 */
void F2C(expcmpc, EXPCMPC)(int *idne, void_int *nmap_ids, void_int *nmap_node_cnts,
                           void_int *nproc_ptrs, void_int *emap_ids, void_int *emap_elem_cnts,
                           void_int *eproc_ptrs, int *ierr)
{
  if ((*ierr = ex_put_cmap_params_cc(*idne, nmap_ids, nmap_node_cnts, nproc_ptrs, emap_ids,
                                     emap_elem_cnts, eproc_ptrs)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write comm map parameters to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the nodal communications map for a single processor
 */
void F2C(exgncm, EXGNCM)(int *idne, entity_id *map_id, void_int *node_ids, void_int *proc_ids,
                         int *processor, int *ierr)
{
  if ((*ierr = ex_get_node_cmap(*idne, *map_id, node_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to read nodal communications map from file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the nodal communications map for a single processor
 */
void F2C(expncm, EXPNCM)(int *idne, entity_id *map_id, void_int *node_ids, void_int *proc_ids,
                         int *processor, int *ierr)
{
  if ((*ierr = ex_put_node_cmap(*idne, *map_id, node_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Error: failed to write nodal communications map to file id %d", *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Read the elemental communications map for a single processor
 */
void F2C(exgecm, EXGECM)(int *idne, entity_id *map_id, void_int *elem_ids, void_int *side_ids,
                         void_int *proc_ids, int *processor, int *ierr)
{
  if ((*ierr = ex_get_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to read elemental comm map from file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}

/*
 * Write the elemental communications map for a single processor
 */
void F2C(expecm, EXPECM)(int *idne, entity_id *map_id, void_int *elem_ids, void_int *side_ids,
                         void_int *proc_ids, int *processor, int *ierr)
{
  if ((*ierr = ex_put_elem_cmap(*idne, *map_id, elem_ids, side_ids, proc_ids, *processor)) != 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "Error: failed to write elemental comm map to file id %d",
             *idne);
    ex_err_fn(*idne, __func__, errmsg, EX_MSG);
  }
}
