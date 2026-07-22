// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

#include "dr_const.h"
#include "dr_hg_readfile.h"
#include "dr_util_const.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define BUF_LEN 1000
#define ERROR(proc,yo,msg,err) \
 {printf("Proc(%d) From(%s): %s\n",proc,yo,msg);return(err);}

static int   readfile  (int, ZOLTAN_FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);

static int nextstr (ZOLTAN_FILE* f, char *string);


/* These routines work with multiple (vector) vertex and hyperedge weights */

/*****************************************************************************/

int HG_readfile (
  int Proc,
  ZOLTAN_FILE* f,
  int *nVtx, int *nEdge, int *nPins,
  int **hindex,  int **hvertex,
  int *vwgt_dim, float **vwgt,
  int *ewgt_dim, float **ewgt,
  int *base)
{
  char string[BUF_LEN];
  const char *yo = "HG_readfile";
  int i, code;

  /* Initialize return values in case of error. */
  *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = *base = 0;
  *hindex = *hvertex = NULL;
  *vwgt = *ewgt = NULL;

  /* Skip comments. */
  while (1) {
    ZOLTAN_FILE_gets (string, BUF_LEN-1, f);

    if (*string != '%')
      break;
  }

  /* Extended format: nVtx nEdge nPins NM where both N&M are [0,9]. N is the
   * number of dimension for vertex weights and M is number of dimensions for
   * the edge weights. */
  if (sscanf(string, "%d %d %d %d", nVtx,nEdge,nPins,&code) != 4)
    ERROR(Proc, yo, "Unrecognized file format.", ZOLTAN_FATAL);

  if (*nVtx <= 1) {                  /* wrong guess if nVtx really is one! */
    i = sscanf(string, "%d %d %d %d %d %d", 
                       base, nVtx, nEdge, nPins, &code, vwgt_dim);
    if (i == 4)
      code = 0;
    else if (i < 4)
      ERROR(Proc, yo, "Error in PaToH file format.", ZOLTAN_FATAL);

    switch (code) {
    case 0:    /* No weights */
      *vwgt_dim = 0;
      *ewgt_dim = 0;
      break;
    case 1: /* Vwgts only; default is 1 vwgt unless 6 ints on input line */
      if (i == 5) *vwgt_dim = 1;
      *ewgt_dim = 0;
      break;
    case 2: /* Ewgts only; PaToH supports only one ewgt. */
      *vwgt_dim = 0;
      *ewgt_dim = 1;
      break;
    case 3: /* Ewgts & vwgts; default is 1 vwgt unless 6 ints on input line */
      if (i == 5) *vwgt_dim = 1;
      *ewgt_dim = 1;
      break;
    default: /* Invalid code */
      ERROR(Proc, yo, "Error in PaToH file format:  Illegal PaToH code.",
            ZOLTAN_FATAL);
      break;
    }

    return readfile(Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim,
                    vwgt, ewgt_dim, ewgt, base);          /* patoh format */
  }
  else {
    *base = 1;
    if (sscanf(string, "%d %d %d %1d%1d",
                        nVtx,nEdge,nPins,vwgt_dim,ewgt_dim) != 5)
      ERROR(Proc, yo, "Error in HG file format.", ZOLTAN_FATAL);
    return readfile(Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim,
                     vwgt, ewgt_dim, ewgt, base);          /* IBM format */
  }
}


/*****************************************************************************/
/* Read IBM & Patoh formats */

static int readfile (
 int Proc,
 ZOLTAN_FILE* f,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
{
  int i, j, pin;
  char string[BUF_LEN];
  float *p;
  const char *yo = "readfile";

  if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
   || !(*vertex = (int*) malloc ((*nPins    * sizeof(int)))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR);

  if (*vwgt_dim > 0
   && !(*vwgt = (float*) malloc (*nVtx * *vwgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)

  if (*ewgt_dim > 0
   && !(*ewgt = (float*) malloc (*nEdge * *ewgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)

  pin = 0;
  p = *ewgt;
  for (i = 0; i < *nEdge; i++) {
     (*index)[i] = pin;
     for (j = 0; j < *ewgt_dim; j++) {
        if (nextstr (f, string) != 1)
           ERROR (Proc, yo, "Wrong number of edge weights", ZOLTAN_FATAL)
        *p++ = atof (string);
     }
     while (nextstr (f,string) == 1)
        (*vertex)[pin++] = atoi(string);
  }
  (*index)[i] = pin;
  if (pin != *nPins)
     ERROR (Proc, yo, "Assertion Error in readfile.", ZOLTAN_FATAL);

  p = *vwgt;
  if (*vwgt_dim > 0)
     for (i = 0; i < *nVtx; i++) {
        for (j = 0; nextstr (f, string) == 1 && j < *vwgt_dim; j++)
           *p++ = atof (string);
        if (j != *vwgt_dim)
	  ERROR (Proc, yo, "Wrong number of vertex weights", ZOLTAN_FATAL)
     }

  if (nextstr (f, string) == 1)
     ERROR (Proc, yo, "Input file is longer than expected", ZOLTAN_FATAL)
  return ZOLTAN_OK;
}


/*****************************************************************************/
/* Uses character reads to read arbitrary long hyperedge definition lines.   */
static int nextstr (ZOLTAN_FILE* f, char *string)
{
  char ch;

  *string = 0;
  while (1) {
    ch = ZOLTAN_FILE_getc(f);
    if (ch == EOF)          return 0;
    if (ch == '\n')         return -1;
    if (isspace(ch) == 0)   break;
  }

  while (1)  {
    if (ch == EOF || isspace(ch)) {
      *string = 0;
      ZOLTAN_FILE_ungetc (ch, f);  /* ??? */
      return 1;
    }
    *string++ = ch;
    ch = ZOLTAN_FILE_getc(f);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
