/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mmio.h"

/* Karen Devine  Feb 28, 2005 */
/* Program to convert a MatrixMarket-formatted matrix to a Chaco graph file. */
/* For now, assumes MatrixMarket matrix is structurally symmetric! */
/* Can add support for non-symmetric matrices later. */

/* to compile:  gcc mm_to_graph.c mmio.c */

/* usage:  a.out matrix.mtx row     -- for row-based decomposition
           a.out matrix.mtx         -- for column-based decomposition */


static void Zoltan_quicksort_pointer_inc_int_int(int *, int *, int *, int, int);


int main (int argc, char *argv[]) {
  MM_typecode mmtype;
  int nz;
  int kk;
  int *I = NULL, *J = NULL, *idx = NULL;
  int *vtxs, *nbors;
  int nvtx;
  int M, N, prevv, v;
  double junkval;
  char filename[50];
  FILE *fp;
  
  /* Read matrix information */

  sprintf(filename, "%s.mtx", argv[1]);
  printf("Reading MM matrix %s\n", filename);

  fp = fopen(filename, "r");
  mm_read_banner(fp, &mmtype);
  mm_read_mtx_crd_size(fp, &M, &N, &nz);

  I = (int *) malloc(3 * nz * sizeof(int));
  J = I + nz;
  idx = J + nz;

  for (kk = 0; kk < nz; kk++) {
    if (mm_is_pattern(mmtype))
      fscanf(fp, "%d %d", &I[kk], &J[kk]);
    else 
      fscanf(fp, "%d %d %lf", &I[kk], &J[kk], &junkval);
    idx[kk] = kk;
  }
  fclose(fp);

  /*  Determine whether partitioning rows or columns. */
  /*  Default is columns (like zdrive) */

  if (argc >= 3 && !strcasecmp(argv[2], "row")) {
    /* vertices == rows */
    sprintf(filename, "%s_row.graph", argv[1]);
    nvtx = M;
    vtxs = I;
    nbors = J;
  }
  else {
    /* vertices == cols:  Sort based on col */
    sprintf(filename, "%s_col.graph", argv[1]);
    nvtx = N;
    vtxs = J;
    nbors = I;
  }

  /* Sort the edges by vertex */
  Zoltan_quicksort_pointer_inc_int_int(idx, vtxs, nbors, 0, nz-1);

  /* Write the graph file. */
  printf("Writing MM graph %s\n", filename);

  fp = fopen(filename, "w");

  fprintf(fp, "%d %d 00\n", nvtx, nz);
  prevv = -1;
  for (kk = 0; kk < nz; kk++) {
    v = vtxs[idx[kk]];
    if (v != prevv) {
      if (prevv != -1) fprintf(fp, "\n");
      prevv = v;
    }
    fprintf(fp, "%d ", nbors[idx[kk]]);
  }
  fprintf(fp, "\n");

  fclose(fp);

  free(I);
  return 0;
}

/*****************************************************************************/
/* Sorting pointers in increasing order. Sort key is int. Sub key is int. */

static void quickpart_pointer_inc_int_int (
int *sorted,
int *val1,
int *val2,
int start,
int end,
int *equal,
int *larger
)
{
int i, next, key1, key1_next, key2, key2_next;

  i = (end + start) / 2;
  key1 = val1 ? val1[sorted[i]] : 1;
  key2 = val2 ? val2[sorted[i]] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++) {
    next = sorted[i];
    key1_next = val1 ? val1[next] : 1;
    key2_next = val2 ? val2[next] : 1;
    if (key1_next < key1 || (key1_next == key1 && key2_next < key2)) {
      sorted[i]           = sorted[*larger];
      sorted[(*larger)++] = sorted[*equal];
      sorted[(*equal)++]  = next;
    }
    else if (key1_next == key1  &&  key2_next == key2) {
      sorted[i]           = sorted[*larger];
      sorted[(*larger)++] = next;
    }
  }
}

/* 
 * Sorts in increasing order with primary key val1 and secondary key val2.
 * The arrays val1 and val2 are not rearranged; rather the index array
 * sorted is rearranged based on values in val1 and val2. 
 */

/*****************************************************************************/
static void Zoltan_quicksort_pointer_inc_int_int (
int *sorted,   /* index array that is rearranged; should be initialized
                  so that sorted[i] == i. */
int* val1,     /* array of primary key values. */
int *val2,     /* array of secondary key values. */
int start,     /* first array position to be considered for sorting. */
int end        /* last array position to be considered for sorting. */
)
{
int  equal, larger;

  if (start < end) {
    quickpart_pointer_inc_int_int (sorted,val1,val2,start,end,&equal,&larger);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, start, equal-1);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, larger, end);
  }
}

