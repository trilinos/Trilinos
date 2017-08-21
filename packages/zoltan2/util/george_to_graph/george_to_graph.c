/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: mm_to_graph.c,v $
 *    $Author: kddevin $
 *    $Date: 2005/03/01 00:43:34 $
 *    $Revision: 1.1 $
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

/* usage:  a.out matrix.mtx     -- for row-based decomposition */


static void Zoltan_quicksort_pointer_inc_int_int(int *, int *, int *, int, int);
int weight_generator(int, int);


int main (int argc, char *argv[]) {
  MM_typecode mmtype;

  int cnt, nz;
  int jj, kk;
  int *I = NULL, *J = NULL, *idx = NULL;
  int *vtxs, *nbors;
  int *vindex = NULL, *symvtxs = NULL, *symnbors = NULL;
  int nvtx;
  int M, N, prevv, v;
  double junkval;
  char filename[50];
  FILE *fp;

  /* Read matrix information */

  int vertex_weights_num = atoi(argv[2]);

  printf("KDDKDD HELLO\n");
  sprintf(filename, "%s", argv[1]);
  printf("Reading MM matrix %s with %d vertex weights\n", filename, vertex_weights_num);

  fp = fopen(filename, "r");
  /* mm_read_banner(fp, &mmtype); */
  /* mm_read_mtx_crd_size(fp, &M, &N, &nz); */
  fscanf(fp, "%d %d", &M, &nz);
  N = M;
  printf("KDDKDD %d %d \n", M, nz);

  I = (int *) malloc(4 * nz * sizeof(int));
  J = I + nz;
  idx = J + nz;

  for (kk = 0; kk < nz; kk++) {
    fscanf(fp, "%d %d", &I[kk], &J[kk]);
    idx[kk] = kk;
  }
  fclose(fp);

  /* vertices == rows
  YE: extracts the name of the graph file, leaving out the .graph term
  YE: strtok is a weird function: https://stackoverflow.com/questions/21097253/how-does-the-strtok-function-in-c-work */

  char * graph_name = strtok(argv[1], ".");

  /* renames to <graph name>>_<# of vertex weights>vw.graph */
  sprintf(filename, "%s_%dvw.graph", graph_name, vertex_weights_num);
  nvtx = M;
  vtxs = I;
  nbors = J;

  /* Sort the edges by vertex */
  Zoltan_quicksort_pointer_inc_int_int(idx, vtxs, nbors, 0, nz-1);

  /* Symmetrize the matrix */
  printf("Symmetrizing MM matrix %s\n", filename);

  /* Create vindex for faster searching */
  vindex = (int *) calloc((2 + nvtx), sizeof(int));
  for (kk = 0; kk < nz; kk++)
    vindex[vtxs[kk]]++;

  cnt = 0;
  for (kk = 0; kk < 2+nvtx; kk++) {
    cnt += vindex[kk];
    vindex[kk] = cnt - vindex[kk];
  }

  symvtxs = (int *) malloc(4 * nz * sizeof(int));
  symnbors = symvtxs + 2 * nz;

  cnt = 0;
  for (kk = 0; kk < nz; kk++) {
    int i, j, k, found;
    k = idx[kk];
    i = symvtxs[cnt] = vtxs[k];
    j = symnbors[cnt] = nbors[k];
    cnt++;
    /* for each pin, see whether the transpose pin exists */
    found = 0;
    for (jj = vindex[j]; !found && jj < vindex[j+1]; jj++) {
      if (nbors[idx[jj]] == i) found = 1;
      if (nbors[idx[jj]] > i) break;  /* nbors are sorted; early out
                                         shortens search */
    }

    if (!found) {
      /* Transpose wasn't found; add it. */
      symvtxs[cnt] = j;
      symnbors[cnt] = i;
      cnt++;
    }
  }

  if (cnt != nz) {
    /* Pins added; need to resort. */
    vtxs = symvtxs;
    nbors = symnbors;
    nz = cnt;
    for (kk = 0; kk < cnt; kk++) idx[kk] = kk;
    Zoltan_quicksort_pointer_inc_int_int(idx, vtxs, nbors, 0, nz-1);
  }

  /* Write the graph file. */
  printf("Writing MM graph %s\n", filename);

  fp = fopen(filename, "w");

  /* Writes the first line as <number of vertices> <number of edges> 010 <number of vertex weights>
  where "010" corresponds to "vertex weights, no edge weights".*/

  fprintf(fp, "%d %d 010 %d\n", nvtx, nz, vertex_weights_num);

  prevv = -1;

  // Iterates through each line of George's graph file.
  for (kk = 0; kk < nz; kk++) {
    v = vtxs[idx[kk]];
    if (v != prevv) {
      if (prevv != -1) fprintf(fp, "\n");
      prevv = v;

      /* Writes the vertex weights by component. You can set your choice of weights in the weight_generator function */
      for(int wgt = 0; wgt < vertex_weights_num; wgt++)
      {
          fprintf(fp, "%d ", weight_generator(wgt, v));
      }
    }
    fprintf(fp, "%d ", nbors[idx[kk]]+1);  // Chaco is one-based; George's is zero-based
  }
  fprintf(fp, "\n");

  fclose(fp);

  free(I);
  free(vindex);
  free(symvtxs);
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

int weight_generator(int wgt_component, int vertex_index)
{
    if(wgt_component == 0)
    {
        return 1;
    }
    else if(wgt_component == 1)
    {
        //Chaco is 1-based, George's is 0-based
        return vertex_index + 1;
    }
    else
    {
        return rand()%100 + 1;
    }
}
