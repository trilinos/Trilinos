/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"

/*****************************************************************************/
void Zoltan_PHG_Plot_2D_Distrib(
  ZZ *zz,
  PHGraph *phg
)
{
/* Routine that produces gnuplot output of 2D data distribution in form of 
 * a matrix.
 * One column for each vertex.
 * One row for each hyperedge.
 * Separate files are produced for each processor.
 * Vertex and edge global node numbers are used for "coordinates" in plotting.
 * No partitioning information is displayed; only the 2D data distribution 
 * is shown.
 */
static int cnt = 0;
char filename[32];
FILE *fp = NULL;
int i, j;
int egno, vgno;

  sprintf(filename, "phg%02d.%02d", cnt, zz->Proc);
  fp = fopen(filename, "w");

  for (i = 0; i < phg->nEdge; i++) {
    egno = EDGE_LNO_TO_GNO(phg, i);
    for (j = phg->hindex[i]; j < phg->hindex[i+1]; j++) {
      vgno = VTX_LNO_TO_GNO(phg, phg->hvertex[j]);
      fprintf(fp, "%d  %d\n", vgno, -egno);
    }
  }
  fclose(fp);
  if (zz->Proc == 0) {
    sprintf(filename, "phg%02d.gnuload", cnt);
    fp = fopen(filename, "w");
    fprintf(fp, "set data style dots\n");
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set xlabel \"vertices\"\n");
    fprintf(fp, "set ylabel \"-hyperedges\"\n");
    fprintf(fp, "plot ");
    for (i = 0; i < zz->Num_Proc; i++) {
      fprintf(fp, "\"phg%02d.%02d\"", cnt, i);
      if (i != zz->Num_Proc-1)
        fprintf(fp, ", ");
      else
        fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Sanity check to ensure Mirror is working correctly */
  /* Don't need to generate both sets of files, but they should differ only 
   * in the order of the points */
  sprintf(filename, "phgmirror%02d.%02d", cnt, zz->Proc);
  fp = fopen(filename, "w");

  for (i = 0; i < phg->nVtx; i++) {
    vgno = VTX_LNO_TO_GNO(phg, i);
    for (j = phg->vindex[i]; j < phg->vindex[i+1]; j++) {
      egno = EDGE_LNO_TO_GNO(phg, phg->vedge[j]);
      fprintf(fp, "%d  %d\n", vgno, -egno);
    }
  }
  fclose(fp);
  if (zz->Proc == 0) {
    sprintf(filename, "phgmirror%02d.gnuload", cnt);
    fp = fopen(filename, "w");
    fprintf(fp, "set data style dots\n");
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set xlabel \"vertices\"\n");
    fprintf(fp, "set ylabel \"-hyperedges\"\n");
    fprintf(fp, "plot ");
    for (i = 0; i < zz->Num_Proc; i++) {
      fprintf(fp, "\"phgmirror%02d.%02d\"", cnt, i);
      if (i != zz->Num_Proc-1)
        fprintf(fp, ", ");
      else
        fprintf(fp, "\n");
    }
    fclose(fp);
  }

  cnt++;
}

/*****************************************************************************/
void Zoltan_PHG_Plot(
  int proc,         /* Processor calling the routine */
  int nvtx,         /* Number of vertices */
  int nparts,       /* Number of partitions; ignored if part == NULL */
  int *vindex,      /* Starting index in vedges of hyperedges for each vertex */
  int *vedge,       /* Hyperedges for each vertex */
  int *part,        /* Partition to which vertices are assigned; if NULL,
                       partition information is not plotted. */
  char *str         /* String included as comment in output files. */
)
{
/* Routine that produces gnuplot output of hypergraph in form of matrix.
 * One column for each vertex.
 * One row for each hyperedge.
 * Entries in matrix entry a[i,j] if vertex i is in hyperedge j.
 * If part == NULL, a single file is produced with all information.
 * (Currently, this capability is serial only.)
 * If part != NULL, a separate file is produced for each partition.
 * (Currently, this capability works in parallel for k >= p.)
 */
static int cnt = 0;
char filename[32];
int i, j, v;
FILE *fp = NULL;
int prev_part = -1;
int *idx = NULL;
int *vtx = NULL;

  idx = (int *) ZOLTAN_MALLOC(2 * nvtx * sizeof(int));
  vtx = idx + nvtx;
  for (i = 0; i < nvtx; i++) vtx[i] = idx[i] = i;
  if (part != NULL) {
    /* sort vertices by partition */
    Zoltan_quicksort_pointer_inc_int_int(idx, part, vtx, 0, nvtx-1);
  }
  else {
    /* write all results to one file */
    sprintf(filename, "hgplot%02d", cnt);
    fp = fopen(filename, "w");
    fprintf(fp, "#%s\n", str);
  }

  for (i = 0; i < nvtx; i++) {
    v = idx[i];
    if (part != NULL) {
      if (part[v] != prev_part) {
        /* open new file for this partition's data */
        if (fp != NULL) fclose(fp);
        sprintf(filename, "hgplot%02d.%02d", cnt, part[v]);
        fp = fopen(filename, "w");
        fprintf(fp, "#%s\n", str);
        prev_part = part[v];
      }
    }
    for (j = vindex[v]; j < vindex[v+1]; j++)
      fprintf(fp, "%d %d\n", vedge[j], -v);
  }

  fclose(fp);
  ZOLTAN_FREE(&idx);
  if (proc == 0) {
    sprintf(filename, "hgplot%02d.gnuload", cnt);
    fp = fopen(filename, "w");
    fprintf(fp, "set data style dots\n");
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set title \"%s\"\n", str);
    fprintf(fp, "set xlabel \"hyperedges\"\n");
    fprintf(fp, "set ylabel \"-vertices\"\n");
    fprintf(fp, "plot ");
    if (part != NULL) {
      for (i = 0; i < nparts; i++) {
        fprintf(fp, "\"hgplot%02d.%02d\"", cnt, i);
        if (i != nparts-1)
          fprintf(fp, ", ");
        else
          fprintf(fp, "\n");
      }
    }
    else {
      fprintf(fp, "\"hgplot%02d\"\n", cnt);
    }
    fclose(fp);
  }
  cnt++;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
