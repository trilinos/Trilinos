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
#include <float.h>
#include <limits.h>


/*****************************************************************************/



void Zoltan_PHG_HGraph_Init(
  PHGraph *hgraph
)
{
  hgraph->info  = 0;
  hgraph->nVtx  = 0;
  hgraph->nEdge = 0;
  hgraph->nNonZero= 0;
  hgraph->nDim  = 0;
  hgraph->EdgeWeightDim   = 0;
  hgraph->VertexWeightDim = 0;
  hgraph->ratio = 0.5;

  hgraph->coor    = NULL;
  hgraph->vwgt    = NULL;
  hgraph->ewgt    = NULL;
  hgraph->hindex  = NULL;
  hgraph->hvertex = NULL;
  hgraph->vindex  = NULL;
  hgraph->vedge   = NULL;
  hgraph->dist_x  = NULL;
  hgraph->dist_y  = NULL;
  hgraph->vmap    = NULL;
}



/*****************************************************************************/

void Zoltan_PHG_Graph_Init(
  PGraph *graph
)
{
  graph->info  = 0;
  graph->nVtx  = 0;
  graph->nEdge = 0;
  graph->nDim  = 0;
  graph->EdgeWeightDim   = 0;
  graph->VertexWeightDim = 0;

  graph->vtxdist = NULL;
  graph->coor    = NULL;
  graph->vwgt    = NULL;
  graph->ewgt    = NULL;
  graph->nindex  = NULL;
  graph->neigh   = NULL;
}



/****************************************************************************/
/* Frees memory associated with a hypergraph; but not the hypergraph itself */

int Zoltan_PHG_HGraph_Free(
  PHGraph *hg
)
{
  if (hg)
    Zoltan_Multifree (__FILE__, __LINE__, 10, &hg->coor, &hg->vwgt, &hg->ewgt,
     &hg->hindex, &hg->hvertex, &hg->vindex, &hg->vedge, &hg->dist_x,
     &hg->dist_y, &hg->vmap);

  return ZOLTAN_OK;
}



/****************************************************************************/
/* Frees all memory associated with a graph; does not free the graph itself. */
int Zoltan_PHG_Graph_Free(
  PGraph *g
)
{
  if (g)
    Zoltan_Multifree (__FILE__, __LINE__, 6, &g->coor, &g->vwgt, &g->ewgt,
     &g->vtxdist, &g->nindex, &g->neigh);

  return ZOLTAN_OK;
}



/****************************************************************************/

int Zoltan_PHG_Info (
  ZZ *zz,
  PHGraph *hg
)
{
  int i, size, size_min, size_max, count;
  float wgt_min, wgt_max, wgt_tot;
  double mean, var, temp;
  char *yo = "Zoltan_PHG_Info";

  ZOLTAN_TRACE_ENTER(zz, yo);

  printf("--------- PHGraph Information (min/ave/max/tot) ------------------\n");
  printf("INFO MAY BE WRONG FOR 2D DATA DISTRIBUTION. KDDKDD \n");
  printf("info:%d |V|=%d |E|=%d |P|=%d \n", hg->info, hg->nVtx, hg->nEdge,
   hg->nNonZero);

  /* print weights */
  if (hg->vwgt) {
    wgt_tot = 0.0;
    wgt_min = FLT_MAX;
    wgt_max = FLT_MIN;
    for (i = 0; i < hg->nVtx; i++) {
      wgt_tot += hg->vwgt[i];
      wgt_min = MIN(wgt_min, hg->vwgt[i]);
      wgt_max = MAX(wgt_max, hg->vwgt[i]);
    }
    printf("Vertex weights   :    %9.2f %9.2f %9.2f %12.2f\n", wgt_min,
     wgt_tot/hg->nVtx, wgt_max, wgt_tot);

    mean = wgt_tot / hg->nVtx;
    var = 0.0;
    for (i = 0; i < hg->nVtx; i++) {
      temp = hg->vwgt[i] - mean;
      var += (temp*temp);
    }
    var = sqrt(var/(hg->nVtx-1));
    printf ("Vertex Stats: stdev %.2f,   Coef of Var %.2f\n", var, var/mean);
    count=0;
    temp=0.0;
    for (i = 0; i < hg->nVtx; i++)
      if (hg->vwgt[i] > (mean + 3.0 * var)) {
        count++;
        temp += hg->vwgt[i];
      }
    printf ("Vertex >3sigma: count %d, rel count x10000 %.1f, "
     "rel weight x100 %.1f\n",
     count, (float)10000*count/(float)hg->nVtx, 100.0*temp/wgt_tot);
  }
  if (hg->ewgt) {
    wgt_tot = 0.0;
    wgt_min = FLT_MAX;
    wgt_max = FLT_MIN;
    for (i = 0; i < hg->nEdge; i++) {
      wgt_tot += hg->ewgt[i];
      wgt_min = MIN(wgt_min, hg->ewgt[i]);
      wgt_max = MAX(wgt_max, hg->ewgt[i]);
    }
    printf("HEdge weights    :    %9.2f %9.2f %9.2f %12.2f\n", wgt_min,
     wgt_tot/hg->nEdge, wgt_max, wgt_tot);

    var = 0.0;
    mean = wgt_tot / hg->nEdge;
    for (i = 0; i < hg->nEdge; i++) {
      temp = hg->ewgt[i] - mean;
      var += (temp*temp);
    }
    var = sqrt(var/(hg->nEdge-1));
    printf ("HEdge Stats: STDV %.2f,   Coef of Var %.2f\n", var, var/mean);
  }

  /* print sizes */
  if (hg->hindex) {
    size_min = INT_MAX;
    size_max = INT_MIN;
    for (i = 0; i < hg->nEdge; i++) {
      size     = hg->hindex[i+1] - hg->hindex[i];
      size_min = MIN(size_min, size);
      size_max = MAX(size_max, size);
    }
    printf("Edge sizes       :    %6d    %9.2f %6d    %9d\n", size_min,
     (float)hg->nNonZero / hg->nEdge, size_max, hg->nNonZero);

    var = 0.0;
    mean = (float)hg->nNonZero / hg->nEdge;
    for (i = 0; i < hg->nEdge; i++) {
      temp = (float)(hg->hindex[i+1]-hg->hindex[i]) - mean;
      var += (temp*temp);
    }
    var = sqrt(var/(hg->nEdge-1));
    printf ("Edge Stats: stdev %.2f,   Coef of Var %.2f\n", var, var/mean);
  }
  if (hg->vindex) {
    size_min = INT_MAX;
    size_max = INT_MIN;
    for (i = 0; i < hg->nVtx; i++) {
      size     = hg->vindex[i+1] - hg->vindex[i];
      size_min = MIN(size_min, size);
      size_max = MAX(size_max, size);
    }
    printf("Vertex sizes     :    %6d    %9.2f %6d    %9d\n", size_min,
     (float)hg->nNonZero / hg->nVtx, size_max, hg->nNonZero);

    var = 0.0;
    mean = (float)hg->nNonZero / hg->nVtx;
    for (i = 0; i < hg->nVtx; i++) {
      temp = (float)(hg->vindex[i+1]-hg->vindex[i]) - mean;
      var += (temp*temp);
    }
    var = sqrt(var/(hg->nVtx-1));
    printf ("Vertex Stats: stdev %.2f,  Coef of Var %.2f\n", var, var/mean);
  }

  printf("-----------------------------------------------------------------\n");

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}



/****************************************************************************/
/**  Given arrays to lookup a vertex for a given hyperedge (hindex, hvertex)
 **  or arrays to lookup a hyperedge given a vertex (vindex, vedge) create the
 **  other set of arrays (mirror).
 **  The mirror (to be filled in) should have NULL pointers in hg on entry!
 **  The calling program is responsible for freeing the mirror's memory  ***/

int Zoltan_PHG_Create_Mirror (
  ZZ *zz,
  PHGraph *hg
)
{
  int i, j;                  /* loop counters */
  int inlength, outlength;   /* input/output array lengths */
  int *index, *data;         /* pointers to input information */
  int *outindex, *outdata;
  char *yo = "Zoltan_PHG_Create_Mirror";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* determine which data to "mirror" and set corresponding data pointers. */
  if (hg &&  hg->hindex && (hg->nNonZero == 0 || hg->hvertex)
   && !hg->vindex && !hg->vedge) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Have hindex; building vindex.");

    inlength  = hg->nEdge;
    outlength = hg->nVtx;
    index     = hg->hindex;
    data      = hg->hvertex;
    outindex  = hg->vindex = (int*) ZOLTAN_MALLOC((hg->nVtx+1) * sizeof(int));
    outdata   = hg->vedge  = (int*) ZOLTAN_MALLOC (hg->nNonZero  * sizeof(int));

    if (outindex == NULL || (hg->nNonZero > 0 && outdata == NULL)) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &hg->vindex, &hg->vedge);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  else if (hg && hg->vindex && (hg->nNonZero == 0 || hg->vedge)
   && !hg->hindex && !hg->hvertex) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Have vindex; building hindex.");

    inlength  = hg->nVtx;
    outlength = hg->nEdge;
    index     = hg->vindex;
    data      = hg->vedge;
    outindex  = hg->hindex  = (int*) ZOLTAN_MALLOC((hg->nEdge+1) * sizeof(int));
    outdata   = hg->hvertex = (int*) ZOLTAN_MALLOC(hg->nNonZero    * sizeof(int));

    if (outindex == NULL || (hg->nNonZero > 0 && outdata == NULL)) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &hg->hindex, &hg->hvertex);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input error.");
    return ZOLTAN_FATAL;  /* unable to proceed */
  }

  /* count number of data objects in the outindex space */
  for (i = 0; i < outlength+1; i++)
    outindex[i] = 0;
  for (i = 0; i < inlength; i++)
    for (j = index[i];  j < index[i+1]; j++)
      (outindex[data[j]+1])++;

  /* compute partial sum */
  for (i = 1; i < outlength; i++)
    outindex [i] += outindex[i-1];

  /* store data in out index location and increment index */
  for (i = 0; i < inlength; i++)
    for (j = index[i]; j < index[i+1]; j++)
      outdata [(outindex[data[j]])++ ] = i;

  /* out index now points past end of sublist. Shift array one position up */
  for (i = outlength; i > 0; i--)
    outindex[i] = outindex[i-1];
  outindex[0] = 0;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}



/****************************************************************************/
/* check that (hindex, hvertex) and (vindex, vedge) are consistant mappings */
/* returns ZOLTAN_WARN if not consistant, returns ZOLTAN_OK otherwise */

int Zoltan_PHG_Check (
  ZZ *zz,
  PHGraph *hg
)
{
  int i;
  int iedge;               /* iedge denotes an hyperedge */
  int j;                   /* j is the index of a vertex */
  int k;                   /* k is an index of hyperedge */
  int *check;
  char str[256];
  int err = ZOLTAN_OK;
  char *yo = "Zoltan_PHG_Check";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (!hg->hindex || !hg->hvertex || !hg->vindex || !hg->vedge) {
    err = ZOLTAN_WARN;
    goto End;
  }

  check = (int*) ZOLTAN_MALLOC(((hg->nEdge > hg->nVtx) ? hg->nEdge : hg->nVtx)
   * sizeof(int));
  if (check == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to allocate memory.");
    err = ZOLTAN_MEMERR;
    goto End;
  }

  for (i = 0; i < hg->nEdge; i++)
    check[i] = -1;
  for (i = 0; i < hg->nVtx; i++)
    for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
      if (check[hg->vedge[j]] < i)
        check[hg->vedge[j]] = i;
      else {
        ZOLTAN_PRINT_WARN(zz->Proc, yo,"Found multiple hedge for same vertex.");
        err = ZOLTAN_WARN;
      }
  for (i = 0; i < hg->nVtx; i++)
    check[i] = -1;
  for (i = 0; i < hg->nEdge; i++)
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
      if (check[hg->hvertex[j]] < i)
        check[hg->hvertex[j]] = i;
      else {
        ZOLTAN_PRINT_WARN(zz->Proc, yo,"Found multiple vertex for same hedge.");
        err =  ZOLTAN_WARN;
      }
  ZOLTAN_FREE ((void**) &check);

  for (i = 0; i < hg->VertexWeightDim * hg->nVtx; i += hg->VertexWeightDim)
    if (hg->vwgt[i] < 0.0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative vertex weight.");
      err = ZOLTAN_WARN;
    }

  for (i = 0; i < hg->EdgeWeightDim * hg->nEdge; i += hg->EdgeWeightDim)
    if (hg->ewgt[i] < 0.0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative edge weight.");
      err = ZOLTAN_WARN;
    }

  if (hg->comm->nProc_x == 1) { 
    /* In 2D distribution, check makes sense only if proc has entire hedge */
    for (i = 0; i < hg->nEdge; i++)
      if ((hg->hindex[i+1] - hg->hindex[i]) < 2) {
        sprintf (str, "Found hedge with less than two vertices: "
         "edge %d has %d vtxs\n", i, (hg->hindex[i+1] - hg->hindex[i]));
        ZOLTAN_PRINT_WARN(zz->Proc, yo, str);
        err = ZOLTAN_WARN;
      }
  }

  for (i = 0; i < hg->nEdge; i++)
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
      if (hg->hvertex[j] < 0  ||  hg->hvertex[j] >= hg->nVtx) {
        ZOLTAN_PRINT_WARN(zz->Proc, yo,"Found vertex out of range in hvertex.");
        err = ZOLTAN_WARN;
      }

  for (i = 0; i < hg->nVtx; i++)
    for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
      if (hg->vedge[j] < 0  ||  hg->vedge[j] >= hg->nEdge) {
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found edge out of range in vedge.");
        err = ZOLTAN_WARN;
      }

  /* starting from (hindex,hvertex), for each edge determine each associated
   * vertex. Then for each vertex lookup associated edges using (vindex, vedge)
   */
  for (iedge = 0; iedge < hg->nEdge; iedge++)
    for (j = hg->hindex[iedge]; j < hg->hindex[iedge+1]; j++) {
      /* for each hyperedge get index to vertices */
      for (k=hg->vindex[hg->hvertex[j]]; k<hg->vindex[hg->hvertex[j]+1]; k++)
        /* for each vertex of current hyperedge get index to hyperedges */
        if (hg->vedge[k] == iedge)    /* does it match with original edge? */
          break;
      if (k == hg->vindex[hg->hvertex[j]+1])  /* if no match was found then */
        return ZOLTAN_WARN;                  /* failure, else keep on */
    }

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}




/****************************************************************************/
/*
 *  Converts a graph g into a hypergraph hg.
 *  One hyperedge is created for each vertex of g.
 *  Hyperedge i consists of vertex i + all edge neighbors of i in graph g. */

int Zoltan_PHG_Graph_to_HGraph(
  ZZ *zz,
  PGraph *g,        /* Input graph */
  PHGraph *hg       /* Ouput hypergraph */
)
{
  int i, j;
  int *hindex = NULL, *hvertex = NULL;  /* temporary array pointers */
  int cnt;
  int err = ZOLTAN_OK;
  char *yo = "Zoltan_PHG_Graph_to_HGraph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  ZOLTAN_PRINT_WARN(zz->Proc, yo, 
               "THIS ROUTINE IS NOT CORRECT FOR 2D DATA DISTRIBUTION KDDKDD");

  Zoltan_PHG_HGraph_Init(hg);
  hg->info  = g->info;
  hg->nVtx  = g->nVtx;
  hg->nEdge = g->nVtx;
  hg->nDim  = g->nDim;
  hg->nNonZero = g->nVtx + g->nEdge;
  hg->VertexWeightDim = g->VertexWeightDim;
  hg->EdgeWeightDim   = g->EdgeWeightDim;

  /* Copy vtxdist array; it is the same for hypergraph as graph. */
  /* KDDKDD NOT FOR 2D DATA DISTRIBUTION */
  if (g->vtxdist) {
    if (!(hg->dist_x = (int *) ZOLTAN_MALLOC((zz->Num_Proc+1)*sizeof(int)))) {
      err = ZOLTAN_MEMERR;
      goto End;
    }
    memcpy(hg->dist_x, g->vtxdist, (zz->Num_Proc+1) * sizeof(int));
  }

  /* Copy coordinates from graph to hypergraph */
  if (g->coor) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Copying coordinates");
    cnt = hg->nVtx * g->nDim;
    if (!(hg->coor = (double *) ZOLTAN_MALLOC(cnt * sizeof(double)))) {
      err = ZOLTAN_MEMERR;
      goto End;
    }
    memcpy(hg->coor, g->coor, cnt * sizeof(double));
  }

  /* Copy vertex weights from graph to hypergraph */
  if (g->vwgt) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Copying vertex weights");
    cnt = hg->nVtx * hg->VertexWeightDim;
    if (!(hg->vwgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float)))) {
      err = ZOLTAN_MEMERR;
      goto End;
    }
    memcpy(hg->vwgt, g->vwgt, cnt * sizeof(float));
  }

  /* Allocate memory for hypergraph edges */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge arrays");
  if (hg->nEdge > 0 && hg->nNonZero > 0)
    if (!(hindex  = hg->hindex =(int*)ZOLTAN_MALLOC((hg->nEdge+1)*sizeof(int)))
     || !(hvertex = hg->hvertex=(int*)ZOLTAN_MALLOC(hg->nNonZero*sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &hindex, &hvertex);
      err = ZOLTAN_MEMERR;
      goto End;
    }

  /* at the moment, the hyperedge weight will be constant !
  if (g->ewgt) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge weights");
    cnt = hg->nEdge * hg->EdgeWeightDim;
    if (!(hg->ewgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float)))) {
      err = ZOLTAN_MEMERR;
      goto End;
    }
  }
  */

  /* KDD  The following section will have to be reworked for parallel HGs. */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Filling edge arrays");
  cnt = 0;
  for (i = 0; i < hg->nEdge; i++) {
    hindex[i] = cnt;
    hvertex[cnt++] = i;
    for (j = g->nindex[i]; j < g->nindex[i+1]; j++) {
      hvertex[cnt++] = g->neigh[j];
      /* at the moment, the hyperedge weight will be constant
      for (k = 0; k < hg->EdgeWeightDim; k++)
        hg->ewgt[i*(hg->EdgeWeightDim)+k] += g->ewgt[j*(hg->EdgeWeightDim)+k];
      */
    }
  }
  hindex[hg->nEdge] = cnt;

  /* Sanity check */
  if (hg->nNonZero != cnt) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Sanity check failed: nNonZero != cnt.");
    err = ZOLTAN_FATAL;
    goto End;
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Creating mirror");
  err = Zoltan_PHG_Create_Mirror(zz, hg);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in building mirror.");

End:
  if (err == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory.");
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    Zoltan_Multifree (__FILE__, __LINE__, 5, &hg->coor, &hg->vwgt, &hg->ewgt,
     &hg->hindex, &hg->hvertex);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}

/****************************************************************************/



/* Routine to print hypergraph weights and edges. Assumes serial execution;
 * put inside Zoltan_Print_Sync_Start/Zoltan_Print_Sync_End for parallel
 * programs. */
void Zoltan_PHG_Print(
  ZZ *zz,
  PHGraph *hg,
  FILE *fp
)
{
  int i, j;
  int num_vwgt;
  int num_ewgt;
  char *yo = "Zoltan_PHG_Print";

  if (hg == NULL)
    return;

  ZOLTAN_TRACE_ENTER(zz, yo);

  num_vwgt = hg->VertexWeightDim;
  num_ewgt = hg->EdgeWeightDim;
  if (hg->vwgt != NULL) {
    fprintf(fp, "Vertices: [weights])\n");
    for (i = 0; i < hg->nVtx; i++) {
      fprintf(fp, "%d:  [", i);
      for (j = 0; j < num_vwgt; j++)
        fprintf(fp, "%f ", hg->vwgt[i*num_vwgt + j]);
      fprintf(fp, "])\n");
    }
  }

  /* Print Hyperedge Info */
  fprintf(fp, "Hyperedges:  (vertices)\n");
  for (i = 0; i < hg->nEdge; i++) {
    fprintf(fp, "%d:  ", i);
    fprintf(fp, "(");
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
      fprintf(fp, "%d ", hg->hvertex[j]);
    fprintf(fp, ")\n");
  }

  if (hg->ewgt != NULL) {
    fprintf(fp, "Hyperedge Weights:  [weights]\n");
    for (i = 0; i < hg->nEdge; i++) {
      fprintf(fp, "%d:  ", i);
      fprintf(fp, "[");
      for (j = 0; j < num_ewgt; j++)
        fprintf(fp, "%f ", hg->ewgt[i*num_ewgt + j]);
      fprintf(fp, "])\n");
    }
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
}

/****************************************************************************/



/*  Random Number generator due to Knuth found in Numerical Recipes in C
   (2nd edition) by Press, Vetterling, Teukolsky, Flannery (Page 284.)
   Currently, for use in the hypergraph library only. Needed because different
   random number implimentations on different machines produced different
   answers! Hence this is a single portable, fast, algorithm with adequate
   random number generation. NOTE: this assumes 32 bit ints, previously
   these variables were unsigned long (as was the return value) which
   gave problems on stratus (which assumed 64 bit longs.) */

static unsigned int idum = 123456789U;

unsigned Zoltan_PHG_Rand (void) 
{
  return idum = (1664525U * idum) + 1013904223U;
}



void Zoltan_PHG_Srand (unsigned int seed) 
{
  idum = seed;
}


/* Randomly permute an array of ints. */
void Zoltan_PHG_Rand_Perm_Int (int *data, int n)
{
  int i, number, temp;

  for (i = n; i > 0; i--) {
    number       = Zoltan_PHG_Rand() % i;
    temp         = data[number];
    data[number] = data[i-1];
    data[i-1]    = temp;
  }
} 



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
