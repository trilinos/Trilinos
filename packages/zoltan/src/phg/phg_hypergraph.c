/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
#include <float.h>
#include <limits.h>
#include "phg_hypergraph.h"

/*****************************************************************************/
void Zoltan_HG_Graph_Init(
  Graph *graph
)
{
  graph->info  = 0;
  graph->nVtx  = 0;
  graph->nEdge = 0;
  graph->nDim  = 0;
  graph->EdgeWeightDim   = 0;
  graph->VtxWeightDim = 0;

  graph->vtxdist = NULL;

  graph->coor    = NULL;
  graph->vwgt    = NULL;
  graph->ewgt    = NULL;
  graph->nindex  = NULL;
  graph->neigh   = NULL;
}

/* Frees all memory associated with a graph; does not free the graph itself. */
int Zoltan_HG_Graph_Free(
  Graph *g
)
{
  if (g)
     Zoltan_Multifree (__FILE__, __LINE__, 6, &g->coor, &g->vwgt, &g->ewgt,
      &g->vtxdist, &g->nindex, &g->neigh);

  return ZOLTAN_OK;
}
/****************************************************************************/


/*************** Hypergraph Function *****************************************/
void Zoltan_HG_HGraph_Init(
  HGraph *phg
)
{
  phg->info  = 0;
  phg->nVtx  = 0;
  phg->nEdge = 0;
  phg->nPins = 0;
  phg->nRepartVtx = 0;
  phg->nRepartEdge = 0;
  phg->nRepartPin = 0;
  phg->nDim = 0;
  phg->EdgeWeightDim = 0;
  phg->VtxWeightDim = 0;
  phg->ratio = 0.5;
  phg->bisec_split = -1;  /* EBEB 6/7/06 */

  phg->comm    = NULL;
  phg->coor    = NULL;
  phg->vwgt    = NULL;
  phg->ewgt    = NULL;
  phg->hindex  = NULL;
  phg->hvertex = NULL;
  phg->vindex  = NULL;
  phg->vedge   = NULL;
  phg->dist_x  = NULL;
  phg->dist_y  = NULL;
  phg->vmap    = NULL;
  phg->fixed   = NULL;
}



/****************************************************************************/
/* Frees memory associated with a hypergraph; but not the hypergraph itself */

int Zoltan_HG_HGraph_Free(
  HGraph *hg
)
{
  if (hg){
/*
    Zoltan_Multifree (__FILE__, __LINE__, 11, 
     &hg->coor, &hg->vwgt, &hg->ewgt,
     &hg->hindex, &hg->hvertex, &hg->vindex, &hg->vedge, &hg->dist_x,
     &hg->dist_y, &hg->vmap, &hg->fixed);
*/

ZOLTAN_FREE(&hg->coor); 
ZOLTAN_FREE(&hg->vwgt); 
ZOLTAN_FREE(&hg->ewgt);
ZOLTAN_FREE(&hg->hindex); 
ZOLTAN_FREE(&hg->hvertex); 
ZOLTAN_FREE(&hg->vindex); 
ZOLTAN_FREE(&hg->vedge); 
ZOLTAN_FREE(&hg->dist_x);
ZOLTAN_FREE(&hg->dist_y); 
ZOLTAN_FREE(&hg->vmap); 
ZOLTAN_FREE(&hg->fixed);
  }
  return ZOLTAN_OK;
}






/****************************************************************************/

int Zoltan_HG_Info (
  ZZ *zz,
  HGraph *hg
)
{
  int i, size, size_min, size_max, count;
  float wgt_min, wgt_max, wgt_tot;
  double mean, var, temp;
  char *yo = "Zoltan_HG_Info";

  ZOLTAN_TRACE_ENTER(zz, yo);

  printf("--------- HGraph Information (min/ave/max/tot) ------------------\n");
  printf("INFO MAY BE WRONG FOR 2D DATA DISTRIBUTION. KDDKDD \n");
  printf("info:%d |V|=%d |E|=%d |P|=%d \n", hg->info, hg->nVtx, hg->nEdge,
   hg->nPins);

  /* print weights */
  if (hg->nVtx && hg->vwgt) {
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

    mean = var = 0.0;
    if (hg->nVtx > 1) {
      mean = wgt_tot / hg->nVtx;
      for (i = 0; i < hg->nVtx; i++) {
        temp = hg->vwgt[i] - mean;
        var += (temp*temp);
      }
      var = sqrt(var/(hg->nVtx-1));
      printf ("Vertex Stats: stdev %.2f,   Coef of Var %.2f\n", var, var/mean);
    }

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
  if (hg->nEdge && hg->ewgt) {
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

    if (hg->nEdge > 1) {
      var = 0.0;
      mean = wgt_tot / hg->nEdge;
      for (i = 0; i < hg->nEdge; i++) {
        temp = hg->ewgt[i] - mean;
        var += (temp*temp);
      }
      var = sqrt(var/(hg->nEdge-1));
      printf ("HEdge Stats: STDV %.2f,   Coef of Var %.2f\n", var, var/mean);
    }
  }

  /* print sizes */
  if (hg->nPins && hg->hindex) {
    size_min = INT_MAX;
    size_max = INT_MIN;
    for (i = 0; i < hg->nEdge; i++) {
      size     = hg->hindex[i+1] - hg->hindex[i];
      size_min = MIN(size_min, size);
      size_max = MAX(size_max, size);
    }
    printf("Edge sizes       :    %6d    %9.2f %6d    %9d\n", size_min,
     (float)hg->nPins / hg->nEdge, size_max, hg->nPins);

    if (hg->nEdge > 1) {
      var = 0.0;
      mean = (float)hg->nPins / hg->nEdge;
      for (i = 0; i < hg->nEdge; i++) {
        temp = (float)(hg->hindex[i+1]-hg->hindex[i]) - mean;
        var += (temp*temp);
      }
      var = sqrt(var/(hg->nEdge-1));
      printf ("Edge Stats: stdev %.2f,   Coef of Var %.2f\n", var, var/mean);
    }
  }
  if (hg->nPins && hg->vindex) {
    size_min = INT_MAX;
    size_max = INT_MIN;
    for (i = 0; i < hg->nVtx; i++) {
      size     = hg->vindex[i+1] - hg->vindex[i];
      size_min = MIN(size_min, size);
      size_max = MAX(size_max, size);
    }
    printf("Vertex sizes     :    %6d    %9.2f %6d    %9d\n", size_min,
     (float)hg->nPins / hg->nVtx, size_max, hg->nPins);

    if (hg->nVtx > 1) {
      var = 0.0;
      mean = (float)hg->nPins / hg->nVtx;
      for (i = 0; i < hg->nVtx; i++) {
        temp = (float)(hg->vindex[i+1]-hg->vindex[i]) - mean;
        var += (temp*temp);
      }
      var = sqrt(var/(hg->nVtx-1));
      printf ("Vertex Stats: stdev %.2f,  Coef of Var %.2f\n", var, var/mean);
    }
  }

  printf("-----------------------------------------------------------------\n");

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}


void Zoltan_HG_Mirror(
  int inlength,     /* Input:  either nVtx or nEdge */
  int *inindex,     /* Input:  index array, either vindex or hindex; 
                               length = inlength+1  */
  int *indata,      /* Input:  non-zeros array, either vedge or hvertex;
                               length = nPins */
  int outlength,    /* Input:  either nEdge or nVtx */
  int *outindex,    /* Output: index array, either hindex or vindex;
                               length = outlength+1 */
  int *outdata      /* Output: non-zeros array, either hvertex or vedge;
                               length = nPins */
)
{
/* Routine to invert arrays describing vertices and edges. 
 * Usually called by Zoltan_HG_Create_Mirror, with arguments 
 * determined by Zoltan_HG_Create_Mirror.
 * Arrays must be allocated to correct size before calling
 * the function.
 */
int i, j;

  /* count number of data objects in the outindex space */
  for (i = 0; i <= outlength; i++)
    outindex[i] = 0;
  for (i = 0; i < inlength; i++)
    for (j = inindex[i];  j < inindex[i+1]; j++)
      (outindex[indata[j]+1])++;

  /* compute partial sum */
  for (i = 1; i < outlength; i++)
    outindex[i] += outindex[i-1];

  /* store indata in outindex location and increment index */
  for (i = 0; i < inlength; i++)
    for (j = inindex[i]; j < inindex[i+1]; j++)
      outdata [(outindex[indata[j]])++ ] = i;

  /* out index now points past end of sublist. Shift array one position up */
  for (i = outlength; i > 0; i--)
    outindex[i] = outindex[i-1];
  outindex[0] = 0;
}


/****************************************************************************/
/**  Given arrays to lookup a vertex for a given hyperedge (hindex, hvertex)
 **  or arrays to lookup a hyperedge given a vertex (vindex, vedge) create the
 **  other set of arrays (mirror).
 **  The mirror (to be filled in) should have NULL pointers in hg on entry!
 **  The calling program is responsible for freeing the mirror's memory  ***/

int Zoltan_HG_Create_Mirror (
  ZZ *zz,
  HGraph *hg
)
{
  int inlength, outlength;   /* input/output array lengths */
  int *index, *data;         /* pointers to input information */
  int *outindex, *outdata;
  char *yo = "Zoltan_HG_Create_Mirror";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* determine which data to "mirror" and set corresponding data pointers. */
  if (hg &&  (hg->nEdge == 0 || hg->hindex) && (hg->nPins == 0 || hg->hvertex)
   && !hg->vindex && !hg->vedge) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Have hindex; building vindex.");

    inlength  = hg->nEdge;
    outlength = hg->nVtx;
    index     = hg->hindex;
    data      = hg->hvertex;
    outindex  = hg->vindex = (int*) ZOLTAN_MALLOC((hg->nVtx+1) * sizeof(int));
    outdata   = hg->vedge  = (int*) ZOLTAN_MALLOC (hg->nPins * sizeof(int));

    if (outindex == NULL || (hg->nPins > 0 && outdata == NULL)) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &hg->vindex, &hg->vedge);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  else if (hg && (hg->nVtx == 0 || hg->vindex) && (hg->nPins == 0 || hg->vedge)
   && !hg->hindex && !hg->hvertex) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Have vindex; building hindex.");

    inlength  = hg->nVtx;
    outlength = hg->nEdge;
    index     = hg->vindex;
    data      = hg->vedge;
    outindex  = hg->hindex  = (int*) ZOLTAN_MALLOC((hg->nEdge+1) * sizeof(int));
    outdata   = hg->hvertex = (int*) ZOLTAN_MALLOC(hg->nPins * sizeof(int));

    if (outindex == NULL || (hg->nPins > 0 && outdata == NULL)) {
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

  Zoltan_HG_Mirror(inlength, index, data, 
                    outlength, outindex, outdata);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}

/*****************************************************************************/



/****************************************************************************/
/* check that (hindex, hvertex) and (vindex, vedge) are consistant mappings */
/* returns ZOLTAN_WARN if not consistant, returns ZOLTAN_OK otherwise */

int Zoltan_HG_Check (
  ZZ *zz,
  HGraph *hg
)
{
  int i;
  int iedge;               /* iedge denotes an hyperedge */
  int j;                   /* j is the index of a vertex */
  int k;                   /* k is an index of hyperedge */
  int *check;
  char str[256];
  int err = ZOLTAN_OK;
  char *yo = "Zoltan_HG_Check";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if ((hg->nEdge && !hg->hindex) ||
      (hg->nVtx && !hg->vindex) || 
      (hg->nPins && (!hg->vedge || !hg->hvertex))) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "NULL arrays found");
    err = ZOLTAN_WARN;
    goto End;
  }

  if (hg->nPins) {
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
  }

  for (i = 0; i < hg->VtxWeightDim * hg->nVtx; i += hg->VtxWeightDim)
    if (hg->vwgt[i] < 0.0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative vertex weight.");
      err = ZOLTAN_WARN;
    }

  for (i = 0; i < hg->EdgeWeightDim * hg->nEdge; i += hg->EdgeWeightDim)
    if (hg->ewgt[i] < 0.0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative edge weight.");
      err = ZOLTAN_WARN;
    }

  if (hg->comm && hg->comm->nProc_x == 1) { 
    /* In 2D distribution, check makes sense only if proc has entire hedge */
    for (i = 0; i < hg->nEdge; i++)
      if ((hg->hindex[i+1] == hg->hindex[i])) {
        sprintf (str, "Found hyperedge with no vertices: "
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
      if (k == hg->vindex[hg->hvertex[j]+1]) { 
        /* if no match was found then failure, else keep on */
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "Inconsistent hvertex/vedge");
        err = ZOLTAN_WARN;                    
        break;
      }
    }

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}




/****************************************************************************/
void Zoltan_HG_Print(
  ZZ *zz,
  HGraph *hg,
  Partition parts,
  FILE *fp,
  char *str
)
{
/* Routine to print hypergraph weights and edges. Assumes serial execution;
 * put inside Zoltan_Print_Sync_Start/Zoltan_Print_Sync_End for parallel
 * programs. 
 */
int i, j;
int num_vwgt;
int num_ewgt;
float *sum;
char *yo = "Zoltan_HG_Print";

  if (hg == NULL)
    return;

  ZOLTAN_TRACE_ENTER(zz, yo);

  num_vwgt = hg->VtxWeightDim;
  num_ewgt = hg->EdgeWeightDim;

  sum = (float *) ZOLTAN_MALLOC(MAX(num_vwgt, num_ewgt) * sizeof(float));

  fprintf(fp, "%s nVtx=%d nEdge=%d nPins=%d vWgt=%d eWgt=%d\n", 
          str, hg->nVtx, hg->nEdge, hg->nPins, 
          hg->VtxWeightDim, hg->EdgeWeightDim);

  /* Print Vertex Info */
  fprintf(fp, "%s Vertices:  (edges)\n", str);
  for (i = 0; i < hg->nVtx; i++) {
    fprintf(fp, "%d (%d) in part %d:  ", 
            i, VTX_LNO_TO_GNO(hg, i), (parts ? parts[i] : -1));
    fprintf(fp, "(");
    for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
      fprintf(fp, "%d ", hg->vedge[j]);
    fprintf(fp, ")\n");
  }

  if (hg->vwgt != NULL) {
    for (j = 0; j < num_vwgt; j++) sum[j] = 0;
    fprintf(fp, "%s Vertices: [weights])\n", str);
    for (i = 0; i < hg->nVtx; i++) {
      fprintf(fp, "%d (%d):  [", i, VTX_LNO_TO_GNO(hg, i));
      for (j = 0; j < num_vwgt; j++) {
        fprintf(fp, "%f ", hg->vwgt[i*num_vwgt + j]);
        sum[j] += hg->vwgt[i*num_vwgt + j];
      }
      fprintf(fp, "])\n");
    }
    fprintf(fp, "Total vertex weight = [");
    for (j = 0; j < num_vwgt; j++) fprintf(fp, "%f  ", sum[j]);
    fprintf(fp, "]\n");
  }

  /* Print Hyperedge Info */
  fprintf(fp, "%s Hyperedges:  (vertices)\n", str);
  for (i = 0; i < hg->nEdge; i++) {
    fprintf(fp, "%d (%d):  ", i, EDGE_LNO_TO_GNO(hg, i));
    fprintf(fp, "(");
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
      fprintf(fp, "%d ", hg->hvertex[j]);
    fprintf(fp, ")\n");
  }

  if (hg->ewgt != NULL) {
    for (j = 0; j < num_ewgt; j++) sum[j] = 0;
    fprintf(fp, "%s Hyperedge Weights:  [weights]\n", str);
    for (i = 0; i < hg->nEdge; i++) {
      fprintf(fp, "%d (%d):  ", i, EDGE_LNO_TO_GNO(hg, i));
      fprintf(fp, "[");
      for (j = 0; j < num_ewgt; j++) {
        fprintf(fp, "%f ", hg->ewgt[i*num_ewgt + j]);
        sum[j] += hg->ewgt[i*num_ewgt + j];
      }
      fprintf(fp, "])\n");
    }
    fprintf(fp, "Total hyperedge weight = [");
    for (j = 0; j < num_ewgt; j++) fprintf(fp, "%f  ", sum[j]);
    fprintf(fp, "]\n");
  }

  ZOLTAN_FREE(&sum);
  ZOLTAN_TRACE_EXIT(zz, yo);
}

/****************************************************************************/

int Zoltan_HG_HGraph_to_Graph(
  ZZ *zz,
  HGraph *hg,
  Graph *g
)
{
int i, j, k, e, roughly_e, pins, *_neigh=NULL, *degrees, current_degree,
 vertex1, vertex2, start, end, empty;
float *w, *_ewgt=NULL, weight;
int err = ZOLTAN_OK;
char  *yo = "Zoltan_HG_HGraph_to_Graph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_HG_Graph_Init(g);
  g->info = hg->info;
  g->nVtx = hg->nVtx;
  g->nDim = hg->nDim;
  g->VtxWeightDim = hg->VtxWeightDim;
  g->EdgeWeightDim = hg->EdgeWeightDim;
  g->vtxdist = NULL;  /* KDD May need to update vtxdist for use in parallel */

   /* copy coordinates */
   /* not necessary at the moment, because matching does not use them */

  /* copy vertex weights */
  if (hg->vwgt) {
     if (!(g->vwgt = (float*) ZOLTAN_MALLOC (g->nVtx * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
        }
     for (i=0; i<g->nVtx; i++)
        g->vwgt[i] = hg->vwgt[i];
     }

  /* calculate the roughly |E| and each degree */
  roughly_e = 0;
  if (!(degrees = (int*) ZOLTAN_CALLOC (hg->nVtx, sizeof(int))))
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto End;
     }
  for (i = 0; i < hg->nEdge; i++) {
     /* if (hg->hindex[i+1] - hg->hindex[i] <= 10) */
     pins = hg->hindex[i+1] - hg->hindex[i];
     roughly_e += (pins * (pins-1));
     pins--;
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        degrees[hg->hvertex[j]] += pins;
     }

  /* Calculate the initial nindex */
  if (!(g->nindex = (int*) ZOLTAN_MALLOC ((hg->nVtx + 1) * sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto End;
     }
  e = 0;
  for (i = 0; i < hg->nVtx; i++) {
     current_degree = degrees[i];
     degrees[i]     = g->nindex[i] = e;
     e += current_degree;
     }
  g->nindex[hg->nVtx] = e;

  /* Calculate the initial neigh and ewgt */
  if (!(_neigh = (int*)   ZOLTAN_MALLOC (roughly_e * sizeof(int)))
   || !(_ewgt  = (float*) ZOLTAN_MALLOC (roughly_e * sizeof(float))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &_neigh, &_ewgt);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
      }
  for (i = 0; i < hg->nEdge; i++) {
     /* OLD: if (hg->hindex[i+1]-hg->hindex[i] <= 10) */
     pins = hg->hindex[i+1] - hg->hindex[i];
     /* OLD: weight = 1.0 / (pins-1); */
     weight = 2.0 / ((pins-1) * pins);
     if (hg->ewgt)
        weight *= hg->ewgt[i];
     for (j  = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
        vertex1 = hg->hvertex[j];
        for (k = j + 1; k < hg->hindex[i+1]; k++) {
           vertex2                   = hg->hvertex[k];
           _neigh[degrees[vertex1]]  = vertex2;
           _neigh[degrees[vertex2]]  = vertex1;
           _ewgt[degrees[vertex1]++] = _ewgt[degrees[vertex2]++] = weight;
           }
        }
     }
  ZOLTAN_FREE ((void**) &degrees);

  /* Compact identical incident edges and their weight */
  if (!(w = (float*) ZOLTAN_CALLOC (hg->nVtx, sizeof(float)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto End;
     }
  empty = 0;
  for (i = 0; i < hg->nVtx; i++)
     {
     start = g->nindex[i];
     end   = g->nindex[i+1];
     g->nindex[i] -= empty;
     for (j = start; j < end; j++)
        if (_ewgt[j] > 0.0)
           w[_neigh[j]] += _ewgt[j];
     for (j = start; j < end; j++)
        if (_ewgt[j] > 0.0  &&  w[_neigh[j]] > 0.0)
           {
           _neigh[j-empty] = _neigh[j];
           _ewgt[j-empty]  = w[_neigh[j]];
           w[_neigh[j]]    = 0.0;
           }
        else
           empty++;
  }
  g->nindex[hg->nVtx] -= empty;
  ZOLTAN_FREE ((void**) &w);

  /* Reallocate to the exact size */
  g->nEdge = g->nindex[g->nVtx];

  if (!(g->neigh = (int*) ZOLTAN_MALLOC (g->nindex[hg->nVtx] * sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto End;
     }
  memcpy(g->neigh, _neigh, g->nEdge * sizeof(int));
  ZOLTAN_FREE ((void**) &_neigh);

  if (!(g->ewgt = (float*)ZOLTAN_MALLOC(g->nindex[hg->nVtx] * sizeof(float)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto End;
     }
  memcpy(g->ewgt, _ewgt, g->nEdge * sizeof(float));
  ZOLTAN_FREE ((void **) &_ewgt);

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}

/****************************************************************************/
/*
 *  Converts a graph g into a hypergraph hg.
 *  One hyperedge is created for each vertex of g.
 *  Hyperedge i consists of vertex i + all edge neighbors of i in graph g. */

int Zoltan_HG_Graph_to_HGraph(
  ZZ *zz,
  Graph *g,        /* Input graph */
  HGraph *hg       /* Ouput hypergraph */
)
{
int i, j;
int *hindex = NULL, *hvertex = NULL;  /* temporary array pointers */
int cnt;
int err = ZOLTAN_OK;
char *yo = "Zoltan_HG_Graph_to_HGraph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_HG_HGraph_Init(hg);
  hg->info  = g->info;
  hg->nVtx  = g->nVtx;
  hg->nEdge = g->nVtx;
  hg->nDim  = g->nDim;
  hg->nPins = g->nVtx + g->nEdge;
  hg->VtxWeightDim = g->VtxWeightDim;
  hg->EdgeWeightDim   = g->EdgeWeightDim;

  /* Copy vtxdist array; it is the same for hypergraph as graph. */
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
     cnt = hg->nVtx * hg->VtxWeightDim;
     if (!(hg->vwgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float)))) {
        err = ZOLTAN_MEMERR;
        goto End;
        }
     memcpy(hg->vwgt, g->vwgt, cnt * sizeof(float));
     }

  /* Allocate memory for hypergraph edges */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge arrays");
  if (hg->nEdge > 0 && hg->nPins > 0)
    if (!(hindex  = hg->hindex =(int*)ZOLTAN_MALLOC((hg->nEdge+1)*sizeof(int)))
     || !(hvertex = hg->hvertex=(int*)ZOLTAN_MALLOC(hg->nPins*sizeof(int))) ) {
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
  if (hg->nPins != cnt) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Sanity check failed: nPins != cnt.");
     err = ZOLTAN_FATAL;
     goto End;
     }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Creating mirror");
  err = Zoltan_HG_Create_Mirror(zz, hg);
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

void Zoltan_HG_HGraph_Print(
  ZZ *zz,          /* the Zoltan data structure */
  ZHG *zoltan_hg,
  HGraph *hg,
  Partition parts,
  FILE *fp
)
{
/* Printing routine. Can be used to print a Zoltan_HGraph or just an HGraph.
 * Set zoltan_hg to NULL if want to print only an HGraph.
 * Lots of output; synchronized across processors, so is a bottleneck.
 */
  int i;
  int num_gid = zz->Num_GID;
  int num_lid = zz->Num_LID;
  char *yo = "Zoltan_PHG_HGraph_Print";

  if (zoltan_hg != NULL  &&  hg != &zoltan_hg->HG) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Input hg != Zoltan HG");
    return;
  }

  Zoltan_Print_Sync_Start (zz->Communicator, 1);

  /* Print Vertex Info */
  fprintf (fp, "%s Proc %d\n", yo, zz->Proc);
  fprintf (fp, "Vertices (GID, LID, index)\n");
  for (i = 0; i < zoltan_hg->nObj; i++) {
    fprintf(fp, "(");
    ZOLTAN_PRINT_GID(zz, &zoltan_hg->GIDs[i * num_gid]);
    fprintf(fp, ", ");
    ZOLTAN_PRINT_LID(zz, &zoltan_hg->LIDs[i * num_lid]);
    fprintf(fp, ", %d)\n", i);
  }
  Zoltan_HG_Print(zz, hg, parts, fp, "Build");
  Zoltan_Print_Sync_End(zz->Communicator, 1);
}

