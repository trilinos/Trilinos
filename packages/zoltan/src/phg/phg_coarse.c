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

#include "zz_sort.h"    
#include "phg.h"
#include "zoltan_comm.h"



#define PLAN_TAG 32010      /* tag for comm plan */ 



/*
  #define _DEBUG1
  #define _DEBUG2
  #define _DEBUG3  
*/

#ifdef _DEBUG
#define _DEBUG1
#endif


    
/* UVC:
   Following quicksort routines are modified from
   Numerical Recipes Software
   these versions of quicksort seems to perform
   better than the ones that exist in Zoltan
*/
    
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define M             7
#define NSTACK        50


static void uqsorti(int n, int *arr)
{
    int         i, ir=n, j, k, l=1;
    int         jstack=0, istack[NSTACK];
    int         a, temp;
    
    --arr;
    for (;;) {
        if (ir-l < M) {
            for (j=l+1;j<=ir;j++) {
                a=arr[j];
                for (i=j-1;i>=1;i--) {
                    if (arr[i] <= a) 
                        break;
                    arr[i+1] = arr[i];
                }
                arr[i+1]=a;
            }
            if (jstack == 0) 
                break;
            ir=istack[jstack--];
            l=istack[jstack--];
        } else {
            k=(l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            if (arr[l+1] > arr[ir]) 
                SWAP(arr[l+1], arr[ir]);
            if (arr[l] > arr[ir]) 
                SWAP(arr[l], arr[ir]);
            if (arr[l+1] > arr[l]) 
                SWAP(arr[l+1], arr[l]);
            i=l+1;
            j=ir;
            a=arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
            }
            arr[l]=arr[j];
            arr[j]=a;
            jstack += 2;
            if (jstack > NSTACK) 
                errexit("uqsort: NSTACK too small in sort.");
            if (ir-i+1 >= j-l) {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            } else {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
}


#undef M
#undef NSTACK
#undef SWAP

static unsigned int hashValue(HGraph *hg, int n, int *ar)
{
    unsigned int l=0;
    int i;

    /* Chris Torek's hash function */
    for (i=0; i<n; ++i) {
        l *= 33;
        l += (unsigned int) VTX_LNO_TO_GNO(hg, ar[i]);
    }
    return l;
}

/*
  identical operator:
  a[i]  b[i]     res[i]
  -1    -1       -1  : -1 means no edge in that proc; identical to everything :)
  -1    y        -1==a[y] ? y : 0   --> UVC note that this is same as x < y
  x     -1       -1==b[x] ? x : 0   --> UVC note that this is same as y < x
  0     y        0   : 0 means not identical to anyone in this proc; hence not identical anyone in all
  x     0        0   
  x     x        x   : they are identical to same net
  x <   y       x==a[y] ? y : 0
  x >   y       y==b[x] ? x : 0
*/

static int *idenOperandBuf=NULL;

static void identicalOperator(void *va, void *vb, int *len, MPI_Datatype *dt)
{
    int *a=(int *)va, *b=(int *)vb; 
    int i, *x=a, *y=b;

    memcpy(idenOperandBuf, b, sizeof(int)*(*len));
    b = idenOperandBuf;
    --a; --b; /* net ids in the a and b are base-1 numbers */
    for (i=0; i < *len; ++i, ++x, ++y) {
        if (*x == -1 && *y == -1)
            ; /* no op *y = *y */
        else if (*x==0 || *y == 0)
            *y = 0;
        else if (*x==*y)
            ; /* no op */
        else if (*x < *y)
            *y = (*x==a[*y]) ? *y : 0;
        else /* *x > *y */ 
            *y = (*y==b[*x]) ? *x : 0;
    }
}


/* Procedure to coarsen a hypergraph based on a matching. All vertices of one
   match are clustered to a single vertex. Currently, we allow more
   than two vertices to be merged/grouped locally on a proc, but
   allow only pairs of two to be merged between processors.
   All hyperedges are kept; identical hyperedges are not collapsed. 
   The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
int Zoltan_PHG_Coarsening
( ZZ     *zz,         /* the Zoltan data structure */
  HGraph *hg,         /* information about hypergraph, weights, etc. */
  int    *match,      /* Matching, Packing or Grouping array */
  HGraph *c_hg,       /* output: the coarse hypergraph */
  int    *LevelMap,
  int    *LevelCnt,
  int    *LevelSndCnt,
  int   **LevelData,  /* information to reverse coarsenings later */
  struct Zoltan_Comm_Obj  **comm_plan,
  PHGPartParams *hgp
    )   
{
  char     *yo = "Zoltan_PHG_Coarsening";
  PHGComm  *hgc = hg->comm;
  int   ierr=ZOLTAN_OK, i, j, count, size, me=hgc->myProc_x, idx, ni;
  int   *vmark=NULL, *listgno=NULL, *listlno=NULL, *listproc=NULL, *msg_size=NULL, *ip;
  int   *ahindex=NULL, *ahvertex=NULL, *hsize=NULL, *hlsize=NULL, *ids=NULL, *iden;
  int   *extmatchsendcnt=NULL, extmatchrecvcnt=0, *extmatchrecvcounts=NULL;
  float *c_ewgt=NULL;
  char  *buffer=NULL, *rbuffer=NULL;
  unsigned int           *hash=NULL, *lhash=NULL;
  struct Zoltan_Comm_Obj *plan=NULL;
  MPI_Op                 idenOp;

  static int timer_merge=-1, timer_shuffle=-1, timer_remove=-1, timer_theend=-1;
  int time_details = (hgp->use_timers > 3);

  if (time_details) {
    if (timer_merge < 0)
      timer_merge = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Merge");
    if (timer_shuffle < 0)
      timer_shuffle = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Shuffle");
    if (timer_remove < 0)
      timer_remove = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Remove");
    if (timer_theend < 0)
      timer_theend = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Finish");
  }

  if (time_details)
    ZOLTAN_TIMER_START(zz->ZTime, timer_merge, hgc->Communicator);

#ifdef _DEBUG1
  int    totiden, totsize1;  
  double t_all, t_coarse, t_redhash, t_redsize, t_userredop, t_suffle, t_sort, t_iden, t_shrink, t_mirror, t_cur;

#ifndef _DEBUG2
  if (!hgc->myProc)
#endif
  uprintf(hgc, "In Coarsening....\n");


  MPI_Barrier(hgc->Communicator);
  t_all=MPI_Wtime();
  t_coarse = -t_all;
#endif

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  c_hg->comm    = hg->comm;         /* set communicators */
  c_hg->info    = hg->info + 1;     /* for debugging */
  c_hg->coor    = NULL;             /* currently we don't use coordinates */
  c_hg->nDim    = hg->nDim;    
  c_hg->vmap    = NULL;             /* only needed by rec-bisec */
  c_hg->redl    = hg->redl;  /* to stop coarsening near desired count */
  c_hg->VtxWeightDim  = hg->VtxWeightDim;
  c_hg->bisec_split = hg->bisec_split;
  c_hg->fixed = (hg->fixed) ? (int*)ZOLTAN_MALLOC(hg->nVtx * sizeof(int)) : NULL;
  
if (c_hg->fixed)  
  for (i = 0; i < hg->nVtx; i++)
    c_hg->fixed[i] = -2;  
  
  
  /* (over) estimate number of external matches that we need to send data to */
  count = 0;
  if (hgp->match_array_type==0) { /* old style */
      for (i = 0; i < hg->nVtx; i++)
          if (match[i] < 0)
              ++count;
      extmatchrecvcnt = count;
  } else { /* new style match array */
      if (!(extmatchsendcnt = (int *) ZOLTAN_CALLOC(hgc->nProc_x, sizeof(int)))
          || !(extmatchrecvcounts = (int *) ZOLTAN_MALLOC(hgc->nProc_x * sizeof(int))))
          MEMORY_ERROR;
      for (i = 0; i < hgc->nProc_x; ++i)
          extmatchrecvcounts[i] = 1;
      for (i = 0; i < hg->nVtx; i++) {
          int px=VTX_TO_PROC_X(hg, match[i]);

          if (px!=me) {
              ++extmatchsendcnt[px];
              ++count;              
          }
      }

      MPI_Reduce_scatter(extmatchsendcnt, &extmatchrecvcnt, extmatchrecvcounts, MPI_INT, MPI_SUM, hgc->row_comm);
      Zoltan_Multifree (__FILE__, __LINE__, 2, &extmatchsendcnt, &extmatchrecvcounts);
  }
 
  if (hg->nVtx > 0 && !(vmark = (int*) ZOLTAN_CALLOC(MAX(hg->nEdge, hg->nVtx),  sizeof(int))))
      MEMORY_ERROR;

  size = MAX(count, hg->nEdge);
  if (size && (
      !(listproc  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))
   || !(msg_size  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))))
      MEMORY_ERROR;

  if (count > 0 && (
      !(listgno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
      || !(listlno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))))
      MEMORY_ERROR;
  if (extmatchrecvcnt && 
      !(*LevelData= (int*) ZOLTAN_MALLOC (extmatchrecvcnt * sizeof(int) * 2)))
      MEMORY_ERROR;

  if (hgp->match_array_type==0) { /* old style */
      /* Assume all rows in a column have the entire (column's) matching info */
      /* Calculate the number of resulting coarse verticeaas. */
      c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
      me = hgc->myProc_x;             /* short name, convenience variable */
      size  = 0;                      /* size (in ints) to communicate */
      count = 0;                      /* number of vertices to communicate */
      for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
          if (match[i] < 0)  {               /* external processor match */
              int gx = -match[i]-1, proc = VTX_TO_PROC_X(hg, gx);
              
              /* rule to determine "ownership" of coarsened vertices between procs */
              proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(proc, me) : MAX(proc, me);
              
              /* prepare to send data to owner */
              if (proc != me)   {             /* another processor owns this vertex */
                  LevelMap[i] = -gx - 1;
                  size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
                  listgno[count]   = gx;                    /* listgno of vtx's to send */
                  listproc[count]  = proc;                  /* proc to send to */
                  listlno[count++] = i;                     /* lno of my match to gno */
              } else   { /* myProc owns the matching across processors */ 
                  LevelMap[i] = c_hg->nVtx++;        /* next available coarse vertex */
              }
          } else if (!vmark[i]) {     /* local matching, packing and groupings */
              int v = i;
              int fixed = -1; 
              while (!vmark[v])  {
                  if (hgp->UseFixedVtx && hg->fixed[v] >= 0)
                      fixed = hg->fixed[v];
                  LevelMap[v] = c_hg->nVtx;         /* next available coarse vertex */      
                  vmark[v] = 1;  /* flag this as done already */
                  v = match[v];  
              }
              if (hgp->UseFixedVtx)
                  c_hg->fixed[c_hg->nVtx] = fixed;      
              ++c_hg->nVtx;
          }
      }
      *LevelSndCnt = count;
  } else {  /* new style match array */
    c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
    me = hgc->myProc_x;             /* short name, convenience variable */
    size  = 0;                      /* size (in ints) to communicate */
    count = 0;                      /* number of vertices to communicate */
    for (i = 0; i < hg->nVtx; ++i)
        if (match[i] == VTX_LNO_TO_GNO(hg, i)) {
            LevelMap[i] = c_hg->nVtx;
            if (c_hg->fixed)
                c_hg->fixed[c_hg->nVtx] = hg->fixed[i];
/*            uprintf(hgc, "match[%d (gno=%d)] = %d   new vtxno=%d\n", i, VTX_LNO_TO_GNO(hg, i), match[i], c_hg->nVtx);*/
            ++c_hg->nVtx;
        }
    
    for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
      if (match[i] != VTX_LNO_TO_GNO(hg, i))  {
              int gx = match[i], proc = VTX_TO_PROC_X(hg, gx);
                            
              if (proc != me)   {             /* owner is external */
                  LevelMap[i] = -gx - 1;     /* prepare to send data to owner */
                  size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
                  listgno[count]   = gx;                    /* listgno of vtx's to send */
                  listproc[count]  = proc;                  /* proc to send to */
                  listlno[count++] = i;                     /* lno of my match to gno */
/*                  uprintf(hgc, "EXTMAT:  match[%d (gno=%d)] = %d \n", i, VTX_LNO_TO_GNO(hg, i), match[i]);*/
      } else   { /* owner is local */
                  LevelMap[i] = LevelMap[VTX_GNO_TO_LNO(hg, gx)];
/*                  uprintf(hgc, "LOCMAT:  match[%d (gno=%d)] = %d   new vtxno=%d\n", i, VTX_LNO_TO_GNO(hg, i), match[i], LevelMap[i]);*/
              }
      }
    }
    *LevelSndCnt = count;
/*      errexit("this type of coarsening is not implemented yet"); */
  }
  
  /* size and allocate the send buffer */
  { 
  int header = (hgp->UseFixedVtx ? 4 : 3);
  size += ((header + hg->VtxWeightDim) * count);
  if (size > 0 && !(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))
      MEMORY_ERROR;
  } 
  
  /* Message is list of <gno, vweights, gno's edge count, list of edge lno's> */
  /* We pack ints and floats into same buffer */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
    int lno=listlno[i], sz=hg->vindex[lno+1] - hg->vindex[lno];
    int *ip_old = ip;
    *ip++ = listlno[i];             /* source lno */
    *ip++ = listgno[i];             /* destination vertex gno */
    if (hg->fixed)
      *ip++ = hg->fixed[lno];
      
    /* UVC Assumes a float is no larger than an int which true in almost all
       platforms except prehistoric 16-bit platforms. */
    memcpy(ip, &hg->vwgt[lno*hg->VtxWeightDim], sizeof(float)*hg->VtxWeightDim);
    ip += hg->VtxWeightDim;
    
    *ip++ = sz;                    /* count */
    memcpy(ip, &hg->vedge[hg->vindex[lno]], sizeof(int)*sz);
    ip +=  sz;
    msg_size[i] = ip - ip_old;     /* save mesg size in #ints */
  }    

  /* Create comm plan. */
  Zoltan_Comm_Create(comm_plan, count, listproc, hgc->row_comm, PLAN_TAG, 
                      &size); /* we'll ignore the size because of resize*/
  
  /* call Comm_Resize since we have variable-size messages */
  Zoltan_Comm_Resize(*comm_plan, msg_size, PLAN_TAG+1, &size); 

  /* Allocate receive buffer. */
  /* size is the size of the received data, measured in #ints */
  if (size  && (
       !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int)))
    || !(ahvertex = (int*) ZOLTAN_MALLOC (size * sizeof(int)))))
      MEMORY_ERROR;
  if (!(ahindex = (int *) ZOLTAN_CALLOC(hg->nEdge+1, sizeof(int))))
      MEMORY_ERROR;
  if (hg->nEdge && (
       !(hlsize  = (int *) ZOLTAN_MALLOC(hg->nEdge*sizeof(int)))
    || !(lhash   = (unsigned int *) ZOLTAN_MALLOC(hg->nEdge*sizeof(unsigned int)))
    || !(hash    = (unsigned int *) ZOLTAN_MALLOC(hg->nEdge*sizeof(unsigned int)))
          ))
      MEMORY_ERROR;

  /* Comm_Do sends personalized messages of variable sizes */
  Zoltan_Comm_Do(*comm_plan, PLAN_TAG+2, buffer, sizeof(int), rbuffer);

  /* Allocate vertex weight array for coarse hgraph */
  if (c_hg->nVtx > 0 && c_hg->VtxWeightDim > 0 &&
      !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx * c_hg->VtxWeightDim,
                                             sizeof(float))))
      MEMORY_ERROR;
  for (i=0; i < hg->nVtx; ++i) {
      int ni=LevelMap[i];
      if (ni>=0)
          for (j=0; j<hg->VtxWeightDim; ++j)
              c_hg->vwgt[ni*hg->VtxWeightDim+j] += hg->vwgt[i*hg->VtxWeightDim+j];
  }
      
  /* index all received data for rapid lookup */
  *LevelCnt   = 0;
  for (ip = (int*) rbuffer, i = 0; i < size; )  {
    int j, sz, source_lno, lno;
    float *pw;

    source_lno              = ip[i++];
    lno = VTX_GNO_TO_LNO (hg, ip[i++]);
       
    if (hgp->UseFixedVtx)  {
      int fixed = ip[i++];
      c_hg->fixed [LevelMap[lno]] = (fixed >= 0) ? fixed : hg->fixed[lno];
    } 
     
    pw=(float*) &ip[i];

    (*LevelData)[(*LevelCnt)++] = source_lno;
    (*LevelData)[(*LevelCnt)++] = lno;              /* to lookup in part[] */

    lno = LevelMap[lno];
    for (j=0; j<hg->VtxWeightDim; ++j)
        c_hg->vwgt[lno*hg->VtxWeightDim+j] += pw[j];
    i += hg->VtxWeightDim;    /* skip vertex weights */
    sz = ip[i++];             /* get sz to a diff var, skip size*/
    for (j=0; j<sz; ++j)      /* count extra vertices per edge */
        ++ahindex[ip[i+j]];
    i += sz;
  }

  for (i=0; i<hg->nEdge; ++i) /* prefix sum over ahindex */
    ahindex[i+1] += ahindex[i];
  /* now prepare ahvertex */
  for (ip = (int*) rbuffer, i = 0; i < size; )  {
    int j, sz, lno=VTX_GNO_TO_LNO (hg, ip[i+1]);

    i += (hgp->UseFixedVtx ? 3 : 2) + hg->VtxWeightDim;  /* skip header + vtx weights */
    sz = ip[i++];             /* get sz to a diff var, skip size*/
    for (j=0; j<sz; ++j)      /* count extra vertices per edge */
        ahvertex[--ahindex[ip[i+j]]] = lno;
    i += sz;
  }
  Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);
  ip = NULL;
  
  c_hg->nPins = hg->nPins + ahindex[hg->nEdge]; /* safe over estimate of nPins */
  c_hg->nEdge = hg->nEdge;

      
  /* Allocate edge weight and index array for coarse hgraph */
  c_hg->EdgeWeightDim = (hg->EdgeWeightDim > 0) ? hg->EdgeWeightDim : 1;
  if (c_hg->nEdge) {
      if (!(c_hg->ewgt =(float*)ZOLTAN_MALLOC(c_hg->nEdge*c_hg->EdgeWeightDim*sizeof(float))))
          MEMORY_ERROR;
      if (hg->EdgeWeightDim)
          for (i = 0; i < hg->nEdge * hg->EdgeWeightDim; ++i)
              c_hg->ewgt[i] = hg->ewgt[i];
      else
          for (i = 0; i < c_hg->nEdge; ++i)
              c_hg->ewgt[i] = 1.0;
  }
  if (!(c_hg->hindex = (int*)ZOLTAN_MALLOC((c_hg->nEdge+1)*sizeof(int))))
      MEMORY_ERROR;

  if (c_hg->nPins>0 && !(c_hg->hvertex = (int*)ZOLTAN_MALLOC (c_hg->nPins*sizeof(int))))
      MEMORY_ERROR;

  memset(vmark, 0xff, sizeof(int)*c_hg->nVtx);
  idx = 0;
  for (i=0; i < hg->nEdge; ++i) { /* loop over edges */
      int sidx=idx;
      c_hg->hindex[i] = idx;
      /* first go over local vertices */
      for (j=hg->hindex[i]; j<hg->hindex[i+1]; ++j) {
          int nvno=LevelMap[hg->hvertex[j]];
          if (nvno>=0 && vmark[nvno]!=i) {
              c_hg->hvertex[idx++] = nvno;
              vmark[nvno] = i;
          }
      }
      /* now go over the received vertices */
      for (j=ahindex[i]; j<ahindex[i+1]; ++j) {
          int nvno=LevelMap[ahvertex[j]];
          if (nvno>=0 && vmark[nvno]!=i) {
              c_hg->hvertex[idx++] = nvno;
              vmark[nvno] = i;
          }
      }          
      /* in qsort start and end indices are inclusive */
      /*    Zoltan_quicksort_list_inc_int(&c_hg->hvertex[sidx], 0, idx-sidx-1); */
      /* UVC my qsort code is a little bit faster :) */
      uqsorti(idx-sidx, &c_hg->hvertex[sidx]);
  }
  c_hg->hindex[hg->nEdge] = c_hg->nPins = idx;

  if (time_details) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_merge, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer_shuffle, hgc->Communicator);
  }
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_coarse += t_cur;
  t_redhash = -t_cur;
#endif
  for (i=0; i < c_hg->nEdge; ++i) { /* compute size and hashvalue */
      hlsize[i] = c_hg->hindex[i+1]-c_hg->hindex[i];
      lhash[i] = hashValue(hg, hlsize[i], &c_hg->hvertex[c_hg->hindex[i]]);
  }

  /* UVC TODO to compute global hash; right now we'll use SUM (UVC TODO: try:bitwise xor);
     we need to check if this is good, if not we need to find a better way */
  if (c_hg->nEdge) 
      MPI_Allreduce(lhash, hash, c_hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);

  Zoltan_Multifree(__FILE__, __LINE__, 3, &vmark, &ahvertex, &ahindex);
  
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur =  MPI_Wtime();
  t_redhash += t_cur;
  t_suffle = -t_cur;
#endif

  for (i=0; i < c_hg->nEdge; ++i)  /* decide where to send */
      listproc[i] = (int) (hash[i] % hgc->nProc_y);
  
  Zoltan_Comm_Create(&plan, c_hg->nEdge, listproc, hgc->col_comm, PLAN_TAG+10, 
                     &size);

  /* send hash values */
  if (size > hg->nEdge) {
      Zoltan_Multifree(__FILE__, __LINE__, 2, &lhash, &listproc);
      if (!(lhash=(unsigned int *) ZOLTAN_MALLOC(size * sizeof(unsigned int)))
          || !(listproc=(int *) ZOLTAN_MALLOC(size * sizeof(int))))
          MEMORY_ERROR;
  }
  Zoltan_Comm_Do(plan, PLAN_TAG+11, (char *) hash, sizeof(unsigned int), (char *) lhash);
  ZOLTAN_FREE(&hash); /* we don't need it anymore */

  /* now local sizes */
  if (!(ahindex = (int *)  ZOLTAN_MALLOC((1+size) * sizeof(int))))
    MEMORY_ERROR;
  if (size && (
       !(ip      = (int *)  ZOLTAN_MALLOC(size * sizeof(int)))
    || !(hsize =   (int *)  ZOLTAN_MALLOC(size * sizeof(int)))       
    || !(c_ewgt  = (float *)ZOLTAN_MALLOC(size * sizeof(float)*c_hg->EdgeWeightDim)))) 
      MEMORY_ERROR;

  Zoltan_Comm_Do(plan, PLAN_TAG+12, (char *) hlsize, sizeof(int), (char *) ip);
  /* now ewgt  */
  Zoltan_Comm_Do(plan, PLAN_TAG+12, (char *) c_hg->ewgt, c_hg->EdgeWeightDim*sizeof(float), (char *) c_ewgt);

  /* now vertices of hyperedges */

  Zoltan_Comm_Resize(plan, hlsize, PLAN_TAG+13, &idx);
  if (idx && !(ahvertex = (int *) ZOLTAN_MALLOC(idx * sizeof(int))))
      MEMORY_ERROR;

  Zoltan_Comm_Do(plan, PLAN_TAG+14, (char *) c_hg->hvertex, sizeof(int), (char *) ahvertex);
  Zoltan_Comm_Destroy (&plan);
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_suffle += t_cur;
  t_redsize = -t_cur;
#endif

  ZOLTAN_FREE(&hlsize);    hlsize=ip;

  if (size) 
      MPI_Allreduce(hlsize, hsize, size, MPI_INT, MPI_SUM, hgc->row_comm);
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_redsize += t_cur;
  t_sort = -t_cur;
#endif

  if (time_details) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_shuffle, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer_remove, hgc->Communicator);
  }

  /* in order to find identical nets; we're going to sort hash values and compare them */
  if (size && !(ids = (int *) ZOLTAN_MALLOC(size * sizeof(int))))
      MEMORY_ERROR;
  ahindex[0] = 0;
  for (i=0; i<size; ++i) {
      ahindex[1+i] = ahindex[i] + hlsize[i];
      ids[i] = i;
  }

  /* lhash is actually global hash */
  Zoltan_quicksort_pointer_inc_int_int (ids, (int *)lhash, hsize, 0, size-1); 
/*  uqsort_ptr_uint_int(size, lhash, hsize, ids); */

#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_sort += t_cur;
  t_iden = -t_cur;
#endif

  
  iden = listproc; /* just better variable name */
  for (j=0; j<size; ++j)
      iden[j] = (hlsize[j]) ? 0 : -1; /* if no local pins iden is -1 */
#ifdef _DEBUG1  
  count = idx = me = 0;
  for (j=0; j<size; ++j)
      if (iden[j]==-1)
          ++me;
#endif
  for (j=0; j<size; ++j) {
      int n1=ids[j];
      if (!iden[n1]) {
          int last=1+n1, minv=last;
              
          for (i = j+1; i<size && lhash[n1] == lhash[ids[i]] && hsize[n1]==hsize[ids[i]]; ++i) {
              int n2=ids[i];
              
#ifdef _DEBUG1
              ++idx;
#endif
              if (!iden[n2] && hlsize[n1]==hlsize[n2]
                  && !memcmp(&ahvertex[ahindex[n1]], &ahvertex[ahindex[n2]],
                             sizeof(int)*hlsize[n1])) {
                  iden[n2] = last; /* n2 is potentially identical to n1 */
                  last = 1+n2;
                  minv = (last<minv) ? last : minv;
#ifdef _DEBUG1
                  ++count;
#endif
              }                  
          }
          /* iden[last] is a link list (in array) now make the
             all identical nets to point the same net with the smallest id;
             it will be needed in identicalOperator; we'll zero(clear) the
             original net (net with the smallest id) after this loop */
          while (last!=1+n1) {
              int prev=iden[last-1];
              iden[last-1] = minv;
              last = prev;
          }
          iden[n1] = minv;
      }
  }

  for (i=0; i<size; ++i)
      if (iden[i]==1+i) /* original net; clear iden */
          iden[i] = 0;
  ZOLTAN_FREE(&ids); 

#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_iden += t_cur;
  t_userredop = -t_cur;

#ifndef _DEBUG2
  if (!hgc->myProc)
#endif
  uprintf(hgc, "#Loc.Iden= %7d   (Computed= %d)   #Comp.PerNet= %.2lf     ElapT= %.3lf\n", count+me, count, (double) idx / (double) size, MPI_Wtime()-t_all);
  count += me;
#endif
  

  ip = (int *) lhash; 
  if (size) {
      MPI_Op_create(identicalOperator, 1, &idenOp);
      if (!(idenOperandBuf = (int *) ZOLTAN_MALLOC(size * sizeof(int))))
          MEMORY_ERROR;
      MPI_Allreduce(iden, ip, size, MPI_INT, idenOp, hgc->row_comm);
      MPI_Op_free(&idenOp);
      ZOLTAN_FREE(&idenOperandBuf);
  }
#ifdef _DEBUG1  
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_userredop +=  t_cur;
  t_shrink = -t_cur;

  me = idx = 0;
#endif
  
  c_hg->nPins = 0;
  c_hg->nEdge = 0;
  for (i=0; i<size; ++i) {
#ifdef _DEBUG1
      if (ip[i]==-1 && hlsize[i])
          errexit("ip[%d]==-1 but hlsize[%d] = %d", i, i, hlsize[i]);
#endif
      if (ip[i]>0) { /* identical net */
          int n1=ip[i]-1; /* i is identical to n1 */
          /* Removing identical net here:  i == n1 */
          for (j=0; j<c_hg->EdgeWeightDim; ++j)
              c_ewgt[n1*c_hg->EdgeWeightDim + j] += c_ewgt[i*c_hg->EdgeWeightDim + j];
#ifdef _DEBUG1
          ++idx;
          if (n1>i)
              errexit("n1(%d) > i(%d)", n1, i);
          if (hsize[n1]>1 && ip[n1]==0)
              errexit("i=%d is pointing net %d but that net is going to be removed", i, n1);
#endif
          ip[i] = 0; /* ignore */
      } else if (hsize[i]>1) { /*  NOT identical and size not 0/1 */          
          c_hg->nPins += hlsize[i];
          ++c_hg->nEdge;
          ip[i] = 1; /* Don't ignore */
      } else {
          /* ignoring net i with size hsize[i] */
          ip[i] = 0; /* ignore size 0/1 nets*/
#ifdef _DEBUG1
          ++me;
#endif          
      }
  }

#ifdef _DEBUG1

  t_cur = MPI_Wtime()-t_all;
#ifdef _DEBUG2
  uprintf(hgc, "#GlobIden= %7d    SuccessRate= %.1lf%%    ElapT= %.3lf HashOpT=%.3lf (%.1lf%%)\n", idx, 100.0 * idx / (double) count, t_cur, t_userredop, 100.0*t_userredop/t_cur);
#endif
  MPI_Allreduce(&idx, &totiden, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  MPI_Allreduce(&me, &totsize1, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  if (!hgc->myProc)
      uprintf(hgc, "Level %d - Orig #Nets=%d    Iden=%d    Size-0/1=%d CurT=%.3lf user-def redopT=%.3lf\n", hg->info, hg->dist_y[hgc->nProc_y], totiden, totsize1, t_cur, t_userredop);  
#endif

  Zoltan_Multifree(__FILE__, __LINE__, 3, &c_hg->hindex, &c_hg->hvertex, &c_hg->ewgt);

  c_hg->hindex = (int*)ZOLTAN_MALLOC((c_hg->nEdge+1)*sizeof(int));
  if (c_hg->nPins &&
      !(c_hg->hvertex=(int*)ZOLTAN_MALLOC(c_hg->nPins*sizeof(int))))
      MEMORY_ERROR;
  if (c_hg->nEdge &&
      !(c_hg->ewgt=(float*)ZOLTAN_MALLOC(c_hg->nEdge*c_hg->EdgeWeightDim*sizeof(float))))
      MEMORY_ERROR;

#ifdef _DEBUG1  
#ifndef _DEBUG2
  if (!hgc->myProc)
#endif
  uprintf(hgc, "Reconstructing coarsen hygr.... ElapT= %.3lf\n", MPI_Wtime()-t_all);
#endif  
  
  for (idx=ni=i=0; i<size; ++i)
      if (ip[i]) {
          c_hg->hindex[ni] = idx;
          memcpy(&c_hg->ewgt[ni], &c_ewgt[i], c_hg->EdgeWeightDim*sizeof(float));
          memcpy(&c_hg->hvertex[idx], &ahvertex[ahindex[i]], hlsize[i]*sizeof(int));
          ++ni;
          idx += hlsize[i];
      }
  c_hg->hindex[c_hg->nEdge] = idx;

  if (time_details) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_remove, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer_theend, hgc->Communicator);
  }

#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);  
  t_cur = MPI_Wtime();
  t_shrink += t_cur;
  t_mirror = -t_cur;


  if (idx!=c_hg->nPins || ni!=c_hg->nEdge)
      errexit("idx(%d)!=c_hg->nPins(%d) || ni(%d)!=c_hg->nEdge(%d)", idx, c_hg->nPins, ni, c_hg->nEdge);
#endif

  /* We need to compute dist_x, dist_y */
  if (!(c_hg->dist_x = (int *) ZOLTAN_CALLOC((hgc->nProc_x+1), sizeof(int)))
	 || !(c_hg->dist_y = (int *) ZOLTAN_CALLOC((hgc->nProc_y+1), sizeof(int))))
      MEMORY_ERROR;

  MPI_Scan(&c_hg->nVtx, c_hg->dist_x, 1, MPI_INT, MPI_SUM, hgc->row_comm);
  MPI_Allgather(c_hg->dist_x, 1, MPI_INT, &(c_hg->dist_x[1]), 1, MPI_INT, hgc->row_comm);
  c_hg->dist_x[0] = 0;
  
  MPI_Scan(&c_hg->nEdge, c_hg->dist_y, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  MPI_Allgather(c_hg->dist_y, 1, MPI_INT, &(c_hg->dist_y[1]), 1, MPI_INT, hgc->col_comm);
  c_hg->dist_y[0] = 0;  
  
if (c_hg->fixed)
  for (i = 0; i < c_hg->nVtx; i++)
    if (c_hg->fixed[i] == -2)
      printf ("RTHRTH BAD COARSENING for FIXED VERTICES\n"); 

  ierr = Zoltan_HG_Create_Mirror(zz, c_hg);
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);  
  t_mirror += MPI_Wtime();
#endif
 End:
#ifdef _DEBUG1
  t_mirror -= MPI_Wtime();
#endif
  Zoltan_Multifree (__FILE__, __LINE__, 14,
                    &listgno, &listlno, &listproc, &msg_size,
                    &buffer, &rbuffer, &ahindex, &ahvertex, &vmark,
                    &hlsize, &hsize, &lhash, &hash, &c_ewgt 
                    );
#ifdef _DEBUG1

  MPI_Barrier(hgc->Communicator);  
  t_cur = MPI_Wtime();
  t_mirror += t_cur;
  t_all = t_cur-t_all;
#ifndef _DEBUG2
  if (!hgc->myProc)
#endif
  uprintf(hgc, "Terminating Coarsening ... ElapT= %.3lf Coar= %.3lf ( %.1lf%% ) RedHash= %.3lf ( %.1lf%% ) RedSize= %.3lf ( %.1lf%% ) RedIden= %.3lf ( %.1lf%% ) Suffle= %.3lf ( %.1lf%% ) Sort= %.3lf ( %.1lf%% ) Iden= %.3lf ( %.1lf%% ) Shrink= %.3lf ( %.1lf%% )  Mirror= %.3lf ( %.1lf%% ) Rest= %.3lf ( %.1lf%% )\n", t_all, t_coarse, 100.0*t_coarse/t_all, t_redhash, 100.0*t_redhash/t_all, t_redsize, 100.0*t_redsize/t_all, t_userredop, 100.0*t_userredop/t_all, t_suffle, 100.0*t_suffle/t_all, t_sort, 100.0*t_sort/t_all, t_iden, 100.0*t_iden/t_all, t_shrink, 100.0*t_shrink/t_all, t_mirror, 100.0*t_mirror/t_all, t_all-(t_coarse+t_redhash+t_redsize+t_userredop+t_suffle+t_sort+t_iden+t_shrink+t_mirror), 100.0*(t_all-(t_coarse+t_redhash+t_redsize+t_userredop+t_suffle+t_sort+t_iden+t_shrink+t_mirror))/t_all);
#endif
  if (time_details)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_theend, hgc->Communicator);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
