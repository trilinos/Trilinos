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
#include "zz_util_const.h"


#define PLAN_TAG 32010      /* tag for comm plan */ 


/*
  #define _DEBUG1
  #define _DEBUG2
  #define _DEBUG3  
*/

#ifdef _DEBUG
#define _DEBUG1
#endif

#define BITSET(data, elem)       ((data)[(elem) >> 5] |=( 1 << ((elem) & 31)))
#define BITCHECK(data, elem)     ((data)[(elem) >> 5] & (1 << ((elem) & 31)))
    
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
    int         temp;
    int a;
    
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


/* Procedure to coarsen a hypergraph based on a matching. All vertices of one
   match are clustered to a single vertex. Currently, we allow more
   than two vertices to be merged/grouped locally on a proc, but
   allow only pairs of two to be merged between processors.
   If a coarse vertex constitutes a fixed fine vertex; coarse vertex is fixed
   in the same part (note that if a coarse vertex constitues of multiple fixed
   vertices; since they have to be fixed in the "same" side of bisection; the
   coarse vertex is fixed in the "first" fixed part).
   Identical hyperedges are identified and collapsed into single one. 
   The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
int Zoltan_PHG_Coarsening
( ZZ     *zz,         /* the Zoltan data structure */
  HGraph *hg,         /* information about hypergraph, weights, etc. */
  ZOLTAN_GNO_TYPE *match,  /* Matching, Packing or Grouping array */
  HGraph *c_hg,       /* output: the coarse hypergraph */
  ZOLTAN_GNO_TYPE *LevelMap,    /* TODO64: ZOLTAN_GNO_TYPE or not? */
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
  int   *vmark=NULL, *listlno=NULL, *listproc=NULL, *ip=NULL;
  int   *ahindex=NULL, *hlsize=NULL, *hsize=NULL, *ids=NULL, *iden;
  int   *emptynets=NULL, emptynetsize, *idennets=NULL, *allemptynets=NULL, *allidennets=NULL,
#ifdef _DEBUG1
      emptynetcnt, 
#endif
      *idennetdest=NULL, rootRank;
  int   *extmatchsendcnt=NULL, extmatchrecvcnt=0, *extmatchrecvcounts=NULL;
  float *c_ewgt=NULL;
  unsigned int *hash=NULL, *lhash=NULL;
  struct Zoltan_Comm_Obj *plan=NULL;
  ZOLTAN_GNO_TYPE *listgno=NULL;
  char *buffer=NULL, *rbuffer=NULL, *b=NULL, *b_end=NULL;
  double *coordbuf=NULL, *coordrecbuf=NULL, *cb=NULL, *cb_end=NULL;
  int *ahvertex=NULL;
  int *msg_size=NULL;     /* TODO64 - would we ever need ZOLTAN_GNO_TYPE for message sizes? */
  double *coorcount = NULL; /* No. of vertices represented by the corresponding c_hg->coor */
  int gno_size, alt_field_count, int_size, float_size;
  ZOLTAN_GNO_TYPE tmp_gno;
  ZOLTAN_GNO_TYPE *gnoptr;
  int *intptr;
  float *floatptr;
  double *doubleptr;
  MPI_Datatype zoltan_gno_mpi_type;
  struct phg_timer_indices *timer = NULL;
  int time_details;

#ifdef _DEBUG1
  int    totiden, totsize1;  
  double t_all, t_coarse, t_redhash, t_redsize, t_userredop, t_suffle, t_sort, t_iden, t_shrink, t_mirror, t_cur;
#endif

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  timer = Zoltan_PHG_LB_Data_timers(zz);
  time_details = (hgp->use_timers > 3);

  if (time_details) {
    if (timer->comerge < 0)
      timer->comerge = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Merge");
    if (timer->coshuffle < 0)
      timer->coshuffle = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Shuffle");
    if (timer->coremove < 0)
      timer->coremove = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Remove");
    if (timer->cotheend < 0)
      timer->cotheend = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsen_Finish");
  }

  if (time_details)
    ZOLTAN_TIMER_START(zz->ZTime, timer->comerge, hgc->Communicator);

#ifdef _DEBUG1

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
  c_hg->coor    = NULL;           
  c_hg->nDim    = hg->nDim;    
  c_hg->vmap    = NULL;             /* only needed by rec-bisec */
  c_hg->redl    = hg->redl;  /* to stop coarsening near desired count */
  c_hg->VtxWeightDim  = hg->VtxWeightDim;
  c_hg->bisec_split = hg->bisec_split;
  c_hg->fixed_part = (hg->fixed_part) ? (int*)ZOLTAN_MALLOC(hg->nVtx * sizeof(int)) : NULL;
  c_hg->pref_part = (hg->pref_part) ? (int*)ZOLTAN_MALLOC(hg->nVtx * sizeof(int)) : NULL;

  if ( ((hg->nVtx && hg->fixed_part) && !c_hg->fixed_part) ||
       ((hg->nVtx && hg->pref_part) && !c_hg->pref_part) ){
    MEMORY_ERROR;
  }

#ifdef _DEBUG1  
  if (c_hg->fixed_part)  
      for (i = 0; i < hg->nVtx; i++)
          c_hg->fixed_part[i] = -2;
  if (c_hg->pref_part)  
      for (i = 0; i < hg->nVtx; i++)
          c_hg->pref_part[i] = -2;
#endif
  
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

  size = MAX(hgc->nProc_x, MAX(count, hg->nEdge));
  if (size && (
      !(listproc  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))
   || !(msg_size  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))))
      MEMORY_ERROR;

  if (count > 0 && (
      !(listgno   = (ZOLTAN_GNO_TYPE*) ZOLTAN_MALLOC (count * sizeof(ZOLTAN_GNO_TYPE)))
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
              ZOLTAN_GNO_TYPE gx = -match[i]-1, proc = VTX_TO_PROC_X(hg, gx);
              
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
                  LevelMap[i] = (ZOLTAN_GNO_TYPE)c_hg->nVtx++;        /* next available coarse vertex */
              }
          } else if (!vmark[i]) {     /* local matching, packing and groupings */
              int v = i;
              int fixed = -1, pref = -1; 
              while (!vmark[v])  {                           /* TODO64 invalid read if v is ZOLTAN_GNO_TYPE*/
                  if (hgp->UseFixedVtx && hg->fixed_part[v] >= 0)
                      fixed = hg->fixed_part[v];
                  if (hgp->UsePrefPart && hg->pref_part[v] >= 0)
                      pref = hg->pref_part[v];
                  LevelMap[v] = (ZOLTAN_GNO_TYPE)c_hg->nVtx;         /* next available coarse vertex */      
                  vmark[v] = 1;  /* flag this as done already */
                  v = (int)match[v];  
              }
              if (hgp->UseFixedVtx) {
                  c_hg->fixed_part[c_hg->nVtx] = fixed;
              }
              if (hgp->UsePrefPart)
                  c_hg->pref_part[c_hg->nVtx] = pref;      
              ++c_hg->nVtx;
          }
      }
      *LevelSndCnt = count;
  } else {  /* new style match array */
    c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
    me = hgc->myProc_x;             /* short name, convenience variable */
    size  = 0;                      /* size (in ints) to communicate */
    count = 0;                      /* number of vertices to communicate */
    for (i = 0; i < hg->nVtx; ++i){
        if (match[i] == VTX_LNO_TO_GNO(hg, i)) {
            LevelMap[i] = (ZOLTAN_GNO_TYPE)c_hg->nVtx;
#ifdef KDDKDD_DEBUG
if (VTX_LNO_TO_GNO(hg, i) == 35 || VTX_LNO_TO_GNO(hg, i) == 65 || VTX_LNO_TO_GNO(hg, i) == 66) printf("%d MATCH %d (%f %f %f) to %d; LevelMap = %d\n", zz->Proc, VTX_LNO_TO_GNO(hg, i), hg->coor[i*3], hg->coor[i*3+1], hg->coor[i*3+2], match[i], LevelMap[i]);
#endif
            if (c_hg->fixed_part)
                c_hg->fixed_part[c_hg->nVtx] = hg->fixed_part[i];
            if (c_hg->pref_part)
                c_hg->pref_part[c_hg->nVtx] = hg->pref_part[i];
/*            uprintf(hgc, "match[%d (gno=%d)] = %zd   new vtxno=%d\n", i, VTX_LNO_TO_GNO(hg, i), match[i], c_hg->nVtx);*/
            ++c_hg->nVtx;
        }
    }
    
    for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertices */
      if (match[i] != VTX_LNO_TO_GNO(hg, i))  {
              ZOLTAN_GNO_TYPE gx = match[i], proc = VTX_TO_PROC_X(hg, gx);
                            
              if (proc != me)   {             /* owner is external */
                  LevelMap[i] = -gx - 1;     /* prepare to send data to owner */
                  size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
                  listgno[count]   = gx;                    /* listgno of vtx's to send */
                  listproc[count]  = proc;                  /* proc to send to */
                  listlno[count++] = i;                     /* lno of my match to gno */
/*                  uprintf(hgc, "EXTMAT:  match[%d (gno=%d)] = %zd \n", i, VTX_LNO_TO_GNO(hg, i), match[i]);*/
      } else   { /* owner is local */
                  LevelMap[i] = LevelMap[VTX_GNO_TO_LNO(hg, gx)];
/*                  uprintf(hgc, "LOCMAT:  match[%d (gno=%zd)] = %zd   new vtxno=%d\n", i, VTX_LNO_TO_GNO(hg, i), match[i], LevelMap[i]);*/
              }
      }
#ifdef KDDKDD_DEBUG
if (VTX_LNO_TO_GNO(hg, i) == 35 || VTX_LNO_TO_GNO(hg, i) == 65 || VTX_LNO_TO_GNO(hg, i) == 66) printf("%d MATCH %d (%f %f %f) to %d; LevelMap = %d\n", zz->Proc, VTX_LNO_TO_GNO(hg, i), hg->coor[i*3], hg->coor[i*3+1], hg->coor[i*3+2], match[i], LevelMap[i]);
#endif
    }
    *LevelSndCnt = count;
/*      errexit("this type of coarsening is not implemented yet"); */
  }
  
  /* Write list of: gno, lno, fixed_part, pref_part, size, vedge * size, weight * VtxWeightDim */

  /* sizes of parts of the messages as multiples of integers */

  gno_size = sizeof(ZOLTAN_GNO_TYPE) / sizeof(int);
  alt_field_count = (hgp->UseFixedVtx ? 1 : 0) + (hgp->UsePrefPart ? 1 : 0);
  int_size = alt_field_count + 2;
  float_size = (hg->VtxWeightDim * sizeof(float)) / sizeof(int);

  size += (count * (int_size + float_size + gno_size));

  if (size > 0 && !(buffer = (char *) ZOLTAN_MALLOC (size * sizeof(int))))
      MEMORY_ERROR;

  b = buffer;
  for (i = 0; i < count; i++)  {

    int lno=listlno[i], sz=hg->vindex[lno+1] - hg->vindex[lno];

    gnoptr = (ZOLTAN_GNO_TYPE *)b;
    intptr = (int *)(gnoptr + 1);
    floatptr = (float *)(intptr + 2 + alt_field_count + sz);

    b = (char *)(floatptr + hg->VtxWeightDim);

    *gnoptr = listgno[i];             /* destination vertex gno */

    *intptr++ = listlno[i];             /* source lno */
    if (hg->fixed_part)
      *intptr++ = hg->fixed_part[lno];
    if (hg->pref_part)
      *intptr++ = hg->pref_part[lno];

    *intptr++ = sz;

    for (j=0; j < sz; j++){
      *intptr++ = hg->vedge[hg->vindex[lno] + j];
    }
      
    if (hg->VtxWeightDim > 0){
      for (j=0; j < hg->VtxWeightDim; j++){
        *floatptr++ = hg->vwgt[lno*hg->VtxWeightDim + j];
      }
    }

    msg_size[i] = (b - (char *)gnoptr) / sizeof(int);
  }    
  if (hg->nDim) { /* Only process coords for coarse graph if fine graph has them */
    if (count &&
	!(coordbuf = (double *) ZOLTAN_MALLOC(count * hg->nDim * sizeof(double))))
      MEMORY_ERROR;

    /* Add local coordinates to their send buffer */
    for (i = 0; i < count; i++){
      int lno = listlno[i];
      for (j = 0; j < hg->nDim; j++)
	coordbuf[i * hg->nDim + j] = hg->coor[lno * hg->nDim + j];
#ifdef KDDKDD_DEBUG
if (VTX_LNO_TO_GNO(hg, lno) == 35 || VTX_LNO_TO_GNO(hg, lno) == 65 || VTX_LNO_TO_GNO(hg, lno) == 66) printf("%d LOADING %d (%f %f %f) (%f %f %f)\n", zz->Proc, VTX_LNO_TO_GNO(hg, lno), hg->coor[lno*hg->nDim], hg->coor[lno*hg->nDim+1], hg->coor[lno*hg->nDim+2], coordbuf[i*hg->nDim], coordbuf[i*hg->nDim+1], coordbuf[i*hg->nDim+2]);
#endif
    }
  }
  /* Create comm plan. */
  Zoltan_Comm_Create(comm_plan, count, listproc, hgc->row_comm, PLAN_TAG, 
                      &size); /* we'll use size for coords and then resize*/
  
  if (hg->nDim) {
    if (size &&
	!(coordrecbuf = (double *) ZOLTAN_MALLOC(size * hg->nDim * sizeof(double))))
      MEMORY_ERROR;

    /* No need for resize yet, since coordinates won't be of variable sizes */
    Zoltan_Comm_Do(*comm_plan, PLAN_TAG+1, (char *)coordbuf, sizeof(double) * hg->nDim,
		   (char *)coordrecbuf);
  
    /* Allocate coordinate array for coarse hgraph */
    if (c_hg->nVtx && (
         !(c_hg->coor = (double *) ZOLTAN_CALLOC(c_hg->nVtx * hg->nDim, sizeof(double)))
      || !(coorcount  = (double *) ZOLTAN_CALLOC(c_hg->nVtx, sizeof(double)))))
      MEMORY_ERROR;

    /* Accumulating on-processor coordinates */
    for (i = 0; i < hg->nVtx; i++) {
      ZOLTAN_GNO_TYPE ni = LevelMap[i];
      if (ni >= 0) {
        double hg_vwgt = hg->vwgt[i*hg->VtxWeightDim];
	for (j = 0; j < hg->nDim; j++)
	  c_hg->coor[ni*hg->nDim + j] += (hg_vwgt * hg->coor[i*hg->nDim + j]);
	coorcount[ni] += hg_vwgt;
#ifdef KDDKDD_DEBUG
if (VTX_LNO_TO_GNO(hg, i) == 35 || VTX_LNO_TO_GNO(hg, i) == 65 || VTX_LNO_TO_GNO(hg, i) == 66) printf("%d SUMMING %d (%f %f %f) into ni %d coorcount %f\n", zz->Proc, VTX_LNO_TO_GNO(hg, i), hg->coor[i*3], hg->coor[i*3+1], hg->coor[i*3+2], ni, coorcount[ni]);
#endif
      }
    }
  }
  
  /* call Comm_Resize since we have variable-size messages */
  Zoltan_Comm_Resize(*comm_plan, msg_size, PLAN_TAG+2, &size); 

  /* Allocate receive buffer. */
  /* size is the size of the received data, measured in #ints */
  if (size  && (
       !(rbuffer = (char *) ZOLTAN_MALLOC (size * sizeof(int)))
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
  Zoltan_Comm_Do(*comm_plan, PLAN_TAG+3, (char *)buffer, sizeof(int), (char *)rbuffer);
     
  /* Allocate vertex weight array for coarse hgraph */
  if (c_hg->nVtx > 0 && c_hg->VtxWeightDim > 0 &&
      !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx * c_hg->VtxWeightDim,
                                             sizeof(float))))
      MEMORY_ERROR;
  for (i=0; i < hg->nVtx; ++i) {
      ZOLTAN_GNO_TYPE ni=LevelMap[i];
      if (ni>=0)
          for (j=0; j<hg->VtxWeightDim; ++j)
              c_hg->vwgt[ni*hg->VtxWeightDim+j] += hg->vwgt[i*hg->VtxWeightDim+j];
  }
      
  /* index all received data for rapid lookup */ 
  *LevelCnt   = 0;
  b = rbuffer;
  b_end = rbuffer + (size * sizeof(int));
  if (hg->nDim)
    doubleptr = (double *)coordrecbuf;

  while (b < b_end){
    int j, sz, source_lno;
    /* ZOLTAN_GNO_TYPE lno; */
    int lno;               /* we get memory errors using ZOLTAN_GNO_TYPE to index int arrays TODO64 */

    gnoptr = (ZOLTAN_GNO_TYPE *)b;
    intptr = (int *)(gnoptr + 1);
    sz = intptr[1 + alt_field_count];
    floatptr = (float *)(intptr + 2 + alt_field_count + sz);
    b = (char *)(floatptr + hg->VtxWeightDim);
    
    source_lno              = *intptr++;
    lno = VTX_GNO_TO_LNO (hg, gnoptr[0]);
       
    if (hgp->UseFixedVtx)  {
      int fixed = *intptr++;
      c_hg->fixed_part [LevelMap[lno]] = (fixed >= 0) ? fixed : hg->fixed_part[lno];
    } 
    if (hgp->UsePrefPart)  {
      int pref = *intptr++;
      c_hg->pref_part [LevelMap[lno]] = (pref >= 0) ? pref : hg->pref_part[lno];
    } 
     
    (*LevelData)[(*LevelCnt)++] = source_lno;
    (*LevelData)[(*LevelCnt)++] = lno;              /* to lookup in part[] */

    lno = (int)LevelMap[lno];
    if (hg->nDim) {
#ifdef KDDKDD_DEBUG
if (gnoptr[0] == 35 || gnoptr[0] == 65 || gnoptr[0] == 66) printf("%d RECEIVED %d (%f %f %f) into lno %d coorcount %f doublptr %x\n", zz->Proc, gnoptr[0], *doubleptr, *(doubleptr+1), *(doubleptr+2), lno, coorcount[lno]+1., doubleptr);
#endif
      for (j = 0; j < hg->nDim; j++){
        /* NOTE:  This code must preceed accumulation of vwgt below so that 
         * floatptr is correct. */
	c_hg->coor[lno * hg->nDim + j] += (*floatptr * *doubleptr++);
      }
      coorcount[lno] += *floatptr;
    }
    
    for (j=0; j<hg->VtxWeightDim; ++j)
        c_hg->vwgt[lno*hg->VtxWeightDim+j] += *floatptr++;
    
    intptr++;                /* skip sz */
    for (j=0; j<sz; ++j)      /* count extra vertices per edge */
        ++ahindex[*intptr++];
  }

  if (c_hg->nDim) {
    /* Average coordinates */
    for (i = 0; i < c_hg->nVtx; i++)
      for (j = 0; j < c_hg->nDim; j++)
	c_hg->coor[i * hg->nDim + j] = c_hg->coor[i * hg->nDim + j] / coorcount[i];
    ZOLTAN_FREE(&coorcount);
  }
  
  for (i=0; i<hg->nEdge; ++i) /* prefix sum over ahindex */
    ahindex[i+1] += ahindex[i];
  /* now prepare ahvertex */

  b = rbuffer;
  b_end = rbuffer + (size * sizeof(int));

  while (b < b_end){
    int lno, sz, j;

    gnoptr = (ZOLTAN_GNO_TYPE *)b;
    intptr = (int *)(gnoptr + 1);
    sz = intptr[1 + alt_field_count];
    floatptr = (float *)(intptr + 2 + alt_field_count + sz);
    b = (char *)(floatptr + hg->VtxWeightDim);

    intptr += (2 + alt_field_count);    

    lno=VTX_GNO_TO_LNO (hg, gnoptr[0]);

    for (j=0; j<sz; ++j)      /* count extra vertices per edge */
        ahvertex[--ahindex[*intptr++]] = lno;
  }

  Zoltan_Multifree (__FILE__, __LINE__, 4, &buffer, &rbuffer, &coordbuf, &coordrecbuf);
  
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
          int nvno=(int)LevelMap[hg->hvertex[j]];
          if (nvno>=0 && vmark[nvno]!=i) {
              c_hg->hvertex[idx++] = (int)nvno;
              vmark[nvno] = i;
          }
      }
      /* now go over the received vertices */
      for (j=ahindex[i]; j<ahindex[i+1]; ++j) {
          int nvno= (int)LevelMap[ahvertex[j]];
          if (nvno>=0 && vmark[nvno]!=i) {
              c_hg->hvertex[idx++] = (int)nvno; 
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
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->comerge, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer->coshuffle, hgc->Communicator);
  }
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_coarse += t_cur;
  t_redhash = -t_cur;
#endif

  for (i=0; i < c_hg->nEdge; ++i) { /* compute size and hashvalue */
    hlsize[i] = c_hg->hindex[i+1]-c_hg->hindex[i];
    lhash[i] = hashValue(hg, hlsize[i], (int *)&c_hg->hvertex[c_hg->hindex[i]]);
  }

  /* UVC TODO to compute global hash; right now we'll use SUM (UVC TODO: try:bitwise xor);
     we need to check if this is good, if not we need to find a better way */
  if (c_hg->nEdge) 
      MPI_Allreduce(lhash, hash, c_hg->nEdge, MPI_UNSIGNED, MPI_SUM, hgc->row_comm);

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
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->coshuffle, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer->coremove, hgc->Communicator);
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
  Zoltan_quicksort_pointer_inc_int_int((int *)ids, (int *)lhash, hsize, 0, size-1); 
/*  uqsort_ptr_uint_int(size, lhash, hsize, ids); */

#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_sort += t_cur;
  t_iden = -t_cur;
#endif

  emptynetsize = (size >> 5)+1;
  emptynets = (int*) ZOLTAN_CALLOC(emptynetsize, sizeof(int));
  idennets = (int *) ZOLTAN_MALLOC(size*sizeof(int));

  if ((emptynetsize && !emptynets) || (size && !idennets)) MEMORY_ERROR;

  iden = listproc; /* just better variable name */
#ifdef _DEBUG1
  emptynetcnt=0;
#endif
  for (j=0; j<size; ++j) {
      if (hlsize[j])
          iden[j] = 0;
      else {
          iden[j] = -1;
          BITSET(emptynets, j);
#ifdef _DEBUG1          
          ++emptynetcnt;
#endif
      }
  }

  count = idx = 0;
  for (j=0; j<size; ++j) {
      int n1=ids[j];
      if (!iden[n1]) {
          int last=1+n1, minv=last;

          for (i = j+1; i<size && lhash[n1] == lhash[ids[i]] && hsize[n1]==hsize[ids[i]]; ++i) {
              int n2=ids[i];
              
              if (!iden[n2] && hlsize[n1]==hlsize[n2]
                  && !memcmp(&ahvertex[ahindex[n1]], &ahvertex[ahindex[n2]], sizeof(int)*hlsize[n1])) {
                  iden[n2] = last; /* n2 is potentially identical to n1 */
                  last = 1+n2;
                  minv = (last<minv) ? last : minv;
                  ++count;
              }                  
          }
          /* iden[last] is a link list (in array) now make the
             all identical nets to point the same net with the smallest id;
             it will be needed in identicalOperator; we'll zero(clear) the
             original net (net with the smallest id) after this loop */
          if (last!=1+n1) {  /* some nets are identical to n1 */
              idennets[idx++] = -minv; /* we need to first write minv */
              if (minv!=1+n1)     /* then we write n1 if it is not minv itself */
                  idennets[idx++] = 1+n1;
          }
          while (last!=1+n1) {
              int prev=iden[last-1];

              if (last!=minv)
                  idennets[idx++] = last;
              iden[last-1] = minv;
              last = prev;
          }
          iden[n1] = minv;
      }
  }
  
  for (i=0; i<size; ++i){
      if (iden[i]==1+i) /* original net; clear iden */
          iden[i] = 0;
  }

  
#ifdef _DEBUG1
  MPI_Barrier(hgc->Communicator);
  t_cur = MPI_Wtime();
  t_iden += t_cur;
  t_userredop = -t_cur;

#ifndef _DEBUG2
  if (!hgc->myProc)
#endif
  uprintf(hgc, "#Loc.Iden= %7d   (Computed= %d)   #Comp.PerNet= %.2lf     ElapT= %.3lf\n", count+emptynetcnt, count, (double) idx / (double) size, MPI_Wtime()-t_all);
#endif


  /* a simple message size opt; proc with largest message is root; so it won't sent it */
  Zoltan_PHG_Find_Root(idx+emptynetsize, hgc->myProc_x, hgc->row_comm, 
                       &i, &rootRank); 
  
  /* communicate empty-nets to row-root */
  if (hgc->myProc_x == rootRank){ 
      allemptynets = (int *) ZOLTAN_MALLOC(emptynetsize*hgc->nProc_x*sizeof(int));
      if ((emptynetsize*hgc->nProc_x > 0) && !allemptynets) MEMORY_ERROR;
  }
  MPI_Gather(emptynets, emptynetsize, MPI_INT, allemptynets, emptynetsize, MPI_INT, rootRank, hgc->row_comm);

  /* communicate identical nets to row-root */
  MPI_Gather(&idx, 1, MPI_INT, msg_size, 1, MPI_INT, rootRank, hgc->row_comm);
  if (hgc->myProc_x == rootRank) {
      int recsize = 0;
      for (i = 0; i < hgc->nProc_x; i++) 
          recsize += msg_size[i];

      allidennets = (int *) ZOLTAN_MALLOC(recsize*sizeof(int));
      idennetdest = (int *) ZOLTAN_MALLOC(hgc->nProc_x*sizeof(int));

      if ((recsize && !allidennets) || (hgc->nProc_x && !idennetdest)) MEMORY_ERROR;
      
      idennetdest[0] = 0;
      for (i = 1; i < hgc->nProc_x; i++)
          idennetdest[i] = idennetdest[i-1] + msg_size[i-1];
  }  
  MPI_Gatherv (idennets, idx, MPI_INT, allidennets, msg_size, idennetdest, MPI_INT, rootRank, hgc->row_comm);

  
  ip = (int *) lhash;
  if (hgc->myProc_x == rootRank) {
      for (j=0; j<hgc->nProc_x; ++j) {
          if (j!=rootRank) {
              int *enets=allemptynets+j*emptynetsize;
              int *inets=allidennets+idennetdest[j], rootnet=-1;
              int *x, *y, *a, *b;

              for (i=0; i<size; ++i) {
                  ip[i] = BITCHECK(enets, i) ? -1 : 0;
              }

              for (i=0; i<msg_size[j]; ++i) {
                  if (inets[i] < 0) 
                      rootnet = -inets[i];
                  else 
                      ip[inets[i]-1] = rootnet;
              }

              /* merge received identical net info with the local one */
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
              
              x = ip;   y = iden;
              memcpy(ids, iden, size*sizeof(int));
              a = ip-1; b=ids-1; /* net ids in the a and b are base-1 numbers */
              for (i=0; i < size; ++i, ++x, ++y) {
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
      }
      ip = iden;
  }
  ZOLTAN_FREE(&ids); 
  MPI_Bcast(ip, size, MPI_INT, rootRank, hgc->row_comm);

  Zoltan_Multifree(__FILE__, __LINE__, 5, &idennets, &emptynets, &allemptynets, &allidennets, &idennetdest);
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
      if (ip[i]==-1 && (hlsize[i]>0 || hsize[i]>0))
          errexit("ip[%d]==-1 but hlsize[%d] = %d hsize[%d]=%d", i, i, hlsize[i], i, hsize[i]);
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
  uprintf(hgc, "#GlobIden= %7d    SuccessRate= %.1lf%%    ElapT= %.3lf HashOpT=%.3lf (%.1lf%%)\n", idx, 100.0 * idx / (double) (count+emptynetcnt), t_cur, t_userredop, 100.0*t_userredop/t_cur);
#endif
  MPI_Allreduce(&idx, &totiden, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  MPI_Allreduce(&me, &totsize1, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  if (!hgc->myProc)
      uprintf(hgc, "Level %d - Orig #Nets=%d    Iden=%d    Size-0/1=%d CurT=%.3lf user-def redopT=%.3lf\n", hg->info, hg->dist_y[hgc->nProc_y], totiden, totsize1, t_cur, t_userredop);  
#endif

  Zoltan_Multifree(__FILE__, __LINE__, 3, &c_hg->hindex, &c_hg->hvertex, &c_hg->ewgt);

  c_hg->hindex = (int*)ZOLTAN_MALLOC((c_hg->nEdge+1)*sizeof(int));
  if (!c_hg->hindex) MEMORY_ERROR;
  if (c_hg->nPins &&
      !(c_hg->hvertex=(int*)ZOLTAN_MALLOC(c_hg->nPins*sizeof(int))))
      MEMORY_ERROR;
  if (c_hg->nEdge &&
      !(c_hg->ewgt=(float*)ZOLTAN_MALLOC(c_hg->nEdge*c_hg->EdgeWeightDim*sizeof(float))))
      MEMORY_ERROR;
  if (c_hg->nEdge &&
      !(c_hg->esize=(ZOLTAN_GNO_TYPE*)ZOLTAN_MALLOC(c_hg->nEdge*sizeof(ZOLTAN_GNO_TYPE))))
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
          c_hg->esize[ni] = hsize[i];
          memcpy(&c_hg->ewgt[ni], &c_ewgt[i], c_hg->EdgeWeightDim*sizeof(float));
          memcpy(&c_hg->hvertex[idx], &ahvertex[ahindex[i]], hlsize[i]*sizeof(int));
          ++ni;
          idx += hlsize[i];
      }
  c_hg->hindex[c_hg->nEdge] = idx;

  if (time_details) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->coremove, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer->cotheend, hgc->Communicator);
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
  if (!(c_hg->dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((hgc->nProc_x+1), sizeof(ZOLTAN_GNO_TYPE)))
	 || !(c_hg->dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((hgc->nProc_y+1), sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;

  tmp_gno = (ZOLTAN_GNO_TYPE)c_hg->nVtx;
  MPI_Scan(&tmp_gno, c_hg->dist_x, 1, zoltan_gno_mpi_type, MPI_SUM, hgc->row_comm);
  MPI_Allgather(c_hg->dist_x, 1, zoltan_gno_mpi_type, &(c_hg->dist_x[1]), 1, zoltan_gno_mpi_type, hgc->row_comm);
  c_hg->dist_x[0] = 0;
  
  tmp_gno = (ZOLTAN_GNO_TYPE)c_hg->nEdge;
  MPI_Scan(&tmp_gno, c_hg->dist_y, 1, zoltan_gno_mpi_type, MPI_SUM, hgc->col_comm);
  MPI_Allgather(c_hg->dist_y, 1, zoltan_gno_mpi_type, &(c_hg->dist_y[1]), 1, zoltan_gno_mpi_type, hgc->col_comm);
  c_hg->dist_y[0] = 0;  

  ierr = Zoltan_HG_Create_Mirror(zz, c_hg);
#ifdef _DEBUG1
  if (c_hg->fixed_part)
      for (i = 0; i < c_hg->nVtx; i++){
          if (c_hg->fixed_part[i] == -2)
              printf ("RTHRTH BAD COARSENING for FIXED VERTICES\n"); 
      }
  if (c_hg->pref_part)
      for (i = 0; i < c_hg->nVtx; i++){
          if (c_hg->pref_part[i] == -2)
              uprintf(hgc, "*******BAD COARSENING for PREF[%d] is unassigned\n", i);
      }
  
  MPI_Barrier(hgc->Communicator);  
  t_mirror += MPI_Wtime();
#endif
 End:
#ifdef _DEBUG1
  t_mirror -= MPI_Wtime();
#endif
  Zoltan_Multifree (__FILE__, __LINE__, 19,
                    &listgno, &listlno, &listproc, &msg_size,
                    &buffer, &rbuffer, &ahindex, &ahvertex, &vmark,
                    &hlsize, &hsize, &lhash, &hash, &c_ewgt,
                    &idennets, &emptynets, &allemptynets, &allidennets, &idennetdest
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
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->cotheend, hgc->Communicator);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
