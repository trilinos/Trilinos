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
#include "zoltan_comm.h"

    
/* #define USE_NEW_COARSENING  */
#define REMOVE_IDENTICAL_NETS
    
    
#define PIN_OVER_ALLOC 1.2  /* Overallocate pins by 20% */
#define PLAN_TAG 32010      /* tag for comm plan */



#ifdef USE_NEW_COARSENING



typedef struct _TagSortItem
{
int id;
unsigned long val;
} SortItem;
    
/* UVC:
   Following quicksort routines are modified from
   Numerical Recipes Software
   these versions of quicksort seems to perform
   better than the ones that exist in Zoltan
*/
    
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define M             7
#define NSTACK        50

static void uqsort(int n, SortItem *arr)
{
    int           i, ir=n, j, k, l=1;
    int           jstack=0, istack[NSTACK];
    unsigned long aval;
    SortItem      a, temp;
    
    --arr;
    for (;;) {
        if (ir-l < M) {
            for (j=l+1;j<=ir;j++) {
                a=arr[j];
                aval = a.val;
                for (i=j-1;i>=1;i--) {
                    if (arr[i].val <= aval) 
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
            if (arr[l+1].val > arr[ir].val) 
                SWAP(arr[l+1],arr[ir]);
            if (arr[l].val > arr[ir].val) 
                SWAP(arr[l],arr[ir]);
            if (arr[l+1].val > arr[l].val) 
                SWAP(arr[l+1],arr[l]);

            i=l+1;
            j=ir;
            a=arr[l];
            aval = a.val;
            for (;;) {
                do i++; while (arr[i].val < aval);
                do j--; while (arr[j].val > aval);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
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

static long hashValue(int n, int *ar)
{
    unsigned long l=0;
    int i;

    /* Chris Torek's hash function */
    for (i=0; i<n; ++i) {
        l *= 33;
        l += ar[i];
    }
    return l;
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
  struct Zoltan_Comm_Obj  **comm_plan
    )   
{
  char *yo = "Zoltan_PHG_Coarsening";
  int  ierr=ZOLTAN_OK, i, j, count, size, me, idx, ni;
  int  *vmark=NULL, *listgno=NULL, *listlno=NULL, *listproc=NULL, *msg_size=NULL, *ip;
  int  *ahindex=NULL, *ahvertex=NULL, *hsize=NULL, *hlsize=NULL;
  float *c_ewgt=NULL;
  unsigned long *hash=NULL, *lhash=NULL;
  char *buffer=NULL, *rbuffer=NULL;
  PHGComm          *hgc = hg->comm;
  struct Zoltan_Comm_Obj  *plan=NULL;
  SortItem *slist=NULL;
  
  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  c_hg->comm    = hg->comm;         /* set communicators */
  c_hg->info    = hg->info + 1;     /* for debugging */
  c_hg->coor    = NULL;             /* currently we don't use coordinates */
  c_hg->nDim    = hg->nDim;    
  c_hg->vmap    = NULL;             /* only needed by rec-bisec */
  c_hg->ratio   = hg->ratio;        /* for "global" recursive bisectioning */
  c_hg->redl    = hg->redl;         /* to stop coarsening near desired count */
  c_hg->VtxWeightDim  = hg->VtxWeightDim;
  
  /* (over) estimate number of external matches that we need to send data to */
  count = 0; 
  for (i = 0; i < hg->nVtx; i++)
    if (match[i] < 0)
      ++count;
 
  if (hg->nVtx > 0 && !(vmark = (int*) ZOLTAN_CALLOC(hg->nVtx,  sizeof(int))))
      MEMORY_ERROR;

  size = MAX(count, hg->nEdge);
  if (size && (
      !(listproc  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))
   || !(msg_size  = (int*) ZOLTAN_MALLOC (size * sizeof(int)))))
      MEMORY_ERROR;

  if (count > 0 && (
      !(listgno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(listlno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(*LevelData= (int*) ZOLTAN_MALLOC (count * sizeof(int) * 2))))
      MEMORY_ERROR;

  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of resulting coarse vertices. */
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* short name, convenience variable */
  size  = 0;                      /* size (in ints) to communicate */
  count = 0;                      /* number of vertices to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int gx = -match[i]-1, proc = VTX_TO_PROC_X(hg,gx);
      
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
      while (!vmark[v])  {
        LevelMap[v] = c_hg->nVtx;         /* next available coarse vertex */      
        vmark[v] = 1;  /* flag this as done already */
        v = match[v];  
      }
      ++c_hg->nVtx;
    }
  }
  *LevelSndCnt = count;
  
  /* size and allocate the send buffer */
  size += ((3 + hg->VtxWeightDim) * count);
  if (size > 0 && !(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))
      MEMORY_ERROR;
  
  
  /* Message is list of <gno, vweights, gno's edge count, list of edge lno's> */
  /* We pack ints and floats into same buffer */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
      int lno=listlno[i], sz=hg->vindex[lno+1] - hg->vindex[lno];
      int *ip_old = ip;
      
    *ip++ = listlno[i];             /* source lno */
    *ip++ = listgno[i];             /* destination vertex gno */
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
    || !(lhash   = (unsigned long *) ZOLTAN_MALLOC(hg->nEdge*sizeof(unsigned long)))
    || !(hash    = (unsigned long *) ZOLTAN_MALLOC(hg->nEdge*sizeof(unsigned long)))
          ))
      MEMORY_ERROR;

  /* Comm_Do sends personalized messages of variable sizes */
  Zoltan_Comm_Do(*comm_plan, PLAN_TAG+2, buffer, sizeof(int), rbuffer);

  /* Allocate vertex weight array for coarse hgraph */
  if (hg->nVtx > 0 && hg->VtxWeightDim > 0 &&
      !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx * hg->VtxWeightDim,
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
    int j, sz, source_lno=ip[i++];
    int lno=VTX_GNO_TO_LNO (hg, ip[i++]);
    float *pw=(float*) &ip[i];

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

    i += 2+hg->VtxWeightDim;  /* skip slno+gno+vertex weights */
    sz = ip[i++];             /* get sz to a diff var, skip size*/
    for (j=0; j<sz; ++j)      /* count extra vertices per edge */
        ahvertex[--ahindex[ip[i+j]]] = lno;
    i += sz;
  }
  Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);
  
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
#ifdef REMOVE_IDENTICAL_NETS
      uqsorti(idx-sidx, &c_hg->hvertex[sidx]);
#endif
  }
  c_hg->hindex[hg->nEdge] = c_hg->nPins = idx;
  Zoltan_Multifree(__FILE__, __LINE__, 3, &vmark, &ahvertex, &ahindex);


#ifdef REMOVE_IDENTICAL_NETS
  for (i=0; i < c_hg->nEdge; ++i) { /* compute size and hashvalue */
      hlsize[i] = c_hg->hindex[i+1]-c_hg->hindex[i];
      lhash[i] = hashValue(hlsize[i], &c_hg->hvertex[c_hg->hindex[i]]);
  }

  /* UVC TODO to compute global hash; right now we'll use bitwise xor;
     we need to check if this is good, if not we need to find a better way */
  MPI_Allreduce(lhash, hash, c_hg->nEdge, MPI_LONG, MPI_PROD, hgc->row_comm);


  for (i=0; i < hg->nEdge; ++i)  /* decide where to send */
      listproc[i] = (int) (hash[i] % hgc->nProc_y);
  
  Zoltan_Comm_Create(&plan, c_hg->nEdge, listproc, hgc->col_comm, PLAN_TAG+10, 
                     &size);

  /* send hash values */
  if (size > hg->nEdge) {
      Zoltan_Multifree(__FILE__, __LINE__, 2, &lhash, &listproc);
      if (!(lhash=(unsigned long *) ZOLTAN_MALLOC(size * sizeof(unsigned long)))
          || !(listproc=(int *) ZOLTAN_MALLOC(size * sizeof(int))))
          MEMORY_ERROR;
  }
  Zoltan_Comm_Do(plan, PLAN_TAG+11, (char *) hash, sizeof(unsigned long), (char *) lhash);
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
  /* now ewgt sizes */
  Zoltan_Comm_Do(plan, PLAN_TAG+12, (char *) c_hg->ewgt, sizeof(int), (char *) c_ewgt);

  /* now vertices of hyperedges */
  Zoltan_Comm_Resize(plan, hlsize, PLAN_TAG+13, &idx);
  if (idx && !(ahvertex = (int *) ZOLTAN_MALLOC(idx * sizeof(int))))
      MEMORY_ERROR; 
  Zoltan_Comm_Do(plan, PLAN_TAG+14, (char *) c_hg->hvertex, sizeof(int), (char *) ahvertex);
  ZOLTAN_FREE(&hlsize);    hlsize=ip;

  Zoltan_Comm_Destroy (&plan);

  /* in order to find identical nets; we're going to sort hash values and compare them */
  if (size && !(slist = (SortItem *) ZOLTAN_MALLOC(size * sizeof(SortItem))))
      MEMORY_ERROR;
  ahindex[0] = 0;
  for (i=0; i<size; ++i) {
      ahindex[1+i] = ahindex[i] + hlsize[i];
      slist[i].id = i;
      slist[i].val = lhash[i];
  }
  uqsort(size, slist);
      
  idx = 0;
  memset(listproc, 0, sizeof(int)*size);
  for (j=0, i=1; i<size; ++i) {
      while (i<size && slist[j].val==slist[i].val) {
          int n1=slist[j].id, n2=slist[i].id;

          if (hlsize[n1]==hlsize[n2] &&
              !memcmp(&ahvertex[ahindex[n1]], &ahvertex[ahindex[n2]], sizeof(int)*hlsize[n1])) {
              listproc[n2] = 1+n1; /* n2 is potentially identical to n1 */
              ++idx;
          }
          ++i;
      }
      j = i;
  }  
  ZOLTAN_FREE(&slist);
  ip = (int *) lhash; /* UVC: should be safe; it is always sizeof(int)<=sizeof(long) */
  MPI_Allreduce(listproc, ip, size, MPI_INT, MPI_LAND, hgc->row_comm);
  MPI_Allreduce(hlsize, hsize, size, MPI_INT, MPI_SUM, hgc->row_comm);
  
  c_hg->nPins = 0;
  c_hg->nEdge = 0;
  for (i=0; i<size; ++i) {
      if (ip[i]) { /* identical net */
          int n1=listproc[i]-1; /* i is identical to n1 */
          for (j=0; j<c_hg->EdgeWeightDim; ++j)
              c_ewgt[n1*c_hg->EdgeWeightDim + j] += c_ewgt[i*c_hg->EdgeWeightDim + j];
          ip[i] = 0; /* ignore */
      } else if (hsize[i]>1) { /* size not 0/1 and NOT identical */          
          c_hg->nPins += hlsize[i];
          ++c_hg->nEdge;
          ip[i] = 1; /* Don't ignore */
      } else 
          ip[i] = 0; /* ignore size 0/1 nets*/
  }
  
  Zoltan_Multifree(__FILE__, __LINE__, 3, &c_hg->hindex, &c_hg->hvertex, &c_hg->ewgt);

  c_hg->hindex = (int*)ZOLTAN_MALLOC((c_hg->nEdge+1)*sizeof(int));
  if (c_hg->nPins &&
      !(c_hg->hvertex=(int*)ZOLTAN_MALLOC(c_hg->nPins*sizeof(int))))
      MEMORY_ERROR;
  if (c_hg->nEdge &&
      !(c_hg->ewgt=(float*)ZOLTAN_MALLOC(c_hg->nEdge*c_hg->EdgeWeightDim*sizeof(float))))
      MEMORY_ERROR;

  
  for (idx=ni=i=0; i<size; ++i)
      if (ip[i]) {
          c_hg->hindex[ni] = idx;
          memcpy(&c_hg->ewgt[ni], &c_ewgt[i], c_hg->EdgeWeightDim*sizeof(float));
          memcpy(&c_hg->hvertex[idx], &ahvertex[ahindex[i]], hlsize[i]*sizeof(int));
          ++ni;
          idx += hlsize[i];
      }
  c_hg->hindex[c_hg->nEdge] = idx;

#ifdef _DEBUG
  if (idx!=c_hg->nPins || ni!=c_hg->nEdge)
      errexit("idx(%d)!=c_hg->nPins(%d) || ni(%d)!=c_hg->nEdge(%d)", idx, c_hg->nPins, ni, c_hg->nEdge);
#endif
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

  ierr = Zoltan_HG_Create_Mirror(zz, c_hg);
 End:
  Zoltan_Multifree (__FILE__, __LINE__, 14,
                    &listgno, &listlno, &listproc, &msg_size,
                    &buffer, &rbuffer, &ahindex, &ahvertex, &vmark,
                    &hlsize, &hsize, &lhash, &hash, &c_ewgt 
                    );
  return ierr;
}

#else
int Zoltan_PHG_Coarsening
( ZZ     *zz,         /* the Zoltan data structure */
  HGraph *hg,         /* information about hypergraph, weights, etc. */
  int    *match,      /* Matching, Packing or Grouping array */
  HGraph *c_hg,       /* output: the coarse hypergraph */
  int    *LevelMap,
  int    *LevelCnt,
  int    *LevelSndCnt,
  int   **LevelData,  /* information to reverse coarsenings later */
  struct Zoltan_Comm_Obj  **comm_plan
  )   
{
  int i, j, vertex, edge, me, size, count, pincnt=hg->nPins;
  int *ip, *ip_old;
  int *cmatch=NULL, *used_edges=NULL, *c_vindex=NULL, *c_vedge=NULL;
  int *listgno=NULL, *listlno=NULL,  *listproc=NULL;
  int *msg_size=NULL;
  float *pwgt;
  char *buffer=NULL, *rbuffer=NULL;
  PHGComm *hgc = hg->comm;
  char *yo = "Zoltan_PHG_Coarsening";

 
  ZOLTAN_TRACE_ENTER (zz, yo);  
  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  
  /* (over) estimate number of external matches that we need to send data to */
  for (count = i = 0; i < hg->nVtx; i++)
    if (match[i] < 0)
      ++count;
 
  if (hg->nVtx > 0 && (
      !(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int))))){
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }        
  if (count > 0 && (
      !(listgno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(listlno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(listproc  = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(msg_size  = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(*LevelData= (int*) ZOLTAN_MALLOC (count * sizeof(int) * 2))))    {   
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }        
  memcpy (cmatch, match, hg->nVtx * sizeof(int));   /* working copy of match[] */
       
  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of resulting coarse vertices. */
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* short name, convenience variable */
  size  = 0;                      /* size (in ints) to communicate */
  count = 0;                      /* number of vertices to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int proc, gx;
      gx = -match[i]-1;
      proc = VTX_TO_PROC_X(hg,gx);
      
      /* rule to determine "ownership" of coarsened vertices between procs */
      proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(proc, me) : MAX(proc, me);
      
      /* prepare to send data to owner */
      if (proc != me)   {             /* another processor owns this vertex */
        size += hg->vindex[i+1] - hg->vindex[i];  /* for send buffer sizing */ 
        listgno [count]   = gx;                   /* listgno of vtx's to send */
        listproc[count]   = proc;                 /* proc to send to */
        listlno [count++] = i;                    /* lno of my match to gno */
      }
      else      
        c_hg->nVtx++;         /* myProc owns the matching across processors */
    }
    else if (cmatch[i] >= 0)  {     /* local matching, packing and groupings */    
      c_hg->nVtx++;
      vertex = i;
      while (cmatch[vertex] >= 0)  {
        cmatch[vertex] = -cmatch[vertex] - 1;  /* flag this as done already */
        vertex         = -cmatch[vertex] - 1;  
      }
    }
  }
  *LevelSndCnt = count;
  
  /* size and allocate the send buffer */
  size += ((3 + hg->VtxWeightDim) * count);
  if (size > 0 && !(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))  {  
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }
        
  /* Message is list of <gno, vweights, gno's edge count, list of edge lno's> */
  /* We pack ints and floats into same buffer */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
    ip_old = ip;
    *ip++ = listlno[i];                            /* source lno */
    *ip++ = listgno[i];                            /* destination vertex gno */
    for (j = 0; j < hg->VtxWeightDim; j++) {
       /* EBEB Assume a float is no larger than an int. 
          This is usually true, but to be safe this trick should be avoided. */
       pwgt = (float*) ip++;                               /* vertex weight */
       *pwgt = hg->vwgt[listlno[i]*hg->VtxWeightDim+j] ;
    }
       
    *ip++ = hg->vindex[listlno[i]+1] - hg->vindex[listlno[i]];      /* count */
    for (j = hg->vindex[listlno[i]]; j < hg->vindex[listlno[i]+1]; j++)
      *ip++ = hg->vedge[j];    

    /* save mesg size in #ints */
    msg_size[i] = ip - ip_old;
  }    

  /* Create comm plan. */
  Zoltan_Comm_Create( comm_plan, count, listproc, hgc->row_comm, PLAN_TAG, 
   &i);

  /* call Comm_Resize since we have variable-size messages */
  Zoltan_Comm_Resize( *comm_plan, msg_size, PLAN_TAG+1, &size); 

  /* Allocate receive buffer. */
  /* size is the size of the received data, measured in #ints */
  if (size > 0 && !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))   {
    ZOLTAN_TRACE_EXIT (zz, yo);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Comm_Do sends personalized messages of variable sizes */
  Zoltan_Comm_Do(*comm_plan, PLAN_TAG+2, buffer, sizeof(int), rbuffer);

  /* index all received data for rapid lookup */
  ip = (int*) rbuffer;
  *LevelCnt = 0;
  for (i = 0; i < size; i++)  {
    int source_lno = ip[i++];
    int lno = VTX_GNO_TO_LNO (hg, ip[i]);
    (*LevelData)[(*LevelCnt)++] = source_lno;
    (*LevelData)[(*LevelCnt)++] = lno;              /* to lookup in part[] */
    
    cmatch [lno] = i;       /* place index into buffer into "match" */
    i++;                    /* skip destination gno */
    i += hg->VtxWeightDim;  /* skip vertex weights */
    i += ip[i];             /* skip hyperedges */
  }
   
  if (hg->nEdge>0 && !(used_edges = (int*)ZOLTAN_CALLOC(hg->nEdge,sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  if (hg->nVtx>0 && !(c_vindex = (int*)ZOLTAN_MALLOC((hg->nVtx+1)*sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  if (hg->nPins>0 && !(c_vedge = (int*)ZOLTAN_MALLOC (hg->nPins*sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }

  /* Allocate vertex weight array for coarse hgraph */
  if (hg->nVtx > 0 && hg->VtxWeightDim > 0) 
     c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx * hg->VtxWeightDim,
      sizeof(float));   
      
  /* Allocate edge weight array for coarse hgraph */
  if (hg->EdgeWeightDim > 0) {
    c_hg->ewgt =(float*)ZOLTAN_MALLOC(hg->nEdge*hg->EdgeWeightDim*sizeof(float));
    if (c_hg->ewgt == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    }    
    for (i = 0; i < hg->nEdge * hg->EdgeWeightDim; i++)
      c_hg->ewgt[i] = hg->ewgt[i];
  }

  /* Construct the LevelMap; match[vertex] is changed back to original value */
  /* Coarsen vertices (create vindex, vedge), sum up coarsened vertex weights */
  /* EBEB This code constructs the vertex-based arrays. It might be
     better to construct the edge-based arrays so we can easily remove edges. */
  c_hg->nPins = 0;                      /* count of coarsened pins */
  c_hg->nVtx  = 0;                      /* count of coarsened vertices */
  for (i = 0; i < hg->nVtx; i++)  {
    if (match[i] < 0 && cmatch[i] < 0)  /* match to external vertex, don't own */
      LevelMap [i] = match[i];         /* negative value => external vtx */         
    else if (match[i] < 0) {            /* match from external vertex, I own */
       LevelMap[i] = c_hg->nVtx;        /* next available coarse vertex */
       c_vindex[c_hg->nVtx] = c_hg->nPins;
       ip =  ((int*) rbuffer) + cmatch[i];      /* point to received data */               
       ip++;                                    /* skip over gno */
      
       for (j = 0; j < hg->VtxWeightDim; j++)  {
          pwgt = (float*) ip++;
          c_hg->vwgt[c_hg->nVtx*hg->VtxWeightDim+j] = *pwgt;
       }  
       count = *ip++;           /* extract edge count, advance to first edge */
       while (count--)  {
          edge = *ip++;           
          used_edges [edge] = i+1;
          if (c_hg->nPins >= pincnt)  {
             pincnt = 1 + PIN_OVER_ALLOC * pincnt;
             c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
          }
          c_vedge[c_hg->nPins++] = edge;
       }
  
      for (j = 0; j < hg->VtxWeightDim; j++)
        c_hg->vwgt[c_hg->nVtx * hg->VtxWeightDim + j]
         += hg->vwgt[i        * hg->VtxWeightDim + j] ;
            
      for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)  {
        if (used_edges [hg->vedge[j]] <= i)   {
          used_edges [hg->vedge[j]] = i+1;  
          if (c_hg->nPins >= pincnt)  {
             pincnt = 1 + PIN_OVER_ALLOC * pincnt;     
             c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
          }                  
          c_vedge[c_hg->nPins++] = hg->vedge[j];
        }      
      }        
      c_hg->nVtx++;          
    }
    else if (match[i] >= 0 && cmatch[i] < 0) { /* match/pack/group local vtx's */
      c_vindex[c_hg->nVtx] = c_hg->nPins;
      vertex = i;
      while (cmatch[vertex] < 0)  {    
        LevelMap[vertex] = c_hg->nVtx; 
        
        for (j = 0; j < hg->VtxWeightDim; j++)
          c_hg->vwgt[c_hg->nVtx * hg->VtxWeightDim + j]
           += hg->vwgt[vertex   * hg->VtxWeightDim + j] ;        
           
        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++)  {
          if (used_edges [hg->vedge[j]] <= i)  {
            used_edges [hg->vedge[j]] = i+1; 
            if (c_hg->nPins >= pincnt)  {
               pincnt = 1 + PIN_OVER_ALLOC * pincnt;            
               c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
            }                             
            c_vedge[c_hg->nPins++] = hg->vedge[j];
          }              
        }                       
        cmatch[vertex] = -cmatch[vertex] - 1;
        vertex         =  cmatch[vertex];
      }
      c_hg->nVtx++;
    }
  }
  c_vindex[c_hg->nVtx] = c_hg->nPins;
  ZOLTAN_FREE ((void**) &used_edges);
  
  /* Update the dist_x array after coarsening: allocate and fill */
  if (!(c_hg->dist_x = (int*)ZOLTAN_MALLOC ((hgc->nProc_x+1) * sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  } 
    
  MPI_Allgather (&(c_hg->nVtx), 1, MPI_INT, &(c_hg->dist_x[1]), 1, MPI_INT, 
                 hgc->row_comm);
  
  /* dist_x is the cumulative sum */
  c_hg->dist_x[0] = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    c_hg->dist_x[i+1] += c_hg->dist_x[i];

  /* Assuming that we do not collapse Edges, dist_y for the coarse hgraph
   * is the same as dist_y for the fine hgraph */
  /* EBEB We should remove hyperedges of size 1; then dist_y will change.  */
  if (!(c_hg->dist_y = (int*)ZOLTAN_MALLOC ((hgc->nProc_y+1) * sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }  
  for (i = 0; i < hgc->nProc_y+1; i++)
    c_hg->dist_y[i] = hg->dist_y[i];

  /* Done if there are no remaining vertices (on this proc) */
  if (c_hg->nVtx == 0)  {
    /* EBEB c_hg->vindex and vedge not properly initialized? */
    if (!(c_hg->vindex = (int*) ZOLTAN_CALLOC (1, sizeof(int))))  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    };
  }
  else  {
    c_hg->vindex  = c_vindex;
    c_hg->vedge   = c_vedge;
  }

  c_hg->vmap = NULL; /* UVC: we don't need vmap in the coarser graphs, it is only
                        needed in recursive bisection; and hence at level 0 */
  
  c_hg->hindex  = NULL;             /* is computed later by HG_Create_Mirror */
  c_hg->hvertex = NULL;             /* is computed later by HG_Create_Mirror */
  c_hg->coor    = NULL;             /* currently we don't use coordinates */
  c_hg->nDim    = hg->nDim;    
  c_hg->nEdge   = hg->nEdge;
  c_hg->comm    = hg->comm;
  c_hg->info    = hg->info + 1;     /* for debugging */
  c_hg->ratio   = hg->ratio;        /* for "global" recursive bisectioning */
  c_hg->redl    = hg->redl;         /* to stop coarsening near desired count */
    
  c_hg->VtxWeightDim  = hg->VtxWeightDim;
  c_hg->EdgeWeightDim = hg->EdgeWeightDim;
  
  Zoltan_Multifree (__FILE__, __LINE__, 7, &buffer, &rbuffer, &listgno, &listlno,
   &cmatch, &listproc, &msg_size);  
  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_HG_Create_Mirror(zz, c_hg);
}

#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
