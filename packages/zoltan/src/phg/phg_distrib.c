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

#include "phg_distrib.h"
    
#define _DEBUG1
#define _DEBUG2
#define _DEBUG3

#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

    
int Zoltan_PHG_Gno_To_Proc_Block(
  int gno,
  int *dist_dim,
  int nProc_dim
)
{
/* Function that locates a given global number gno within a distribution 
 * vector dist.
 * Works for both vtx and edges.
 * Takes an initial guess based on equal distribution of gno's.
 * Modifies guess based on actual distribution.
 */

int idx;
int maxgno = dist_dim[nProc_dim];

  idx = gno * nProc_dim / maxgno;

  while (gno < dist_dim[idx]) idx--;
  while (gno >= dist_dim[idx+1]) idx++;

  return idx;
}

int Zoltan_PHG_Redistribute(
  ZZ *zz, 
  HGraph  *ohg,           /* Input: Local part of distributed hypergraph */
  int     lo, int hi,     /* Input: range of proc ranks (inclusive)
                             to be included in new communicator: ncomm */
  PHGComm *ncomm,         /* Output: Communicators of new distribution */
  HGraph  *nhg,           /* Output: Newly redistributed hypergraph */
  int     **vmap          /* Output: allocated with the size nhg->nVtx and
                             vertex map from nhg to ohg's GLOBAL vertex number*/
    )   
{
    char * yo = "Zoltan_PHG_Redistribute";
    PHGComm *ocomm = ohg->comm;
    int     *v2Col, *n2Row, ierr=ZOLTAN_OK, tmp, i, *ranks;
    float   frac;
    MPI_Group allgrp, newgrp;

    if (ocomm->nProc==1){
        errexit("Zoltan_PHG_Redistribute: ocomm->nProc==1");
        return ZOLTAN_FATAL;
    }

    /* create a new communicator for procs[lo..hi] */
    MPI_Comm_group(ocomm->Communicator, &allgrp);
    ranks = (int *) ZOLTAN_MALLOC(ocomm->nProc * sizeof(int));
    for (i=lo; i<=hi; ++i)
        ranks[i-lo] = i;
    
    MPI_Group_incl(allgrp, hi-lo+1, ranks, &newgrp);
    MPI_Comm_create(ocomm->Communicator, newgrp, &ncomm->Communicator);
    MPI_Group_free(&newgrp);
    MPI_Group_free(&allgrp);   
    ZOLTAN_FREE(ranks);
    
    /* fill ncomm */
    ncomm->nProc = hi-lo+1;
    tmp = (int) sqrt((double)ncomm->nProc+0.1);
    while (ncomm->nProc % tmp)
        --tmp;
    ncomm->nProc_y = tmp;
    ncomm->nProc_x = ncomm->nProc / tmp;

    /* if new communicator is not NULL; this process is in that group
       so compute the rest of the stuff that it will need */
    if (ncomm->Communicator!=MPI_COMM_NULL) {
        MPI_Comm_rank(ncomm->Communicator, &ncomm->myProc);
#ifdef _DEBUG1
        if (ncomm->myProc != ocomm->myProc-lo)
            errexit("Zoltan_PHG_Redistribute: ncomm->myProc(%d) != ocomm->myProc(%d)-lo(%d)", ncomm->myProc, ocomm->myProc, lo);
#endif
        ncomm->myProc_x = ncomm->myProc % ncomm->nProc_x;
        ncomm->myProc_y = ncomm->myProc / ncomm->nProc_x;

        if ((MPI_Comm_split(ncomm->Communicator, ncomm->myProc_x, ncomm->myProc_y, 
                            &ncomm->col_comm) != MPI_SUCCESS)
            || (MPI_Comm_split(ncomm->Communicator, ncomm->myProc_y, ncomm->myProc_x, 
                               &ncomm->row_comm) != MPI_SUCCESS)) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "MPI_Comm_Split failed");
            return ZOLTAN_FATAL;
        }        
    } else {
        ncomm->myProc = ncomm->myProc_x = ncomm->myProc_y = -1;
    }
    
    v2Col = (int *) ZOLTAN_MALLOC(ohg->nVtx * sizeof(int));    
    n2Row = (int *) ZOLTAN_MALLOC(ohg->nEdge * sizeof(int));

    /* UVC: TODO very simple straight forward partitioning right now;
       later we can implemet a more "load balanced", or smarter
       mechanisms */
    frac = (float) ohg->nVtx / (float) ncomm->nProc_x;
    for (i=0; i<ohg->nVtx; ++i) 
        v2Col[i] = (int) ((float) i / frac);
    frac = (float) ohg->nEdge / (float) ncomm->nProc_y;
    for (i=0; i<ohg->nEdge; ++i) 
        n2Row[i] = (int) ((float) i / frac);
    
    ierr = Zoltan_PHG_Redistribute_Hypergraph(zz, ohg, v2Col, n2Row, ncomm, nhg, vmap);
    Zoltan_Multifree(__FILE__, __LINE__, 2,
                     &v2Col, &n2Row);
    
    return ierr;
}

    
    
int Zoltan_PHG_Redistribute_Hypergraph(
  ZZ *zz, 
  HGraph  *ohg,           /* Input:  Local part of distributed hypergraph */
  int     *v2Col,         /* Input:  Vertex to processor Column Mapping */
  int     *n2Row,         /* Input:  Net to processor Row Mapping */
  PHGComm *ncomm,         /* Input:  communicators of new distribution */
  HGraph  *nhg,           /* Output: Newly redistributed hypergraph */
  int     **vmap          /* Output: allocated with the size nhg->nVtx and
                             vertex map from nhg to ohg's GLOBAL vertex number*/
    )
{
    char * yo = "Zoltan_PHG_Redistribute_Hypergraph";
    PHGComm *ocomm = ohg->comm;
    int ierr=ZOLTAN_OK;
    int i, v, n, p=0, nsend, elemsz, nVtx, nEdge;
    int msg_tag = 90000;
    int *proclist=NULL, *sendbuf=NULL;
    int *vno=NULL, *nno=NULL, *dist_x=NULL, *dist_y=NULL,
        *vsn=NULL, *nsn=NULL, *pins=NULL, *cnt=NULL;
    ZOLTAN_COMM_OBJ *plan;    

    Zoltan_HG_HGraph_Init (nhg);
    nhg->comm = ncomm;
    /* UVCUVC: ADD memory checks for mallocs */
    nhg->dist_x = (int *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(int));
    nhg->dist_y = (int *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(int));
    dist_x = (int *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(int));
    dist_y = (int *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(int));
    vsn = (int *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(int));
    nsn = (int *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(int));
    vno = (int *) ZOLTAN_MALLOC(ohg->nVtx * sizeof(int));
    nno = (int *) ZOLTAN_MALLOC(ohg->nEdge * sizeof(int));

    if (!nhg->dist_x || !nhg->dist_y || !dist_x || !dist_y ||
        !vsn || !nsn || !vno || !nno )
        MEMORY_ERROR;
      
    for (v = 0; v < ohg->nVtx; ++v)
        ++dist_x[v2Col[v]];
    for (n = 0; n < ohg->nEdge; ++n)
        ++dist_y[n2Row[n]];

    /* compute prefix sum to find new vertex start numbers; for each processor */
    MPI_Scan(dist_x, vsn, ncomm->nProc_x, MPI_INT, MPI_SUM, ocomm->row_comm);
    /* All reduce to compute how many each processor will have */ 
    MPI_Allreduce(dist_x, &(nhg->dist_x[1]), ncomm->nProc_x, MPI_INT, MPI_SUM, ocomm->row_comm);
    nhg->dist_x[0] = 0;
    
    for (i=1; i<=ncomm->nProc_x; ++i) 
        nhg->dist_x[i] += nhg->dist_x[i-1];
    MPI_Scan(dist_y, nsn, ncomm->nProc_y, MPI_INT, MPI_SUM, ocomm->col_comm);
    MPI_Allreduce(dist_y, &(nhg->dist_y[1]), ncomm->nProc_y, MPI_INT, MPI_SUM, ocomm->col_comm);
    nhg->dist_y[0] = 0;
    for (i=1; i<=ncomm->nProc_y; ++i)
        nhg->dist_y[i] += nhg->dist_y[i-1];

    /* find mapping of current LOCAL vertex no (in my node)
       to "new" vertex no LOCAL to dest node*/
    for (v = ohg->nVtx-1; v>=0; --v)
        vno[v] = --vsn[v2Col[v]];
    for (n = ohg->nEdge-1; n>=0; --n)
        nno[n] = --nsn[n2Row[n]];

    nsend = MAX(MAX(ohg->nPins, ohg->nVtx), ohg->nEdge);
    elemsz = MAX(MAX(2, ohg->VtxWeightDim), ohg->EdgeWeightDim);
    
    proclist = (int *) ZOLTAN_MALLOC(nsend * sizeof(int));
    sendbuf = (int *) ZOLTAN_MALLOC(nsend * elemsz * sizeof(int));
    if (sizeof(float)>sizeof(int))
        errexit("Zoltan_PHG_Redistribute_Hypergraph: this code assumes sizeof(float)(%d)<=sizeof(int)(%d)", sizeof(float), sizeof(int));


    /* first communicate pins */    
    for (v = 0; v < ohg->nVtx; ++v) { 
        for (i = ohg->vindex[v]; i < ohg->vindex[v+1]; ++i) {
            proclist[p]   = n2Row[ohg->vedge[i]] * ncomm->nProc_x + v2Col[v];
            sendbuf[2*p]  = vno[v];
            sendbuf[2*p+1]= nno[ohg->vedge[i]];
            ++p; 
        }
    }
#ifdef _DEBUG1
    if (p!=ohg->nPins) {
        uprintf(ocomm, "sanity check failed p(%d)!=hg->nPins(%d)\n", p, ohg->nPins);
        errexit("terminating");
    }
#endif

    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nPins, proclist, ocomm->Communicator,
                               msg_tag, &p);

#ifdef _DEBUG1
    if (ncomm->myProc==-1 && p>1) { /* this processor is not in new comm but receiving data?*/
        uprintf(ocomm, "Something wrong; why I'm receiving data p=%d\n", p);
        errexit("terminating");
    }
#endif
    
    if (p && (pins = (int *) ZOLTAN_MALLOC(p * 2 * sizeof(int)))==NULL) 
        MEMORY_ERROR;

    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                   (char *) pins);

    Zoltan_Comm_Destroy(&plan);

#ifdef _DEBUG1
    /* now communicate local vertex no's wrt dest proc; for debugging */
    for (v = 0; v < ohg->nVtx; ++v) {
        proclist[v] = v2Col[v];
        sendbuf[v] = vno[v];
    }
    
    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nVtx, proclist, ocomm->row_comm,
                               msg_tag, &nVtx);

    if (nVtx && (cnt = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int)))==NULL) 
        MEMORY_ERROR;
    
    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                   (char *) cnt);
    

    if (ncomm->myProc==-1 && nVtx>1) { /* this processor is not in new comm but receiving data?*/
        uprintf(ocomm, "Something wrong; why I'm receiving data nVtx=%d\n", nVtx);
        errexit("terminating");
    } else if (ncomm->myProc!=-1) {
        for (v=0; v<nVtx; ++v)
            if (cnt[v]!=v)
                errexit("numbering error; for vertex %d received %d", v, cnt[v]);
    }
#endif 
    /* now communicate vertex map */ 
    for (v = 0; v < ohg->nVtx; ++v) { 
        proclist[v] = v2Col[v];
        sendbuf[v] = VTX_LNO_TO_GNO(ohg, v);
    }
    
    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nVtx, proclist, ocomm->row_comm,
                               msg_tag, &nVtx);

#ifdef _DEBUG1
    if (ncomm->myProc==-1 && nVtx>1) { /* this processor is not in new comm but receiving data?*/
        uprintf(ocomm, "Something wrong; why I'm receiving data nVtx=%d\n", nVtx);
        errexit("terminating");
    }
#endif
    
    if (nVtx && (*vmap = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int)))==NULL) 
        MEMORY_ERROR;
    
    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                   (char *) *vmap);


    /* now communicate vertex weights */
    if (nVtx && ohg->vwgt && ohg->VtxWeightDim)
        nhg->vwgt = (float*) ZOLTAN_MALLOC(nVtx*ohg->VtxWeightDim*sizeof(float));
    
    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->vwgt, ohg->VtxWeightDim*sizeof(float),
                   (char *) nhg->vwgt);
    
    Zoltan_Comm_Destroy(&plan);
    

    /* now communicate edge weights */ 
    for (n = 0; n < ohg->nEdge; ++n) 
        proclist[v] = n2Row[n];
    
    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nEdge, proclist, ocomm->col_comm,
                               msg_tag, &nEdge);

#ifdef _DEBUG1
    if (ncomm->myProc==-1 && nEdge>1) { /* this processor is not in new comm but receiving data?*/
        uprintf(ocomm, "Something wrong; why I'm receiving data nEdge=%d\n", nEdge);
        errexit("terminating");
    }
#endif
    
    if (nEdge && ohg->ewgt && ohg->EdgeWeightDim)
        nhg->ewgt = (float*) ZOLTAN_MALLOC(nEdge*ohg->EdgeWeightDim*sizeof(float));
    
    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->ewgt, ohg->EdgeWeightDim*sizeof(float),
                   (char *) nhg->ewgt);    
    Zoltan_Comm_Destroy(&plan);
    


    if (ncomm->myProc==-1) {
        nhg->nEdge = nhg->nVtx = nhg->nPins = 0;
    } else {
        nhg->nEdge = nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y];
        nhg->nVtx = nhg->dist_x[ncomm->myProc_x+1] - nhg->dist_x[ncomm->myProc_x];
        nhg->nPins = p;
    
        /* Unpack the pins received. */
        cnt = (int *) ZOLTAN_CALLOC(nhg->nVtx + 1, sizeof(int));
        nhg->vindex = (int *) ZOLTAN_CALLOC(nhg->nVtx + 1, sizeof(int));
        nhg->vedge = (int *) ZOLTAN_MALLOC(p * sizeof(int));
        
        if (!cnt || !nhg->vindex || (p && !nhg->vedge))
            MEMORY_ERROR;
        
        /* Count the number of pins per vertex */
        for (i = 0; i < p; ++i)
            ++cnt[pins[2*i]];
        
        
        /* Compute prefix sum to represent hindex correctly. */
        for (i = 0; i < nhg->nVtx; ++i)  {
            nhg->vindex[i+1] = nhg->vindex[i] + cnt[i];
            cnt[i] = nhg->vindex[i];
        }
        
        for (i = 0; i < p; ++i) 
            nhg->vedge[cnt[pins[2*i]++]] = pins[2*i+1];
        
        nhg->info               = ohg->info;
        nhg->VtxWeightDim       = ohg->VtxWeightDim;
        nhg->EdgeWeightDim      = ohg->EdgeWeightDim;
        
        
        ierr = Zoltan_HG_Create_Mirror(zz, nhg);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
            MEMORY_ERROR;
    }

 End:
    Zoltan_Multifree(__FILE__, __LINE__, 10,
                     &proclist, &sendbuf, &pins, &cnt,
                     &vno, &nno, &dist_x, &dist_y, &vsn, &nsn
        );
    
    return ierr;
}
    


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
