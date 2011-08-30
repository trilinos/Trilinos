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
#include "zz_util_const.h"

    /*
#define _DEBUG1
#define _DEBUG2
#define _DEBUG3
    */

    
int Zoltan_PHG_Gno_To_Proc_Block(
  ZOLTAN_GNO_TYPE gno,
  ZOLTAN_GNO_TYPE *dist_dim,
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
double fidx;
ZOLTAN_GNO_TYPE maxgno = dist_dim[nProc_dim];

  fidx = (double)gno * (double)nProc_dim / (double)maxgno;
  idx = (int) fidx;

  while (gno < dist_dim[idx]) idx--;
  while (gno >= dist_dim[idx+1]) idx++;

  return idx;
}



#ifdef _DEBUG1
static void PrintArr(PHGComm *hgc, char *st, int *ar, int n)
{
    int i;
    Zoltan_Print_Sync_Start(hgc->Communicator, TRUE);
    uprintf(hgc, "%s = [", st);
    for (i=0; i<n; ++i)
        printf("%d, ", ar[i]);
    printf("]\n");
    Zoltan_Print_Sync_End(hgc->Communicator, TRUE);
}
#endif


static int Zoltan_PHG_Redistribute_Hypergraph(
    ZZ *zz, 
    PHGPartParams *hgp,     /* Input:  parameters; used only for UseFixedVtx */
    HGraph  *ohg,           /* Input:  Local part of distributed hypergraph */
    int     firstproc,      /* Input:  rank (in ocomm) of the first proc of 
                                       the ncomm*/
    int     *v2Col,         /* Input:  Vertex to processor Column Mapping */
    int     *n2Row,         /* Input:  Net to processor Row Mapping */
    PHGComm *ncomm,         /* Input:  communicators of new distribution */
    HGraph  *nhg,           /* Output: Newly redistributed hypergraph */
    int     **vmap,         /* Output: allocated with the size nhg->nVtx and
                               vertex map from nhg to ohg's local vertex number*/
    int     **vdest         /* Output: allocated with the size nhg->nVtx and
                               stores dest proc in ocomm */
    )
{
    char * yo = "Zoltan_PHG_Redistribute_Hypergraph";
    PHGComm *ocomm = ohg->comm;
    int ierr=ZOLTAN_OK;
    int i, v, n, nPins, nsend, elemsz, nVtx, nEdge;
    int msg_tag = 9999;
    int *proclist=NULL, *cnt=NULL; 
    int *intbuf;
    ZOLTAN_GNO_TYPE *sendbuf=NULL;
    ZOLTAN_GNO_TYPE *vno=NULL, *nno=NULL, *dist_x=NULL, *dist_y=NULL,
        *vsn=NULL, *nsn=NULL, *pins=NULL;
    ZOLTAN_COMM_OBJ *plan;    
    MPI_Datatype zoltan_gno_mpi_type;

    zoltan_gno_mpi_type = Zoltan_mpi_gno_type();
    
    Zoltan_HG_HGraph_Init (nhg);
    nhg->comm = ncomm;
    
    nhg->dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(ZOLTAN_GNO_TYPE));
    nhg->dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(ZOLTAN_GNO_TYPE));
    dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(ZOLTAN_GNO_TYPE));
    dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(ZOLTAN_GNO_TYPE));
    vsn = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_x+1, sizeof(ZOLTAN_GNO_TYPE));
    nsn = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(ncomm->nProc_y+1, sizeof(ZOLTAN_GNO_TYPE));
    vno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(ohg->nVtx * sizeof(ZOLTAN_GNO_TYPE));
    nno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(ohg->nEdge * sizeof(ZOLTAN_GNO_TYPE));

    if (!nhg->dist_x || !nhg->dist_y || !dist_x || !dist_y ||
        !vsn || !nsn || (ohg->nVtx && !vno) || (ohg->nEdge && !nno) ) {
        uprintf(ocomm, " new comm nProcx=%d nProcy=%d nvtx=%d nedge=%d", ncomm->nProc_x, ncomm->nProc_y, ohg->nVtx, ohg->nEdge);
        MEMORY_ERROR;
    }
      
    for (v = 0; v < ohg->nVtx; ++v)
        ++dist_x[v2Col[v]];
    for (n = 0; n < ohg->nEdge; ++n)
        ++dist_y[n2Row[n]];

    /* UVCUVC: CHECK ASSUMPTION
       This code assumes that the objects in the receive buffer of
       Zoltan_Comm_Do function are
         1- in the increasing processor order,
         2- order of the items send by a processor is preserved.
     */
    

    /* compute prefix sum to find new vertex start numbers; for each processor */
    MPI_Scan(dist_x, vsn, ncomm->nProc_x, zoltan_gno_mpi_type, MPI_SUM, ocomm->row_comm);
    /* All reduce to compute how many each processor will have */ 
    MPI_Allreduce(dist_x, &(nhg->dist_x[1]), ncomm->nProc_x, zoltan_gno_mpi_type, MPI_SUM, 
                  ocomm->row_comm);
    nhg->dist_x[0] = 0;    
    for (i=1; i<=ncomm->nProc_x; ++i) 
        nhg->dist_x[i] += nhg->dist_x[i-1];
    
    MPI_Scan(dist_y, nsn, ncomm->nProc_y, zoltan_gno_mpi_type, MPI_SUM, ocomm->col_comm);

    MPI_Allreduce(dist_y, &(nhg->dist_y[1]), ncomm->nProc_y, zoltan_gno_mpi_type, MPI_SUM, ocomm->col_comm);
    nhg->dist_y[0] = 0;
    for (i=1; i<=ncomm->nProc_y; ++i)
        nhg->dist_y[i] += nhg->dist_y[i-1];

#ifdef _DEBUG1
    PrintArr(ocomm, "vsn", vsn, ncomm->nProc_x);
    PrintArr(ocomm, "nsn", nsn, ncomm->nProc_y);
#endif
    
    /* find mapping of current LOCAL vertex no (in my node)
       to "new" vertex no LOCAL to dest node*/
    for (v = ohg->nVtx-1; v>=0; --v)
        vno[v] = --vsn[v2Col[v]];
    for (n = ohg->nEdge-1; n>=0; --n)
        nno[n] = --nsn[n2Row[n]];

    nsend = MAX(MAX(ohg->nPins, ohg->nVtx), ohg->nEdge);
    elemsz = MAX(MAX(2, ohg->VtxWeightDim), ohg->EdgeWeightDim);
    elemsz = (sizeof(float)>sizeof(ZOLTAN_GNO_TYPE)) ? sizeof(float)*elemsz : sizeof(ZOLTAN_GNO_TYPE)*elemsz;

    proclist = (int *) ZOLTAN_MALLOC(nsend * sizeof(int));
    sendbuf = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nsend * elemsz);

    if (nsend && (!proclist || !sendbuf)) MEMORY_ERROR;

    /* first communicate pins */
    nPins = 0;
    for (v = 0; v < ohg->nVtx; ++v) { 
        for (i = ohg->vindex[v]; i < ohg->vindex[v+1]; ++i) {
#ifdef _DEBUG1
            if ((n2Row[ohg->vedge[i]] * ncomm->nProc_x + v2Col[v])<0 ||
                (n2Row[ohg->vedge[i]] * ncomm->nProc_x + v2Col[v])>=ocomm->nProc)
                errexit("vertex %d vedge[%d]=%d n2Row=%d #Proc_x=%d v2Col=%d", i, ohg->vedge[i], n2Row[ohg->vedge[i]], ncomm->nProc_x , v2Col[v]);
#endif
            proclist[nPins]   = firstproc + n2Row[ohg->vedge[i]] * ncomm->nProc_x + v2Col[v];
            sendbuf[2*nPins]  = vno[v];
            sendbuf[2*nPins+1]= nno[ohg->vedge[i]];
            ++nPins; 
        }
    }
#ifdef _DEBUG1
    if (nPins!=ohg->nPins) {
        uprintf(ocomm, "sanity check failed nPins(%d)!=hg->nPins(%d)\n", nPins, ohg->nPins);
        errexit("terminating");
    }
#endif

    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nPins, proclist, ocomm->Communicator, msg_tag, &nPins);

#ifdef _DEBUG1
    if (ncomm->myProc==-1 && nPins>1) { /* this processor is not in new comm but receiving data?*/
        uprintf(ocomm, "Something wrong; why I'm receiving data nPins=%d\n", nPins);
        errexit("terminating");
    }
#endif
    
    if (nPins && (pins = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nPins * 2 * sizeof(ZOLTAN_GNO_TYPE)))==NULL) 
        MEMORY_ERROR;

    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(ZOLTAN_GNO_TYPE), (char *) pins);
    Zoltan_Comm_Destroy(&plan);

    /* now communicate vertex map */
    nsend = 0;
    intbuf = (int *)sendbuf;
    if (!ocomm->myProc_y) { /* only first row sends to the first row of ncomm */
        for (v = 0; v < ohg->nVtx; ++v) { 
            proclist[nsend] = firstproc+v2Col[v];
            intbuf[nsend++] = ohg->vmap[v];
        }
    }
        
    --msg_tag; 
    ierr |= Zoltan_Comm_Create(&plan, nsend, proclist, ocomm->Communicator, msg_tag, &nVtx); 

#ifdef _DEBUG1
    if (ncomm->myProc==-1 && nVtx>1) { /* this processor is not in new comm but receiving data?*/ 
        uprintf(ocomm, "Something wrong; why I'm receiving data nVtx=%d\n", nVtx);
        errexit("terminating");
    }
#endif

    /* those are only needed in the first row of ncomm */
    *vmap = *vdest = NULL;  
    if (!ncomm->myProc_y && nVtx &&
        (!(*vmap = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int))) ||
         !(*vdest = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int)))))
        MEMORY_ERROR;
    

    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, sizeof(int), (char *) *vmap);

    intbuf = (int *)sendbuf;
    if (!ocomm->myProc_y) { /* only first row sends to the first row of ncomm */
        for (v = 0; v < ohg->nVtx; ++v) 
            intbuf[v] = ocomm->myProc;
    }
    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, sizeof(int), (char *) *vdest);
        
    if (ncomm->myProc!=-1) { /* I'm in the new comm */
        /* ncomm's first row now bcast to other rows */
        MPI_Bcast(&nVtx, 1, MPI_INT, 0, ncomm->col_comm);
#ifdef _DEBUG1
        if (nVtx!=(int)(nhg->dist_x[ncomm->myProc_x+1] - nhg->dist_x[ncomm->myProc_x]))
            errexit("nVtx(%d)!= nhg->dist_x[ncomm->myProc_x+1] - nhg->dist_x[ncomm->myProc_x](%d)", nVtx, nhg->dist_x[ncomm->myProc_x+1] - nhg->dist_x[ncomm->myProc_x]);
#endif
        if (nVtx && (nhg->vmap = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int)))==NULL) 
            MEMORY_ERROR;
        for (i=0; i<nVtx; ++i)
            nhg->vmap[i] = i;
    }


    /* now communicate vertex weights */
    if (ohg->VtxWeightDim) {
        if (nVtx){
            nhg->vwgt = (float*) ZOLTAN_MALLOC(nVtx*ohg->VtxWeightDim*sizeof(float));
            if (!nhg->vwgt) MEMORY_ERROR;
        }
    
        --msg_tag;
        Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->vwgt,
                       ohg->VtxWeightDim*sizeof(float), (char *) nhg->vwgt);
        if (ncomm->myProc!=-1)  /* ncomm's first row now bcast to other rows */
            MPI_Bcast(nhg->vwgt, nVtx*ohg->VtxWeightDim, MPI_FLOAT, 0, ncomm->col_comm);
    }

    /* communicate coordinates */
    if (ohg->nDim) {
      if (nVtx) {
	nhg->coor = (double *) ZOLTAN_MALLOC(nVtx * ohg->nDim * sizeof(double));
	if (!nhg->coor) MEMORY_ERROR;
      }

      --msg_tag;
      Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->coor, ohg->nDim * sizeof(double),
		     (char *) nhg->coor);
      if (ncomm->myProc != -1) /* ncomm's first row bcast to other rows */
	MPI_Bcast(nhg->coor, nVtx * ohg->nDim, MPI_DOUBLE, 0, ncomm->col_comm);
    }
    
    /* communicate fixed vertices, if any */
    if (hgp->UseFixedVtx) {
        if (nVtx){
            nhg->fixed_part = (int *) ZOLTAN_MALLOC(nVtx*sizeof(int));
            if (!nhg->fixed_part) MEMORY_ERROR;
        }
        --msg_tag;
        Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->fixed_part,
                       sizeof(int), (char *) nhg->fixed_part);
        if (ncomm->myProc!=-1)  /* ncomm's first row now bcast to other rows */
            MPI_Bcast(nhg->fixed_part, nVtx, MPI_INT, 0, ncomm->col_comm);
    }    
    /* communicate pref parts, if any */
    if (hgp->UsePrefPart) {
        if (nVtx){
            nhg->pref_part = (int *) ZOLTAN_MALLOC(nVtx*sizeof(int));
            if (!nhg->pref_part) MEMORY_ERROR;
        }
        --msg_tag;
        Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->pref_part,
                       sizeof(int), (char *) nhg->pref_part);
        if (ncomm->myProc!=-1)  /* ncomm's first row now bcast to other rows */
            MPI_Bcast(nhg->pref_part, nVtx, MPI_INT, 0, ncomm->col_comm);
    }    

    /* this comm plan is no longer needed. */
    Zoltan_Comm_Destroy(&plan);

    
    if (ohg->EdgeWeightDim) { /* now communicate edge weights */
        nsend = 0;
        if (!ocomm->myProc_x)  /* only first column sends to first column of ncomm */
            for (n = 0; n < ohg->nEdge; ++n) 
                proclist[nsend++] = firstproc + n2Row[n]*ncomm->nProc_x;
    
        --msg_tag;
        ierr |= Zoltan_Comm_Create(&plan, nsend, proclist, ocomm->Communicator,
                                   msg_tag, &nEdge);

#ifdef _DEBUG1
        if (ncomm->myProc==-1 && nEdge>1) { /* this processor is not in new comm but receiving data?*/
            uprintf(ocomm, "Something wrong; why I'm receiving data nEdge=%d\n", nEdge);
            errexit("terminating");
        }
#endif
        if (ncomm->myProc!=-1) { /* if we're in the new comm */
            /* ncomm's first column now bcast to other columns */
            MPI_Bcast(&nEdge, 1, MPI_INT, 0, ncomm->row_comm);
#ifdef _DEBUG1
            if (nEdge != (nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y]))
            errexit("nEdge(%d)!=nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y](%d)", nEdge, nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y]);
#endif
        }
        
        if (nEdge){
            nhg->ewgt = (float*) ZOLTAN_MALLOC(nEdge*ohg->EdgeWeightDim*sizeof(float));
            if (!nhg->ewgt) MEMORY_ERROR;
        }
    
        --msg_tag;
        Zoltan_Comm_Do(plan, msg_tag, (char *) ohg->ewgt,
                       ohg->EdgeWeightDim*sizeof(float), (char *) nhg->ewgt);
        if (ncomm->myProc!=-1) { /* if we're in the new comm */
            /* ncomm's first column now bcast to other columns */
            if (nEdge) 
                MPI_Bcast(nhg->ewgt, nEdge*ohg->EdgeWeightDim, MPI_FLOAT, 0, 
                          ncomm->row_comm);
        }

        Zoltan_Comm_Destroy(&plan);
    } else 
        nEdge = (ncomm->myProc==-1) 
                ? 0 
                : (int)(nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y]);
    

    if (ncomm->myProc==-1) {
#ifdef _DEBUG1
        if (nPins || nVtx || nEdge)
            errexit("I should not have any data: hey nPins=%d  nVtx=%d  nEdge=%d\n", nPins, nVtx, nEdge);
#endif
        nhg->nEdge = nhg->nVtx = nhg->nPins = 0;
    } else {
        nhg->nEdge = (int)(nhg->dist_y[ncomm->myProc_y+1] - nhg->dist_y[ncomm->myProc_y]);
        nhg->nVtx = (int)(nhg->dist_x[ncomm->myProc_x+1] - nhg->dist_x[ncomm->myProc_x]);
        nhg->nPins = nPins;
    
        /* Unpack the pins received. */
        cnt = (int *) ZOLTAN_CALLOC(nhg->nVtx + 1, sizeof(int));
        nhg->vindex = (int *) ZOLTAN_CALLOC(nhg->nVtx + 1, sizeof(int));
        nhg->vedge = (int *) ZOLTAN_MALLOC(nhg->nPins * sizeof(int));
        
        if (!cnt || !nhg->vindex || (nPins && !nhg->vedge))
            MEMORY_ERROR;

        /* Count the number of pins per vertex */
        for (i = 0; i < nPins; ++i)
            ++cnt[pins[2*i]];
        
        /* Compute prefix sum to represent hindex correctly. */
        for (i = 0; i < nhg->nVtx; ++i)  {
            nhg->vindex[i+1] = nhg->vindex[i] + cnt[i];
            cnt[i] = nhg->vindex[i];
        }

        for (i = 0; i < nPins; ++i) 
            nhg->vedge[cnt[pins[2*i]]++] = pins[2*i+1];
        
        nhg->info               = ohg->info;
	nhg->nDim               = ohg->nDim;
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
    


/**********************************************************************/
int Zoltan_PHG_Redistribute(
  ZZ *zz, 
  PHGPartParams *hgp,     /* Input: parameters; used only for user's
                             request of nProc_x and nProc_y */
  HGraph  *ohg,           /* Input: Local part of distributed hypergraph */
  int     lo, int hi,     /* Input: range of proc ranks (inclusive)
                             to be included in new communicator: ncomm */
  PHGComm *ncomm,         /* Output: Communicators of new distribution */
  HGraph  *nhg,           /* Output: Newly redistributed hypergraph */
  int     **vmap,         /* Output: allocated with the size nhg->nVtx and
                             vertex map from nhg to ohg's local vertex number*/
  int     **vdest         /* Output: allocated with the size nhg->nVtx and
                             stores dest proc in ocomm */
    )   
{
    char * yo = "Zoltan_PHG_Redistribute";
    PHGComm *ocomm = ohg->comm;
    int     *v2Col, *n2Row, ierr=ZOLTAN_OK, i, *ranks;
    int     reqx=hgp->nProc_x_req, reqy=hgp->nProc_y_req;
    float   frac;
    MPI_Group allgrp, newgrp;
    MPI_Comm  nmpicomm;

    if (ocomm->nProc==1){
        errexit("%s: ocomm->nProc==1", yo);
        return ZOLTAN_FATAL;
    }

    /* create a new communicator for procs[lo..hi] */
    MPI_Comm_group(ocomm->Communicator, &allgrp);
    ranks = (int *) ZOLTAN_MALLOC(ocomm->nProc * sizeof(int));
    if (!ranks) MEMORY_ERROR;

    for (i=lo; i<=hi; ++i)
        ranks[i-lo] = i;
    
    MPI_Group_incl(allgrp, hi-lo+1, ranks, &newgrp);
    MPI_Comm_create(ocomm->Communicator, newgrp, &nmpicomm);
    MPI_Group_free(&newgrp);
    MPI_Group_free(&allgrp);   
    ZOLTAN_FREE(&ranks);

    if (reqx==1 || reqy==1)
        ;
    else
        reqx = reqy = -1;
    
    /* fill ncomm */
    ierr = Zoltan_PHG_Set_2D_Proc_Distrib(ocomm->zz, nmpicomm, 
                                          ocomm->myProc-lo, hi-lo+1, 
                                          reqx, reqy, ncomm);
    
    v2Col = (int *) ZOLTAN_MALLOC(ohg->nVtx * sizeof(int));    
    n2Row = (int *) ZOLTAN_MALLOC(ohg->nEdge * sizeof(int));

    if ( (ohg->nVtx && !v2Col) || (ohg->nEdge && !n2Row)) MEMORY_ERROR;

    /* UVC: TODO very simple straight forward partitioning right now;
       later we can implement a more "load balanced", or smarter
       mechanisms */
    /* KDDKDD 5/11/07:  Round-off error in the computation of v2Col
     * and n2Row can lead to different answers on different platforms.
     * Vertices or edges get sent to different processors during the 
     * split, resulting in different matchings and, thus, different
     * answers.
     * Problem was observed on hg_cage10, zdrive.inp.phg.ipm.nproc_vertex1
     * and zdrive.inp.phg.ipm.nproc_edge1;
     * solaris machine seamus and linux machine patches give different
     * results due to differences in n2Row and v2Col, respectively.  
     * Neither answer is wrong,
     * but the linux results result in FAILED test in test_zoltan.
     */
    frac = (float) ohg->nVtx / (float) ncomm->nProc_x;
    for (i=0; i<ohg->nVtx; ++i) 
        v2Col[i] = (int) ((float) i / frac);
    frac = (float) ohg->nEdge / (float) ncomm->nProc_y;
    for (i=0; i<ohg->nEdge; ++i) 
        n2Row[i] = (int) ((float) i / frac);

    ierr |= Zoltan_PHG_Redistribute_Hypergraph(zz, hgp, ohg, lo, 
                                               v2Col, n2Row, ncomm, 
                                               nhg, vmap, vdest);
    Zoltan_Multifree(__FILE__, __LINE__, 2,
                     &v2Col, &n2Row);

End:
    
    return ierr;
}

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
