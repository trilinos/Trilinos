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

    
    
int Zoltan_PHG_Redistribute_Hypergraph(
  ZZ *zz, 
  HGraph  *ohg,           /* Input:  Local part of distributed hypergraph */
  int     *v2Col,         /* Input:  Vertex to processor Column Mapping */
  int     *n2Row,         /* Input:  Net to processor Row Mapping */
  PHGComm *ncomm,         /* Input:  communicators of new distribution */
  HGraph  *nhg            /* Output: Newly redistributed hypergraph */ 
    )
{
    char * yo = "Zoltan_PHG_Redistribute_Hypergraph";
    PHGComm *ocomm = ohg->comm;
    int ierr=ZOLTAN_OK;
    int i, v, n, p=0;
    int msg_tag = 90000;
    int *proclist=NULL, *sendbuf=NULL;
    int *vno=NULL, *nno=NULL, *dist_x=NULL, *dist_y=NULL,
        *vsn=NULL, *nsn=NULL, *pins=NULL, *cnt;
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

    if (!nhg->dist_x || !nhg->dist_y || !vno || !nno) MEMORY_ERROR;
      
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

#ifdef _DEBUG1    
    for (i=0; i<=ncomm->nProc_x; ++i)
        if (vsn[i])
            errexit("hey vsn[%d]=%d", i, vsn[i]); 
    for (i=0; i<=ncomm->nProc_y; ++i)
        if (nsn[i])
            errexit("hey nsn[%d]=%d", i, nsn[i]);   
#endif

    proclist = (int *) ZOLTAN_MALLOC(ohg->nPins * sizeof(int));
    sendbuf = (int *) ZOLTAN_MALLOC(ohg->nPins * 2 * sizeof(int));

    for (v = 0; v < ohg->nVtx; ++v) {
        for (i = ohg->vindex[v]; i < ohg->vindex[v+1]; ++i) {
            proclist[p]   = n2Row[ohg->vedge[i]] * ncomm->nProc_x + v2Col[v];
            sendbuf[2*p]  = vno[v];
            sendbuf[2*p+1]= nno[ohg->vedge[i]];
            ++p;
        }
    }

    if (p!=ohg->nPins) {
        uprintf(ocomm, "sanity check failed p(%d)!=hg->nPins(%d)\n", p, ohg->nPins);
        ierr = ZOLTAN_FATAL;
        goto End;
    }
    
    --msg_tag;
    ierr |= Zoltan_Comm_Create(&plan, ohg->nPins, proclist, ocomm->Communicator,
                               msg_tag, &p);
    
    if (p && (pins = (int *) ZOLTAN_MALLOC(p * 2 * sizeof(int)))==NULL) 
        MEMORY_ERROR;
    

    --msg_tag;
    Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                   (char *) pins);

    Zoltan_Comm_Destroy(&plan);
    
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

  
    ierr = Zoltan_HG_Create_Mirror(zz, nhg);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
        MEMORY_ERROR;

    /* UVCUVC TODO CHECK: communicate vertex&net weights!!!
       make sure to send old vertex # so that partitioning
       results can be sent to appropriate node */

 End:
    Zoltan_Multifree(__FILE__, __LINE__, 10,
                     &proclist, &sendbuf, &pins, &cnt,
                     &vno, &nno, &dist_x, &dist_y, &vsn, &nsn
        );

    return ierr;
}
    

int Zoltan_PHG_Redistribute(
  ZZ *zz, 
  HGraph  *ohg,           /* Input:  Local part of distributed hypergraph */
  PHGComm *ncomm,         /* Output: Communicators of new distribution */
  HGraph  *nhg            /* Output: Newly redistributed hypergraph */
    )
    
{
    char * yo = "Zoltan_PHG_Redistribute";
    PHGComm *ocomm = ohg->comm;
    int     *v2Col, *n2Row, ierr=ZOLTAN_OK, tmp, i;
    float   frac;

    if (ocomm->nProc==1){
        errexit("Zoltan_PHG_Redistribute: ocomm->nProc==1");
        return ZOLTAN_FATAL;
    }

    /* UVCUVC: TODO BUGBUG; we probably don't want to divide machine from middle??
       we should probably doing something similar that we do in rdivide */
    ncomm->nProc = ocomm->nProc/2;
    ncomm->myProc = ocomm->myProc; /* this is clearly wrong; FIX IT! */
    ncomm->Communicator = ocomm->Communicator; /* so does this line */

       /* Compute default */
    tmp = (int) sqrt((double)ncomm->nProc+0.1);
    while (ncomm->nProc % tmp)
        --tmp;
    ncomm->nProc_y = tmp;
    ncomm->nProc_x = ncomm->nProc / tmp;
    
    ncomm->myProc_x = ncomm->myProc % ncomm->nProc_x;
    ncomm->myProc_y = ncomm->myProc / ncomm->nProc_x;

    
    if ((MPI_Comm_split(ncomm->Communicator, ncomm->myProc_x, ncomm->myProc_y, 
                        &ncomm->col_comm) != MPI_SUCCESS)
        || (MPI_Comm_split(ncomm->Communicator, ncomm->myProc_y, ncomm->myProc_x, 
                           &ncomm->row_comm) != MPI_SUCCESS)) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "MPI_Comm_Split failed");
        return ZOLTAN_FATAL;
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
    
    ierr = Zoltan_PHG_Redistribute_Hypergraph(zz, ohg, v2Col, n2Row, ncomm, nhg);
    Zoltan_Multifree(__FILE__, __LINE__, 2,
                     &v2Col, &n2Row);
    
    return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
