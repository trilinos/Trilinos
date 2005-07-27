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


#include <limits.h>
#include <ctype.h>
#include "zoltan_mem.h"    
#include "zz_const.h"
#include "coloring.h"
#include "params_const.h"
#include "zz_util_const.h"
#include "parmetis_jostle.h"

#define SWAP(a,b) tmp=(a);(a)=(b);(b)=tmp;

#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}
    
/* function prototypes */
int ReorderGraph(ZZ *, int, int **, int **, int **, int *, int *, int *, int *);
int pickColor(ZZ *, char, int, int, int *, int *);
int IntColoring(ZZ *, int *, int, int, int *, int *, int *, int *, int *, int *, char);
int ParallelColoring(ZZ *, int, int *, int *, int *, int *, int *, int, int *, int *, int **, int *, char, char, int *, int *, MPI_Request *, MPI_Request *, MPI_Status *);
int Conflict(ZZ *, int, int *, int *, int *, int *, int *, int *, int *, int *);
int D1coloring(ZZ *, char, char, char, int, int, int, int *, int *, int **, int **, int *, int *);
    
/*****************************************************************************/
/*  Parameters structure for Color method.  Used in  */
/*  Zoltan_Color_Set_Param and Zoltan_Color.         */
static PARAM_VARS Color_params[] = {
                  { "DISTANCE", NULL, "INT", 0 },
                  { "SUPERSTEP_SIZE", NULL, "INT", 0},
                  { "COMM_PATTERN", NULL, "CHAR", 0 },
                  { "COLOR_ORDER", NULL, "CHAR", 0 },
                  { "COLORING_METHOD", NULL, "CHAR", 0},
                  { NULL, NULL, NULL, 0 } };
    
/*---------------------------------------------------------------------------*/

int Zoltan_Color_Set_Param(
char *name,			/* name of variable */
char *val			/* value of variable */
)
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, Color_params, &result, &index);

    return(status);
}

    
/**********************************************************/
/* Interface routine for Graph Coloring                   */
/**********************************************************/    
int Zoltan_Color(
    ZZ *zz,                   /* Zoltan structure */
    int *num_gid_entries,     /* # of entries for a global id */
    int *num_lid_entries,     /* # of entries for a local id */
    int num_obj,              /* Input: number of objects */
    ZOLTAN_ID_PTR global_ids, /* Input: global ids of the vertices */
                              /* The application must allocate enough space */    
    ZOLTAN_ID_PTR local_ids,  /* Input: local ids of the vertices */
                              /* The application must allocate enough space */
    int *color_exp            /* Output: Colors assigned to local vertices */
                              /* The application must allocate enough space */
) 
{
  int distance;         /* Input: which coloring to perform;
                           currently only supports D1 and D2 coloring */
  int ss;               /* Superstep size: detemines how many vertices are
                           locally colored before the next color
                           information exchange */
  char comm_pattern;    /* (A) asynchronous (S) synchronous supersteps*/
  char color_order;     /* (I) interior vertices first
                           (B) boundary vertices first (U) interleaved */
  char color_method;    /* Coloring method. (F) First fit
                           (S) staggered first fit (L) load balancing */

  static char *yo = "color_fn";
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt, *partvec;
  int *adjproc, nvtx=num_obj, gVtx;
  int *input_parts;                   /* Initial partitions for objects. */
  float *ewgts, *float_vwgt;
  int *color;
  int graph_type, check_graph;
  int i, j;
  int obj_wgt_dim, edge_wgt_dim;
  int ierr = ZOLTAN_OK;

      
  Zoltan_Bind_Param(Color_params, "DISTANCE", (void *) &distance);
  Zoltan_Bind_Param(Color_params, "SUPERSTEP_SIZE", (void *) &ss);
  Zoltan_Bind_Param(Color_params, "COMM_PATTERN", (void *) &comm_pattern);
  Zoltan_Bind_Param(Color_params, "COLOR_ORDER", (void *) &color_order);
  Zoltan_Bind_Param(Color_params, "COLORING_METHOD", (void *) &color_method);
  
  /* Set default values */
  distance = 1;
  ss = 100;
  color_method = 'F';
  color_order = 'I';
  comm_pattern = 'A';

  Zoltan_Assign_Param_Vals(zz->Params, Color_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  /* Compute Max number of array entries per ID over all processors.
     This is a sanity-maintaining step; we don't want different
     processors to have different values for these numbers. */
/*    comm[0] = zz->Num_GID;
    comm[1] = zz->Num_LID;
    MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
    zz->Num_GID = *num_gid_entries = gcomm[0];
    zz->Num_LID = *num_lid_entries = gcomm[1];
*/
   /* Return if this processor is not in the Zoltan structure's
      communicator. */
  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
      return ZOLTAN_OK;
  }
  /* Check that the user has allocated space for the return args. */
  if (!color_exp)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Output argument is NULL. Please allocate all required arrays before calling this routine.");

  /* Construct the heterogenous machine description. */
  /*  ierr = Zoltan_Build_Machine_Desc(zz);
      if (ierr == ZOLTAN_FATAL)
      ZOLTAN_COLOR_ERROR(ierr, "Error in constructing heterogeneous machine description.");
  */
  
  /* Initialize all local pointers to NULL. This is necessary
     because we free all non-NULL pointers upon errors. */
  vtxdist = xadj = adjncy = vwgt = adjwgt = adjproc = NULL;
  ewgts = float_vwgt = NULL;
  input_parts = NULL;

  /* Default graph type is GLOBAL. */
  graph_type = GLOBAL_GRAPH;
  check_graph = 1;
  obj_wgt_dim = 0; /* We do not use weights */
  edge_wgt_dim = 0;

  /* Get object ids and part information */
  ierr = Zoltan_Get_Obj_List(zz, &nvtx, &global_ids, &local_ids,
                             obj_wgt_dim, &float_vwgt, &input_parts);
  if (ierr) { /* Return error */      
      ZOLTAN_COLOR_ERROR(ierr, "Get_Obj_List returned error.");
  }

#if 0
  printf("after Get_Obj_list\n");
  for (i=0; i<nvtx; ++i)
      printf("%d: gid=%d lid=%d part=%d\n", i, global_ids[i], local_ids[i], input_parts[i]);
#endif
  
  /* Build ParMetis data structures, or just get vtxdist. */
  ierr = Zoltan_Build_Graph(zz, graph_type, check_graph, nvtx,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_COLOR_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }

  /* Vertex IDs start from 0 */
  for (i=0; i<nvtx; i++)
      global_ids[i] = global_ids[i]-1;
  /* Determine the global number of vertices */
  MPI_Allreduce(&nvtx, &gVtx, 1, MPI_INT, MPI_SUM, zz->Communicator);

  /* create part vector */
  partvec = (int *) ZOLTAN_MALLOC(gVtx*sizeof(int));
  if (!partvec)
      MEMORY_ERROR;
  for (i=0; i<xadj[nvtx]; i++)
      partvec[adjncy[i]] = adjproc[i];

#if 0
  printf("GRAPH: nvtx:%d Proc:%d\n", nvtx, zz->Proc);
  for (i=0; i < nvtx; i++) {
      printf("%d :: ", global_ids[i]);
      for (j=xadj[i]; j < xadj[i+1]; j++)
          printf("%d ", adjncy[j]);
      printf("\n");
  }
#endif
  
  /* Allocate color array. Local and D1 neighbor colors will be stored here. */ 
  if (!(color = (int *) ZOLTAN_CALLOC(gVtx, sizeof(int))))
      MEMORY_ERROR;

  if (distance == 1)  
      D1coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, gVtx, global_ids, vtxdist, &xadj, &adjncy, partvec, color);
  else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Coloring with requested distance is not implemented. Using Distance-1 coloring.");
      D1coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, gVtx, global_ids, vtxdist, &xadj, &adjncy, partvec, color);
  }

  /* fill the return array */
  for (i=0; i<nvtx; i++)
      color_exp[i] = color[global_ids[i]]; 
  
#if 0
  printf("P%d: Vtx colors: ", zz->Proc);
  for (i=0; i<nvtx; i++) {
      int gu = global_ids[i];
      printf("%d->%d  ", gu, color[gu]);
  }
  printf("\n");
#endif
  
  /* Check if there is an error in coloring */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) { /* DB: ZOLTAN_DEBUG_ALL)*/
      for (i=0; i<nvtx; i++) {
          int gu = global_ids[i];
          for (j = xadj[i]; j < xadj[i+1]; j++) {
              int gv = adjncy[j];
              if (color[gu] == color[gv])
                  printf("Error in coloring! u:%d, v:%d, cu:%d, cv:%d\n", gu, gv, color[gu], color[gv]);
          }
      }
  }  

 End:
  Zoltan_Multifree(__FILE__,__LINE__, 7, &vtxdist, &xadj, &adjncy, &adjproc,
                   &input_parts, &partvec, &color);
  
  return ierr;
}





int D1coloring(
    ZZ *zz,
    char color_order,
    char color_method,
    char comm_pattern,
    int ss,
    int nVtx,
    int gVtx,
    int *global_ids,
    int *vdist,
    int **pxadj,
    int **padj,
    int *partvec,
    int *color
)
{
    static char *yo = "D1coloring";
    int i;
    int nColor = 0;
    int nConflict;
    int nTotConflict;
    int nRound;
    int confCont;     
    int *rand_key = NULL;
    int *xadj = *pxadj;
    int *adj = *padj;
    int *xbadj = NULL;
    int nbound;
    int *isbound=NULL;
    int *visit = NULL;
    int *mark = NULL;
    int *conflicts = NULL;
    int *replies = NULL;
    MPI_Status *stats = NULL;
    MPI_Request *rreqs = NULL;
    MPI_Request *sreqs = NULL;
    int *rreqfrom = NULL;
    int **newcolored = NULL;
    double times[6];
    int get_times;
    int ierr;

    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int));
    if (!isbound || !visit)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    
    /* Start timer */
    get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
    if (get_times){
        MPI_Barrier(zz->Communicator);
        times[0] = Zoltan_Time(zz->Timer);
    }    
    ReorderGraph(zz, nVtx, &xadj, &xbadj, &adj, partvec, &nbound, isbound, visit);   
    if (get_times) times[1] = Zoltan_Time(zz->Timer);

#ifdef _DEBUG2    
    printf("[%d]: nbound:%d. Visit order: ", zz->Proc, nbound);
    for (i=0; i<nVtx; i++)
        printf("%d ", visit[i]);
    printf("\n");
#endif
    
    /* Memory allocation */
    mark = (int *) ZOLTAN_CALLOC(gVtx+1, sizeof(int));
    conflicts = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int));
    replies = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    stats = (MPI_Status *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Status));
    rreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    rreqfrom = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    newcolored = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    if (!mark || !conflicts || !replies || !stats || !rreqs || !sreqs || !rreqfrom || !newcolored)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    for (i=0; i<zz->Num_Proc; ++i) {
        newcolored[i] = (int *) ZOLTAN_MALLOC(2 * ss * sizeof(int));
        if (!newcolored[i])
            ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }

    /* Generate random numbers associated with vertices, with fixed seed */
    srand48(0); 
    rand_key= (int *) ZOLTAN_MALLOC(sizeof(int) * gVtx);
    if (!rand_key)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    for(i=0; i<gVtx; i++)
        rand_key[i] = drand48()*1000000;

    /* Color internal vertices and determine the visit order */
    if (get_times) {
        MPI_Barrier(zz->Communicator);
        times[2] = Zoltan_Time(zz->Timer);
    }
    if (color_order == 'U') { 
        nConflict = nVtx;
        for (i=0; i<nVtx; i++)
            visit[i] = i;
        if (zz->Num_Proc==1)
            IntColoring(zz, &nColor, nVtx, gVtx, visit, xadj, adj, global_ids, color, mark, color_method);
    }
    else if (color_order == 'I') {
        IntColoring(zz, &nColor, nVtx - nbound, gVtx, visit + nbound, xadj, adj, global_ids, color, mark, color_method);
        nConflict = nbound;
    }
    else if (color_order=='B')
        nConflict = nbound;

    if (get_times) times[3] = Zoltan_Time(zz->Timer);
        
    nRound = 0; 
    nTotConflict = 0;
            
    /* Color boundary vertices */
    if (zz->Num_Proc >= 2) {
        do {
            int *tp = visit;
            memset(mark, 0xff, (gVtx+1) * sizeof(int));
            ParallelColoring(zz, nConflict, global_ids, visit, xadj, adj, isbound, ss,
                             &nColor, color, newcolored, mark, color_method, comm_pattern,
                             rreqfrom, replies, sreqs, rreqs, stats);
            nConflict = Conflict(zz, nConflict, global_ids, visit, xadj, xbadj, adj,
                                 color, conflicts, rand_key);
            /* swap conflicts list with visit list so that if there are conflicts,
               next coloring will color them */
            visit = conflicts;
            conflicts = tp;
            confCont = 0;
            MPI_Allreduce(&nConflict, &confCont, 1, MPI_INT, MPI_SUM, zz->Communicator);
            nTotConflict += confCont;
            ++nRound;
        } while (confCont);
    }

    if (get_times) times[4] = Zoltan_Time(zz->Timer);
    /* Color internal vertices after boundaries if boundary first ordering */
    if (color_order == 'B')
        IntColoring(zz, &nColor, nVtx - nbound, gVtx, visit + nbound, xadj, adj, global_ids, color, mark, color_method);
    if (get_times) {
        MPI_Barrier(zz->Communicator);
        times[5] = Zoltan_Time(zz->Timer);
    }
    	    
    /* Output timing results if desired */
    if (get_times){
        if (zz->Proc == zz->Debug_Proc) printf("\nZOLTAN timing statistics:\n");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, nColor, 
                           " Number of Colors                           ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, nTotConflict, 
                           " Number of Conflicts                        ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, nRound, 
                           " Number of Rounds                           ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[1]-times[0], 
                           " Graph reordering time                      ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[3]-times[2]+times[5]-times[4], 
                           " Internal vertex coloring time              ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[4]-times[3], 
                           " Parallel boundary vertex coloring time     ");
        Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[5]-times[0], 
                           " Total coloring time (including reordering) ");
        if (zz->Proc==zz->Debug_Proc) printf("\n");
    }
    
    *pxadj = xadj;
    *padj = adj;   

    ierr = ZOLTAN_OK;
    
 End:

    ZOLTAN_FREE(&mark);
    ZOLTAN_FREE(&conflicts);
    ZOLTAN_FREE(&replies);
    ZOLTAN_FREE(&stats);
    ZOLTAN_FREE(&rreqs);
    ZOLTAN_FREE(&sreqs);
    ZOLTAN_FREE(&rreqfrom);
    for (i=0; i<zz->Num_Proc; i++)
        ZOLTAN_FREE(&newcolored[i]);
    ZOLTAN_FREE(&newcolored);
    ZOLTAN_FREE(&isbound);
    ZOLTAN_FREE(&visit);
    ZOLTAN_FREE(&xbadj);
    ZOLTAN_FREE(&rand_key);
    
    return ierr;
}


/* Reorder graph */
int ReorderGraph(
    ZZ *zz,           
    int nVtx,         /* In: Number of objects */
    int **pxadj,      /* In/Out: Pointer to xadj */
    int **pxbadj,     /* Out: Pointer to start of cut edges in the
                         adj lists of boundary vertices*/
    int **padj,       /* In/Out: Pointer to adj */
    int *partvec,     /* In: Part vector */
    int *nbound,      /* Out: Number of boundary vertices */
    int *isbound,     /* Out: Indicates if a vertex is on the boundary */
    int *visit        /* Out: Visit order */
)
{
    static char *yo = "ReorderGraph";
    int *oxadj = *pxadj; /* xadj before reordering            */
    int *oadj = *padj;   /* adj before reordering             */
    int *xadj = NULL;      /* xadj after reordering             */
    int *xbadj = NULL;     /* pointer to start of cut edges in the
                            adj lists of boundary vertices    */
    int *adj = NULL;       /* adj after reordering              */
    int i, j, idx, nboundIdx;
    int ierr;
    
    /* Memory allocation */
    xadj = (int *) ZOLTAN_CALLOC(nVtx+1, sizeof(int));
    xbadj = (int *) ZOLTAN_CALLOC(nVtx+1, sizeof(int));
    adj = (int *) ZOLTAN_CALLOC(oxadj[nVtx], sizeof(int));

    if (!xadj || !xbadj || !adj)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");

    /* Determine the boundary vertices and fill isbound */
    *nbound = 0;
    for (i=0; i<nVtx; ++i) {
        isbound[i] = 0;        
        for (j=oxadj[i]; j<oxadj[i+1]; ++j) {
            int v = oadj[j];
            if (partvec[v] != zz->Proc) {
                isbound[i] = 1;
                ++(*nbound);
                break;
            }
        }
    }

    /* Order vertices: boundaries first, internals last */
    idx = 0;
    nboundIdx = *nbound;
    for (i=0; i<nVtx; ++i) {
        if (isbound[i])
            visit[idx++] = i;
        else
            visit[nboundIdx++] = i;
    }

    /* xadj and adj after reordering */
    /* Currently there is no reordering so just copy the original lists */
    for (i=0; i<nVtx; i++) {
        xadj[i] = oxadj[i];
        for (j=oxadj[i]; j<oxadj[i+1]; j++)
            adj[j] = oadj[j];
    }
    xadj[nVtx] = oxadj[nVtx];
    
    /* move cut edges to the beginning of adj lists for boundary vertices */
    for (i=0; i<nVtx; ++i) {
        int j, b, tmp;
        b = xadj[i+1] - 1;
        j = xadj[i];
        while (b > j) {
            for ( ; j<b && zz->Proc != partvec[adj[j]]; ++j);
            for ( ; j<b && zz->Proc == partvec[adj[b]]; --b);
            if (j < b) {
                SWAP(adj[j], adj[b]);
                ++j; --b;
            }
        }
        xbadj[i] = j;
        if ((b==j) && (zz->Proc != partvec[adj[j]]))
            ++xbadj[i];
    }

    *pxadj = xadj;
    *pxbadj = xbadj;
    *padj = adj;

    ierr = ZOLTAN_OK;
    
 End:

    ZOLTAN_FREE(&oxadj);
    ZOLTAN_FREE(&oadj);
    
    return ierr;
}


/* Pick a color for a vertex based on the coloring method */
int pickColor(ZZ *zz, char color_method, int u, int oldColor, int *nColor, int *mark)
{
    int c, pickedColor;
    static char *yo="pickColor";
    
    switch(color_method){    
    case 'F':
    default: 
        if (color_method !='F')
            ZOLTAN_PRINT_WARN(zz->Proc, yo, "Unknown coloring method. Using First Fit (F).");	      
        for (c = 1; (c <= *nColor) && (mark[c] == u); ++c) ;
        if (c > *nColor)  /* no available color with # less than nColor */
            c = ++(*nColor);
        pickedColor = c;
    }
    
    return pickedColor;   
}

/* Color internal vertices */
int IntColoring(ZZ *zz, int *nColor, int nVtx, int gVtx, int *visit, int * xadj, int *adj, int *global_ids, int *color, int *mark, char color_method)
{
    int i, j;

    memset(mark, 0xff, (gVtx+1) * sizeof(int));

    for (i=0; i<nVtx; ++i) {
        int u = visit[i], gu;        
        for (j = xadj[u]; j < xadj[u+1]; ++j) {
            int gv = adj[j], c;            
            if ((c = color[gv]) != 0) 
                mark[c] = u;
        }
        gu = global_ids[u];
	color[gu] = pickColor(zz, color_method, u, color[gu], nColor, mark);
    }
    return ZOLTAN_OK;
}


/* Parallel coloring of boundary vertices */
int ParallelColoring(ZZ *zz, int nVtx, int *global_ids, int *visit, int *xadj, int *adj, int *isbound, int ss, int *nColor, int *color, int **newcolored, int *mark, char color_method, char comm_pattern, int *rreqfrom, int *replies, MPI_Request *sreqs, MPI_Request *rreqs, MPI_Status *stats)
{
    static char *yo="ParallelColoring";
    int colortag=1001, i, j, p, q, l;
    int *colored, n=0;
    int rreqcnt=0, sreqcnt=0, repcount;
    int ierr;
    
    colored = newcolored[zz->Proc];
    
    /*** Issue async recvs *******************************/
    for (rreqcnt = p = 0; p < zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            rreqfrom[rreqcnt] = p;
            if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[rreqcnt]))
                ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
            ++rreqcnt;
        }
    
    /*** Coloring *************************************/
    for (i=0; i<nVtx; ++i) {
        int u = visit[i], gu;        
        for (j=xadj[u]; j<xadj[u+1]; ++j) {
            int gv = adj[j], c;            
            if ((c = color[gv]) != 0) 
                mark[c] = u;
        }
        gu = global_ids[u];
        color[gu] = pickColor(zz, color_method, u, color[gu], nColor, mark);

        if (!isbound[u]) /* send only boundary vertices */
            continue;
        colored[n++] = gu;
        colored[n++] = color[gu];
        
	/*** If superstep is finished, communicate ******************************/
        if (n >= 2*ss) {
            for (sreqcnt = p = 0; p < zz->Num_Proc; ++p)
                if (p != zz->Proc) {
                    MPI_Isend(colored, 2*ss, MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
                    ++sreqcnt;
                }
            
            /* wait color lists from other processors */
            if (comm_pattern == 'S') { /* wait all results if sync communication*/
                MPI_Waitall(rreqcnt, rreqs, stats); 
                repcount = rreqcnt;
            }
            else /* wait some results if async communication*/
                MPI_Waitsome(rreqcnt, rreqs, &repcount, replies, stats); 

            for (l = repcount-1; l >= 0; --l) {
                if (comm_pattern == 'S') 
                    q = l;                 
                else
                    q = replies[l];                
                p = rreqfrom[q];

                /* Read received color list from p */
                for (j = 0; j < 2*ss; ) {
                    int gv  = newcolored[p][j++], c;                    
                    if (gv < 0) 
                        break;
                    c = newcolored[p][j++];
		    color[gv] = c;                    
                }
                /* If p hasn't finished coloring, issue new color request */
                if (j >= 2*ss) {
                    if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[q]))
                        ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");                    
                }
                else { 
                    rreqs[q] = rreqs[--rreqcnt];     /* we don't need rreqs[q]; delete it by */
                    rreqfrom[q] = rreqfrom[rreqcnt]; /* overwriting with the last one */
                } 
            }  
            n = 0;
            MPI_Waitall(sreqcnt, sreqs, stats); /* wait all to be sent */
        }
    }

    /* Mark that coloring is finished on zz->Proc and send the last color list*/
    colored[n++] = -1;
    colored[n++] = -1;
    for (sreqcnt = p = 0; p < zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            MPI_Isend(colored, 2*ss, MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
            ++sreqcnt;
        }
    
    /* Receive the remaining color lists from other processors */
    while (rreqcnt) {
        if (comm_pattern == 'S') { /* wait all results if sync communication*/
            MPI_Waitall(rreqcnt, rreqs, stats); /* wait all results */
            repcount = rreqcnt;
        }
        else
            MPI_Waitsome(rreqcnt, rreqs, &repcount, replies, stats); /* wait some results */

        for (l=repcount-1; l>=0; --l) {
            if (comm_pattern == 'S') /* wait all results if sync communication*/
                q = l;                
            else
                q = replies[l];
            p = rreqfrom[q];
            
            /* Read received color list from p */
            for (j = 0; j < 2*ss; ) {
                int gv  = newcolored[p][j++], c;                
                if (gv < 0) 
                    break;
                c = newcolored[p][j++];
		color[gv] = c;
            }
            
            /* If p hasn't finished coloring, issue new color request */
            if (j >= 2*ss) {
                if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[q]))
                    ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
            }
            else { 
                rreqs[q] = rreqs[--rreqcnt];     /* we don't need rreqs[q]; delete it by */
                rreqfrom[q] = rreqfrom[rreqcnt]; /* overwriting with the last one */
            } 
        }  
    }
    MPI_Waitall(sreqcnt, sreqs, stats); /* wait all sends to complete */

    ierr = ZOLTAN_OK;
    
 End:
    
    return ierr;
}


/* Detect Conflicts */ 
int Conflict(ZZ *zz, int nConflict, int *global_ids, int *visit, int *xadj, int *xbadj, int *adj, int *color, int *conflicts, int *rand_key)
{
    int i, j, conflict = 0;
    
    for (i = 0; i < nConflict; ++i) {
        int u = visit[i], gu;
        gu = global_ids[u];

        for (j = xadj[u]; j < xbadj[u]; ++j) {
            int gv = adj[j];

            if ((color[gu] == color[gv]) && (rand_key[gu] < rand_key[gv])){
		conflicts[conflict++] = u;
		break;
            }
        }
    }
    return conflict;
}


#ifdef __cplusplus
}
#endif
