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
#include "g2l_hash.h"
#include "params_const.h"
#include "zz_util_const.h"
#include "parmetis_jostle.h"
#include "all_allo_const.h"

    
/* Function prototypes */
static int ReorderGraph(ZZ *, int, int *, int **, int *, int *, int *, int *, int *);
static int PickColor(ZZ *, char, int, int, int *, int *);
static int InternalColoring(ZZ *zz, int *nColor, int nVtx, int *visit, int * xadj, int *adj, int *color, int *mark, int mark_size,char color_method);
static int ParallelColoring (ZZ *zz, int nvtx, int *visit, int *xadj, int *adj, int *isbound, int ss, int *nColor, int *color, int **newcolored, int *mark, int gmaxdeg, G2LHash *hash, char color_method, char comm_pattern, int *rreqfrom, int *replies, MPI_Request *sreqs, MPI_Request *rreqs, MPI_Status *stats);
static int DetectConflicts(ZZ *, int, int *, int *, int *, int *, int *, int *, int *);
static int D1coloring(ZZ *zz, char color_order, char color_method, char comm_pattern, int ss, int nVtx, G2LHash *hash, int *xadj, int *adj, int *adjproc, int *color);
static int GenPrime(int stop, int *prime_num);
    
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
    
/*****************************************************************************/
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
    
/*****************************************************************************/
/* Interface routine for Graph Coloring */

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
  char color_method;    /* Coloring method. (F) First fit */


  static char *yo = "color_fn";
  idxtype *vtxdist, *xadj, *adjncy; /* arrays to store the graph structure */
  int *adjproc;                     
  int *input_parts;                 /* Initial partitions for objects. */
  int nvtx = num_obj;               /* number of vertices */
  idxtype *vwgt, *adjwgt;           /* weights - not used */
  float *ewgts, *float_vwgt;        /* weights - not used */
  int obj_wgt_dim, edge_wgt_dim;    /* weight dimensions - not used */
  int *color;                       /* array to store colors of local and D1
                                       neighbor vertices */
  int graph_type, check_graph;
  int i, j;
  int lastlno;                      /* total number of local and D1 neighbor vertices */
  G2LHash hash;                     /* hash to map global ids of local and D1 neighbor
                                       vertices to consecutive local ids */
  int hsize;                        /* hash size */
  int ierr = ZOLTAN_OK;


  /* PARAMETER SETTINGS */
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
  /* comm[0] = zz->Num_GID;
    comm[1] = zz->Num_LID;
    MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
    zz->Num_GID = *num_gid_entries = gcomm[0];
    zz->Num_LID = *num_lid_entries = gcomm[1];
  */
  /* Return if this processor is not in the Zoltan structure's
     communicator. */
  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz))
      return ZOLTAN_OK;

  /* Construct the heterogenous machine description. */
  /*  ierr = Zoltan_Build_Machine_Desc(zz);
      if (ierr == ZOLTAN_FATAL)
      ZOLTAN_COLOR_ERROR(ierr, "Error in constructing heterogeneous machine description.");
  */
  

  /* BUILD THE GRAPH */
  /* Check that the user has allocated space for the return args. */
  if (!color_exp)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Output argument is NULL. Please allocate all required arrays before calling this routine.");
  
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

  /* Build ParMetis data structures, or just get vtxdist. */
  ierr = Zoltan_Build_Graph(zz, graph_type, check_graph, nvtx,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_COLOR_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }

  /* CREATE THE HASH TABLE */
  /* Allocate hash table */
  /* UVCUVC: TODO we can allocate smaller hash; check this later */
  if (GenPrime(2*xadj[nvtx], &hsize)==ZOLTAN_MEMERR)
      MEMORY_ERROR;
  if (Zoltan_G2LHash_Create(&hash, hsize)==ZOLTAN_MEMERR)
      MEMORY_ERROR;

  /* Insert the ids of the local vertices into the hash table */
  /* By inserting local vertices first, it is ensured that the
     local ids of local vertices start from 0 and are consecutive*/
  for (j=0; j<nvtx; j++) 
      Zoltan_G2LHash_G2L(&hash, vtxdist[zz->Proc] + j);


  /* Add ids of the d1 neighbors into the hash table*/
  for (i=0; i<xadj[nvtx]; ++i)
      adjncy[i] = Zoltan_G2LHash_G2L(&hash, adjncy[i]);
  /* lastno is the total number of local and d1 neighbors */
  lastlno = hash.lastlno; 

#if 0
  printf("[%d] GRAPH after hashing: nvtx:%d\n", zz->Proc, nvtx);
  for (i=0; i < nvtx; i++) {
      printf("%d :: ", i);
      for (j=xadj[i]; j < xadj[i+1]; j++)
          printf("%d ", adjncy[j]);
      printf("\n");
  }
#endif
  
  /* Allocate color array. Local and D1 neighbor colors will be stored here. */ 
  if (!(color = (int *) ZOLTAN_CALLOC(lastlno, sizeof(int))))
      MEMORY_ERROR;

  /* SELECT COLORING ALGORITHM AND PERFORM COLORING*/
  if (distance == 1)
      D1coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, &hash, xadj, adjncy, adjproc, color);
  else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Coloring with requested distance is not implemented. Using Distance-1 coloring.");
      D1coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, &hash, xadj, adjncy, adjproc, color);
  }
  
  /* FILL THE RETURN ARRAY */
  for (i=0; i<nvtx; i++) 
      color_exp[i] = color[i];
  
  /* Check if there is an error in coloring */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) { 
      for (i=0; i<nvtx; i++) {
          for (j = xadj[i]; j < xadj[i+1]; j++) {
              int v = adjncy[j];
              if (color[i] == color[v])
                  printf("Error in coloring! u:%d, v:%d, cu:%d, cv:%d\n", i, v, color[i], color[v]);
          }
      }
  }

 End:
  
  Zoltan_Multifree(__FILE__,__LINE__, 6, &vtxdist, &xadj, &adjncy, &input_parts, &adjproc, &color);
  Zoltan_G2LHash_Destroy(&hash);
  
  return ierr;
}

/*****************************************************************************/
/* Distance-1 coloring */

static int D1coloring(
    ZZ *zz,
    char color_order,  /* (I) interior vertices first
                          (B) boundary vertices first (U) interleaved */
    char color_method, /* Coloring method. (F) First fit
                          (S) staggered first fit (L) load balancing */
    char comm_pattern, /* (A) asynchronous (S) synchronous supersteps */
    int ss,            /* Superstep size: detemines how many vertices are
                          locally colored before the next color
                          information exchange */
    int nvtx,          /* number of vertices in the graph */
    G2LHash *hash,     /* hash to map global ids of local and D1 neighbor
                          vertices to consecutive local ids */
    int *xadj,         /* arrays that store the graph structure */
    int *adj,
    int *adjproc,
    int *color         /* return array to store colors of local and D1
                          neighbor vertices */
)
{
    static char *yo = "D1coloring";
    int i;
    int nColor = 0;              /* Number of colors */
    int nConflict;               /* Number of local vertices to be recolored
                                    in the next round */
    int confCont;                /* Total number of vertices to be recolored
                                    in the next round */
    int nTotConflict;            /* Total number of conflicts over all rounds
                                    on all processors */ 
    int nRound;                  /* Number of parallel coloring rounds */
    int *rand_key = NULL;        /* Array of random numbers associated with
                                    global numbers of vertices */
    int nbound;                  /* Number of local boundary vertices */
    int *xbadj = NULL;           /* Pointer to end of cut edges in adj lists
                                    of local vertices */
    int *isbound=NULL;           /* Indicates whether a local vertex is a boundary
                                    vertex */
    int *visit = NULL;           /* Visit (coloring) order of local vertices */
    int *mark = NULL;            /* Array to mark forbidden colors for a vertex */
    int gmaxdeg = 0;             /* Maximum vertex degree in the graph */
    int lmaxdeg = 0;             /* Maximum vertex degree for the local vertices */
    int lastlno = hash->lastlno; /* Total number of local and D1 neighbor vertices */
    int *conflicts = NULL;       /* List of vertices to be recolored in the next
                                    round */
    int *replies = NULL;         /* Arrays used for MPI communication */
    MPI_Status *stats = NULL;
    MPI_Request *rreqs = NULL;
    MPI_Request *sreqs = NULL;
    int *rreqfrom = NULL;
    int **newcolored = NULL;     /* Array used for communicating boundary vertex
                                    colors at the end of supersteps. Global number
                                    of the vertex - color of vertex pairs are filled
                                    in this array */
    double times[6];             /* Used for timing measurements */
    int get_times;               /* (1) Measure timings (0) Don't */ 
    int ierr;

    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    if (!isbound || !visit)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    
    /* Start timer */
    get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
    if (get_times){
        MPI_Barrier(zz->Communicator);
        times[0] = Zoltan_Time(zz->Timer);
    }    
    ierr = ReorderGraph(zz, nvtx, xadj, &xbadj, adj, adjproc, &nbound, isbound, visit);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
        ZOLTAN_COLOR_ERROR(ierr, "Error in ReorderGraph");
    if (get_times) times[1] = Zoltan_Time(zz->Timer);

#if 0
    printf("After reordering: nvtx:%d Proc:%d\n", nvtx, zz->Proc);
    for (i=0; i < nvtx; i++) {
        int j;
        printf("%d :: ", i);
        for (j=xadj[i]; j < xadj[i+1]; j++)
            printf("%d ", adj[j]);
        printf("\n");
    }
#endif

    /* Calculate the maximum degree of the graph */
    lmaxdeg = 0;
    for (i=0; i<nvtx; i++)
        if (lmaxdeg < xadj[i+1] - xadj[i])
            lmaxdeg = xadj[i+1] - xadj[i];    
    MPI_Allreduce(&lmaxdeg, &gmaxdeg, 1, MPI_INT, MPI_MAX, zz->Communicator);        
    /* gmaxdeg+1 is the upper bound for #colors */
    ++gmaxdeg; 
    
    /* Memory allocation */
    mark = (int *) ZOLTAN_CALLOC(gmaxdeg, sizeof(int)); 
    conflicts = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
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

    /* Generate random numbers associated with global numbers of the vertices */
    /* All processors generate the same random number corresponding
       to the same global vertex number */
    rand_key = (int *) ZOLTAN_MALLOC(sizeof(int) * lastlno);
    if (!rand_key)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    for(i=0; i<lastlno; i++) {
        srand48(Zoltan_G2LHash_L2G(hash, i));
        rand_key[i] = drand48()*1000000;
    }

    /* Color internal vertices and determine the visit order */
    if (get_times) {
        MPI_Barrier(zz->Communicator);
        times[2] = Zoltan_Time(zz->Timer);
    }
    if (color_order == 'U') { 
        nConflict = nvtx;
        for (i=0; i<nvtx; i++)
            visit[i] = i;
        if (zz->Num_Proc==1)
            InternalColoring(zz, &nColor, nvtx, visit, xadj, adj, color, mark, gmaxdeg, color_method);
    }
    else if (color_order == 'I') {
        InternalColoring(zz, &nColor, nvtx - nbound, visit + nbound, xadj, adj, color, mark, gmaxdeg, color_method);
        nConflict = nbound;
    }
    else if (color_order=='B')
        nConflict = nbound;

    if (get_times) times[3] = Zoltan_Time(zz->Timer);
        
    /* Color boundary vertices */
    nRound = 0; 
    nTotConflict = 0;
    if (zz->Num_Proc >= 2) {
        do {
            int *tp = visit;
            memset(mark, 0xff, gmaxdeg * sizeof(int));
            ierr = ParallelColoring(zz, nConflict, visit, xadj, adj, isbound, ss,
                                    &nColor, color, newcolored, mark, gmaxdeg, hash,
                                    color_method, comm_pattern, rreqfrom, replies,
                                    sreqs, rreqs, stats);
            if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                ZOLTAN_COLOR_ERROR(ierr, "Error in ParallelColoring");
            nConflict = DetectConflicts(zz, nConflict, visit, xadj, xbadj, adj,
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
        InternalColoring(zz, &nColor, nvtx-nbound, visit + nbound, xadj, adj, color, mark, gmaxdeg, color_method);

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


/*****************************************************************************/
/* Reorder graph adjacency list */

static int ReorderGraph(
    ZZ *zz,           
    int nvtx,         /* In: Number of objects */
    int *xadj,        /* In: Pointer to xadj */
    int **pxbadj,     /* Out: Pointer to start of cut edges in the
                         adj lists of boundary vertices*/
    int *adj,         /* In/Out: Pointer to adj */
    int *adjproc,     /* In: Part vector */
    int *nbound,      /* Out: Number of boundary vertices */
    int *isbound,     /* Out: Indicates if a vertex is on the boundary */
    int *visit        /* Out: Visit order */
)
{
    static char *yo = "ReorderGraph";
    int *xbadj = NULL;     /* pointer to start of cut edges in the
                              adj lists of boundary vertices    */  
    int i, j, idx, nboundIdx;
    int ierr;
    
    /* Memory allocation */
    xbadj = (int *) ZOLTAN_CALLOC(nvtx+1, sizeof(int));
    if (!xbadj)
        ZOLTAN_COLOR_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    
    /* Determine the boundary vertices and fill isbound */
    *nbound = 0;
    for (i=0; i<nvtx; ++i) {
        isbound[i] = 0;        
        for (j=xadj[i]; j<xadj[i+1]; ++j) {
            if (adjproc[j] != zz->Proc) {
                isbound[i] = 1;
                ++(*nbound);
                break;
            }
        }
    }

    /* Order vertices: boundaries first, internals last */
    idx = 0;
    nboundIdx = *nbound;
    for (i=0; i<nvtx; ++i) {
        if (isbound[i])
            visit[idx++] = i;
        else
            visit[nboundIdx++] = i;
    }

    /* move cut edges to the beginning of adj lists for boundary vertices */
    for (i=0; i<nvtx; ++i) {
        int j, b, tmp;
        b = xadj[i+1] - 1;
        j = xadj[i];
        while (b > j) {
            for ( ; j<b && zz->Proc != adjproc[j]; ++j);
            for ( ; j<b && zz->Proc == adjproc[b]; --b);
            if (j < b) {
                SWAP(adj[j], adj[b]);
                SWAP(adjproc[j], adjproc[b]);
                ++j; --b;
            }
        }
        xbadj[i] = j;
        if ((b==j) && (zz->Proc != adjproc[j]))
            ++xbadj[i];
    }

    *pxbadj = xbadj;

    ierr = ZOLTAN_OK;
    
 End:

    return ierr;
}


/*****************************************************************************/
/* Pick a color for a vertex based on the coloring method */

static int PickColor(
    ZZ *zz,
    char color_method,
    int u,
    int oldColor,
    int *nColor,
    int *mark
)
{
    int c, pickedColor;
    static char *yo="PickColor";
    
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

/*****************************************************************************/
/* Color internal vertices */

static int InternalColoring(
    ZZ *zz,
    int *nColor,
    int nvtx,
    int *visit,
    int * xadj,
    int *adj,
    int *color,
    int *mark,
    int gmaxdeg,
    char color_method
)
{
    int i, j;

    memset(mark, 0xff, gmaxdeg * sizeof(int));

    for (i=0; i<nvtx; ++i) {
        int u = visit[i];        
        for (j = xadj[u]; j < xadj[u+1]; ++j) {
            int v = adj[j], c;            
            if ((c = color[v]) != 0) 
                mark[c] = u;
        }
	color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);
    }
    return ZOLTAN_OK;
}

/*****************************************************************************/
/* Parallel coloring of boundary vertices */       

static int ParallelColoring (
    ZZ *zz,
    int nvtx,         /* number of vertices in the graph */
    int *visit,       /* Visit (coloring) order of local vertices */
    int *xadj,        /* arrays that store the graph structure */
    int *adj,
    int *isbound,     /* Indicates whether a local vertex is a boundary
                         vertex */
    int ss,           /* Superstep size: detemines how many vertices are
                         locally colored before the next color
                         information exchange */
    int *nColor,      /* Number of colors */
    int *color,       /* return array to store colors of local and D1
                         neighbor vertices */
    int **newcolored, /* Array used for communicating boundary vertex
                         colors at the end of supersteps. Global number
                         of the vertex - color of vertex pairs are filled
                         in this array */
    int *mark,        /* Array to mark forbidden colors for a vertex */
    int gmaxdeg,      /* Maximum vertex degree in the graph */
    G2LHash *hash,    /* hash to map global ids of local and D1 neighbor
                         vertices to consecutive local ids */ 
    char color_method, /* Coloring method. (F) First fit
                          (S) staggered first fit (L) load balancing */
    char comm_pattern, /* (A) asynchronous (S) synchronous supersteps */
    int *rreqfrom,     /* Arrays used for MPI communication */
    int *replies,
    MPI_Request *sreqs,
    MPI_Request *rreqs,
    MPI_Status *stats
)
{
    static char *yo="ParallelColoring";
    int colortag=1001, i, j, p, q, l;
    int *colored, n=0;
    int rreqcnt=0, sreqcnt=0, repcount;
    int ierr;
    
    colored = newcolored[zz->Proc];
    
    /* Issue async recvs */
    for (rreqcnt = p = 0; p < zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            rreqfrom[rreqcnt] = p;
            if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[rreqcnt]))
                ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
            ++rreqcnt;
        }
    
    /* Coloring */
    for (i=0; i<nvtx; ++i) {
        int u = visit[i];        
        for (j=xadj[u]; j<xadj[u+1]; ++j) {
            int gv = adj[j], c;            
            if ((c = color[gv]) != 0) 
                mark[c] = u;
        }
        color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);

        if (!isbound[u]) /* send only boundary vertices */
            continue;
        colored[n++] = Zoltan_G2LHash_L2G(hash, u);
        colored[n++] = color[u];
        
	/* If superstep is finished, communicate */
        if (n >= 2*ss) {
            for (sreqcnt = p = 0; p < zz->Num_Proc; ++p)
                if (p != zz->Proc) {
                    MPI_Isend(colored, 2*ss, MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
                    ++sreqcnt;
                }
            
            /* Wait color lists from other processors */
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
                    int v  = newcolored[p][j++], c, hv;                    
                    if (v < 0) 
                        break;
                    hv = Zoltan_G2LHash_G2L(hash, v);
                    c = newcolored[p][j++];
		    color[hv] = c;                    
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
                int v  = newcolored[p][j++], c, hv;                
                if (v < 0) 
                    break;
                hv = Zoltan_G2LHash_G2L(hash, v);
                c = newcolored[p][j++];
		color[hv] = c;
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

/*****************************************************************************/
/* Detect conflicts */

static int DetectConflicts(
    ZZ *zz,
    int nConflict,
    int *visit,
    int *xadj,
    int *xbadj,
    int *adj,
    int *color,
    int *conflicts,
    int *rand_key
)
{
    int i, j, conflict = 0;
    
    for (i = 0; i < nConflict; ++i) {
        int u = visit[i];
        for (j = xadj[u]; j < xbadj[u]; ++j) {
            int v = adj[j];
            if ((color[u] == color[v]) && (rand_key[u] < rand_key[v])){
		conflicts[conflict++] = u;
		break;
            }
        }
    }
    return conflict;
}

/*****************************************************************************/
/* Returns the prime number closest to (and smaller than) stop */

static int GenPrime(
    int stop,
    int *prime_num
)
{
    int nap;
    int num, c, i;
    int *prime;
    
    prime = (int *) ZOLTAN_MALLOC((stop/2) * sizeof(int));
    if (!prime)
        return ZOLTAN_MEMERR;
    
    prime[0] = 2; prime[1] = 3;
    c = 2; /* initial primes */
    
    /* only have to check odd numbers */
    for (num=5; num < stop; num = num + 2) {
        nap = 0;  /* set not-a-prime false */        
        /* cycle through list of known primes */
        for (i=0; i < c; i++) { 
            /* check if a previous prime divides evenly */
            /* if so the number is not a prime */
            if ((num % prime[i]) == 0) {
                nap = 1;
                break;
            }            
            /* stop if prime squared is bigger than the number */
            if ((prime[i] * prime[i]) > num)
                break;
        }        
        /* if not-a-prime, then we found a prime */
        if (nap != 1) {
           /* add prime to list of known primes */
            prime[c] = num;
            c++;
        }
    }
    *prime_num = prime[c-1];

    ZOLTAN_FREE(&prime);
    return ZOLTAN_OK;
}




