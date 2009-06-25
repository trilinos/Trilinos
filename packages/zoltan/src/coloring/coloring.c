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
#include "third_library_const.h"
#include "all_allo_const.h"
#include "zz_rand.h"

/* when sending new colored vertices to processors,
   sent only the "relevant" ones; i.e., send the color info if the processors has
   D1 neighbour */
#define RELEVANT_COLORS
    
#define COLORTAG     1001
#define SNTAG        1002
#define XFORBIDTAG   1003
#define FORBIDTAG    1004
#define RECOLORTAG   1005
    
/* Function prototypes */
static int D1coloring(ZZ *zz, char color_order, char color_method, char comm_pattern, int ss,
                      int nVtx, G2LHash *hash, int *xadj, int *adj, int *adjproc, int *color);
static int D2coloring(ZZ *zz, char color_order, char color_method, char comm_pattern, int ss,
                      int nVtx, G2LHash *hash, int *xadj, int *adj, int *adjproc, int *color);

static int ReorderGraph(ZZ *, int, int, int *, int **, int *,
                        int *, int *, int *, int *);
static int PickColor(ZZ *, char, int, int, int *, int *);
static int InternalColoring(ZZ *zz, int distance, int *nColor,
                            int nVtx, int *visit, int * xadj, int *adj,
                            int *color, int *mark, int mark_size, char color_method);

#ifdef RELEVANT_COLORS    
static int D1ParallelColoring (ZZ *zz, int nvtx, int *visit, int *xadj, int *adj,
                               int *isbound, int ss, int *nColor, int *color,
                               int **newcolored, int *mark, int gmaxdeg, G2LHash *hash,
                               char color_method, char comm_pattern, int *rreqfrom,
                               int *replies, MPI_Request *sreqs, MPI_Request *rreqs,
                               MPI_Status *stats,
                               int *xadjproc, int *adjproc, int **persSbuf, int *Ssize, int plstcnt, int *plst);
#else
static int D1ParallelColoring (ZZ *zz, int nvtx, int *visit, int *xadj, int *adj,
                               int *isbound, int ss, int *nColor, int *color,
                               int **newcolored, int *mark, int gmaxdeg, G2LHash *hash,
                               char color_method, char comm_pattern, int *rreqfrom,
                               int *replies, MPI_Request *sreqs, MPI_Request *rreqs,
                               MPI_Status *stats);
#endif
    
static int sendNextStepsForbiddenColors(ZZ *zz, G2LHash *hash, int nlvtx, int p, int **srp, int *xadj, int *adj, int *xadjnl, int *adjnl, int *forbsizeS, int **xforbiddenS, int **forbiddenS, int *nColor, int *color, int *mark, int *confChk, int sncnt, int *wset, int *wsize, MPI_Request *sreqsFx, MPI_Request *sreqsF);
static int waitPtrAndForbiddenColors(ZZ* zz, int rreqcntFx, MPI_Request *rreqsFx, int *rreqfromFx, MPI_Status *stats, int *forbsize, int **xforbidden, int **forbidden, int **xfp, MPI_Request *rreqsF);
static int D2ParallelColoring (ZZ *zz, int nvtx, int nlvtx, int *visit, int *xadj, int *adj, int *xbadj, int *xadjnl, int *adjnl, int *adjproc, int *isbound, int ss, int *nColor, int *color, int **newcolored, int *mark, int gmaxdeg, G2LHash *hash, int **forbidden, int **xforbidden, int **forbiddenS, int **xforbiddenS, int *forbsize, int *forbsizeS, char color_method, char comm_pattern, int *rreqfromC, int *repliesC, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromF, int *repliesF, MPI_Request *sreqsF, MPI_Request *rreqsF, int *rreqfromFx, int *repliesFx, MPI_Request *sreqsFx, MPI_Request *rreqsFx, MPI_Status *stats, int *confChk, int *ssendsize, int *srecsize, int **ssendbuf, int **srecbuf, int ** ssp, int **srp, int **xfp, int *xfpMark, int *wset, int *wsize);
static int D1DetectConflicts(ZZ *, int, int *, int *, int *, int *, int *, int *, int *, G2LHash *);
static int D2DetectConflicts(ZZ *zz, G2LHash *hash, int nlvtx, int *wset, int wsize, int *xadj, int *adj, int *adjproc, int *nColor, int *color, int *conflicts, int *rand_key, int *vmark, int *seen, int *where, int *pwhere, int **rcsendbuf, int **rcrecbuf, int *rcsendsize, int *rcrecsize, int **srp, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromC, MPI_Status *stats, int *nconflict);
    
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
    

#if 0
static void PrintGraph(ZZ *zz, char *name, int base, int nvtx, int *xadj, int *adj, int *adjproc)
{
    int i, j;
    FILE *outf;
    char fname[256];
    static int callcnt=0;
   

    sprintf(fname, "umit-dbg.%d.txt", zz->Proc);
    outf = fopen(fname, (!callcnt) ? "w" : "a");
    fprintf(outf, "%s\n", name);
    for (i=0; i<nvtx; ++i) {
        fprintf(outf, "%d [%d] : ", base+i, i);
        for (j=xadj[i]; j<xadj[i+1]; ++j)
            fprintf(outf, "%d (%d), ", adj[j], adjproc[j]);
        fprintf(outf, "\n");
    }
    fclose(outf);
    ++callcnt;
}
#endif

    
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
  indextype *vtxdist, *xadj, *adjncy; /* arrays to store the graph structure */
  int *adjproc;                     
  int *input_parts;                 /* Initial partitions for objects. */
  int nvtx = num_obj;               /* number of vertices */
  int gvtx;                         /* number of global vertices */
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
  int comm[2],gcomm[2]; 

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
  comm_pattern = 'S';

  Zoltan_Assign_Param_Vals(zz->Params, Color_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);
  
  /* Check validity of parameters */
  if (distance != 1 && distance != 2) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Coloring with requested distance is not implemented. Using Distance-1 coloring.");
      distance = 1;
  }
  if (ss == 0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid superstep size. Using default value 100.");
      ss = 100;
  }
  if (comm_pattern != 'S' && comm_pattern != 'A') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid communication pattern. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (comm_pattern == 'A' && distance == 2) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Asynchronous communication pattern is not implemented for distance-2 coloring. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (color_order != 'I' && color_order != 'B' && color_order != 'U') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring order. Using internal first coloring order (I).");
      color_order = 'I';
  }  
  if (color_order == 'U' && distance == 2) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Interleaved coloring order is not implemented for distance-2 coloring. Using internal first coloring order (I).");
      color_order = 'I';
  }
  if (color_method !='F') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring method. Using first fit method (F).");
      color_method = 'F';
  }
  

  /* Compute Max number of array entries per ID over all processors.
     This is a sanity-maintaining step; we don't want different
     processors to have different values for these numbers. */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = *num_gid_entries = gcomm[0];
  zz->Num_LID = *num_lid_entries = gcomm[1];
  
  /* Return if this processor is not in the Zoltan structure's
     communicator. */
  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz))
      return ZOLTAN_OK;


  /* BUILD THE GRAPH */
  /* Check that the user has allocated space for the return args. */
  if (!color_exp)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Output argument is NULL. Please allocate all required arrays before calling this routine.");
  
  /* Initialize all local pointers to NULL. This is necessary
     because we free all non-NULL pointers upon errors. */
  vtxdist = xadj = adjncy = adjproc = NULL;
  ewgts = float_vwgt = NULL;
  input_parts = NULL;
  
  /* Default graph type is GLOBAL. */
  SET_GLOBAL_GRAPH(&graph_type);
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
  ierr = Zoltan_Build_Graph(zz, &graph_type, check_graph, nvtx,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_COLOR_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }

  /* CREATE THE HASH TABLE */
  /* Determine hash size and allocate hash table */
  /* UVCUVC: TODO we can allocate smaller hash; check this later */
  MPI_Allreduce(&nvtx, &gvtx, 1, MPI_INT, MPI_SUM, zz->Communicator);        
  i = gvtx > xadj[nvtx] ? 2*xadj[nvtx] : 2*gvtx; /* i is the minimum hash size */
  i = i > nvtx ? i : 2*nvtx;   
  if (Zoltan_GenPrime(i , &hsize)==ZOLTAN_MEMERR)
      MEMORY_ERROR;
  if (Zoltan_G2LHash_Create(&hash, hsize)==ZOLTAN_MEMERR)
      MEMORY_ERROR;

  /* Insert the ids of the local vertices into the hash table */
  /* By inserting local vertices first, it is ensured that the
     local ids of local vertices start from 0 and are consecutive*/

#if 0
  PrintGraph(zz, "Before Global-2-Local", vtxdist[zz->Proc], nvtx, xadj, adjncy, adjproc);
#endif
  
  for (j=0; j<nvtx; j++) 
      Zoltan_G2LHash_Insert(&hash, vtxdist[zz->Proc] + j);

  /* Add ids of the d1 neighbors into the hash table*/
  for (i=0; i<xadj[nvtx]; ++i)
      adjncy[i] = Zoltan_G2LHash_Insert(&hash, adjncy[i]);
  /* lastno is the total number of local and d1 neighbors */
  lastlno = hash.lastlno; 

#if 0
  PrintGraph(zz, "After Global-2-Local", vtxdist[zz->Proc], nvtx, xadj, adjncy, adjproc);
#endif

  /* Allocate color array. Local and D1 neighbor colors will be stored here. */ 
  if (!(color = (int *) ZOLTAN_CALLOC(lastlno, sizeof(int))))
      MEMORY_ERROR;

  /* SELECT COLORING ALGORITHM AND PERFORM COLORING*/
  if (distance == 1)      
      D1coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, &hash, xadj, adjncy, adjproc, color); 
  else if (distance == 2)
      D2coloring(zz, color_order, color_method, comm_pattern, ss, nvtx, &hash, xadj, adjncy, adjproc, color);
    
  /* FILL THE RETURN ARRAY */
  for (i=0; i<nvtx; i++) 
      color_exp[i] = color[i];

#if 0
  /* Check if there is an error in coloring */
  if (distance == 1 && zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) { 
      for (i=0; i<nvtx; i++) {
          for (j = xadj[i]; j < xadj[i+1]; j++) {
              int v = adjncy[j];
              if (color[i] == color[v])
                  printf("Error in coloring! u:%d(%d), v:%d(%d), cu:%d, cv:%d\n", i, Zoltan_G2LHash_L2G(&hash, i), v, Zoltan_G2LHash_L2G(&hash, v), color[i], color[v]);
          }
      }
  }
#endif
  
 End:  
  Zoltan_Multifree(__FILE__,__LINE__, 6, &vtxdist, &xadj, &adjncy, &input_parts, &adjproc, &color);
  Zoltan_G2LHash_Destroy(&hash);
  
  return ierr;
}

/*****************************************************************************/
/* Distance-1 coloring. No two adjacent vertices get the same color. */

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
    int distance = 1;
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
    int *visitIntern = NULL;           
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
#ifdef RELEVANT_COLORS
    int j, k;
    int *pmark=NULL;        
    int **persSbuf=NULL; /* personalized send buffers */
    int *Ssize=NULL; /* send buffer sizes */
    int *relproc=NULL;
    int *xrelproc=NULL;
    int *plst=NULL, plstcnt=0;
#endif

    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    if (!isbound || !visit)
        MEMORY_ERROR;
    
    /* Start timer */
    get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
    if (get_times){
        MPI_Barrier(zz->Communicator);
        times[0] = Zoltan_Time(zz->Timer);
    }

#if 0
    printf("Before reordering: nvtx:%d lastlno=%d Proc:%d\n", nvtx, lastlno, zz->Proc);
    for (i=0; i < nvtx; i++) {
        int j;
        printf("%d [%d] :: ", Zoltan_G2LHash_L2G(hash, i), i);
        for (j=xadj[i]; j < xadj[i+1]; j++) {
            printf("%d [%d] (%d) ",   Zoltan_G2LHash_L2G(hash, adj[j]), adj[j], adjproc[j]);
        }
        printf("\n");
    }
#endif
    
    ierr = ReorderGraph(zz, 1, nvtx, xadj, &xbadj, adj, adjproc, &nbound, isbound, visit);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
        ZOLTAN_COLOR_ERROR(ierr, "Error in ReorderGraph");
    if (get_times) times[1] = Zoltan_Time(zz->Timer);

#if 0  
    printf("After reordering: nvtx:%d lastlno=%d Proc:%d\n", nvtx, lastlno, zz->Proc);
    for (i=0; i < nvtx; i++) {
        int j;
        printf("%d [%d] :: ", Zoltan_G2LHash_L2G(hash, i), i);
        for (j=xadj[i]; j < xadj[i+1]; j++) {
            printf("%d [%d] (%d) ",   Zoltan_G2LHash_L2G(hash, adj[j]), adj[j], adjproc[j]);
        }
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
    mark = (int *) ZOLTAN_CALLOC(gmaxdeg > 0 ? gmaxdeg : 1, sizeof(int)); 
    conflicts = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    replies = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    stats = (MPI_Status *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Status));
    rreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    rreqfrom = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    newcolored = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    if (!mark || !conflicts || !replies || !stats || !rreqs || !sreqs || !rreqfrom || !newcolored)
        MEMORY_ERROR;
    for (i=0; i<zz->Num_Proc; ++i) {
        newcolored[i] = (int *) ZOLTAN_MALLOC(2 * ss * sizeof(int));
        if (!newcolored[i])
            MEMORY_ERROR;
	memset(newcolored[i], 0, 2 * ss * sizeof(int));
    }
#ifdef RELEVANT_COLORS
    persSbuf = (int **) ZOLTAN_MALLOC(sizeof(int *) * zz->Num_Proc);
    Ssize = (int *) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
    plst = (int *) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);    
    for (i=0; i<zz->Num_Proc; i++) 
        persSbuf[i] = (int *) ZOLTAN_MALLOC(sizeof(int) * 2*ss);
#endif
                                            
    /* Generate random numbers associated with global numbers of the vertices */
    /* All processors generate the same random number corresponding
       to the same global vertex number */
    rand_key = (int *) ZOLTAN_MALLOC(sizeof(int) * lastlno);
    if (!rand_key)
        MEMORY_ERROR;
    for(i=0; i<lastlno; i++) {
        Zoltan_Srand(Zoltan_G2LHash_L2G(hash, i), NULL);
        rand_key[i] = ((double)Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000;
    }
    
    /* Color internal vertices and determine the visit order */
    if (get_times) {
        MPI_Barrier(zz->Communicator);
        times[2] = Zoltan_Time(zz->Timer);
    }
    visitIntern = visit + nbound; /* Start of internal vertex visit order. Used with I and B options below */
    if (color_order == 'U') { 
        nConflict = nvtx;
        for (i=0; i<nvtx; i++)
            visit[i] = i;
        if (zz->Num_Proc==1)
            InternalColoring(zz, distance, &nColor, nvtx, visit, xadj, adj, color, mark, gmaxdeg, color_method);
    }
    else if (color_order == 'I') {
        InternalColoring(zz, distance, &nColor, nvtx - nbound, visitIntern, xadj, adj, color, mark, gmaxdeg, color_method);
        nConflict = nbound;
    }
    else if (color_order == 'B')
        nConflict = nbound;

    if (get_times) times[3] = Zoltan_Time(zz->Timer);

    /* Color boundary vertices */
    nRound = 0; 
    nTotConflict = 0;

    if (zz->Num_Proc >= 2) {

#ifdef RELEVANT_COLORS
        /*** Determine which processors each vertex is connected to ***/
        /*** Updated color info of a vertex will only be sent to the connected processors ***/
        pmark = (int *) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
        memset(pmark, 0xff, sizeof(int) * zz->Num_Proc);
        relproc = (int *) ZOLTAN_MALLOC(sizeof(int) * xadj[nvtx]); /* UVC: TODO check this might be too large */
        /* more than or equal to the requirement, careful about the "unordered" case if optimizing the above allocation. */
        xrelproc = (int *) ZOLTAN_MALLOC(sizeof(int) * (nConflict+1));
        
        xrelproc[0] = 0;
        k = 0;
        for (i=0; i<nConflict; i++) {
            int u = visit[i];
            isbound[u] = 1+i;  /* we use isbound for indirection for boundary vertices to xrelproc array */
            pmark[zz->Proc] = u; /* mark myself so I don't send message to myself */
            for (j=xadj[u]; j<xbadj[u]; j++) {
                int pv = adjproc[j];
                if (pmark[pv]!=u) {
                    relproc[k++] = pv;
                    pmark[pv] = u;
                }
            }
            xrelproc[i+1] = k;
        }
        memset(pmark, 0, sizeof(int) * zz->Num_Proc);
        pmark[zz->Proc] = 1;
        for (i=0; i<nConflict; i++) {
            int u = visit[i];

            for (j=xadj[u]; j<xbadj[u]; j++) {
                int pv = adjproc[j];
                if (!pmark[pv]) {
                    pmark[pv] = 1;
                    plst[plstcnt++] = pv;
                }
            }
        }

#if 0
        printf("[%d] #of relevalnt procs are %d and procs are: ", zz->Proc, plstcnt);
        for (i=0; i<plstcnt; ++i)
            printf("%d ", plst[i]);
        printf("for vertices related procs are:\n");
        for (i=0; i<nConflict; i++) {
            int u = visit[i];
            printf("%d: (", u);
            for (j=xrelproc[isbound[u]-1]; j<xrelproc[isbound[u]]; j++)
                printf("%d ", relproc[j]);
            printf(") \n full adjproc: ");
            for (j=xadj[u]; j<xadj[u+1]; j++)
                printf("%d ", adjproc[j]);
            printf("\n");            
        }
        printf("\n");
#endif
                    
        
#endif
        
        do {
            int *tp = visit;
            memset(mark, 0xff, (1+nColor) * sizeof(int)); /* reset dirty entries */
#ifdef RELEVANT_COLORS
            ierr = D1ParallelColoring(zz, nConflict, visit, xadj, adj, isbound, ss,
                                      &nColor, color, newcolored, mark, gmaxdeg, hash,
                                      color_method, comm_pattern, rreqfrom, replies,
                                      sreqs, rreqs, stats,
                                      xrelproc, relproc, persSbuf, Ssize, plstcnt, plst);
#else
            ierr = D1ParallelColoring(zz, nConflict, visit, xadj, adj, isbound, ss,
                                      &nColor, color, newcolored, mark, gmaxdeg, hash,
                                      color_method, comm_pattern, rreqfrom, replies,
                                      sreqs, rreqs, stats);
#endif
            if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                ZOLTAN_COLOR_ERROR(ierr, "Error in D1ParallelColoring");
            nConflict = D1DetectConflicts(zz, nConflict, visit, xadj, xbadj, adj,
                                          color, conflicts, rand_key, hash);

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
        InternalColoring(zz, distance, &nColor, nvtx-nbound, visitIntern, xadj, adj, color, mark, gmaxdeg, color_method);
    
#if 0
    printf("[%d] vtx(gno)->color: ", zz->Proc);
    for (i=0; i<nvtx; i++)
        printf("%d(%d)->%d ", i, Zoltan_G2LHash_L2G(hash, i), color[i]);
    printf("\n");
#endif
    
    
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
                                        

#ifdef RELEVANT_COLORS
    for (i=0; i<zz->Num_Proc; i++) 
        ZOLTAN_FREE(&persSbuf[i]);
    ZOLTAN_FREE(&persSbuf);
    ZOLTAN_FREE(&Ssize);

    ZOLTAN_FREE(&plst);    
    ZOLTAN_FREE(&pmark);
    ZOLTAN_FREE(&relproc);
    ZOLTAN_FREE(&xrelproc);
#endif

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
/* Distance-2 Coloring */

static int D2coloring(
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
    static char *yo = "D2coloring";
    int i, j, p;
    int distance = 2;
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
    int nbound;                  /* Number of local boundary1 vertices */
    int *xbadj = NULL;           /* Pointer to end of cut edges in adj lists
                                    of local vertices */
    int *isbound=NULL;           /* Indicates whether a local vertex is a boundary
                                    vertex */
    int *visit = NULL;           /* Visit (coloring) order of local vertices */
    int *visitIntern = NULL;       
    int *mark = NULL;            /* Array to mark forbidden colors for a vertex */
    int gmaxcolor = 0;           /* Maximum #colors */
    int lmaxdeg = 0;             /* Maximum vertex degree for the local vertices */
    int lastlno = hash->lastlno; /* Total number of local and D1 neighbor vertices */
    int *conflicts = NULL;       /* List of vertices to be recolored in the next
                                    round */
    /* Arrays used for MPI communication */
    int *repliesF = NULL, *repliesFx = NULL, *repliesC = NULL;
    MPI_Status *stats = NULL;
    MPI_Request *rreqsF = NULL, *rreqsFx = NULL, *rreqsC = NULL;
    MPI_Request *sreqsF = NULL, *sreqsFx = NULL, *sreqsC = NULL;
    int *rreqfromF = NULL, *rreqfromFx = NULL, *rreqfromC = NULL;
    int rreqcntC, sreqcntC;
    int **newcolored = NULL;     /* Array used for communicating boundary vertex
                                    colors at the end of supersteps. Global number
                                    of the vertex - color of vertex pairs are filled
                                    in this array */
    double times[6];             /* Used for timing measurements */
    int get_times;               /* (1) Measure timings (0) Don't */ 
    int ierr;
    /* conflict detection */
    int *confChk=NULL, *seen=NULL, *where=NULL, *pwhere=NULL, *vmark=NULL;
    int *wset, wsize=0;    
    /* forbidden color information exchange */
    int **xfp=NULL, *xfpMark;
    int **forbidden=NULL, *forbsize = NULL, **xforbidden=NULL; 
    int **forbiddenS=NULL, *forbsizeS = NULL, **xforbiddenS=NULL;
    /* storing partial adjacency lists of nonlocal vertices */
    int nnl=0, *xadjnl=NULL, *adjnl=NULL, *xadjnltemp=NULL;
    /* superstep size info exchange */
    int **srecbuf = NULL, **ssendbuf=NULL;
    int *srecsize=NULL, *ssendsize=NULL;
    int **ssp, **srp;
    /* vertex recolor info exchange */
    int **rcsendbuf = NULL, **rcrecbuf=NULL;
    int *rcsendsize=NULL, *rcrecsize=NULL;
    
    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    if (!isbound || !visit)
        MEMORY_ERROR;
    
    /* Start timer */
    get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
    if (get_times){
        MPI_Barrier(zz->Communicator);
        times[0] = Zoltan_Time(zz->Timer);
    }    
    ierr = ReorderGraph(zz, 2, nvtx, xadj, &xbadj, adj, adjproc, &nbound, isbound, visit);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
        ZOLTAN_COLOR_ERROR(ierr, "Error in ReorderGraph");
    if (get_times) times[1] = Zoltan_Time(zz->Timer);

#if 0
    printf("After reordering: nvtx:%d lastlno=%d Proc:%d\n", nvtx, lastlno, zz->Proc);
    for (i=0; i < nvtx; i++) {
        printf("%d (%d) :: ", Zoltan_G2LHash_L2G(hash, i), i);
        for (j=xadj[i]; j < xadj[i+1]; j++) {
            if (j == xbadj[i])
                printf(" | ");
            printf("%d (%d) ", Zoltan_G2LHash_L2G(hash, adj[j]), adj[j]);
        }
        printf("\n");
    }
#endif
    
    /* Calculate the maximum degree of the graph */
    lmaxdeg = 0;
    for (i=0; i<nvtx; i++)
        if (lmaxdeg < xadj[i+1] - xadj[i])
            lmaxdeg = xadj[i+1] - xadj[i];    
    MPI_Allreduce(&lmaxdeg, &gmaxcolor, 1, MPI_INT, MPI_MAX, zz->Communicator);        
    /* gmaxdeg^2+1 is the upper bound for #colors */
    gmaxcolor = gmaxcolor*gmaxcolor+1; 


    /***** CONSTRUCT PARTIAL ADJ LIST OF NON-LOCAL VERTICES *****/ 
    /* Count the number of times a non-local vertex appears in adj list */
    /* Id's of nonlocal vertices start from nvtx and ends at lastlno - 1 due to hashing*/
    nnl = lastlno - nvtx;
    if (!(xadjnl = (int *) ZOLTAN_CALLOC(nnl + 1, sizeof(int))))
        MEMORY_ERROR;
    if (!(xadjnltemp = (int *) ZOLTAN_MALLOC((nnl + 1) * sizeof(int))))
        MEMORY_ERROR;
    xadjnl[0] = 0;    
    for (i=0; i<nvtx; i++) {
        /* consider the cut edges only */
        for (j=xadj[i]; j<xbadj[i]; j++) 
            ++xadjnl[adj[j]-nvtx+1];
    }
    for (i=0; i<nnl; i++)
        xadjnl[i+1] += xadjnl[i];

    /* Construct the partial adjancency list of non-local vertices */
    /* Only local vertices appear in the partial adj list */
    memcpy(xadjnltemp, xadjnl, (nnl + 1) * sizeof(int));
    if (!(adjnl = (int *) ZOLTAN_CALLOC(xadjnl[nnl] > 0 ? xadjnl[nnl] : 1, sizeof(int))))
        MEMORY_ERROR;    
    for (i=0; i<nvtx; i++) {
        /* consider the cut edges only */
        for (j=xadj[i]; j<xbadj[i]; j++) {
            int idx = xadjnltemp[adj[j]-nvtx]++;
            adjnl[idx] = i;
        }
    }       

#if 0
    printf("[%d] Partial adj list of non-local vertices: #nonlocal: %d\n", zz->Proc, nnl);
    for (i=0; i<nnl; i++) {
        printf("%d :: ", i+nvtx);
        for (j=xadjnl[i]; j<xadjnl[i+1]; j++)
            printf("%d ", adjnl[j]);
        printf("\n");
    }
#endif
    
    
    /***** MEMORY ALLOCATON *****/
    stats = (MPI_Status *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Status));
    rreqsC = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqsC = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    rreqfromC = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    repliesC = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    srecsize = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    ssendsize = (int *) ZOLTAN_CALLOC(zz->Num_Proc, sizeof(int));
    rcsendsize = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    rcrecsize = (int *) ZOLTAN_CALLOC(zz->Num_Proc, sizeof(int));
    if (!stats || !repliesC || !rreqsC || !sreqsC || !rreqfromC || !srecsize || !ssendsize || !rcsendsize || !rcrecsize)
        MEMORY_ERROR;    
    /* Issue async recvs for superstep info size */
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p) 
        if (p!=zz->Proc)
            MPI_Irecv(&srecsize[p], 1, MPI_INT, p, SNTAG, zz->Communicator, &rreqsC[rreqcntC++]);
    /* Calculate size of superstep info to be sent to each processor*/
    for (i=0; i<nvtx; ++i) {
        /* consider cut edges only */
        for (j=xadj[i]; j<xbadj[i]; j++)
            ++ssendsize[adjproc[j]]; 
    }    
    for (i=0; i<zz->Num_Proc; i++)
        if (i != zz->Proc)
            ssendsize[i] += ceil((double)nbound / (double)ss) + 1;
    /* Send the superstep info size so that other processors will allocate enough space */
    for (sreqcntC=p=0; p<zz->Num_Proc; ++p) 
        if (p != zz->Proc) 
            MPI_Isend(&ssendsize[p], 1, MPI_INT, p, SNTAG, zz->Communicator, &sreqsC[sreqcntC++]);
    
    mark = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int)); 
    vmark = (int *) ZOLTAN_CALLOC(lastlno, sizeof(int));
    conflicts = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    if (!mark || !conflicts || !vmark)
        MEMORY_ERROR;
    repliesF = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    repliesFx = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    rreqsF = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqsF = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    rreqsFx = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqsFx = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc  *sizeof(MPI_Request));
    rreqfromF = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    rreqfromFx = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    if (!repliesF || !rreqsF || !sreqsF || !rreqfromF || 
        !repliesFx || !rreqsFx || !sreqsFx || !rreqfromFx)
        MEMORY_ERROR;
    forbsize = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    forbsizeS = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    xfp = (int**) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));    
    xfpMark = (int*) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));    
    newcolored = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    forbidden = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    xforbidden = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    forbiddenS = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    xforbiddenS = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    if (!forbsize || !forbsizeS || !xfp || !xfpMark || !newcolored || !forbidden || 
        !xforbidden || !forbiddenS || !xforbiddenS)
        MEMORY_ERROR;
    for (i=0; i<zz->Num_Proc; ++i) {
        newcolored[i] = (int *) ZOLTAN_MALLOC(2 * ss * sizeof(int)); /* DB: can be reduced */
        xforbidden[i] = (int *) ZOLTAN_MALLOC((ss+1) * sizeof(int));
        forbsize[i] = (int)(ceil((double) xadj[nvtx] / (double)nvtx) * ss); /* size of forbidden color allocation - init to avgDeg*ss */
        forbidden[i] = (int *) ZOLTAN_MALLOC((forbsize[i] > 0 ? forbsize[i] : 1) * sizeof(int));
        xforbiddenS[i] = (int *) ZOLTAN_MALLOC((ss+1) * sizeof(int));
        forbsizeS[i] = (int)(ceil((double) xadj[nvtx] / (double) nvtx) * ss); /* size of forbidden color allocation for send buffer - init to avgDeg*ss */
        forbiddenS[i] = (int *) ZOLTAN_MALLOC((forbsizeS[i] > 0 ? forbsizeS[i] : 1) * sizeof(int));
        if (!newcolored[i] || !xforbidden[i] || !forbidden[i] || !xforbiddenS[i] || !forbiddenS[i])
            MEMORY_ERROR;
    }
    confChk = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int)); /*the vertices to be checked in conflict detection */
    wset = (int *) ZOLTAN_MALLOC((nbound > 0 ? nbound : 1) * sizeof(int));
    seen = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    where = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    pwhere = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    if (!confChk || !wset || !seen || !where || !pwhere)
        MEMORY_ERROR;
    
    /* Wait for superstep size communication to end */
    MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all superstep numbers to be received */
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all superstep numbers to be sent */ 

    /* Allocate buffers for superstep size and recolor info to be sent and received */
    ssendbuf = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    srecbuf = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    rcrecbuf = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    rcsendbuf = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    ssp = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    srp = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    if (!ssendbuf || !srecbuf || !ssp || !srp || !rcsendbuf || !rcrecbuf)
        MEMORY_ERROR;
    for (p=0; p<zz->Num_Proc; p++) 
        if (p != zz->Proc) {
            ssendbuf[p] = (int *) ZOLTAN_MALLOC((ssendsize[p] > 0 ? ssendsize[p] : 1) * sizeof(int));
            srecbuf[p] = (int *) ZOLTAN_MALLOC((srecsize[p] > 0 ? srecsize[p] : 1) * sizeof(int));
            rcrecsize[p] = 2 * ssendsize[p]; /* recolor info exchange buffer sizes. */ 
            rcsendsize[p] = 2 * srecsize[p]; /* better estimates or realloc in D2DetectConflicts can be used */
            rcrecbuf[p] = (int *) ZOLTAN_MALLOC((rcrecsize[p] > 0 ? rcrecsize[p] : 1) * sizeof(int));
            rcsendbuf[p] = (int *) ZOLTAN_MALLOC((rcsendsize[p] > 0 ? rcsendsize[p] : 1) * sizeof(int));
            if (!srecbuf[p] || !ssendbuf[p] || !rcsendbuf[p] || !rcrecbuf[p])
                MEMORY_ERROR;
        } else {
            ssendbuf[p] = NULL;
            srecbuf[p] = NULL;
            rcrecbuf[p] = NULL;
            rcsendbuf[p] = NULL;
        }
    
    /* Generate random numbers associated with global numbers of the vertices */
    /* All processors generate the same random number corresponding
       to the same global vertex number */
    rand_key = (int *) ZOLTAN_MALLOC(sizeof(int) * lastlno);
    if (!rand_key)
        MEMORY_ERROR;
    for(i=0; i<lastlno; i++) {
        Zoltan_Srand(Zoltan_G2LHash_L2G(hash, i), NULL);
        rand_key[i] = ((double)Zoltan_Rand(NULL)/(double)ZOLTAN_RAND_MAX)*1000000;
    }   

    /* Color internal vertices and determine the visit order */
    if (get_times) {
        MPI_Barrier(zz->Communicator);
        times[2] = Zoltan_Time(zz->Timer);
    }
    visitIntern = visit + nbound; /* Start of internal vertex visit order. Used with I and B options below */
    if (color_order == 'I') {
        InternalColoring(zz, distance, &nColor, nvtx - nbound, visitIntern, xadj, adj, color, mark, gmaxcolor, color_method);
        nConflict = nbound;
    }
    else if (color_order=='B') 
        nConflict = nbound;
/*  
    else if (color_order == 'U') { ** not implemented **
        nConflict = nvtx;
        for (i=0; i<nvtx; i++)
            visit[i] = i;
        if (zz->Num_Proc==1)
            InternalColoring(zz, distance, &nColor, nvtx, visit, xadj, adj, color, mark, gmaxcolor, color_method);
    }
*/
    
    if (get_times) times[3] = Zoltan_Time(zz->Timer);

    /* Color boundary vertices */
    nTotConflict = 0;
    nRound = 0;    
    if (zz->Num_Proc >= 2) {
        do {
            int *tp = visit, lColor;
            wsize = 0;
            memset(mark, 0xff, (nColor+1) * sizeof(int));
            memset(xfpMark, 0xff, zz->Num_Proc * sizeof(int));
            memset(confChk, 0xff, nvtx * sizeof(int)); /* UVCUVC: Check if we can do better;
                                                          i.e. not to re-set again in every round */
            ierr = D2ParallelColoring(zz, nConflict, nvtx, visit, xadj, adj, xbadj, xadjnl, adjnl, adjproc, isbound, ss,
                                      &nColor, color, newcolored, mark, gmaxcolor, hash,
                                      forbidden, xforbidden, forbiddenS, xforbiddenS, forbsize, forbsizeS,
                                      color_method, comm_pattern,
                                      rreqfromC, repliesC, sreqsC, rreqsC,
                                      rreqfromF, repliesF, sreqsF, rreqsF,
                                      rreqfromFx, repliesFx, sreqsFx, rreqsFx,
                                      stats, confChk, ssendsize, srecsize, ssendbuf, srecbuf, ssp, srp, xfp, xfpMark, wset, &wsize);
            if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                ZOLTAN_COLOR_ERROR(ierr, "Error in D2ParallelColoring");
            /* Allreduce below avoids insufficient memset in conflict detection. 
            If a received color of a non-local vertex is greater than nColor on this processor,
            then arrays accessed with color number will cause seg fault if this is not used. We may do better */
            lColor = nColor; 
            MPI_Allreduce(&lColor, &nColor, 1, MPI_INT, MPI_MAX, zz->Communicator);        

            ierr = D2DetectConflicts(zz, hash, nvtx, wset, wsize, xadj, adj, adjproc, &nColor, color, conflicts, rand_key,
                              vmark, seen, where, pwhere, rcsendbuf, rcrecbuf, rcsendsize, rcrecsize, srp, sreqsC, rreqsC,
                              rreqfromC, stats, &nConflict);
            if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                ZOLTAN_COLOR_ERROR(ierr, "Error in D2DetectConflicts");
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
        InternalColoring(zz, distance, &nColor, nvtx-nbound, visitIntern, xadj, adj, color, mark, gmaxcolor, color_method);
    
    
#if 0
    printf("[%d] vtx(gno)->color: ", zz->Proc);
    for (i=0; i<nvtx; i++)
        printf("%d(%d)->%d ", i, Zoltan_G2LHash_L2G(hash, i), color[i]);
    printf("\n");
#endif


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
    
    /* Free */
    for (i=0; i<zz->Num_Proc; ++i) 
        Zoltan_Multifree (__FILE__, __LINE__, 9, &newcolored[i], &forbidden[i],
                          &xforbidden[i], &forbiddenS[i], &xforbiddenS[i],
                          &srecbuf[i], &ssendbuf[i], &rcsendbuf[i], &rcrecbuf[i]);

    Zoltan_Multifree (__FILE__, __LINE__, 47, &isbound, &visit, &conflicts,
                      &mark, &vmark, &newcolored, &repliesF, &repliesC, &repliesFx, &stats,
                      &rreqsC, &rreqsF, &rreqsFx, &sreqsC,
                      &sreqsF, &sreqsFx, &rreqfromC, &rreqfromF,
                      &rreqfromFx, &xbadj, &srecsize, &ssendsize, &rcsendsize, &rcrecsize,
                      &srecbuf, &ssendbuf, &rcsendbuf, &rcrecbuf,
                      &forbidden, &xforbidden, &forbiddenS,
                      &xforbiddenS, &forbsize, &forbsizeS,  &confChk, &seen,
                      &where, &pwhere, &xfp, &xfpMark, &rand_key, &wset, &xadjnl, &adjnl, &xadjnltemp, &ssp, &srp);
    
    return ierr;
}


/*****************************************************************************/
/* Reorder graph adjacency list */

static int ReorderGraph(
    ZZ *zz,           
    int distance,     /* Coloring type */
    int nvtx,         /* In: Number of objects */
    int *xadj,        /* In: Pointer to xadj */
    int **pxbadj,     /* Out: Pointer to start of cut edges in the
                         adj lists of boundary vertices*/
    int *adj,         /* In/Out: Pointer to adj */
    int *adjproc,     /* In/Out: parallel to adj; holds the processor owns the adj vertex */
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
        MEMORY_ERROR;
    
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
    /* Construct vertex visit order: boundaries first, internals last */
    idx = 0;
    nboundIdx = *nbound;
    for (i=0; i<nvtx; ++i) {
        if (isbound[i] == 1)
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
         
    switch(color_method){    
    case 'F':
    default: 
        for (c = 1; (c <= *nColor) && (mark[c] == u); ++c) ;
        if (c > *nColor) { /* no available color with # less than nColor */
            c = ++(*nColor);
            mark[c] = -1;
        }
        pickedColor = c;
    }
    
    return pickedColor;   
}

/*****************************************************************************/
/* Color internal vertices */

static int InternalColoring(
    ZZ *zz,
    int distance,
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
    int i, j, k, c, u, v, w;
    static char *yo = "InternalColoring";
    int ierr = ZOLTAN_OK;
 
    memset(mark, 0xff, gmaxdeg * sizeof(int));
    if (distance == 1) {
        for (i=0; i<nvtx; ++i) {
            u = visit[i];        
            for (j = xadj[u]; j < xadj[u+1]; ++j) {
                v = adj[j];            
                if ((c = color[v]) != 0) 
                    mark[c] = u;
            }
            color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);
        }
    } else if (distance == 2) {        
        for (i=0; i<nvtx; ++i) {
            u = visit[i];        
            for (j = xadj[u]; j < xadj[u+1]; ++j) {
                v = adj[j];              
                if ((c = color[v]) != 0) 
                    mark[c] = u;
                if (v < nvtx) /* neighbor vertex v is internal as well - not sure if this is true for boundary first d2 coloring DBDBcheck*/ 
                    for (k = xadj[v]; k < xadj[v+1]; ++k) {	  
                        w = adj[k];
                        if ((c = color[w]) != 0) 
                            mark[c] = u;	      
                    }                
            }
            color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);
        }
    } else
        ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Unknown coloring distance");

 End:
    
    return ierr;
}

/*****************************************************************************/
/* Parallel coloring of boundary vertices */       

#ifdef RELEVANT_COLORS    
static int D1ParallelColoring (ZZ *zz, int nvtx, int *visit, int *xadj, int *adj,
                               int *isbound, int ss, int *nColor, int *color,
                               int **newcolored, int *mark, int gmaxdeg, G2LHash *hash,
                               char color_method, char comm_pattern, int *rreqfrom,
                               int *replies, MPI_Request *sreqs, MPI_Request *rreqs,
                               MPI_Status *stats,
                               int *xrelproc, int *relproc, int **persSbuf, int *Ssize, int plstcnt, int *plst)
#else
static int D1ParallelColoring (
    ZZ *zz,
    int nvtx,         /* number of vertices to be colored in this round */
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
    int *rreqfrom,    /* Arrays used for MPI communication */
    int *replies,
    MPI_Request *sreqs,
    MPI_Request *rreqs,
    MPI_Status *stats
)
#endif
{
    static char *yo="D1ParallelColoring";
    int colortag=1001, i, j, p, q, l;
    int *colored, n=0;
    int rreqcnt=0, sreqcnt=0, repcount;
    int ierr;

#ifdef RELEVANT_COLORS
    memset(Ssize, 0, sizeof(int) * zz->Num_Proc);
#endif
    
    colored = newcolored[zz->Proc];
    
    /* Issue async recvs */
#ifdef RELEVANT_COLORS
    for (rreqcnt = i = 0; i < plstcnt; ++i) {
        p = plst[i];
        rreqfrom[rreqcnt] = p;
        if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[rreqcnt]))
            ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
        ++rreqcnt;
    }
#else        
    for (rreqcnt = p = 0; p < zz->Num_Proc; ++p) 
        if (p != zz->Proc) {
            rreqfrom[rreqcnt] = p;
            if (MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, colortag, zz->Communicator, &rreqs[rreqcnt]))
                ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
            ++rreqcnt;
        }
#endif
    
    /* Coloring */
    for (i=0; i<nvtx; ++i) {
        int u = visit[i], gu;        
        for (j=xadj[u]; j<xadj[u+1]; ++j) {
            int gv = adj[j], c;
            if ((c = color[gv]) != 0) {
                while (*nColor < c) /* if nColor < c, set nColor to c and reset intermediate mark entries. May be done better. */
                    mark[++(*nColor)] = -1;
                mark[c] = u;
            }
        }
        color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);

        if (!isbound[u]) /* send only boundary vertices */
            continue;

        gu = Zoltan_G2LHash_L2G(hash, u);
#ifdef RELEVANT_COLORS
        for (j=xrelproc[isbound[u]-1]; j<xrelproc[isbound[u]]; ++j) {
            int ap = relproc[j];
            persSbuf[ap][Ssize[ap]++] = gu;
            persSbuf[ap][Ssize[ap]++] = color[u];
        }
        n += 2;
#else
            
        colored[n++] = gu;
        colored[n++] = color[u];
#endif
        
	/* If superstep is finished, communicate */
        if (n >= 2*ss) {
#ifdef RELEVANT_COLORS
            for (sreqcnt = j = 0; j<plstcnt; ++j) {
                p = plst[j];
                if (Ssize[p] < 2*ss) {
                    persSbuf[p][Ssize[p]++] = -2;
                    persSbuf[p][Ssize[p]++] = -2;
                }
                MPI_Isend(persSbuf[p], Ssize[p], MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
                Ssize[p] = 0;
            }
#else            
            for (sreqcnt = p = 0; p < zz->Num_Proc; ++p)
                if (p != zz->Proc) {
                    MPI_Isend(colored, 2*ss, MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
                    ++sreqcnt;
                }
#endif
            
            /* Wait color lists from other processors */
            if (comm_pattern == 'S') { /* wait all results if sync communication*/
                MPI_Waitall(rreqcnt, rreqs, stats); 
                repcount = rreqcnt;
            }
            else /* if async communication: check if any arrivied */
                MPI_Waitsome(rreqcnt, rreqs, &repcount, replies, stats); 

            for (l = repcount-1; l >= 0; --l) {
                int v=0;
                
                if (comm_pattern == 'S') 
                    q = l;                 
                else
                    q = replies[l];                
                p = rreqfrom[q];

                /* Read received color list from p */
                for (j = 0; j < 2*ss; ) {
                    int c, hv;
                    v  = newcolored[p][j++];
                    if (v < 0) 
                        break;                    
                    if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
                        c = newcolored[p][j++];
                        color[hv] = c;
                    } else
                        ++j;                    
                }
                /* If p hasn't finished coloring, issue new color request */
                if (v!=-1) {
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
#ifdef RELEVANT_COLORS        
    for (j=0; j<plstcnt; ++j) {
        p = plst[j];
        persSbuf[p][Ssize[p]++] = -1;
        persSbuf[p][Ssize[p]++] = -1;
    }
#else    
    colored[n++] = -1;
    colored[n++] = -1;
#endif

#ifdef RELEVANT_COLORS
    for (sreqcnt = j = 0; j<plstcnt; ++j) {
        p = plst[j];
        MPI_Isend(persSbuf[p], Ssize[p], MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
    }
#else            
    for (sreqcnt = p = 0; p < zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            MPI_Isend(colored, 2*ss, MPI_INT, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
            ++sreqcnt;
        }
#endif

    /* Receive the remaining color lists from other processors */
    while (rreqcnt) {
        if (comm_pattern == 'S') { /* wait all results if sync communication*/
            MPI_Waitall(rreqcnt, rreqs, stats); /* wait all results */
            repcount = rreqcnt;
        }
        else
            MPI_Waitsome(rreqcnt, rreqs, &repcount, replies, stats); /* wait some results */

        for (l=repcount-1; l>=0; --l) {
            int v=0;
            if (comm_pattern == 'S') /* wait all results if sync communication*/
                q = l;                
            else
                q = replies[l];
            p = rreqfrom[q];
            
            /* Read received color list from p */
            for (j = 0; j < 2*ss; ) {
                int c, hv;
                v  = newcolored[p][j++];
                if (v < 0) 
                    break;
                if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
                    c = newcolored[p][j++];
                    color[hv] = c;
                } else
                    ++j;
            }
            
            /* If p hasn't finished coloring, issue new color request */
            if (v!=-1) {
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
/* Send forbidden colors for the vertices to be colored in the next superstep */       

static int sendNextStepsForbiddenColors(ZZ *zz, G2LHash *hash, int nlvtx, int p, int **srp, int *xadj, int *adj, int *xadjnl, int *adjnl,
                                 int *forbsizeS, int **xforbiddenS, int **forbiddenS,
                                 int *nColor, int *color, int *mark, int *confChk, int sncnt, int *wset, int *wsize, MPI_Request *sreqsFx, MPI_Request *sreqsF)
{
    int u, m1 = 0, m2 = 0, m3, m4;
    int cnt;

    /*** calculate forbidden colors for each vertex u
         that proc p will color in its current superstep ***/
    while (*srp[p] >= 0) {
        u = Zoltan_G2LHash_G2L(hash, *srp[p]);
        cnt = 0; /* count the max number of forbidden colors for u */
        /* put the start of u's forbidden colors into xforbiddenS */
        xforbiddenS[p][m1++] = m2;
        for (m3 = xadjnl[u-nlvtx]; m3 < xadjnl[u-nlvtx+1]; m3++) { /*mark d2 forbidden colors of u*/
            int v = adjnl[m3]; /* v is a local neighbor of u*/
            /* mark v for conflict detection if it has more than 1 neighbors colored in this superstep */
            if (confChk[v] == sncnt) {
                confChk[v] = -2; /*mark for conflict check*/
                wset[(*wsize)++] = v;
            }
            else if (confChk[v] != -2)
                confChk[v] = sncnt; /*mark v the first time for this superstep*/
            /*resume calculating forbidden colors*/
            for (m4 = xadj[v]; m4 < xadj[v+1]; m4++) {
                int w = adj[m4], c; /* w is a d2 neighbor of u*/
                if ((c = color[w]) != 0) {
                    mark[c] = u;
                    cnt++;
                }
            }
        }
        if (forbsizeS[p] < m2 + cnt) { /*if send buffer is not large enough, increase buffer by 20% */
            forbsizeS[p] = (int)((double)(m2+cnt) * 1.2); 
            forbiddenS[p] = (int *) ZOLTAN_REALLOC(forbiddenS[p], sizeof(int) * forbsizeS[p]); 
        }
        for (m3 = 1; m3 <= *nColor; m3++) {
            if (mark[m3] == u) {
                forbiddenS[p][m2++] = m3;
            }
        }
        ++srp[p];
    }
    if (*srp[p] != -2)
        ++srp[p];
    xforbiddenS[p][m1++] = xforbiddenS[p][0] = m2;
    
    /*** send forbidden colors to p for its next superstep***/
    MPI_Isend(xforbiddenS[p], m1, MPI_INT, p, XFORBIDTAG, zz->Communicator, &sreqsFx[p]);
    MPI_Isend(forbiddenS[p], xforbiddenS[p][0], MPI_INT, p, FORBIDTAG, zz->Communicator, &sreqsF[p]);
    
    return 0;                    
}

/*****************************************************************************/
/* Wait until receiving all forbidden color pointers and then forbidden colors */       

static int waitPtrAndForbiddenColors(ZZ *zz, int rreqcntFx, MPI_Request *rreqsFx, int *rreqfromFx, MPI_Status *stats, int *forbsize, int **xforbidden, int **forbidden, int **xfp, MPI_Request *rreqsF)
{
    int rreqcntF=0, l, p;
        
    MPI_Waitall(rreqcntFx, rreqsFx, stats); /* wait all forbidden color pointers */ 
    for (l=rreqcntFx-1; l>=0; --l) { /* process the received messages */
        p = rreqfromFx[l];	  
        if (forbsize[p] < xforbidden[p][0]) { /*xforbidden[p][0] is the buffersize for msg from p*/
            forbidden[p] = (int *) ZOLTAN_REALLOC(forbidden[p], sizeof(int) * xforbidden[p][0]); 
            forbsize[p] = xforbidden[p][0];
        }
        /*** Issue async recvs for corresponding forbidden colors ***/
        MPI_Irecv(forbidden[p], xforbidden[p][0], MPI_INT, p, FORBIDTAG, zz->Communicator, &rreqsF[rreqcntF++]);
    }        
    MPI_Waitall(rreqcntF, rreqsF, stats); /* wait all forbidden color messages to be received */
    for (p = 0; p < zz->Num_Proc; p++) { /*reset the xforbidden vertex pointers to point the first vertices*/
        if (p != zz->Proc) {
            xforbidden[p][0] = 0;
            xfp[p] = xforbidden[p];
        }
    }
    return rreqcntF;
}



/*****************************************************************************/
/* Parallel coloring of boundary vertices */       

static int D2ParallelColoring (
    ZZ *zz,
    int nvtx,         /* number of vertices to be colored in this round */
    int nlvtx,
    int *visit,       /* Visit (coloring) order of local vertices */
    int *xadj,        /* arrays that store the graph structure */
    int *adj,
    int *xbadj,
    int *xadjnl,
    int *adjnl,
    int *adjproc,
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
    int **forbidden,  /* Arrays used for forbidden color exchange */
    int **xforbidden,
    int **forbiddenS,
    int **xforbiddenS,
    int *forbsize,
    int *forbsizeS,
    char color_method, /* Coloring method. (F) First fit
                          (S) staggered first fit (L) load balancing */
    char comm_pattern, /* (A) asynchronous (S) synchronous supersteps */
    int *rreqfromC,     /* Arrays used for MPI communication */
    int *repliesC,
    MPI_Request *sreqsC,
    MPI_Request *rreqsC,
    int *rreqfromF,     /* Arrays used for MPI communication */
    int *repliesF,
    MPI_Request *sreqsF,
    MPI_Request *rreqsF,
    int *rreqfromFx,     /* Arrays used for MPI communication */
    int *repliesFx,
    MPI_Request *sreqsFx,
    MPI_Request *rreqsFx,
    MPI_Status *stats,
    int *confChk,        /* Arrays used for conflict detection */
    int *ssendsize,
    int *srecsize,
    int **ssendbuf,
    int **srecbuf,
    int **ssp,
    int **srp, 
    int **xfp,
    int *xfpMark,
    int *wset,
    int *wsize
)
{    
    int i, j, k, p, l;
    int *colored = newcolored[zz->Proc];
    int n=0, m = 0;
    int rreqcntC=0, sreqcntC=0; /*keeps track of updated D-1 neighbor colors comm.*/
    int rreqcntFx=0; /*keeps track of forbidden color pointer communication*/
    int sncnt = 0; /*superstep number counter */
    int fp = 0; /*index to forbidProc*/
    int ierr = ZOLTAN_OK;

    /* DETEMINE THE SUPERSTEP NUMBER OF LOCAL VERTICES AND COMMUNICATE THEM */
    /* Issue async recvs for superstep info */
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
        if (p!=zz->Proc) 
            MPI_Irecv(srecbuf[p], srecsize[p], MPI_INT, p, SNTAG, zz->Communicator, &rreqsC[rreqcntC++]);


    /* Calculate superstep info to be sent */
    /* Initialize pointers for filling the buffers */
    for (p=0; p<zz->Num_Proc; p++)
        if (p != zz->Proc)
            ssp[p] = ssendbuf[p];
    /* Determine which local vertex will be colored in which superstep
       and put this info to the buffer to be sent to processors holding
       the neighbors of the vertex */
    n = 1;
    for (i=0; i<nvtx; ++i) {
        int u = visit[i], gu;
        gu = Zoltan_G2LHash_L2G(hash, u);
        for (j=xadj[u]; j<xbadj[u]; j++) {
            if (ssp[adjproc[j]] != ssendbuf[adjproc[j]]) {
                if (*(ssp[adjproc[j]]-1) != gu) {
                    *(ssp[adjproc[j]]++) = gu;
                }
            } else 
                *(ssp[adjproc[j]]++) = gu;
        }
        if (n == ss) { /* mark the end of superstep */
            n = 0;
            for (p=0; p<zz->Num_Proc; p++)
                if (p!= zz->Proc)
                    *(ssp[p]++) = -1; 
        }        
        ++n;
    }
    /* Mark end of superstep info by -2 */
    for (p=0; p<zz->Num_Proc; p++) {
        if (p != zz->Proc) {
            if (ssp[p] == ssendbuf[p]) 
                *ssp[p] = -2;
            else if (*(ssp[p]-1) == -1) 
                *(--ssp[p]) = -2;
            else 
                *ssp[p] = -2;
        }
    }

    /* Issue async sends for superstep info */
    for (sreqcntC=p=0; p<zz->Num_Proc; ++p) 
        if (p != zz->Proc) 
            MPI_Isend(ssendbuf[p], ssp[p] - ssendbuf[p] + 1, MPI_INT, p, SNTAG, zz->Communicator, &sreqsC[sreqcntC++]);
    

    MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all superstep info to be received */
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all superstep info to be sent */ 
    
    /***** COLORING *****/
    /* Issue async recvs for updated d1 colors */
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            rreqfromC[rreqcntC] = p;
            MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, COLORTAG, zz->Communicator, &rreqsC[rreqcntC++]);
        }

    /* Initialize superstep info buffer pointers */
    for (p=0; p<zz->Num_Proc; p++)
        if (p != zz->Proc) {
            ssp[p] = ssendbuf[p];
            srp[p] = srecbuf[p];
        }

    /* Initial forbidden color communication and confChk marking for the first ss ****/
    for (rreqcntFx=p=0; p<zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            sendNextStepsForbiddenColors(zz, hash, nlvtx, p, srp, xadj, adj, xadjnl, adjnl, forbsizeS, xforbiddenS, forbiddenS, nColor, color, mark, confChk, sncnt, wset, wsize, sreqsFx, sreqsF);
            rreqfromFx[rreqcntFx] = p;
            MPI_Irecv(xforbidden[p], ss + 1, MPI_INT, p, XFORBIDTAG, zz->Communicator, &rreqsFx[rreqcntFx++]);
        } else
            xforbidden[p][0] = 0;
    waitPtrAndForbiddenColors(zz, rreqcntFx, rreqsFx, rreqfromFx, stats, forbsize, xforbidden, forbidden, xfp, rreqsF);                

    /* Coloring */
    n = 0; 
    for (i = 0; i < nvtx; ++i) {
        int u = visit[i];
        if (confChk[u] == sncnt) { /* add boundary vertices of this round into confChk */
            confChk[u] = -2;
            wset[(*wsize)++] = u;
        }
        for (j = xadj[u]; j < xadj[u+1]; ++j) {
            int v =  adj[j], c;            
            int vp = adjproc[j]; /* owner of vertex v */
            if ((c = color[v]) != 0) {
                while (*nColor < c) /* if nColor < c, set nColor to c and reset intermediate mark entries. May be done better. */
                    mark[++(*nColor)] = -1;
                mark[c] = u;
            }
            if (vp == zz->Proc) { /* if v is local as well */
                if (isbound[v] && confChk[v] == sncnt) { /* add boundary vertices of this round into confChk */
                    confChk[v] = -2;
                    wset[(*wsize)++] = v;
                }
                for (k = xadj[v]; k < xadj[v+1]; ++k) {	  
                    int w = adj[k];
                    if ((c = color[w]) != 0) {
                        while (*nColor < c) /* if nColor < c, set nColor to c and reset intermediate mark entries. May be done better. */
                            mark[++(*nColor)] = -1;
                        mark[c] = u;
                    }
                }
            }
            else if (xfpMark[vp]!=u) { /* if v is non-local, read forbidden colors sent by the owner of v */                
                                                                            /* but do it only once for vp */
                for (m = *xfp[vp]; m < *(xfp[vp] + 1); m++) {
                    c = forbidden[vp][m];
                    while (*nColor < c) /* if nColor < c, set nColor to c and reset intermediate mark entries. May be done better. */
                        mark[++(*nColor)] = -1;
                    mark[c] = u;
                }
                xfpMark[vp] = u;
                ++(xfp[vp]);
            }
        }
        color[u] = PickColor(zz, color_method, u, color[u], nColor, mark);
        colored[n++] = Zoltan_G2LHash_L2G(hash, u);
        colored[n++] = color[u];
        
        /* If superstep is finished, communicate the updated colors and forbidden colors */
        if (n >= 2*ss) {
            ++sncnt;
            /* Issue new async recvs for forbidden color pointers */
            for (rreqcntFx=p=0; p<zz->Num_Proc; ++p)
                if (p!=zz->Proc) {
                    rreqfromFx[rreqcntFx] = p;
                    MPI_Irecv(xforbidden[p], ss + 1, MPI_INT, p, XFORBIDTAG, zz->Communicator, &rreqsFx[rreqcntFx++]);
                }                
            /* send local color updates */
            for (sreqcntC=p=0; p<zz->Num_Proc; ++p) 
                if (p!=zz->Proc) 
                    MPI_Isend(colored, 2*ss, MPI_INT, p, COLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);
            
            /* process the received color update messages */
            MPI_Waitall(rreqcntC, rreqsC, stats);
            fp = 0; /*number of procesors requesting forbidden colors for their next ss*/
            for (l = 0; l < rreqcntC; ++l) {
                p = rreqfromC[l];	
                /* update non-local vertex colors according to the received color update message content */  
                for (j = 0; j < 2*ss; ) { 
                    int v = newcolored[p][j++], c, hv;                   
                    if (v < 0)
                        break;
                    if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
                        c = newcolored[p][j++];
                        color[hv] = c;
                    } else
                        ++j;
                }
                /* if round is not finished on proc p */
                if (j >= 2*ss) 
                    rreqfromC[fp++] = p;                    
            }
            
            /* send forbidden colors to requesting processors */
	    rreqcntC = fp;
            for (l = 0; l < rreqcntC; l++) { 
                p = rreqfromC[l];
                sendNextStepsForbiddenColors(zz, hash, nlvtx, p, srp, xadj, adj, xadjnl, adjnl, forbsizeS, xforbiddenS, forbiddenS, nColor, color, mark, confChk, sncnt, wset, wsize, sreqsFx, sreqsF);
                /* Issue new async receive for color update */
                MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, COLORTAG, zz->Communicator, &rreqsC[l]);
            }            
            /* read forbidden color pointer messages and issue corresponding async recvs for forbidden colors*/
            waitPtrAndForbiddenColors(zz, rreqcntFx, rreqsFx, rreqfromFx, stats, forbsize, xforbidden, forbidden, xfp, rreqsF);                
            n = 0;
            MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all color updates to be sent */            
        }
    }
    
    /* round is finished on zz->Proc, make the last send */
    colored[n++] = -1;
    colored[n++] = -1;
    for (sreqcntC=p=0; p<zz->Num_Proc; ++p)
        if (p!=zz->Proc) 
            MPI_Isend(colored, n, MPI_INT, p, COLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);

    /* process the remaining color update messages */
    while (rreqcntC) { /*until all color update messages in this superstep are received*/
        MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all color updates to be received */
        fp = 0;
        ++sncnt;
        for (l = 0; l < rreqcntC; ++l) {            
            p = rreqfromC[l];	
            /* update non-local vertex colors according to the received color update message content */  
            for (j = 0; j < 2*ss; ) { 
                int v = newcolored[p][j++], c, hv;                   
                if (v < 0)
                    break;
                if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
                    c = newcolored[p][j++];
                    color[hv] = c;
                } else
                    ++j;
            }
	    /* if round is not finished on proc p */
	    if (j >= 2*ss)
                rreqfromC[fp++] = p;
	}
	/* send forbidden colors to requesting processors */
	rreqcntC = fp;
	for (l = 0; l < rreqcntC; l++) {
            p = rreqfromC[l];
            sendNextStepsForbiddenColors(zz, hash, nlvtx, p, srp, xadj, adj, xadjnl, adjnl, forbsizeS, xforbiddenS, forbiddenS, nColor, color, mark, confChk, sncnt, wset, wsize, sreqsFx, sreqsF);
            /* Issue new async receive for color update */
            MPI_Irecv(newcolored[p], 2*ss, MPI_INT, p, COLORTAG, zz->Communicator, &rreqsC[l]);
	}
    }
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait last color update to be sent */            

    return ierr;
}

/*****************************************************************************/
/* Detect D1 conflicts */

static int D1DetectConflicts(
    ZZ *zz,
    int nConflict,
    int *visit,
    int *xadj,
    int *xbadj,
    int *adj,
    int *color,
    int *conflicts,
    int *rand_key,
    G2LHash *hash
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
            } else if ((color[u] == color[v]) && (rand_key[u] == rand_key[v]) && (Zoltan_G2LHash_L2G(hash, u) < Zoltan_G2LHash_L2G(hash, v))){
                conflicts[conflict++] = u;
		break;
            }
        }
    }
    return conflict;
}

/*****************************************************************************/
/* Detect D2 conflicts */

static int D2DetectConflicts(ZZ *zz, G2LHash *hash, int nlvtx, int *wset, int wsize, int *xadj, int *adj, int *adjproc, int *nColor, int *color, int *conflicts, int *rand_key, int *vmark, int *seen, int *where, int *pwhere, int **rcsendbuf, int **rcrecbuf, int *rcsendsize, int *rcrecsize, int **srp, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromC, MPI_Status *stats, int *nconflict)
{
    static char *yo="D2DetectConflicts";
    int w, i, j, p;
    int rreqcntC=0, sreqcntC=0;
    int *q;
    int ierr = ZOLTAN_OK;

    for (i=0; i<zz->Num_Proc; ++i)
        if (i != zz->Proc)
            srp[i] = rcsendbuf[i];
        else
            srp[i] = conflicts;

    /*** Issue async recv for recolor info ***/
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
        if (p != zz->Proc) {
            rreqfromC[rreqcntC] = p;
            MPI_Irecv(rcrecbuf[p], rcrecsize[p], MPI_INT, p, RECOLORTAG, zz->Communicator, &rreqsC[rreqcntC++]);
        }

    /*** Detect conflicts and mark vertices to be recolored ***/
    memset(seen, 0xff, sizeof(int) * (*nColor+1));
    memset(where, 0xff, sizeof(int) * (*nColor+1));
    memset(pwhere, 0xff, sizeof(int) * (*nColor+1));

    for (i = 0; i < wsize; i++) {
        int x, v, cw, cx, px, pv, gv, gx;
        w = wset[i];
        cw = color[w];
        seen[cw] = w;
        where[cw] = w;
        pwhere[cw] = zz->Proc;
        for (j = xadj[w]; j < xadj[w+1]; j++) {
            x = adj[j];
            px = adjproc[j];
            cx = color[x];
            gx = Zoltan_G2LHash_L2G(hash, x);
            if (seen[cx] == w) {
                v = where[cx];
                pv = pwhere[cx];
                gv = Zoltan_G2LHash_L2G(hash, v);
                if (rand_key[v] <= rand_key[x]) {
                    if (!vmark[x]) {
                        *srp[px]++ = gx;
                        vmark[x] = 1;
                    }
                } else {
                    if (!vmark[v]) {
                        *srp[pv]++ = gv;
                        vmark[v] = 1;
                    }
                    where[cx] = x;
                    pwhere[cx] = px;
                }
            } else {
                seen [cx] = w;
                where[cx] = x;
                pwhere [cx] = px;
            }
        }
    }
    /*** Construct recolor messages ***/
    for (p = 0; p < zz->Num_Proc; p++) 
        if (p != zz->Proc) {
            *srp[p]++ = -1;
            if ((srp[p] - rcsendbuf[p]) > rcsendsize[p])
                ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Insufficient memory allocation, use larger sized buffers.");
            MPI_Isend(rcsendbuf[p], srp[p] - rcsendbuf[p] , MPI_INT, p, RECOLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);
            for (q = rcsendbuf[p]; *q!=-1; ++q) {
                int u = Zoltan_G2LHash_G2L(hash, *q);
                vmark[u] = 0;
            }
	    /* *q = 0; UVC: BUGBUG CHECK */
        }

    /*** Build the local vertices to be recolored ***/                    
    MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all recolor info to be received */
    for (p = 0; p < zz->Num_Proc; p++)
        if (p != zz->Proc) {
            j = 0;
            while (rcrecbuf[p][j] != -1) {
                int gu = rcrecbuf[p][j++], u;
                u = Zoltan_G2LHash_G2L(hash, gu);
                if (!vmark[u]) {
                    *srp[zz->Proc]++ = gu;
                    vmark[u] = 1;
                }
            }
        }
    *srp[zz->Proc] = -1;
    *nconflict = srp[zz->Proc] - conflicts;
    for (q = conflicts; *q != -1; ++q) {
        *q = Zoltan_G2LHash_G2L(hash, *q);
        vmark[*q] = 0;
    }
    /* *q = 0; UVC: this is OK to keep but it is useless; therefore it is commented out*/	
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all recolor info to be sent */
    
 End:
    
    return ierr;
}



