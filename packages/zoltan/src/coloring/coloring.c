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

#include <assert.h>
#include <ctype.h>
#include "zoltan_mem.h"
#include "zz_const.h"
#include "coloring.h"
#include "g2l_hash.h"
#include "params_const.h"
#include "zz_util_const.h"
#include "graph.h"
#include "all_allo_const.h"
#include "zz_rand.h"
#include "bucket.h"

/* when sending new colored vertices to processors,
   sent only the "relevant" ones; i.e., send the color info if the processors has
   D1 neighbour */

#define COLORTAG     1001
#define SNTAG        1002
#define XFORBIDTAG   1003
#define FORBIDTAG    1004
#define RECOLORTAG   1005
#define FORWARD        11
#define REVERSE        12
#define NONDECREASING  13
#define NONINCREASING  14
#define SYNCHRONOUS    21
#define ASYNCHRONOUS   22

/* Function prototypes */
static int D1coloring(ZZ *zz, char coloring_problem, char coloring_order, char coloring_method, char comm_pattern, int ss,
		      int nVtx, G2LHash *hash, int *xadj, int *adj, int *adjproc, int *color, int recoloring_permutation,
                      int recoloring_type, int recoloring_num_of_iterations);
static int D2coloring(ZZ *zz, char coloring_problem, char coloring_order, char coloring_method, char comm_pattern, int ss,
		      int nVtx, G2LHash *hash, int *xadj, int *adj, int *adjproc, int *color, int *partialD2);

static int ReorderGraph(ZZ *, char, int, int *, int **, int *,
			int *, int *, int *, int *, int *partialD2, int *nintvisit, int *nboundvisit);
static int PickColor(ZZ *, char, int, int, int *, int *);
static int InternalColoring(ZZ *zz, char coloring_problem, int *nColor,
			    int nVtx, int *visit, int * xadj, int *adj,
			    int *color, int *mark, int mark_size, char coloring_method);

static int D1ParallelColoring (ZZ *zz, int nvtx, int *visit, int *xadj, int *adj,
			       int *isbound, int ss, int *nColor, int *color,
			       ZOLTAN_GNO_TYPE **newcolored, int *mark, int gmaxdeg, G2LHash *hash,
			       char coloring_method, char comm_pattern, int *rreqfrom,
			       int *replies, MPI_Request *sreqs, MPI_Request *rreqs,
			       MPI_Status *stats,
			       int *xadjproc, int *adjproc, ZOLTAN_GNO_TYPE **persSbuf, int *Ssize, int plstcnt, int *plst);

static int sendNextStepsForbiddenColors(ZZ *zz, G2LHash *hash, int nlvtx, int p, ZOLTAN_GNO_TYPE **srp, int *xadj, int *adj, int *xadjnl, int *adjnl, int *forbsizeS, int **xforbiddenS, int **forbiddenS, int *nColor, int *color, int *mark, int *confChk, int sncnt, int *wset, int *wsize, MPI_Request *sreqsFx, MPI_Request *sreqsF);
static int waitPtrAndForbiddenColors(ZZ* zz, int rreqcntFx, MPI_Request *rreqsFx, int *rreqfromFx, MPI_Status *stats, int *forbsize, int **xforbidden, int **forbidden, int **xfp, MPI_Request *rreqsF);
static int D2ParallelColoring (ZZ *zz, int nvtx, int nlvtx, int *visit, int *xadj, int *adj, int *xbadj, int *xadjnl, int *adjnl, int *adjproc, int *isbound, int ss, int *nColor, int *color, ZOLTAN_GNO_TYPE **newcolored, int *mark, int gmaxdeg, G2LHash *hash, int **forbidden, int **xforbidden, int **forbiddenS, int **xforbiddenS, int *forbsize, int *forbsizeS, char coloring_problem, char coloring_method, char comm_pattern, int *rreqfromC, int *repliesC, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromF, int *repliesF, MPI_Request *sreqsF, MPI_Request *rreqsF, int *rreqfromFx, int *repliesFx, MPI_Request *sreqsFx, MPI_Request *rreqsFx, MPI_Status *stats, int *confChk, int *ssendsize, int *srecsize, ZOLTAN_GNO_TYPE **ssendbuf, ZOLTAN_GNO_TYPE **srecbuf, ZOLTAN_GNO_TYPE  ** ssp,  ZOLTAN_GNO_TYPE **srp, int **xfp, int *xfpMark, int *wset, int *wsize);
static int D1DetectConflicts(ZZ *, int, int *, int *, int *, int *, int *, int *, int *, G2LHash *);
static int D2DetectConflicts(ZZ *zz, char coloring_problem, G2LHash *hash, int nlvtx, int *wset, int wsize, int *xadj, int *adj, int *adjproc, int *nColor, int *color, int *conflicts, int *rand_key, int *vmark, int *seen, int *where, int *pwhere, ZOLTAN_GNO_TYPE **rcsendbuf, ZOLTAN_GNO_TYPE  **rcrecbuf, int *rcsendsize, int *rcrecsize, ZOLTAN_GNO_TYPE **srp, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromC, MPI_Status *stats, int *nconflict);
static int Recoloring (ZZ *zz, int recoloring_permutation, int recoloring_type,
                       int recoloring_num_of_iterations, int nvtx, int carrierbufsize, int *nConflict,
                       int *confCont, int *nTotConflict, int *nRound,int lastlno, int globmaxnvtx,
                       int *nColor, int *color, int *visit,int *xadj, int *adj, int *xbadj,
                       int *conflicts, int *rand_key, G2LHash *hash, int *isbound, 
                       ZOLTAN_GNO_TYPE **newcolored, int *mark, int gmaxdeg,
                       char coloring_method, char comm_pattern, int *rreqfrom, int *replies,
                       MPI_Request *sreqs, MPI_Request *rreqs, MPI_Status *stats,
                       int *xrelproc, int *relproc, ZOLTAN_GNO_TYPE **persSbuf,
                       int *Ssize, int plstcnt, int *plst);
static int pairofintsup_ind_nd (const void *a, const void *b);
static int pairofintsup_ind_ni (const void *a, const void *b);

/*****************************************************************************/
/*  Parameters structure for Color method.  Used in  */
/*  Zoltan_Color_Set_Param and Zoltan_Color.         */
static PARAM_VARS Color_params[] = {
		  { "COLORING_PROBLEM", NULL, "STRING", 0 },
		  { "SUPERSTEP_SIZE",   NULL, "INT", 0},
		  { "COMM_PATTERN",     NULL, "CHAR", 0 },
		  { "VERTEX_VISIT_ORDER",   NULL, "CHAR", 0 },
		  { "COLORING_METHOD",  NULL, "CHAR", 0},
                  { "RECOLORING_TYPE", NULL, "STRING", 0},
                  { "RECOLORING_PERMUTATION", NULL, "STRING", 0},
                  { "RECOLORING_NUM_OF_ITERATIONS", NULL, "INT", 0},
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
static void PrintGraph(ZZ *zz, char *name, int base, int nvtx, int *xadj, ZOLTAN_GNO_TYPE *adj, int *adjproc)
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
	    fprintf(outf, "%ld (%d), ", adj[j], adjproc[j]);
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
    int num_gid_entries,      /* # of entries for a global id */
    int num_req_objs,         /* Input: number of requested objects (vtx) */
    ZOLTAN_ID_PTR req_objs,   /* Input: global ids of the objects for
                                 which the colors should be returned.
                                 (May be different from local objects.) */
			      /* The application must allocate enough space */
    int *color_exp            /* Output: Colors assigned to objects 
                                 given by req_objs (may be non-local) */
			      /* The application must allocate enough space */
)
{
  char coloring_problem;   /* Input: which coloring to perform;
			   currently only supports D1, D2 coloring and partial D2 */
  char coloring_problemStr[MAX_PARAM_STRING_LEN]; /* string version coloring problem name */
  int ss;               /* Superstep size: detemines how many vertices are
			   locally colored before the next color
			   information exchange */
  char comm_pattern;    /* (A) asynchronous (S) synchronous supersteps*/
  char coloring_order;     /* (I) interior vertices first
			   (B) boundary vertices first (U) interleaved = (N) natural,
                           (L) largest degree first, (S) smallest degree last */
  char coloring_method;    /* Coloring method. (F) First fit */
  char recoloring_permutationStr[MAX_PARAM_STRING_LEN];
  int recoloring_permutation; /* FORWARD, REVERSE, NONDECREASING, NONINCREASING, default is NONDECREASING */
  char recoloring_typeStr[MAX_PARAM_STRING_LEN];
  int recoloring_type; /* ASYNCHRONOUS, SYNCHRONOUS, default is SYNCHRONOUS */
  int recoloring_num_of_iterations; /* must be a positive number, default is 1 if not specified */
  static char *yo = "Zoltan_Color";
  ZOLTAN_GNO_TYPE *vtxdist=NULL, *adjncy=NULL;
  int *itmp, *xadj=NULL;
  int *adjproc=NULL;
  int nvtx;                         /* number of local vertices */
  ZOLTAN_GNO_TYPE gvtx;             /* number of global vertices */

  int *color=NULL;                  /* array to store colors of local and D1
				       neighbor vertices */
  int i, lno;
  int lastlno;                      /* total number of local and D1 neighbor vertices */
  G2LHash hash;                     /* hash to map global ids of local and D1 neighbor
				       vertices to consecutive local ids */
  int ierr = ZOLTAN_OK;
  int comm[2],gcomm[2];
  ZOLTAN_ID_PTR my_global_ids= NULL;  /* gids local to this proc */
  struct Zoltan_DD_Struct *dd_color;  /* DDirectory for colors */
  ZOLTAN_GNO_TYPE *requested_GNOs = NULL;  /* Return GNOs of 
                                              the requested GIDs.  */


#ifdef _DEBUG_TIMES  
  double times[6]={0.,0.,0.,0.,0.,0.}, rtimes[6]; /* Used for timing measurements */
  double gtimes[6]={0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
  char *timenames[6]= {"", "setup", "graph build", "renumber", "color", "clean up"};
#endif
  int *partialD2 = NULL;       /* binary array showing which vertices to be colored */
  ZOLTAN_GNO_TYPE *local_GNOs = NULL;
  ZG graph;

  memset (&graph, 0, sizeof(ZG));
  memset(&hash, 0 , sizeof(G2LHash)); /* To allow a correct free */

#ifdef _DEBUG_TIMES    
  MPI_Barrier(zz->Communicator);
  times[0] = Zoltan_Time(zz->Timer);
#endif

  /* PARAMETER SETTINGS */

  Zoltan_Bind_Param(Color_params, "COLORING_PROBLEM",   (void *) &coloring_problemStr);
  Zoltan_Bind_Param(Color_params, "SUPERSTEP_SIZE",     (void *) &ss);
  Zoltan_Bind_Param(Color_params, "COMM_PATTERN",       (void *) &comm_pattern);
  Zoltan_Bind_Param(Color_params, "VERTEX_VISIT_ORDER", (void *) &coloring_order);
  Zoltan_Bind_Param(Color_params, "COLORING_METHOD",    (void *) &coloring_method);
  Zoltan_Bind_Param(Color_params, "RECOLORING_PERMUTATION", (void *) &recoloring_permutationStr);
  Zoltan_Bind_Param(Color_params, "RECOLORING_TYPE",        (void *) &recoloring_typeStr);
  Zoltan_Bind_Param(Color_params, "RECOLORING_NUM_OF_ITERATIONS", (void *) &recoloring_num_of_iterations);

  /* Set default values */
  strncpy(coloring_problemStr, "distance-1", MAX_PARAM_STRING_LEN);
  coloring_problem = '1';
  ss = 100;
  coloring_method = 'F';
  coloring_order = 'I';
  comm_pattern = 'S';
  strncpy(recoloring_permutationStr, "NONDECREASING", MAX_PARAM_STRING_LEN);
  strncpy(recoloring_typeStr, "SYNCHRONOUS", MAX_PARAM_STRING_LEN);
  recoloring_permutation = NONDECREASING;
  recoloring_type = SYNCHRONOUS;
  recoloring_num_of_iterations = 0;

  Zoltan_Assign_Param_Vals(zz->Params, Color_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);

  /* Check validity of parameters */
  if (!strcasecmp(coloring_problemStr, "distance-1"))
      coloring_problem = '1';
  else if (!strcasecmp(coloring_problemStr, "distance-2"))
      coloring_problem = '2';
  else if (!strcasecmp(coloring_problemStr, "partial-distance-2")
      || !strcasecmp(coloring_problemStr, "bipartite"))
      coloring_problem = 'P';
  else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Unknown coloring requested. Using Distance-1 coloring.");
      coloring_problem = '1';
  }
  if (ss == 0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid superstep size. Using default value 100.");
      ss = 100;
  }
  if (comm_pattern != 'S' && comm_pattern != 'A') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid communication pattern. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (comm_pattern == 'A' && (coloring_problem == '2' || coloring_problem == 'P')) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Asynchronous communication pattern is not implemented for distance-2 coloring and its variants. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (coloring_order != 'I' && coloring_order != 'B' && coloring_order != 'U' && coloring_order != 'N' && coloring_order != 'L' && coloring_order != 'S') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring order. Using internal first coloring order (I).");
      coloring_order = 'I';
  }
  if (coloring_order == 'U' && (coloring_problem == '2' || coloring_problem == 'P')) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Interleaved coloring order is not implemented for distance-2 coloring and its variants. Using internal first coloring order (I).");
      coloring_order = 'I';
  }
  if (coloring_method !='F') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring method. Using first fit method (F).");
      coloring_method = 'F';
  }
  if (recoloring_num_of_iterations > 0) {
      if (!strcasecmp(recoloring_permutationStr, "FORWARD"))
          recoloring_permutation = FORWARD;
      else if (!strcasecmp(recoloring_permutationStr, "REVERSE"))
          recoloring_permutation = REVERSE;
      else if (!strcasecmp(recoloring_permutationStr, "NONDECREASING"))
          recoloring_permutation = NONDECREASING;
      else if (!strcasecmp(recoloring_permutationStr, "NONINCREASING"))
          recoloring_permutation = NONINCREASING;
      else {
          ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid recoloring permutation. Using NONDECREASING  recoloring permutation.");
          recoloring_permutation = NONDECREASING;
      }
      if (zz->Num_Proc > 1) {
          if (!strcasecmp(recoloring_typeStr, "SYNCHRONOUS"))
              recoloring_type = SYNCHRONOUS;
          else if (!strcasecmp(recoloring_typeStr, "ASYNCHRONOUS"))
              recoloring_type = ASYNCHRONOUS;
          else {
              ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid recoloring type. Using SYNCHRONOUS recoloring type ");
              recoloring_type = SYNCHRONOUS;
          }
      }
  }
  else if (recoloring_num_of_iterations < 0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid recoloring number of iterations. Using 1 iteration for recoloring.");
      recoloring_num_of_iterations = 0;
  }


  /* Compute Max number of array entries per ID over all processors.
     This is a sanity-maintaining step; we don't want different
     processors to have different values for these numbers. */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = gcomm[0];
  zz->Num_LID = gcomm[1];

  if (num_gid_entries != zz->Num_GID)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "num_gid_entries is not consistent with the queries.");

  /* Return if this processor is not in the Zoltan structure's
     communicator. */
  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz))
      return ZOLTAN_OK;

  /* BUILD THE GRAPH */
  /* TODO: Allow req_objs==NULL as special case for local vertices. */
  /* Check that the user has allocated space for the return args. */
  if (num_req_objs && !color_exp)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Output argument color_exp is NULL. Please allocate all required arrays before calling this routine.");

#ifdef _DEBUG_TIMES
  times[1] = Zoltan_Time(zz->Timer);
#endif

  requested_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(num_req_objs *
                                                       sizeof(ZOLTAN_GNO_TYPE));

  ierr =  Zoltan_ZG_Build (zz, &graph, 0, 1, num_req_objs,
                           req_objs, requested_GNOs);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Cannot construct graph.");

  ierr = Zoltan_ZG_Export (zz, &graph,
		    &gvtx, &nvtx, NULL, NULL, 
                    &vtxdist, &xadj, &adjncy, &adjproc,
		     NULL, NULL);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Cannot construct graph (2).");
#ifdef _DEBUG_TIMES
  times[2] = Zoltan_Time(zz->Timer);
#endif

  if (nvtx && !(local_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nvtx * sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;
  for (i=0; i<nvtx; ++i)
      local_GNOs[i] = i+vtxdist[zz->Proc];


  /* CREATE THE HASH TABLE */
  /* Determine hash size and allocate hash table */
  i = xadj[nvtx]; /* i is the minimum hash size */
  if (Zoltan_G2LHash_Create(&hash, i, vtxdist[zz->Proc], nvtx)==ZOLTAN_MEMERR)
      MEMORY_ERROR;

  /* Add global ids of the d1 neighbors into the hash table,
   *    create a "local ID" for each neighbor if it's not mine */

  if (sizeof(ZOLTAN_GNO_TYPE) ==          /* size of global id */
      sizeof(int)) {                       /* size of local id */

    for (i=0; i<xadj[nvtx]; ++i)
        adjncy[i] = (ZOLTAN_GNO_TYPE)Zoltan_G2LHash_Insert(&hash, adjncy[i]);
  }
  else {
    itmp = (int *)adjncy;
    for (i=0; i<xadj[nvtx]; ++i){
        lno = Zoltan_G2LHash_Insert(&hash, adjncy[i]);
        itmp[i] = lno;
    }
  }

  /* lastlno is the total number of local and d1 neighbors */
  lastlno = nvtx+hash.size;

  /* Allocate color array. Local and D1 neighbor colors will be stored here. */
  if (lastlno && !(color = (int *) ZOLTAN_CALLOC(lastlno, sizeof(int))))
      MEMORY_ERROR;

  if (coloring_problem == 'P') {
      int sz=(nvtx>num_req_objs) ? nvtx : num_req_objs;

      if (sz && !(partialD2 = (int *) ZOLTAN_CALLOC(sz, sizeof(int))))
	  MEMORY_ERROR;

      ierr = Zoltan_DD_Create (&dd_color, zz->Communicator, 
                               sizeof(ZOLTAN_GNO_TYPE)/sizeof(ZOLTAN_ID_TYPE), 0, 0, 0, 0);
      if (ierr != ZOLTAN_OK)
          ZOLTAN_COLOR_ERROR(ierr, "Cannot construct DDirectory.");
      /* Put req obs with 1 but first inialize the rest with 0 */
      ierr = Zoltan_DD_Update (dd_color, (ZOLTAN_ID_TYPE *)local_GNOs, NULL,
                               NULL, partialD2, nvtx);
      if (ierr != ZOLTAN_OK)
          ZOLTAN_COLOR_ERROR(ierr, "Cannot update DDirectory.");
      for (i=0; i<num_req_objs; ++i)
          partialD2[i] = 1;
      ierr = Zoltan_DD_Update (dd_color, (ZOLTAN_ID_TYPE *)requested_GNOs, NULL,
                               NULL, partialD2, num_req_objs);
      if (ierr != ZOLTAN_OK)
          ZOLTAN_COLOR_ERROR(ierr, "Cannot update DDirectory.");
      /* Get requested colors from the DDirectory. */
      ierr = Zoltan_DD_Find (dd_color, (ZOLTAN_ID_TYPE *)local_GNOs, NULL, NULL,
                             partialD2, nvtx, NULL);
      if (ierr != ZOLTAN_OK)
          ZOLTAN_COLOR_ERROR(ierr, "Cannot find object in DDirectory.");
      /* Free DDirectory */
      Zoltan_DD_Destroy(&dd_color);
  }

#ifdef _DEBUG_TIMES
  times[3] = Zoltan_Time(zz->Timer);
#endif
  /* Select Coloring algorithm and perform the coloring */
  if (coloring_problem == '1')
      D1coloring(zz, coloring_problem, coloring_order, coloring_method, comm_pattern, ss, nvtx, &hash, xadj, (int *)adjncy, adjproc, color,
		 recoloring_permutation, recoloring_type, recoloring_num_of_iterations);
  else if (coloring_problem == '2' || coloring_problem == 'P')
      D2coloring(zz, coloring_problem, coloring_order, coloring_method, comm_pattern, ss, nvtx, &hash, xadj, (int *)adjncy, adjproc, color, partialD2);
#ifdef _DEBUG_TIMES    
  times[4] = Zoltan_Time(zz->Timer);
#endif

   ierr = Zoltan_DD_Create (&dd_color, zz->Communicator, 
                            sizeof(ZOLTAN_GNO_TYPE)/sizeof(ZOLTAN_ID_TYPE), 0, 0, 0, 0);
   if (ierr != ZOLTAN_OK)
       ZOLTAN_COLOR_ERROR(ierr, "Cannot construct DDirectory.");
   /* Put color in part field. */
   ierr = Zoltan_DD_Update (dd_color, (ZOLTAN_ID_TYPE *)local_GNOs, NULL,
                            NULL, color, nvtx);
   if (ierr != ZOLTAN_OK)
       ZOLTAN_COLOR_ERROR(ierr, "Cannot update DDirectory.");

   /* Get requested colors from the DDirectory. */
   ierr = Zoltan_DD_Find (dd_color, (ZOLTAN_ID_TYPE *)requested_GNOs, NULL, NULL,
            color_exp, num_req_objs, NULL);
   if (ierr != ZOLTAN_OK)
     ZOLTAN_COLOR_ERROR(ierr, "Cannot find object in DDirectory.");

   /* Free DDirectory */
   Zoltan_DD_Destroy(&dd_color);
   ZOLTAN_FREE(&my_global_ids); 

#ifdef _DEBUG_TIMES    
  MPI_Barrier(zz->Communicator);
  times[5] = Zoltan_Time(zz->Timer);
  rtimes[0] = times[5]-times[0];
  for (i=1; i<6; ++i) 
      rtimes[i] = times[i]-times[i-1]; 
  MPI_Reduce(rtimes, gtimes, 6, MPI_DOUBLE, MPI_MAX, 0, zz->Communicator);
  if (!zz->Proc) {
      int i;
      for (i=1; i<6; ++i)
          printf("Zoltan_Color %-13s in Proc-0: %8.2lf  Max: %8.2lf\n", timenames[i], times[i]-times[i-1], gtimes[i]);
      printf("Zoltan_Color %-13s in Proc-0: %8.2lf  Max: %8.2lf\n", "Total Time", times[5]-times[0], gtimes[0]);
  }
#endif
 End:
  /* First, free graph */
  Zoltan_ZG_Free (&graph);
  ZOLTAN_FREE(&requested_GNOs);
  ZOLTAN_FREE(&adjproc);
  ZOLTAN_FREE(&color);
  ZOLTAN_FREE(&partialD2);
  ZOLTAN_FREE(&local_GNOs);
  Zoltan_G2LHash_Destroy(&hash);

  return ierr;
}


/* fills the visit array with the n first vertices of xadj using the
   Largest Degree First ordering. The algorithm used to compute this
   ordering is a stable count sort. */
static void LargestDegreeFirstOrdering(
    ZZ  *zz, 
    int *visit, /*Out*/
    int *xadj,
    int n,
    int max_degree) 
{	
    static char* yo = "LargestDegreeFirstOrdering";
    int ierr = ZOLTAN_OK;
    int i;
    int *cnt;
    
    cnt = (int*) ZOLTAN_MALLOC (sizeof(int)*(max_degree+1));
    if (!cnt)
        MEMORY_ERROR;
    
    memset(cnt, 0, sizeof(int)*(max_degree+1));
    for (i=0; i<n; ++i) {
        int degree_of_i = xadj[i+1]-xadj[i];
        ++cnt[degree_of_i];
    }
    /* now, cnt[i] is the number of vertice which degree is i*/

    /* this loop performs a prefix sum array computation on cnt*/
    for (i=1; i<max_degree+1; ++i) 
        cnt[i] = cnt[i] + cnt[i-1];

  
    /* cnt[x]-1 is the index of the next node of degree x*/
    for (i=0; i<n; ++i) {
        int degree_of_i = xadj[i+1]-xadj[i];
        cnt[degree_of_i]--;
        visit[n - 1 - cnt[degree_of_i]] = i;
    }
    
End:
    ZOLTAN_FREE(&cnt);
}


static void SmallestDegreeLastOrdering(
    ZZ  *zz, 
    int *visit, /*Out*/
    int *xadj,
    int *adj,
    int n,
    int max_degree) 
{	
  static char* yo = "SmallestDegreeLastOrdering";
  int ierr = ZOLTAN_OK;
  int i;
  Bucket bs;
  
  bs = Zoltan_Bucket_Initialize(max_degree+1,n);
  if (bs.buckets == NULL) 
      MEMORY_ERROR;

  for (i=0; i<n; i++) {
      int degree_of_i = xadj[i+1] - xadj[i];
      Zoltan_Bucket_Insert(&bs, i, degree_of_i);
  }

  for (i=0; i<n; i++) {
      int u = Zoltan_Bucket_PopMin(&bs);
      int j;
      
      visit[n - 1 - i] = u;
      for (j = xadj[u]; j < xadj[u+1]; ++j) { /* decrease the degree of the neighboors*/
	  if (adj[j] < n)  /* Otherwise it is not a local vertex */
              Zoltan_Bucket_DecVal(&bs, adj[j]);
      }    
  }
End:
  Zoltan_Bucket_Free(&bs);		
}




/*****************************************************************************/
/* Distance-1 coloring. No two adjacent vertices get the same color. */

static int D1coloring(
    ZZ *zz,
    char coloring_problem,/* Coloring problem. '1' in this case */
    char coloring_order,  /* (I) interior vertices first (B) boundary
			  vertices first (U) interleaved, (N) Natural
			  Ordering, (L) Largest Degree First, (S)
			  Smallest Degree Last (dynamic) */
    char coloring_method, /* Coloring method. (F) First fit
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
    int *color,         /* return array to store colors of local and D1
			  neighbor vertices */
    int recoloring_permutation, /* recoloring permutation type; FORWARD, 
                                     REVERSE, NONDECREASING, NONINCREASING */
    int recoloring_type, /* recoloring type; SYNCHRONOUS or ASYNCHRONOUS */
    int recoloring_num_of_iterations
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
    int *visitIntern = NULL;
    int nintvisit, nboundvisit;  /* NOTUSED: Number of internal and boundary vertices to be colored */
    int *mark = NULL;            /* Array to mark forbidden colors for a vertex */
    int gmaxdeg = 0;             /* Maximum vertex degree in the graph */
    int lmaxdeg = 0;             /* Maximum vertex degree for the local vertices */
    int lastlno = nvtx+hash->size; /* Total number of local and D1 neighbor vertices */
    int *conflicts = NULL;       /* List of vertices to be recolored in the next
				    round */
    int *replies = NULL;         /* Arrays used for MPI communication */
    MPI_Status *stats = NULL;
    MPI_Request *rreqs = NULL;
    MPI_Request *sreqs = NULL;
    int *rreqfrom = NULL;
    ZOLTAN_GNO_TYPE **newcolored = NULL;     /* Array used for communicating boundary vertex
				    colors at the end of supersteps. Global number
				    of the vertex - color of vertex pairs are filled
				    in this array */
    double times[6]={0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
    int get_times;               /* (1) Measure timings (0) Don't */
    int ierr;
    int j, k;
    int *pmark=NULL;
    ZOLTAN_GNO_TYPE **persSbuf=NULL; /* personalized send buffers */
    int *Ssize=NULL; /* send buffer sizes */
    int *relproc=NULL;
    int *xrelproc=NULL;
    int *plst=NULL, plstcnt=0;
    int globmaxnvtx = 0; /* maximum number of vertices in all the procs */
    int carrierbufsize = ss; /* size of the buffer to transfer data between procs,
                                it is either superstep size or globmaxnvtx */

    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    if (nvtx && (!isbound || !visit))
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
	printf(ZOLTAN_GNO_SPEC " [%d] :: ", Zoltan_G2LHash_L2G(hash, i), i);
	for (j=xadj[i]; j < xadj[i+1]; j++) {
	    printf(ZOLTAN_GNO_SPEC " [%d] (%d) ",   Zoltan_G2LHash_L2G(hash, adj[j]), adj[j], adjproc[j]);
	}
	printf("\n");
    }
#endif
    
    ierr = ReorderGraph(zz, coloring_problem, nvtx, xadj, &xbadj, adj, adjproc, &nbound, isbound, visit, NULL, &nintvisit, &nboundvisit);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
	ZOLTAN_COLOR_ERROR(ierr, "Error in ReorderGraph");
    if (get_times) times[1] = Zoltan_Time(zz->Timer);
    
#if 0
    printf("After reordering: nvtx:%d lastlno=%d Proc:%d\n", nvtx, lastlno, zz->Proc);
    for (i=0; i < nvtx; i++) {
	int j;
	printf((ZOLTAN_GNO_SPEC " [%d] :: ", Zoltan_G2LHash_L2G(hash, i), i);
	for (j=xadj[i]; j < xadj[i+1]; j++) {
	    printf((ZOLTAN_GNO_SPEC " [%d] (%d) ",   Zoltan_G2LHash_L2G(hash, adj[j]), adj[j], adjproc[j]);
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
    /* gmaxdeg+1 is the upper bound for #colors and colors start at one */
    gmaxdeg += 2;
    /* if recoloring is enabled, maximum number of vertices is computed and
       carrierbufsize is set accordingly for allocation purposes */
    if (recoloring_num_of_iterations > 0) {
        MPI_Allreduce(&nvtx, &globmaxnvtx, 1, MPI_INT, MPI_MAX, zz->Communicator);
        carrierbufsize = globmaxnvtx;
    }

    /* Memory allocation */
    mark = (int *) ZOLTAN_CALLOC(gmaxdeg, sizeof(int));
    conflicts = (int *) ZOLTAN_MALLOC((1+nvtx) * sizeof(int));
    replies = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    stats = (MPI_Status *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Status));
    rreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    sreqs = (MPI_Request *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(MPI_Request));
    rreqfrom = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    newcolored = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    if (!mark || !conflicts || !replies || !stats || !rreqs || !sreqs || !rreqfrom || !newcolored)
	MEMORY_ERROR;
    for (i=0; i<zz->Num_Proc; ++i) {
	newcolored[i] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(2 * carrierbufsize * sizeof(ZOLTAN_GNO_TYPE));
	if (!newcolored[i])
	    MEMORY_ERROR;
	memset(newcolored[i], 0, 2 * carrierbufsize * sizeof(ZOLTAN_GNO_TYPE));
    }
    persSbuf = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE *) * zz->Num_Proc);
    Ssize = (int *) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
    plst = (int *) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
    for (i=0; i<zz->Num_Proc; i++)
	persSbuf[i] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * 2 * carrierbufsize);

    /* Generate random numbers associated with global numbers of the vertices */
    /* All processors generate the same random number corresponding
       to the same global vertex number */
    rand_key = (int *) ZOLTAN_MALLOC(sizeof(int) * lastlno);
    if (lastlno && !rand_key)
	MEMORY_ERROR;
    for (i=0; i<lastlno; i++) {

	/* TODO64 OK if 8byte int is cast as 4byte int here? */
	Zoltan_Srand((unsigned int)Zoltan_G2LHash_L2G(hash, i), NULL);

	rand_key[i] = (int) (((double)Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
    }

    /* Color internal vertices and determine the visit order */
    if (get_times) {
	MPI_Barrier(zz->Communicator);
	times[2] = Zoltan_Time(zz->Timer);
    }
    visitIntern = visit + nbound; /* Start of internal vertex visit order. Used with I and B options below */
    if (coloring_order == 'U' || coloring_order == 'N' || coloring_order == 'L' || coloring_order == 'S') {
	nConflict = nvtx;

	/* change the order of the vertices */
	/*The ReorderGraph function also set the visit array. But I believe it is only useful for 'B' and 'I' */
	if (coloring_order == 'N' || coloring_order == 'U') /*natural order*/
	  {
	    for (i=0; i<nvtx; i++)
	      visit[i] = i;
	  }
	else if (coloring_order == 'L') /*largest first*/
	    LargestDegreeFirstOrdering(zz, visit, xadj, nvtx, lmaxdeg);
	else if (coloring_order == 'S') /*smallest last*/
	    SmallestDegreeLastOrdering(zz, visit, xadj, adj, nvtx, lmaxdeg);


	if (zz->Num_Proc==1)
	    InternalColoring(zz, coloring_problem, &nColor, nvtx, visit, xadj, adj, color, mark, gmaxdeg, coloring_method);
    }
    else if (coloring_order == 'I') {
	InternalColoring(zz, coloring_problem, &nColor, nvtx - nbound, visitIntern, xadj, adj, color, mark, gmaxdeg, coloring_method);
	nConflict = nbound;
    }
    else if (coloring_order == 'B')
	nConflict = nbound;


    if (get_times) times[3] = Zoltan_Time(zz->Timer);

    /* Color boundary vertices */
    nRound = 0;
    nTotConflict = 0;

    if (zz->Num_Proc >= 2) {

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
	printf("[%d] #of relevalnt procs is %d and procs are: ", zz->Proc, plstcnt);
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

	if (recoloring_num_of_iterations > 0)
	    carrierbufsize = ss;
	do {
	    int *tp = visit;

	    memset(mark, 0xff, (1+nColor) * sizeof(int)); /* reset dirty entries */
	    ierr = D1ParallelColoring(zz, nConflict, visit, xadj, adj, isbound, carrierbufsize,
				      &nColor, color, newcolored, mark, gmaxdeg, hash,
				      coloring_method, comm_pattern, rreqfrom, replies,
				      sreqs, rreqs, stats,
				      xrelproc, relproc, persSbuf, Ssize, plstcnt, plst);
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
    if (coloring_order == 'B')
	InternalColoring(zz, coloring_problem, &nColor, nvtx-nbound, visitIntern, xadj, adj, color, mark, gmaxdeg, coloring_method);

    /* Recoloring is run if enabled,i.e. num of iters is greater than 0 */
    if (recoloring_num_of_iterations > 0) {
        ierr = Recoloring(zz, recoloring_permutation, recoloring_type,
                          recoloring_num_of_iterations, nvtx, carrierbufsize, &nConflict,
                          &confCont, &nTotConflict, &nRound, lastlno, globmaxnvtx,
                          &nColor, color, visit, xadj, adj, xbadj,
                          conflicts, rand_key, hash, isbound,
                          newcolored, mark, gmaxdeg,
                          coloring_method, comm_pattern, rreqfrom, replies,
                          sreqs, rreqs, stats,
                          xrelproc, relproc, persSbuf,
                          Ssize, plstcnt, plst);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
            ZOLTAN_COLOR_ERROR(ierr, "Error in Recoloring");
    }

#if 0
    printf("[%d] vtx(gno)->color: ", zz->Proc);
    for (i=0; i<nvtx; i++)
	printf("%d(" ZOLTAN_GNO_SPEC ")->%d ", i, Zoltan_G2LHash_L2G(hash, i), color[i]);
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

    for (i=0; i<zz->Num_Proc; i++)
	ZOLTAN_FREE(&persSbuf[i]);
    ZOLTAN_FREE(&persSbuf);
    ZOLTAN_FREE(&Ssize);

    ZOLTAN_FREE(&plst);
    ZOLTAN_FREE(&pmark);
    ZOLTAN_FREE(&relproc);
    ZOLTAN_FREE(&xrelproc);
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
    char coloring_problem,/* Coloring problem. '2' or 'P'  in this case */
    char coloring_order,  /* (I) interior vertices first (B) boundary vertices first 
		             (U) interleaved (N) Natural Ordering (L) Largest First (S) Smallest Last
 		 	     Note that LF and SL orderings are same with the ones in distance-1 function.*/
    char coloring_method, /* Coloring method. (F) First fit
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
    int *color,        /* return array to store colors of local and D1
			  neighbor vertices */
    int *partialD2     /* binary array showing which vertices will be colored */
)
{
    static char *yo = "D2coloring";
    int i, j, p;
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
    int nintvisit, nboundvisit;  /* Number of internal and boundary vertices to be colored */
    int *mark = NULL;            /* Array to mark forbidden colors for a vertex */
    int gmaxcolor = 0;           /* Maximum #colors */
    int lmaxdeg = 0;             /* Maximum vertex degree for the local vertices */
    int lastlno = nvtx+hash->size; /* Total number of local and D1 neighbor vertices */
    int *conflicts = NULL;       /* List of vertices to be recolored in the next
				    round */
    /* Arrays used for MPI communication */
    int *repliesF = NULL, *repliesFx = NULL, *repliesC = NULL;
    MPI_Status *stats = NULL;
    MPI_Request *rreqsF = NULL, *rreqsFx = NULL, *rreqsC = NULL;
    MPI_Request *sreqsF = NULL, *sreqsFx = NULL, *sreqsC = NULL;
    int *rreqfromF = NULL, *rreqfromFx = NULL, *rreqfromC = NULL;
    int rreqcntC, sreqcntC;
    ZOLTAN_GNO_TYPE **newcolored = NULL;     /* Array used for communicating boundary vertex
				    colors at the end of supersteps. Global number
				    of the vertex - color of vertex pairs are filled
				    in this array */
    double times[6]={0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
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
    int *srecsize=NULL, *ssendsize=NULL;
    ZOLTAN_GNO_TYPE **srecbuf = NULL, **ssendbuf=NULL;
    ZOLTAN_GNO_TYPE **ssp, **srp;
    /* vertex recolor info exchange */
    ZOLTAN_GNO_TYPE **rcsendbuf = NULL, **rcrecbuf=NULL;
    int *rcsendsize=NULL, *rcrecsize=NULL;

    /* Memory allocation */
    isbound = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int));
    visit = (int *) ZOLTAN_MALLOC((1+nvtx) * sizeof(int));
    if ((nvtx && !isbound) || !visit)
	MEMORY_ERROR;

    /* Start timer */
    get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
    if (get_times){
	MPI_Barrier(zz->Communicator);
	times[0] = Zoltan_Time(zz->Timer);
    }
    ierr = ReorderGraph(zz, coloring_problem, nvtx, xadj, &xbadj, adj, adjproc, &nbound, isbound, visit, partialD2, &nintvisit, &nboundvisit);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
	ZOLTAN_COLOR_ERROR(ierr, "Error in ReorderGraph");
    if (get_times) times[1] = Zoltan_Time(zz->Timer);

#if 0
    printf("After reordering: nvtx:%d lastlno=%d Proc:%d\n", nvtx, lastlno, zz->Proc);
    for (i=0; i < nvtx; i++) {
	printf((ZOLTAN_GNO_SPEC " (%d) :: ", Zoltan_G2LHash_L2G(hash, i), i);
	for (j=xadj[i]; j < xadj[i+1]; j++) {
	    if (j == xbadj[i])
		printf(" | ");
	    printf((ZOLTAN_GNO_SPEC " (%d) ", Zoltan_G2LHash_L2G(hash, adj[j]), adj[j]);
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
    /* gmaxdeg^2+1 is the upper bound for #colors and colors start at 1*/
    gmaxcolor = gmaxcolor*gmaxcolor+2;


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
	    ssendsize[i] += (int) (ceil((double)nvtx / (double)ss) + 1);
          
    /* Send the superstep info size so that other processors will allocate enough space */
    for (sreqcntC=p=0; p<zz->Num_Proc; ++p)
	if (p != zz->Proc)
	    MPI_Isend(&ssendsize[p], 1, MPI_INT, p, SNTAG, zz->Communicator, &sreqsC[sreqcntC++]);

    mark = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    vmark = (int *) ZOLTAN_CALLOC(lastlno, sizeof(int));
    conflicts = (int *) ZOLTAN_MALLOC((1+nvtx) * sizeof(int));
    if ((gmaxcolor && !mark) || !conflicts || (lastlno && !vmark))
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
    newcolored = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    forbidden = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    xforbidden = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    forbiddenS = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    xforbiddenS = (int **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int *));
    if (!forbsize || !forbsizeS || !xfp || !xfpMark || !newcolored || !forbidden ||
	!xforbidden || !forbiddenS || !xforbiddenS)
	MEMORY_ERROR;
    for (i=0; i<zz->Num_Proc; ++i) {
	newcolored[i] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(2 * ss * sizeof(ZOLTAN_GNO_TYPE)); /* DB: can be reduced */
	xforbidden[i] = (int *) ZOLTAN_MALLOC((ss+1) * sizeof(int));
	forbsize[i] = (int)(ceil((double) xadj[nvtx] / (double)nvtx) * ss); /* size of forbidden color allocation - init to avgDeg*ss */
	forbidden[i] = (int *) ZOLTAN_MALLOC((forbsize[i] > 0 ? forbsize[i] : 1) * sizeof(int));
	xforbiddenS[i] = (int *) ZOLTAN_MALLOC((ss+1) * sizeof(int));
	forbsizeS[i] = (int)(ceil((double) xadj[nvtx] / (double) nvtx) * ss); /* size of forbidden color allocation for send buffer - init to avgDeg*ss */
	forbiddenS[i] = (int *) ZOLTAN_MALLOC((forbsizeS[i] > 0 ? forbsizeS[i] : 1) * sizeof(int));
	if ((ss && !newcolored[i]) || !xforbidden[i] || !forbidden[i] || !xforbiddenS[i] || !forbiddenS[i])
	    MEMORY_ERROR;
    }
    confChk = (int *) ZOLTAN_MALLOC(nvtx * sizeof(int)); /*the vertices to be checked in conflict detection */
    wset = (int *) ZOLTAN_MALLOC((nbound > 0 ? nbound : 1) * sizeof(int));
    seen = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    where = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    pwhere = (int *) ZOLTAN_MALLOC(gmaxcolor * sizeof(int));
    if ((nvtx && !confChk) || !wset || (gmaxcolor && (!seen || !where || !pwhere)))
	MEMORY_ERROR;

    /* Wait for superstep size communication to end */
    MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all superstep numbers to be received */
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all superstep numbers to be sent */

    /* Allocate buffers for superstep size and recolor info to be sent and received */
    ssendbuf = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    srecbuf = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    rcrecbuf = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    rcsendbuf = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    ssp = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    srp = (ZOLTAN_GNO_TYPE **) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(ZOLTAN_GNO_TYPE *));
    if (!ssendbuf || !srecbuf || !ssp || !srp || !rcsendbuf || !rcrecbuf)
	MEMORY_ERROR;
    for (p=0; p<zz->Num_Proc; p++){
	if (p != zz->Proc) {
	    ssendbuf[p] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((ssendsize[p] > 0 ? ssendsize[p] : 1) * sizeof(ZOLTAN_GNO_TYPE));
	    srecbuf[p] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((srecsize[p] > 0 ? srecsize[p] : 1) * sizeof(ZOLTAN_GNO_TYPE));
	    rcrecsize[p] = 2 * ssendsize[p]; /* recolor info exchange buffer sizes. */
	    rcsendsize[p] = 2 * srecsize[p]; /* better estimates or realloc in D2DetectConflicts can be used */
	    rcrecbuf[p] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((rcrecsize[p] > 0 ? rcrecsize[p] : 1) * sizeof(ZOLTAN_GNO_TYPE));
	    rcsendbuf[p] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((rcsendsize[p] > 0 ? rcsendsize[p] : 1) * sizeof(ZOLTAN_GNO_TYPE));
	    if (!srecbuf[p] || !ssendbuf[p] || !rcsendbuf[p] || !rcrecbuf[p])
		MEMORY_ERROR;
	} else {
	    ssendbuf[p] = NULL;
	    srecbuf[p] = NULL;
	    rcrecbuf[p] = NULL;
            /* We need this buffer in the case where sizeof(ZOLTAN_GNO_TYPE) > sizeof(int) because the
             *  the processes store global IDs for local vertices here in D2DetectConflicts.
             */
	    rcsendbuf[p] = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((lastlno > 0 ? lastlno : 1) * sizeof(ZOLTAN_GNO_TYPE));
	    if (!rcsendbuf[p])
		MEMORY_ERROR;
	}
    }

    /* Generate random numbers associated with global numbers of the vertices */
    /* All processors generate the same random number corresponding
       to the same global vertex number */
    rand_key = (int *) ZOLTAN_MALLOC(sizeof(int) * lastlno);
    if (lastlno && !rand_key)
	MEMORY_ERROR;
    for (i=0; i<lastlno; i++) {
	/* TODO64 OK if 8byte int is cast as 4byte int here? */
	Zoltan_Srand((unsigned int)Zoltan_G2LHash_L2G(hash, i), NULL);
	rand_key[i] = (int) (((double)Zoltan_Rand(NULL)/(double)ZOLTAN_RAND_MAX)*1000000);
    }

    /* Color internal vertices and determine the visit order */
    if (get_times) {
	MPI_Barrier(zz->Communicator);
	times[2] = Zoltan_Time(zz->Timer);
    }
    visitIntern = visit + nbound; /* Start of internal vertex visit order. Used with I and B options below */
    if (coloring_order == 'I') {
	InternalColoring(zz, coloring_problem, &nColor, nintvisit, visitIntern, xadj, adj, color, mark, gmaxcolor, coloring_method);
	nConflict = nboundvisit;
    }
    else if (coloring_order=='B')
	nConflict = nboundvisit;
    else if (coloring_order=='U' || coloring_order=='N' || coloring_order=='L' || coloring_order=='S') {
        nConflict = nvtx;
        if (coloring_order=='U' || coloring_order=='N') {
            for (i=0; i<nvtx; i++)
                visit[i] = i;
        }   
        else if (coloring_order=='L') {
            LargestDegreeFirstOrdering(zz, visit, xadj, nvtx, lmaxdeg);
        }
        else if (coloring_order=='S') {            
	    SmallestDegreeLastOrdering(zz, visit, xadj, adj, nvtx, lmaxdeg);  
        }
        if (zz->Num_Proc==1)
            InternalColoring(zz, coloring_problem, &nColor, nvtx, visit, xadj, adj, color, mark, gmaxcolor, coloring_method);
    }


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
				      coloring_problem, coloring_method, comm_pattern,
				      rreqfromC, repliesC, sreqsC, rreqsC,
				      rreqfromF, repliesF, sreqsF, rreqsF,
				      rreqfromFx, repliesFx, sreqsFx, rreqsFx,
				      stats, confChk, ssendsize, srecsize, ssendbuf, srecbuf, ssp, srp, xfp, xfpMark, 
                                      wset, &wsize);
	    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
		ZOLTAN_COLOR_ERROR(ierr, "Error in D2ParallelColoring");
	    /* Allreduce below avoids insufficient memset in conflict detection.
	    If a received color of a non-local vertex is greater than nColor on this processor,
	    then arrays accessed with color number will cause seg fault if this is not used. We may do better */
	    lColor = nColor;
	    MPI_Allreduce(&lColor, &nColor, 1, MPI_INT, MPI_MAX, zz->Communicator);

	    ierr = D2DetectConflicts(zz, coloring_problem, hash, nvtx, wset, wsize, xadj, adj, adjproc, &nColor, color, conflicts, rand_key,
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
    if (coloring_order == 'B')
	InternalColoring(zz, coloring_problem, &nColor, nintvisit, visitIntern, xadj, adj, color, mark, gmaxcolor, coloring_method);


#if 0
    printf("[%d] vtx(gno)->color: ", zz->Proc);
    for (i=0; i<nvtx; i++)
	printf("%d(" ZOLTAN_GNO_SPEC ")->%d ", i, Zoltan_G2LHash_L2G(hash, i), color[i]);
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
    char color_problem, /* Coloring problem */
    int nvtx,         /* In: Number of objects */
    int *xadj,        /* In: Pointer to xadj */
    int **pxbadj,     /* Out: Pointer to start of cut edges in the
			 adj lists of boundary vertices*/
    int *adj,         /* In/Out: Pointer to adj */
    int *adjproc,     /* In: Part vector */
    int *nbound,      /* Out: Number of boundary vertices */
    int *isbound,     /* Out: Indicates if a vertex is on the boundary */
    int *visit,       /* Out: Visit order */
    int *partialD2,   /* In: binary array showing which vertices will be colored */
    int *nintvisit,   /* Out: Number of internal vertices to be colored */
    int *nboundvisit  /* Out: Number of boundary vertices to be colored */
)
{
    static char *yo = "ReorderGraph";
    int *xbadj = NULL;     /* pointer to start of cut edges in the
			      adj lists of boundary vertices    */
    int i, j, iidx, bidx;
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
    iidx = *nbound;
    bidx = 0;
    for (i=0; i<nvtx; ++i) {
	if (color_problem != 'P' || partialD2[i] == 1) {
	    if (isbound[i] == 1)
		visit[bidx++] = i;
	    else {
		visit[iidx++] = i;
	    }
	}
    }
    *nintvisit = iidx-*nbound;
    *nboundvisit = bidx;


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
    char coloring_method,
    int u,
    int oldColor,
    int *nColor,
    int *mark
)
{
    int c, pickedColor;

    switch(coloring_method){
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
    char coloring_problem,
    int *nColor,
    int nvtx,
    int *visit,
    int * xadj,
    int *adj,
    int *color,
    int *mark,
    int gmaxdeg,
    char coloring_method
)
{
    int i, j, k, c, u, v, w;
    static char *yo = "InternalColoring";
    int ierr = ZOLTAN_OK;

    memset(mark, 0xff, gmaxdeg * sizeof(int));
    if (coloring_problem == '1') {
	for (i=0; i<nvtx; ++i) {
	    u = visit[i];
	    for (j = xadj[u]; j < xadj[u+1]; ++j) {
		v = adj[j];
		if ((c = color[v]) != 0)
		    mark[c] = u;
	    }
	    color[u] = PickColor(zz, coloring_method, u, color[u], nColor, mark);
	}
    } else if (coloring_problem == '2' || coloring_problem == 'P') {
	for (i=0; i<nvtx; ++i) {
	    u = visit[i];
	    for (j = xadj[u]; j < xadj[u+1]; ++j) {
		v = adj[j];
		if ((c = color[v]) != 0)
		    mark[c] = u;
		for (k = xadj[v]; k < xadj[v+1]; ++k) {
		    w = adj[k];
		    if ((c = color[w]) != 0)
			mark[c] = u;
		}
	    }
	    color[u] = PickColor(zz, coloring_method, u, color[u], nColor, mark);
	}
    } else
	ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Unknown coloring problem");

 End:

    return ierr;
}

/*****************************************************************************/
/* Parallel coloring of boundary vertices */
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
    ZOLTAN_GNO_TYPE **newcolored, /* Array used for communicating boundary vertex
			 colors at the end of supersteps. Global number
			 of the vertex - color of vertex pairs are filled
			 in this array */
    int *mark,        /* Array to mark forbidden colors for a vertex */
    int gmaxdeg,      /* Maximum vertex degree in the graph */
    G2LHash *hash,    /* hash to map global ids of local and D1 neighbor
			 vertices to consecutive local ids */
    char coloring_method, /* Coloring method. (F) First fit
			  (S) staggered first fit (L) load balancing */
    char comm_pattern, /* (A) asynchronous (S) synchronous supersteps */
    int *rreqfrom,    /* Arrays used for MPI communication */
    int *replies,
    MPI_Request *sreqs,
    MPI_Request *rreqs,
    MPI_Status *stats,
    int *xrelproc, int *relproc, ZOLTAN_GNO_TYPE **persSbuf, int *Ssize, int plstcnt, int *plst
)
{
    static char *yo="D1ParallelColoring";
    int colortag=1001, i, j, p, q, l;
    int n=0;
    int rreqcnt=0, sreqcnt=0, repcount;
    int ierr;
    MPI_Datatype gno_mpi_type;
    ZOLTAN_GNO_TYPE *colored=NULL;

    int flag = 0; /* if set to one, will color all the vertices with color 1 */
    gno_mpi_type = Zoltan_mpi_gno_type();

    memset(Ssize, 0, sizeof(int) * zz->Num_Proc);

    colored = newcolored[zz->Proc];

    /* Issue async recvs */
    for (rreqcnt = i = 0; i < plstcnt; ++i) {
	p = plst[i];
	rreqfrom[rreqcnt] = p;
	if (MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, colortag, zz->Communicator, &rreqs[rreqcnt]))
	    ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "MPI_Irecv failed.");
	++rreqcnt;
    }

    /* at first iteration of recoloring nColor is set to -1 to enable flag so that color 1 is manually assigned */
    if (*nColor < 0) {
        flag = 1;
        *nColor = 0;
    }
    /* Coloring */
    for (i=0; i<nvtx; ++i) {
        int u = visit[i];
        ZOLTAN_GNO_TYPE gu;
        if (!flag) {
            for (j=xadj[u]; j<xadj[u+1]; ++j) {
                int gv = adj[j], c;
                if ((c = color[gv]) != 0) {
                    while (*nColor < c) /* if nColor < c, set nColor to c and reset intermediate mark entries. May be done better. */
                        mark[++(*nColor)] = -1;
                    mark[c] = u;
                }
            }
            color[u] = PickColor(zz, coloring_method, u, color[u], nColor, mark);
        }
        else { /* at first iteration of recoloring, all vertices in that color class is colored with 1, in order not to lose time */
            color[u] = 1;
            *nColor = 1;
        }

	if (!isbound[u]) /* send only boundary vertices */
	    continue;

	gu = Zoltan_G2LHash_L2G(hash, u);
	for (j=xrelproc[isbound[u]-1]; j<xrelproc[isbound[u]]; ++j) {
	    int ap = relproc[j];
	    persSbuf[ap][Ssize[ap]++] = gu;
	    persSbuf[ap][Ssize[ap]++] = (ZOLTAN_GNO_TYPE)color[u];
	}
	n += 2;

	/* If superstep is finished, communicate */
	if (n >= 2*ss) {
	    for (sreqcnt = j = 0; j<plstcnt; ++j) {
		p = plst[j];
                if (Ssize[p] || (comm_pattern == 'S')) {
                    if (Ssize[p] < 2*ss) {
                        persSbuf[p][Ssize[p]++] = -2;
                        persSbuf[p][Ssize[p]++] = -2;
                    }
                    MPI_Isend(persSbuf[p], Ssize[p], gno_mpi_type, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
                    sreqcnt++;
                    Ssize[p] = 0;
                    }
	    }

	    /* Wait color lists from other processors */
	    if (comm_pattern == 'S') { /* wait all results if sync communication*/
		MPI_Waitall(rreqcnt, rreqs, stats);
		repcount = rreqcnt;
	    }
	    else /* if async communication: check if any arrived */
		MPI_Waitsome(rreqcnt, rreqs, &repcount, replies, stats);

	    for (l = repcount-1; l >= 0; --l) {
		ZOLTAN_GNO_TYPE v=0;

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
			c = (int)newcolored[p][j++];
			color[hv] = c;
		    } else
			++j;
		}
		/* If p hasn't finished coloring, issue new color request */
		if (v!=-1) {
		    if (MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, colortag, zz->Communicator, &rreqs[q]))
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
    for (j=0; j<plstcnt; ++j) {
	p = plst[j];
	persSbuf[p][Ssize[p]++] = -1;
	persSbuf[p][Ssize[p]++] = -1;
    }

    for (sreqcnt = j = 0; j<plstcnt; ++j) {
	p = plst[j];
	MPI_Isend(persSbuf[p], Ssize[p], gno_mpi_type, p, colortag, zz->Communicator, &sreqs[sreqcnt]);
        sreqcnt++;
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
	    ZOLTAN_GNO_TYPE v=0;
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
		    c = (int)newcolored[p][j++];
		    color[hv] = c;
		} else
		    ++j;
	    }

	    /* If p hasn't finished coloring, issue new color request */
	    if (v!=-1) {
		if (MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, colortag, zz->Communicator, &rreqs[q]))
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

static int sendNextStepsForbiddenColors(ZZ *zz, G2LHash *hash, int nlvtx, int p, ZOLTAN_GNO_TYPE **srp, 
                                 int *xadj, int *adj, int *xadjnl, int *adjnl,
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
    ZOLTAN_GNO_TYPE **newcolored, /* Array used for communicating boundary vertex
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
    char coloring_problem, /* D2 or partial D2 */
    char coloring_method, /* Coloring method. (F) First fit
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
    ZOLTAN_GNO_TYPE **ssendbuf, /* the size of ssendbuf[i] is ssendsize[i] */
    ZOLTAN_GNO_TYPE **srecbuf,
    ZOLTAN_GNO_TYPE **ssp, /* ssp[i] is a pointer inside ssendbuf[i] memory */
    ZOLTAN_GNO_TYPE **srp,
    int **xfp,
    int *xfpMark,
    int *wset,
    int *wsize
)
{
    int i, j, k, p, l;
    MPI_Datatype gno_mpi_type;
    ZOLTAN_GNO_TYPE *colored = newcolored[zz->Proc];
    int n=0, m = 0;
    int rreqcntC=0, sreqcntC=0; /*keeps track of updated D-1 neighbor colors comm.*/
    int rreqcntFx=0; /*keeps track of forbidden color pointer communication*/
    int sncnt = 0; /*superstep number counter */
    int fp = 0; /*index to forbidProc*/
    int ierr = ZOLTAN_OK;

    gno_mpi_type = Zoltan_mpi_gno_type();

    /* DETERMINE THE SUPERSTEP NUMBER OF LOCAL VERTICES AND COMMUNICATE THEM */
    /* Issue async recvs for superstep info */
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
	if (p!=zz->Proc)
	    MPI_Irecv(srecbuf[p], srecsize[p], gno_mpi_type, p, SNTAG, zz->Communicator, &rreqsC[rreqcntC++]);


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
        ZOLTAN_GNO_TYPE gu;
	int u = (int)visit[i];
	gu = Zoltan_G2LHash_L2G(hash, u);
	for (j=xadj[u]; j<xbadj[u]; j++) {
	    if (ssp[adjproc[j]] != ssendbuf[adjproc[j]]) {
		if (*(ssp[adjproc[j]]-1) != gu) {
		    *(ssp[adjproc[j]]++) = gu;
		}
	    } else {
		*(ssp[adjproc[j]]++) = gu;
            }
	}
	if (n == ss) { /* mark the end of superstep */
	    n = 0;
	    for (p=0; p<zz->Num_Proc; p++)
		if (p != zz->Proc) {
		    *(ssp[p]++) = -1;
		}
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
	    MPI_Isend(ssendbuf[p], ssp[p] - ssendbuf[p] + 1, gno_mpi_type, p, SNTAG, zz->Communicator, &sreqsC[sreqcntC++]);


    MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all superstep info to be received */
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all superstep info to be sent */

    /***** COLORING *****/
    /* Issue async recvs for updated d1 colors */
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
	if (p != zz->Proc) {
	    rreqfromC[rreqcntC] = p;
	    MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, COLORTAG, zz->Communicator, &rreqsC[rreqcntC++]);
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
	int u = (int)visit[i];
	if (/*coloring_problem != 'P' && */confChk[u] == sncnt) { /* add boundary vertices of this round into confChk */
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
	color[u] = PickColor(zz, coloring_method, u, color[u], nColor, mark);
	colored[n++] = Zoltan_G2LHash_L2G(hash, u);
	colored[n++] = (ZOLTAN_GNO_TYPE)color[u];

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
		    MPI_Isend(colored, 2*ss, gno_mpi_type, p, COLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);

	    /* process the received color update messages */
	    MPI_Waitall(rreqcntC, rreqsC, stats);
	    fp = 0; /*number of procesors requesting forbidden colors for their next ss*/
	    for (l = 0; l < rreqcntC; ++l) {
		p = rreqfromC[l];
		/* update non-local vertex colors according to the received color update message content */
		for (j = 0; j < 2*ss; ) {
		    ZOLTAN_GNO_TYPE v = newcolored[p][j++];
		    int c, hv;
		    if (v < 0)
			break;
		    if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
			c = (int)newcolored[p][j++];
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
		MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, COLORTAG, zz->Communicator, &rreqsC[l]);
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
	    MPI_Isend(colored, n, gno_mpi_type, p, COLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);

    /* process the remaining color update messages */
    while (rreqcntC) { /*until all color update messages in this superstep are received*/
	MPI_Waitall(rreqcntC, rreqsC, stats); /* wait all color updates to be received */
	fp = 0;
	++sncnt;
	for (l = 0; l < rreqcntC; ++l) {
	    p = rreqfromC[l];
	    /* update non-local vertex colors according to the received color update message content */
	    for (j = 0; j < 2*ss; ) {
		ZOLTAN_GNO_TYPE v = newcolored[p][j++];
		int c, hv;
		if (v < 0)
		    break;
		if ((hv = Zoltan_G2LHash_G2L(hash, v)) != -1) {
		    c = (int)newcolored[p][j++];
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
	    MPI_Irecv(newcolored[p], 2*ss, gno_mpi_type, p, COLORTAG, zz->Communicator, &rreqsC[l]);
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

static int D2DetectConflicts(ZZ *zz, char coloring_problem, G2LHash *hash, int nlvtx, int *wset, int wsize, int *xadj, int *adj, int *adjproc, int *nColor, int *color, int *conflicts, int *rand_key, int *vmark, int *seen, int *where, int *pwhere, ZOLTAN_GNO_TYPE **rcsendbuf, ZOLTAN_GNO_TYPE **rcrecbuf, int *rcsendsize, int *rcrecsize, ZOLTAN_GNO_TYPE **srp, MPI_Request *sreqsC, MPI_Request *rreqsC, int *rreqfromC, MPI_Status *stats, int *nconflict)
{
    static char *yo="D2DetectConflicts";
    int w, i, j, p;
    int rreqcntC=0, sreqcntC=0;
    ZOLTAN_GNO_TYPE *q;
    MPI_Datatype gno_mpi_type;
    int ierr = ZOLTAN_OK;

    gno_mpi_type = Zoltan_mpi_gno_type();

    for (i=0; i<zz->Num_Proc; ++i)
      srp[i] = rcsendbuf[i];

    /*** Issue async recv for recolor info ***/
    for (rreqcntC=p=0; p<zz->Num_Proc; ++p)
	if (p != zz->Proc) {
	    rreqfromC[rreqcntC] = p;
	    MPI_Irecv(rcrecbuf[p], rcrecsize[p], gno_mpi_type, p, RECOLORTAG, zz->Communicator, &rreqsC[rreqcntC++]);
	}

    /*** Detect conflicts and mark vertices to be recolored ***/
    memset(seen, 0xff, sizeof(int) * (*nColor+1));
    memset(where, 0xff, sizeof(int) * (*nColor+1));
    memset(pwhere, 0xff, sizeof(int) * (*nColor+1));

    for (i = 0; i < wsize; i++) {
        ZOLTAN_GNO_TYPE gx, gv;
	int x, v, cw, cx, px, pv;
	w = wset[i];
/*	if (coloring_problem != 'P') { */
	    cw = color[w];
	    seen[cw] = w;
	    where[cw] = w;
	    pwhere[cw] = zz->Proc;
/*	}*/
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
	    MPI_Isend(rcsendbuf[p], srp[p] - rcsendbuf[p] , gno_mpi_type, p, RECOLORTAG, zz->Communicator, &sreqsC[sreqcntC++]);
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
		ZOLTAN_GNO_TYPE gu = rcrecbuf[p][j++];
		int u;
		u = Zoltan_G2LHash_G2L(hash, gu);
		if (!vmark[u]) {
		    *srp[zz->Proc]++ = gu;
		    vmark[u] = 1;
		}
	    }
	}
    *srp[zz->Proc] = -1;
    *nconflict = srp[zz->Proc] - rcsendbuf[zz->Proc];
    for (q = rcsendbuf[zz->Proc]; *q != -1; ++q) {
	*conflicts = Zoltan_G2LHash_G2L(hash, *q);
	vmark[*conflicts] = 0;
        conflicts++;
    }
    /* *q = 0; UVC: this is OK to keep but it is useless; therefore it is commented out*/
    MPI_Waitall(sreqcntC, sreqsC, stats); /* wait all recolor info to be sent */

 End:

    return ierr;
}

/* sorting functions for NONDECREASING and NONINCREASING permutations */
static int pairofintsup_ind_nd(const void* a, const void* b)
{
  const struct foo{int x,y;} *pa = a, *pb = b;
  return pa->y - pb->y;
}

static int pairofintsup_ind_ni(const void* a, const void* b)
{
  const struct foo{int x,y;} *pa = a, *pb = b;
  return pb->y - pa->y;
}

static int Recoloring(ZZ *zz, int recoloring_permutation, int recoloring_type, 
		      int recoloring_num_of_iterations, int nvtx, int carrierbufsize, int *nConflict,
		      int *confCont, int *nTotConflict, int *nRound, int lastlno, int globmaxnvtx,
		      int *nColor, int *color, int *visit,int *xadj, int *adj, int *xbadj,
		      int *conflicts, int *rand_key, G2LHash *hash, int *isbound,
		      ZOLTAN_GNO_TYPE **newcolored, int *mark, int gmaxdeg,
		      char coloring_method, char comm_pattern, int *rreqfrom, int *replies,
		      MPI_Request *sreqs, MPI_Request *rreqs, MPI_Status *stats, 
		      int *xrelproc, int *relproc, ZOLTAN_GNO_TYPE **persSbuf,
		      int *Ssize, int plstcnt, int *plst)
{
  int i, j, isSequential;
  int ierr = ZOLTAN_OK;
  static char *yo = "Recoloring";
  int globnumofcolor = 0;
  int numofcolor = 0; 
  int marksize, nStart, nEnd, length, color_of_i;
  int *howmanynodeinthatcolor = NULL; /* number of vertices in the color classes of a proc*/
  int *howmanynodeinthatcolorglobal = NULL; /* number of vertices in the color classes of all procs */
  int *sorted_color = NULL; /* sorted colors wrt to number of vertices at them */
  int *permutation = NULL; /* permutation of colors to be used in recoloring permuatations */
  int *dummy = NULL; 
  int *prefixcount = NULL; /* count of vertices at each color classes in the specif ird recoloring permutation order  */
  int *colorindex = NULL; /* starting indexes of color classes */
  int *tp = NULL;
  int *dummyvisit = NULL;
  
  if (zz->Num_Proc == 1)
      isSequential = 1;
  else
      isSequential = 0;
  

  /* number of colors in a proc is calculated */
  for (i=0; i<nvtx; i++) {
      if (color[i] > numofcolor)
          numofcolor = color[i];
  }
  /* total number of colors in all procs is calculated */
  if (isSequential)
      globnumofcolor = numofcolor;
  else
      MPI_Allreduce(&numofcolor, &globnumofcolor, 1, MPI_INT, MPI_MAX, zz->Communicator);
  
  howmanynodeinthatcolorglobal = (int *)ZOLTAN_MALLOC(sizeof(int) * (globnumofcolor + 1));
  if (!isSequential)
      howmanynodeinthatcolor = (int *)ZOLTAN_MALLOC(sizeof(int)*(globnumofcolor+1));

  sorted_color = (int *)ZOLTAN_MALLOC(sizeof(int) * 2 * (globnumofcolor));
  permutation = (int *)ZOLTAN_MALLOC(sizeof(int)*(globnumofcolor+1));
  dummy = (int *)ZOLTAN_MALLOC(sizeof(int)*(globnumofcolor+1));
  prefixcount = (int *)ZOLTAN_MALLOC(sizeof(int)*(globnumofcolor+1));
  colorindex = (int *)ZOLTAN_MALLOC(sizeof(int)*(globnumofcolor+1));

  if (recoloring_permutation == REVERSE)
      dummyvisit = (int *)ZOLTAN_MALLOC(nvtx * sizeof(int));
  
  for (j=0; j<recoloring_num_of_iterations; j++) {
      if (j != 0) { /*already done for iteration 1 */
          /* number of colors in a proc is calculated */
          for (i=0; i<nvtx; i++) {
	      if (color[i] > numofcolor)
	          numofcolor = color[i];
          }

          /* total number of colors in all procs is calculated */
          if (isSequential)
	      globnumofcolor = numofcolor;
          else
	      MPI_Allreduce(&numofcolor, &globnumofcolor, 1, MPI_INT, MPI_MAX, zz->Communicator);
      }

      memset(howmanynodeinthatcolorglobal, 0, sizeof(int) * (globnumofcolor + 1));
      
      /* number of vertices at each color classes in a proc is calculated */
      if (isSequential) {
          for (i=0; i<nvtx; i++)
	      howmanynodeinthatcolorglobal[color[i]]++;
      }
      else {
          memset(howmanynodeinthatcolor, 0, sizeof(int) * (globnumofcolor+1));
          for (i=0; i<nvtx; i++)
              howmanynodeinthatcolor[color[i]]++;
      }
      
      /* color classes are sorted wrt specif ied recoloring permuatation */
      if ((recoloring_permutation == NONDECREASING) || (recoloring_permutation == NONINCREASING)) {
          /* number of vertices at each color classes in all procs is calculated */
          if (!isSequential)
              MPI_Allreduce(howmanynodeinthatcolor, howmanynodeinthatcolorglobal, globnumofcolor+1, MPI_INT, MPI_SUM, zz->Communicator);
          
          
          /* colors are sorted based on their counts */
          for (i=0; i < 2*globnumofcolor; i+=2) {
              sorted_color[i] = i/2+1;
              sorted_color[i+1] = howmanynodeinthatcolorglobal[i/2+1];
          }
          
          if (recoloring_permutation == NONDECREASING)
              qsort (sorted_color, globnumofcolor, 2 * sizeof(int), pairofintsup_ind_nd);
          else
              qsort (sorted_color, globnumofcolor, 2 * sizeof(int), pairofintsup_ind_ni);
          
          for (i=0; i<2*globnumofcolor; i+=2)
	      permutation[sorted_color[i]] = i/2+1;
      
          memset(dummy, 0, sizeof(int) * (globnumofcolor+1));
          /* color ids are changed acc to new permutation calculated */
          for (i=0; i<lastlno; i++)
	      color[i] = permutation[color[i]];
          if (isSequential) {
	      for (i=1; i<globnumofcolor+1; i++)
	          dummy[permutation[i]] = howmanynodeinthatcolorglobal[i];
	      for (i=1; i<globnumofcolor+1; i++)
                  howmanynodeinthatcolorglobal[i] = dummy[i];
          }
          else {
              for (i=1; i<globnumofcolor+1; i++)
                  dummy[permutation[i]] = howmanynodeinthatcolor[i];
              memset(howmanynodeinthatcolor, 0, sizeof(int) * (globnumofcolor+1));
              for (i=1; i<globnumofcolor+1; i++)
                  howmanynodeinthatcolor[i] = dummy[i];
          }  
      }   
      
      /* color indexes are stored to keep track of each color classes */
      memset(prefixcount, 0, sizeof(int) * (globnumofcolor+1));
      memset(colorindex, 0, sizeof(int) * (globnumofcolor+1));
      colorindex[0] = prefixcount[0]=0;
      for (i=1; i<globnumofcolor+1; ++i) {
          if (isSequential)
              colorindex[i] = prefixcount[i]=prefixcount[i-1]+howmanynodeinthatcolorglobal[i];
          else
              colorindex[i] = prefixcount[i]=prefixcount[i-1]+howmanynodeinthatcolor[i];
      }
      for (i=0; i<nvtx; ++i) {
          color_of_i = color[i];
          prefixcount[color_of_i]--;
          visit[prefixcount[color_of_i]] = i;
      }                   
      memset(color, 0, sizeof(int) * lastlno);
      
      if (!isSequential) {
          tp = visit;
          /* for  synchronous coloring, recoloring is done for  each color classes separately,i.e. relevant neigs wait each other while  coloring same colored vertices*/
          if (recoloring_type == SYNCHRONOUS) {
              carrierbufsize = globmaxnvtx;
              marksize = *nColor;
              *nColor = -1;
              for (i=1; i<globnumofcolor+1; i++) {
                  if (recoloring_permutation == REVERSE) {
                      nStart = colorindex[globnumofcolor-i];
                      nEnd = colorindex[globnumofcolor-i+1];
                  }
                  else {
                      nStart = colorindex[i-1];
                      nEnd = colorindex[i];
                  }	  
                  length = nEnd-nStart;
                  memset(mark, 0xff, (1+*nColor) * sizeof(int));
                  ierr = D1ParallelColoring(zz, length, visit+nStart, xadj, adj, isbound, carrierbufsize,
                                            nColor, color, newcolored, mark, gmaxdeg, hash,
                                            coloring_method, comm_pattern, rreqfrom, replies,
                                            sreqs, rreqs, stats, xrelproc, relproc, persSbuf,
                                            Ssize, plstcnt, plst);
                  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                      ZOLTAN_COLOR_ERROR(ierr, "Error in D1ParallelColoring");
              }
          }
          /* for  asynchronous coloring, a new visiting order is obtained based on previous coloring and recoloring is done just like previous coloring */
          else if (recoloring_type == ASYNCHRONOUS) {
              if (recoloring_permutation == REVERSE) {
                  for (i=0; i<nvtx; i++)
                      dummyvisit[nvtx-i-1] = visit[i];
                  for (i=0; i<nvtx; i++)
                      visit[i] = dummyvisit[i];
              }
              
              *nConflict = nvtx;
              *confCont = 1; /* dummy to put asynch into while  loop */
              while (*confCont) {
                  tp = visit;
                  memset(mark, 0xff, (1+*nColor) * sizeof(int)); /* reset dirty entries */
                  ierr = D1ParallelColoring(zz, *nConflict, visit, xadj, adj, isbound, carrierbufsize,
                                            nColor, color, newcolored, mark, gmaxdeg, hash,
                                            coloring_method, comm_pattern, rreqfrom, replies,
                                            sreqs, rreqs, stats, xrelproc, relproc, persSbuf,
                                            Ssize, plstcnt, plst);
                  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
                      ZOLTAN_COLOR_ERROR(ierr, "Error in D1ParallelColoring");
                  
                  *nConflict = D1DetectConflicts(zz, *nConflict, visit, xadj, xbadj, adj, color, conflicts, rand_key, hash);
                  visit = conflicts;
                  conflicts = tp;
                  *confCont = 0;
                  MPI_Allreduce(nConflict, confCont, 1, MPI_INT, MPI_SUM, zz->Communicator);
                  *nTotConflict += *confCont;
                  ++(*nRound);
              }
          }
      }
      else {
	  *nColor = 0;
          if (recoloring_permutation == REVERSE) {
              for (i=0; i<nvtx; i++)
                  dummyvisit[nvtx-i-1] = visit[i];
              for (i=0; i<nvtx; i++)
                  visit[i] = dummyvisit[i];
          }  
          InternalColoring(zz, '1', nColor, nvtx, visit, xadj, adj, color, mark, gmaxdeg, coloring_method);
      }
  }
  
  ierr = ZOLTAN_OK;
  
End:
  
  ZOLTAN_FREE(&howmanynodeinthatcolorglobal);
  ZOLTAN_FREE(&sorted_color);
  if (permutation)
      ZOLTAN_FREE(&permutation);
  if (dummy)
      ZOLTAN_FREE(&dummy);
  if (!isSequential)
      ZOLTAN_FREE(&howmanynodeinthatcolor);
  ZOLTAN_FREE(&prefixcount);
  ZOLTAN_FREE(&colorindex);
  if (dummyvisit)
      ZOLTAN_FREE(&dummyvisit);
  
  return ierr;
}
