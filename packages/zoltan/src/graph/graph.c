// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_matrix.h"
#include "graph.h"
#include "graph_const.h"
#include "params_const.h"
#include "graph_util.h"
#include "graph_params.h"

/* #define CC_TIMERS */

#define AFFECT_NOT_NULL(ptr, src) do { if ((ptr) != NULL) (*(ptr)) = (src); } while (0)

/* This function needs a distribution : rows then cols to work properly */
int
Zoltan_ZG_Build (ZZ* zz, ZG* graph, int local,
  int request_GNOs,                /* Input:  Flag indicating calling code 
                                              needs translation of extra GIDs
                                              to GNOs; partial 2D coloring
                                              needs this feature. */
  int num_requested,               /* Input:  Local # of GIDs needing 
                                              translation to GNOs. */
  ZOLTAN_ID_PTR requested_GIDs,    /* Input:  Calling code requests the 
                                              GNOs for these GIDs */
  ZOLTAN_GNO_TYPE *requested_GNOs  /* Output: Return GNOs of 
                                              the requested GIDs.  */
)
{
  static char *yo = "Zoltan_ZG_Build";
  int ierr = ZOLTAN_OK;
  int diag;
  int *diagarray=NULL;
  Zoltan_matrix_options opt;
  char symmetrization[MAX_PARAM_STRING_LEN];
  char bipartite_type[MAX_PARAM_STRING_LEN];
  char weigth_type[MAX_PARAM_STRING_LEN];
  char matrix_build_type[MAX_PARAM_STRING_LEN];
  int graph_fast_build_base;
  int bipartite = 0;
#ifdef CC_TIMERS
  double times[9]={0.,0.,0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
  double gtimes[9]={0.,0.,0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
  char *timenames[9]= {"", "setup", "matrix build", "diag", "symmetrize", "dist lin", "2D dist", "complete", "clean up"};

  MPI_Barrier(zz->Communicator);
  times[0] = Zoltan_Time(zz->Timer);
#endif /* CC_TIMERS */

  ZOLTAN_TRACE_ENTER(zz, yo);
  memset (graph, 0, sizeof(ZG));

  /* Read graph build parameters */
  Zoltan_Bind_Param(ZG_params, "GRAPH_SYMMETRIZE", (void *) &symmetrization);
  Zoltan_Bind_Param(ZG_params, "GRAPH_SYM_WEIGHT", (void *) &weigth_type);
  Zoltan_Bind_Param(ZG_params, "GRAPH_BIPARTITE_TYPE", (void *) &bipartite_type);
  Zoltan_Bind_Param(ZG_params, "GRAPH_BUILD_TYPE", (void*) &matrix_build_type);
  Zoltan_Bind_Param(ZG_params, "GRAPH_FAST_BUILD_BASE", (void*) &graph_fast_build_base);

  /* Set default values */
  strncpy(symmetrization, "NONE", MAX_PARAM_STRING_LEN);
  strncpy(bipartite_type, "OBJ", MAX_PARAM_STRING_LEN);
  strncpy(weigth_type, "ADD", MAX_PARAM_STRING_LEN);
  strncpy(matrix_build_type, "NORMAL", MAX_PARAM_STRING_LEN);
  graph_fast_build_base = 0;

  Zoltan_Assign_Param_Vals(zz->Params, ZG_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);

  Zoltan_Matrix2d_Init(&graph->mtx);

  graph->mtx.comm = (PHGComm*)ZOLTAN_MALLOC (sizeof(PHGComm));
  if (graph->mtx.comm == NULL) MEMORY_ERROR;
  Zoltan_PHGComm_Init (graph->mtx.comm);

  memset(&opt, 0, sizeof(Zoltan_matrix_options));
  opt.enforceSquare = 1;      /* We want a graph: square matrix */
  if (!strcasecmp(weigth_type, "ADD"))
    opt.pinwgtop = ADD_WEIGHT;
  else if (!strcasecmp(weigth_type, "MAX"))
    opt.pinwgtop = MAX_WEIGHT;
  else if (!strcasecmp(weigth_type, "CMP"))
    opt.pinwgtop = MAX_WEIGHT;
  opt.pinwgt = 1;
  opt.randomize = 0;
  opt.local = local;
  opt.keep_distribution = 1;
  if (strcasecmp(symmetrization, "NONE")) {
    opt.symmetrize = 1;
  }
  if (!strcasecmp(matrix_build_type, "FAST"))
    opt.speed = MATRIX_FAST;
  else if (!strcasecmp(matrix_build_type, "FAST_NO_DUP"))
    opt.speed = MATRIX_NO_REDIST;
  else
    opt.speed = MATRIX_FULL_DD;
  opt.fast_build_base = graph_fast_build_base;

#ifdef CC_TIMERS
  times[1] = Zoltan_Time(zz->Timer);
#endif

  ierr = Zoltan_Matrix_Build(zz, &opt, &graph->mtx.mtx, request_GNOs,
                             num_requested, requested_GIDs, requested_GNOs);
  CHECK_IERR;

#ifdef CC_TIMERS
  times[2] = Zoltan_Time(zz->Timer);
#endif

  ierr = Zoltan_Matrix_Mark_Diag (zz, &graph->mtx.mtx, &diag, &diagarray);
  CHECK_IERR;
  if (diag) { /* Some Diagonal Terms have to be removed */
    ierr = Zoltan_Matrix_Delete_nnz(zz, &graph->mtx.mtx, diag, diagarray);
    ZOLTAN_FREE(&diagarray);
    CHECK_IERR;
  }

#ifdef CC_TIMERS
  times[3] = Zoltan_Time(zz->Timer);
#endif

  if (opt.symmetrize) {
    if (!strcasecmp(symmetrization, "BIPARTITE"))
      bipartite = 1;
    ierr = Zoltan_Matrix_Sym(zz, &graph->mtx.mtx, bipartite);
    CHECK_IERR;
  }

#ifdef CC_TIMERS
  times[4] = Zoltan_Time(zz->Timer);
#endif

  ierr = Zoltan_Distribute_LinearY(zz, graph->mtx.comm);
  CHECK_IERR;

#ifdef CC_TIMERS
  times[5] = Zoltan_Time(zz->Timer);
  MPI_Barrier(zz->Communicator);
#endif
  ierr = Zoltan_Matrix2d_Distribute (zz, graph->mtx.mtx, &graph->mtx, 0);
  CHECK_IERR;

#ifdef CC_TIMERS
  times[6] = Zoltan_Time(zz->Timer);
#endif
  ierr = Zoltan_Matrix_Complete(zz, &graph->mtx.mtx);

#ifdef CC_TIMERS
  times[7] = Zoltan_Time(zz->Timer);
#endif

  if (bipartite) {

/*    int vertlno; */
/*    int limit; */
/*    int offset; */
    graph->bipartite = 1;
    graph->fixed_vertices = graph->mtx.mtx.ybipart;
/*     graph->fixed_vertices = (int*) ZOLTAN_MALLOC(graph->mtx.mtx.nY*sizeof(int)); */
/*     if (graph->mtx.mtx.nY && graph->fixed_vertices == NULL) MEMORY_ERROR; */
/*     limit = graph->mtx.mtx.offsetY; */
/*     /\* What kind of vertices do we want to keep ? *\/ */
/*     graph->fixObj = !strcasecmp(bipartite_type, "OBJ"); /\* Non-zero value means "objects" *\/ */

/*     offset = graph->mtx.mtx.offsetY - graph->mtx.dist_y[graph->mtx.comm->myProc_y]; */
/*     if (graph->fixObj) /\* What kind of vertices do we want to keep ? *\/ */
/*       for (vertlno = 0 ; vertlno < graph->mtx.mtx.nY ; ++ vertlno) */
/* 	graph->fixed_vertices[vertlno] = (vertlno < offset); */
/*     else */
/*       for (vertlno = 0 ; vertlno < graph->mtx.mtx.nY ; ++ vertlno) */
/* 	graph->fixed_vertices[vertlno] = (vertlno >= offset); */
  }

#ifdef CC_TIMERS
  MPI_Barrier(zz->Communicator);
  times[8] = Zoltan_Time(zz->Timer);

  MPI_Reduce(times, gtimes, 9, MPI_DOUBLE, MPI_MAX, 0, zz->Communicator);
  if (!zz->Proc) {
      int i;
      printf("Total Build Time in Proc-0: %.2lf    Max: %.2lf\n", times[8]-times[0], gtimes[8]-times[0]);
      for (i=1; i<9; ++i)
          printf("%-13s in Proc-0: %8.2lf  Max: %8.2lf\n", timenames[i],  times[i]-times[i-1], gtimes[i]-gtimes[i-1]);
  }
#endif

 End:
  ZOLTAN_FREE(&diagarray);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

int
Zoltan_ZG_Export (ZZ* zz, const ZG* const graph, ZOLTAN_GNO_TYPE *gvtx, int *nvtx,
		  int *obj_wgt_dim, int *edge_wgt_dim,
		  ZOLTAN_GNO_TYPE **vtxdist, int **xadj, ZOLTAN_GNO_TYPE **adjncy, int **adjproc,
		  float **ewgt, int **partialD2)
{
  AFFECT_NOT_NULL(gvtx, graph->mtx.mtx.globalY);
  AFFECT_NOT_NULL(nvtx, graph->mtx.mtx.nY);
  AFFECT_NOT_NULL(vtxdist, graph->mtx.dist_y);
  AFFECT_NOT_NULL(xadj, graph->mtx.mtx.ystart); 
  AFFECT_NOT_NULL(adjncy, graph->mtx.mtx.pinGNO);
  AFFECT_NOT_NULL(partialD2, graph->fixed_vertices);

  AFFECT_NOT_NULL(obj_wgt_dim, graph->mtx.mtx.ywgtdim);
  AFFECT_NOT_NULL(edge_wgt_dim, graph->mtx.mtx.pinwgtdim);

  AFFECT_NOT_NULL(ewgt, graph->mtx.mtx.pinwgt);

  return Zoltan_Matrix2d_adjproc(zz, &graph->mtx, adjproc);
}

int
Zoltan_ZG_Vertex_Info(ZZ* zz, const ZG *const graph,
		      ZOLTAN_ID_PTR *pgid, ZOLTAN_ID_PTR *plid, float **pwwgt, int **pinput_part) {
  static char *yo = "Zoltan_ZG_Vertex_Info";
  int ierr = ZOLTAN_OK;
  float *wgt = NULL;
  int *input_part = NULL;
  ZOLTAN_ID_PTR lid = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  AFFECT_NOT_NULL(pgid, graph->mtx.mtx.yGID);
  if (pwwgt != NULL) {
    wgt = *pwwgt = (float*) ZOLTAN_MALLOC(graph->mtx.mtx.nY*zz->Obj_Weight_Dim*sizeof(float));
    if (graph->mtx.mtx.nY >0 && zz->Obj_Weight_Dim > 0 && *pwwgt == NULL) MEMORY_ERROR;
  }
  if (pinput_part != NULL) {
    input_part = *pinput_part = (int *) ZOLTAN_MALLOC(graph->mtx.mtx.nY*sizeof(int));
    if (graph->mtx.mtx.nY > 0 && *pinput_part == NULL) MEMORY_ERROR;
  }
  if (plid != NULL) {
    lid = *plid = ZOLTAN_MALLOC_LID_ARRAY(zz, graph->mtx.mtx.nY);
    if (graph->mtx.mtx.nY >0 && zz->Num_LID >0 && *plid == NULL)
      MEMORY_ERROR;
  }
  ierr = Zoltan_Matrix_Vertex_Info(zz, &graph->mtx.mtx, lid,
				   wgt, input_part);

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


  /* This function may work on any distribution of the bipartite graph */
int
Zoltan_ZG_Register(ZZ* zz, ZG* graph, int* properties)
{
  static char *yo = "Zoltan_ZG_Register";
  int ierr = ZOLTAN_OK;
  int *props;
  struct Zoltan_DD_Struct *dd;
  int size;
  ZOLTAN_ID_PTR GID;

  ZOLTAN_TRACE_ENTER(zz, yo);
  size = graph->mtx.mtx.nY;
  dd = graph->mtx.mtx.ddY;

  if (graph->bipartite) { /* Need to construct another properties array with only the fixed elements ! */
    int vertlno;

    if (graph->fixObj) {
      dd = graph->mtx.mtx.ddX;
    }
    props = (int*)ZOLTAN_MALLOC(sizeof(int)*size);
    if (graph->mtx.mtx.nY  && props == NULL) MEMORY_ERROR;
    GID = ZOLTAN_MALLOC_GID_ARRAY(zz, size);
    if (size && GID == NULL) MEMORY_ERROR;
    for (size = 0, vertlno = 0 ; vertlno < graph->mtx.mtx.nY ; ++vertlno) {
      if (graph->fixed_vertices[vertlno]) {
	props[size] = properties[vertlno];
	ZOLTAN_SET_GID(zz, GID+ size*zz->Num_GID,
		       graph->mtx.mtx.yGID+vertlno*zz->Num_GID);
	size ++;
      }
    }
  }
  else {
    props = properties;
    GID = graph->mtx.mtx.yGID;
    if (graph->mtx.mtx.ddY == NULL) {
      ierr = Zoltan_DD_Create (&graph->mtx.mtx.ddY, zz->Communicator, 1, zz->Num_GID,
			       sizeof(ZOLTAN_ID_TYPE), graph->mtx.mtx.globalY/zz->Num_Proc, 0);
      CHECK_IERR;
      /* Hope a linear assignment will help a little */
      if (graph->mtx.mtx.globalX/zz->Num_Proc)
        Zoltan_DD_Set_Neighbor_Hash_Fn1(graph->mtx.mtx.ddY,
                                        graph->mtx.mtx.globalX/zz->Num_Proc);
    }
    dd = graph->mtx.mtx.ddY;
  }
  /* Make our new numbering public */
  ierr = Zoltan_DD_Update (dd, GID, NULL, NULL, props, size);
  CHECK_IERR;

  End:
  if (graph->bipartite) {
    ZOLTAN_FREE(&props);
    ZOLTAN_FREE(&GID);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


/* This function may work with any distribution of the bipartite graph */
int
Zoltan_ZG_Query (ZZ* zz, const ZG* const graph,
	  ZOLTAN_ID_PTR GID, int GID_length, int* properties)
{
  struct Zoltan_DD_Struct *dd;

  dd = graph->mtx.mtx.ddY;
/*   if (graph->bipartite && graph->fixObj) */
/*     dd = graph->mtx.mtx.ddX; */

  return Zoltan_DD_Find(dd, GID, NULL, NULL, properties, GID_length, NULL);
}

void
Zoltan_ZG_Free(ZG *graph){
/*   if (graph->bipartite) */
/*     ZOLTAN_FREE(&graph->fixed_vertices); */

  Zoltan_Matrix2d_Free(&graph->mtx);
}


int Zoltan_ZG_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int  index;
  PARAM_UTYPE result;

  return Zoltan_Check_Param(name, val, ZG_params, &result, &index);
}



#ifdef __cplusplus
}
#endif
