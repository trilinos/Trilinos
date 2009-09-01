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

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "matrix.h"
#include "graph.h"

#define CHECK_IERR do {   if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) \
    goto End;  } while (0)

#define AFFECT_NOT_NULL(ptr, src) do { if ((ptr) != NULL) (*(ptr)) = (src); } while (0)

/* At this time function is in parmetis directory but it will change soon */
extern int
Zoltan_Verify_Graph(MPI_Comm comm, int *vtxdist, int *xadj,
		    int *adjncy, int *vwgt, int *adjwgt,
		    int vwgt_dim, int ewgt_dim,
		    int graph_type, int check_graph, int output_level);

/* This function needs a distribution : rows then cols to work properly */
int
Zoltan_ZG_Build (ZZ* zz, ZG* graph, int bipartite, int fixObj, int local)
{
  static char *yo = "ZG_Build";
  int ierr = ZOLTAN_OK;
  int diag;
  int *diagarray=NULL;
  Zoltan_matrix_options opt;

  ZOLTAN_TRACE_ENTER(zz, yo);
  memset (graph, 0, sizeof(ZG));

  graph->mtx.comm = (PHGComm*)ZOLTAN_MALLOC (sizeof(PHGComm));
  if (graph->mtx.comm == NULL) MEMORY_ERROR;

  memset(&opt, 0, sizeof(Zoltan_matrix_options));
  opt.enforceSquare = 1;
  opt.pinwgtop = ADD_WEIGHT;
  opt.pinwgt = 1;
  opt.randomize = 0;
  opt.local = local;
  opt.keep_distribution = 1;

  ierr = Zoltan_Matrix_Build(zz, &opt, &graph->mtx.mtx);
  CHECK_IERR;
  ierr = Zoltan_Matrix_Mark_Diag (zz, &graph->mtx.mtx, &diag, &diagarray);
  CHECK_IERR;
  if (diag) { /* Some Diagonal Terms have to be removed */
    ierr = Zoltan_Matrix_Delete_nnz(zz, &graph->mtx.mtx, diag, diagarray);
    ZOLTAN_FREE(&diagarray);
    CHECK_IERR;
  }

/*   ierr = Zoltan_Matrix_Sym(zz, &graph->mtx.mtx, bipartite); */
/*   CHECK_IERR; */

  ierr = Zoltan_Distribute_LinearY(zz, graph->mtx.comm);
  CHECK_IERR;
  ierr = Zoltan_Matrix2d_Distribute (zz, graph->mtx.mtx, &graph->mtx, 0);
  CHECK_IERR;

  ierr = Zoltan_Matrix_Complete(zz, &graph->mtx.mtx);

  if (bipartite) {
    int vertlno;
    int limit;
    int offset;

    graph->bipartite = 1;
    graph->fixed_vertices = (int*) ZOLTAN_MALLOC(graph->mtx.mtx.nY*sizeof(int));
    if (graph->mtx.mtx.nY && graph->fixed_vertices == NULL) MEMORY_ERROR;
    limit = graph->mtx.mtx.offsetY;
    graph->fixObj = fixObj;

    offset = graph->mtx.mtx.offsetY - graph->mtx.dist_y[graph->mtx.comm->myProc_y];
    if (fixObj)
      for (vertlno = 0 ; vertlno < graph->mtx.mtx.nY ; ++ vertlno)
	graph->fixed_vertices[vertlno] = (vertlno < offset);
    else
      for (vertlno = 0 ; vertlno < graph->mtx.mtx.nY ; ++ vertlno)
	graph->fixed_vertices[vertlno] = (vertlno >= offset);
  }

 End:
  ZOLTAN_FREE(&diagarray);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

int
Zoltan_ZG_Export (ZZ* zz, const ZG* const graph, int *gvtx, int *nvtx,
		  int *obj_wgt_dim, int *edge_wgt_dim,
		  int **vtxdist, int **xadj, int **adjncy, int **adjproc,
		  float **xwgt, float **ewgt, int **partialD2)
{
  int ierr;

  AFFECT_NOT_NULL(gvtx, graph->mtx.mtx.globalY);
  AFFECT_NOT_NULL(nvtx, graph->mtx.mtx.nY);
  AFFECT_NOT_NULL(vtxdist, graph->mtx.dist_y);
  AFFECT_NOT_NULL(xadj, graph->mtx.mtx.ystart);
  AFFECT_NOT_NULL(adjncy, graph->mtx.mtx.pinGNO);
  AFFECT_NOT_NULL(partialD2, graph->fixed_vertices);
  /* I have to convert from float to int */
  AFFECT_NOT_NULL(obj_wgt_dim, graph->mtx.mtx.ywgtdim);
  AFFECT_NOT_NULL(edge_wgt_dim, graph->mtx.mtx.pinwgtdim);
  AFFECT_NOT_NULL(xwgt, graph->mtx.mtx.ywgt);
  AFFECT_NOT_NULL(ewgt, graph->mtx.mtx.pinwgt);

  /* TODO: convert wgt to int to be able to call Zoltan_Verify_Graph */
  ierr = Zoltan_Verify_Graph(zz->Communicator, *vtxdist, *xadj,
			     *adjncy, NULL, NULL,
			     0, 0,
			     0, 2, 2);

  return Zoltan_Matrix2d_adjproc(zz, &graph->mtx, adjproc);
}

int
Zoltan_ZG_Vertex_Info(ZZ* zz, const ZG *const graph,
		      ZOLTAN_ID_PTR *pgid, int **pinput_part) {
  static char *yo = "Zoltan_ZG_Vertex_Info";
  int ierr = ZOLTAN_OK;
  int nX;
  ZOLTAN_ID_PTR gid = NULL;
  int *input_part = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  nX = graph->mtx.mtx.nY;
  gid = ZOLTAN_MALLOC_GID_ARRAY(zz, nX);
  if (nX && gid == NULL) MEMORY_ERROR;
  (*pgid) = gid;
  if (pinput_part != NULL) {
    input_part = (int*) ZOLTAN_MALLOC(nX * sizeof(int));
    if (nX && input_part == NULL) MEMORY_ERROR;
    (*pinput_part) = input_part;
  }
  else {
    input_part = NULL;
  }

  ierr = Zoltan_DD_Find(graph->mtx.mtx.ddX, (ZOLTAN_ID_PTR)graph->mtx.mtx.yGNO,
			gid, NULL, input_part, nX, NULL);
 End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


  /* This function may work on any distribution of the bipartite graph */
int
Zoltan_ZG_Register(ZZ* zz, ZG* graph, int* properties)
{
  static char *yo = "ZG_Register";
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
			       1, graph->mtx.mtx.globalY/zz->Num_Proc, 0);
      /* Hope a linear assignment will help a little */
      Zoltan_DD_Set_Neighbor_Hash_Fn1(graph->mtx.mtx.ddY, graph->mtx.mtx.globalX/zz->Num_Proc);
    }
    dd = graph->mtx.mtx.ddY;
  }
  /* Make our new numbering public */
  ierr = Zoltan_DD_Update (dd, GID, NULL, NULL, props, size);

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
  if (graph->bipartite && graph->fixObj)
    dd = graph->mtx.mtx.ddX;

  return Zoltan_DD_Find(dd, GID, NULL, NULL, properties, GID_length, NULL);
}

void
Zoltan_ZG_Free(ZG *graph){
  /* TODO : free the communicators properly */

  if (graph->bipartite)
    ZOLTAN_FREE(&graph->fixed_vertices);


/*   Zoltan_Matrix_Free(zz, &graph->mtx.mtx); */

  Zoltan_Matrix2d_Free(&graph->mtx);
  ZOLTAN_FREE(&graph->mtx.comm);
}





#ifdef __cplusplus
}
#endif
