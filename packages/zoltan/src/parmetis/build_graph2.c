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

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_dd.h"
#include "all_allo_const.h"
#include "third_library_const.h"
#include "third_library_tools.h"
#include "params_const.h"

/* Parameters for how to build the graph */
static PARAM_VARS Graph_params[] = {
	{ "GRAPH_UNSYMMETRIC", NULL, "INT", 0 },
	{ "GRAPH_FAST_BUILD", NULL, "STRING", 0 },
	{ NULL, NULL, NULL, 0 } };

extern int Zoltan_Build_Graph_Erik(ZZ *zz, int graph_type, int check_graph,
       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
       int obj_wgt_dim, int edge_wgt_dim,
       indextype **vtxdist, indextype **xadj, indextype **adjncy, float **ewgts,
       int **adjproc);

int Zoltan_Build_Graph_Cedric(ZZ *zz, int graph_type, int check_graph,
       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
       int obj_wgt_dim, int edge_wgt_dim,
       indextype **vtxdist, indextype **xadj, indextype **adjncy, float **ewgts,
       int **adjproc, int *success);

int Zoltan_Build_Graph_NoComm(
    ZZ *zz, int graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int edge_wgt_dim,
    indextype **vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc, int *success);



static int
Zoltan_Sort_Edges(ZZ *zz, int num_obj, int edge_wgt_dim, const int * const xadj,
		  int *adjncy, float *ewgts, int* adjproc);

int Zoltan_Build_Graph(
    ZZ *zz, int *graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int* edge_wgt_dim,
    indextype **vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc)
{
  int success =0;
  int ierr=ZOLTAN_OK;

  int unsymmetric = 0;
  int symmetrize = 0;
  char try_fast_str[MAX_PARAM_STRING_LEN+1];

  strcpy(try_fast_str, "NO");  /* default */

  Zoltan_Bind_Param(Graph_params, "GRAPH_UNSYMMETRIC", (void *) &unsymmetric);
  Zoltan_Bind_Param(Graph_params, "GRAPH_FAST_BUILD", (void *) try_fast_str);
  Zoltan_Assign_Param_Vals(zz->Params, Graph_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);

  (*graph_type) &= (~(1<<FORCE_FAST));
  if (!strcasecmp(try_fast_str, "NO")) {
    (*graph_type) &= (~(1<<TRY_FAST));
  }
  else {
    if (!strcasecmp(try_fast_str, "FORCE"))
      (*graph_type) |= (1 << FORCE_FAST);
    (*graph_type) |= (1<< TRY_FAST);
  }

  if (unsymmetric)
    (*graph_type) |= (1 << UNSYMMETRIC);
  else
    (*graph_type) &= (~(1<<UNSYMMETRIC));

  if (symmetrize)
    (*graph_type) |= (1 << SYMMETRIZE);
  else
    (*graph_type) &= (~(1<<SYMMETRIZE));

  if ((*graph_type)&(1<<TRY_FAST)) /* Try to construct the graph without any communication */
    ierr = Zoltan_Build_Graph_NoComm(zz, *graph_type, check_graph, num_obj,
				     global_ids, local_ids,
				     obj_wgt_dim, *edge_wgt_dim,
				     vtxdist, xadj, adjncy, ewgts, adjproc, &success);

  if ((*graph_type)&(1<<UNSYMMETRIC)) {
    if (!success) /* If the graph is not built yet */
      ierr = Zoltan_Build_Graph_Cedric(zz, *graph_type, check_graph, num_obj,
				       global_ids, local_ids,
				       obj_wgt_dim, *edge_wgt_dim,
				       vtxdist, xadj, adjncy, ewgts, adjproc, &success);

    if ((success) && ((*graph_type)&(1<<SYMMETRIZE))) /* Symmetrize the graph, A+At way */
      ierr = Zoltan_Symmetrize_Graph(zz, *graph_type, check_graph, num_obj,
				     global_ids, local_ids,
				     obj_wgt_dim, edge_wgt_dim,
				     vtxdist, xadj, adjncy, ewgts, adjproc);
  }
  else if (!success)
      ierr = Zoltan_Build_Graph_Erik(zz, *graph_type, check_graph, num_obj,
				     global_ids, local_ids,
				     obj_wgt_dim, *edge_wgt_dim,
				     vtxdist, xadj, adjncy, ewgts, adjproc);

  /* CC: in order to compare to the new implementation I need to sort edges */
  if ((ierr != ZOLTAN_OK) || (ierr != ZOLTAN_WARN)) {
    ierr = Zoltan_Sort_Edges(zz, num_obj, *edge_wgt_dim, *xadj, *adjncy, *ewgts, *adjproc);
  }


  return (ierr);
}


/*
 * Build a graph in ParMetis format from Zoltan query functions.
 * The dynamic arrays are allocated in this function
 * and should be freed by the calling function after use.
 * Geometric methods may use this function to only
 * compute vtxdist (and nothing else) by setting graph_type&NO_GRAPH.
 * The global and local ids are assumed to be available.
 * Also, the vertex weights are not computed here since
 * they are typically obtained together with the gids.
 *
 * If graph_type&GLOBAL_GRAPH, construct the global graph,
 * otherwise (LOCAL_GRAPH) construct a local subgraph on each proc
 * (and discard all inter-proc edges).
 *
 * This version will be slower than Zoltan_Build_Graph, but it
 * can handle unsymmetric graphs.
 */

int Zoltan_Build_Graph_Cedric(
    ZZ *zz, int graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int edge_wgt_dim,
    indextype **vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc, int *success)
{
  /* Local variables */
  int  num_edges, max_edges;
  int *nbors_proc = NULL;
  int *edges_per_obj = NULL;
  int nself;
  int ierr = ZOLTAN_OK, i, i99, tmp;
  int max_obj;
  float *tmp_ewgts;
  ZOLTAN_ID_PTR lid;
  ZOLTAN_ID_PTR nbors_global;
  ZOLTAN_ID_PTR proc_list_nbor; /* Global IDs of neighbors of proc_list
				       entries.  This array is separate from
				       proc_list to prevent individual mallocs
				       for nbor global IDs.   */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  char msg[256];


  char *yo = "Zoltan_Build_Graph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set pointers to NULL */
  *vtxdist = *xadj = *adjncy = NULL;
  *adjproc = NULL;
  *ewgts = tmp_ewgts = NULL;
  nbors_global = proc_list_nbor = lid = NULL;
  nbors_proc = NULL;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: num_obj =%d\n", zz->Proc, num_obj);

  if (!(graph_type&NO_GRAPH)){
      if ((zz->Get_Num_Edges == NULL && zz->Get_Num_Edges_Multi == NULL) ||
	  (zz->Get_Edge_List == NULL && zz->Get_Edge_List_Multi == NULL))
	ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
	      "A graph query function is not registered.\n");
  }

  *vtxdist = (indextype *)ZOLTAN_MALLOC((zz->Num_Proc+1)* sizeof(indextype));
  if (num_obj>0){
    if (!(*vtxdist)){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }

  /* Construct *vtxdist[i] = the number of objects on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&num_obj, *vtxdist, 1, IDX_DATATYPE, MPI_SUM, zz->Communicator);
  /* Gather data from all procs */
  MPI_Allgather (&((*vtxdist)[0]), 1, IDX_DATATYPE,
		 &((*vtxdist)[1]), 1, IDX_DATATYPE, zz->Communicator);
  (*vtxdist)[0] = 0;

  max_obj = num_obj;
  for (i = 0 ; i < zz->Num_Proc ; ++i) {
    max_obj = MAX(max_obj, ((*vtxdist)[i+1] - (*vtxdist)[i]));
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: vtxdist = ", zz->Proc);
    for (i99=0; i99<=zz->Num_Proc; i99++)
      printf("%d ", (*vtxdist)[i99]);
    printf("\n");
  }

  if (zz->Debug_Level){
     if ((zz->Proc ==0) && ((*vtxdist)[zz->Num_Proc]==0))
	ZOLTAN_PRINT_WARN(zz->Proc, yo, "No objects to balance.");
  }

  if (!(graph_type&(1<<NO_GRAPH))){
    int vertex;
    int *vertexid;
    int offset;
    ZOLTAN_ID_PTR edgetab;
    struct Zoltan_DD_Struct *dd;

    /* Get edge data */
    Zoltan_Get_Num_Edges_Per_Obj(zz, num_obj, global_ids, local_ids,
				 &edges_per_obj, &max_edges, &num_edges);
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: num_edges = %d\n", zz->Proc, num_edges);

    /* Allocate space for ParMETIS data structs */
    *xadj   = (indextype *)ZOLTAN_MALLOC((num_obj+1) * sizeof(indextype));
    *adjncy = (indextype *)ZOLTAN_MALLOC(num_edges * sizeof(indextype));
    *adjproc = (int *)ZOLTAN_MALLOC(num_edges * sizeof(int));

    if (!(*xadj) || (num_edges && (!(*adjncy) || !(*adjproc)))){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Successfully allocated ParMetis space\n", zz->Proc);

    ierr = Zoltan_DD_Create (&dd, zz->Communicator, num_gid_entries, 0, sizeof(indextype)/sizeof(int), num_obj, 0);
    if (ierr != ZOLTAN_OK) {
      ZOLTAN_PARMETIS_ERROR(ierr, "Cannot create data directory.");
    }

    /* Optimize hash function */
    if (num_gid_entries == 1)   /* Assume gids are consecutive numbers */
    {
      int *proc;

      proc = (int*)ZOLTAN_MALLOC(zz->Num_Proc*sizeof(int));
      for (i= 0 ; i < zz->Num_Proc ; ++i) {
	proc[i] = i;
      }

      Zoltan_DD_Set_Neighbor_Hash_Fn2 (dd, proc, (*vtxdist), (*vtxdist)+1, zz->Num_Proc);
      ZOLTAN_FREE(&proc);
    }


    vertexid = (indextype*)ZOLTAN_MALLOC(num_obj*sizeof(indextype));
    for (vertex = 0 ; vertex < num_obj ; ++ vertex) {
      if (IS_GLOBAL_GRAPH(graph_type))
        vertexid[vertex] = (*vtxdist)[zz->Proc]+vertex; /* Make this a global number */
      else /* graph_type is LOCAL_GRAPH */
        vertexid[vertex] = vertex;                      /* Make this a local number */
    }

    /* Make our new numbering public */
    Zoltan_DD_Update (dd, global_ids, NULL, (ZOLTAN_ID_PTR) vertexid, NULL, num_obj);
    ZOLTAN_FREE(&vertexid);

    edgetab = ZOLTAN_MALLOC_GID_ARRAY(zz, num_edges);
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, num_edges);
    nbors_proc = (int *)ZOLTAN_MALLOC(num_edges * sizeof(int));
    if (edge_wgt_dim && max_edges){
      tmp_ewgts = (float *)ZOLTAN_MALLOC(edge_wgt_dim * num_edges * sizeof(float));
      *ewgts = (float *)ZOLTAN_MALLOC(edge_wgt_dim * num_edges * sizeof(float));
    }

    if (((!nbors_global) || (!nbors_proc) ||
                       (edge_wgt_dim && !(*ewgts)) ||
                       (edge_wgt_dim && !tmp_ewgts)) || (!edgetab)){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }


    if (zz->Get_Edge_List_Multi) {
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
				num_gid_entries, num_lid_entries,
				num_obj, global_ids,
				local_ids,
				edges_per_obj,
				nbors_global, nbors_proc, edge_wgt_dim,
				tmp_ewgts, &ierr);
    }
    else {
      int edge;
      for (vertex = 0, edge = 0 ; vertex < num_obj ; ++vertex) {
	zz->Get_Edge_List(zz->Get_Edge_List_Data,
                          num_gid_entries, num_lid_entries,
                          &(global_ids[vertex*num_gid_entries]), &(local_ids[vertex*num_lid_entries]),
                          nbors_global+edge*num_gid_entries, nbors_proc+edge, edge_wgt_dim,
                          tmp_ewgts+edge*edge_wgt_dim, &ierr);
	edge += edges_per_obj[vertex];
      }
    }

    nself = 0; (*xadj)[0] = 0;
    for (vertex = 0, offset=0 ; vertex < num_obj ; ++vertex) {
      int edge;
      int trueedge;

      trueedge = (*xadj)[vertex];

      for (edge = offset ; edge < offset+edges_per_obj[vertex] ; ++edge) {
	if (IS_LOCAL_GRAPH(graph_type) && (nbors_proc[edge] != zz->Proc))
	  continue;                  /* Keep only local edges */
	if (ZOLTAN_EQ_GID(zz, &(global_ids[vertex*num_gid_entries]),
			  &(nbors_global[edge*num_gid_entries]))) {
	  nself ++;                  /* Remove self edges */
	  continue;
	}
	ZOLTAN_SET_GID (zz, &(edgetab[trueedge*num_gid_entries]),
			&(nbors_global[edge*num_gid_entries]));
	if (edge_wgt_dim > 0)
	  *(ewgts)[trueedge] = tmp_ewgts[edge];
	trueedge ++;
      }

      (*xadj)[vertex+1] = trueedge;
      offset+=edges_per_obj[vertex];
    }

    ZOLTAN_FREE(&nbors_global);
    ZOLTAN_FREE(&nbors_proc);
    ZOLTAN_FREE(&edges_per_obj);

    Zoltan_DD_Find (dd, edgetab, NULL, (ZOLTAN_ID_PTR)((*adjncy)), NULL,
		    num_edges, (*adjproc));

    Zoltan_DD_Destroy (&dd);

    ZOLTAN_FREE(&edgetab);

    /* Warn if we removed any self-edges */
    if (check_graph >= 1){
      if (nself>0) ierr = ZOLTAN_WARN;
      MPI_Reduce(&nself, &tmp, 1, MPI_INT, MPI_SUM, 0, zz->Communicator);
      if ((zz->Proc==0) && (tmp>0) && (zz->Debug_Level>0)){
	  sprintf(msg, "Found and removed %d self edges in the graph.\n",
		  tmp);
	  ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      }
    }

    *success = 1;

  }

 End:

  ZOLTAN_FREE(&tmp_ewgts);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


/*
 * Build a graph in ParMetis format from Zoltan query functions.
 * The dynamic arrays are allocated in this function
 * and should be freed by the calling function after use.
 * Geometric methods may use this function to only
 * compute vtxdist (and nothing else) by setting graph_type = NO_GRAPH.
 * The global and local ids are assumed to be available.
 * Also, the vertex weights are not computed here since
 * they are typically obtained together with the gids.
 *
 * If graph_type = GLOBAL_GRAPH, construct the global graph,
 * otherwise (LOCAL_GRAPH) construct a local subgraph on each proc
 * (and discard all inter-proc edges).
 *
 * This version will be slower than Zoltan_Build_Graph, but it
 * can handle unsymmetric graphs.
 */

int Zoltan_Build_Graph_NoComm(
    ZZ *zz, int graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int edge_wgt_dim,
    indextype **vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc, int* success)
{
  /* Local variables */
  int  num_edges, max_edges;
  int *edges_per_obj = NULL;
  int nself;
  int ierr, i, tmp;
  ZOLTAN_ID_PTR lid;
  ZOLTAN_ID_PTR nbors_global;
  ZOLTAN_ID_PTR proc_list_nbor; /* Global IDs of neighbors of proc_list
				       entries.  This array is separate from
				       proc_list to prevent individual mallocs
				       for nbor global IDs.   */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  char msg[256];
  int offset = 0;


  char *yo = "Zoltan_Build_Graph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (num_gid_entries != 1) {
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_WARN,
			  "Unhandled size of GID\n");
  }

  /* Set pointers to NULL */
  *vtxdist = *xadj = *adjncy = NULL;
  *adjproc = NULL;
  *ewgts =  NULL;
  nbors_global = proc_list_nbor = lid = NULL;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: num_obj =%d\n", zz->Proc, num_obj);

  if (!(graph_type&(1<<NO_GRAPH))){
      if ((zz->Get_Num_Edges == NULL && zz->Get_Num_Edges_Multi == NULL) ||
	  (zz->Get_Edge_List == NULL && zz->Get_Edge_List_Multi == NULL))
	ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
	      "A graph query function is not registered.\n");
  }

  *vtxdist = (indextype *)ZOLTAN_MALLOC((zz->Num_Proc+1)* sizeof(indextype));
  if (num_obj>0){
    if (!(*vtxdist)){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }

  /* Construct *vtxdist[i] = the number of objects on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&num_obj, *vtxdist, 1, IDX_DATATYPE, MPI_SUM, zz->Communicator);
  /* Gather data from all procs */
  MPI_Allgather (&((*vtxdist)[0]), 1, IDX_DATATYPE,
		 &((*vtxdist)[1]), 1, IDX_DATATYPE, zz->Communicator);
  (*vtxdist)[0] = 0;

  if (num_obj > 0) {
    offset = (*vtxdist)[zz->Proc] - global_ids[0] ;
  }
  else
    offset = 0;

  if (!(graph_type&(1<<FORCE_FAST))) {
    int flag;
    int sndoffset[3];
    int rcvoffset[3];

    flag = 0;
    for (i = 0; i < num_obj ; ++i) {
      if (offset != (*vtxdist)[zz->Proc] + i - global_ids[i])
	flag =1;
    }
    sndoffset[0] = flag;
    sndoffset[1] = offset;
    sndoffset[2] = -offset;

    MPI_Allreduce (sndoffset, rcvoffset, 3, MPI_INT, MPI_MAX, zz->Communicator);

    if ((rcvoffset[0] != 0) || (rcvoffset[1] != -rcvoffset[2])) {
      ZOLTAN_FREE(vtxdist);
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_WARN,
			    "Cannot quickly handle GID\n");
    }
    offset = rcvoffset[1];
  }
/*   else { */
/*     ierr = ZOLTAN_WARN; */
/*     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Assume FAST graph construction\n"); */
/*   } */

  if (!(graph_type&(1<<NO_GRAPH))){
    int vertex;
    int begin;

    /* Get edge data */
    Zoltan_Get_Num_Edges_Per_Obj(zz, num_obj, global_ids, local_ids,
				 &edges_per_obj, &max_edges, &num_edges);
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: num_edges = %d\n", zz->Proc, num_edges);

    /* Allocate space for ParMETIS data structs */
    /* Trick, use 2 times the same array : do clean in place */
    *xadj   = (indextype *)ZOLTAN_MALLOC((num_obj+1) * sizeof(indextype));
    *adjncy = (indextype *)ZOLTAN_MALLOC(num_edges * sizeof(indextype));
    *adjproc = (int *)ZOLTAN_MALLOC(num_edges * sizeof(int));

    if (!(*xadj) || (num_edges && (!(*adjncy) || !(*adjproc)))){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Successfully allocated ParMetis space\n", zz->Proc);


    if (edge_wgt_dim && max_edges){
      *ewgts = (float *)ZOLTAN_MALLOC(edge_wgt_dim * num_edges * sizeof(float));
    }

    if (edge_wgt_dim && !(*ewgts)) {
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }


    if (zz->Get_Edge_List_Multi) {
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
			      num_gid_entries, num_lid_entries,
			      num_obj, global_ids,
			      local_ids,
			      edges_per_obj,
			      (ZOLTAN_ID_PTR)(*adjncy), (*adjproc), edge_wgt_dim,
			      (*ewgts), &ierr);
    }
    else {
      int edge;
      for (vertex = 0, edge = 0 ; vertex < num_obj ; ++vertex) {
	zz->Get_Edge_List(zz->Get_Edge_List_Data,
                          num_gid_entries, num_lid_entries,
                          &(global_ids[vertex*num_gid_entries]), &(local_ids[vertex*num_lid_entries]),
                          (ZOLTAN_ID_PTR)((*adjncy)+edge), ((*adjproc)+edge), edge_wgt_dim,
                          (*ewgts)+edge, &ierr);
	edge += edges_per_obj[vertex];
      }
    }

    nself = 0; (*xadj)[0] = 0;
    for (vertex = 0, begin=0 ; vertex < num_obj ; ++vertex) {
      int edge;
      int trueedge;

      trueedge = (*xadj)[vertex];

      for (edge = begin ; edge < begin+edges_per_obj[vertex] ; ++edge) {
	if (IS_LOCAL_GRAPH(graph_type) && ((*adjproc)[edge] != zz->Proc))
	  continue;                  /* Keep only local edges */
	if (vertex + (*vtxdist)[zz->Proc] == (*adjncy)[edge] + offset) {
	  nself ++;                  /* Remove self edges */
	  continue;
	}
	(*adjncy)[trueedge] = (*adjncy)[edge] + offset;
	if (edge_wgt_dim > 0)
	  *(ewgts)[trueedge] = (*ewgts)[edge];
	trueedge ++;
      }

      (*xadj)[vertex+1] = trueedge;
      begin+=edges_per_obj[vertex];
    }

    ZOLTAN_FREE (&edges_per_obj);

    /* Perhaps helps to free some unnecessary memory */
    (*adjncy) = (indextype*) ZOLTAN_REALLOC((*adjncy), (*xadj)[num_obj] * sizeof(indextype));
    (*adjproc) = (int*) ZOLTAN_REALLOC((*adjproc), (*xadj)[num_obj] * sizeof(int));
    if (edge_wgt_dim > 0)
      (*ewgts) = (float*) ZOLTAN_REALLOC((*ewgts), (*xadj)[num_obj] * sizeof(float));

    /* Warn if we removed any self-edges */
    if (check_graph >= 1){
      if (nself>0) ierr = ZOLTAN_WARN;
      MPI_Reduce(&nself, &tmp, 1, MPI_INT, MPI_SUM, 0, zz->Communicator);
      if ((zz->Proc==0) && (tmp>0) && (zz->Debug_Level>0)){
	  sprintf(msg, "Found and removed %d self edges in the graph.\n",
		  tmp);
	  ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      }
    }

  }

  *success = 1;


 End:

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

static int
Zoltan_Sort_Edges(ZZ *zz, int num_obj, int edge_wgt_dim, const int * const xadj,
		  int *adjncy, float *ewgts, int* adjproc)
{
  int i;
  int *position;
  float *tmp_ewgts;
  int * tmp_adjproc;

  tmp_ewgts = (float*) ZOLTAN_MALLOC(xadj[num_obj]*sizeof(float)*edge_wgt_dim);
  tmp_adjproc = (int*) ZOLTAN_MALLOC(xadj[num_obj]*sizeof(int));
  position = (int*)ZOLTAN_MALLOC(xadj[num_obj]*sizeof(int));

  memcpy(tmp_ewgts, ewgts, xadj[num_obj]*sizeof(float)*edge_wgt_dim);
  memcpy(tmp_adjproc, adjproc, xadj[num_obj]*sizeof(int));

  for (i = 0 ; i<num_obj ; ++i) {
    int degree = xadj[i+1] - xadj[i];
    int j;

    for (j = 0 ; j < degree ; ++j)
      position[j] = xadj[i]+j;
    Zoltan_Comm_Sort_Ints( adjncy + xadj[i], position , degree);
    for (j = 0 ; j < degree ; ++j) {
      adjproc[position[j]] = tmp_adjproc[xadj[i]+j];
      memcpy (ewgts+position[j]*edge_wgt_dim, tmp_ewgts+(xadj[i]+j)*edge_wgt_dim, edge_wgt_dim*sizeof(float));
    }
  }
  return ZOLTAN_OK;
}


int Zoltan_Build_Graph_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int  index;
  PARAM_UTYPE result;

  return Zoltan_Check_Param(name, val, Graph_params, &result, &index);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
