#include "lb_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_jostle_const.h"

/*
 * Scatter a ParMetis-style graph to all processors such that each proc
 * gets an even chunk of vertices and corresponding edges. 
 * The entire graph structure is shuffled around; no attempt is made
 * to keep data local. Processor i receives nodes i*n/p through (i+1)*n/p. 
 * This routine is only useful if the distribution is highly imbalanced!
 * Weights are currently ignored.
 *
 * The input graph structure is lost and fresh memory is allocated
 * for the new distributed graph.
 *
 * This routine can also handle geometry data. The graph data may be NULL.
 *
 * The vertex communication plan is returned so it may be used to compute
 * the export lists after partitioning. 
 */

int LB_scatter_graph(
  int     have_graph,		/* do I have graph data, or only the geometry? */
  idxtype **vtxdist,
  idxtype **xadj,
  idxtype **adjncy,
  idxtype **vwgt,
  idxtype **adjwgt,
  float   **xyz,
  int     ndims,		/* # dimensions of xyz geometry data */
  LB      *lb,
  struct Comm_Obj **plan
)
{
  static char *yo = "LB_scatter_graph";
  idxtype *old_vtxdist, *old_xadj, *old_adjncy, *old_vwgt, *old_adjwgt;
  float   *old_xyz;
  int *ptr, *proclist, *proclist2;
  int i, j, num_obj, old_num_obj, num_edges, nrecv;
  int vwgt_dim= lb->Obj_Weight_Dim, ewgt_dim= lb->Comm_Weight_Dim;
  struct Comm_Obj *plan2;

  LB_TRACE_ENTER(lb, yo);

  /* Save pointers to "old" data distribution */
  old_vtxdist = *vtxdist;
  old_xadj = *xadj;
  old_adjncy = *adjncy;
  old_vwgt = *vwgt;
  old_adjwgt = *adjwgt;
  old_xyz = *xyz;

  old_num_obj = old_vtxdist[lb->Proc+1] - old_vtxdist[lb->Proc]; 
  if (lb->Debug_Level >= LB_DEBUG_ALL) 
    printf("[%1d] Debug: Old number of objects = %d\n", lb->Proc, old_num_obj);

  /* Reset all data pointers to NULL for now */
  *vtxdist = *xadj = *adjncy = *vwgt = *adjwgt = NULL;
  *xyz = NULL;

  /* Compute new distribution, *vtxdist */
  (*vtxdist) = (idxtype *)LB_MALLOC((lb->Num_Proc+1)* sizeof(idxtype));
  for (i=0; i<=lb->Num_Proc; i++){
    (*vtxdist)[i] = (i*old_vtxdist[lb->Num_Proc])/lb->Num_Proc;
  }

  /* Convert the xdj array so that it contains the degree of each vertex */
  if (have_graph){
    for (i=0; i<old_num_obj; i++){
      old_xadj[i] = old_xadj[i+1] - old_xadj[i];
    }
  }

  /* Allocate new space for vertex data */
  num_obj = (*vtxdist)[lb->Proc+1] - (*vtxdist)[lb->Proc];
  if (lb->Debug_Level >= LB_DEBUG_ALL) 
    printf("[%1d] Debug: New number of objects = %d\n", lb->Proc, num_obj);
  if (have_graph)
    *xadj = (idxtype *) LB_MALLOC((num_obj+1)*sizeof(idxtype));
  if (vwgt_dim)
    *vwgt = (idxtype *) LB_MALLOC(vwgt_dim*num_obj*sizeof(idxtype));
  if (ndims)
    *xyz = (float *) LB_MALLOC(ndims*num_obj*sizeof(float));

  /* Set up the communication plan for the vertex data */
  proclist = (int *) LB_MALLOC(old_num_obj * sizeof(int));
  /* Let j be the new owner of vertex old_vtxdist[lb->Proc]+i */
  j = 0;
  while (old_vtxdist[lb->Proc] >= (*vtxdist)[j+1]) j++;
  for (i=0; i<old_num_obj; i++){
    if (old_vtxdist[lb->Proc]+i >= (*vtxdist)[j+1]) j++;
    proclist[i] = j;
  }

  LB_Comm_Create( plan, old_num_obj, proclist, lb->Communicator, TAG1, lb->Deterministic, &nrecv);

  if (nrecv != num_obj){
    fprintf(stderr,"Zoltan internal error in %s: Proc %d received %d object but expected %d.\n",
      yo, lb->Proc, nrecv, num_obj);
    /* Free data */
    LB_FREE(&proclist);
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
  }

  /* Do the communication. To save memory, we do not pack all the data into
   * a buffer, but send directly from the old arrays to the new arrays. 
   * We use the vertex communication plan for all the vertex-based arrays
   * and the edge communication plan for all the edge-based arrays.
   */

  if (lb->Debug_Level >= LB_DEBUG_ALL) 
    printf("[%1d] Debug: Starting vertex-based communication.\n", lb->Proc);

  if (have_graph){
    LB_Comm_Do( *plan, TAG2, (char *) old_xadj, sizeof(idxtype), (char *) *xadj);
  }
  if (vwgt_dim){
    LB_Comm_Do( *plan, TAG3, (char *) old_vwgt, vwgt_dim*sizeof(idxtype), (char *) *vwgt);
  }
  if (ndims){
    LB_Comm_Do( *plan, TAG4, (char *) old_xyz, ndims*sizeof(idxtype), (char *) *xyz);
  }
  if (lb->Debug_Level >= LB_DEBUG_ALL) 
    printf("[%1d] Debug: Finished vertex-based communication.\n", lb->Proc);

  if (have_graph){

    /* Rebuild xadj from degrees */
    for (i=1; i<num_obj; i++)
      (*xadj)[i] += (*xadj)[i-1];
    for (i=num_obj; i>0; i--)
      (*xadj)[i] = (*xadj)[i-1];
    (*xadj)[0] = 0;
  
    /* Allocate space for new edge data structures */
    num_edges = (*xadj)[num_obj];
    *adjncy = (idxtype *) LB_MALLOC(num_edges*sizeof(idxtype));
  
    if (ewgt_dim)
      *adjwgt = (idxtype *) LB_MALLOC(ewgt_dim*num_edges*sizeof(idxtype));
  
    /* Set up the communication plan for the edge data. */
    ptr = proclist2 = (int *) LB_MALLOC(old_xadj[old_num_obj] * sizeof(int));
    for (i=0; i<old_num_obj; i++)
      for (j=0; j<old_xadj[i]; j++)
        *ptr++ = proclist[i];
  
    LB_Comm_Create( &plan2, old_xadj[old_num_obj], proclist2, lb->Communicator, TAG1, 
      lb->Deterministic, &nrecv);
  
    if (nrecv != num_edges){
      fprintf(stderr,"Zoltan internal error in %s: Proc %d received %d edges but expected %d.\n",
        yo, lb->Proc, nrecv, num_edges);
      /* Free data */
      LB_FREE(&proclist);
      LB_FREE(&proclist2);
      LB_TRACE_EXIT(lb, yo);
      return LB_FATAL;
    }
  
    if (lb->Debug_Level >= LB_DEBUG_ALL) 
      printf("[%1d] Debug: Starting edge-based communication.\n", lb->Proc);
  
    /* Do the communication. */
    LB_Comm_Do( plan2, TAG2, (char *) old_adjncy, sizeof(idxtype), (char *) *adjncy);
    if (ewgt_dim){
      LB_Comm_Do( plan2, TAG3, (char *) old_adjwgt, ewgt_dim*sizeof(idxtype), (char *) *adjwgt);
    }
  
    if (lb->Debug_Level >= LB_DEBUG_ALL) 
      printf("[%1d] Debug: Finished edge-based communication.\n", lb->Proc);
  
    /* Free the comm. plan for edge data */
    LB_Comm_Destroy(&plan2);

  } /* end of have_graph */

  /* Free data structures */
  LB_FREE(&proclist);
  LB_FREE(&proclist2);
  LB_FREE(&old_vtxdist);
  LB_FREE(&old_xadj);
  LB_FREE(&old_adjncy);
  LB_FREE(&old_vwgt);
  LB_FREE(&old_adjwgt);
  LB_FREE(&old_xyz);

  LB_TRACE_EXIT(lb, yo);
  return LB_OK;
}
