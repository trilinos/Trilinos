/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "comm_const.h"
#include "parmetis_jostle.h"

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

int Zoltan_Scatter_Graph(
  int     have_graph,		/* do I have graph data, or only the geometry? */
  idxtype **vtxdist,
  idxtype **xadj,
  idxtype **adjncy,
  idxtype **vwgt,
  idxtype **adjwgt,
  float   **xyz,
  int     ndims,		/* # dimensions of xyz geometry data */
  LB      *lb,
  ZOLTAN_COMM_OBJ **plan
)
{
  static char *yo = "Zoltan_Scatter_Graph";
  char     msg[256];
  idxtype *old_vtxdist, *old_xadj, *old_adjncy, *old_vwgt, *old_adjwgt;
  float   *old_xyz;
  int *ptr, *proclist = NULL, *proclist2;
  int i, j, num_obj, old_num_obj, num_edges, nrecv;
  int vwgt_dim= lb->Obj_Weight_Dim, ewgt_dim= lb->Comm_Weight_Dim;
  ZOLTAN_COMM_OBJ *plan2;

  ZOLTAN_TRACE_ENTER(lb, yo);

  /* Save pointers to "old" data distribution */
  old_vtxdist = *vtxdist;
  old_xadj = *xadj;
  old_adjncy = *adjncy;
  old_vwgt = *vwgt;
  old_adjwgt = *adjwgt;
  old_xyz = *xyz;

  old_num_obj = old_vtxdist[lb->Proc+1] - old_vtxdist[lb->Proc]; 
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
    printf("[%1d] Debug: Old number of objects = %d\n", lb->Proc, old_num_obj);

  /* Reset all data pointers to NULL for now */
  *vtxdist = *xadj = *adjncy = *vwgt = *adjwgt = NULL;
  *xyz = NULL;

  /* Compute new distribution, *vtxdist */
  (*vtxdist) = (idxtype *)ZOLTAN_MALLOC((lb->Num_Proc+1)* sizeof(idxtype));
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
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
    printf("[%1d] Debug: New number of objects = %d\n", lb->Proc, num_obj);
  if (have_graph)
    *xadj = (idxtype *) ZOLTAN_MALLOC((num_obj+1)*sizeof(idxtype));
  if (vwgt_dim)
    *vwgt = (idxtype *) ZOLTAN_MALLOC(vwgt_dim*num_obj*sizeof(idxtype));
  if (ndims)
    *xyz = (float *) ZOLTAN_MALLOC(ndims*num_obj*sizeof(float));

  if (old_num_obj > 0) {
    /* Set up the communication plan for the vertex data */
    proclist = (int *) ZOLTAN_MALLOC(old_num_obj * sizeof(int));
    /* Let j be the new owner of vertex old_vtxdist[lb->Proc]+i */
    j = 0;
    while (old_vtxdist[lb->Proc] >= (*vtxdist)[j+1]) j++;
    for (i=0; i<old_num_obj; i++){
      if (old_vtxdist[lb->Proc]+i >= (*vtxdist)[j+1]) j++;
      proclist[i] = j;
    }
  }

  Zoltan_Comm_Create(plan, old_num_obj, proclist, lb->Communicator, TAG1, &nrecv);

  if (nrecv != num_obj){
    sprintf(msg,"Proc %d received %d object but expected %d.",
      lb->Proc, nrecv, num_obj);
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, msg);
    /* Free data */
    ZOLTAN_FREE(&proclist);
    ZOLTAN_TRACE_EXIT(lb, yo);
    return ZOLTAN_FATAL;
  }

  /* Do the communication. To save memory, we do not pack all the data into
   * a buffer, but send directly from the old arrays to the new arrays. 
   * We use the vertex communication plan for all the vertex-based arrays
   * and the edge communication plan for all the edge-based arrays.
   */

  if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
    printf("[%1d] Debug: Starting vertex-based communication.\n", lb->Proc);

  if (have_graph){
    Zoltan_Comm_Do( *plan, TAG2, (char *) old_xadj, sizeof(idxtype), (char *) *xadj);
  }
  if (vwgt_dim){
    Zoltan_Comm_Do( *plan, TAG3, (char *) old_vwgt, vwgt_dim*sizeof(idxtype), (char *) *vwgt);
  }
  if (ndims){
    Zoltan_Comm_Do( *plan, TAG4, (char *) old_xyz, ndims*sizeof(idxtype), (char *) *xyz);
  }
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
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
    *adjncy = (idxtype *) ZOLTAN_MALLOC(num_edges*sizeof(idxtype));
  
    if (ewgt_dim)
      *adjwgt = (idxtype *) ZOLTAN_MALLOC(ewgt_dim*num_edges*sizeof(idxtype));
  
    /* Set up the communication plan for the edge data. */
    ptr = proclist2 = (int *) ZOLTAN_MALLOC(old_xadj[old_num_obj] * sizeof(int));
    for (i=0; i<old_num_obj; i++)
      for (j=0; j<old_xadj[i]; j++)
        *ptr++ = proclist[i];
  
    Zoltan_Comm_Create(&plan2, old_xadj[old_num_obj], proclist2, lb->Communicator, 
                   TAG1, &nrecv);
  
    if (nrecv != num_edges){
      sprintf(msg,"Proc %d received %d edges but expected %d.",
        lb->Proc, nrecv, num_edges);
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, msg);
      /* Free data */
      ZOLTAN_FREE(&proclist);
      ZOLTAN_FREE(&proclist2);
      ZOLTAN_TRACE_EXIT(lb, yo);
      return ZOLTAN_FATAL;
    }
  
    if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
      printf("[%1d] Debug: Starting edge-based communication.\n", lb->Proc);
  
    /* Do the communication. */
    Zoltan_Comm_Do( plan2, TAG2, (char *) old_adjncy, sizeof(idxtype), (char *) *adjncy);
    if (ewgt_dim){
      Zoltan_Comm_Do( plan2, TAG3, (char *) old_adjwgt, ewgt_dim*sizeof(idxtype), (char *) *adjwgt);
    }
  
    if (lb->Debug_Level >= ZOLTAN_DEBUG_ALL) 
      printf("[%1d] Debug: Finished edge-based communication.\n", lb->Proc);
  
    /* Free the comm. plan for edge data */
    Zoltan_Comm_Destroy(&plan2);

  } /* end of have_graph */

  /* Free data structures */
  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&proclist2);
  ZOLTAN_FREE(&old_vtxdist);
  ZOLTAN_FREE(&old_xadj);
  ZOLTAN_FREE(&old_adjncy);
  ZOLTAN_FREE(&old_vwgt);
  ZOLTAN_FREE(&old_adjwgt);
  ZOLTAN_FREE(&old_xyz);

  ZOLTAN_TRACE_EXIT(lb, yo);
  return ZOLTAN_OK;
}
