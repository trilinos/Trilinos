/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_parmetis_part_id = "$Id$";
#endif

/* Interface routine to ParMETIS. */

#include "lb_const.h"
#include "all_allo_const.h"
#include "parmetis_const.h"

void LB_ParMETIS_part(LB *lb, int *num_imp, LB_GID** imp_gids, 
                       LB_LID** imp_lids, int **imp_procs)
{
  int i, j, ierr, size, flag, wgtflag, offset, hi;
  int num_obj, nedges, sum_edges, max_edges, cross_edges, edgecut;
  int nsend, nrecv;
  int options[4], *destproc;
  idxtype *xadj, *adjncy, *vtxdist, *vwgt, *adjwgt, *part;
  struct LB_vtx_list  *proc_list, *ptr, new;
  struct LB_hash_node **hashtab, *hash_nodes;
  LB_LID *local_ids;
  LB_GID *global_ids;
  char *sendbuf, *recvbuf;
  struct Comm_Obj *plan;
  MPI_Request *request;
  MPI_Status *status;

  /* Get options */

  /* Set up ParMETIS data structures */
  num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr){
    /* Return error code */
  }

  global_ids = (LB_GID *) LB_array_alloc(__FILE__, __LINE__, 
                1, num_obj, sizeof(LB_GID) );
  local_ids = (LB_LID *) LB_array_alloc(__FILE__, __LINE__, 
                1, num_obj, sizeof(LB_LID) );
  if (!global_ids || !local_ids){
    /* Return not-enough-memory error code */
  }
  if (lb->Get_Obj_List != NULL){
    lb->Get_Obj_List(lb->Get_Obj_List_Data, global_ids, local_ids, &ierr);
  } else {
    /* Use Iterator functions to loop through list */
  }
  if (ierr){
    /* Return error code */
  }

  sum_edges = 0;
  max_edges = 0;
  for (i=0; i< num_obj; i++){
    nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
             local_ids[i], &ierr);
    sum_edges += nedges;
    if (nedges>max_edges) max_edges = nedges;
  }

  /* Allocate space for ParMETIS data structs */
  vtxdist= (idxtype *)LB_array_alloc(__FILE__, __LINE__,
            1, nprocs+1, sizeof(idxtype));
  xadj   = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
            1, num_obj+1, sizeof(idxtype));
  adjncy = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
            1, sum_edges+1, sizeof(idxtype));
  if (!vtxdist || !xadj || !adjncy){
    /* Return not-enough-memory error code */
  }

  /* EB TODO: allocate space for vertex and edge weights.
        Note that weights are of type idxtype (=int) in Metis */

  /* Construct ParMETIS graph */
  /* First compute a global dense numbering of the objects/vertices */

  /* Scan over all procs to determine the number range for each proc */
  MPI_Scan (&num_obj, vtxdist, 1, MPI_IDXTYPE, MPI_SUM, lb->Communicator);
  MPI_Allgather (&vtxdist[0], 1, MPI_IDXTYPE, 
                 &vtxdist[1], 1, MPI_IDXTYPE, lb->Communicator);
  vtxdist[0] = 0;

  /* Construct local hash table */
  hash_nodes = (struct LB_hash_node *)LB_array_alloc(__FILE__, __LINE__,
      1, num_obj, sizeof(struct LB_hash_node));
  hashtab = (struct LB_hash_node **) LB_array_alloc(__FILE__, __LINE__,
      1, num_obj, sizeof(struct LB_hash_node *) );
  if ((!hash_nodes) || (!hashtab)){
    /* Return not-enough-memory error code */
  }
  
  for (i=0; i< num_obj; i++){
    hashtab[i] = NULL;
  }

  for (i=0; i< num_obj; i++){
    /* insert hashf(global_ids[i]) into hash table */
    j = LB_hashf(global_ids[i], num_obj);
    hash_nodes[i]->next = hashtab[j];
    hashtab[j] = hash_nodes[i];
    /* assign it the global number vtxdist[proc]+i */
    hashtab[j]->gid = global_ids[i];
    hashtab[j]->gno = vtxdist[lb->Proc]+i;
  }

  
  /* Construct edge list */
  nbors_global = (LB_GID *)LB_array_alloc(__FILE__, __LINE__,
            1, max_edges, sizeof(LB_GID));
  nbors_proc = (int *)LB_array_alloc(__FILE__, __LINE__,
            1, max_edges, sizeof(int));
  if (using_ewgts)
    nbors_ewgts = (double *)LB_array_alloc(__FILE__, __LINE__,
            1, max_edges, sizeof(double));
  proc_list = (struct LB_vtx_list **) LB_array_alloc(__FILE__, __LINE__,
                1, lb->Num_Proc, sizeof(struct LB_vtx_list *) );
  if (!nbors_global || !nbors_proc || (using_ewgts && !nbors_ewgts)
      || (!proc_list)){
    /* Return not-enough-memory error code */
  }

  /* Initialize pointers */
  for (i=0; i< lb->Num_Proc; i++){
    proc_list[i] = NULL;
  }

  sum_edges = 0;
  max_edges = 0;
  adjptr = adjcny;
  xadj[0] = 0;

  for (i=0; i< num_obj; i++){
    nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
             local_ids[i], &ierr);
    xadj[i+1] = xadj[i] + nedges;
    lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
        nbors_global, nbors_proc, nbors_ewgts, &ierr);
    if (ierr){
      /* Return error code */
    }
    /* Separate inter-processor edges from the local ones */
    cross_edges = 0;
    nsend = 0;
    for (j=0; j<nedges; j++){
      if (nbors_proc[j] == lb->proc){
        /* local edge */
        *adjptr++ = LB_hash_lookup(hashtab, num_obj, nbors_global[j]);
      } else {
        /* Inter-processor edge */
        cross_edges++;
        /* Add it to list */
        new = (struct LB_vtx_list *) LB_SMALLOC (sizeof((struct LB_vtx_list));
        new->next = proc_list[nbors_proc[j]];
        if (new->next == NULL){
          new->length = 1;
          nsend ++;
        } else {
          new->length = (new->next)->length+1;
        }
        new->my_gid = global_ids[i];
        new->my_gno = LB_hash_lookup(hashtab, num_obj, global_ids[i]);
        new->nbor_gid = nbors_global[j];
        new->adj = adjptr++;
        proc_list[nbors_proc[j]] = new;
      }
    }
    sum_edges += cross_edges;
    if (cross_edges>max_edges) max_edges = cross_edges;
  }

  /* Exchange info between processors to resolve global_number */
  /* Allocate buffers */
  size = sizeof(LB_GID) + sizeof(int);
  sendbuf = (char *) LB_SMALLOC(max_edges*size);
  recvbuf = (char *) LB_SMALLOC(sum_edges*size);
  request = (MPI_Request *) LB_SMALLOC(nsend*sizeof(MPI_Request));
  status  = (MPI_Status *) LB_SMALLOC(nsend*sizeof(MPI_Status));
  /* Issue the recvs */
  offset = 0;
  j = 0;
  for (i=0; i<lb->Num_Proc; i++){
    if (proc_list[i] != NULL){
      MPI_Irecv(&recvbuf[offset], proc_list[i]->length *size,
        MPI_BYTE, i, 1, lb->Communicator, &request[j]);
      offset += proc_list[i]->length * size);
      j++;
    }
  }
  /* Barrier */
  MPI_Barrier(lb->Communicator);
  /* Issue the sends */
  for (i=0; i<lb->Num_Proc; i++){
    if (proc_list[i] != NULL){
      /* Pack data to send to proc i */
      offset = 0;
      for (ptr = proc_list[i]; ptr != NULL; ptr = ptr->next){
        memcpy(&sendbuf[offset], &(ptr->my_gid), sizeof(LB_GID)); 
        offset += sizeof(LB_GID);
        memcpy(&sendbuf[offset], &(ptr->my_gno), sizeof(int)); 
        offset += sizeof(int);
      }
      MPI_Rsend(sendbuf, proc_list[i]->length *size,
        MPI_BYTE, i, 1, lb->Communicator);
    }
  }
  /* Wait for all */
  MPI_Waitall(nsend, request, status);

  /* Unpack data into ParMETIS struct */
  /* Resolve off-proc global_ids. */
  for (i=0; i<lb->Num_Proc; i++){
    offset = 0;
    for (ptr = proc_list[i]; ptr != NULL; ptr = ptr->next){
      /* Look for matching global_id in recvbuf */
      flag = 0;
      hi = offset + (proc_list[i]->length)*size;
      for (j=offset; j<hi; j += size){
        if (LB_EQ_GID((LB_GID *)&recvbuf[j], ptr->gid)){
          /* Found match. Amend adjcny array. */
          flag = 1;
          j += sizeof(LB_GID);
          ptr->adj = *((int *)&recvbuf[j]);
          break;
        }
      }
      if (!flag){
         /* Error in algorithm! */
      }
      offset = hi;
    }
  }
  /* Free space for tmp data structures */
  LB_safe_free((void **) &hashtab);
  LB_safe_free((void **) &hash_nodes);
  LB_safe_free((void **) &sendbuf);
  LB_safe_free((void **) &recvbuf);
  LB_safe_free((void **) &request);
  LB_safe_free((void **) &status);

  /* Call ParMETIS */
  options[0] = 0; /* No options for now */
  wgtflag = 0; /* weights not impl. yet */
  part = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
          1, num_obj, sizeof(idxtype));
  
  /* Only PartKway is supported for now */
  ParMETIS_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 0,
    lb->Num_Proc, options, &edgecut, part, lb->Communicator);

  /* Construct send/recv data from ParMETIS output */
  destproc = (int *)LB_array_alloc(__FILE__, __LINE__,
             1, num_obj, sizeof(int));
  j = 0;
  for (i=0; i<num_obj; i++){
    if (part[i] != lb->Proc){
      /* Need to send this data to another proc */
      destproc[j++] = part[i];
    }
  }
  nsend = j;
  size = sizeof(LB_GID) + sizeof(LB_LID) + sizeof(int);
  sendbuf = (char *)LB_array_alloc(__FILE__, __LINE__,
             1, nsend, size);
  j = 0;
  for (i=0; i<num_obj; i++){
    if (part[i] != lb->Proc){
      /* Need to send this data to another proc */
      memcpy(&sendbuf[j], &global_ids[i], sizeof(LB_GID));
      j += sizeof(LB_GID);
      memcpy(&sendbuf[j], &local_ids[i], sizeof(LB_LID));
      j += sizeof(LB_LID);
      memcpy(&sendbuf[j], &part[i], sizeof(int));
      j += sizeof(int);
    }
  }

  /* Create a communication plan */
  plan = LB_comm_create(nsend, destproc, lb->Communicator, &nrecv);

  /* Allocate enough space for receive buffer */
  recvbuf = (char *)LB_array_alloc(__FILE__, __LINE__,
             1, nrecv, size);

  /* Do the communication */
  LB_comm_do(plan, sendbuf, nrecv, recvbuf);

  /* Unpack received data into proper places */
  *num_imp = nrecv;
  *imp_gids = (LB_GID *)LB_array_alloc(__FILE__, __LINE__,
               1, nrecv, sizeof(LB_GID));
  *imp_lids = (LB_LID *)LB_array_alloc(__FILE__, __LINE__,
               1, nrecv, sizeof(LB_LID));
  *imp_procs = (int *)LB_array_alloc(__FILE__, __LINE__,
               1, nrecv, sizeof(int));
  j = 0;
  for (i=0; i<nrecv; i++){
      memcpy(&imp_gids[i], &recvbuf[j], sizeof(LB_GID));
      j += sizeof(LB_GID);
      memcpy(&imp_lids[i], &recvbuf[j], sizeof(LB_LID));
      j += sizeof(LB_LID);
      memcpy(&imp_procs[i], &recvbuf[j], sizeof(int));
      j += sizeof(int);
  }


  /* Free space */
  LB_safe_free((void **) &sendbuf);
  LB_safe_free((void **) &recvbuf);
  LB_safe_free((void **) &vtxdist);
  LB_safe_free((void **) &xadj);
  LB_safe_free((void **) &adjncy);
  LB_safe_free((void **) &vwgt);
  LB_safe_free((void **) &adjwgt);
  LB_safe_free((void **) &part);
  LB_safe_free((void **) &local_ids);
  LB_safe_free((void **) &global_ids);

}

