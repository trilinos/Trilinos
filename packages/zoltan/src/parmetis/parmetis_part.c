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

/* Interface routine between Zoltan and ParMetis. */

#define LB_DEBUG /* turn on debug print statements? */

#include <math.h>
#include <strings.h>
#include "lb_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_const.h"
#include "ParMetis/parmetis.h"

void LB_ParMetis_Part(
  LB *lb,             /* load balancing object */
  int *num_imp,       /* number of objects to be imported */
  LB_GID** imp_gids,  /* global ids of objects to be imported */
  LB_LID** imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs     /* list of processors to import from */
)
{
  int i, j, ierr, size, flag, offset, hi, ndims;
  int num_obj, nedges, sum_edges, max_edges, edgecut;
  int options[4], *destproc, *nbors_proc;
  int nsend, nrecv, use_ewgts, use_vwgts, wgtflag, numflag;
  int get_graph_data, get_geom_data;
  idxtype *vtxdist, *xadj, *adjncy, *adjptr, *vwgt, *adjwgt, *part;
  float  max_wgt, *float_vwgt, *xyz;
  double geom_vec[6];
  struct LB_vtx_list  **proc_list, *ptr, *new;
  struct LB_hash_node **hashtab, *hash_nodes;
  LB_LID *local_ids;
  LB_GID *global_ids, *nbors_global;
  char *sendbuf, *recvbuf, *alg;
  struct Comm_Obj *plan;
  MPI_Request *request;
  MPI_Status *status;

#ifdef LB_DEBUG
  int i99, *p99, myproc;
  MPI_Comm_rank(lb->Communicator, &myproc);
  printf("[%1d] Debug: Entering ParMetis_Part()\n", myproc);
#endif

  /* Set default return values (in case of early exit) */
  *num_imp = 0;
  *imp_gids = NULL;
  *imp_lids = NULL;
  *imp_procs = NULL;

  /* Get options. For now, use default values.  */
  alg = "PartKway";
  use_ewgts = 0;
  use_vwgts = 0;
  get_graph_data = 1;
  get_geom_data = 0;

  /* Some algorithms use geometry data */
  if (strcasecmp(alg, "PartGeomKway") == 0){
    get_graph_data = 1;
    get_geom_data = 1;
  }
  else if (strcasecmp(alg, "PartGeom") == 0){
    get_graph_data = 0;
    get_geom_data = 1;
  }

  /* Set up ParMETIS data structures */
  num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr){
    /* Return error code */
  }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: num_obj =%d\n", myproc, num_obj);
#endif
  
  global_ids = (LB_GID *) LB_array_alloc(__FILE__, __LINE__, 
                1, num_obj, sizeof(LB_GID) );
  local_ids = (LB_LID *) LB_array_alloc(__FILE__, __LINE__, 
                1, num_obj, sizeof(LB_LID) );
  if (use_vwgts)
    float_vwgt = (float *)LB_array_alloc(__FILE__, __LINE__,
                  1, num_obj, sizeof(float));
  else
    float_vwgt = NULL;

  if (!global_ids || !local_ids || (use_vwgts && !float_vwgt)){
    /* Return not-enough-memory error code */
  }
  LB_Get_Obj_List(lb, global_ids, local_ids, use_vwgts, float_vwgt, &ierr);
  if (ierr){
    /* Return error code */
  }
  
  if (get_graph_data){

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
              1, lb->Num_Proc+1, sizeof(idxtype));
    xadj   = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
              1, num_obj+1, sizeof(idxtype));
    adjncy = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
              1, sum_edges, sizeof(idxtype));
    if (use_ewgts) 
      adjwgt = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
                1, sum_edges, sizeof(idxtype));
    else
      adjwgt = NULL;
  
    if (!vtxdist || !xadj || !adjncy || (use_ewgts && !adjwgt)){
      /* Return not-enough-memory error code */
    }
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Scan over all procs to determine the number range for each proc */
    MPI_Scan (&num_obj, vtxdist, 1, IDX_DATATYPE, MPI_SUM, lb->Communicator);
    MPI_Allgather (&vtxdist[0], 1, IDX_DATATYPE, 
                   &vtxdist[1], 1, IDX_DATATYPE, lb->Communicator);
    vtxdist[0] = 0;
  
#ifdef LB_DEBUG
    if (myproc<2){
      printf("[%1d] Debug: vtxdist = ", myproc);
      for (i99=0; i99<=lb->Num_Proc; i99++)
        printf("%3d", vtxdist[i99]);
      printf("\n");
    }
#endif
  
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
      hash_nodes[i].gid = global_ids[i];
      hash_nodes[i].gno = vtxdist[lb->Proc]+i;
    }
  
    for (i=0; i< num_obj; i++){
      /* insert hashed elements into hash table */
      j = LB_hashf(global_ids[i], num_obj);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
    }
  
    
    /* Construct edge list */
    nbors_global = (LB_GID *)LB_array_alloc(__FILE__, __LINE__,
              1, max_edges, sizeof(LB_GID));
    nbors_proc = (int *)LB_array_alloc(__FILE__, __LINE__,
              1, max_edges, sizeof(int));
    proc_list = (struct LB_vtx_list **) LB_array_alloc(__FILE__, __LINE__,
                  1, lb->Num_Proc, sizeof(struct LB_vtx_list *) );
    if ((!nbors_global) || (!nbors_proc) || (!proc_list)){
      /* Return not-enough-memory error code */
    }
  
    /* Initialize pointers */
    for (i=0; i< lb->Num_Proc; i++){
      proc_list[i] = NULL;
    }
  
    sum_edges = 0;
    max_edges = 0;
    nsend = 0;
    adjptr = adjncy;
    xadj[0] = 0;
  
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      xadj[i+1] = xadj[i] + nedges;
      lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
          nbors_global, nbors_proc, use_ewgts, adjwgt, &ierr);
      if (ierr){
        /* Return error code */
      }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: local no =%d, nedges =%d\n", myproc, i, nedges);
#endif
  
      /* Separate inter-processor edges from the local ones */
      for (j=0; j<nedges; j++){
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: edge %d is to %d on proc %d\n", myproc,
      j, nbors_global[j], nbors_proc[j]);
#endif
  
        if (nbors_proc[j] == lb->Proc){
          /* local edge */
          *adjptr++ = LB_hash_lookup(hashtab, num_obj, nbors_global[j]);
#ifdef LB_DEBUG
          adjptr--; /* Go back a step */
          printf("[%1d] Debug: found local edge to %d\n", myproc, *adjptr++);
#endif
        } else {
          /* Inter-processor edge */
          /* Add it to beginning of the list */
          new = (struct LB_vtx_list *) LB_SMALLOC (sizeof(struct LB_vtx_list));
          new->next = proc_list[nbors_proc[j]];
          if (new->next == NULL){
            new->length = 1;
            nsend++;
#ifdef LB_DEBUG
          printf("[%1d] Debug: creating new list, nsend =%d\n", myproc, nsend);
#endif
          } else {
            new->length = new->next->length + 1;
#ifdef LB_DEBUG
          printf("[%1d] Debug: appending to old list, new length =%d\n", 
          myproc, new->length);
#endif
          }
          new->my_gid = global_ids[i];
          new->my_gno = LB_hash_lookup(hashtab, num_obj, global_ids[i]);
          new->nbor_gid = nbors_global[j];
          new->adj = adjptr++;
          proc_list[nbors_proc[j]] = new;
#ifdef LB_DEBUG
          printf("[%1d] Debug: found edge to proc %d\n", myproc, nbors_proc[j]);
          printf("[%1d] Debug: new edge = (%d,%d,%d)\n", 
            myproc, new->my_gid, new->my_gno, new->nbor_gid);
#endif
        }
      }
    }
    /* Self test */
    if (((int)adjptr - (int)adjncy)/sizeof(int) != xadj[num_obj]){
      printf("Warning: Internal error in LB_Parmetis_Part, incorrect pointer\n");
      printf("adjptr-adjncy =%d, #edges =%d\n", ((int)adjptr - (int)adjncy)/sizeof(int), xadj[num_obj]);
    }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: nsend =%d\n", myproc, nsend);
#endif
  
    if (nsend>0){
      /* Exchange info between processors to resolve global_number */
  
      /* Determine required buffer sizes */
      max_edges = 0;
      sum_edges = 0;
      for (i=0; i<lb->Num_Proc; i++){
        if (proc_list[i] != NULL){
           sum_edges += proc_list[i]->length;
           if (proc_list[i]->length > max_edges) 
             max_edges = proc_list[i]->length;
        }
      }
  
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
#ifdef LB_DEBUG
          printf("[%1d] Debug: Receive from proc %d\n", myproc, i);
#endif
          MPI_Irecv(&recvbuf[offset], proc_list[i]->length * size,
            MPI_BYTE, i, 1, lb->Communicator, &request[j]);
          offset += proc_list[i]->length * size;
          j++;
        }
      }
#ifdef LB_DEBUG
      printf("[%1d] Debug: Finished issuing the recvs\n", myproc);
#endif
      /* Barrier */
      MPI_Barrier(lb->Communicator);
      /* Issue the sends */
      for (i=0; i<lb->Num_Proc; i++){
        if (proc_list[i] != NULL){
          /* Pack data to send to proc i */
          offset = 0;
          for (ptr = proc_list[i]; ptr != NULL; ptr = ptr->next){
#ifdef LB_DEBUG
            printf("[%1d] Debug: Sending (%d,%d) to proc %d\n", myproc, 
              ptr->my_gid, ptr->my_gno, i);
#endif
            memcpy(&sendbuf[offset], &(ptr->my_gid), sizeof(LB_GID)); 
            offset += sizeof(LB_GID);
            memcpy(&sendbuf[offset], &(ptr->my_gno), sizeof(int)); 
            offset += sizeof(int);
          }
          MPI_Rsend(sendbuf, proc_list[i]->length *size,
            MPI_BYTE, i, 1, lb->Communicator);
        }
      }
#ifdef LB_DEBUG
      printf("[%1d] Debug: Finished issuing the sends\n", myproc);
#endif
      /* Wait for all */
      MPI_Waitall(nsend, request, status);
#ifdef LB_DEBUG
      printf("[%1d] Debug: received %d pairs of data from %d procs.\n", 
        myproc, sum_edges, nsend);
      printf("[%1d] Debug: received data: ", myproc);
      p99 = (int *)recvbuf;
      for (i99=0; i99<sum_edges; i99++){
        printf("%d %d ", p99[2*i99], p99[2*i99+1]);
      }
      printf("\n");
#endif
    
      /* Unpack data into ParMETIS struct */
      /* Resolve off-proc global_ids. */
      size = sizeof(LB_GID) + sizeof(int);
      hi = 0;
      for (i=0; i<lb->Num_Proc; i++){
        offset = hi;
        for (ptr = proc_list[i]; ptr != NULL; ){
#ifdef LB_DEBUG
          printf("[%1d] Debug: Matching data from proc %d, offset=%d\n", 
            myproc, i, offset);
#endif
          /* Look for matching global_id in recvbuf */
          /* The sought gid should be in recvbuf between offset and hi */
          flag = 0;
          hi = offset + (proc_list[i]->length)*size;
          for (j=offset; j<hi; j += size){
#ifdef LB_DEBUG
            printf("[%1d] Debug: Comparing GIDs %d and %d\n", myproc, 
            *((LB_GID *)&recvbuf[j]), ptr->nbor_gid);
#endif
            if (LB_EQ_GID(*((LB_GID *)&recvbuf[j]), ptr->nbor_gid)){
#ifdef LB_DEBUG
              printf("[%1d] Debug: Match!\n", myproc);
#endif
              /* Found match. Amend adjncy array. */
              flag = 1;
              /* Insert the global number into adjncy vector */
              *(ptr->adj) = *((int *)&recvbuf[j+sizeof(LB_GID)]); 
              /* Free the node in the list we don't need anymore */
              new = ptr;
              ptr = ptr->next;
              LB_safe_free((void **)&new);
              break;
            }
          }
          if (!flag){
             /* Error in algorithm! */
             printf("WARNING: Internal error in LB_ParMetis_Part, could not resolve off-proc ID.\n");
          }
        }
      }
      /* Free space for communication data */
      LB_safe_free((void **) &sendbuf);
      LB_safe_free((void **) &recvbuf);
      LB_safe_free((void **) &request);
      LB_safe_free((void **) &status);
    } /* end if (nsend>0) */
  
    /* Free space for temp data structures */
    LB_safe_free((void **) &hash_nodes);
    LB_safe_free((void **) &hashtab);
    LB_safe_free((void **) &proc_list);
    LB_safe_free((void **) &nbors_global);
    LB_safe_free((void **) &nbors_proc);
  
    /* Get vertex weights if needed */
    if (use_vwgts){
#ifdef LB_DEBUG
      printf("[%1d] Debug: Converting vertex weights...\n", myproc);
#endif
      vwgt = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
                1, num_obj, sizeof(idxtype));
      max_wgt = 0;
      for (i=0; i<num_obj; i++){
        if (float_vwgt[i]>max_wgt) max_wgt = float_vwgt[i];
      }
      /* Convert weights to integers between 1 and 100 */
      for (i=0; i<num_obj; i++){
        vwgt[i] = ceil(float_vwgt[i]*100/max_wgt);
      }
      LB_safe_free((void **) &float_vwgt);
    }

  } /* end get_graph_data */

  if (get_geom_data){
    /* Determine how many dimensions the data have */
    ndims = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
    /* Allocate space for the geometry data */
    xyz = (float *) LB_array_alloc(__FILE__, __LINE__, 
                  1, ndims*num_obj, sizeof(float) );
    if (!xyz){
      /* Not enough space */
    }
    /* Get the geometry data */
    for (i=0; i<num_obj; i++){
      lb->Get_Geom(lb->Get_Geom_Data, global_ids[i], local_ids[i], 
        geom_vec, &ierr);
      if (ierr) {
        /* error */
      }
      for (j=0; j<ndims; j++)
        xyz[i*ndims+j] = geom_vec[j];
    }
  }

  /* Call ParMETIS */
  options[0] = 0; /* No options for now */
  wgtflag = 2*use_vwgts + use_ewgts;
  numflag = 0;
  part = (idxtype *)LB_array_alloc(__FILE__, __LINE__,
          1, num_obj, sizeof(idxtype));
  
  /* Select the desired ParMetis function */
#ifdef LB_DEBUG
    printf("[%1d] Debug: Calling ParMETIS partitioner ...\n", myproc);
#endif
  if (strcasecmp(alg, "PartKway") == 0){
    ParMETIS_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, &(lb->Num_Proc), options, &edgecut, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "PartGeomKway") == 0){
    ParMETIS_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &(lb->Num_Proc), options, &edgecut, 
      part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "PartGeom") == 0){
    ParMETIS_PartGeom (vtxdist, &ndims, xyz, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "RepartLDiffusion") == 0){
    ParMETIS_RepartLDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "RepartGDiffusion") == 0){
    ParMETIS_RepartGDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "RepartRemap") == 0){
    ParMETIS_RepartRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "RepartMLRemap") == 0){
    ParMETIS_RepartMLRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &(lb->Communicator));
  }
  else if (strcasecmp(alg, "RefineKway") == 0){
    ParMETIS_RefineKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &(lb->Communicator));
  }
  else {
    printf("Error: Unknown ParMetis algorithm %s\n", alg);
    /* return DLB_FATAL; */
  }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Returned from ParMETIS partitioner with edgecut= %d\n", myproc, edgecut);
#endif

  /* Free weights; they are no longer needed */
  LB_safe_free((void **) &vwgt);
  LB_safe_free((void **) &adjwgt);

  /* Construct send/recv data from ParMETIS output */
  destproc = (int *)LB_array_alloc(__FILE__, __LINE__,
             1, num_obj, sizeof(int));
  nsend = 0;
  for (i=0; i<num_obj; i++){
    if (part[i] != lb->Proc){
      /* Need to send this data to another proc */
      destproc[nsend++] = part[i];
    }
  }

#ifdef LB_DEBUG
  printf("[%1d] Debug: nsend =%d\n", myproc, nsend);
  if (nsend>0){
    printf("[%1d] Debug: destproc = ", myproc);
    for (i99=0; i99<nsend; i99++)
      printf("%d ", destproc[i99]);
    printf("\n");
  }
#endif

  /* Check if any proc needs to communicate */
  MPI_Allreduce(&nsend, &flag, 1, MPI_INT, MPI_SUM, lb->Communicator);

  if (flag>0){

    size = sizeof(LB_GID) + sizeof(LB_LID) + sizeof(int);
    sendbuf = (char *)LB_SMALLOC(nsend*size);
#ifdef LB_DEBUG
  printf("[%1d] Debug: copying data to sendbuf.\n", myproc);
#endif
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
#ifdef LB_DEBUG
  printf("[%1d] Debug: copied data to sendbuf, j=%d\n", myproc, j);
#endif
  
    /* Create a communication plan */
    plan = LB_comm_create(nsend, destproc, lb->Communicator, &nrecv);
  
    /* Allocate enough space for receive buffer */
    recvbuf = (char *)LB_SMALLOC(nrecv*size);
  
    /* Do the communication */
#ifdef LB_DEBUG
  printf("[%1d] Debug: calling LB_comm_do.\n", myproc);
#endif
    LB_comm_do(plan, sendbuf, size, recvbuf);
  
    /* Unpack received data into proper places */
    *num_imp = nrecv;
    *imp_gids = (LB_GID *)LB_SMALLOC(nrecv * sizeof(LB_GID));
    *imp_lids = (LB_LID *)LB_SMALLOC(nrecv * sizeof(LB_LID));
    *imp_procs = (int *)LB_SMALLOC(nrecv * sizeof(int));
    if (!(*imp_gids) || !(*imp_lids) || !(*imp_procs)){
      /* Not enough memory */
    }

#ifdef LB_DEBUG
  printf("[%1d] Debug: copying data into output parameters. nrecv =%d\n", 
    myproc, nrecv);
#endif
    j = 0;
    for (i=0; i<nrecv; i++){
        memcpy(&((*imp_gids)[i]), &recvbuf[j], sizeof(LB_GID));
        j += sizeof(LB_GID);
        memcpy(&((*imp_lids)[i]), &recvbuf[j], sizeof(LB_LID));
        j += sizeof(LB_LID);
        memcpy(&((*imp_procs)[i]), &recvbuf[j], sizeof(int));
        j += sizeof(int);
    }
    /* Free buffers */
    LB_safe_free((void **) &sendbuf);
    LB_safe_free((void **) &recvbuf);

  } /* end if (flag>0) */

#ifdef LB_DEBUG
  printf("[%1d] Debug: almost finished with Parmetis_Part, freeing space...\n", myproc);
#endif
  /* Free space */
  LB_safe_free((void **) &part);
  LB_safe_free((void **) &local_ids);
  LB_safe_free((void **) &global_ids);
  if (get_graph_data){
    LB_safe_free((void **) &destproc);
    LB_safe_free((void **) &vtxdist);
    LB_safe_free((void **) &xadj);
    LB_safe_free((void **) &adjncy);
  }
  if (get_geom_data){
    LB_safe_free((void **) &xyz);
  }

}


