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

/* #define LB_DEBUG */  /* turn on debug print statements? */

#include <math.h>
#include <strings.h>
#include "lb_const.h"
#include "lb_util_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_const.h"
#include "params_const.h"

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define MAX_OPTIONS             4


/******** Interface routine between Zoltan and ParMetis. ************/

int LB_ParMetis_Part(
  LB *lb,             /* load balancing object */
  int *num_imp,       /* number of objects to be imported */
  LB_GID **imp_gids,  /* global ids of objects to be imported */
  LB_LID **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs     /* list of processors to import from */
)
{
#ifdef LB_NO_PARMETIS
  fprintf(stderr, "Error: ParMetis requested but not compiled into library.\n");
  return LB_FATAL;

#else /* !LB_NO_PARMETIS */
  int i, j, ierr, size, flag, offset, hi, ndims;
  int num_obj, nedges, sum_edges, max_edges, edgecut;
  int options[MAX_OPTIONS], *destproc, *nbors_proc;
  int nsend, nrecv, ewgt_dim, vwgt_dim, wgtflag, numflag;
  int get_graph_data, get_geom_data;
  idxtype *vtxdist, *xadj, *adjncy, *adjptr, *vwgt, *adjwgt, *part;
  float  max_wgt, *float_vwgt, *xyz;
  double geom_vec[6];
  struct LB_vtx_list  **proc_list, *ptr, *new;
  struct LB_hash_node **hashtab, *hash_nodes;
  LB_LID *local_ids;
  LB_GID *global_ids, *nbors_global;
  char *sendbuf, *recvbuf, alg[MAX_PARAM_STRING_LEN+1];
  struct Comm_Obj *plan;
  MPI_Request *request;
  MPI_Status *status;

  PARAM_VARS parmetis_params[] = {
        { "PARMETIS_METHOD", NULL, "STRING" },
        { "PARMETIS_VWGT_DIM", NULL, "INT" },
        { "PARMETIS_EWGT_DIM", NULL, "INT" },
        { "PARMETIS_COARSE_ALG", NULL, "INT" },
        { "PARMETIS_FOLD", NULL, "INT" },
        { "PARMETIS_OUTPUT_LEVEL", NULL, "INT" },
        { NULL, NULL, NULL } };
  
#ifdef LB_DEBUG
  int i99, *p99;
  printf("[%1d] Debug: Entering ParMetis_Part()\n", lb->Proc);
#endif

  /* Set default return values (in case of early exit) */
  /* Unnecessary because this was done in LB_Balance.
  *num_imp = 0;
  *imp_gids = NULL;
  *imp_lids = NULL;
  *imp_procs = NULL;
  */

  /* Set parameters */
  strcpy(alg, "PartKway");
  vwgt_dim = 0;
  ewgt_dim = 0;
  for (i=0; i<MAX_OPTIONS; i++)
    options[i] = -1;
  parmetis_params[0].ptr = (void *) alg;
  parmetis_params[1].ptr = (void *) &vwgt_dim;
  parmetis_params[2].ptr = (void *) &ewgt_dim;
  parmetis_params[3].ptr = (void *) &(options[OPTION_IPART]);
  parmetis_params[4].ptr = (void *) &(options[OPTION_FOLDF]);
  parmetis_params[5].ptr = (void *) &(options[OPTION_DBGLVL]);

  LB_Assign_Param_Vals(lb->Params, parmetis_params);

  /* Set options[0] to 1 if any of the low level ParMetis options were set,
     or 0 otherwise. This is required by ParMetis. */
  for (i=1; i<MAX_OPTIONS; i++)
    if (options[i]>0) options[0] = 1;
  if (options[0] == -1) 
    options[0] = 0;
  else {
    /* If one option was set, fill in the others. We need to do this
     * because ParMetis requires all or nothing! 
     * The default values below are consistent with ParMetis 2.0 
     */
    if (options[OPTION_IPART] == -1)
      options[OPTION_IPART] = 2; 
    if (options[OPTION_FOLDF] == -1)
      options[OPTION_FOLDF] = 150; 
    if (options[OPTION_DBGLVL] == -1)
      options[OPTION_DBGLVL] = 0; 
  }

#ifdef LB_DEBUG
    printf("[%1d] Debug: alg=%s, vwgt_dim=%d, ewgt_dim=%d\n", lb->Proc, 
      alg, vwgt_dim, ewgt_dim);
    printf("[%1d] Debug: ParMetis options = %d, %d, %d, %d\n", lb->Proc,
      options[0], options[1], options[2], options[3]);
#endif

  /* Most ParMetis methods use only graph data */
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
    printf("[%1d] Debug: num_obj =%d\n", lb->Proc, num_obj);
#endif
  
  global_ids = (LB_GID *) LB_MALLOC(num_obj * sizeof(LB_GID) );
  local_ids = (LB_LID *) LB_MALLOC(num_obj * sizeof(LB_LID) );
  if (vwgt_dim)
    float_vwgt = (float *)LB_MALLOC(vwgt_dim*num_obj * sizeof(float));
  else {
    float_vwgt = NULL;
    vwgt = NULL;
  }

  if (!global_ids || !local_ids || (vwgt_dim && !float_vwgt)){
    /* Not enough memory */
    return LB_MEMERR;
  }
  LB_Get_Obj_List(lb, global_ids, local_ids, vwgt_dim, float_vwgt, &ierr);
  if (ierr){
    /* Return error code ? */
#ifdef LB_DEBUG
    printf("[%1d] Error: LB_Get_Obj_List failed!\n", lb->Proc);
#endif
  }

#ifdef LB_DEBUG
    printf("[%1d] Debug: Global ids = ", lb->Proc);
    for (i99=0; i99<num_obj; i99++) printf("%d ", global_ids[i99]);
    printf("\n");
#endif
  
  if (get_graph_data){

    sum_edges = 0;
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      if (ierr){
      }
      sum_edges += nedges;
      if (nedges>max_edges) max_edges = nedges;
    }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Sum_edges = %d\n", lb->Proc, sum_edges);
#endif
  
    /* Allocate space for ParMETIS data structs */
    vtxdist= (idxtype *)LB_MALLOC((lb->Num_Proc+1)* sizeof(idxtype));
    xadj   = (idxtype *)LB_MALLOC((num_obj+1)* sizeof(idxtype));
    adjncy = (idxtype *)LB_MALLOC(sum_edges * sizeof(idxtype));
    if (ewgt_dim) 
      adjwgt = (idxtype *)LB_MALLOC(ewgt_dim*sum_edges * sizeof(idxtype));
    else
      adjwgt = NULL;
  
    if (!vtxdist || !xadj || !adjncy || (ewgt_dim && !adjwgt)){
      /* Not enough memory */
      return LB_MEMERR;
    }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Successfully allocated ParMetis space\n", lb->Proc);
#endif
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Scan over all procs to determine the number range for each proc */
    MPI_Scan (&num_obj, vtxdist, 1, IDX_DATATYPE, MPI_SUM, lb->Communicator);
    MPI_Allgather (&vtxdist[0], 1, IDX_DATATYPE, 
                   &vtxdist[1], 1, IDX_DATATYPE, lb->Communicator);
    vtxdist[0] = 0;
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: vtxdist = ", lb->Proc);
    for (i99=0; i99<=lb->Num_Proc; i99++)
      printf("%d ", vtxdist[i99]);
    printf("\n");
#endif
  
    /* Construct local hash table */
    hash_nodes = (struct LB_hash_node *)LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node));
    hashtab = (struct LB_hash_node **) LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node *) );
    if ((!hash_nodes) || (!hashtab)){
      /* Not enough memory */
      return LB_MEMERR;
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
    nbors_global = (LB_GID *)LB_MALLOC(max_edges * sizeof(LB_GID));
    nbors_proc = (int *)LB_MALLOC(max_edges * sizeof(int));
    proc_list = (struct LB_vtx_list **) LB_MALLOC(lb->Num_Proc *
      sizeof(struct LB_vtx_list *) );
    if ((!nbors_global) || (!nbors_proc) || (!proc_list)){
      /* Not enough memory */
      return LB_MEMERR;
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
          nbors_global, nbors_proc, ewgt_dim, adjwgt, &ierr);
      if (ierr){
        /* Return error code */
      }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: i=%d, gid=%d, lid=%d, nedges=%d\n", lb->Proc, i, 
      global_ids[i], local_ids[i], nedges);
#endif
  
      /* Separate inter-processor edges from the local ones */
      for (j=0; j<nedges; j++){
  
        if (nbors_proc[j] == lb->Proc){
          /* local edge */
          *adjptr++ = LB_hash_lookup(hashtab, nbors_global[j], num_obj);
        } else {
          /* Inter-processor edge; add it to beginning of the list. */

          /* Note: Here we allocate elements in the list one-by-one,
           * which may be inefficient. Alternatively we could use
           * the num_border_obj query function to obtain
           * the number of cross-edges in advance, but this
           * query function might not be available.
           */
          new = (struct LB_vtx_list *) LB_MALLOC (sizeof(struct LB_vtx_list));
          new->next = proc_list[nbors_proc[j]];
          if (new->next == NULL){
            new->length = 1;
            nsend++;
#ifdef LB_DEBUG
          printf("[%1d] Debug: creating new list, nsend =%d\n", lb->Proc, nsend);
#endif
          } else {
            new->length = new->next->length + 1;
#ifdef LB_DEBUG
          printf("[%1d] Debug: appending to old list, new length =%d\n", 
          lb->Proc, new->length);
#endif
          }
          new->my_gid = global_ids[i];
          new->my_gno = LB_hash_lookup(hashtab, global_ids[i], num_obj);
          new->nbor_gid = nbors_global[j];
          new->adj = adjptr++;
          proc_list[nbors_proc[j]] = new;
        }
      }
    }
    /* Sanity check */
    if (((int)adjptr - (int)adjncy)/sizeof(int) != xadj[num_obj]){
      printf("Warning: Internal error in LB_Parmetis_Part, incorrect pointer\n");
      printf("adjptr-adjncy =%d, #edges =%d\n", ((int)adjptr - (int)adjncy)/sizeof(int), xadj[num_obj]);
    }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: nsend =%d\n", lb->Proc, nsend);
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
      sendbuf = (char *) LB_MALLOC(max_edges*size);
      recvbuf = (char *) LB_MALLOC(sum_edges*size);
      request = (MPI_Request *) LB_MALLOC(nsend*sizeof(MPI_Request));
      status  = (MPI_Status *) LB_MALLOC(nsend*sizeof(MPI_Status));
      /* Issue the recvs */
      offset = 0;
      j = 0;
      for (i=0; i<lb->Num_Proc; i++){
        if (proc_list[i] != NULL){
#ifdef LB_DEBUG
          printf("[%1d] Debug: Receive from proc %d\n", lb->Proc, i);
#endif
          MPI_Irecv(&recvbuf[offset], proc_list[i]->length * size,
            MPI_BYTE, i, 1, lb->Communicator, &request[j]);
          offset += proc_list[i]->length * size;
          j++;
        }
      }
#ifdef LB_DEBUG
      printf("[%1d] Debug: Finished issuing the recvs\n", lb->Proc);
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
            printf("[%1d] Debug: Sending (%d,%d) to proc %d\n", lb->Proc, 
              ptr->my_gid, ptr->my_gno, i);
#endif
            memcpy(&sendbuf[offset], (char *) &(ptr->my_gid), sizeof(LB_GID)); 
            offset += sizeof(LB_GID);
            memcpy(&sendbuf[offset], (char *) &(ptr->my_gno), sizeof(int)); 
            offset += sizeof(int);
          }
          MPI_Rsend(sendbuf, proc_list[i]->length *size,
            MPI_BYTE, i, 1, lb->Communicator);
        }
      }
#ifdef LB_DEBUG
      printf("[%1d] Debug: Finished issuing the sends\n", lb->Proc);
#endif
      /* Wait for all */
      MPI_Waitall(nsend, request, status);
#ifdef LB_DEBUG
      printf("[%1d] Debug: received %d pairs of data from %d procs.\n", 
        lb->Proc, sum_edges, nsend);
      printf("[%1d] Debug: received data: ", lb->Proc);
      p99 = (int *)recvbuf;
      for (i99=0; i99<sum_edges; i99++){
        printf("(%d %d) ", p99[2*i99], p99[2*i99+1]);
      }
      printf("\n");
#endif
    
      /* Unpack data into ParMETIS struct */
      /* Resolve off-proc global_ids. */
      size = sizeof(LB_GID) + sizeof(int);
      hi = 0;
      for (i=0; i<lb->Num_Proc; i++){
        if (proc_list[i] != NULL){
          offset = hi;
          hi = offset + (proc_list[i]->length)*size;
          for (ptr = proc_list[i]; ptr != NULL; ){
#ifdef LB_DEBUG
            printf("[%1d] Debug: Matching data from proc %d, offset=%d, hi=%d\n", 
              lb->Proc, i, offset, hi);
#endif
            /* Look for matching global_id in recvbuf */
            /* The sought gid should be in recvbuf between offset and hi */
            flag = 0;
            for (j=offset; j<hi; j += size){
              if (LB_EQ_GID(*((LB_GID *)&recvbuf[j]), ptr->nbor_gid)){
#ifdef LB_DEBUG
                printf("[%1d] Debug: Matched %d\n", lb->Proc, ptr->nbor_gid);
#endif
                /* Found match. Amend adjncy array. */
                flag = 1;
                /* Insert the global number into adjncy vector */
                *(ptr->adj) = *((int *)&recvbuf[j+sizeof(LB_GID)]); 
                /* Free the node in the list we don't need anymore */
                new = ptr;
                ptr = ptr->next;
                LB_Free((void **)&new);
                break;
              }
            }
            if (!flag){
               /* Error in algorithm! */
               printf("WARNING: Internal error in LB_ParMetis_Part, could not resolve off-proc ID.\n");
            }
          }
        }
      }
      /* Free space for communication data */
      LB_Free((void **) &sendbuf);
      LB_Free((void **) &recvbuf);
      LB_Free((void **) &request);
      LB_Free((void **) &status);
    } /* end if (nsend>0) */
  
    /* Free space for temp data structures */
    LB_Free((void **) &hash_nodes);
    LB_Free((void **) &hashtab);
    LB_Free((void **) &proc_list);
    LB_Free((void **) &nbors_global);
    LB_Free((void **) &nbors_proc);
  
    /* Get vertex weights if needed */
    if (vwgt_dim){
#ifdef LB_DEBUG
      printf("[%1d] Debug: Converting vertex weights...\n", lb->Proc);
#endif
      vwgt = (idxtype *)LB_MALLOC(vwgt_dim*num_obj * sizeof(idxtype));
      max_wgt = 0;
      for (i=0; i<num_obj; i++){
        if (float_vwgt[i]>max_wgt) max_wgt = float_vwgt[i];
      }
      /* Convert weights to integers between 1 and 100 */
      for (i=0; i<vwgt_dim*num_obj; i++){
        vwgt[i] = ceil(float_vwgt[i]*100/max_wgt);
      }
      LB_Free((void **) &float_vwgt);
    }

  } /* end get_graph_data */

  if (get_geom_data){
    /* Determine how many dimensions the data have */
    ndims = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
    /* Allocate space for the geometry data */
    xyz = (float *) LB_MALLOC(ndims*num_obj * sizeof(float));
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
  wgtflag = 2*(vwgt_dim>0) + (ewgt_dim>0); /* Multidim wgts not supported yet */
  numflag = 0;
  part = (idxtype *)LB_MALLOC(num_obj * sizeof(idxtype));
  
  /* Select the desired ParMetis function */
#ifdef LB_DEBUG
    printf("[%1d] Debug: Calling ParMETIS partitioner ...\n", lb->Proc);
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
    return LB_FATAL;
  }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Returned from ParMETIS partitioner with edgecut= %d\n", lb->Proc, edgecut);
#endif

  /* Free weights; they are no longer needed */
  if (vwgt_dim) LB_Free((void **) &vwgt);
  if (ewgt_dim) LB_Free((void **) &adjwgt);

  /* Construct send/recv data from ParMETIS output */
  destproc = (int *)LB_MALLOC(num_obj * sizeof(int));
  nsend = 0;
  for (i=0; i<num_obj; i++){
    if (part[i] != lb->Proc){
      /* Need to send this data to another proc */
      destproc[nsend++] = part[i];
    }
  }

#ifdef LB_DEBUG
  printf("[%1d] Debug: nsend =%d\n", lb->Proc, nsend);
  if (nsend>0){
    printf("[%1d] Debug: destproc = ", lb->Proc);
    for (i99=0; i99<nsend; i99++)
      printf("%d ", destproc[i99]);
    printf("\n");
  }
#endif

  /* Check if any proc needs to communicate */
  MPI_Allreduce(&nsend, &flag, 1, MPI_INT, MPI_SUM, lb->Communicator);

  if (flag>0){

    size = sizeof(LB_GID) + sizeof(LB_LID) + sizeof(int);
    sendbuf = (char *)LB_MALLOC(nsend*size);
#ifdef LB_DEBUG
  printf("[%1d] Debug: copying data to sendbuf.\n", lb->Proc);
#endif
    j = 0;
    for (i=0; i<num_obj; i++){
      if (part[i] != lb->Proc){
        /* Need to send this data to another proc */
        memcpy(&sendbuf[j], (char *) &global_ids[i], sizeof(LB_GID));
        j += sizeof(LB_GID);
        memcpy(&sendbuf[j], (char *) &local_ids[i], sizeof(LB_LID));
        j += sizeof(LB_LID);
        memcpy(&sendbuf[j], (char *) &lb->Proc, sizeof(int));
        j += sizeof(int);
      }
    }
#ifdef LB_DEBUG
  printf("[%1d] Debug: copied data to sendbuf, j=%d\n", lb->Proc, j);
#endif
  
    /* Create a communication plan */
    plan = LB_comm_create(nsend, destproc, lb->Communicator, &nrecv);
  
    /* Allocate enough space for receive buffer */
    recvbuf = (char *)LB_MALLOC(nrecv*size);
  
    /* Do the communication */
#ifdef LB_DEBUG
  printf("[%1d] Debug: calling LB_comm_do.\n", lb->Proc);
#endif
    LB_comm_do(plan, sendbuf, size, recvbuf);
  
    /* Unpack received data into proper places */
    *num_imp = nrecv;
    *imp_gids = (LB_GID *)LB_MALLOC(nrecv * sizeof(LB_GID));
    *imp_lids = (LB_LID *)LB_MALLOC(nrecv * sizeof(LB_LID));
    *imp_procs = (int *)LB_MALLOC(nrecv * sizeof(int));
    if (!(*imp_gids) || !(*imp_lids) || !(*imp_procs)){
      /* Not enough memory */
      return LB_MEMERR;
    }

#ifdef LB_DEBUG
  printf("[%1d] Debug: copying data into output parameters. nrecv =%d\n", 
    lb->Proc, nrecv);
#endif
    j = 0;
    for (i=0; i<nrecv; i++){
        memcpy(&((*imp_gids)[i]), (char *) &recvbuf[j], sizeof(LB_GID));
        j += sizeof(LB_GID);
        memcpy(&((*imp_lids)[i]), (char *) &recvbuf[j], sizeof(LB_LID));
        j += sizeof(LB_LID);
        memcpy(&((*imp_procs)[i]), (char *) &recvbuf[j], sizeof(int));
        j += sizeof(int);
    }
#ifdef LB_DEBUG
    printf("[%1d] Debug: import data (gid,proc) is\n", lb->Proc);
    for (i99=0; i99<nrecv; i99++){
      printf(" (%2d,%2d) ", (*imp_gids)[i99], (*imp_procs)[i99]);
    }
    printf("\n");
#endif

    /* Free buffers */
    LB_Free((void **) &sendbuf);
    LB_Free((void **) &recvbuf);

  } /* end if (flag>0) */

  /* Free space */
  LB_Free((void **) &part);
  LB_Free((void **) &local_ids);
  LB_Free((void **) &global_ids);
  if (get_graph_data){
    LB_Free((void **) &destproc);
    LB_Free((void **) &vtxdist);
    LB_Free((void **) &xadj);
    LB_Free((void **) &adjncy);
  }
  if (get_geom_data){
    LB_Free((void **) &xyz);
  }
#ifdef LB_DEBUG
  printf("[%1d] Debug: exiting ParMetis_Part\n", lb->Proc);
#endif
  return LB_OK;
#endif /* LB_NO_PARMETIS */
}

/*********************************************************************/

int LB_Set_ParMetis_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status, i;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    PARAM_VARS parmetis_params[] = {
        { "PARMETIS_METHOD", NULL, "STRING" },
        { "PARMETIS_VWGT_DIM", NULL, "INT" },
        { "PARMETIS_EWGT_DIM", NULL, "INT" },
        { "PARMETIS_COARSE_ALG", NULL, "INT" },
        { "PARMETIS_FOLD", NULL, "INT" },
        { "PARMETIS_OUTPUT_LEVEL", NULL, "INT" },
        { NULL, NULL, NULL } };
    char *valid_methods[] = {
        "PartKway", "PartGeomKway", "PartGeom", 
        "RepartLDiffusion", "RepartGDiffusion",
        "RepartRemap", "RepartMLRemap",
        "RefineKway",
         NULL };

    status = LB_Check_Param(name, val, parmetis_params, &result, &index);

    if (status == 0){
      /* OK so far, do sanity check of parameter values */

      if (index == 0){
        status = 2;
        for (i=0; valid_methods[i] != NULL; i++){
          if (strcasecmp(val, valid_methods[i]) == 0){
            status = 0; 
            break;
          }
        }
      }
      else{ /* index > 0 */
        if (result.ival < 0)
          status = 2; /* all integer parameters must be non-negative */
      }
    }

    return(status);
}
