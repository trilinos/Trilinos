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

/* #define LB_DEBUG */   /* turn on debug print statements? */

#include "lb_const.h"
#include "lb_util_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_const.h"
#include "params_const.h"
#include "timer_const.h"
#include <string.h>

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define MAX_OPTIONS             4

/* Misc. local constants */
#define CHUNKSIZE 10  /* Number of list nodes to allocate in one chunk. */

/* Macro to free all allocated memory */
#define FREE_MY_MEMORY \
  { int i; \
  LB_FREE(&vtxdist); LB_FREE(&xadj); LB_FREE(&adjncy); \
  LB_FREE(&vwgt); LB_FREE(&adjwgt); LB_FREE(&part); \
  LB_FREE(&float_vwgt); LB_FREE(&xyz); \
  LB_FREE(&sendbuf); LB_FREE(&recvbuf); LB_FREE(&request); \
  LB_FREE(&status); LB_FREE(&hash_nodes); LB_FREE(&hashtab); \
  LB_FREE(&nbors_proc); LB_FREE(&nbors_global); \
  LB_FREE(&local_ids); LB_FREE(&global_ids); \
  if (proc_nodes !=NULL) \
    for (i=0; proc_nodes[i] != NULL; i++) \
      LB_FREE(&(proc_nodes[i])); \
  LB_FREE(&proc_nodes); LB_FREE(&proc_list); \
  }


/******** Interface routine between Zoltan and ParMetis. ************/

int LB_ParMetis(
  LB *lb,             /* load balancing object */
  int *num_imp,       /* number of objects to be imported */
  LB_GID **imp_gids,  /* global ids of objects to be imported */
  LB_LID **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,    /* list of processors to import from */
  int *num_exp,       /* number of objects to be exported */
  LB_GID **exp_gids,  /* global ids of objects to be exported */
  LB_LID **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs     /* list of processors to export to */
)
{
#ifdef LB_NO_PARMETIS
  fprintf(stderr, "Error: ParMetis requested but not compiled into library.\n");
  return LB_FATAL;

#else /* !LB_NO_PARMETIS */
  int i, j, ierr, size, flag, offset, hi, row, ndims;
  int num_obj, nedges, sum_edges, max_edges, edgecut;
  int options[MAX_OPTIONS], *nbors_proc;
  int nsend, ewgt_dim, vwgt_dim, wgtflag, numflag;
  int get_graph_data, get_geom_data, get_times; 
  idxtype *vtxdist, *xadj, *adjncy, *adjptr, *vwgt, *adjwgt, *part;
  float  max_wgt, *float_vwgt, *xyz;
  double geom_vec[6];
  struct LB_vtx_list  **proc_list, **proc_nodes, *ptr;
  struct LB_hash_node **hashtab, *hash_nodes;
  LB_LID *local_ids;
  LB_GID *global_ids, *nbors_global;
  char *sendbuf, *recvbuf, alg[MAX_PARAM_STRING_LEN+1];
  MPI_Request *request;
  MPI_Status *status;
  double times[5];

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
  printf("[%1d] Debug: Entering LB_ParMetis()\n", lb->Proc);
#endif

  /* Set default return values (in case of early exit) */
  *num_exp = 0;
  *num_imp = -1; /* No import data */

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  nbors_proc = NULL;
  vtxdist = xadj = adjncy = adjptr = vwgt = adjwgt = part = NULL;
  float_vwgt = xyz = NULL;
  proc_list = proc_nodes = NULL;
  ptr = NULL;
  hashtab = NULL;
  hash_nodes = NULL;
  sendbuf = recvbuf = NULL;
  request = NULL;
  status = NULL;

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

  /* Start timer as soon as options have been set */
  get_times = (options[OPTION_DBGLVL]>0);
  if (get_times){
    MPI_Barrier(lb->Communicator);
    times[0] = MPI_Wtime();
  }

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
    printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
    FREE_MY_MEMORY;
    return (ierr);
  }
  
#ifdef LB_DEBUG
    printf("[%1d] Debug: num_obj =%d\n", lb->Proc, num_obj);
#endif
  
  vtxdist = (idxtype *)LB_MALLOC((lb->Num_Proc+1)* sizeof(idxtype));
  global_ids = (LB_GID *) LB_MALLOC(num_obj * sizeof(LB_GID) );
  local_ids = (LB_LID *) LB_MALLOC(num_obj * sizeof(LB_LID) );
  if (vwgt_dim)
    float_vwgt = (float *)LB_MALLOC(vwgt_dim*num_obj * sizeof(float));
  else {
    float_vwgt = NULL;
    vwgt = NULL;
  }
  if (!vtxdist || !global_ids || !local_ids || (vwgt_dim && !float_vwgt)){
    /* Not enough memory */
    printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
    FREE_MY_MEMORY;
    return LB_MEMERR;
  }
  LB_Get_Obj_List(lb, global_ids, local_ids, vwgt_dim, float_vwgt, &ierr);
  if (ierr){
    /* Return error */
    printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
    FREE_MY_MEMORY;
    return LB_FATAL;
  }

#ifdef LB_DEBUG
    printf("[%1d] Debug: Global ids = ", lb->Proc);
    for (i99=0; i99<num_obj; i99++) printf("%d ", global_ids[i99]);
    printf("\n");
#endif
  
  /* The vtxdist array is required by all ParMetis routines */
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
  
  if (get_graph_data){

    sum_edges = 0;
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      if (ierr){
        /* Return error */
        printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
        FREE_MY_MEMORY;
        return (ierr);
      }
      sum_edges += nedges;
      if (nedges>max_edges) max_edges = nedges;
    }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Sum_edges = %d\n", lb->Proc, sum_edges);
#endif
  
    /* Allocate space for ParMETIS data structs */
    xadj   = (idxtype *)LB_MALLOC((num_obj+1)* sizeof(idxtype));
    adjncy = (idxtype *)LB_MALLOC(sum_edges * sizeof(idxtype));
    if (ewgt_dim) 
      adjwgt = (idxtype *)LB_MALLOC(ewgt_dim*sum_edges * sizeof(idxtype));
  
    if (!xadj || !adjncy || (ewgt_dim && !adjwgt)){
      /* Not enough memory */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
      return LB_MEMERR;
    }
#ifdef LB_DEBUG
    printf("[%1d] Debug: Successfully allocated ParMetis space\n", lb->Proc);
#endif
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Construct local hash table */
    hash_nodes = (struct LB_hash_node *)LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node));
    hashtab = (struct LB_hash_node **) LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node *) );
    if ((!hash_nodes) || (!hashtab)){
      /* Not enough memory */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
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
  
    
    /* Allocate edge list data */
    nbors_global = (LB_GID *)LB_MALLOC(max_edges * sizeof(LB_GID));
    nbors_proc = (int *)LB_MALLOC(max_edges * sizeof(int));
    proc_list = (struct LB_vtx_list **) LB_MALLOC(lb->Num_Proc *
      sizeof(struct LB_vtx_list *) );
    proc_nodes = (struct LB_vtx_list **) LB_MALLOC((sum_edges/CHUNKSIZE+1) *
      sizeof(struct LB_vtx_list *) );
    if ((!nbors_global) || (!nbors_proc) || (!proc_list) || (!proc_nodes)){
      /* Not enough memory */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
      return LB_MEMERR;
    }
  
    /* proc_list[i] will contain a pointer to a linked list of data 
     * to send to processor i. proc_nodes is a list of pointers
     * to chunks of list elements we have MALLOC'ed. The motivation
     * for this design is to reduce the number of calls to
     * MALLOC.  We could have reduced it into a single MALLOC 
     * if we were willing to call the query function get_edge_list 
     * twice for every object, but that is probably
     * even more expensive.
     */
    for (i=0; i< lb->Num_Proc; i++){
      proc_list[i] = NULL;
    }
    for (i=0; i< sum_edges/CHUNKSIZE+1; i++){
      proc_nodes[i] = NULL;
    }
  
    row = 0;        /* Index to current row in proc_nodes */
    offset = 0;     /* Index of next available node in proc_nodes[row] */
    nsend = 0;      /* # procs to send to */
    adjptr = adjncy;
    xadj[0] = 0;
  
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      xadj[i+1] = xadj[i] + nedges;
      lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
          nbors_global, nbors_proc, ewgt_dim, adjwgt, &ierr);
      if (ierr){
        /* Return error */
        printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
        FREE_MY_MEMORY;
        return (ierr);
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
          /* Check whether we need to allocate more space */
          if (offset == 0){
            proc_nodes[row] = (struct LB_vtx_list *) LB_MALLOC(CHUNKSIZE*sizeof(struct LB_vtx_list));
#ifdef LB_DEBUG
          printf("[%1d] Debug: Allocating more list space, row = %d\n", 
          lb->Proc, row);
#endif
            if (!proc_nodes[row]){
              /* Not enough memory */
              printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
              FREE_MY_MEMORY;
              return LB_MEMERR;
            }
          } 
          /* Set ptr to point to next available node */
          ptr = &(proc_nodes[row][offset]);
          ptr->next = proc_list[nbors_proc[j]];
          if (ptr->next == NULL){
            ptr->length = 1;
            nsend++;
#ifdef LB_DEBUG
          printf("[%1d] Debug: creating new list, nsend =%d\n", lb->Proc, nsend);
#endif
          } else {
            ptr->length = ptr->next->length + 1;
#ifdef LB_DEBUG
          printf("[%1d] Debug: appending to proc_list[%1d], new length =%d\n", 
          lb->Proc, nbors_proc[j], ptr->length);
#endif
          }
          ptr->my_gid = global_ids[i];
          ptr->my_gno = LB_hash_lookup(hashtab, global_ids[i], num_obj);
          ptr->nbor_gid = nbors_global[j];
          ptr->adj = adjptr++;
          proc_list[nbors_proc[j]] = ptr;

          offset++;
          if (offset == CHUNKSIZE){
            offset = 0;
            row++;
          }
        }
      }
    }
    /* Sanity check */
    if (((int)adjptr - (int)adjncy)/sizeof(int) != xadj[num_obj]){
      printf("Warning: Internal error in LB_ParMetis, incorrect pointer\n");
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
#ifdef LB_DEBUG
   printf("[%1d] Debug: Proc_list[%1d] has length %d, first elem = (%d,%d)\n", 
   lb->Proc, i, proc_list[i]->length, proc_list[i]->my_gid, proc_list[i]->my_gno);
#endif
        }
      }
  
      /* Allocate buffers */
      size = sizeof(LB_GID) + sizeof(int);
      sendbuf = (char *) LB_MALLOC(max_edges*size);
      recvbuf = (char *) LB_MALLOC(sum_edges*size);
      request = (MPI_Request *) LB_MALLOC(nsend*sizeof(MPI_Request));
      status  = (MPI_Status *) LB_MALLOC(nsend*sizeof(MPI_Status));
      if ((!sendbuf) || (!recvbuf) || (!request) || (!status)){
        /* Not enough space */
        printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
        FREE_MY_MEMORY;
        return LB_MEMERR;
      }

      /* Issue the recvs */
      offset = 0;
      j = 0;
      for (i=0; i<lb->Num_Proc; i++){
        if (proc_list[i] != NULL){
#ifdef LB_DEBUG
          printf("[%1d] Debug: Ready to receive from proc %d\n", lb->Proc, i);
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
      /* Barrier is required when we use Rsend */
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
                /* Found match. */
                flag = 1;
                /* Insert the global number into adjncy vector */
                *(ptr->adj) = *((int *)&recvbuf[j+sizeof(LB_GID)]); 
                /* Advance pointer to next unidentified GID */
                ptr = ptr->next;
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
      LB_FREE(&sendbuf);
      LB_FREE(&recvbuf);
      LB_FREE(&request);
      LB_FREE(&status);
    } /* end if (nsend>0) */
  
    /* Free space for temp data structures */
    for (i=0; proc_nodes[i] != NULL; i++){
      LB_FREE(&(proc_nodes[i]));
    }
    LB_FREE(&proc_nodes);
    LB_FREE(&proc_list);
    LB_FREE(&nbors_global);
    LB_FREE(&nbors_proc);
    LB_FREE(&hash_nodes);
    LB_FREE(&hashtab);
  
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
      LB_FREE(&float_vwgt);
    }

  } /* end get_graph_data */

  if (get_geom_data){
    /* Determine how many dimensions the data have */
    ndims = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
    if (ierr){
      /* Return error */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
      return (ierr);
    }
    /* Allocate space for the geometry data */
    xyz = (float *) LB_MALLOC(ndims*num_obj * sizeof(float));
    if (!xyz){
      /* Not enough space */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
      return LB_MEMERR;
    }
    /* Get the geometry data */
    for (i=0; i<num_obj; i++){
      lb->Get_Geom(lb->Get_Geom_Data, global_ids[i], local_ids[i], 
        geom_vec, &ierr);
      if (ierr) {
        /* Return error code */
        printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
        FREE_MY_MEMORY;
        return (ierr);
      }
      for (j=0; j<ndims; j++)
        xyz[i*ndims+j] = geom_vec[j];
    }
  }

  /* Get ready to call ParMETIS */
  wgtflag = 2*(vwgt_dim>0) + (ewgt_dim>0); /* Multidim wgts not supported yet */
  numflag = 0;
  part = (idxtype *)LB_MALLOC(num_obj * sizeof(idxtype));
  if (!part){
    /* Not enough memory */
    printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
    FREE_MY_MEMORY;
    return LB_MEMERR;
  }
  
  /* Get a time here */
  if (get_times) times[1] = MPI_Wtime();

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
    printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
    FREE_MY_MEMORY;
    return LB_FATAL;
  }

  /* Get a time here */
  if (get_times) times[2] = MPI_Wtime();

#ifdef LB_DEBUG 
    printf("[%1d] Debug: Returned from ParMETIS partitioner with edgecut= %d\n", lb->Proc, edgecut);
#endif

  /* Free weights; they are no longer needed */
  if (vwgt_dim) LB_FREE(&vwgt);
  if (ewgt_dim) LB_FREE(&adjwgt);

  /* Determine number of objects to export */
  nsend = 0;
  for (i=0; i<num_obj; i++)
    if (part[i] != lb->Proc) nsend++;
  (*num_exp) = nsend;

  /* Create export lists */
  if (nsend>0){
    *exp_gids = (LB_GID *)LB_MALLOC(nsend * sizeof(LB_GID));
    *exp_lids = (LB_LID *)LB_MALLOC(nsend * sizeof(LB_LID));
    *exp_procs = (int *)LB_MALLOC(nsend * sizeof(int));
    if (!(*exp_gids) || !(*exp_lids) || !(*exp_procs)){
      /* Not enough memory */
      printf("[%1d] Error on line %d in %s\n", lb->Proc, __LINE__, __FILE__);
      FREE_MY_MEMORY;
      return LB_MEMERR;
    }
    j = 0;
    for (i=0; i<num_obj; i++){
      if (part[i] != lb->Proc){
        (*exp_gids)[j] = global_ids[i];
        (*exp_lids)[j] = local_ids[i];
        (*exp_procs)[j] = part[i];
        j++;
      }
    }
  }

  /* Free space */
  LB_FREE(&part);
  LB_FREE(&local_ids);
  LB_FREE(&global_ids);
  LB_FREE(&vtxdist);
  LB_FREE(&xadj);
  LB_FREE(&adjncy);
  LB_FREE(&xyz);

  /* Get a time here */
  if (get_times) times[3] = MPI_Wtime();

  /* Output timing results if desired */
  if (get_times){
    if (lb->Proc==0) printf("\nZoltan/ParMETIS timing statistics:\n");
    LB_Print_Stats(lb, times[1]-times[0], "ParMETIS  Pre-processing time  ");
    LB_Print_Stats(lb, times[2]-times[1], "ParMETIS  Library time         ");
/*  LB_Print_Stats(lb, times[3]-times[2], "ParMETIS  Post-processing time "); */
    LB_Print_Stats(lb, times[3]-times[0], "ParMETIS  Total time           ");
    if (lb->Proc==0) printf("\n");
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
