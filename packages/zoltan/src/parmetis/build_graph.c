/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *   *****************************************************************************/
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
#include "all_allo_const.h"
#include "parmetis_jostle.h"
#include "params_const.h"


/* Prototypes */
static int hash_lookup (ZZ *zz, struct Hash_Node **hashtab, ZOLTAN_ID_PTR key,
                        int n);
static int process_edge_list(ZZ *, int, ZOLTAN_ID_PTR, int, ZOLTAN_ID_PTR, 
  int *, float *, struct Hash_Node **, int, int, int, int, int *, int *, int *, 
  float *, int *, int *, int *, int *, int *, int *, int *, struct Edge_Info **,
  ZOLTAN_ID_PTR *);

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
 * The graph should be symmetric. BUG: Currently, a graceful
 * exit is not guaranteed if the graph is non-symmetric.
 * Future versions should consider automatically
 * symmetrizing the graph (at least as an option).
 */

int Zoltan_Build_Graph(
    ZZ *zz, int graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int edge_wgt_dim,
    idxtype **vtxdist, idxtype **xadj, idxtype **adjncy, 
    float **ewgts)
{
  /* Local variables */
  int  num_edges, cross_edges, max_edges;
  int *nbors_proc, *plist = NULL;
  int *edges_per_obj = NULL;
  int nsend, nrecv, nself, num_border, max_proc_list_len;
  int i, i99, j, jj, k, ierr, packet_size, offset, tmp;
  float *tmp_ewgts;
  char *sendbuf = NULL, *recvbuf = NULL;
  struct Hash_Node **hashtab = NULL, *hash_nodes = NULL;
  struct Edge_Info *ptr;
  struct Edge_Info *proc_list;   /* Edge information; contains global IDs
                                 of objects with off-processor nbors. */
  ZOLTAN_ID_PTR lid;
  ZOLTAN_ID_PTR nbors_global;
  ZOLTAN_ID_PTR proc_list_nbor; /* Global IDs of neighbors of proc_list
                                       entries.  This array is separate from
                                       proc_list to prevent individual mallocs
                                       for nbor global IDs.   */
  MPI_Comm comm = zz->Communicator; 
  ZOLTAN_COMM_OBJ *comm_plan;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_off, lid_off;
  int sum, edge_list_size;
  int gid_size = num_gid_entries * sizeof(ZOLTAN_ID_TYPE);
  char msg[256];


  char *yo = "Zoltan_Build_Graph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set pointers to NULL */
  *vtxdist = *xadj = *adjncy = NULL;
  *ewgts = tmp_ewgts = NULL;
  nbors_global = proc_list_nbor = lid = NULL;
  proc_list = NULL;
  nbors_proc = NULL;
  
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: num_obj =%d\n", zz->Proc, num_obj);

  if (graph_type != NO_GRAPH){
      if ((zz->Get_Num_Edges == NULL && zz->Get_Num_Edges_Multi == NULL) || 
          (zz->Get_Edge_List == NULL && zz->Get_Edge_List_Multi == NULL))
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, 
              "A graph query function is not registered.\n");
  }
  
  *vtxdist = (idxtype *)ZOLTAN_MALLOC((zz->Num_Proc+1)* sizeof(idxtype));
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

  if (graph_type != NO_GRAPH){
    /* Get edge data */
    Zoltan_Get_Num_Edges_Per_Obj(zz, num_obj, global_ids, local_ids, 
                                 &edges_per_obj, &max_edges, &num_edges);
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: num_edges = %d\n", zz->Proc, num_edges);

    /* Allocate space for ParMETIS data structs */
    *xadj   = (idxtype *)ZOLTAN_MALLOC((num_obj+1) * sizeof(idxtype));
    *adjncy = (idxtype *)ZOLTAN_MALLOC(num_edges * sizeof(idxtype));
  
    if (!(*xadj) || (num_edges && !(*adjncy))){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Successfully allocated ParMetis space\n", zz->Proc);
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Construct local hash table */
    hash_nodes = (struct Hash_Node *)ZOLTAN_MALLOC(num_obj *
      sizeof(struct Hash_Node));
    hashtab = (struct Hash_Node **) ZOLTAN_MALLOC(num_obj *
      sizeof(struct Hash_Node *) );
    if (num_obj && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    
    /* Assign consecutive numbers based on the order of the ids */
    for (i=0; i< num_obj; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      if (graph_type == GLOBAL_GRAPH)
        /* Make this a global number */
        hash_nodes[i].gno = (*vtxdist)[zz->Proc]+i;
      else /* graph_type == LOCAL_GRAPH */
        /* Make this a local number */
        hash_nodes[i].gno = i;
    }

    for (i=0; i< num_obj; i++){
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&(global_ids[i*num_gid_entries]), num_gid_entries, (unsigned int)num_obj);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
    }
  
    /* Estimate the number of inter-proc edges. 
     * First estimate the number of border objects 
     * based on a 2d regular grid. (We could alternatively
     * use Get_Num_Border_Obj but this could be expensive
     * and the function may not be reqistered. )
     */
    max_proc_list_len = 0;

    if (graph_type == GLOBAL_GRAPH){
      num_border = 4*sqrt((double) num_obj);
      if (num_border > num_obj) num_border = num_obj;
       
      /* Assume that the edges are approx. evenly distributed among the objs. */
      if (num_obj>0){
         max_proc_list_len = num_edges * num_border / num_obj;
         if (max_proc_list_len < CHUNKSIZE)
            max_proc_list_len = CHUNKSIZE;
      }
    }
    
    /* Allocate edge list data */
    if (zz->Get_Edge_List_Multi) 
      edge_list_size = num_edges;   /* Get all edges at once */
    else
      edge_list_size = max_edges;   /* Get one object's edges at a time */
   
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, edge_list_size);
    nbors_proc = (int *)ZOLTAN_MALLOC(edge_list_size * sizeof(int));
    plist = (int *)ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    if (edge_wgt_dim && max_edges){
      tmp_ewgts = (float *)ZOLTAN_MALLOC(edge_wgt_dim * edge_list_size 
                                                      * sizeof(float));
      *ewgts = (float *)ZOLTAN_MALLOC(edge_wgt_dim * num_edges * sizeof(float));
    }

    if ((edge_list_size && ((!nbors_global) || (!nbors_proc) ||
                       (edge_wgt_dim && !(*ewgts)) || 
                       (edge_wgt_dim && !tmp_ewgts))) || (!plist)){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    for (i=0; i<zz->Num_Proc; i++)
      plist[i] = -1;
  
    if (max_proc_list_len){
      /* Allocate space for processor list */
      while ((proc_list==NULL || proc_list_nbor == NULL)
          && (max_proc_list_len>=CHUNKSIZE)){
        proc_list = (struct Edge_Info *) ZOLTAN_MALLOC(max_proc_list_len *
          sizeof(struct Edge_Info) );
        proc_list_nbor = ZOLTAN_MALLOC_GID_ARRAY(zz, max_proc_list_len);
        if (!proc_list || !proc_list_nbor){
          /* Not enough memory, try shorter list */
          ZOLTAN_FREE(&proc_list);
          ZOLTAN_FREE(&proc_list_nbor);
          if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
            printf("[%1d] Debug: Could not allocate %d list nodes, "
                   "trying %d instead.\n", zz->Proc,
                   max_proc_list_len, max_proc_list_len/2);
          }
          max_proc_list_len /= 2;
        }
      }
      if (!proc_list || !proc_list_nbor){
        /* Not enough memory */
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
      }
    }

    /* proc_list[i] will contain a struct with data to send
     * to another processor.  proc_list_nbor[i*num_gid_entries] contains
     * the global ID of the neighboring object.
     * We don't know yet the total number of inter-proc edges so we may 
     * have to adjust the size of proc_list and proc_list_nbor through
     * REALLOC.  We increase the size by a factor REALLOC_FACTOR each time.
     * The motivation for this design is to reduce the number of calls
     * to REALLOC.  We could have reduced it into a single MALLOC 
     * if we were willing to call the query function get_edge_list 
     * twice for every object, but that is probably
     * even more expensive.
     *
     * Note that proc_list may contain duplicate elements.
     * These are eliminated before the communication phase.
     */
    
    nsend = 0;      /* Number of objects we need to send info about */
    offset = 0;     /* Index of next available node in proc_list */
    ptr = proc_list;
    (*xadj)[0] = 0;
    jj = 0;          /* Index into (*xadj) */
    nself = 0;       /* Number of self-edges in the graph */
    cross_edges = 0; /* Number of edges that cross over to a different proc. */
  
    if (zz->Get_Edge_List_Multi) {
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data, 
                              num_gid_entries, num_lid_entries,
                              num_obj, global_ids, local_ids, edges_per_obj,
                              nbors_global, nbors_proc, edge_wgt_dim, 
                              tmp_ewgts, &ierr);

      sum = 0;
      for (i = 0; i < num_obj; i++) {
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: i=%d, gid=", zz->Proc, i);
          ZOLTAN_PRINT_GID(zz, &(global_ids[i*num_gid_entries]));
          printf("lid=");
          ZOLTAN_PRINT_LID(zz, 
                 (num_lid_entries ? &local_ids[i*num_lid_entries] : NULL));
          printf("nedges=%d\n", edges_per_obj[i]);
        }

        ierr = process_edge_list(zz, i, &(global_ids[i*num_gid_entries]), 
                                 edges_per_obj[i], 
                                 &(nbors_global[sum]), &(nbors_proc[sum]), 
                                 &(tmp_ewgts[sum*edge_wgt_dim]),
                                 hashtab, graph_type, num_gid_entries, 
                                 num_obj, edge_wgt_dim,
                                 *vtxdist, *xadj, *adjncy, *ewgts, plist,
                                 &jj, &nself, &cross_edges, &offset, &nsend,
                                 &max_proc_list_len, &proc_list,
                                 &proc_list_nbor);
        if (ierr) {
          ZOLTAN_PARMETIS_ERROR(ierr, "Error in process_edge_list");
        }
        sum += edges_per_obj[i];
      }
    }
    else {
      for (i=0; i< num_obj; i++){
        gid_off = i * num_gid_entries;
        lid_off = i * num_lid_entries;
        lid = (num_lid_entries ? &(local_ids[lid_off]) : NULL);
  
        zz->Get_Edge_List(zz->Get_Edge_List_Data,
                          num_gid_entries, num_lid_entries,
                          &(global_ids[gid_off]), lid, 
                          nbors_global, nbors_proc, edge_wgt_dim, 
                          tmp_ewgts, &ierr);
        if (ierr){
          /* Return error */
          ZOLTAN_PARMETIS_ERROR(ierr, "Error in Get_Edge_List.");
        }
  
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: i=%d, gid=", zz->Proc, i);
          ZOLTAN_PRINT_GID(zz, &(global_ids[gid_off]));
          printf("lid=");
          ZOLTAN_PRINT_LID(zz, lid);
          printf("nedges=%d\n", edges_per_obj[i]);
        }

        ierr = process_edge_list(zz, i, &(global_ids[gid_off]),
                                 edges_per_obj[i],
                                 nbors_global, nbors_proc, tmp_ewgts,
                                 hashtab, graph_type, num_gid_entries, 
                                 num_obj, edge_wgt_dim,
                                 *vtxdist, *xadj, *adjncy, *ewgts, plist,
                                 &jj, &nself, &cross_edges, &offset, &nsend,
                                 &max_proc_list_len, &proc_list,
                                 &proc_list_nbor);
        if (ierr) {
          ZOLTAN_PARMETIS_ERROR(ierr, "Error in process_edge_list");
        }
      } /* end object loop */
    } /* End zz->Get_Edge_List */
    /* Remark: cross_edges == offset when graph_type==GLOBAL_GRAPH */

    ZOLTAN_FREE(&plist);
    ZOLTAN_FREE(&nbors_global);
    ZOLTAN_FREE(&nbors_proc);
    ZOLTAN_FREE(&edges_per_obj);
    if (tmp_ewgts) ZOLTAN_FREE(&tmp_ewgts);

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

    /* Sanity check for edges. */
    k = (graph_type == GLOBAL_GRAPH ? num_edges : num_edges - cross_edges);
    if ((check_graph >= 1) && ((*xadj)[num_obj] + nself != k)){
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
         printf("[%3d] xadj[num_obj=%1d] = %d, nself = %d, num_edges = %d, cross_edges = %d, k =%d\n", zz->Proc, num_obj, (*xadj)[num_obj],  nself, num_edges, cross_edges, k);
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "Incorrect edge count. Please check your graph query functions.");
    }

    /* Set cross_edges to zero if we threw them out. */
    if (offset==0)
      cross_edges = 0;

    /* Exchange info between processors to resolve global number 
     * for objects that are off-proc.
     */

    /* Allocate send buffer */
    packet_size = gid_size + sizeof(int);
    sendbuf = (char *) ZOLTAN_MALLOC(nsend * packet_size);
    plist = (int *) ZOLTAN_MALLOC(nsend * sizeof(int));

    if (nsend && (!sendbuf || !plist) ){
      /* Not enough space */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }

    /* Pack the data to send */
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: %d messages to send.\n", zz->Proc, nsend);

    offset = 0;
    j = 0;
    for (i=0, ptr=proc_list; i<cross_edges; i++, ptr++){
      if (ptr->nbor_proc >= 0){
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: Sending gid = ", zz->Proc);
          ZOLTAN_PRINT_GID(zz, ptr->my_gid);
          printf(", gno = %d to proc %d\n", ptr->my_gno, ptr->nbor_proc);
        }
        memcpy(&sendbuf[offset], (char *) (ptr->my_gid), gid_size); 
        offset += gid_size;
        memcpy(&sendbuf[offset], (char *) &(ptr->my_gno), sizeof(int)); 
        offset += sizeof(int);
        plist[j++] = ptr->nbor_proc;
      }
    }

    /* Create the communication plan */
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Calling Zoltan_Comm_Create with %d packets to send.\n",
             zz->Proc, nsend);

    ierr = Zoltan_Comm_Create(&comm_plan, nsend, plist, comm, TAG1, &nrecv);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Comm_Create returned error.");
    }
    /* nrecv should be <= cross_edges. */
    /* fatal error if nrecv==0, cross_edges>0. */
    if (nrecv==0 && cross_edges>0){
      /* Return error code and msg  */
      sprintf(msg, "cross_edges=%d but received no edge data from comm plan. "
        "Possible error in the graph query functions.\n", cross_edges);
      ZOLTAN_PARMETIS_ERROR(ierr, msg);
    }

    /* Allocate recv buffer */
    recvbuf = (char *) ZOLTAN_MALLOC(nrecv * packet_size);
    if (nrecv && (!recvbuf) ){
      /* Not enough space */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Ready to receive %d packets.\n", 
        zz->Proc, nrecv);

    /* Do the communication */
    ierr = Zoltan_Comm_Do( comm_plan, TAG2, sendbuf, packet_size, recvbuf);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Comm_Do returned error.");
    }

    /* Destroy the comm. plan */
    Zoltan_Comm_Destroy( &comm_plan);
  
    /* Unpack data into the ParMETIS/Jostle struct.
     * Resolve off-proc global numbers by:
     * 1) Insert global ids from recvbuf into hash table;
     * 2) Look up missing references in proc_list 
     *    using the hash table.
     */

    /* Change hash table to contain only border objects */
    /* that we received from other procs.               */
    hash_nodes = (struct Hash_Node *)ZOLTAN_REALLOC(hash_nodes,
      nrecv * sizeof(struct Hash_Node));
    hashtab = (struct Hash_Node **) ZOLTAN_REALLOC(hashtab,
      nrecv * sizeof(struct Hash_Node *) );
    if (nrecv && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    
    /* Copy data from recvbuf into hash table nodes */
    for (i=0; i< nrecv; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = (ZOLTAN_ID_PTR) &(recvbuf[i*packet_size]);
      hash_nodes[i].gno = *((int *)&recvbuf[i*packet_size+gid_size]);
      /* Do we need to pad for byte alignment? */
    }
  
    /* Insert nodes into hash table */
    for (i=0; i< nrecv; i++){
      j = Zoltan_Hash(hash_nodes[i].gid, num_gid_entries, (unsigned int)nrecv);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
        printf("[%1d] Debug: Hashed GID ", zz->Proc);
        ZOLTAN_PRINT_GID(zz, hash_nodes[i].gid);
        printf(" to %d, gno = %d\n", j, hash_nodes[i].gno);
      }
    }

    for (i=0; i<cross_edges; i++){
      /* Look up unresolved global_ids */
      if ((tmp=hash_lookup(zz, hashtab, &(proc_list_nbor[i*num_gid_entries]), 
                           nrecv)) <0){
        /* Error: Global ID is not in hash table. 
           This only happens if the graph is invalid. */
        sprintf(msg, "Debug info: cross_edges = %d, global id not found: ",
           cross_edges); 
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        ZOLTAN_PRINT_GID(zz, &(proc_list_nbor[i*num_gid_entries]));
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "Invalid or nonsymmetric graph. "
           "Possible error in the graph query functions.\n");
      }
      else{
        /* Insert the global number into adjncy vector */
        *(proc_list[i].adj) = tmp;
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: GID ", zz->Proc);
          ZOLTAN_PRINT_GID(zz, &(proc_list_nbor[i*num_gid_entries]));
          printf(" has global number %d\n", tmp);
        }
      }
      
    }
  }

  /* Successful finish */
  ierr = ZOLTAN_OK; 

End:
  /* Free all local arrays */
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&recvbuf);
  ZOLTAN_FREE(&proc_list);
  ZOLTAN_FREE(&proc_list_nbor);
  ZOLTAN_FREE(&plist);
  ZOLTAN_FREE(&hash_nodes);
  ZOLTAN_FREE(&hashtab);
  ZOLTAN_FREE(&edges_per_obj);

  /* Free tmp_ewgts if they haven't been freed already (error occurred) */
  if (tmp_ewgts) ZOLTAN_FREE(ewgts); 

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

/**************************************************************************/

int Zoltan_Get_Num_Edges_Per_Obj(
  ZZ *zz,
  int num_obj,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int **edges_per_obj,
  int *max_edges,
  int *num_edges
)
{
/* Calls ZOLTAN_NUM_EDGE_FN or ZOLTAN_NUM_EDGE_MULTI_FN to obtain number
 * of edges per object.
 * Returns number of edges per object in array edges_per_obj.
 * Computes max edges per obj and total edges per obj.
 */
char *yo = "Zoltan_Get_Num_Edges_Per_Obj";
int ierr = ZOLTAN_OK;
int i;
int nedges;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
ZOLTAN_ID_PTR lid;

  *max_edges = *num_edges = 0;
  if (num_obj) {

    *edges_per_obj = (int *) ZOLTAN_MALLOC(num_obj * sizeof(int));
    if (*edges_per_obj == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (zz->Get_Num_Edges_Multi) {
      zz->Get_Num_Edges_Multi(zz->Get_Num_Edges_Multi_Data,
                              num_gid_entries, num_lid_entries, num_obj,
                              global_ids, local_ids, *edges_per_obj, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Get_Num_Edges_Multi.");
        goto End;
      }

      for (i = 0; i < num_obj; i++) {
        nedges = (*edges_per_obj)[i];
        *num_edges += nedges;
        if (nedges > *max_edges) *max_edges = nedges;
      }
    }
    else {
      for (i=0; i< num_obj; i++) {
        lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
        nedges = zz->Get_Num_Edges(zz->Get_Num_Edges_Data,
                                   num_gid_entries, num_lid_entries,
                                   &(global_ids[i*num_gid_entries]),
                                   lid, &ierr);
        if (ierr) {
          /* Return error */
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Get_Num_Edges.");
          goto End;
        }
        *num_edges += nedges;
        if (nedges > *max_edges) *max_edges = nedges;
        (*edges_per_obj)[i] = nedges;
      }
    }
  }

End:
  return ierr;
}


/*******************************************************************
 * hash_lookup uses Zoltan_Hash to lookup a key 
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type ZOLTAN_ID_PTR (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   the global number of the element with key key,
 *   or -1 if the key is not in the hash table
 *
 *******************************************************************/

static int hash_lookup (ZZ *zz, struct Hash_Node **hashtab, ZOLTAN_ID_PTR key,
                        int n)
{
  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int)n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}


static int process_edge_list(
  ZZ *zz, 
  int i,                        /* Object being processed */
  ZOLTAN_ID_PTR gid,            /* GID of object being processed */
  int nedges,                   /* Number of edges this GID has */
  ZOLTAN_ID_PTR nbors_global,   /* GIDs of this object's neighboring vertices */
  int *nbors_proc,              /* Procs owning the GIDs in nbors_global */
  float *tmp_ewgts,             /* Edge weights for this GID's edges */
  struct Hash_Node **hashtab,   /* Hash table for global dense vertex #ing */
  int graph_type,               /* Global or local? */
  int num_gid_entries,          /* # of unsigned ints in a GID */
  int num_obj,                  /* Number of objects on this processor */
  int edge_wgt_dim,             /* Number of edge weights per edge */
  int *vtxdist,                 /* Distrib across procs of global dense #ing */
  int *xadj,                    /* Index of each obj's first edge in adjncy */
  int *adjncy,                  /* Edges for all objects. */
  float *ewgts,                 /* Edge weights for all edges */
  int *plist,                   /* Procs with nbors of this object */
  int *jj,                      /* Position to store non-self edge in adjncy. */
  int *nself,                   /* # of self-edges */
  int *cross_edges,             /* # of edges crossing processor boundaries */
  int *offset,      
  int *nsend,                   /* # of edges to be sent to other procs. */
  int *max_proc_list_len,       /* Max # of entries alloc'ed in proc_list */
  struct Edge_Info **proc_list, /* Edge info to be sent to other procs */
  ZOLTAN_ID_PTR *proc_list_nbor /* Global IDs of neighbors of proc_list
                                   entries.  This array is separate from
                                   proc_list to prevent individual mallocs
                                   for nbor global IDs.   */
)
{
/* Perform processing on list of edges for a single object. 
 * Add edges to graph data structure.
 * Add off-processor edges to communication lists.
 * This routine has lots of side effects; it is a subroutine only
 * to support both single-object and multi-object edge list callbacks.
 */
char *yo = "process_edge_list";
int j, k;
int tmp, flag;
struct Edge_Info *ptr;
int ierr = ZOLTAN_OK;

  /* Separate inter-processor edges from the local ones */
  /* Also remove any self-edges */
  for (j=0; j<nedges; j++){
    if (nbors_proc[j] == zz->Proc){
      /* Local edge */
      tmp = hash_lookup(zz, hashtab, 
                        &(nbors_global[j*num_gid_entries]), num_obj);
      if (tmp == i+(vtxdist)[zz->Proc]){
        /* Self-edge! */
        (*nself)++;
      }
      else{
        /* Copy over edge weights. */
        for (k=0; k<edge_wgt_dim; k++)
          ewgts[(*jj)*edge_wgt_dim+k] = tmp_ewgts[j*edge_wgt_dim+k];
        /* Put the global number into the adjacency array */
        adjncy[(*jj)++] = tmp;
      }
    } else {
      /* Inter-processor edge. */
      (*cross_edges)++;
      /* Skip this edge if local graphs have been requested. */
      if (graph_type == LOCAL_GRAPH)
        /* do nothing */ ;
      else {
        /* Check if we already have gid[i] in proc_list with */
        /* the same destination.                             */
        flag = 0;
        if (plist[nbors_proc[j]] < i){
          /* We need to send info about this edge */
          flag = 1;
          (*nsend)++; 
          plist[nbors_proc[j]] = i;
        }

        /* Check if we need to allocate more space for proc_list.*/
        if ((*offset) == (*max_proc_list_len)){
          if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
            printf("[%1d] Debug: Allocating more list space, "
                   "max_proc_list_len = %d, increasing by factor %f\n", 
                   zz->Proc, *max_proc_list_len, REALLOC_FACTOR);

          (*max_proc_list_len) *= REALLOC_FACTOR;
          (*proc_list) = (struct Edge_Info *) ZOLTAN_REALLOC(*proc_list,
                       *max_proc_list_len*sizeof(struct Edge_Info));
          (*proc_list_nbor) = ZOLTAN_REALLOC_GID_ARRAY(zz, *proc_list_nbor,
                            *max_proc_list_len);
          if (!(*proc_list)){
            /* Not enough memory */
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
            ierr = ZOLTAN_MEMERR;
            goto End;
          }
        }
        ptr = &(*proc_list)[(*offset)];
        ptr->my_gid = gid;
        ptr->my_gno = hash_lookup(zz, hashtab, gid, num_obj);
        ZOLTAN_SET_GID(zz, &((*proc_list_nbor)[(*offset)*num_gid_entries]),
                       &(nbors_global[j*num_gid_entries]));
        if (flag)
          ptr->nbor_proc = nbors_proc[j];
        else
          ptr->nbor_proc = -1;
        ptr->adj = &adjncy[(*jj)];

        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: proc_list[%1d] my_gid=", zz->Proc, (*offset));
          ZOLTAN_PRINT_GID(zz, ptr->my_gid);
          printf(", my_gno=%d, nbor_proc=%d\n", ptr->my_gno, ptr->nbor_proc);
        }

        /* Copy over edge weights. */
        for (k=0; k<edge_wgt_dim; k++)
          ewgts[(*jj)*edge_wgt_dim+k] = tmp_ewgts[j*edge_wgt_dim+k];

        /* Still don't know the global number, need to come back here later */
        adjncy[(*jj)++] = -1; 

        (*offset)++;
      }
    } /* end inter-proc edge */
  } /* end edge loop */

  xadj[i+1] = *jj; /* NB: We do not count self-edges. */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%3d] In %s: xadj[%3d] = %3d\n", zz->Proc, yo, i+1, *jj);

End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
