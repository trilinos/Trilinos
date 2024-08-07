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


#include "zz_util_const.h"
#include "third_library_const.h"
#include "third_library_tools.h"
#include "graph_util.h"

/* comparison routine for bsearch */
static int Zoltan_Compare_Indextypes(const void *key, const void *arg);

#include "zz_sort.h"

/*********************************************************************/
/* Verify ParMetis graph structure.                                  */
/*                                                                   */
/* Input parameter graph_type specifies the type of graph:           */
/*   LOCAL_GRAPH  - each proc owns a local METIS style graph         */
/*   GLOBAL_GRAPH - the procs share a global ParMETIS graph          */
/*                                                                   */
/* Input parameter check_graph specifies the level of verification:  */
/*   0 - perform no checks at all                                    */
/*   1 - verify on-processor part of graph                           */
/*   2 - verify complete graph (requires communication)              */
/*                                                                   */
/* Input parameter output_level specifies the level of verbosity:    */
/*   0 - suppress warnings                                           */
/*   1 - print summary of warnings                                   */
/*   2 - detailed output for each warning                            */
/*                                                                   */
/* Output: an error code (the same on all procs)                     */
/*                                                                   */
/* Fatal error :                                                     */
/*   - non-symmetric graph                                           */
/*   - incorrect vertex number                                       */
/*   - negative vertex or edge weight                                */
/*                                                                   */
/* Warning :                                                         */
/*   - no vertices or no edges in graph                              */
/*   - zero sum of vertex or edge weights                            */
/*   - self-edge                                                     */
/*   - multiple edges between a pair of vertices                     */
/*   - singletons (vertices with no edges)                           */
/*   - non-symmetric edge weights                                    */
/*                                                                   */
/*********************************************************************/

int Zoltan_Verify_Graph(MPI_Comm comm, indextype *vtxdist, indextype *xadj, 
       indextype *adjncy, weighttype *vwgt, weighttype *adjwgt, 
       int vwgt_dim, int ewgt_dim, 
       int graph_type, int check_graph, int output_level)
{
  int flag, cross_edges = 0, mesg_size, sum;
  int ierr;  /* Error in the graph; always return the worst ierr found */
  int runerr;/* Error running this routine */
  int nprocs, proc, *proclist, errors, global_errors;
  int *perm=NULL; 
  int free_adjncy_sort=0;
  ZOLTAN_COMM_OBJ *comm_plan;
  static char *yo = "Zoltan_Verify_Graph";
  char msg[256];
  ZOLTAN_GNO_TYPE num_obj;
  int nrecv;
  indextype *ptr, *ptr1, *ptr2;
  indextype global_i, global_j;
  indextype *sendgno=NULL, *recvgno=NULL, *adjncy_sort=NULL;
  weighttype *sendwgt, *recvwgt;
  ZOLTAN_GNO_TYPE num_duplicates, num_singletons;
  ZOLTAN_GNO_TYPE num_selfs, nedges, global_sum, num_zeros;
  ZOLTAN_GNO_TYPE i, j, ii, k;
  MPI_Datatype zoltan_gno_mpi_type;

  ierr = ZOLTAN_OK;
  runerr = ZOLTAN_OK;
  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  /* Make sure all procs have same value of check_graph. */
  MPI_Allreduce(&check_graph, &i, 1, MPI_INT, MPI_MAX, comm);
  check_graph = i;

  if (check_graph == 0) /* perform no error checking at all */
     return ierr;

  /* Get number of procs and my rank */
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &proc);

  /* Check number of vertices (objects) */
  num_obj = (ZOLTAN_GNO_TYPE)(vtxdist[proc+1] - vtxdist[proc]);
  MPI_Reduce(&num_obj, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum==0)){
    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
    if (output_level>0)
      ZOLTAN_PRINT_WARN(proc, yo, "No vertices in graph.");
  }

  /* Verify that vertex weights are non-negative */
  num_zeros = 0;
  if (vwgt_dim){
    for (i=0; i<num_obj; i++){
       sum = 0;
       for (k=0; k<vwgt_dim; k++){
         if (vwgt[i*vwgt_dim+k] < 0) {
            sprintf(msg, "Negative object weight of " TPL_WGT_SPEC " for object " ZOLTAN_GNO_SPEC ".", 
                    vwgt[i*vwgt_dim+k], i);
            ZOLTAN_PRINT_ERROR(proc, yo, msg);
            ierr = ZOLTAN_FATAL;
         }
         sum += vwgt[i*vwgt_dim+k];
       }
       if (sum == 0){
          num_zeros++;
          if (output_level>1) {
            sprintf(msg, "Zero vertex (object) weights for object " ZOLTAN_GNO_SPEC ".", i);
            ZOLTAN_PRINT_WARN(proc, yo, msg);
          }
          if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
       }
    }
    MPI_Reduce(&num_zeros, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
    if ((proc==0) && (global_sum>0)){
      if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
      if (output_level>0){
        sprintf(msg,  ZOLTAN_GNO_SPEC " objects have zero weights.", global_sum);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
  }

  /* Check number of edges */
  nedges = (ZOLTAN_GNO_TYPE)xadj[num_obj];
  MPI_Reduce(&nedges, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum==0)){
    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
    if (output_level>0)
      ZOLTAN_PRINT_WARN(proc, yo, "No edges in graph.");
  }

  /* Verify that edge weights are non-negative */
  num_zeros = 0;
  if (ewgt_dim){
    for (j=0; j<nedges; j++){
      sum = 0;
      for (k=0; k<ewgt_dim; k++){
        if (adjwgt[j*ewgt_dim+k] < 0) {
          sprintf(msg, "Negative edge weight of " TPL_WGT_SPEC " in edge " ZOLTAN_GNO_SPEC ".", 
                  adjwgt[j*ewgt_dim+k], j);
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
          ierr = ZOLTAN_FATAL;
        }
        sum += adjwgt[j*ewgt_dim+k];
      }
      if (sum == 0){
        num_zeros++;
        if (output_level>1) {
          sprintf(msg, "Zero edge (communication) weights for edge " ZOLTAN_GNO_SPEC ".", j);
          ZOLTAN_PRINT_WARN(proc, yo, msg);
        }
        if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
      }
    }

    MPI_Reduce(&num_zeros, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
    if ((proc==0) && (global_sum>0)){
      if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
      if (output_level>0){
        sprintf(msg,  ZOLTAN_GNO_SPEC " edges have zero weights.", global_sum);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
  }

  /* Verify that the graph is symmetric (edge weights, too) */
  /* Also check for self-edges and multi-edges */

  /* Pre-processing: Check if edge lists are sorted. If not, */
  /* make a copy and sort so we can save time in the lookups. */
  flag = 0; /* Assume sorted. */
  for (i=0; (i<num_obj) && (flag==0); i++){
    for (ii=xadj[i]; ii<xadj[i+1]-1; ii++){
      if (adjncy[ii] > adjncy[ii+1]){
        flag = 1; /* Not sorted. */
        break; 
      }
    }
  }
  if (flag){ /* Need to sort. */
    adjncy_sort = (indextype *) ZOLTAN_MALLOC(nedges*sizeof(indextype));
    perm = (int *) ZOLTAN_MALLOC(nedges*sizeof(int));
    free_adjncy_sort = 1;
    if (nedges && (!adjncy_sort || !perm)){
      /* Out of memory. */
      ZOLTAN_PRINT_ERROR(proc, yo, "Out of memory.");
      runerr = ZOLTAN_MEMERR;
    }
    else {
      for (k=0; k<nedges; k++){
        adjncy_sort[k] = adjncy[k];
        perm[k] = k;
      }
      if (sizeof(indextype) == sizeof(short)){
        for (i=0; i<num_obj; i++) 
          Zoltan_quicksort_list_inc_short((short *)adjncy_sort, perm, 
                                          (int)xadj[i], (int)xadj[i+1]-1);
      }
      else if (sizeof(indextype) == sizeof(int)){
        for (i=0; i<num_obj; i++) 
          Zoltan_quicksort_list_inc_int((int *)adjncy_sort, perm, 
                                        (int)xadj[i], (int)xadj[i+1]-1);
      }
      else if (sizeof(indextype) == sizeof(long)){
        for (i=0; i<num_obj; i++) 
          Zoltan_quicksort_list_inc_long((long *)adjncy_sort, perm, 
                                         (int)xadj[i], (int)xadj[i+1]-1);
      }
      else if (sizeof(indextype) == sizeof(int64_t)){
        for (i=0; i<num_obj; i++) 
          Zoltan_quicksort_list_inc_long_long((int64_t*)adjncy_sort, perm, 
                                              (int)xadj[i], (int)xadj[i+1]-1);
      }
      else{
        ZOLTAN_PRINT_ERROR(proc, yo, 
                           "Error in third party library data type support.");
        ierr = ZOLTAN_FATAL;
      }
    }
  }
  else { /* Already sorted. */
    adjncy_sort = adjncy;
  }

  /* First pass: Check on-processor edges and count # off-proc edges */
  cross_edges = 0;
  num_selfs = 0;
  num_duplicates = 0;
  num_singletons = 0;
  for (i=0; i<num_obj; i++){
    if (IS_GLOBAL_GRAPH(graph_type)){
      global_i = vtxdist[proc]+i;
    }
    else{ /* graph_type == LOCAL_GRAPH */
      global_i = i; /* A bit confusingly, global_i = i for local graphs */
    }
    /* Singleton? */
    if (xadj[i] == xadj[i+1]){
      num_singletons++;
      if (output_level>1){
        sprintf(msg, "Vertex " TPL_IDX_SPEC " has no edges.", global_i);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
    for (ii=xadj[i]; ii<xadj[i+1]; ii++){
      global_j = adjncy_sort[ii];
      /* Valid vertex number? */
      if ((IS_GLOBAL_GRAPH(graph_type) &&
           ((global_j < vtxdist[0]) || (global_j >= vtxdist[nprocs])))
          || (IS_LOCAL_GRAPH(graph_type) && 
           ((global_j < 0) || (global_j >= num_obj)))){
        sprintf(msg, "Edge to invalid vertex " TPL_IDX_SPEC " detected.", global_j);
        ZOLTAN_PRINT_ERROR(proc, yo, msg);
        ierr = ZOLTAN_FATAL;
      }
      /* Self edge? */
      if (global_j == global_i){
        num_selfs++;
        if (output_level>1){
          sprintf(msg, "Self edge for vertex " TPL_IDX_SPEC " detected.", global_i);
          ZOLTAN_PRINT_WARN(proc, yo, msg);
        }
      }
      /* Duplicate edge? */
      if ((ii+1<xadj[i+1]) && (adjncy_sort[ii]==adjncy_sort[ii+1])){
        num_duplicates++;
        if (output_level>1){
          sprintf(msg, "Duplicate edge (" TPL_IDX_SPEC "," TPL_IDX_SPEC ") detected.", global_i, global_j);
          ZOLTAN_PRINT_WARN(proc, yo, msg);
        }
      }
      /* Is global_j a local vertex? */
      if (IS_LOCAL_GRAPH(graph_type) || (IS_GLOBAL_GRAPH(graph_type) &&
          (global_j >= vtxdist[proc]) && (global_j < vtxdist[proc+1]))){
        /* Check if (global_j, global_i) is an edge */
        if (IS_GLOBAL_GRAPH(graph_type))
          j = global_j - vtxdist[proc];
        else /* graph_type == LOCAL_GRAPH */
          j = global_j;
        /* Binary search for edge (global_j, global_i) */
        ptr = (indextype *)bsearch(&global_i, &adjncy_sort[xadj[j]], (int)(xadj[j+1]-xadj[j]),
              sizeof(indextype), Zoltan_Compare_Indextypes);
        if (ptr){
          /* OK, found edge (global_j, global_i) */
          if ((adjncy_sort==adjncy) && ewgt_dim){
            /* Compare weights */
            /* EBEB For now, don't compare weights if we sorted edge lists. */
            flag = 0;
            for (k=0; k<ewgt_dim; k++){
              if (adjwgt[(ptr-adjncy_sort)*ewgt_dim+k] != 
                                           adjwgt[ii*ewgt_dim+k]){
                 /* Numerically nonsymmetric */
                 flag = -1;
                 if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
              }
            }
            if (flag<0 && output_level>0){
              sprintf(msg, "Graph is numerically nonsymmetric "
                "in edge (" TPL_IDX_SPEC "," TPL_IDX_SPEC ")", global_i, global_j);
              ZOLTAN_PRINT_WARN(proc, yo, msg);
            }
          }
        }
        else { /* bsearch failed */
          sprintf(msg, "Graph is not symmetric. "
                  "Edge (" TPL_IDX_SPEC "," TPL_IDX_SPEC ") exists, but no edge (" TPL_IDX_SPEC "," TPL_IDX_SPEC ").", 
                  global_i, global_j, global_j, global_i);
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
          ierr = ZOLTAN_FATAL;
        }
      }
      else {
        cross_edges++;
      }
    }
  }

  /* Sum up warnings so far. */
  MPI_Reduce(&num_selfs, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg,  ZOLTAN_GNO_SPEC " self-edges in graph.", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  MPI_Reduce(&num_duplicates, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg,  ZOLTAN_GNO_SPEC " duplicate edges in graph.", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  MPI_Reduce(&num_singletons, &global_sum, 1, zoltan_gno_mpi_type, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg,  ZOLTAN_GNO_SPEC " vertices in the graph are singletons (have no edges).", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  
  
  /* Check if any processor has encountered an error so far */
  errors = 0;
  if (ierr == ZOLTAN_WARN)
    errors |= 1;
  if (ierr == ZOLTAN_MEMERR)
    errors |= 2;
  if (ierr == ZOLTAN_FATAL)
    errors |= 4;

  MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_BOR, comm);

  if (global_errors & 4){
    /* Fatal error: return now */
    if (free_adjncy_sort) ZOLTAN_FREE(&adjncy_sort);
    if (free_adjncy_sort) ZOLTAN_FREE(&perm);
    return ZOLTAN_FATAL;
  }

  if ((IS_GLOBAL_GRAPH(graph_type)) && (check_graph >= 2)) {
    /* Test for consistency across processors. */

    /* Allocate space for off-proc data */
    mesg_size = (2*sizeof(indextype)) + (ewgt_dim * sizeof(weighttype));
    sendgno = (indextype *) ZOLTAN_MALLOC(cross_edges*mesg_size);
    recvgno = (indextype  *) ZOLTAN_MALLOC(cross_edges*mesg_size);
    proclist = (int *) ZOLTAN_MALLOC(cross_edges*sizeof(int));

    if (cross_edges && !(sendgno && recvgno && proclist)){
       ZOLTAN_PRINT_ERROR(proc, yo, "Out of memory.");
       runerr = ZOLTAN_MEMERR;
    }
    else {

      /* Second pass: Copy data to send buffer */
      nedges = 0;
      ptr1 = (indextype *) sendgno;
      for (i=0; i<num_obj; i++){
        global_i = vtxdist[proc]+i;
        for (ii=xadj[i]; ii<xadj[i+1]; ii++){
          global_j = adjncy[ii];
          /* Is global_j off-proc? */
          if ((global_j < vtxdist[proc]) || (global_j >= vtxdist[proc+1])){
             /* Add to list */
             k=0; 
             while (global_j >= vtxdist[k+1]) k++;
             proclist[nedges++] = k;
             /* Copy (global_i,global_j) and corresponding weights to sendgno */
             *ptr1++ = global_i;
             *ptr1++ = global_j;
             for (k=0; k<ewgt_dim; k++){
               *ptr1++ = adjwgt[ii*ewgt_dim+k];
             }
          }
        }
      }

      /* Do the irregular communication */
      runerr = Zoltan_Comm_Create(&comm_plan, cross_edges, proclist, comm,
                                  TAG1, &nrecv);
      if (runerr != ZOLTAN_OK && runerr != ZOLTAN_WARN) {
        sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
                (runerr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
        Zoltan_Comm_Destroy(&comm_plan);
        ZOLTAN_PRINT_ERROR(proc, yo, msg);
      }
      else {
        if (nrecv != cross_edges){
          sprintf(msg, "Incorrect number of edges to/from proc %d.", proc);
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
          ierr = ZOLTAN_FATAL;
        }
  
        runerr = Zoltan_Comm_Do(comm_plan, TAG2, (char *)sendgno, mesg_size,
                               (char *)recvgno);
        Zoltan_Comm_Destroy(&comm_plan);
        if (runerr != ZOLTAN_OK && runerr != ZOLTAN_WARN) {
          sprintf(msg, "Error %s returned from Zoltan_Comm_Do.",
                  (runerr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
        }
        else {
  
          /* Third pass: Compare on-proc data to off-proc data we received */
          /* sendgno and recvgno should contain the same data except (i,j) is */
          /* (j,i)                                                            */
          for (i=0, ptr1=sendgno; i<cross_edges; i++){
            flag = 0;
            sendwgt = (weighttype *)(ptr1 + 2);
            for (j=0, ptr2=recvgno; j<cross_edges; j++){
              recvwgt = (weighttype *)(ptr2 + 2);
              if ((ptr2[0] == ptr1[1]) && (ptr2[1] == ptr1[0])){
                /* Found matching edge */
                flag = 1;
                /* Check weights */
                for (k=0; k<ewgt_dim; k++){
                  if (sendwgt[k] != recvwgt[k]){
                    flag = -1;
                    if (ierr == ZOLTAN_OK) ierr = ZOLTAN_WARN;
                  }
                }
                if (flag<0 && output_level>0){
                    sprintf(msg, "Edge weight (" TPL_IDX_SPEC "," TPL_IDX_SPEC ") is not symmetric",
                            ptr1[0], ptr1[1]);
                    ZOLTAN_PRINT_WARN(proc, yo, msg);
                }
              }
              ptr2 = (indextype *)(recvwgt + ewgt_dim);
            }
            if (!flag){
              sprintf(msg, "Graph is not symmetric.  "
                      "Edge (" TPL_IDX_SPEC "," TPL_IDX_SPEC ") exists, but not (" TPL_IDX_SPEC "," TPL_IDX_SPEC ").", 
                      ptr1[0], ptr1[1], ptr1[1], ptr1[0]);
              ZOLTAN_PRINT_ERROR(proc, yo, msg);
              ierr = ZOLTAN_FATAL;
            }
            ptr1 = (indextype *)(sendwgt + ewgt_dim);
          }
        }
      }
    }

    /* Free memory */
    ZOLTAN_FREE(&sendgno);
    ZOLTAN_FREE(&recvgno);
    ZOLTAN_FREE(&proclist);
  }

  if (free_adjncy_sort) ZOLTAN_FREE(&adjncy_sort);
  if (free_adjncy_sort) ZOLTAN_FREE(&perm);

  /* Compute global error code */
  errors = 0;
  if (ierr == ZOLTAN_WARN)
    errors |= 1;
  if (ierr == ZOLTAN_MEMERR)
    errors |= 2;
  if (ierr == ZOLTAN_FATAL)
    errors |= 4;

  MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_BOR, comm);

  if (global_errors & 4){
    return ZOLTAN_FATAL;
  }
  else if (global_errors & 2){
    return ZOLTAN_MEMERR;
  }
  else if (global_errors & 1){
    return ZOLTAN_WARN;
  }
  else {
    if (proc==0 && output_level>0){
      printf("ZOLTAN %s: The graph is valid with check_graph = %1d\n", 
             yo, check_graph);
    }
    return ZOLTAN_OK;
  }

}

/* comparison routine for bsearch */
static int Zoltan_Compare_Indextypes(const void *key, const void *arg)
{
   if ( *(indextype*) key > (*(indextype*) arg))  return  1;
   if ( *(indextype*) key < (*(indextype*) arg))  return -1;

   return 0;  /* equal */
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
