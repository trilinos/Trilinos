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

#include "zz_const.h"
#include "parmetis_jostle.h"

/*********************************************************************/
/* Verify ParMetis graph structure.                                  */
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
/*   - non-symmetric graph (or edge weights)                         */
/*   - incorrect vertex number                                       */
/*   - negative vertex or edge weight                                */
/*                                                                   */
/* Warning :                                                         */
/*   - no vertices or no edges in graph                              */
/*   - zero sum of vertex or edge weights                            */
/*   - self-edge                                                     */
/*   - multiple edges between a pair of vertices                     */
/*   - singletons (vertices with no edges)                           */
/*                                                                   */
/*********************************************************************/

int Zoltan_Verify_Graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj, 
       idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
       int vwgt_dim, int ewgt_dim, int check_graph, int output_level)
{
  int i, j, ii, jj, k, kk, num_obj, nedges, ierr;
  int flag, cross_edges, mesg_size, sum, global_sum;
  int nprocs, proc, *proclist, errors, global_errors;
  int num_zeros, num_selfs, num_duplicates, num_singletons;
  idxtype global_i, global_j;
  idxtype *ptr1, *ptr2;
  char *sendbuf, *recvbuf;
  ZOLTAN_COMM_OBJ *comm_plan;
  static char *yo = "Zoltan_Verify_Graph";
  char msg[256];

  ierr = ZOLTAN_OK;
  if (check_graph == 0) /* perform no error checking at all */
     return ierr;

  /* Get number of procs and my rank */
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &proc);

  /* Check number of vertices (objects) */
  num_obj = vtxdist[proc+1] - vtxdist[proc];
  MPI_Reduce(&num_obj, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum==0)){
    ierr = ZOLTAN_WARN;
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
            sprintf(msg, "Negative object weight of %d for object %d.", 
                    vwgt[i*vwgt_dim+k], i);
            ZOLTAN_PRINT_ERROR(proc, yo, msg);
            ierr = ZOLTAN_FATAL;
            goto barrier1;
         }
         sum += vwgt[i*vwgt_dim+k];
       }
       if (sum == 0){
          num_zeros++;
          if (output_level>1) {
            sprintf(msg, "Zero vertex (object) weights for object %d.", i);
            ZOLTAN_PRINT_WARN(proc, yo, msg);
          }
          ierr = ZOLTAN_WARN;
       }
    }
    MPI_Reduce(&num_zeros, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
    if ((proc==0) && (global_sum>0)){
      ierr = ZOLTAN_WARN;
      if (output_level>0){
        sprintf(msg, "%d objects have zero weights.", global_sum);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
  }

  /* Check number of edges */
  nedges = xadj[num_obj];
  MPI_Reduce(&nedges, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum==0)){
    ierr = ZOLTAN_WARN;
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
          sprintf(msg, "Negative communication weight of %d in edge %d.", 
                  adjwgt[j*ewgt_dim+k], j);
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
          ierr = ZOLTAN_FATAL;
          goto barrier1;
        }
        sum += adjwgt[j*ewgt_dim+k];
      }
      if (sum == 0){
        num_zeros++;
        if (output_level>1) {
          sprintf(msg, "Zero edge (communication) weights for edge %d.", j);
          ZOLTAN_PRINT_WARN(proc, yo, msg);
        }
        ierr = ZOLTAN_WARN;
      }
    }

    MPI_Reduce(&num_zeros, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
    if ((proc==0) && (global_sum>0)){
      ierr = ZOLTAN_WARN;
      if (output_level>0){
        sprintf(msg, "%d edges have zero weights.", global_sum);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
  }

  /* Verify that the graph is symmetric (edge weights, too) */
  /* Also check for self-edges and multi-edges */

  /* First pass: Check on-processor edges and count # off-proc edges */
  cross_edges = 0;
  num_selfs = 0;
  num_duplicates = 0;
  num_singletons = 0;
  for (i=0; i<num_obj; i++){
    global_i = vtxdist[proc]+i;
    /* Singleton? */
    if (xadj[i] == xadj[i+1]){
      num_singletons++;
      if (output_level>1){
        sprintf(msg, "Vertex %d has no edges.", global_i);
        ZOLTAN_PRINT_WARN(proc, yo, msg);
      }
    }
    for (ii=xadj[i]; ii<xadj[i+1]; ii++){
      global_j = adjncy[ii];
      /* Valid vertex number? */
      if ((global_j < vtxdist[0]) || (global_j >= vtxdist[nprocs])){
        sprintf(msg, "Edge to invalid vertex %d detected.", global_j);
        ZOLTAN_PRINT_ERROR(proc, yo, msg);
        ierr = ZOLTAN_FATAL;
        goto barrier1;
      }
      /* Self edge? */
      if (global_j == global_i){
        num_selfs++;
        if (output_level>1){
          sprintf(msg, "Self edge for vertex %d detected.", global_i);
          ZOLTAN_PRINT_WARN(proc, yo, msg);
        }
      }
      /* Duplicate edge? */
      for (kk=xadj[i]; kk<xadj[i+1]; kk++){
        if ((kk != ii) && (adjncy[kk] == adjncy[ii])){
          num_duplicates++;
          if (output_level>1){
            sprintf(msg, "Duplicate edge (%d,%d) detected.", global_i, global_j);
            ZOLTAN_PRINT_WARN(proc, yo, msg);
          }
        }
      }
      /* Is global_j a local vertex? */
      if ((global_j >= vtxdist[proc]) && (global_j < vtxdist[proc+1])){
        /* Check if (global_j, global_i) is an edge */
        j = global_j - vtxdist[proc];
        flag = 0;
        for (jj=xadj[j]; jj<xadj[j+1]; jj++){
          if (adjncy[jj] == global_i) {
            flag = 1;
            if (ewgt_dim){
              /* Compare weights */
              for (k=0; k<ewgt_dim; k++){
                if (adjwgt[jj*ewgt_dim+k] != adjwgt[ii*ewgt_dim+k])
                   ierr = ZOLTAN_FATAL;
              }
            }
            break;
          }
        }
        if (!flag) {
          sprintf(msg, "Graph is not symmetric. "
                  "Edge (%d,%d) exists, but no edge (%d,%d).", 
                  global_i, global_j, global_j, global_i);
          ZOLTAN_PRINT_ERROR(proc, yo, msg);
          ierr = ZOLTAN_FATAL;
          goto barrier1;
        }
      }
      else {
        cross_edges++;
      }
    }
  }

barrier1:
  /* Sum up warnings so far. */
  MPI_Reduce(&num_selfs, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg, "%d self-edges in graph.", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  MPI_Reduce(&num_duplicates, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg, "%d duplicate edges in graph.", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  MPI_Reduce(&num_singletons, &global_sum, 1, MPI_INT, MPI_SUM, 0, comm);
  if ((proc==0) && (global_sum>0)){
    ierr = ZOLTAN_WARN;
    if (output_level>0){
      sprintf(msg, "%d vertices in the graph are singletons (have no edges).", global_sum);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
    }
  }
  
  
  /* Check if any processor has encountered a fatal error so far */
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
    return ZOLTAN_FATAL;
  }

  if (check_graph >= 2) {
    /* Allocate space for off-proc data */
    mesg_size = (2+ewgt_dim)*sizeof(idxtype);
    sendbuf = (char *) ZOLTAN_MALLOC(cross_edges*mesg_size);
    recvbuf = (char *) ZOLTAN_MALLOC(cross_edges*mesg_size);
    proclist = (int *) ZOLTAN_MALLOC(cross_edges*sizeof(int));

    if (cross_edges && !(sendbuf && recvbuf && proclist)){
       ZOLTAN_PRINT_ERROR(proc, yo, "Out of memory.");
       ierr = ZOLTAN_MEMERR;
    }

    /* Second pass: Copy data to send buffer */
    nedges = 0;
    ptr1 = (idxtype *) sendbuf;
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
           /* Copy (global_i, global_j) and corresponding weights to sendbuf */
           *ptr1++ = global_i;
           *ptr1++ = global_j;
           for (k=0; k<ewgt_dim; k++){
             *ptr1++ = adjwgt[ii*ewgt_dim+k];
           }
        }
      }
    }

    /* Do the irregular communication */
    ierr = Zoltan_Comm_Create(&comm_plan, cross_edges, proclist, comm, 
                          TAG1, &k);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
              (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
      ZOLTAN_PRINT_ERROR(proc, yo, msg);
    }
    else {
      if (k != cross_edges){
        sprintf(msg, "Incorrect number of edges to/from proc %d.", proc);
        ZOLTAN_PRINT_ERROR(proc, yo, msg);
        ierr = ZOLTAN_FATAL;
      }

      ierr = Zoltan_Comm_Do(comm_plan, TAG2, sendbuf, mesg_size, recvbuf);
      Zoltan_Comm_Destroy(&comm_plan);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        sprintf(msg, "Error %s returned from Zoltan_Comm_Do.",
                (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
        ZOLTAN_PRINT_ERROR(proc, yo, msg);
      }
      else {



        /* Third pass: Compare on-proc data to the off-proc data we received */
        /* sendbuf and recvbuf should contain the same data except (i,j) is  */
        /* (j,i)                                                             */
        for (i=0, ptr1=(idxtype *)sendbuf; i<cross_edges; 
             i++, ptr1 += (2+ewgt_dim)){
          flag = 0;
          for (j=0, ptr2=(idxtype *)recvbuf; j<cross_edges; 
               j++, ptr2 += (2+ewgt_dim)){
            if ((ptr2[0] == ptr1[1]) && (ptr2[1] == ptr1[0])){
              /* Found matching edge */
              flag = 1;
              /* Check weights */
              for (k=0; k<ewgt_dim; k++){
                if (ptr1[2+k] != ptr2[2+k]){
                  sprintf(msg, "Edge weight (%d,%d) is not symmetric: %d != %d",
                           ptr1[0], ptr1[1], ptr1[2+k], ptr2[2+k]);
                  ZOLTAN_PRINT_ERROR(proc, yo, msg);
                  ierr = ZOLTAN_FATAL;
                }
              }
            }
          }
          if (!flag){
            sprintf(msg, "Graph is not symmetric.  "
                    "Edge (%d,%d) exists, but not (%d,%d).", 
                    ptr1[0], ptr1[1], ptr1[1], ptr1[0]);
            ZOLTAN_PRINT_ERROR(proc, yo, msg);
            ierr = ZOLTAN_FATAL;
          }
        }
      }
    }

    /* Free memory */
    ZOLTAN_FREE(&sendbuf);
    ZOLTAN_FREE(&recvbuf);
    ZOLTAN_FREE(&proclist);
  }

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

