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
#include "lb_util_const.h"
#include "comm_const.h"
#include "parmetis_jostle_const.h"

/*********************************************************************/
/* Verify ParMetis graph structure.                                  */
/*                                                                   */
/* Input parameter check_graph specifies the level of verification:  */
/*   0 - perform no checks at all                                    */
/*   1 - verify on-processor part of graph                           */
/*   2 - verify complete graph (requires communication)              */
/*                                                                   */
/* Output: an error code (the same on all procs)                     */
/*                                                                   */
/* Fatal error :                                                     */
/*   - non-symmetric graph                                           */
/*   - incorrect vertex number                                       */
/*   - negative vertex or edge weight                                */
/*                                                                   */
/* Warning :                                                         */
/*   - zero sum of vertex or edge weights                            */
/*   - self-edge                                                     */
/*   - multiple edges between a pair of vertices                     */
/*                                                                   */
/*********************************************************************/

int LB_Verify_Graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj, 
       idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
       int vwgt_dim, int ewgt_dim, int check_graph, int debug_level)
{
  int i, j, ii, jj, k, kk, num_obj, nedges, ierr;
  int flag, cross_edges, mesg_size, sum;
  int nprocs, proc, *proclist, errors, global_errors;
  idxtype global_i, global_j;
  idxtype *ptr1, *ptr2;
  char *sendbuf, *recvbuf;
  struct Comm_Obj *comm_plan;
  static char *yo = "LB_Verify_Graph";
  char msg[256];

  ierr = LB_OK;
  if (check_graph == 0) /* perform no error checking at all */
     return ierr;

  /* Get number of procs and my rank */
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &proc);

  num_obj = vtxdist[proc+1] - vtxdist[proc];

  /* Verify that vertex weights are positive */
  if (vwgt_dim){
    for (i=0; i<num_obj; i++){
       sum = 0;
       for (k=0; k<vwgt_dim; k++){
         if (vwgt[i*vwgt_dim+k] < 0) {
            sprintf(msg, "Negative object weight of %d.", vwgt[i*vwgt_dim+k]);
            LB_PRINT_ERROR(proc, yo, msg);
            ierr = LB_FATAL;
            goto barrier1;
         }
         sum += vwgt[i*vwgt_dim+k];
       }
       if (sum == 0){
          sprintf(msg, "Zero weight sum for object %d (after scaling).", i);
          LB_PRINT_WARN(proc, yo, msg);
          ierr = LB_WARN;
       }
    }
  }

  /* Verify that edge weights are positive */
  nedges = xadj[num_obj];
  if (ewgt_dim){
    for (j=0; j<nedges; j++){
      sum = 0;
      for (k=0; k<ewgt_dim; k++){
        if (adjwgt[j*ewgt_dim+k] < 0) {
          sprintf(msg, "Negative communication weight of %d.", vwgt[j]);
          LB_PRINT_ERROR(proc, yo, msg);
          ierr = LB_FATAL;
          goto barrier1;
        }
        sum += adjwgt[j*ewgt_dim+k];
      }
      if (sum == 0){
         LB_PRINT_WARN(proc, yo, "Zero communication weight.");
         ierr = LB_WARN;
      }
    }
  }

  /* Verify that the graph is symmetric (edge weights, too) */
  /* Also check for self-edges and multi-edges */

  /* First pass: Check on-processor edges and count # off-proc edges */
  cross_edges = 0;
  for (i=0; i<num_obj; i++){
    global_i = vtxdist[proc]+i;
    for (ii=xadj[i]; ii<xadj[i+1]; ii++){
      global_j = adjncy[ii];
      /* Reasonable vertex value? */
      if ((global_j < vtxdist[0]) || (global_j >= vtxdist[nprocs])){
        sprintf(msg, "Edge to invalid vertex %d detected.", global_j);
        LB_PRINT_ERROR(proc, yo, msg);
        ierr = LB_FATAL;
        goto barrier1;
      }
      /* Self edge? */
      if (global_j == global_i){
        sprintf(msg, "Self edge for vertex %d detected.", global_i);
        LB_PRINT_WARN(proc, yo, msg);
        ierr = LB_WARN;
      }
      /* Duplicate edge? */
      for (kk=xadj[i]; kk<xadj[i+1]; kk++){
        if ((kk != ii) && (adjncy[kk] == adjncy[ii])){
          sprintf(msg, "Duplicate edge (%d,%d) detected.", global_i, global_j);
          LB_PRINT_WARN(proc, yo, msg);
          ierr = LB_WARN;
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
                   ierr = LB_FATAL;
              }
            }
            break;
          }
        }
        if (!flag) {
          sprintf(msg, "Graph is not symmetric. "
                  "Edge (%d,%d) exists, but no edge (%d,%d).", 
                  global_i, global_j, global_j, global_i);
          LB_PRINT_ERROR(proc, yo, msg);
          ierr = LB_FATAL;
          goto barrier1;
        }
      }
      else {
        cross_edges++;
      }
    }
  }

barrier1:
  /* Check if any processor has encountered a fatal error so far */
  errors = 0;
  if (ierr == LB_WARN)
    errors |= 1;
  if (ierr == LB_MEMERR)
    errors |= 2;
  if (ierr == LB_FATAL)
    errors |= 4;

  MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_BOR, comm);

  if (global_errors & 4){
    /* Fatal error: return now */
    return LB_FATAL;
  }

  if ((check_graph >= 2) && cross_edges) {
    /* Allocate space for off-proc data */
    mesg_size = (2+ewgt_dim)*sizeof(idxtype);
    sendbuf = (char *) LB_MALLOC(cross_edges*mesg_size);
    recvbuf = (char *) LB_MALLOC(cross_edges*mesg_size);
    proclist = (int *) LB_MALLOC(cross_edges*sizeof(int));

    if (!(sendbuf && recvbuf && proclist)){
       LB_PRINT_ERROR(proc, yo, "Out of memory.");
       ierr = LB_MEMERR;
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
    ierr = LB_Comm_Create(&comm_plan, cross_edges, proclist, comm, TAG1, 
                          TRUE, &k);
    if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
      sprintf(msg, "Error %s returned from LB_Comm_Create.", 
              (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
      LB_PRINT_ERROR(proc, yo, msg);
      ierr = (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
    }
    else {
      if (k != cross_edges){
        sprintf(msg, "Incorrect number of edges to/from proc %d.", proc);
        LB_PRINT_ERROR(proc, yo, msg);
        ierr = LB_FATAL;
      }

      ierr = LB_Comm_Do(comm_plan, TAG2, sendbuf, mesg_size, recvbuf);
      LB_Comm_Destroy(&comm_plan);
      if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
        sprintf(msg, "Error %s returned from LB_Comm_Do.",
                (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
        LB_PRINT_ERROR(proc, yo, msg);
        ierr = (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
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
                  LB_PRINT_ERROR(proc, yo, msg);
                  ierr = LB_FATAL;
                }
              }
            }
          }
          if (!flag){
            sprintf(msg, "Graph is not symmetric.  "
                    "Edge (%d,%d) exists, but not (%d,%d).", 
                    ptr1[0], ptr1[1], ptr1[1], ptr1[0]);
            LB_PRINT_ERROR(proc, yo, msg);
            ierr = LB_FATAL;
          }
        }
      }
    }

    /* Free memory */
    LB_FREE(&sendbuf);
    LB_FREE(&recvbuf);
    LB_FREE(&proclist);
  }

  /* Compute global error code */
  errors = 0;
  if (ierr == LB_WARN)
    errors |= 1;
  if (ierr == LB_MEMERR)
    errors |= 2;
  if (ierr == LB_FATAL)
    errors |= 4;

  MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_BOR, comm);

  if (global_errors & 4){
    return LB_FATAL;
  }
  else if (global_errors & 2){
    return LB_MEMERR;
  }
  else if (global_errors & 1){
    return LB_WARN;
  }
  else {
    if (debug_level && (proc==0)){
      if (check_graph == 1)
        printf("ZOLTAN %s: graph is OK locally on each processor\n", yo);
      else /* check_graph >= 2 */
        printf("ZOLTAN %s: graph is OK\n", yo);
    }
    return LB_OK;
  }

}

