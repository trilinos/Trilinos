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

#include "lb_const.h"
#include "lb_util_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_jostle_const.h"
#include "params_const.h"

/*********************************************************************/
/* Verify ParMetis graph structure                                   */
/* Only 1-dimensional weights are allowed.                           */
/*********************************************************************/
int LB_verify_graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
              idxtype *vwgt, idxtype *adjwgt, int wgtflag, int check_graph, char *yo)
{
  int i, j, num_obj, nedges, ierr;
  int nprocs, proc;

  if (check_graph == 0) /* perform no error checking at all */
     return LB_OK;

  /* Get number of procs and my rank */
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &proc);

  num_obj = vtxdist[proc+1] - vtxdist[proc];

  /* Verify that vertex weights are positive */
  if (wgtflag & 2){
    for (i=0; i<num_obj; i++){
       if (vwgt[i] < 0) {
          fprintf(stderr, "Zoltan error: Negative object weight of %g in %s\n", 
                  vwgt[i], yo);
          ierr = LB_FATAL;
       }
       else if (vwgt[i] == 0){
          fprintf(stderr, "Zoltan warning: Zero object weight in %s (after scaling)\n", 
                  yo);
          ierr = LB_WARN;
       }
    }
  }

  /* Verify that edge weights are positive */
  nedges = xadj[num_obj];
  if (wgtflag & 1){
    for (j=0; j<nedges; j++){
       if (adjwgt[j] < 0) {
          fprintf(stderr, "Zoltan error: Negative communication weight of %g in %s\n", 
          vwgt[j], yo);
          ierr = LB_FATAL;
       }
       else if (adjwgt[j] == 0){
          fprintf(stderr, "Zoltan warning: Zero communication weight in %s\n", yo);
          ierr = LB_WARN;
       }
    }
  }

  /* Verify that the graph is symmetric (edge weights, too) */
  /* Also check for self-edges and multi-edges */

  return ierr;
}

