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
static char *cvs_dr_evalc_id = "$Id$";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "dr_const.h"
#include "dr_eval_const.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Utility functions for evaluating a partition.
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void driver_eval()
{
/*
 * Function to evaluate a partition.  Largely duplicates functionality
 * of LB_Eval, but provides sanity checking.
 */
int i;
int proc;
int cuts;
int gsumcuts, gmaxcuts, gmincuts;
int gsumelems, gmaxelems, gminelems;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  cuts = 0;
  for (i = 0; i < Mesh.necmap; i++) {
    cuts += Mesh.ecmap_cnt[i];
  }

  MPI_Allreduce(&cuts, &gsumcuts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmaxcuts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmincuts, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  MPI_Allreduce(&(Mesh.num_elems), &gsumelems, 1, MPI_INT, MPI_SUM, 
                MPI_COMM_WORLD);
  MPI_Allreduce(&(Mesh.num_elems), &gmaxelems, 1, MPI_INT, MPI_MAX, 
                MPI_COMM_WORLD);
  MPI_Allreduce(&(Mesh.num_elems), &gminelems, 1, MPI_INT, MPI_MIN, 
                MPI_COMM_WORLD);

  if (proc == 0) {
    printf("DRIVER EVAL:  load:  max %d  min %d  sum %d\n", 
           gmaxelems, gminelems, gsumelems);
    printf("DRIVER EVAL:  cuts:  max %d  min %d  sum %d\n",
           gmaxcuts, gmincuts, gsumcuts);
  }
}
