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

void driver_eval(MESH_INFO_PTR mesh)
{
/*
 * Function to evaluate a partition.  Largely duplicates functionality
 * of Zoltan_LB_Eval, but provides sanity checking.
 *
 * Currently uses only the first cpu weight.
 */
int i;
int proc;
int cuts = 0;
float load = 0.;
int gsumcuts, gmaxcuts, gmincuts;
int gsumelems, gmaxelems, gminelems;
float gsumload, gmaxload, gminload;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  for (i = 0; i < mesh->necmap; i++) {
    cuts += mesh->ecmap_cnt[i];
  }
  
  for (i = 0; i < mesh->num_elems; i++) {
    load += mesh->elements[i].cpu_wgt[0];
  }

  MPI_Allreduce(&cuts, &gsumcuts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmaxcuts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmincuts, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  MPI_Allreduce(&(mesh->num_elems), &gsumelems, 1, MPI_INT, MPI_SUM, 
                MPI_COMM_WORLD);
  MPI_Allreduce(&(mesh->num_elems), &gmaxelems, 1, MPI_INT, MPI_MAX, 
                MPI_COMM_WORLD);
  MPI_Allreduce(&(mesh->num_elems), &gminelems, 1, MPI_INT, MPI_MIN, 
                MPI_COMM_WORLD);

  MPI_Allreduce(&load, &gsumload, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&load, &gmaxload, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&load, &gminload, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  if (proc == 0) {
    printf("DRIVER EVAL:  load:  max %f  min %f  sum %f\n", 
           gmaxload, gminload, gsumload);
    printf("DRIVER EVAL:  objs:  max %d  min %d  sum %d\n", 
           gmaxelems, gminelems, gsumelems);
    printf("DRIVER EVAL:  cuts:  max %d  min %d  sum %d\n",
           gmaxcuts, gmincuts, gsumcuts);
  }
}
