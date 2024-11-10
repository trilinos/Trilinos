// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*#define DEBUG_ME */

#define NUM_LEVELS 4

MPI_Comm comm[NUM_LEVELS] = {
zoltan_get_global_comm(),
MPI_COMM_NODE,
MPI_COMM_SOCKET,
MPI_COMM_CACHE
};

char *commName[NUM_LEVELS] = { "world", "node", "socket", "cache"};

int num_significant_levels=0;
int level_part_range[NUM_LEVELS][2];
int level_number[NUM_LEVELS];
int my_part_number;

int *buf=NULL;

void level_down(int level, int *classes, int *me, int *nprocs);
void new_part_range(int *classes, int *rangeCurrent, int *rangeNext, int me, int nprocs);

int main(int argc, char *argv[])
{
int i, j, k, lval, gval, last_level, numClasses;
int nextPart, rc, len;
int me[NUM_LEVELS], nprocs[NUM_LEVELS];
int *classes=NULL;
char *errorstr;

  MPI_Init(&argc, &argv);

  for (i=0; i < NUM_LEVELS; i++){
    rc = MPI_Comm_size(comm[i], nprocs + i);
    if (rc != MPI_SUCCESS){
      MPI_Error_string(rc, errorstr, &len);
      fprintf(stderr,"(%d) MPI_Comm_size %s : %s\n",me[0],commName[i],errorstr); 
    }
    rc = MPI_Comm_rank(comm[i], me + i);
    if (rc != MPI_SUCCESS){
      MPI_Error_string(rc, errorstr, &len);
      fprintf(stderr,"(%d) MPI_Comm_size %s : %s\n",me[0],commName[i],errorstr); 
    }
#ifdef DEBUG_ME
    MPI_Barrier(zoltan_get_global_comm());
    printf("(%d) %s communicator, size %d, my rank %d\n",me[0],commName[i],nprocs[i],me[i]);
    MPI_Barrier(zoltan_get_global_comm());
#endif
  }

  buf = (int *)malloc(sizeof(int) * nprocs[0]);
  classes = (int *)malloc(sizeof(int) * nprocs[0]);

  level_part_range[0][0] = 0;
  level_part_range[0][1] = nprocs[0]-1;

  for (i=0; i < NUM_LEVELS-1; i++){
    /*
     * classes[k] contains the rank (in level i) of the rank 0 element of element k's subcommunicator
     */
    level_down(i, classes, me, nprocs); 

    /*
     * my sub communicator will create which parts in the final partitioning?
     */
    new_part_range(classes, level_part_range[i], level_part_range[i+1], me[i], nprocs[i]);
  }

#ifdef DEBUG_ME
  MPI_Barrier(zoltan_get_global_comm());
  printf("(%d) ranges: %s (%d %d iam %d) %s (%d %d iam %d) %s (%d %d iam %d) %s (%d %d iam %d) \n",
  me[0],
  commName[0], level_part_range[0][0],  level_part_range[0][1], level_part_range[0][0] + me[0],
  commName[1], level_part_range[1][0],  level_part_range[1][1], level_part_range[1][0] + me[1],
  commName[2], level_part_range[2][0],  level_part_range[2][1], level_part_range[2][0] + me[2],
  commName[3], level_part_range[3][0],  level_part_range[3][1], level_part_range[3][0] + me[3]);
  MPI_Barrier(zoltan_get_global_comm());
#endif

  /* Figure out which levels are significant */

#ifdef COLLAPSE_LEVELS
  num_significant_levels = 1;
  level_number[0] = 0;
  gval = 0;

  for (i=0; i < NUM_LEVELS-1; i++){
    if ((nprocs[i+1] == nprocs[i]) || (nprocs[i+1] == 1)){
      lval=1;   /* insignificant */
    }
    else{
      lval=0;   /* meaningful level in hierarchy */
    }

    MPI_Allreduce(&lval, &gval, 1, MPI_INT, MPI_SUM, zoltan_get_global_comm());

    if (gval < nprocs[0]){
      /* next level in hierarchy is significant for at least some processes */
      level_number[num_significant_levels++] = i + 1;
    }
  }
#else
  num_significant_levels = NUM_LEVELS;
  for (i=0; i < NUM_LEVELS; i++){
    level_number[i] = i;
  }
#endif

#ifdef DEBUG_ME
  MPI_Barrier(zoltan_get_global_comm());
  printf("(%d) %d significant levels %d %d %d %d \n",me[0], num_significant_levels,
  level_number[0], level_number[1], level_number[2], level_number[3]);
  MPI_Barrier(zoltan_get_global_comm());
#endif

  last_level = level_number[num_significant_levels-1];

  my_part_number = level_part_range[last_level][0] + me[last_level];

  /* Visualize hierarchy */

  MPI_Gather(&my_part_number, 1, MPI_INT, buf, 1, MPI_INT, 0, zoltan_get_global_comm());

  if (me[0] == 0){

    printf("\nPart numbers in process rank order:\n    ");
    for (i=0; i < nprocs[0]; i++){
      printf("%d ",buf[i],i);
      if (i && (i % 25 == 0)) printf("\n  ");
    }
    printf("\n");

    printf("Parts in each subcommunicator:\n");

    for (i=0; i < num_significant_levels; i++){
      printf("Level %s:\n", commName[level_number[i]]);
      MPI_Gather(&(level_part_range[level_number[i]][0]), 1, MPI_INT, buf, 1, MPI_INT, 0, zoltan_get_global_comm());
      memset(classes, 0, sizeof(int) * nprocs[0]);
      for (j=0; j < nprocs[0]; j++){
        classes[buf[j]]++;
      }
      memset(buf, 0, sizeof(int) * nprocs[0]);
      for (j=0, numClasses=0; j < nprocs[0]; j++){
        if (classes[j] > 0){
          buf[numClasses++] = classes[j];   /* number of parts/procs in sub communicator */
        }
      }
      for (j=0, nextPart=0; j < numClasses; j++){
        printf("    [ ");
        for (k=0; k < buf[j]; k++){
          printf("%d ",nextPart++);
        }
        printf("]\n");
      }
      printf("\n");
    }
  }
  else{
    for (i=0; i < num_significant_levels; i++){
      MPI_Gather(&(level_part_range[level_number[i]][0]), 1, MPI_INT, buf, 1, MPI_INT, 0, zoltan_get_global_comm());
    }
  }

  free(buf);
  free(classes);

  MPI_Finalize();
}

/*
 * Determine which processes in my communicator are in which communicators in
 * the next level down. (Communicators are equivalent to topological entities.)
 */
void level_down(int level, int *classes, int *me, int *nprocs)
{
int i;

  MPI_Group current, next;

  MPI_Comm_group(comm[level], &current);
  MPI_Comm_group(comm[level+1], &next);

  /* Rank zero in each sub communicator tells Rank zero at the current level
   * which processes are in its communicator.
   */

  if (me[level+1] > 0){
    for (i=0; i < nprocs[level]; i++){
      buf[i] = -1;
    }
  }
  else{
    for (i=0; i < nprocs[level+1]; i++){
      buf[i] = i;
    }
    MPI_Group_translate_ranks(next, nprocs[level+1], buf, current, classes);

#ifdef DEBUG_ME
    int j;
    for (j=0; j < nprocs[level+1]; j++){
      printf("(%d) Next level's (%s) rank %d is current level's (%s) rank %d\n",me[0],commName[level],buf[j],commName[level+1],classes[j]);
    }
#endif

    for (i=0; i < nprocs[level]; i++){
      buf[i] = -1;
    }

    for (i=0; i < nprocs[level+1]; i++){
      buf[classes[i]] = me[level];   /* mark procs in my sub communicator */
    }
  }

  MPI_Allreduce(buf, classes, nprocs[level], MPI_INT, MPI_MAX, comm[level]);
}

void new_part_range(int *classes, int *rangeCurrent, int *rangeNext, int me, int nprocs)
{
int i, numClasses, myClass, myClassLeader, rangeStart;

  /* how many equivalence classes? how large are they? */

  numClasses = 0;
  memset(buf, 0, sizeof(int) * nprocs);

  for (i=0; i < nprocs; i++){
    buf[classes[i]]++;
    if (buf[classes[i]] == 1) numClasses++;
  }

  /* which class am I in? what is the narrower range that my part number will be in? */

  myClassLeader = classes[me];
  rangeStart = 0;

  for (i=0; i < nprocs; i++){
    if ((buf[i] > 0) && (i == myClassLeader)) break;
    rangeStart += buf[i];
  } 

  rangeNext[0] = rangeCurrent[0] + rangeStart;
  rangeNext[1] = rangeNext[0] + buf[myClassLeader] - 1;
}

