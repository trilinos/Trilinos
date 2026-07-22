// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/* 
 * Discover hardware topology using hwloc.
 * Tested with hwloc 1.0.2: http://www.open-mpi.org/software/hwloc
 *
 * This function returns hierarchy on a node or on a machine.
 * A dual socket quad core might return the following:
 *
 * depth = Zoltan_Get_Toplogy(branch_degrees, num_cpus, names, tmem, lmem)
 *
 * depth = 3
 * degrees = {2, 4, 0}
 * num_cpus = {8, 4, 1}
 * names = {"MACHINE", "SOCKET", "PU"}
 * tmem = {34359738368, 17179869184, 4294967296}
 * lmem = {34359738368, 17179869184, 4294967296}
 *
 * Because a machine (a node to us) contains 2 sockets with 8 total cpus.
 * A socket contains 4 cores with 4 total cpus.
 * A processing unit is always the "leaf node" of the topology.
 *
 * tmem is the total memory of that object and its children
 * lmem is the total memory of that object
 *
 * If memory information is not available it will be zero.
 *
 * If any parameter is NULL it will be ignored.
 */

#include <hwloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <mpi.h>

#define MEMORY_ERROR(buf) { \
  if (buf == NULL){ \
    fprintf(stderr,"Memory error\n"); \
    exit(1); \
  } \
}

/* 
 * Define BUILD_MAIN to create a test application
 */
#define BUILD_MAIN 

#define MAX_NAME_LEN 64           /* what should this be? */
#define MAX_NUM_LEVELS 12          /* ditto */

int Zoltan_Get_Topology(int **branch_degree, 
                        int **num_cpus, 
                        char ***level_name, 
                        uint64_t **total_memory, 
                        uint64_t **local_memory)
{
hwloc_topology_t topology;
hwloc_obj_t rootobj, obj, prev_obj;
char description[4096];
int size, rank, depth;
int i, j;

int real_depth, next, from, to, num;
int real_level[MAX_NUM_LEVELS];

char type_name[MAX_NUM_LEVELS][MAX_NAME_LEN];    /* "machine", "node", "socket", ... */
int  type_size[MAX_NUM_LEVELS];                  /* how many (dual node is 2, quad socket is 4) */
hwloc_cpuset_t ancestor_cpuset[MAX_NUM_LEVELS];  /* mask representing cpus in this object */
uint64_t lmem[MAX_NUM_LEVELS];           /* memory owned by this object */
uint64_t tmem[MAX_NUM_LEVELS];           /* owned by this object and all its children */

  MPI_Comm_size(zoltan_get_global_comm(), &size);
  MPI_Comm_rank(zoltan_get_global_comm(), &rank);

  /* allocate & initialize topology object */

  hwloc_topology_init(&topology);

#if 0
  /* Configure detection to skip levels with only one child or which are only children.
   *
   * If we only want to know branch degrees, we can skip the levels in the hierarchy
   * that don't branch, but for now it's interesting to see what they are.
   */

  hwloc_topology_ignore_all_keep_structure(topology);
#endif

  /* detect topology */

  hwloc_topology_load(topology);

  /* 
   * Get the depth of topology and the root.
   */

  depth = hwloc_topology_get_depth(topology);

  for (i=0; i < depth; i++){
    ancestor_cpuset[i] = hwloc_cpuset_alloc();
  }

  /* 
   * Get topology information.  We'll use sibling 0 at each level of the topology.
   *
   * The root object typically is a node on a multi-node machine, not the 
   * collection of nodes in the application.
   */

  rootobj = hwloc_get_root_obj(topology);

  hwloc_obj_type_snprintf(type_name[0], MAX_NAME_LEN, rootobj, 1);
  hwloc_cpuset_copy(ancestor_cpuset[0], rootobj->cpuset); 
  type_size[0] = 1;
  lmem[0] = rootobj->memory.local_memory;  /*0 if not available */
  tmem[0] = rootobj->memory.total_memory;

  real_level[0] = 0;
  next = 1;

  prev_obj = rootobj;

  for (i=1; i < depth; i++){

    type_size[i] = prev_obj->arity;      /* parent's number of children */

    if ( (type_size[i] > 1) || (i == depth-1)){
      /* 
       * Significant levels are processing units (the leaf nodes) and levels
       * that are genuine branches - not the only child of parent.
       */
      real_level[next++] = i;
    }
  
    obj = prev_obj->children[0];
    hwloc_obj_type_snprintf(type_name[i], MAX_NAME_LEN, obj, 1);
    hwloc_cpuset_copy(ancestor_cpuset[i], obj->cpuset);
    lmem[i] = obj->memory.local_memory;
    tmem[i] = obj->memory.total_memory;

    prev_obj = obj;
  }

  real_depth = next;

  if (branch_degree){
    *branch_degree = (int *)malloc(sizeof(int) * real_depth);
    MEMORY_ERROR(*branch_degree);
  }
  if (num_cpus){
    *num_cpus = (int *)malloc(sizeof(int) * real_depth);
    MEMORY_ERROR(*num_cpus);
  }
  if (level_name){
    *level_name = (char **)malloc(sizeof(char*) * real_depth);
  MEMORY_ERROR(*level_name);
  }
  if (total_memory){
    *total_memory = (uint64_t *)malloc(sizeof(uint64_t) * real_depth);
    MEMORY_ERROR(*total_memory);
  }
  if (local_memory){
    *local_memory = (uint64_t *)malloc(sizeof(uint64_t) * real_depth);
    MEMORY_ERROR(*local_memory);
  }


  for (i=0; i < real_depth; i++){
    from = real_level[i];

    if (i < real_depth-1) /* branch node */
      to = real_level[i+1];
    else
      to = 0;             /* leaf node */
 
    /* Total number of CPUs within this object
     */

    if (num_cpus){
      (*num_cpus)[i] = hwloc_cpuset_weight(ancestor_cpuset[from]);
    }

    /* number of next level down objects in one of these objects
     */

    if (branch_degree){
      if (to > 0){
        (*branch_degree)[i] = type_size[to];
      }
      else{
        (*branch_degree)[i] = 0;
      }
    }

    /* 
     * name of this level in the topology and the succession of children
     * level that have only one parent.
     */

    if (level_name){
      description[0] = 0;
      if (to){
        num = type_size[from];
        for (j=from; j < to; j++){
          strcat(description, type_name[j]);
          if (num > 1){
             strcat(description,"s");
          }
          if (j < (to - 1)){
             strcat(description, "/");
          }
        }
      }
      else{
        strcat(description, type_name[from]);
      }
  
      (*level_name)[i] = (char *)malloc(strlen(description) +1);
      MEMORY_ERROR((*level_name)[i]);
      strcpy((*level_name)[i], description);
    }

    if (local_memory){
      (*local_memory)[i] = lmem[from];
      if ((lmem[from] == 0) && (to > from)){
        for (j=from+1; j < to; j++){
          if (lmem[j] > 0){
            (*local_memory)[j] = lmem[from];
            break;
          }
        }
      }
    }

    if (total_memory){
      (*total_memory)[i] = tmem[from];
      if ((tmem[from]) == 0 && (to > from)){
        for (j=from+1; j < to; j++){
          if (tmem[j] > 0){
            (*total_memory)[j] = tmem[from];
            break;
          }
        }
      }
    }

  }

  for (i=0; i < depth; i++){
    hwloc_cpuset_free(ancestor_cpuset[i]);
  }

  return real_depth;
}

#ifdef BUILD_MAIN

int main(int argc, char *argv[])
{
int rank, depth;
int i, go;
int *branching_degree=NULL, *num_cpus=NULL;
uint64_t *local_memory=NULL, *total_memory=NULL;
char **type_name=NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(zoltan_get_global_comm(), &rank);

  depth = Zoltan_Get_Topology(&branching_degree, &num_cpus, &type_name, &total_memory, &local_memory);

  MPI_Barrier(zoltan_get_global_comm());

  if (rank == 0){


    printf("A %s with %d CPUs has:\n",type_name[0],num_cpus[0]);

    for (i=1; i < depth; i++){
      if (i < depth-1){
        printf("  %d %s (%d total CPUs) %swith\n",
              branching_degree[i-1], type_name[i],num_cpus[i],
              ((num_cpus[i] > 1) ? "each " : ""));
      } 
      else{
        printf("  %d %s\n", branching_degree[i-1], type_name[i]);
      }
    }

    printf("\n");

    go = 0;
    for (i=0; i <depth; i++){
      if ((total_memory[i] > 0) || (local_memory[i] > 0 )){
        go = 1;
        break;
      }
    }

    if (go){
      printf("Total memory (for all children), Local memory (for object at level)\n");
      for (i=0; i < depth; i++){
        printf("Memory at level %s: (%f10.0KB, %f10.0KB)\n",
         type_name[i],(float)total_memory[i]/1024.0,(float)local_memory[i]/1024.0);
      } 
    }
    else{
        printf("Memory available at each level is not available\n");
    }
  }

  /* Allocated cpusets can be freed with hwloc_cpuset_free() */

  MPI_Barrier(zoltan_get_global_comm());
  MPI_Finalize();

  return 0;
}
#endif
