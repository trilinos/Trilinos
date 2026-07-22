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
 * Print out topology and which MPI processes are where.
 * Printing out the MPI ranks may not be possible 
 * if hwloc_get_cpubind() doesn't work or if it
 * works but gives us a cpuset containing more than one cpu.
 *
 * Here we are assuming homogeneous topologies.  At any level,
 * each child of that level has identical topology.
 *
 * When running this code make sure the path to the hwloc shared library
 * is in your environment.
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

#define MAX_NAME_LEN 64           /* what should this be? */
#define MAX_NUM_LEVELS 12          /* ditto */

static char *which_mpi();

int main(int argc, char *argv[])
{
int have_my_cpu=1;                /* We know the one cpu that each MPI process is on */
char *recvbuf=NULL;
char type_name[MAX_NAME_LEN];    /* "machine", "node", "socket", ... */
char mask[MAX_NAME_LEN];
hwloc_topology_t topology;
const struct hwloc_topology_support *support=NULL;
hwloc_obj_t obj;
hwloc_cpuset_t *cpuset=NULL;
hwloc_cpuset_t binding;
uint64_t local_memory;
uint64_t total_memory;
int size, rank, depth;
int i, j, p;
int rc, num;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(zoltan_get_global_comm(), &size);
  MPI_Comm_rank(zoltan_get_global_comm(), &rank);

  /* allocate & initialize topology object */

  hwloc_topology_init(&topology);

#if 0
  /* configure detection - skip levels with only one child or which are only children */
  hwloc_topology_ignore_all_keep_structure(topology);
#endif

  /* detect topology */

  hwloc_topology_load(topology);

  /* capabilities of hwloc on this machine */

  support = hwloc_topology_get_support((struct hwloc_topology *)topology);

  /* 
   * Get the depth of topology and the root.
   * The root typically is a node on a multi-node machine, not the collection of nodes
   * in the application.
   */

  depth = hwloc_topology_get_depth(topology);

  /* 
   * Which cpu am I running on?
   *
   * HWLOC_CPUBIND_STRICT - says assume each process is running on one processor and won't be moved
   *   (may ensure we get a cpuset containing no more than one cpu - doesn't work on glory)
   */

  binding = hwloc_cpuset_alloc();
  hwloc_cpuset_zero(binding);

  if (support->cpubind->get_thisproc_cpubind != 1){
    have_my_cpu = 0;
  }
  else{
    rc = hwloc_get_cpubind(topology, binding, HWLOC_CPUBIND_STRICT);
  
    if ((rc < 0) || (hwloc_cpuset_weight(binding) > 1)){
      have_my_cpu = 0;
    }
  }

  MPI_Barrier(zoltan_get_global_comm());

  if (!have_my_cpu && (rank == 0)){
    printf("Warning: Unable to identify each MPI process with its CPU\n");
  }

  if (rank == 0){
    recvbuf = (char *)calloc(size , MAX_NAME_LEN);
  }
  else{
    recvbuf = NULL;
  }

  if (have_my_cpu){

    hwloc_cpuset_snprintf(mask, MAX_NAME_LEN-1, binding);
    MPI_Gather(mask, MAX_NAME_LEN, MPI_CHAR, recvbuf, MAX_NAME_LEN, MPI_CHAR, 0, zoltan_get_global_comm());
 
    if (rank == 0){
      cpuset = (hwloc_cpuset_t*)malloc(sizeof(hwloc_cpuset_t) * size);

      for (p=0; p < size; p++){
        cpuset[p] = hwloc_cpuset_alloc();
        hwloc_cpuset_from_string(cpuset[p], recvbuf + (p * MAX_NAME_LEN));
      }
    }

  }

  if (rank == 0){

    /* Get topology information.
     * Which of the MPI processes are each socket, cache, etc?
     */

    for (i=0; i < depth; i++){

      num = hwloc_get_nbobjs_by_depth(topology, i);
   
      for (j = 0; j < num; j++){

        obj = hwloc_get_obj_by_depth(topology, i, j);

        if (j==0){
          hwloc_obj_type_snprintf(type_name, MAX_NAME_LEN-1, obj, 1);
          printf("\n%d %s%s:\n",num,type_name, ((num> 1) ? "s" : "")); 
        }

        hwloc_cpuset_snprintf(mask, MAX_NAME_LEN - 1, obj->cpuset); 
        local_memory = obj->memory.local_memory;
        total_memory= obj->memory.total_memory;

        printf("  cpu set mask %s,  total memory %10.0fKB, local memory %10.0fKB",
              mask, (float)total_memory/1024.0, (float)local_memory/1024.0);

        if (have_my_cpu){
          printf(", ( MPI ranks ");
          for (p=0; p < size; p++){
            if (hwloc_cpuset_isincluded(cpuset[p], obj->cpuset)){
              printf("%d ",p);
            }
          }
        }
        printf(")\n");
      }
    }

    if (cpuset){
      for (p=0; p < size; p++){
        hwloc_cpuset_free(cpuset[p]);
      }
      free(cpuset);
    }
    free(recvbuf);
    printf("\n%s\n",which_mpi());
  }

  hwloc_cpuset_free(binding);

  MPI_Finalize();

  return 0;
}

/* Curious - do different mpi implementations lay out processes differently
 * on tiered multicore architectures?
 */
static char *which_mpi()
{
char *info=NULL;

  info = (char *)calloc(64, 1);

#ifdef OMPI_MAJOR_VERSION
  sprintf(info,"Open MPI %d %d %d", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION);
#endif

#ifdef LAM_MAJOR_VERSION
  sprintf(info,"LAM/MPI %d %d %d", LAM_MAJOR_VERSION, LAM_MINOR_VERSION, LAM_RELEASE_VERSION);
#endif

#ifdef MPICH2_VERSION
  sprintf(info,"MVAPICH2 or MPICH2 %s", MPICH2_VERSION);
#endif

#ifdef MPICH_VERSION
  sprintf(info,"MVAPICH or MPICH %s", MPICH_VERSION);
#endif

  if (info[0] == 0){
    sprintf(info,"unidentified");
  }

  return info;
}

