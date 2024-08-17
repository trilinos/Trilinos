// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "comm.h"
#include <pthread.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif
static pthread_mutex_t zoltan_global_mpi_lock;
static MPI_Comm Zoltan_Global_MPI_Comm = MPI_COMM_WORLD; // CHECK: ALLOW MPI_COMM_WORLD

/* Function to set the default communicator */
void zoltan_initialize_global_comm(MPI_Comm comm) {
  pthread_mutex_lock(&zoltan_global_mpi_lock);
  Zoltan_Global_MPI_Comm = comm;
  pthread_mutex_unlock(&zoltan_global_mpi_lock);
}

/* Function to get the default communicator */
MPI_Comm zoltan_get_global_comm() {
  pthread_mutex_lock(&zoltan_global_mpi_lock);
  MPI_Comm comm = Zoltan_Global_MPI_Comm;
  pthread_mutex_unlock(&zoltan_global_mpi_lock);
  return comm;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
