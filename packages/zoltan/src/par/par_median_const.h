// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __PAR_MEDIAN_CONST_H
#define __PAR_MEDIAN_CONST_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int Zoltan_RB_find_median(int, double *, double *, double, int *,
  int, int, double*, MPI_Comm,
  double *, int, int, int,
  int, int, int, double, double, double,
  double *, double *, int *, int, int);

extern int Zoltan_RB_find_median_randomized(int, double *, double *, double, int *,
  int, int, double*, MPI_Comm,
  double *, int, int, int,
  int, int, int, double, double, double,
  double *, double *, int *, int, int);

/* Prototype for function used with TFLOPS_SPECIAL */
extern void Zoltan_RB_reduce(int, int, int, void*, void*,
                             int, int*, MPI_Datatype, MPI_Comm, 
                             MPI_User_function);

extern void par_median_accumulate_counts(int nprocs, int num_procs, int rank, int count);
extern void par_median_print_counts(MPI_Comm comm, int print_proc);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
