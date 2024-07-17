/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef __ADELUS_MYTIME_HPP__
#define __ADELUS_MYTIME_HPP__

#include <stdio.h>
#include <mpi.h>
#include "Adelus_defines.h"

namespace Adelus {

inline
double get_seconds(double start)
{
  double time;		/* total seconds */
  time = MPI_Wtime();
  time = time - start;

  return (time);
}

// Exchange and calculate max, min, and average timing information
inline
void showtime(int comm_id, MPI_Comm comm, int me, int nprocs_cube, const char *label, double *value)
{
  double avgtime;

  struct {
    double val;
    int proc;
  } max_in, max_out, min_in, min_out;
  max_in.val = *value;
  max_in.proc = me;
  MPI_Allreduce(&max_in,&max_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,comm);
  min_in.val = *value;
  min_in.proc = me;
  MPI_Allreduce(&min_in,&min_out,1,MPI_DOUBLE_INT,MPI_MINLOC,comm);

  MPI_Allreduce(value,&avgtime,1,MPI_DOUBLE,MPI_SUM,comm);

  avgtime /= nprocs_cube;

  if (me == 0) {
    fprintf(stderr, "Communicator %d -- %s = %.4f (min, on proc %d), %.4f (avg), %.4f (max, on proc %d).\n",
      comm_id,label,min_out.val,min_out.proc,avgtime, max_out.val,max_out.proc);
  }
}

}//namespace Adelus

#endif
