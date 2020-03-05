/*
//@HEADER
// ************************************************************************
//
//                        Adelus v. 1.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of NTESS nor the names of the contributors may be
// used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL NTESS OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Vinh Dang (vqdang@sandia.gov)
//                    Joseph Kotulski (jdkotul@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef __ADELUS_MYTIME_HPP__
#define __ADELUS_MYTIME_HPP__

#include <stdio.h>
#include <mpi.h>
#include "Adelus_defines.h"

namespace Adelus {

double get_seconds(double start)
{
  double time;		/* total seconds */
  time = MPI_Wtime();
  time = time - start;
  
  return (time);
}

// Exchange and calculate max, min, and average timing information

void showtime(const char *label, double *value)
{
  extern int me;		// current processor rank
  extern int nprocs_cube;
  
  double avgtime;
  
  struct {
    double val;
    int proc;
  } max_in, max_out, min_in, min_out;
  max_in.val = *value;
  max_in.proc = me;
  MPI_Allreduce(&max_in,&max_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
  min_in.val = *value;
  min_in.proc = me;
  MPI_Allreduce(&min_in,&min_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
  
  MPI_Allreduce(value,&avgtime,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  avgtime /= nprocs_cube;
  
  if (me == 0) {
    fprintf(stderr, "%s = %.4f (min, on proc %d), %.4f (avg), %.4f (max, on proc %d).\n",
      label,min_out.val,min_out.proc,avgtime, max_out.val,max_out.proc);
  }
}

}//namespace Adelus

#endif
