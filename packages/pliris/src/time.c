/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdio.h>
#include <mpi.h>
#include "defines.h"
double  seconds(double );
double  timing(double , int );
#include "mytime.h"

double
seconds(double start)
{
    double time;		/* total seconds */
    time = MPI_Wtime();
    time = time - start;

    return (time);
}

/*
** Exchange and calculate average timing information
**
** secs:        number of ticks for this processor
** type:        type of message to collect
*/

double
timing(double secs, int type)
{

    extern int me;		/* current processor number */
    extern int nprocs_cube;
    double avgtime;


    struct {
      double val;
      int proc;
    } max_in, max_out;

    max_in.val = secs;
    max_in.proc = me;
    MPI_Allreduce(&max_in,&max_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);

    MPI_Allreduce(&secs,&avgtime,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    avgtime /= nprocs_cube;

    if (me == 0) {
	fprintf(stderr, "%.4f (avg), %.4f (max on processor %d).\n",
		avgtime, max_out.val, max_out.proc);
    }

    return avgtime;
}

void showtime(char *label, double *value)
{

    extern int me;		/* current processor number */
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
	fprintf(stderr, "%s = %.4f (on proc %d), %.4f, %.4f (on proc %d).\n",
		label,min_out.val,min_out.proc,avgtime, max_out.val,max_out.proc);
    }
  }

