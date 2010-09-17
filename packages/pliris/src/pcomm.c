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

#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "defines.h"
#include "macros.h"

#define DEBUG1 0   
/*  define variables to avoid compiler error    */

int one = 1;
double d_one = 1.;

int ringnext,ringprev,hbit,rmbit,my_col_id,my_row_id;
int ringnex2,ringpre2,ringnex3,ringpre3,ringnex4,ringpre4;
typedef struct {
  DATA_TYPE entry;
  DATA_TYPE current;
  int row;
} pivot_type;

void initcomm(){
  extern int nprocs_col, nprocs_row, me, hbit, my_col_id, my_row_id, rmbit;
  extern int ringnext,ringprev,ringnex2,ringpre2,ringnex3,ringpre3,ringnex4,ringpre4;
  int col_id,bit;

  my_col_id = mesh_col(me);
  my_row_id = mesh_row(me);

  
  col_id = my_col_id + 1;
  if (col_id >= nprocs_row) col_id = 0;
  ringnext = proc_num(my_row_id,col_id);

  col_id = my_col_id + 2;
  if (col_id >= nprocs_row) col_id -= nprocs_row;
  ringnex2 = proc_num(my_row_id,col_id);

  col_id = my_col_id + 3;
  if (col_id >= nprocs_row) col_id -= nprocs_row;
  ringnex3 = proc_num(my_row_id,col_id);

  col_id = my_col_id + 4;
  if (col_id >= nprocs_row) col_id -= nprocs_row;
  ringnex4 = proc_num(my_row_id,col_id);

  col_id = my_col_id - 1;
  if (col_id < 0) col_id = nprocs_row - 1;
  ringprev = proc_num(my_row_id,col_id);

  col_id = my_col_id - 2;
  if (col_id < 0) col_id += nprocs_row;
  ringpre2 = proc_num(my_row_id,col_id);

  col_id = my_col_id - 3;
  if (col_id < 0) col_id += nprocs_row;
  ringpre3 = proc_num(my_row_id,col_id);

  col_id = my_col_id - 4;
  if (col_id < 0) col_id += nprocs_row;
  ringpre4 = proc_num(my_row_id,col_id);

  /* calculate first power of two bigger or equal to the number of rows,
     and low order one bit in own name*/

  for (hbit = 1; nprocs_col > hbit ; hbit = hbit << 1);

  rmbit = 0;
  for (bit = 1; bit < hbit; bit = bit << 1) {
    if ((my_row_id & bit) == bit) {
      rmbit = bit; break;}
  }

#if (DEBUG1 > 0)
  printf("In initcomm, node %d: my_col_id = %d, my_row_id = %d, hbit = %d, rmbit = %d, ringnext = %d, ringprev = %d\n",me,my_col_id,my_row_id,hbit,rmbit,ringnext,ringprev);
#endif
}

double 
max_all(double buf, int type)
{

    double maxval;

    maxval = buf;

    MPI_Allreduce(&buf,&maxval,1,MPI_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD);

    return maxval;
}

