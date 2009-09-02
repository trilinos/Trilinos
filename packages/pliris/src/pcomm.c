/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
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

