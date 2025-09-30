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

Authors:

Brian Driessen
Sandia National Labs
(505)-844-9297
bjdries@sandia.gov

Joseph D. Kotulski
Sandia National Labs
(505)-845-7955
jdkotul@sandia.gov

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "BLAS_prototypes.h"
double  timing(double secs, int type);
#include "x_factor.h"
#include "mpi.h"
#include "extern_vars.h"
#include "macros.h"
#include "block.h"
#include "factor.h"
#include "pcomm.h"

#define PERMTYPE ((1 << 5) + (1 << 4))


void X_FACTOR_ (DATA_TYPE *matrix,int *matrixsize,
   int *num_procsr, int *permute, double *secs)
{

  DATA_TYPE *mat;
  int *permutations;
  double run_secs;              /* time (in secs) during which the prog ran */
  double seconds(double);       /* function to generate timings */
  double tsecs;                 /* intermediate storage of timing info */

  int totmem1;
/*
   Determine who I am (me ) and the total number of nodes (nprocs_cube)
                                                                        */

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs_cube);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  permutations = permute;
  mat = matrix;

  matrix_size =  *matrixsize;
  nrows_matrix = *matrixsize;
  ncols_matrix = *matrixsize;
  nprocs_row = *num_procsr;

  totmem1=0;                      /* Initialize the total memory used */
  nprocs_col = nprocs_cube/nprocs_row;
  max_procs = (nprocs_row < nprocs_col) ? nprocs_col : nprocs_row;

    /* set up communicators for rows and columns */

    myrow = mesh_row(me);
    mycol = mesh_col(me);
    MPI_Comm_split(MPI_COMM_WORLD,myrow,mycol,&row_comm);
    MPI_Comm_split(MPI_COMM_WORLD,mycol,myrow,&col_comm);

    {int checkcol,checkrow;
     MPI_Comm_rank(col_comm, &checkrow) ;
     MPI_Comm_rank(row_comm, &checkcol) ;
     if (myrow != checkrow) {
       printf("Node %d: my row = %d but rank in col = %d\n",me,myrow,checkrow);     if (mycol != checkcol)
       printf("Node %d: my col = %d but rank in row = %d\n",me,mycol,checkcol);
     }
    }

  /* Distribution for the matrix on me */

  my_first_col = mesh_col(me);
  my_first_row = mesh_row(me);

  my_rows = nrows_matrix / nprocs_col;
  if (my_first_row < nrows_matrix % nprocs_col)
    ++my_rows;
  my_cols = ncols_matrix / nprocs_row;
  if (my_first_col < ncols_matrix % nprocs_row)
    ++my_cols;

  /* blksz paramter must be set */

  blksz = DEFBLKSZ;


  /* allocate arrays for factor/solve */


  pivot_vec = (int *) malloc(my_cols * sizeof(int));
  totmem1 += my_cols * sizeof(int);
  if (pivot_vec == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  row3 = (DATA_TYPE *) malloc((my_cols +1+ blksz + nrhs) * sizeof(DATA_TYPE));
  totmem1 += (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row3 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }


  row2 = (DATA_TYPE *) malloc((my_cols + blksz+10 + nrhs) * sizeof(DATA_TYPE));
  totmem1 += (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  row1_stride = my_cols+blksz+1;
  row1 = (DATA_TYPE *) malloc(blksz*(my_cols+blksz+nrhs+3)*sizeof(DATA_TYPE));
  totmem1 += blksz * (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row1 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  col2 = (DATA_TYPE *) malloc((my_rows + 10) * sizeof(DATA_TYPE));
  totmem1 += (my_rows + 1) * sizeof(DATA_TYPE);
  if (col2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  col1_stride = my_rows;
  col1 = (DATA_TYPE *) malloc(blksz * (my_rows + 10) * sizeof(DATA_TYPE));
  totmem1 += blksz * (my_rows + 1) * sizeof(DATA_TYPE);
  if (col1 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  mat_stride = my_rows;

  /* Factor and Solve the system */

  tsecs = seconds(0.0);
  /* Initialize Communication  */

  initcomm();
  factor(mat);

  tsecs = seconds(tsecs);

  run_secs = (double) tsecs;

  /* Solve time secs */

  *secs = run_secs;

  free(row2);

}


