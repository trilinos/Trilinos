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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "defines.h"
#include "BLAS_prototypes.h"
#include "macros.h"
#include "pcomm.h"
#include "permute_mat.h"
#include "permute_rhs.h"
#include "forward.h"
#include "exchange_pivots.h"

extern int myrow;
extern int mycol;

extern int me;	               /* processor id information */
extern int nprocs_row;         /* num of procs to which a row is assigned */
extern int nprocs_col;         /* num of procs to which a col is assigned */
extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */
extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */
extern int my_rhs;             /* num of right hand side I own */

extern DATA_TYPE *col1;           /* ptrs to col used for updating a col */
extern DATA_TYPE *row1;           /* ptr to diagonal row */
extern DATA_TYPE *row2;           /* ptr to pivot row */
extern DATA_TYPE *row3;           /* ptr to temporary vector for rows */
extern int *pivot_vec;         /* ptr to vector storing list of pivot rows */
extern int mat_stride,col1_stride,row1_stride;  /* strides for 2nd dim of 2d mats */

extern int blksz;              /* block size for BLAS 3 operations */
extern int ringnext,ringprev,ringnex2,ringpre2,ringnex3,ringpre3,ringnex4,ringpre4;
#define LUSTATUSINT 64

extern MPI_Comm col_comm;
extern MPI_Comm row_comm;

#define LUPIVOTTYPE (1<<19)
#define LUCOLTYPE (1<<20)
#define LUROWTYPE (1<<21)
#define LUPIVROWTYPE ((1<<21) + (1<<20))
#define LUSENDTYPE (1<<22)

#define rowplus(I) proc_num(mesh_row(me),(mesh_col(me)+(I) < nprocs_row) ? mesh_col(me)+(I) : mesh_col(me)+(I)-nprocs_row)
#define rowminus(I) rowplus(nprocs_row - (I))
#define MAXDIST 1

void
forwardall(DATA_TYPE *mat, int *permutations, DATA_TYPE *rhs_copy, 
        DATA_TYPE *rhs)
{
  int colcnt; /* number of columns stored for BLAS 3 ops */
 
  int col_len,row_len,rows_used,cols_used;
  int *sav_pivot_ptr;
  

  DATA_TYPE *sav_col_ptr,*sav_row_ptr,*sav_piv_row_ptr;
  DATA_TYPE *cur_col1_row_ptr;
  DATA_TYPE *temp_row_ptr;
  DATA_TYPE *act_col_ptr,*act_row_ptr;
 
  
  int numprocs;

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
  if ( (numprocs/nprocs_row) * nprocs_row != numprocs )
  {
     if (me == 0)
     {
       printf("nprocs_row must go into numprocs perfectly!\n");
       printf("Try a different value of nprocs_row.\n");
     }
     MPI_Barrier(MPI_COMM_WORLD);
     exit(0);
  }
  colcnt = 0;           /* number of column's currently saved for update */
  col_len = my_rows;    /* length of column in remaining local matrix */
  sav_col_ptr = col1;   /* location to store next active column */
  act_col_ptr = col1;   /* location of matrix of columns being saved for dgemm update */

  row_len = my_cols + my_rhs;  /* length of row in local matrix including 
			          rhs's*/

  rows_used = 0;      /* haven't used any local rows yet */
  cols_used = 0;
  cur_col1_row_ptr = col1;  /* location of first row in col1 matrix */
  act_row_ptr = row1; /* location of matrix of rows being saved for dgemm 
			 update */

  sav_piv_row_ptr = row1; /* location for next row being saved for dgemm 
	                 update */

  temp_row_ptr = row3; /* location for pivot row while being sent and 
		         before transposing */
 
  sav_row_ptr = row2;  /* location to save current row and send to 
			 owner of pivot row */

  sav_pivot_ptr = pivot_vec; /* location to store name of pivot row */

 
  exchange_pivots(permutations);
  permute_mat(mat,permutations);
  permute_rhs(rhs,permutations);
  forward(mat, rhs);
  MPI_Barrier(MPI_COMM_WORLD);

}/* End of function forwardall */
