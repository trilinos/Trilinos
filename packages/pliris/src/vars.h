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

int   me;                       /* processor id information */

int   nprocs_cube;		/* num of procs in the allocated cube */
int   nprocs_row;		/* num of procs to which a row is assigned */
int   nprocs_col;		/* num of procs to which a col is assigned */
int   max_procs;		/* max num of procs in any dimension */

int   nrows_matrix;		/* number of rows in the matrix */
int   ncols_matrix;		/* number of cols in the matrix */
int   matrix_size;		/* order of matrix=nrows_matrix=ncols_matrix */

int   my_first_row;		/* proc position in a row */
int   my_first_col;		/* proc position in a col */

int   my_rows;			/* num of rows I own */
int   my_cols;			/* num of cols I own */

int   my_cols_last;
int   my_cols_seg;              /* vars for compatibility with init */
int   ncols_last;

int   nrhs;                     /* number of right hand sides in the matrix */
int   my_rhs;                   /* number of right hand sides that I own */

int   mat_stride;               /* stride to second dimension of mat */
int   col1_stride;              /* stride to second dimension of col1 */
int   row1_stride;              /* stride to second dimension of row1 */

DATA_TYPE *mat;			/* incore storage for col being factored */
DATA_TYPE *rhs;                 /* storage for right hand sides */

DATA_TYPE *col1;         	/* ptrs to col used for updating a col */
DATA_TYPE *col2;		/* ptr to col received in message buf */
DATA_TYPE *row1;		/* ptr to diagonal row */
DATA_TYPE *row2;		/* ptr to pivot row */
DATA_TYPE *row3;
DATA_TYPE *rhs_temp;
                                /* ptr to row used for pivots */
int *pivot_vec;                 /* stores pivot information */

int   blksz;			/* block size for BLAS 3 operations */
int   rhs_blksz;                /* agglomeration block size for backsolve */
int   colcnt;			/* number of columns stored for BLAS 3 ops */


int  myrow,mycol;
MPI_Comm row_comm,col_comm;


/* volatile int   MSPLIT;           ZGEMM splitting parameter */

