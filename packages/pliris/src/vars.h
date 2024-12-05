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

#ifndef __VARSH__
#define __VARSH__

#if defined(Pliris_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Pliris package is deprecated"
#endif
#endif

#include "defines.h"

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

#endif
