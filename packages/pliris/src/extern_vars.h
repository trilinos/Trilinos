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

extern int    me;                       /* processor id information */

extern int    nprocs_cube;		/* num of procs in the allocated cube */
extern int    nprocs_row;		/* num of procs to which a row is assigned */
extern int    nprocs_col;		/* num of procs to which a col is assigned */
extern int    max_procs;		/* max num of procs in any dimension */

extern int    nrows_matrix;		/* number of rows in the matrix */
extern int    ncols_matrix;		/* number of cols in the matrix */
extern int    matrix_size;		/* order of matrix=nrows_matrix=ncols_matrix */

extern int    my_first_row;		/* proc position in a row */
extern int    my_first_col;		/* proc position in a col */

extern int    my_rows;			/* num of rows I own */
extern int    my_cols;			/* num of cols I own */

extern int    my_cols_last;
extern int    my_cols_seg;              /* vars for compatibility with init */
extern int    ncols_last;

extern int    nrhs;                     /* number of right hand sides in the matrix */
extern int    my_rhs;                   /* number of right hand sides that I own */

extern int    mat_stride;               /* stride to second dimension of mat */
extern int    col1_stride;              /* stride to second dimension of col1 */
extern int    row1_stride;              /* stride to second dimension of row1 */

extern DATA_TYPE *mat;			/* incore storage for col being factored */
extern DATA_TYPE *rhs;                 /* storage for right hand sides */

extern DATA_TYPE *col1;         	/* ptrs to col used for updating a col */
extern DATA_TYPE *col2;		/* ptr to col received in message buf */
extern DATA_TYPE *row1;		/* ptr to diagonal row */
extern DATA_TYPE *row2;		/* ptr to pivot row */
extern DATA_TYPE *row3;
extern DATA_TYPE *rhs_temp;
                                /* ptr to row used for pivots */
extern int  *pivot_vec;                 /* stores pivot information */

extern int    blksz;			/* block size for BLAS 3 operations */
extern int    rhs_blksz;                /* agglomeration block size for backsolve */
extern int    colcnt;			/* number of columns stored for BLAS 3 ops */


extern int   myrow,mycol;
extern MPI_Comm row_comm,col_comm;





#if defined(Pliris_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Pliris package is deprecated"
#endif
#endif
