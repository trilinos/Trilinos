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




