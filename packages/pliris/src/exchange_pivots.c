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
#include <mpi.h>

#include "defines.h"
#include "BLAS_prototypes.h"
#include "macros.h"
#include "pcomm.h"


extern int myrow;
extern int mycol;

extern int me;	               /* processor id information */
extern int nprocs_row;         /* num of procs to which a row is assigned */
extern int nprocs_col;         /* num of procs to which a col is assigned */
extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */
extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */

extern int *pivot_vec;         /* ptr to vector storing list of pivot rows */
extern int mat_stride,col1_stride,row1_stride;  /* strides for 2nd dim of 2d mats */


extern MPI_Comm col_comm;
extern MPI_Comm row_comm;


void
exchange_pivots(int *permutations)

{ int j,k;                        /* loop counter */


  


  int colcnt;        /* number of columns stored for BLAS 3 ops */
  int col_len;
  
 

  MPI_Status msgstatus;


  int i,rank_row,k_row,pivot_col;

  colcnt = 0;           /* number of column's currently saved for update */
  col_len = my_rows;    /* length of column in remaining local matrix */

/*  First gather the permutation vector to processor 0 in row_comm   */
 
  if (myrow == 0 || mycol == 0)
  {
	for (k=0;k<=nrows_matrix-1;k++)
	{
		pivot_col = k%nprocs_row;
		k_row = k%nprocs_col;
		rank_row = k_row * nprocs_row;
		if (me == pivot_col)
		{
			j=k/nprocs_row;
			MPI_Send(&pivot_vec[j],1,MPI_INT,
				rank_row,0,MPI_COMM_WORLD);
		}
		if (me == rank_row)
		{
			i=k/nprocs_col;
			MPI_Recv(&permutations[i],1,MPI_INT,pivot_col,0,
				MPI_COMM_WORLD,&msgstatus);
		}
	}
  }	
  MPI_Barrier(MPI_COMM_WORLD);
  /*   Broadcast to the rest of the processors  in row_comm     */
  MPI_Bcast(permutations,my_rows,MPI_INT,0,row_comm);
}/* End of function exchange_pivots*/
