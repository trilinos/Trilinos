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
#   include <mpi.h>
#include "defines.h"
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
extern int my_rhs;             /* num of right hand side I own */

extern MPI_Comm col_comm;
extern MPI_Comm row_comm;



void
permute_rhs(DATA_TYPE *rhs, int *permutations)
{ int k;                        /* loop counter */


  int col_len,row_len;
 

  MPI_Status msgstatus;

  int pivot_row;

  DATA_TYPE tmpv;      /* Temporary variable for receive operation */

  int rhs_col,k_row;

  
  col_len = my_rows;    /* length of column in remaining local matrix */


  row_len = my_cols + my_rhs;  /* length of row in local matrix including 
			          rhs's*/
 
  rhs_col = 0;
  for (k=0; k<=nrows_matrix-2; k++)
  {
    k_row = k%nprocs_col;
    if (mycol == rhs_col)
    {
	if (myrow  == k_row) 
        {
	  pivot_row = permutations[k/nprocs_col]; 
        }
        MPI_Bcast(&pivot_row,1,MPI_INT,k_row,col_comm);
        if (k != pivot_row)
	{
	  /* Perform necessary sends to do the swap: */
	  if (myrow == k_row) 
          {
	    MPI_Send((char *)(&rhs[k/nprocs_col]),sizeof(DATA_TYPE),MPI_CHAR,
		pivot_row%nprocs_col,2,col_comm);
	  }
	  if (myrow == pivot_row%nprocs_col) 
          {
   	    MPI_Send((char *)(&rhs[pivot_row/nprocs_col]),sizeof(DATA_TYPE),
		MPI_CHAR,k_row,3,col_comm);
	  }
	  /* Perform Receives: */
	  if (myrow == k_row)
	  {
	    MPI_Recv((char *)(&tmpv),sizeof(DATA_TYPE),MPI_CHAR,
		pivot_row%nprocs_col,3,col_comm,&msgstatus);      
    	    rhs[k/nprocs_col] = tmpv;
          }
	  if (myrow == pivot_row%nprocs_col)
	  {
            MPI_Recv((char *)(&tmpv),sizeof(DATA_TYPE),
		MPI_CHAR,k_row,2,col_comm,&msgstatus);
    	    rhs[pivot_row/nprocs_col] = tmpv;
          } 
       }
    }
  }

}/* End of function permute_rhs */
