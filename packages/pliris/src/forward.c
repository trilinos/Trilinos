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
#include <stdlib.h>
#include <mpi.h>
#include "defines.h"
#include "BLAS_prototypes.h"
#include "solve.h"
#include "macros.h"
#include "pcomm.h"


extern int myrow;
extern int mycol;

extern int me;

extern int nrhs;	       /* number of rhs's stored in last col's of 
					matrix */
extern MPI_Comm row_comm;
extern MPI_Comm col_comm;

extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */

extern int nprocs_col;		/* num of procs to which a col is assigned */
extern int nprocs_row;		/* num of procs to which a row is assigned */

extern int my_rows;		/* num of rows I own */
extern int my_cols;            /* num of cols I own */

extern int nrhs;                /* number of right hand sides */
extern int my_rhs;              /* number of right hand sides that I own */


#define SOSTATUSINT 32768

void
forward(DATA_TYPE *mat, DATA_TYPE *rhs)
{
  int i, k;               /* loop counters */
  int rhs_col;			/* torus-wrap column containing the rhs */
  int k_row;                    /* torus-wrap row corresponding to kth
                                   global row */
  int k_col;			/* torus-wrap column corresponding to kth 
				   global col */


  int istart;			/* Starting row index for pivot column */
		
  int index; 			
  int count_row;		/* dummy index */ 
  DATA_TYPE *piv_col;		/* portion of pivot column I am sending */
  DATA_TYPE ck;			/* rhs corresponding to current column
				   of the backsubstitution */
 
  MPI_Status msgstatus;

#ifdef COMPLEX
  double tmpr, tmpi;
#endif

  piv_col = (DATA_TYPE *) malloc(my_rows * sizeof(DATA_TYPE));
  /* Perform the Forward Substitution: */
  rhs_col = 0; 
  for (k=0; k<= nrows_matrix-2; k++)
  {
    k_row=k%nprocs_col;
    k_col=k%nprocs_row;
    istart = (k+1-myrow)/nprocs_col;
    if (istart * nprocs_col < k+1-myrow)
	istart++;
    count_row = 0;
    for (i=istart;i<=my_rows-1;i++)
    {
 	index=(k/nprocs_row)*my_rows+i;
 	piv_col[count_row]=mat[index];
        count_row++;
    }
    if (mycol == rhs_col && myrow == k_row)
    {
	ck = rhs[k/nprocs_col];
    }
    if (mycol == k_col)
    {
	MPI_Send((char *)piv_col,count_row*sizeof(DATA_TYPE),
		MPI_CHAR,rhs_col,0,row_comm);
    }
    if (mycol == rhs_col)
    {
	MPI_Recv((char *)piv_col,count_row*sizeof(DATA_TYPE),
		MPI_CHAR,k_col,0,row_comm,&msgstatus);
    } 
    if (mycol == rhs_col)
    {
	MPI_Bcast((char *)(&ck),sizeof(DATA_TYPE),MPI_CHAR,k_row,col_comm);
	count_row=0;
#ifndef COMPLEX
	  for (i=istart;i<=my_rows-1;i++)
          {
		rhs[i] = rhs[i] - piv_col[count_row] * ck;
		count_row++;
	  }
#endif
#ifdef COMPLEX
          for (i=istart;i<=my_rows-1;i++)
          {
		tmpr =  ((piv_col[count_row]).r) * ck.r -
			((piv_col[count_row]).i) * ck.i;
		tmpi =  ((piv_col[count_row]).r) * ck.i + 
			((piv_col[count_row]).i) * ck.r;
		(rhs[i]).r = (rhs[i]).r - tmpr;
		(rhs[i]).i = (rhs[i]).i - tmpi;	
/*                rhs[i] = rhs[i] - piv_col[count_row] * ck; */
                count_row++;
          }
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }/* end of for (k=0; k<= nrows_matrix-2; k++) */
  free(piv_col); 
	

}/* End of function forward */
