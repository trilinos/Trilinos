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

/*  Code to distribute the matrix for a parallel solve  */

/*
Author:

Joseph D. Kotulski
Sandia National Labs
(505)-845-7955
jdkotul@sandia.gov

*/
/*  Variables  INPUT
               nprocsr --- number of processors assigned to a row
               ncols   --- number of columns(=rows) for the matrix
               nrhs    --- number of right hand sides

               OUTPUT
               my_rows  --- number of rows of the total matrix I own
               my_cols  --- number of columns of the total matrix I own
               my_first_row --- global number of my first row
               my_first_col --- global number of my first column
               my_rhs   --- number of right hand sides that I own
               my_row   --- my subblock of the matrix
               my_col   --- my subblock of the matrix
                                                                     */

#include <mpi.h>
#include "distribute.h"

void distmat_( 
                int *nprocsr,
                int *ncols,
                int *nrhs,
                int *my_rows,
                int *my_cols,
                int *my_first_row,
                int *my_first_col,
                int *my_rhs,
                int *my_row,
                int *my_col) 
{

    int me,nprocs;
    int nprocs_col, nrows;
    int nprocs_row;

/*  Determine who I am and the number of processors that are being used    */

    MPI_Comm_rank(MPI_COMM_WORLD, &me) ;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    nrows = *ncols;
    
    nprocs_row = *nprocsr;

    nprocs_col = nprocs/(*nprocsr) ;

    /* Distribute the rows and columns   */

    *my_row = me/(*nprocsr);
    *my_col = me %(nprocs_row);  


    *my_rows = nrows / nprocs_col;

    *my_first_row = (*my_row)*(*my_rows) +1;

    *my_first_row = ((*my_row) > (nrows%nprocs_col)) ? *my_first_row + (nrows%nprocs_col) :
       *my_first_row + (*my_row);

    if (*my_row < nrows % nprocs_col)
        ++(*my_rows);

    *my_cols = nrows / nprocs_row;

    *my_first_col = (*my_col)*(*my_cols) + 1;

    *my_first_col = ((*my_col) > (nrows%nprocs_row)) ? *my_first_col + (nrows%nprocs_row) :
       *my_first_col + (*my_col);

    *my_cols = *ncols / *nprocsr;
    if (*my_col < *ncols % (*nprocsr))
        ++(*my_cols);

    /* Distribute the RHS per processor */

    *my_rhs = *nrhs / *nprocsr;
    if (*my_first_col < *nrhs % (*nprocsr)) ++(*my_rhs);




} 




