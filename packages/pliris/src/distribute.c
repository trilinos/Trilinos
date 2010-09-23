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
    if (*my_col < *nrhs % (*nprocsr)) ++(*my_rhs);




} 




