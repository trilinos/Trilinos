/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

//  Code to distribute the matrix for a parallel solve

//Author:
//
//Joseph D. Kotulski
//Sandia National Labs
//(505)-845-7955
//jdkotul@sandia.gov

//  Variables  INPUT
//             comm    --- communicator that Adelus is running on
//             nprocsr --- number of processors assigned to a row
//             ncols   --- number of columns(=rows) for the matrix
//             nrhs    --- number of right hand sides
//
//             OUTPUT
//             my_rows  --- number of rows of the total matrix I own
//             my_cols  --- number of columns of the total matrix I own
//             my_first_row --- global number of my first row
//             my_first_col --- global number of my first column
//             my_rhs   --- number of right hand sides that I own
//             my_row   --- my subblock of the matrix
//             my_col   --- my subblock of the matrix

#include <mpi.h>
#include "Adelus_distribute.hpp"

namespace Adelus {

void distmat_( MPI_Comm comm,
               const int nprocsr,
               const int ncols,
               const int nrhs,
               int& my_rows,
               int& my_cols,
               int& my_first_row,
               int& my_first_col,
               int& my_rhs,
               int& my_row,
               int& my_col )
{

    int rank, nprocs, nprocs_col, nprocs_row, nrows;

    //  Determine who I am and the number of processors that are being used
    MPI_Comm_rank(comm, &rank) ;

    MPI_Comm_size(comm, &nprocs);

    nrows = ncols;

    nprocs_row = nprocsr;

    nprocs_col = nprocs/nprocsr;

    // Distribute the rows and columns

    my_row = rank / nprocsr;
    my_col = rank % nprocsr;

    //
    my_rows = nrows / nprocs_col;

    my_first_row = my_row * my_rows + 1;

    my_first_row = (my_row > (nrows%nprocs_col)) ? my_first_row + (nrows%nprocs_col) :
       my_first_row + my_row;

    if (my_row < (nrows%nprocs_col)) ++my_rows;

    //
    my_cols = ncols / nprocs_row;

    my_first_col = my_col * my_cols + 1;

    my_first_col = (my_col > (ncols%nprocs_row)) ? my_first_col + (ncols%nprocs_row) :
       my_first_col + my_col;

    if (my_col < (ncols%nprocs_row)) ++my_cols;

    // Distribute the RHS per processor

    my_rhs = nrhs / nprocs_row;
    if (my_col < (nrhs%nprocs_row)) ++my_rhs;

}

}//namespace Adelus
