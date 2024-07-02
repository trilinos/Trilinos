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

#ifndef __ADELUS_DISTRIBUTE_HPP__
#define __ADELUS_DISTRIBUTE_HPP__

#include <mpi.h>

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
               int& my_col );

}//namespace Adelus

#endif
