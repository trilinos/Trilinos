/*
//@HEADER
// ************************************************************************
//
//                        Adelus v. 1.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// 3. Neither the name of NTESS nor the names of the contributors may be
// used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL NTESS OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Vinh Dang (vqdang@sandia.gov)
//                    Joseph Kotulski (jdkotul@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
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

void distmat_( int *nprocsr,
               int *ncols,
               int *nrhs_,
               int *my_rows_,
               int *my_cols_,
               int *my_first_row_,
               int *my_first_col_,
               int *my_rhs_,
               int *my_row,
               int *my_col )
{

    int rank,nprocs;
    int nprocs_col_, nrows;
    int nprocs_row_;

    //  Determine who I am and the number of processors that are being used

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    nrows = *ncols;

    nprocs_row_ = *nprocsr;

    nprocs_col_ = nprocs/(*nprocsr) ;

    // Distribute the rows and columns

    *my_row = rank/(*nprocsr);
    *my_col = rank %(nprocs_row_);


    *my_rows_ = nrows / nprocs_col_;

    *my_first_row_ = (*my_row)*(*my_rows_) +1;

    *my_first_row_ = ((*my_row) > (nrows%nprocs_col_)) ? *my_first_row_ + (nrows%nprocs_col_) :
       *my_first_row_ + (*my_row);

    if (*my_row < nrows % nprocs_col_)
        ++(*my_rows_);

    *my_cols_ = nrows / nprocs_row_;

    *my_first_col_ = (*my_col)*(*my_cols_) + 1;

    *my_first_col_ = ((*my_col) > (nrows%nprocs_row_)) ? *my_first_col_ + (nrows%nprocs_row_) :
       *my_first_col_ + (*my_col);

    *my_cols_ = *ncols / *nprocsr;
    if (*my_col < *ncols % (*nprocsr))
        ++(*my_cols_);

    // Distribute the RHS per processor

    *my_rhs_ = *nrhs_ / *nprocsr;
    if (*my_col < *nrhs_ % (*nprocsr)) ++(*my_rhs_);

}

}//namespace Adelus
