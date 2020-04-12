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

#ifndef __ADELUS_DISTRIBUTE_HPP__
#define __ADELUS_DISTRIBUTE_HPP__

#include <mpi.h>

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
               int *my_col );

}//namespace Adelus

#endif
