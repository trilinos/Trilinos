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

#ifndef __ADELUS_PERMRHS_HPP__
#define __ADELUS_PERMRHS_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"

//extern int me;	               // processor id information
//extern int nprocs_row;         // num of procs to which a row is assigned
//extern int nprocs_col;         // num of procs to which a col is assigned
//extern int nrows_matrix;       // number of rows in the matrix
//extern int ncols_matrix;       // number of cols in the matrix
//extern int my_rows;            // num of rows I own
//extern int my_cols;            // num of cols I own
//extern int my_rhs;             // num of right hand side I own
//extern int myrow;
//extern int mycol;
//extern MPI_Comm col_comm;

namespace Adelus {
  
  template<class ZViewType, class PViewType>
  inline
  void permute_rhs(ZViewType& RHS, PViewType& permute) {
    //NOTE: Currently assume that a single RHS resides in host memory
    using value_type  = typename ZViewType::value_type;

    MPI_Status msgstatus;
  
    int pivot_row, k_row, rhs_col;
    value_type tmpr, tmps;

#ifdef GET_TIMING
   double permuterhstime,t1;

   t1 = MPI_Wtime();
#endif

    rhs_col = 0;
    for (int k=0;k<=nrows_matrix-2;k++) {
      k_row=k%nprocs_col;
      if (mycol == rhs_col) {
        if (myrow==k_row)
          pivot_row=permute(k/nprocs_col);
        MPI_Bcast(&pivot_row,1,MPI_INT,k_row,col_comm);
        if (k != pivot_row) {
          if (myrow == k_row) {
            tmps = RHS(k/nprocs_col,0);
            MPI_Send((char *)(&tmps),sizeof(value_type),MPI_CHAR,pivot_row%nprocs_col,2,col_comm);
          }
          if (myrow == pivot_row%nprocs_col) {
            tmps = RHS(pivot_row/nprocs_col,0);
            MPI_Send((char *)(&tmps),sizeof(value_type),MPI_CHAR,k_row,3,col_comm);
          }
          if (myrow == k_row) {
            MPI_Recv((char *)(&tmpr),sizeof(value_type),MPI_CHAR,pivot_row%nprocs_col,3,col_comm,&msgstatus);
            RHS(k/nprocs_col,0) = tmpr;
          }
          if (myrow == pivot_row%nprocs_col) {
            MPI_Recv((char *)(&tmpr),sizeof(value_type),MPI_CHAR,k_row,2,col_comm,&msgstatus);
            RHS(pivot_row/nprocs_col,0)  = tmpr;
          }
        }// End of if (k != pivot_row)
      }
    }// End of for (k=0;k<=nrows_matrix-2;k++)

#ifdef GET_TIMING
    permuterhstime = MPI_Wtime()-t1;

    showtime("Time to permute rhs",&permuterhstime);    
#endif
  }// End of function permute_rhs

}//namespace Adelus

#endif
