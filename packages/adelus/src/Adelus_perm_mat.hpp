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

#ifndef __ADELUS_PERMMAT_HPP__
#define __ADELUS_PERMMAT_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"

//extern int me;	             // processor id information
//extern int nprocs_row;         // num of procs to which a row is assigned
//extern int nprocs_col;         // num of procs to which a col is assigned
//extern int nrows_matrix;       // number of rows in the matrix
//extern int ncols_matrix;       // number of cols in the matrix
//extern int my_rows;            // num of rows I own
//extern int my_cols;            // num of cols I own
//extern int myrow;
//extern int mycol;
//extern MPI_Comm col_comm;

namespace Adelus {

  template<class PViewType>
  inline 
  void exchange_pivots(PViewType& lpiv_view, PViewType& permute) {
  
    MPI_Status msgstatus;
    int rank_row,k_row,pivot_col;

    //  First gather the permutation vector to processor 0 in row_comm
    if (myrow == 0 || mycol == 0) {
      for (int k=0;k<=nrows_matrix-1;k++) {
        pivot_col = k%nprocs_row;
        k_row = k%nprocs_col;
        rank_row = k_row*nprocs_row;
        if (me == pivot_col) {
          int j=k/nprocs_row;
          MPI_Send(lpiv_view.data()+j,1,MPI_INT,rank_row,0,MPI_COMM_WORLD);
        }
        if (me == rank_row) {
          int i=k/nprocs_col;
          MPI_Recv(permute.data()+i,1,MPI_INT,pivot_col,0,MPI_COMM_WORLD,&msgstatus);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Broadcast to the rest of the processors in row_comm
    MPI_Bcast(permute.data(),my_rows,MPI_INT,0,row_comm);

  }// End of function exchange_pivots
  
  template<class ZViewType, class PViewType>
  inline
  void permute_mat(ZViewType& Z, PViewType& lpiv_view, PViewType& permute) {
    //TODO: add host pinned memory support
    using value_type  = typename ZViewType::value_type;
#ifndef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    using execution_space = typename ZViewType::device_type::execution_space ;
    using memory_space    = typename ZViewType::device_type::memory_space ;
    using ViewVectorType  = Kokkos::View<value_type*, Kokkos::LayoutLeft, memory_space>;
#ifdef PRINT_STATUS
  printf("Rank %i -- permute_mat() Begin permute mat with myrow %d, mycol %d, nprocs_row %d, nprocs_col %d, nrows_matrix %d, ncols_matrix %d, my_rows %d, my_cols %d, my_rhs %d, nrhs %d, value_type %s, execution_space %s, memory_space %s\n", me, myrow, mycol, nprocs_row, nprocs_col, nrows_matrix, ncols_matrix, my_rows, my_cols, my_rhs, nrhs, typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif
#endif

    MPI_Status msgstatus;
  
    int pivot_row, k_row;
#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    value_type tmpr, tmps;
#else
    ViewVectorType tmpr( "tmpr", Z.extent(1) );
    ViewVectorType tmps( "tmps", Z.extent(1) );
#endif

#ifdef GET_TIMING
   double exchpivtime,permutemattime,t1;

   t1 = MPI_Wtime();
#endif

    exchange_pivots(lpiv_view, permute);

#ifdef GET_TIMING
    exchpivtime = MPI_Wtime()-t1;

    t1 = MPI_Wtime();
#endif

#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    for (int j=0;j<=my_cols-1;j++) {
      int J=j*nprocs_row+mycol; // global column index
      for (int k=J+1;k<=nrows_matrix-1;k++) {
        k_row=k%nprocs_col;
        if (myrow==k_row)
          pivot_row=permute(k/nprocs_col);
        MPI_Bcast(&pivot_row,1,MPI_INT,k_row,col_comm);
        if (k != pivot_row) {
          if (myrow == k_row) {
            tmps = Z(k/nprocs_col, J/nprocs_row);
            MPI_Send((char *)(&tmps),sizeof(value_type),MPI_CHAR,pivot_row%nprocs_col,2,col_comm);
          }
          if (myrow == pivot_row%nprocs_col) {
            tmps = Z(pivot_row/nprocs_col, J/nprocs_row);
            MPI_Send((char *)(&tmps),sizeof(value_type),MPI_CHAR,k_row,3,col_comm);
          }
          if (myrow == k_row) {
            MPI_Recv((char *)(&tmpr),sizeof(value_type),MPI_CHAR,pivot_row%nprocs_col,3,col_comm,&msgstatus);
            Z(k/nprocs_col, J/nprocs_row) = tmpr;
          }
          if (myrow == pivot_row%nprocs_col) {
            MPI_Recv((char *)(&tmpr),sizeof(value_type),MPI_CHAR,k_row,2,col_comm,&msgstatus);
            Z(pivot_row/nprocs_col, J/nprocs_row)  = tmpr;
          }
        }// End of if (k != pivot_row)
      }// End of for (k=J+1;k<=nrows_matrix-1;k++)
    }// End of for (j=0;j<=my_cols-1;j++)
#else
    for (int k = 1 + mycol; k <= nrows_matrix - 1; k++) {
      int max_gcol_k=k-1; // max. global column index in the k row
      int max_lcol_k=0;   // max. local column index in the k row
      k_row=k%nprocs_col; // mesh row id (in the MPI process mesh) of the process that holds k

      if (myrow==k_row) pivot_row = permute(k/nprocs_col);
      MPI_Bcast(&pivot_row,1,MPI_INT,k_row,col_comm);

      int max_gcol_pivot=pivot_row-1;          // max. global column index in the pivot row
      int max_lcol_pivot=0;                    // max. local column index in the pivot row
      int pivot_row_pid = pivot_row%nprocs_col;// mesh row id (in the MPI process mesh) of the process that holds pivot_row

      //Find max. local column index in the k row that covers the lower triangular part
      if ( mycol <= max_gcol_k%nprocs_row)
        max_lcol_k = max_gcol_k/nprocs_row;
      else
        max_lcol_k = max_gcol_k/nprocs_row - 1;//one element less

      //Find max. local column index in the pivot row that covers the lower triangular part
      if ( mycol <= max_gcol_pivot%nprocs_row)
        max_lcol_pivot = max_gcol_pivot/nprocs_row;
      else
        max_lcol_pivot = max_gcol_pivot/nprocs_row - 1;//one element less

      //Find the number of columns needs to be exchanged
      int min_len = std::min(max_lcol_k,max_lcol_pivot) + 1;

      if (k != pivot_row) {//k row is differrent from pivot_row, i.e. needs permutation
        if (k_row == pivot_row_pid) {//pivot row is in the same rank
          if (myrow == k_row) {//I am the right process to do permutation
            int curr_lrid = k/nprocs_col;
            int piv_lrid  = pivot_row/nprocs_col;
            Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,min_len), KOKKOS_LAMBDA (const int i) {
              value_type tmp = Z(curr_lrid,i);
              Z(curr_lrid,i) = Z(piv_lrid,i);
              Z(piv_lrid,i)  = tmp;
            });
            Kokkos::fence();
          }
        }
        else {//k row and pivot row are in different processes (rank)
          if (myrow == k_row) {//I am holding k row
            int curr_lrid = k/nprocs_col;
            Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,min_len), KOKKOS_LAMBDA (const int i) {
              tmps(i) = Z(curr_lrid,i);
            });
            Kokkos::fence();

            MPI_Send(reinterpret_cast<char *>(tmps.data()),min_len*sizeof(value_type),MPI_CHAR,pivot_row_pid,2,col_comm);
            MPI_Recv(reinterpret_cast<char *>(tmpr.data()),min_len*sizeof(value_type),MPI_CHAR,pivot_row_pid,3,col_comm,&msgstatus);

            Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,min_len), KOKKOS_LAMBDA (const int i) {
              Z(curr_lrid,i) = tmpr(i);
            });
            Kokkos::fence();
          }
          if (myrow == pivot_row_pid) {//I am holding the pivot row
            int piv_lrid = pivot_row/nprocs_col;
            Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,min_len), KOKKOS_LAMBDA (const int i) {
              tmps(i) = Z(piv_lrid,i);
            });
            Kokkos::fence();

            MPI_Recv(reinterpret_cast<char *>(tmpr.data()),min_len*sizeof(value_type),MPI_CHAR,k_row,2,col_comm,&msgstatus);
            MPI_Send(reinterpret_cast<char *>(tmps.data()),min_len*sizeof(value_type),MPI_CHAR,k_row,3,col_comm);

            Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,min_len), KOKKOS_LAMBDA (const int i) {
              Z(piv_lrid,i) = tmpr(i);
            });
            Kokkos::fence();
          }
        }//End of k row and pivot row are in different processes (rank)
      }// End of if (k != pivot_row)
    }// End of for (int k=1+mycol;k<=nrows_matrix-1;k++) {
#endif

#ifdef GET_TIMING
    permutemattime = MPI_Wtime()-t1;

    showtime("Time to exchange pivot information",&exchpivtime);
    showtime("Time to permute matrix",&permutemattime);    
#endif
  }// End of function permute_mat

}//namespace Adelus

#endif
