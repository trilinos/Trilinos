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

namespace Adelus {

  template<class HandleType, class ZViewType, class PViewType>
  inline
  void permute_rhs(HandleType& ahandle, ZViewType& RHS, PViewType& permute) {
    using value_type      = typename ZViewType::value_type;
    using execution_space = typename ZViewType::device_type::execution_space ;
    using memory_space    = typename ZViewType::device_type::memory_space ;
    using ViewVectorType  = Kokkos::View<value_type*, Kokkos::LayoutLeft, memory_space>;
#ifdef ADELUS_HOST_PINNED_MEM_MPI
  #if defined(KOKKOS_ENABLE_CUDA)
    using ViewVectorHostPinnType = Kokkos::View<value_type*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace>;//CudaHostPinnedSpace
  #elif defined(KOKKOS_ENABLE_HIP)
    using ViewVectorHostPinnType = Kokkos::View<value_type*, Kokkos::LayoutLeft, Kokkos::HIPHostPinnedSpace>;//HIPHostPinnedSpace
  #endif
#endif

    MPI_Comm col_comm = ahandle.get_col_comm();
    int myrow         = ahandle.get_myrow();
    int nprocs_col    = ahandle.get_nprocs_col();
    int nrows_matrix  = ahandle.get_nrows_matrix();

    int pivot_row, k_row;
    ViewVectorType tmpr( "tmpr", RHS.extent(1) );
    ViewVectorType tmps( "tmps", RHS.extent(1) );
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    ViewVectorHostPinnType h_tmpr( "h_tmpr", RHS.extent(1) );
    ViewVectorHostPinnType h_tmps( "h_tmps", RHS.extent(1) );
#endif

    MPI_Status msgstatus;

    //TODO: try this later
    //MPI_Datatype strided_vec_type;
    //int strided_vec_nblocks  = RHS.extent(1);
    //int strided_vec_blocklen = 1;
    //int strided_vec_stride   = RHS.extent(0);
    //MPI_Type_vector( strided_vec_nblocks, strided_vec_blocklen, strided_vec_stride,
    //                 ADELUS_MPI_DATA_TYPE, &strided_vec_type);
    //MPI_Type_commit(&strided_vec_type);

#ifdef GET_TIMING
   double permuterhstime,t1;

   t1 = MPI_Wtime();
#endif

    for (int k=0;k<=nrows_matrix-2;k++) {
      k_row=k%nprocs_col;

      if (ahandle.get_my_rhs() > 0) {
        if (myrow==k_row) pivot_row = static_cast<int>(permute(k/nprocs_col));
        MPI_Bcast(&pivot_row,1,MPI_INT,k_row,col_comm);
        int pivot_row_pid = pivot_row%nprocs_col;

        if (k != pivot_row) {
          if (k_row == pivot_row_pid) {//pivot row is in the same rank
            if (myrow == k_row) {
              int curr_lrid = k/nprocs_col;
              int piv_lrid  = pivot_row/nprocs_col;
              Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
                value_type tmp   = RHS(curr_lrid,i);
                RHS(curr_lrid,i) = RHS(piv_lrid,i);
                RHS(piv_lrid,i)  = tmp;
              });
              Kokkos::fence();
            }
          }
          else {//pivot row is is a different rank
            if (myrow == k_row) {
              int curr_lrid = k/nprocs_col;
              Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
                tmps(i) = RHS(curr_lrid,i);
              });

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
              Kokkos::deep_copy(h_tmps,tmps);
              MPI_Send(reinterpret_cast<char *>(h_tmps.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,pivot_row_pid,2,col_comm);
              MPI_Recv(reinterpret_cast<char *>(h_tmpr.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,pivot_row_pid,3,col_comm,&msgstatus);
              Kokkos::deep_copy(tmpr,h_tmpr);
#else //GPU-aware MPI
              Kokkos::fence();

              MPI_Send(reinterpret_cast<char *>(tmps.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,pivot_row_pid,2,col_comm);
              MPI_Recv(reinterpret_cast<char *>(tmpr.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,pivot_row_pid,3,col_comm,&msgstatus);
#endif

              Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
                RHS(curr_lrid,i) = tmpr(i);
              });
              Kokkos::fence();
            }
            if (myrow == pivot_row_pid) {
              int piv_lrid = pivot_row/nprocs_col;
              Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
                tmps(i) = RHS(piv_lrid,i);
              });

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
              Kokkos::deep_copy(h_tmps,tmps);
              MPI_Recv(reinterpret_cast<char *>(h_tmpr.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,k_row,2,col_comm,&msgstatus);
              MPI_Send(reinterpret_cast<char *>(h_tmps.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,k_row,3,col_comm);
              Kokkos::deep_copy(tmpr,h_tmpr);
#else // GPU-aware MPI
              Kokkos::fence();

              MPI_Recv(reinterpret_cast<char *>(tmpr.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,k_row,2,col_comm,&msgstatus);
              MPI_Send(reinterpret_cast<char *>(tmps.data()),RHS.extent(1)*sizeof(value_type),MPI_CHAR,k_row,3,col_comm);
#endif

              Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
                RHS(piv_lrid,i) = tmpr(i);
              });
              Kokkos::fence();
            }
          }//End of pivot row is is a different rank
        }// End of if (k != pivot_row)

      }// End of if (my_num_rhs > 0)

    }// End of for (k=0;k<=nrows_matrix-2;k++)

#ifdef GET_TIMING
    permuterhstime = MPI_Wtime()-t1;

    showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
             "Time to permute rhs", &permuterhstime);
#endif
  }// End of function permute_rhs

}//namespace Adelus

#endif
