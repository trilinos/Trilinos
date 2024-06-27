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

#ifndef __ADELUS_FORWARD_HPP__
#define __ADELUS_FORWARD_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosBlas3_gemm.hpp"

namespace Adelus {

template<class HandleType, class ZViewType, class RHSViewType>
inline
void forward(HandleType& ahandle, ZViewType& Z, RHSViewType& RHS)
{
  using value_type      = typename ZViewType::value_type ;
  using execution_space = typename ZViewType::device_type::execution_space ;
  using memory_space    = typename ZViewType::device_type::memory_space ;
  using ViewMatrixType  =  Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space>;
#ifdef ADELUS_HOST_PINNED_MEM_MPI
  #if defined(KOKKOS_ENABLE_CUDA)
    using ViewMatrixHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace>;//CudaHostPinnedSpace
  #elif defined(KOKKOS_ENABLE_HIP)
    using ViewMatrixHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::HIPHostPinnedSpace>;//HIPHostPinnedSpace
  #endif
#endif

  MPI_Comm row_comm = ahandle.get_row_comm();
  MPI_Comm col_comm = ahandle.get_col_comm();
  int myrow         = ahandle.get_myrow();
  int nprocs_row    = ahandle.get_nprocs_row();
  int nprocs_col    = ahandle.get_nprocs_col();
  int nrows_matrix  = ahandle.get_nrows_matrix();
  int my_rows       = ahandle.get_my_rows();

  int k_row;       // torus-wrap row corresponding to kth global row
  int k_col;       // torus-wrap column corresponding to kth global col
  int istart;      // Starting row index for pivot column
  int count_row;   // dummy index
  value_type d_one     = static_cast<value_type>( 1.0);
  value_type d_min_one = static_cast<value_type>(-1.0);

  ViewMatrixType piv_col( "piv_col", my_rows, 1 ); // portion of pivot column I am sending
  ViewMatrixType ck( "ck", 1, RHS.extent(1) ); // rhs corresponding to current column of the backsubstitution
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    ViewMatrixHostPinnType h_piv_col( "h_piv_col", my_rows, 1 );
    ViewMatrixHostPinnType h_ck( "h_ck", 1, RHS.extent(1) );
#endif

#ifdef PRINT_STATUS
  printf("Rank %i -- forward() Begin forward solve with myrow %d, nprocs_row %d, nprocs_col %d, nrows_matrix %d, ncols_matrix %d, my_rows %d, my_cols %d, my_rhs %d, nrhs %d, value_type %s, execution_space %s, memory_space %s\n", ahandle.get_myrank(), myrow, nprocs_row, nprocs_col, nrows_matrix, ahandle.get_ncols_matrix(), my_rows, ahandle.get_my_cols(), ahandle.get_my_rhs(), ahandle.get_nrhs(), typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

#ifdef GET_TIMING
  double t1, fwdsolvetime;
  t1 = MPI_Wtime();
#endif

  // Perform the Forward Substitution
  for (int k=0; k<= nrows_matrix-2; k++) {
    k_row=k%nprocs_col;
    k_col=k%nprocs_row;
    istart = (k+1-myrow)/nprocs_col;
    if (istart * nprocs_col < k+1-myrow) istart++;

    if (istart < my_rows) {
      Kokkos::deep_copy( subview(piv_col, Kokkos::make_pair(0, my_rows - istart), 0),
                         subview(Z, Kokkos::make_pair(istart, my_rows), k/nprocs_row) );
    }
    count_row = my_rows - istart;

    //Note: replace MPI_Send/MPI_Irecv with MPI_Bcast
    //      Rank k_col broadcasts the pivot_col to all
    //      other ranks in the row_comm
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    Kokkos::deep_copy(h_piv_col,piv_col);
    MPI_Bcast(reinterpret_cast<char *>(h_piv_col.data()), count_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_col, row_comm);
    Kokkos::deep_copy(piv_col,h_piv_col);
#else //GPU-aware MPI
    MPI_Bcast(reinterpret_cast<char *>(piv_col.data()), count_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_col, row_comm);
#endif

    if (ahandle.get_my_rhs() > 0) {
      //ck = RHS(k/nprocs_col,0);
      //MPI_Bcast((char *)(&ck),sizeof(ADELUS_DATA_TYPE),MPI_CHAR,k_row,col_comm);
      //count_row=0;
      //printf("Point 2: k %d, istart %d, my_rows %d\n", k, istart, my_rows);
      //for (int i=istart;i<=my_rows-1;i++) {
      //  RHS(i,0) = RHS(i,0) - piv_col(count_row) * ck;
      //  count_row++;
      //}
      int curr_lrid = k/nprocs_col;//note: nprocs_col (global var) cannot be read in a device function
      Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,RHS.extent(1)), KOKKOS_LAMBDA (const int i) {
        ck(0,i) = RHS(curr_lrid,i);
      });

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
      Kokkos::deep_copy(h_ck,ck);
      MPI_Bcast(reinterpret_cast<char *>(h_ck.data()), RHS.extent(1)*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_row, col_comm);
      Kokkos::deep_copy(ck,h_ck);
#else //GPU-aware MPI
      Kokkos::fence();
      MPI_Bcast(reinterpret_cast<char *>(ck.data()), RHS.extent(1)*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_row, col_comm);
#endif

      auto sub_pivot_col = subview(piv_col, Kokkos::make_pair(0, my_rows - istart), Kokkos::ALL());
      auto sub_rhs       = subview(RHS, Kokkos::make_pair(istart, my_rows), Kokkos::ALL());
      if (istart < my_rows) {
        KokkosBlas::gemm("N", "N", d_min_one, sub_pivot_col, ck, d_one, sub_rhs);
      }
    }
    MPI_Barrier(ahandle.get_comm());
  }// end of for (k=0; k<= nrows_matrix-2; k++)

#ifdef GET_TIMING
  fwdsolvetime = MPI_Wtime() - t1;
  showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
           "Total time in forward solve", &fwdsolvetime);
#endif
}

}//namespace Adelus

#endif
