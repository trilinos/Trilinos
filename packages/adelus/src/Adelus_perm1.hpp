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

#ifndef __ADELUS_PERM1_HPP__
#define __ADELUS_PERM1_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"

#define PERMTYPE ((1 << 5) + (1 << 4))

#define IBM_MPI_WRKAROUND

namespace Adelus {

#ifndef IBM_MPI_WRKAROUND

  //  Customized copy
  template<class XView, class YView>
  void zcopy_wr_local_index(int N, XView& X, YView& Y, int lidx) {
    Kokkos::parallel_for(Kokkos::RangePolicy<typename XView::device_type::execution_space>(0,N+1), KOKKOS_LAMBDA (const int i) {
      if (i<N)
        Y(i) = X(i);
      else
        Y(i) = (double)lidx + 0.1;
    });
  }

  //  Customized zcopy
  template<class XView, class YView>
  void zcopy_ld_local_index(int N, XView& X, YView& Y) {
    Kokkos::parallel_for(Kokkos::RangePolicy<typename XView::device_type::execution_space>(0,N), KOKKOS_LAMBDA (const int i) {
  #ifdef ADELUS_COMPLEX
      int lidx = (int)(X(N).real());
  #else
      int lidx = (int)(X(N));
  #endif
      Y(lidx,i) = X(i);
    });
  }

  //  Permutes -- unwraps the torus-wrap for the solution
  //              using the communication buffer
  template<class HandleType, class ZDView>
  inline
  void perm1_(HandleType& ahandle, ZDView& ZV) {

    MPI_Comm comm     = ahandle.get_comm();
    int me            = ahandle.get_myrank();
    int my_rhs_       = ahandle.get_my_rhs();
    int my_rows       = ahandle.get_my_rows();
    int nprocs_row    = ahandle.get_nprocs_row();
    int nprocs_col    = ahandle.get_nprocs_col();
    int nrows_matrix  = ahandle.get_nrows_matrix();
    int ncols_matrix  = ahandle.get_ncols_matrix();
    int my_first_row  = ahandle.get_my_first_row();
    int my_first_col  = ahandle.get_my_first_col();

    int i;

    int bytes;
    int dest;
    int type;

    int global_index;
    int local_index;


    int col_offset, row_offset;
    int ncols_proc1, ncols_proc2, nprocs_row1;
    int nrows_proc1, nrows_proc2, nprocs_col1;

    int change_nosend;
    int change_send;
    int next, next_s;
    int inc;

  #ifdef GET_TIMING
    double t2;
    double totalpermtime;
  #endif

    MPI_Request msgrequest;
    MPI_Status msgstatus;

    int ptr1_idx;

  #ifdef GET_TIMING
    t2 = MPI_Wtime();
  #endif

    typedef typename ZDView::value_type value_type;
    typedef typename ZDView::device_type::execution_space execution_space;
    typedef typename ZDView::device_type::memory_space memory_space;
    typedef Kokkos::View<value_type*, Kokkos::LayoutLeft, memory_space>  ViewVectorType;

    ViewVectorType temp_s     ( "temp_s", (my_rhs_ + 1) );
    ViewVectorType rhs_temp   ( "rhs_temp", (my_rows*(my_rhs_+1)) );
    ViewVectorType my_rhs_temp( "my_rhs_temp", (my_rows*(my_rhs_+1)) );

    if (my_rhs_ > 0) {
      ncols_proc1 = ncols_matrix/nprocs_row;
      ncols_proc2 = ncols_proc1;
      if (ncols_matrix%nprocs_row > 0) ncols_proc1++;
      nprocs_row1 = (ncols_matrix%nprocs_row);
      row_offset = ncols_proc1 * nprocs_row1;

      nrows_proc1 = nrows_matrix/nprocs_col;
      nrows_proc2 = nrows_proc1;
      if (nrows_matrix%nprocs_col > 0) nrows_proc1++;
      nprocs_col1 = (nrows_matrix%nprocs_col);
      col_offset = nrows_proc1 * nprocs_col1;

      ptr1_idx = 0;
      //change_count = 0;
      change_nosend = 0;
      next = 0;
      change_send = 0;
      next_s = 0;

  #ifdef PRINT_STATUS
      printf("Rank %i -- perm1_() Begin permutation, execution_space %s, memory_space %s\n",me,typeid(execution_space).name(),typeid(memory_space).name());
  #endif

      for (i=0; i<my_rows; i++) {
        global_index = my_first_row + i*nprocs_col;
        // break down global index using torus wrap in row direction
        dest = global_index%nprocs_row;
        local_index = global_index/nprocs_row;

        // rebuild global index using block in row direction
        if (dest < nprocs_row1) {
          global_index = dest*ncols_proc1 + local_index;
        } else {
          global_index = row_offset + (dest-nprocs_row1)*ncols_proc2
                         + local_index;
        }

        // break down global index using blocks in the column direction
        if (global_index < col_offset) {
          dest = global_index/nrows_proc1;
          local_index = global_index%nrows_proc1;
        } else {
          dest = (global_index - col_offset)/nrows_proc2 + nprocs_col1;
          local_index = (global_index - col_offset)%nrows_proc2;
        }

        dest = dest*nprocs_row + my_first_col;


        if ((local_index != i) || (dest != me)) {

          // Check if I need to send the data or just change position

          if( dest == me ) {

           auto sub_ZV          = subview(ZV,          ptr1_idx, Kokkos::ALL());
           auto sub_my_rhs_temp = subview(my_rhs_temp, Kokkos::make_pair(next, next+(my_rhs_ + 1)));
           zcopy_wr_local_index(my_rhs_, sub_ZV, sub_my_rhs_temp, local_index);

            change_nosend++;

            next = change_nosend * (my_rhs_ + 1);

          }

          if( dest != me ) {

            bytes = (my_rhs_ + 1)*sizeof(ADELUS_DATA_TYPE);

            MPI_Irecv( (char *)(reinterpret_cast<ADELUS_DATA_TYPE *>(rhs_temp.data())+next_s),bytes,MPI_CHAR,MPI_ANY_SOURCE,
                  MPI_ANY_TAG,comm,&msgrequest);

           auto sub_ZV = subview(ZV, ptr1_idx, Kokkos::ALL());
           zcopy_wr_local_index(my_rhs_, sub_ZV, temp_s, local_index);

           type = PERMTYPE+change_send;
           MPI_Send((char *)(reinterpret_cast<ADELUS_DATA_TYPE *>(temp_s.data())),bytes,MPI_CHAR,dest,
                   type,comm);
           change_send++;

           next_s = change_send * (my_rhs_+1);

           MPI_Wait(&msgrequest,&msgstatus);

          }

        }
        ptr1_idx++;
      }

      // Unpack changes from other processors

      next_s = 0;
      inc = 0;
      for (i = 0; i < change_send; i++) {
        auto sub_rhs_temp = subview(rhs_temp, Kokkos::make_pair(next_s, next_s+(my_rhs_ + 1)));
        zcopy_ld_local_index(my_rhs_, sub_rhs_temp, ZV);

        inc++;
        next_s = inc * (my_rhs_+1);
      }

      // Unpack my changes
      next = 0;
      inc = 0;
      for (i = 0; i < change_nosend; i++) {
        auto sub_my_rhs_temp = subview(my_rhs_temp, Kokkos::make_pair(next, next+(my_rhs_ + 1)));
        zcopy_ld_local_index(my_rhs_, sub_my_rhs_temp, ZV);
        inc++;
        next = inc * (my_rhs_+1);
      }

    }

  #ifdef GET_TIMING
    totalpermtime = MPI_Wtime() - t2;
  #endif
  #ifdef GET_TIMING
    showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Total time in perm", &totalpermtime);
  #endif
  }

#else//define IBM_MPI_WRKAROUND

  //  Customized copy
  template<class XView, class YView>
  void zcopy_wr_global_index(int N, XView& X, YView& Y, int idx) {
    Kokkos::parallel_for(Kokkos::RangePolicy<typename XView::device_type::execution_space>(0,N), KOKKOS_LAMBDA (const int i) {
      Y(idx,i) = X(i);
    });
  }

  //  Permutes -- unwraps the torus-wrap for the solution
  //              using the communication buffer
  template<class HandleType, class ZDView>
  inline
  void perm1_(HandleType& ahandle, ZDView& ZV) {

    MPI_Comm col_comm = ahandle.get_col_comm();
    int myrow         = ahandle.get_myrow();
    int my_rhs_       = ahandle.get_my_rhs();
    int my_rows       = ahandle.get_my_rows();
    int nprocs_row    = ahandle.get_nprocs_row();
    int nprocs_col    = ahandle.get_nprocs_col();
    int nrows_matrix  = ahandle.get_nrows_matrix();
    int ncols_matrix  = ahandle.get_ncols_matrix();
    int my_first_row  = ahandle.get_my_first_row();

    int i;

    int dest, global_index, local_index;

    int row_offset;
    int ncols_proc1, ncols_proc2, nprocs_row1;
    int ptr1_idx, myfirstrow;

  #ifdef GET_TIMING
  #if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    double t1, copyhostpinnedtime;
  #endif
    double t2;
    double totalpermtime;
  #endif

  #ifdef GET_TIMING
    t2 = MPI_Wtime();
  #endif

    typedef typename ZDView::value_type value_type;
  #ifdef PRINT_STATUS
    typedef typename ZDView::device_type::execution_space execution_space;
  #endif
    typedef typename ZDView::device_type::memory_space memory_space;
    typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewMatrixType;

  #ifdef ADELUS_HOST_PINNED_MEM_MPI
  #if defined(KOKKOS_ENABLE_CUDA)
    typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> View2DHostPinnType;//CudaHostPinnedSpace
  #elif defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::HIPHostPinnedSpace> View2DHostPinnType;//HIPHostPinnedSpace
  #endif
  #endif

    if (my_rhs_ > 0) {

      ViewMatrixType rhs_temp ( "rhs_temp", nrows_matrix, my_rhs_ );//allocate full-size RHS vectors
  #if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
      View2DHostPinnType h_rhs_temp( "h_rhs_temp", nrows_matrix, my_rhs_ );
  #endif

      Kokkos::deep_copy(rhs_temp, 0);//initialize with 0s

      ncols_proc1 = ncols_matrix/nprocs_row;
      ncols_proc2 = ncols_proc1;
      if (ncols_matrix%nprocs_row > 0) ncols_proc1++;
      nprocs_row1 = (ncols_matrix%nprocs_row);
      row_offset = ncols_proc1 * nprocs_row1;

      myfirstrow = myrow * (nrows_matrix / nprocs_col) + 1;
      myfirstrow = ( myrow > (nrows_matrix%nprocs_col) ) ? myfirstrow + (nrows_matrix%nprocs_col) :
                                                           myfirstrow + myrow;

      ptr1_idx = 0;

  #ifdef PRINT_STATUS
      printf("Rank %i -- perm1_() Begin permutation, execution_space %s, memory_space %s\n",ahandle.get_myrank(),typeid(execution_space).name(),typeid(memory_space).name());
  #endif

      for (i=0; i<my_rows; i++) {
        global_index = my_first_row + i*nprocs_col;

        /* break down global index using torus wrap in row direction */
        dest = global_index%nprocs_row;
        local_index = global_index/nprocs_row;

        /* rebuild global index using block in row direction */
        if (dest < nprocs_row1) {
          global_index = dest*ncols_proc1 + local_index;
        } else {
          global_index = row_offset + (dest-nprocs_row1)*ncols_proc2 + local_index;
        }

        auto sub_ZV = subview(ZV, ptr1_idx, Kokkos::ALL());
        zcopy_wr_global_index(my_rhs_, sub_ZV, rhs_temp, global_index);

        ptr1_idx++;
      }

  //Workaround for GPU-aware MPI issue on HIP: fence here to make sure "zcopy_wr_global_index"
  //                                           completes before MPI_Allreduce
  #if defined(KOKKOS_ENABLE_HIP) && !defined(ADELUS_HOST_PINNED_MEM_MPI)
      Kokkos::fence();
  #endif

  #if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  #ifdef GET_TIMING
      t1 = MPI_Wtime();
  #endif
      Kokkos::deep_copy(h_rhs_temp,rhs_temp);
  #ifdef GET_TIMING
      copyhostpinnedtime = (MPI_Wtime()-t1);
  #endif

      MPI_Allreduce( MPI_IN_PLACE, h_rhs_temp.data(), nrows_matrix*my_rhs_, ADELUS_MPI_DATA_TYPE, MPI_SUM, col_comm);

      Kokkos::deep_copy( subview(ZV, Kokkos::ALL(), Kokkos::make_pair(0, my_rhs_)),
                         subview(h_rhs_temp, Kokkos::make_pair(myfirstrow-1, myfirstrow-1+my_rows), Kokkos::ALL()) );
  #else //GPU-aware MPI
      MPI_Allreduce( MPI_IN_PLACE, rhs_temp.data(), nrows_matrix*my_rhs_, ADELUS_MPI_DATA_TYPE, MPI_SUM, col_comm);

      Kokkos::deep_copy( subview(ZV, Kokkos::ALL(), Kokkos::make_pair(0, my_rhs_)),
                         subview(rhs_temp, Kokkos::make_pair(myfirstrow-1, myfirstrow-1+my_rows), Kokkos::ALL()) );
  #endif
    }

  #ifdef GET_TIMING
    totalpermtime = MPI_Wtime() - t2;
  #endif

  #ifdef GET_TIMING
  #if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
             "Time to copy dev mem --> host pinned mem", &copyhostpinnedtime);
  #endif
    showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
             "Total time in perm", &totalpermtime);
  #endif
  }

#endif

}//namespace Adelus

#endif
